"""Worker function for the sweep (Tab 2), runs in separate processes via
ProcessPoolExecutor. Needs to live in an importable module (not inside the
Streamlit page script) to be picklable by multiprocessing.

Only calls scf.* — all physics lives there (R1)."""

import sys
from pathlib import Path

_DASHBOARD_DIR = Path(__file__).resolve().parent
_SCF_DIR = _DASHBOARD_DIR.parent / "scf"
for _p in (str(_SCF_DIR), str(_DASHBOARD_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# Certified methodology (docs/teoria.md Sec 6.2b-c): the self-consistent
# toroidal branch needs continuation from K=0 (cold start does not converge
# at the K values that matter, e.g. K>=5e-3) and a domain sized until the
# star's polar radius comfortably clears the box edge (a purely-toroidal
# field deforms PROLATE, so it is the POLE that overflows first -- checking
# only the equator, as an oblate/rotation-only check would, misses this
# family entirely). frac_pol<=0.2 was the tightened criterion adopted after
# frac_pol=0.4 was found still VE-failing at frac_pol=0.83 (the criterion
# used before this fix).
_FRAC_MAX = 0.2
_INITIAL_DOMAIN_MULT = 10.0
_MAX_DOMAIN_GROWTHS = 4
_K_STEP = 1e-3  # increment validated by the Sec 6.2c continuation studies
# Delta r / R_guess validated in Sec 6.2 (nr=986 at domain=10xR_guess).
# NOT derived from the UI's Nr slider (default 129, dr/R_guess~0.08) --
# that resolution was never validated for this branch and silently gives
# an uncertified VE (found while porting this methodology: dr/R_guess~0.08
# gave VE~4.6e-3 at a point independently confirmed certified at ~4e-4
# with the finer, validated ratio below).
_DR_OVER_RGUESS = 3.824e5 / 3.766e7


def _solve_toroidal_certified(rho_c, R_guess, K_tor, m_tor_sc, rotation, mu_e,
                               Nr_base, Ntheta, lmax, tol, max_iter):
    """Continuation in K from 0 to K_tor, growing the domain (Delta r held
    fixed, Nr scaled up with it -- Experiment B, docs/teoria.md Sec 6.2b)
    until frac_pol/frac_eq <= _FRAC_MAX or the growth budget runs out.

    Returns (result, r, theta, overflow) -- result is None if the SCF
    itself failed to converge at some step along the continuation path
    (a real non-convergence, not a domain-sizing problem, reported as a
    failure same as before); overflow is the domain_overflow_check() dict
    for the final attempt either way (best-effort even if _FRAC_MAX was
    never reached, so the caller can see how far off it was rather than
    silently trusting an under-sized box)."""
    import numpy as np
    import scf as scf_mod
    import diagnostics as diag
    from terms.toroidal_sc import ToroidalSC

    n_steps = max(1, int(np.ceil(K_tor / _K_STEP))) if K_tor > 0 else 0
    k_path = list(np.linspace(0.0, K_tor, n_steps + 1))

    dr_target = _DR_OVER_RGUESS * R_guess
    domain_mult = _INITIAL_DOMAIN_MULT
    result = overflow = r = theta = None

    for _attempt in range(_MAX_DOMAIN_GROWTHS + 1):
        domain = domain_mult * R_guess
        Nr = max(int(round(domain / dr_target)) + 1, Nr_base)
        r = np.linspace(0, domain, Nr)
        theta = np.linspace(0, np.pi, Ntheta)
        rho_seed = scf_mod.initial_guess(r, theta, rho_c, R_guess)

        converged_path = True
        for K in k_path:
            toroidal = ToroidalSC(K=K, m=m_tor_sc) if K > 0 else None
            result = scf_mod.hachisu_scf(rho_seed, r, theta, rho_c, rotation=rotation,
                                          toroidal=toroidal, mu_e=mu_e, lmax=lmax,
                                          tol=tol, max_iter=int(max_iter))
            if not result["converged"]:
                converged_path = False
                break
            rho_seed = result["rho"]

        if not converged_path:
            return None, r, theta, None

        R_eq, R_pol = diag.equatorial_polar_radii(result["H"], r, theta)
        overflow = diag.domain_overflow_check(R_eq, R_pol, r[-1], tol=0.1)
        if max(overflow["frac_eq"], overflow["frac_pol"]) <= _FRAC_MAX:
            return result, r, theta, overflow
        domain_mult *= 2.0

    return result, r, theta, overflow


def run_one(params: dict) -> dict:
    """Runs one SCF point (rho_c + whatever field/rotation params are
    active) and returns scalars + fields (to save with store.save_run).
    Does not converge -> returns converged=False, no fields.

    field_mode: "none" | "poloidal" | "toroidal (self-consistent)"
    rotation_mode: "none" | "rigid" | "differential"
    Same mode strings as dashboard/pages/1_equilibrium.py, so a params
    dict built by either page means the same thing.
    """
    import numpy as np
    import scf as scf_mod
    import diagnostics as diag
    import toroidal as tor
    import units
    import eos
    from terms.poloidal import Poloidal
    from terms.rotation import Rotation
    from terms.toroidal_sc import ToroidalSC

    rho_c = params["rho_c"]
    k0 = params.get("k0", 0.0)
    K_tor = params.get("K_tor", 0.0)
    m_tor_sc = params.get("m_tor_sc", 1.0)
    Omega_c = params.get("Omega_c", 0.0)
    A_over_Req = params.get("A_over_Req", 0.0)
    field_mode = params.get("field_mode", "poloidal" if k0 != 0.0 else "none")
    rotation_mode = params.get("rotation_mode", "none")
    mu_e = params.get("mu_e", 2.0)
    R_guess = params["R_guess"]
    Nr = params.get("Nr", 129)
    Ntheta = params.get("Ntheta", 129)
    lmax = params.get("lmax", 16)
    tol = params.get("tol", 1e-6)
    max_iter = params.get("max_iter", 200)

    rotation = None
    if rotation_mode == "rigid":
        rotation = Rotation(Omega_c=Omega_c, A=float("inf"))
    elif rotation_mode == "differential":
        rotation = Rotation(Omega_c=Omega_c, A=A_over_Req * R_guess)

    if field_mode == "toroidal (self-consistent)":
        result, r, theta, overflow = _solve_toroidal_certified(
            rho_c, R_guess, K_tor, m_tor_sc, rotation, mu_e, Nr, Ntheta, lmax, tol, max_iter)
        if result is None:
            return {"converged": False, "rho_c": rho_c, "k0": k0, "K_tor": K_tor, "Omega_c": Omega_c}
        poloidal = None
        toroidal_sc = ToroidalSC(K=K_tor, m=m_tor_sc) if K_tor > 0 else None
    else:
        r = np.linspace(0, 1.3 * R_guess, Nr)
        theta = np.linspace(0, np.pi, Ntheta)
        rho0 = scf_mod.initial_guess(r, theta, rho_c, R_guess)
        poloidal = Poloidal(k0=k0, lmax=lmax) if field_mode == "poloidal" and k0 != 0.0 else None
        toroidal_sc = None

        result = scf_mod.hachisu_scf(rho0, r, theta, rho_c, rotation=rotation, poloidal=poloidal,
                                      toroidal=toroidal_sc, mu_e=mu_e, lmax=lmax,
                                      tol=tol, max_iter=int(max_iter))
        if not result["converged"]:
            return {"converged": False, "rho_c": rho_c, "k0": k0, "K_tor": K_tor, "Omega_c": Omega_c}

        # Domain sizing here is still the historical 1.3xR_guess default,
        # NOT the certified frac_pol<=0.2 methodology (that was only
        # validated for the self-consistent toroidal branch, docs/teoria.md
        # Sec 6.2 -- porting it to poloidal/rotation-only points is
        # unvalidated scope creep, not done here). The check below is
        # reported for visibility only; a flagged overflow on this path is
        # NOT auto-corrected the way the toroidal branch above is.
        R_eq, R_pol = diag.equatorial_polar_radii(result["H"], r, theta)
        overflow = diag.domain_overflow_check(R_eq, R_pol, r[-1], tol=0.1)

    rho, Phi, u, H = result["rho"], result["Phi"], result["u"], result["H"]
    ve = diag.virial_error_terms(rho, Phi, H, r, theta, mu_e,
                                  rotation=rotation, poloidal=poloidal, toroidal=toroidal_sc)
    Br, Bth, Bphi = ve["Br"], ve["Btheta"], ve["Bphi"]
    VE, W, T = ve["VE"], ve["W"], ve["T"]
    M = scf_mod.total_mass(rho, r, theta)
    R_eq, R_pol = diag.equatorial_polar_radii(H, r, theta)

    Bpol_grid = np.sqrt(Br**2 + Bth**2)
    B_pol_max = float(np.max(Bpol_grid))
    B_tor_max = float(np.max(np.abs(Bphi)))
    ratio_energy, ratio_amp = tor.bt_bp_ratios(Br, Bth, Bphi, r, theta)
    T_over_W = T / abs(W) if W != 0 else float("nan")
    mass_loss_ratio = diag.equatorial_mass_loss_ratio(Phi, rotation, r, theta, R_eq)

    # dipolarity only means something with a poloidal field (Bpol==0
    # identically otherwise, e.g. pure toroidal or unmagnetized)
    if poloidal is not None:
        dip = diag.surface_dipolarity(Bpol_grid, H, r, theta)
        B_polo, B_eq_surf, dipolarity = dip["B_pole"], dip["B_eq"], dip["dipolarity"]
    else:
        B_polo = B_eq_surf = dipolarity = float("nan")

    # neutronization validity gate (item 3) — flagged, not filtered: the
    # dashboard plots this point distinctly instead of silently dropping it
    rho_c_neutronization = eos.neutronization_threshold_rho_c(mu_e)
    valid_rho_c = rho_c < rho_c_neutronization

    scalars = {
        "M/M_sun": units.g_to_msun(M),
        "R_eq (km)": units.cm_to_km(R_eq),
        "R_pol (km)": units.cm_to_km(R_pol),
        "R_pol/R_eq": R_pol / R_eq if R_eq > 0 else float("nan"),
        "W (erg)": W,
        "E_mag (erg)": ve["E_mag"],
        "E_mag/|W|": ve["E_mag"] / abs(W) if W != 0 else float("nan"),
        "B_pol,max (G)": B_pol_max,
        "B_tor,max (G)": B_tor_max,
        "B_polo (G)": B_polo,
        "B_eq (G)": B_eq_surf,
        "dipolarity": dipolarity,
        "Bt/Bp (energy)": ratio_energy,
        "Bt/Bp (amplitude)": ratio_amp,
        "T (erg)": T,
        "T/|W|": T_over_W,
        "Omega_c (rad/s)": Omega_c,
        "equatorial mass-loss ratio": mass_loss_ratio,
        "rho_c_valid": bool(valid_rho_c),
        "VE": VE,
        "frac_eq": overflow["frac_eq"],
        "frac_pol": overflow["frac_pol"],
        "domain_overflow": bool(overflow["overflow"]),
    }
    fields = {"rho": rho, "Phi": Phi, "u": u, "H": H, "Bphi": Bphi, "r": r, "theta": theta}
    return {"converged": True, "rho_c": rho_c, "k0": k0, "K_tor": K_tor, "Omega_c": Omega_c,
            "scalars": scalars, "fields": fields}
