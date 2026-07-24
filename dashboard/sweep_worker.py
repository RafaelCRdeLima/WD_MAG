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

    r = np.linspace(0, 1.3 * R_guess, Nr)
    theta = np.linspace(0, np.pi, Ntheta)
    rho0 = scf_mod.initial_guess(r, theta, rho_c, R_guess)

    poloidal = Poloidal(k0=k0, lmax=lmax) if field_mode == "poloidal" and k0 != 0.0 else None
    toroidal_sc = (ToroidalSC(K=K_tor, m=m_tor_sc)
                   if field_mode == "toroidal (self-consistent)" and K_tor > 0 else None)
    rotation = None
    if rotation_mode == "rigid":
        rotation = Rotation(Omega_c=Omega_c, A=float("inf"))
    elif rotation_mode == "differential":
        rotation = Rotation(Omega_c=Omega_c, A=A_over_Req * R_guess)

    result = scf_mod.hachisu_scf(rho0, r, theta, rho_c, rotation=rotation, poloidal=poloidal,
                                  toroidal=toroidal_sc, mu_e=mu_e, lmax=lmax,
                                  tol=tol, max_iter=int(max_iter))

    if not result["converged"]:
        return {"converged": False, "rho_c": rho_c, "k0": k0, "K_tor": K_tor, "Omega_c": Omega_c}

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
    }
    fields = {"rho": rho, "Phi": Phi, "u": u, "H": H, "Bphi": Bphi, "r": r, "theta": theta}
    return {"converged": True, "rho_c": rho_c, "k0": k0, "K_tor": K_tor, "Omega_c": Omega_c,
            "scalars": scalars, "fields": fields}
