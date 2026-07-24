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
    """Runs one SCF for (rho_c, k0) and returns scalars + fields (to save
    with store.save_run). Does not converge -> returns converged=False, no fields."""
    import numpy as np
    import scf as scf_mod
    import diagnostics as diag
    import units
    from terms.poloidal import Poloidal

    rho_c = params["rho_c"]
    k0 = params["k0"]
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
    poloidal = Poloidal(k0=k0, lmax=lmax) if k0 != 0.0 else None
    result = scf_mod.hachisu_scf(rho0, r, theta, rho_c, poloidal=poloidal, mu_e=mu_e,
                                  lmax=lmax, tol=tol, max_iter=int(max_iter))

    if not result["converged"]:
        return {"converged": False, "rho_c": rho_c, "k0": k0}

    rho, Phi, u, H = result["rho"], result["Phi"], result["u"], result["H"]
    ve = diag.virial_error_terms(rho, Phi, H, r, theta, mu_e, poloidal=poloidal)
    Br, Bth, VE, W, E_mag = ve["Br"], ve["Btheta"], ve["VE"], ve["W"], ve["E_mag"]
    M = scf_mod.total_mass(rho, r, theta)
    R_eq, R_pol = diag.equatorial_polar_radii(H, r, theta)

    Bpol_grid = np.sqrt(Br**2 + Bth**2)
    B_pol_max = np.max(Bpol_grid)
    # dipolarity (B_pole/B_eq at the stellar surface, via the same per-theta
    # surface_radius interpolation find_uc() uses) — exactly 2 for a pure
    # dipole, deviation flags multipole content. See diagnostics.surface_dipolarity.
    dip = diag.surface_dipolarity(Bpol_grid, H, r, theta)
    scalars = {
        "M/M_sun": units.g_to_msun(M),
        "R_eq (km)": units.cm_to_km(R_eq),
        "R_pol (km)": units.cm_to_km(R_pol),
        "R_pol/R_eq": R_pol / R_eq if R_eq > 0 else float("nan"),
        "W (erg)": W,
        "E_mag (erg)": E_mag,
        "E_mag/|W|": E_mag / abs(W) if W != 0 else float("nan"),
        "B_pol,max (G)": B_pol_max,
        "B_polo (G)": dip["B_pole"],
        "B_eq (G)": dip["B_eq"],
        "dipolarity": dip["dipolarity"],
        "VE": VE,
    }
    Bphi = np.zeros_like(rho)  # self-consistent toroidal not exposed in the sweep grid (yet); D6 imposition happens in Tab 3
    fields = {"rho": rho, "Phi": Phi, "u": u, "H": H, "Bphi": Bphi, "r": r, "theta": theta}
    return {"converged": True, "rho_c": rho_c, "k0": k0, "scalars": scalars, "fields": fields}
