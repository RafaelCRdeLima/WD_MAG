"""Funcao de trabalho para a varredura (Aba 2), roda em processos separados
via ProcessPoolExecutor. Precisa estar num modulo importavel (nao dentro do
script da pagina Streamlit) para ser picklable pelo multiprocessing.

So' chama scf.* — toda fisica mora la (R1)."""

import sys
from pathlib import Path

_SCF_DIR = Path(__file__).resolve().parent.parent / "scf"
if str(_SCF_DIR) not in sys.path:
    sys.path.insert(0, str(_SCF_DIR))


def run_one(params: dict) -> dict:
    """Roda um SCF para (rho_c, k0) e retorna escalares + campos (para salvar
    com store.save_run). Nao converge -> retorna converged=False, sem campos."""
    import numpy as np
    import scf as scf_mod
    import diagnostics as diag

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
    result = scf_mod.hachisu_scf(rho0, r, theta, rho_c, k0=k0, mu_e=mu_e,
                                  lmax=lmax, tol=tol, max_iter=int(max_iter))

    if not result["converged"]:
        return {"converged": False, "rho_c": rho_c, "k0": k0}

    rho, Phi, u, H = result["rho"], result["Phi"], result["u"], result["H"]
    Br, Bth = diag.poloidal_field(u, r, theta)
    Bphi = np.zeros_like(rho)
    VE, W, Pi, E_mag = diag.virial_error(rho, Phi, H, Br, Bth, Bphi, r, theta, mu_e)
    M = scf_mod.total_mass(rho, r, theta)
    R_eq, R_pol = diag.equatorial_polar_radii(rho, r, theta)

    scalars = {
        "M/M☉": M / 1.989e33,
        "R_eq (km)": R_eq / 1.0e5,
        "R_pol (km)": R_pol / 1.0e5,
        "R_pol/R_eq": R_pol / R_eq if R_eq > 0 else float("nan"),
        "W (erg)": W,
        "E_mag (erg)": E_mag,
        "E_mag/|W|": E_mag / abs(W) if W != 0 else float("nan"),
        "VE": VE,
    }
    fields = {"rho": rho, "Phi": Phi, "u": u, "H": H, "Bphi": Bphi, "r": r, "theta": theta}
    return {"converged": True, "rho_c": rho_c, "k0": k0, "scalars": scalars, "fields": fields}
