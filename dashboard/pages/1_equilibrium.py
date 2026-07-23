"""Tab 1 — Equilibrium: single SCF run, inspection. R1: physics only via scf.*"""

import json
import sys
from pathlib import Path

import numpy as np
import streamlit as st

_DASHBOARD_DIR = Path(__file__).resolve().parent.parent
_REPO_ROOT = _DASHBOARD_DIR.parent
sys.path.insert(0, str(_REPO_ROOT / "scf"))
sys.path.insert(0, str(_DASHBOARD_DIR))

import scf as scf_mod
import diagnostics as diag
import toroidal as tor
import units
import store
import plots
import seed

st.set_page_config(page_title="Equilibrium — wd-magnetizada", layout="wide")
st.title("Tab 1 — Equilibrium")

K0_RANGE_CACHE = _DASHBOARD_DIR / "k0_range_cache.json"


def _load_k0_cache():
    if K0_RANGE_CACHE.exists():
        return json.loads(K0_RANGE_CACHE.read_text())
    return {}


def _save_k0_cache(cache):
    K0_RANGE_CACHE.write_text(json.dumps(cache, indent=2))


def _estimate_k0_max(rho_c, mu_e, R_guess):
    """Raises k0 geometrically (coarse grid, continuation) until VE>1e-3
    or the SCF stops converging. Only calls scf.* (R1)."""
    nr, ntheta = 81, 33
    r = np.linspace(0, 1.3 * R_guess, nr)
    theta = np.linspace(0, np.pi, ntheta)
    rho_seed = scf_mod.initial_guess(r, theta, rho_c, R_guess)
    k0 = 1e-20
    k0_ok = 0.0
    for _ in range(60):
        result = scf_mod.hachisu_scf(rho_seed, r, theta, rho_c, k0=k0, mu_e=mu_e,
                                      lmax=16, tol=1e-7, max_iter=150)
        if not result["converged"]:
            break
        rho_seed = result["rho"]
        Br, Bth = diag.poloidal_field(result["u"], r, theta)
        VE, _, _, _ = diag.virial_error(result["rho"], result["Phi"], result["H"],
                                         Br, Bth, np.zeros_like(Br), r, theta, mu_e)
        if VE > 1e-3:
            break
        k0_ok = k0
        k0 *= 2.5
    return max(k0_ok, 1e-20)


# ---------------- reload a run (Tab 4 -> "reload in Tab 1") ----------------
_reload = st.session_state.pop("reload_run_params", None)
if _reload:
    st.info(f"Parameters loaded from a saved run (rho_c={_reload['rho_c']:.3e}, "
            f"k0={_reload['k0']:.3e}).")

# ---------------- sidebar ----------------
st.sidebar.header("Physical parameters")
_rho_c_options = [10 ** e for e in np.arange(6, 12.01, 0.1)]
# default high enough to reproduce Chandrasekhar within <1% at k0=0 (V1) —
# see tests/test_scf_v1.py and dashboard/tests/test_smoke.py
_rho_c_target = _reload["rho_c"] if _reload else 1e12
_rho_c_default = min(_rho_c_options, key=lambda x: abs(x - _rho_c_target))
rho_c = st.sidebar.select_slider(
    "rho_c (g/cm³)", options=_rho_c_options, value=_rho_c_default,
    format_func=lambda x: f"{x:.2e}",
)
mu_e = st.sidebar.number_input("mu_e", min_value=1.0, max_value=2.5,
                                value=_reload["mu_e"] if _reload else 2.0, step=0.1)

R_guess = seed.r_guess(rho_c)
cache_key = f"{rho_c:.3e}_{mu_e:.2f}"
k0_cache = _load_k0_cache()

st.sidebar.markdown("**Poloidal field (k0)**")
_k0_reload = _reload["k0"] if _reload else 0.0
field_on = st.sidebar.checkbox("poloidal field on", value=(_k0_reload != 0.0))
k0 = 0.0
if field_on:
    if st.sidebar.button("find useful k0 range (empirical)"):
        with st.spinner("probing k0 (coarse grid, ~a few seconds)..."):
            k0_max = _estimate_k0_max(rho_c, mu_e, R_guess)
        k0_cache[cache_key] = k0_max
        _save_k0_cache(k0_cache)
        st.sidebar.success(f"k0_max ≈ {k0_max:.3e} (VE crosses 1e-3 here)")

    k0_max_known = k0_cache.get(cache_key, 1e-12)
    st.sidebar.caption(f"known range (cache): up to {k0_max_known:.2e}. "
                        "Not known a priori — see plan, D6.")
    use_slider = st.sidebar.checkbox("use log slider", value=True)
    _sign_default = "-" if _k0_reload < 0 else "+"
    sign = st.sidebar.radio("k0 sign", ["+", "-"], horizontal=True,
                             index=(0 if _sign_default == "+" else 1))
    _log_k0_default = float(np.log10(abs(_k0_reload))) if _k0_reload != 0.0 else -14.0
    if use_slider:
        log_k0 = st.sidebar.slider("log10(|k0|)", -20.0, float(np.log10(k0_max_known * 3)),
                                    _log_k0_default, 0.1)
        k0 = (1 if sign == "+" else -1) * 10 ** log_k0
    else:
        k0 = st.sidebar.number_input("k0 (free entry)",
                                      value=abs(_k0_reload) if _k0_reload else 1e-14, format="%.3e")
        if sign == "-":
            k0 = -abs(k0)

st.sidebar.markdown("**Toroidal (D6)**")
zeta_target_ratio = st.sidebar.slider("Bt/Bp target (energy)", 0.0, 2.0,
                                       _reload["zeta_target_ratio"] if _reload else 0.0, 0.05)
m_tor = st.sidebar.slider("m_tor", 1, 4, _reload["m_tor"] if _reload else 1)

st.sidebar.header("Numerical parameters")
_Nr_opts = [65, 129, 161, 257]
_Ntheta_opts = [33, 65, 129]
Nr = st.sidebar.select_slider(
    "Nr", options=_Nr_opts,
    value=min(_Nr_opts, key=lambda x: abs(x - _reload["Nr"])) if _reload else 129)
Ntheta = st.sidebar.select_slider(
    "Ntheta", options=_Ntheta_opts,
    value=min(_Ntheta_opts, key=lambda x: abs(x - _reload["Ntheta"])) if _reload else 129)
lmax = st.sidebar.slider("l_max", 4, 32, _reload["lmax"] if _reload else 16)
_tol_opts = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
tol = st.sidebar.select_slider(
    "tol", options=_tol_opts,
    value=min(_tol_opts, key=lambda x: abs(x - _reload["tol"])) if _reload else 1e-6)
max_iter = st.sidebar.number_input("max_iter", value=_reload["max_iter"] if _reload else 200, step=50)

st.sidebar.caption(
    "theta covers the full [0, pi] range (no equatorial symmetry) — "
    "deliberate, protects against spurious asymmetric modes (m=1)."
)
st.sidebar.caption(
    "Note: the (rho_c, k0) formulation used here does not need "
    "sub-relaxation (omega) — the original plan's two-surface-condition "
    "method needed it, but is unstable for this EOS (see scf/scf.py). "
    "There is no omega control because it doesn't correspond to anything "
    "in the real solver."
)

# ---------------- run the SCF ----------------
params = {
    "rho_c": rho_c, "mu_e": mu_e, "k0": k0, "zeta_target_ratio": zeta_target_ratio,
    "m_tor": m_tor, "Nr": Nr, "Ntheta": Ntheta, "lmax": lmax, "tol": tol,
    "max_iter": int(max_iter),
}

r = np.linspace(0, 1.3 * R_guess, Nr)
theta = np.linspace(0, np.pi, Ntheta)
rho0 = scf_mod.initial_guess(r, theta, rho_c, R_guess)

with st.spinner("running SCF..."):
    result = scf_mod.hachisu_scf(rho0, r, theta, rho_c, k0=k0, mu_e=mu_e, lmax=lmax,
                                  tol=tol, max_iter=int(max_iter), track_virial=True)

if not result["converged"]:
    st.error(f"SCF did not converge in {result['iterations']} iterations "
             f"(last delta_rho/rho_c = {result['history'][-1]:.3e}).")
    st.stop()

rho, Phi, u, H = result["rho"], result["Phi"], result["u"], result["H"]

# toroidal (D6), if target > 0
Bphi = np.zeros_like(rho)
u_c = None
zeta_used = 0.0
if zeta_target_ratio > 0 and k0 != 0.0:
    try:
        Bphi, zeta_used, u_c = tor.solve_zeta_for_energy_ratio(
            u, rho, r, theta, zeta_target_ratio, m_tor=m_tor)
    except ValueError as e:
        st.warning(f"Toroidal not imposed: {e}")

Br, Bth = diag.poloidal_field(u, r, theta)
VE, W, Pi, E_mag = diag.virial_error(rho, Phi, H, Br, Bth, Bphi, r, theta, mu_e)
E_pol, E_tor, _ = diag.magnetic_energies(Br, Bth, Bphi, r, theta)
M = scf_mod.total_mass(rho, r, theta)
R_eq, R_pol = diag.equatorial_polar_radii(rho, r, theta)
rho_mean = M / (4.0 / 3.0 * np.pi * ((R_eq**2 * R_pol) ** (1.0 / 3.0)) ** 3) if R_eq > 0 else float("nan")

# ---------------- main panel ----------------
col1, col2 = st.columns(2)
with col1:
    st.pyplot(plots.plot_convergence(result["history"]))
with col2:
    if result["ve_history"]:
        st.pyplot(plots.plot_virial_history(result["ve_history"]))
    else:
        st.info("VE history not available")

if VE < 1e-3:
    st.success(f"VE = {VE:.3e}  (< 1e-3, plan's V3 criterion ✓)")
else:
    st.error(f"VE = {VE:.3e}  (>= 1e-3, plan's V3 criterion ✗ — equilibrium not reliable)")

st.subheader("Scalars")
B_pol_max_gauss = np.max(np.sqrt(Br**2 + Bth**2))
B_tor_max_gauss = np.max(np.abs(Bphi))
frac_torus = tor.closed_torus_volume_fraction(u, rho, r, theta, u_c) if u_c is not None else 0.0

scalars_display = {
    "M/M_sun": M / units.M_SUN,
    "R_eq (km)": units.cm_to_km(R_eq),
    "R_pol (km)": units.cm_to_km(R_pol),
    "R_pol/R_eq": R_pol / R_eq if R_eq > 0 else float("nan"),
    "rho_c confirmed (g/cm³)": rho[0, 0],
    "mean rho (g/cm³)": rho_mean,
    "W (erg)": W,
    "E_int = ∫P dV (erg)": Pi,
    "E_mag (erg)": E_mag,
    "E_pol (erg)": E_pol,
    "E_tor (erg)": E_tor,
    "E_mag/|W|": E_mag / abs(W) if W != 0 else float("nan"),
    "B_pol,max (G)": B_pol_max_gauss,
    # Br,Bth at r=0 are zeroed by construction (poloidal_field, coordinate
    # singularity); use the first grid point with r>0 as a proxy for the center
    "B_central (G)": np.sqrt(Br[1, 0] ** 2 + Bth[1, 0] ** 2),
    "B_tor,max (G)": B_tor_max_gauss,
    "torus volume fraction": frac_torus,
    "VE": VE,
}


def _format_scalar(key, value):
    """Formatting rule (R4): field in gauss -> scientific notation; radii
    in km -> fixed decimal places; everything else -> generic scientific
    notation. Values in scalars_display are already in display units
    (km, G) — this only decides the STRING, via units.py (single source
    of truth)."""
    if not isinstance(value, (int, float, np.floating)):
        return value
    if key.endswith("(G)"):
        return units.format_gauss(value)
    if key.endswith("(km)"):
        return units.format_km_value(value)
    return f"{value:.4e}"


st.table({"quantity": list(scalars_display.keys()),
          "value": [_format_scalar(k, v) for k, v in scalars_display.items()]})

st.subheader("Bt/Bp — two definitions")
ratio_energy, ratio_amp = tor.bt_bp_ratios(Br, Bth, Bphi, r, theta)
c1, c2 = st.columns(2)
c1.metric("Bt/Bp (energy) = E_tor/E_pol", f"{ratio_energy:.4f}")
c2.metric("Bt/Bp (amplitude) = max|Bphi|/max|Bpol|", f"{ratio_amp:.4f}")
st.caption(
    "The two differ by orders of magnitude because the toroidal field is "
    "confined to a small volume (D6) — the literature is often careless "
    "about which one it uses."
)

st.subheader("Figures (meridional plane)")
f1, f2, f3 = st.columns(3)
with f1:
    st.pyplot(plots.plot_density(rho, r, theta, H=H))
with f2:
    st.pyplot(plots.plot_flux_contours(u, r, theta, u_c=u_c))
with f3:
    if np.any(Bphi != 0):
        st.pyplot(plots.plot_toroidal(Bphi, r, theta))
    else:
        st.info("no toroidal field (Bt/Bp target = 0)")

# ---------------- persistence (R2) ----------------
st.divider()
if st.button("save this run"):
    scalars_json = {k: float(v) for k, v in scalars_display.items()}
    fields = {"rho": rho, "Phi": Phi, "u": u, "H": H, "Bphi": Bphi, "r": r, "theta": theta}
    h = store.save_run(params, scalars_json, fields)
    st.success(f"run saved: {h}")
