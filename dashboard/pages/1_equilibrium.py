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
from terms.poloidal import Poloidal
from terms.rotation import Rotation
from terms.toroidal_sc import ToroidalSC

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
        poloidal = Poloidal(k0=k0, lmax=16)
        result = scf_mod.hachisu_scf(rho_seed, r, theta, rho_c, poloidal=poloidal, mu_e=mu_e,
                                      lmax=16, tol=1e-7, max_iter=150)
        if not result["converged"]:
            break
        rho_seed = result["rho"]
        ve = diag.virial_error_terms(result["rho"], result["Phi"], result["H"],
                                      r, theta, mu_e, poloidal=poloidal)
        if ve["VE"] > 1e-3:
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

st.sidebar.markdown("**Rotation**")
_Omega_c_reload = _reload.get("Omega_c", 0.0) if _reload else 0.0
_A_over_Req_reload = _reload.get("A_over_Req", 0.0) if _reload else 0.0
_rot_mode_default = "none" if _Omega_c_reload == 0.0 else ("rigid" if _A_over_Req_reload == 0.0 else "differential")
rotation_mode = st.sidebar.radio("rotation", ["none", "rigid", "differential"], horizontal=True,
                                  index=["none", "rigid", "differential"].index(_rot_mode_default))
Omega_c = 0.0
A_over_Req = 0.0
if rotation_mode != "none":
    Omega_c = st.sidebar.slider("Omega_c (rad/s)", 0.0, 45.0,
                                 _Omega_c_reload if _Omega_c_reload > 0 else 1.0, 0.1)
    if rotation_mode == "differential":
        A_over_Req = st.sidebar.slider(
            "A / R_eq (decay scale of the j-constant profile)", 0.05, 3.0,
            _A_over_Req_reload if _A_over_Req_reload > 0 else 0.3, 0.05)
        st.sidebar.caption(
            "j-constant law, Omega(varpi) = Omega_c A^2/(A^2+varpi^2). Small "
            "A/R_eq = strongly differential (core spins much faster than the "
            "envelope) — this is what reaches masses well above the rigid-rotation "
            "ceiling (~1.5 Msun) toward ~2.2 Msun (Yoon & Langer 2005), because the "
            "envelope stays far from equatorial breakup even as the core spins fast."
        )
    st.sidebar.caption(
        "Watch the mass-loss ratio and T/|W| scalars below — rotation has its own "
        "failure modes (equatorial mass shedding, secular/dynamical instability), "
        "separate from the VE gate."
    )

st.sidebar.markdown("**Magnetic field**")
_k0_reload = _reload["k0"] if _reload else 0.0
_K_tor_reload = _reload.get("K_tor", 0.0) if _reload else 0.0
_field_mode_default = "poloidal" if _k0_reload != 0.0 else ("toroidal (self-consistent)" if _K_tor_reload != 0.0 else "none")
field_mode = st.sidebar.radio(
    "field", ["none", "poloidal", "toroidal (self-consistent)"],
    index=["none", "poloidal", "toroidal (self-consistent)"].index(_field_mode_default),
    help="poloidal: f(u)=k0, optionally with a D6 twisted-torus imposed on top "
         "after convergence. toroidal (self-consistent): purely toroidal, "
         "B_phi=K*rho^m*varpi^(2m-1), resolved INSIDE the SCF loop (mutually "
         "exclusive with poloidal — see scf/terms/toroidal_sc.py)."
)
field_on = field_mode == "poloidal"
k0 = 0.0
zeta_target_ratio = 0.0
m_tor = 1
K_tor = 0.0
m_tor_sc = 1.0

if field_mode == "poloidal":
    st.sidebar.markdown("*Poloidal field (k0)*")
    if st.sidebar.button("find useful k0 range (empirical)"):
        with st.spinner("probing k0 (coarse grid, ~a few seconds)..."):
            k0_max = _estimate_k0_max(rho_c, mu_e, R_guess)
        k0_cache[cache_key] = k0_max
        _save_k0_cache(k0_cache)
        st.sidebar.success(f"k0_max ≈ {k0_max:.3e} (VE crosses 1e-3 here)")

    k0_max_known = k0_cache.get(cache_key, 1e-12)
    st.sidebar.caption(
        f"empirically known VE<1e-3 range (cache): up to {k0_max_known:.2e}. "
        "Not known a priori — see plan, D6. The slider goes well past this "
        "on purpose: beyond it VE will exceed 1e-3 (equilibrium unreliable) "
        "and eventually the SCF stops converging entirely (overflow in the "
        "EOS/Grad-Shafranov solve) — both are informative failure modes, "
        "not bugs, and mark the end of the equilibrium sequence (see "
        "docs/teoria.md §6). Each run is seeded from the previous converged "
        "one at the same (rho_c, mu_e, mesh) — dragging the slider up "
        "gradually reaches noticeably higher k0 than jumping straight to a "
        "large value from a cold start (tested: ~30% further at rho_c=1e12), "
        "because the crude spherical guess used for a cold start is too far "
        "from the true, already-oblate equilibrium for the Picard iteration "
        "to reach it directly."
    )
    use_slider = st.sidebar.checkbox("use log slider", value=True)
    _sign_default = "-" if _k0_reload < 0 else "+"
    sign = st.sidebar.radio("k0 sign", ["+", "-"], horizontal=True,
                             index=(0 if _sign_default == "+" else 1))
    _log_k0_default = float(np.log10(abs(_k0_reload))) if _k0_reload != 0.0 else -14.0
    if use_slider:
        # fixed, generous upper bound (not tied to k0_max_known) so the
        # slider can explore well past the known-safe range up to and
        # through where the SCF breaks down
        log_k0 = st.sidebar.slider("log10(|k0|)", -20.0, -8.0, _log_k0_default, 0.1)
        k0 = (1 if sign == "+" else -1) * 10 ** log_k0
    else:
        k0 = st.sidebar.number_input("k0 (free entry)",
                                      value=abs(_k0_reload) if _k0_reload else 1e-14, format="%.3e")
        if sign == "-":
            k0 = -abs(k0)

    st.sidebar.markdown("*Toroidal (D6, imposed on top after convergence)*")
    _zeta_reload = _reload["zeta_target_ratio"] if _reload else 0.0
    _large_bt_bp = st.sidebar.checkbox("large Bt/Bp (>2, free entry)", value=(_zeta_reload > 2.0))
    if _large_bt_bp:
        zeta_target_ratio = st.sidebar.number_input(
            "Bt/Bp target (energy)", min_value=0.0, max_value=200.0,
            value=_zeta_reload if _zeta_reload > 2.0 else 5.0, step=1.0)
        st.sidebar.caption(
            "solve_zeta_for_energy_ratio() has no upper bound on the target "
            "ratio, but the toroidal field is imposed after the SCF converges, "
            "not solved self-consistently (D6) — a strongly dominant toroidal "
            "field breaks the virial balance (VE) and will likely be blocked "
            "at export (R5, VE >= 1e-3). That's physical, not a bug."
        )
    else:
        zeta_target_ratio = st.sidebar.slider("Bt/Bp target (energy)", 0.0, 2.0,
                                               _zeta_reload if _zeta_reload <= 2.0 else 0.0, 0.05)
    m_tor = st.sidebar.slider("m_tor", 1, 4, _reload["m_tor"] if _reload else 1)

elif field_mode == "toroidal (self-consistent)":
    st.sidebar.markdown("*Self-consistent toroidal field*")
    st.sidebar.caption(
        "B_phi = K rho^m varpi^(2m-1), resolved INSIDE the SCF (implicit rho "
        "inversion, scf.py::_solve_rho_implicit) — NOT the D6 twisted torus "
        "(that needs a poloidal field first; see scf/terms/toroidal_sc.py for "
        "the distinction). Produces PROLATE deformation (R_pol > R_eq), the "
        "opposite of the poloidal/oblate case — verified in "
        "scf/tests/test_toroidal_sc.py (V-R3)."
    )
    _log_K_default = float(np.log10(_K_tor_reload)) if _K_tor_reload > 0 else -3.0
    log_K = st.sidebar.slider("log10(K)", -8.0, -2.0, _log_K_default, 0.1)
    K_tor = 10 ** log_K
    m_tor_sc = st.sidebar.slider("m (toroidal power law, B_phi ~ rho^m)", 1.0, 3.0,
                                  float(_reload.get("m_tor_sc", 1.0)) if _reload else 1.0, 0.5)

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
    "m_tor": m_tor, "Omega_c": Omega_c, "A_over_Req": A_over_Req,
    "K_tor": K_tor, "m_tor_sc": m_tor_sc,
    "Nr": Nr, "Ntheta": Ntheta, "lmax": lmax, "tol": tol,
    "max_iter": int(max_iter),
}

r = np.linspace(0, 1.3 * R_guess, Nr)
theta = np.linspace(0, np.pi, Ntheta)
rho0_cold = scf_mod.initial_guess(r, theta, rho_c, R_guess)

# A is specified in units of R_eq (per the physics prompt), but R_eq is an
# OUTPUT of the solve, not known in advance -- R_guess (the same seed
# estimate already used to size the mesh) stands in for it, matching the
# approximation already used throughout this page. Rotation does change
# the star's size somewhat (see scf/tests/test_differential_rotation.py),
# so this is a zeroth-order reference scale, not an exact R_eq.
rotation = None
if rotation_mode == "rigid":
    rotation = Rotation(Omega_c=Omega_c, A=float("inf"))
elif rotation_mode == "differential":
    rotation = Rotation(Omega_c=Omega_c, A=A_over_Req * R_guess)

# Continuation (warm start): seed from the last converged solution at the
# same (rho_c, mu_e, mesh) instead of always restarting from the crude
# spherical n=3 guess. The Picard iteration (rho -> Phi -> H -> rho_new,
# direct substitution, no sub-relaxation — see scf/scf.py) has a genuinely
# limited basin of attraction for large k0: adaptive under-relaxation was
# tested and does NOT extend it (still diverges to NaN even damped to
# omega~0.02), because the failure is the cold guess being too far from
# the true (already oblate) equilibrium, not a too-large step size.
# Continuation is the standard fix for tracing equilibrium sequences near
# this kind of limit (same technique already used internally by
# _estimate_k0_max above) — empirically it pushed the reachable k0 by
# >30% at rho_c=1e12 (0.83->0.72 in R_pol/R_eq) before hitting what looks
# like the genuine end of the sequence (see docs/teoria.md §6).
_seed_key = (rho_c, mu_e, Nr, Ntheta, rotation_mode, field_mode)
_warm = st.session_state.get("scf_warm_seed")
use_warm = (field_on or rotation_mode != "none") and _warm is not None and _warm["key"] == _seed_key
rho0 = _warm["rho"] if use_warm else rho0_cold

poloidal = Poloidal(k0=k0, lmax=lmax) if field_mode == "poloidal" and k0 != 0.0 else None
toroidal_sc = ToroidalSC(K=K_tor, m=m_tor_sc) if field_mode == "toroidal (self-consistent)" and K_tor > 0 else None

with st.spinner("running SCF..."):
    result = scf_mod.hachisu_scf(rho0, r, theta, rho_c, rotation=rotation, poloidal=poloidal,
                                  toroidal=toroidal_sc, mu_e=mu_e, lmax=lmax,
                                  tol=tol, max_iter=int(max_iter), track_virial=True)

if not result["converged"] and use_warm:
    st.caption("warm-started run did not converge — retrying from a cold (spherical) guess...")
    with st.spinner("running SCF (cold restart)..."):
        result = scf_mod.hachisu_scf(rho0_cold, r, theta, rho_c, rotation=rotation, poloidal=poloidal,
                                      toroidal=toroidal_sc, mu_e=mu_e, lmax=lmax,
                                      tol=tol, max_iter=int(max_iter), track_virial=True)

if not result["converged"]:
    st.error(f"SCF did not converge in {result['iterations']} iterations "
             f"(last delta_rho/rho_c = {result['history'][-1]:.3e}).")
    st.session_state.pop("scf_warm_seed", None)
    st.stop()

st.session_state["scf_warm_seed"] = {"key": _seed_key, "rho": result["rho"]}

rho, Phi, u, H = result["rho"], result["Phi"], result["u"], result["H"]

# toroidal (D6), if target > 0 -- only valid together with a poloidal
# field (see terms/toroidal_sc.py distinction table); mutually exclusive
# with field_mode=="toroidal (self-consistent)" by construction (that mode
# never sets `poloidal`, so u==0 there and this block is skipped by the
# k0!=0.0 guard)
Bphi_d6 = np.zeros_like(rho)
u_c = None
zeta_used = 0.0
if field_mode == "poloidal" and zeta_target_ratio > 0 and k0 != 0.0:
    try:
        Bphi_d6, zeta_used, u_c = tor.solve_zeta_for_energy_ratio(
            u, H, r, theta, zeta_target_ratio, m_tor=m_tor)
    except ValueError as e:
        st.warning(f"Toroidal not imposed: {e}")

# self-consistent terms' contributions (T from rotation, Br/Bth from
# poloidal, Bphi from the self-consistent toroidal branch — zero for
# whichever terms are None)
ve_terms = diag.virial_error_terms(rho, Phi, H, r, theta, mu_e,
                                    rotation=rotation, poloidal=poloidal, toroidal=toroidal_sc)
T = ve_terms["T"]
Br, Bth = ve_terms["Br"], ve_terms["Btheta"]
# D6-imposed and self-consistent toroidal are mutually exclusive in
# practice (poloidal vs toroidal-self-consistent field_mode), so summing
# is safe -- exactly one of the two is ever nonzero
Bphi = Bphi_d6 + ve_terms["Bphi"]

# recompute VE/E_mag with the COMBINED Bphi (ve_terms doesn't know about
# the D6 imposition, which isn't a "term") — single source of truth for
# the residual formula (diag.virial_error), not reimplemented here
VE, W, Pi, E_mag, _T2 = diag.virial_error(rho, Phi, H, Br, Bth, Bphi, r, theta, mu_e, T=T)
E_pol, E_tor, _ = diag.magnetic_energies(Br, Bth, Bphi, r, theta)
M = scf_mod.total_mass(rho, r, theta)
R_eq, R_pol = diag.equatorial_polar_radii(H, r, theta)
rho_mean = M / (4.0 / 3.0 * np.pi * ((R_eq**2 * R_pol) ** (1.0 / 3.0)) ** 3) if R_eq > 0 else float("nan")
T_over_W = T / abs(W) if W != 0 else float("nan")
mass_loss_ratio = diag.equatorial_mass_loss_ratio(Phi, rotation, r, theta, R_eq)
rotation_period_s = (2 * np.pi / Omega_c) if Omega_c > 0 else float("inf")

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

if rotation is not None:
    # T/|W| stability gate (Ostriker & Bodenheimer 1973): >= 0.14 is the
    # secular (non-axisymmetric, viscosity/GR-driven) instability
    # threshold; >= 0.27 is the dynamical (bar-mode) threshold. Same
    # traffic-light convention as VE, and blocks export (Tab 3) like VE
    # does (R5) once >= 0.14.
    if T_over_W < 0.14:
        st.success(f"T/|W| = {T_over_W:.3e}  (< 0.14, secular stability limit ✓)")
    elif T_over_W < 0.27:
        st.warning(f"T/|W| = {T_over_W:.3e}  (>= 0.14 — secular non-axisymmetric "
                    "instability threshold, Ostriker & Bodenheimer 1973 ⚠)")
    else:
        st.error(f"T/|W| = {T_over_W:.3e}  (>= 0.27 — dynamical bar-mode instability "
                  "threshold ✗)")
    if mass_loss_ratio >= 1.0:
        st.error(f"equatorial mass-loss ratio = {mass_loss_ratio:.3f} >= 1 — this "
                  "configuration is past Keplerian breakup and should not be trusted "
                  "(effective gravity has already vanished at the equator).")
    elif mass_loss_ratio >= 0.9:
        st.warning(f"equatorial mass-loss ratio = {mass_loss_ratio:.3f} — close to "
                   "Keplerian breakup (ratio -> 1).")

st.subheader("Scalars")
Bpol_grid = np.sqrt(Br**2 + Bth**2)
B_pol_max_gauss = np.max(Bpol_grid)
B_tor_max_gauss = np.max(np.abs(Bphi))
frac_torus = tor.closed_torus_volume_fraction(u, rho, r, theta, u_c) if u_c is not None else 0.0

# Surface field (D6/observational quantity): B AT the stellar surface
# (rho=0 boundary, via the same per-theta interpolation toroidal.find_uc()
# uses — diagnostics.surface_radius), not the interior max. This is what
# would actually be compared to a real magnetized white dwarf's measured
# (polarimetric) field. Generally << B_pol,max, since the poloidal field
# peaks well inside the star for this k0=const source, not at the surface.
# B_pole/B_eq ("dipolarity") is exactly 2 for a pure dipole; deviation
# flags multipole content. B_pole sits on the symmetry axis (same
# 1/sin(theta) singularity fixed in poloidal_field()) — cross-checked
# against an independent quadratic extrapolation from off-axis points; see
# diag.surface_dipolarity() and scf/tests/test_diagnostics.py (verified
# convergent across ntheta=65/129, agreement ~0.03%/0.01%).
_dip = diag.surface_dipolarity(Bpol_grid, H, r, theta)
B_central_gauss = float(np.sqrt(Br[1, 0] ** 2 + Bth[1, 0] ** 2))
_extrap_rel_diff = (abs(_dip["B_pole"] - _dip["B_pole_extrapolated"]) / _dip["B_pole"]
                     if _dip["B_pole"] != 0 else 0.0)
if poloidal is not None:
    # these checks are specific to the poloidal field (B_r on-axis
    # regularization, interior-vs-surface field) -- meaningless (and
    # trivially 0>=0) when there is no poloidal field at all
    if _dip["B_surf_max"] >= B_central_gauss:
        st.warning(
            f"Sanity check failed: B_surf_max ({_dip['B_surf_max']:.3e} G) >= "
            f"B_central ({B_central_gauss:.3e} G) — surface field should not "
            "exceed the interior field for this k0=const source. Investigate "
            "before trusting these numbers."
        )
    if _extrap_rel_diff > 0.05:
        st.warning(
            f"B_pole direct ({_dip['B_pole']:.3e} G) vs quadratic extrapolation "
            f"({_dip['B_pole_extrapolated']:.3e} G) disagree by {_extrap_rel_diff:.1%} "
            "— the on-axis field estimate may not be resolved at this mesh; "
            "try a higher Ntheta."
        )

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
    # B_r at r>0 on-axis is a real, generally nonzero on-axis field value
    # (see diagnostics.poloidal_field docstring — L'Hopital limit, not the
    # spurious 0 the old coordinate-singularity mask used to return); r_idx=1
    # (not r=0) sidesteps the separate 1/r^2 singularity at the true origin.
    "B_central (G)": B_central_gauss,
    "B_polo (G)": _dip["B_pole"],
    "B_eq (G)": _dip["B_eq"],
    "B_surf,max (G)": _dip["B_surf_max"],
    "B_polo/B_eq (dipolarity)": _dip["dipolarity"],
    "B_tor,max (G)": B_tor_max_gauss,
    "torus volume fraction": frac_torus,
    "VE": VE,
    "T (erg)": T,
    "T/|W|": T_over_W,
    "Omega_c (rad/s)": Omega_c,
    "rotation period (s)": rotation_period_s,
    "equatorial mass-loss ratio": mass_loss_ratio,
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
    st.pyplot(plots.plot_flux_contours(u, r, theta, rho=rho, u_c=u_c, H=H))
    if u_c is not None:
        # numerical tangency check (D6) — NOT a visual read of the plot: the
        # u=u_c and H=0 contours can appear offset by up to ~1 grid cell
        # purely from matplotlib's 2D contour interpolation on the curved
        # (r,theta)->(varpi,z) grid, even when mathematically exactly
        # tangent. See toroidal.check_uc_tangency() / scf/tests/test_toroidal.py.
        _tc = tor.check_uc_tangency(u, rho, H, r, theta, u_c)
        _deg = np.degrees(_tc["theta_tangent"])
        if _tc["unique"] and not _tc["vacuum_leak"]:
            st.caption(
                f"torus boundary verified tangent to the surface at θ≈{_deg:.1f}° "
                f"(margin {_tc['margin']:.1e}, exact by construction — the shaded "
                "region never extends past the physical surface)."
            )
        else:
            st.warning(
                f"torus boundary tangency check failed: unique={_tc['unique']}, "
                f"vacuum_leak={_tc['vacuum_leak']} — see toroidal.check_uc_tangency()."
            )
with f3:
    if np.any(Bphi != 0):
        st.pyplot(plots.plot_toroidal(Bphi, r, theta))
    else:
        st.info("no toroidal field (Bt/Bp target = 0)")

st.subheader("Density profile")
st.pyplot(plots.plot_density_profile(rho, r, theta))

# ---------------- persistence (R2) ----------------
st.divider()
if st.button("save this run"):
    scalars_json = {k: float(v) for k, v in scalars_display.items()}
    fields = {"rho": rho, "Phi": Phi, "u": u, "H": H, "Bphi": Bphi, "r": r, "theta": theta}
    h = store.save_run(params, scalars_json, fields)
    st.success(f"run saved: {h}")
