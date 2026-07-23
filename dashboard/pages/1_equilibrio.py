"""Aba 1 — Equilibrio: execucao unica do SCF, inspecao. R1: fisica so' via scf.*"""

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

st.set_page_config(page_title="Equilíbrio — wd-magnetizada", layout="wide")
st.title("Aba 1 — Equilíbrio")

K0_RANGE_CACHE = _DASHBOARD_DIR / "k0_range_cache.json"


def _load_k0_cache():
    if K0_RANGE_CACHE.exists():
        return json.loads(K0_RANGE_CACHE.read_text())
    return {}


def _save_k0_cache(cache):
    K0_RANGE_CACHE.write_text(json.dumps(cache, indent=2))


def _estimate_k0_max(rho_c, mu_e, R_guess):
    """Sobe k0 geometricamente (malha grosseira, continuacao) ate' VE>1e-3
    ou a SCF parar de convergir. So' chama scf.* (R1)."""
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


# ---------------- reload de uma corrida (Aba 4 -> "recarregar na Aba 1") ----------------
_reload = st.session_state.pop("reload_run_params", None)
if _reload:
    st.info(f"Parâmetros carregados de uma corrida salva (ρc={_reload['rho_c']:.3e}, "
            f"k0={_reload['k0']:.3e}).")

# ---------------- barra lateral ----------------
st.sidebar.header("Parâmetros físicos")
_rho_c_options = [10 ** e for e in np.arange(6, 12.01, 0.1)]
# default alto o suficiente p/ reproduzir Chandrasekhar a <1% com k0=0 (V1) —
# ver tests/test_scf_v1.py e dashboard/tests/test_smoke.py
_rho_c_target = _reload["rho_c"] if _reload else 1e12
_rho_c_default = min(_rho_c_options, key=lambda x: abs(x - _rho_c_target))
rho_c = st.sidebar.select_slider(
    "ρc (g/cm³)", options=_rho_c_options, value=_rho_c_default,
    format_func=lambda x: f"{x:.2e}",
)
mu_e = st.sidebar.number_input("μₑ", min_value=1.0, max_value=2.5,
                                value=_reload["mu_e"] if _reload else 2.0, step=0.1)

R_guess = seed.r_guess(rho_c)
cache_key = f"{rho_c:.3e}_{mu_e:.2f}"
k0_cache = _load_k0_cache()

st.sidebar.markdown("**Campo poloidal (k0)**")
_k0_reload = _reload["k0"] if _reload else 0.0
campo_ligado = st.sidebar.checkbox("campo poloidal ligado", value=(_k0_reload != 0.0))
k0 = 0.0
if campo_ligado:
    if st.sidebar.button("descobrir faixa útil de k0 (empírico)"):
        with st.spinner("sondando k0 (malha grosseira, ~poucos segundos)..."):
            k0_max = _estimate_k0_max(rho_c, mu_e, R_guess)
        k0_cache[cache_key] = k0_max
        _save_k0_cache(k0_cache)
        st.sidebar.success(f"k0_max ≈ {k0_max:.3e} (VE cruza 1e-3 aqui)")

    k0_max_known = k0_cache.get(cache_key, 1e-12)
    st.sidebar.caption(f"faixa conhecida (cache): até {k0_max_known:.2e}. "
                        "Não conhecida a priori — ver plano, D6.")
    use_slider = st.sidebar.checkbox("usar slider log", value=True)
    _sign_default = "-" if _k0_reload < 0 else "+"
    sign = st.sidebar.radio("sinal de k0", ["+", "-"], horizontal=True,
                             index=(0 if _sign_default == "+" else 1))
    _log_k0_default = float(np.log10(abs(_k0_reload))) if _k0_reload != 0.0 else -14.0
    if use_slider:
        log_k0 = st.sidebar.slider("log10(|k0|)", -20.0, float(np.log10(k0_max_known * 3)),
                                    _log_k0_default, 0.1)
        k0 = (1 if sign == "+" else -1) * 10 ** log_k0
    else:
        k0 = st.sidebar.number_input("k0 (entrada livre)",
                                      value=abs(_k0_reload) if _k0_reload else 1e-14, format="%.3e")
        if sign == "-":
            k0 = -abs(k0)

st.sidebar.markdown("**Toroidal (D6)**")
zeta_target_ratio = st.sidebar.slider("Bt/Bp alvo (energia)", 0.0, 2.0,
                                       _reload["zeta_target_ratio"] if _reload else 0.0, 0.05)
m_tor = st.sidebar.slider("m_tor", 1, 4, _reload["m_tor"] if _reload else 1)

st.sidebar.header("Parâmetros numéricos")
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
    "θ cobre [0, π] inteiro (sem simetria equatorial) — deliberado, protege "
    "contra soluções assimétricas espúrias (m=1)."
)
st.sidebar.caption(
    "Nota: a formulação (ρc, k0) usada aqui não precisa de sub-relaxação "
    "(ω) — o método de duas condições de superfície do plano original "
    "precisava, mas é instável para esta EOS (ver scf/scf.py). Não há "
    "controle de ω porque não corresponde a nada no solver real."
)

# ---------------- executa o SCF ----------------
params = {
    "rho_c": rho_c, "mu_e": mu_e, "k0": k0, "zeta_target_ratio": zeta_target_ratio,
    "m_tor": m_tor, "Nr": Nr, "Ntheta": Ntheta, "lmax": lmax, "tol": tol,
    "max_iter": int(max_iter),
}

r = np.linspace(0, 1.3 * R_guess, Nr)
theta = np.linspace(0, np.pi, Ntheta)
rho0 = scf_mod.initial_guess(r, theta, rho_c, R_guess)

with st.spinner("rodando SCF..."):
    result = scf_mod.hachisu_scf(rho0, r, theta, rho_c, k0=k0, mu_e=mu_e, lmax=lmax,
                                  tol=tol, max_iter=int(max_iter), track_virial=True)

if not result["converged"]:
    st.error(f"SCF não convergiu em {result['iterations']} iterações "
             f"(último Δρ/ρc = {result['history'][-1]:.3e}).")
    st.stop()

rho, Phi, u, H = result["rho"], result["Phi"], result["u"], result["H"]

# toroidal (D6), se alvo > 0
Bphi = np.zeros_like(rho)
u_c = None
zeta_used = 0.0
if zeta_target_ratio > 0 and k0 != 0.0:
    try:
        Bphi, zeta_used, u_c = tor.solve_zeta_for_energy_ratio(
            u, rho, r, theta, zeta_target_ratio, m_tor=m_tor)
    except ValueError as e:
        st.warning(f"Toroidal não imposto: {e}")

Br, Bth = diag.poloidal_field(u, r, theta)
VE, W, Pi, E_mag = diag.virial_error(rho, Phi, H, Br, Bth, Bphi, r, theta, mu_e)
E_pol, E_tor, _ = diag.magnetic_energies(Br, Bth, Bphi, r, theta)
M = scf_mod.total_mass(rho, r, theta)
R_eq, R_pol = diag.equatorial_polar_radii(rho, r, theta)
rho_mean = M / (4.0 / 3.0 * np.pi * ((R_eq**2 * R_pol) ** (1.0 / 3.0)) ** 3) if R_eq > 0 else float("nan")

# ---------------- painel principal ----------------
col1, col2 = st.columns(2)
with col1:
    st.pyplot(plots.plot_convergence(result["history"]))
with col2:
    if result["ve_history"]:
        st.pyplot(plots.plot_virial_history(result["ve_history"]))
    else:
        st.info("histórico de VE não disponível")

if VE < 1e-3:
    st.success(f"VE = {VE:.3e}  (< 1e-3, critério V3 do plano ✓)")
else:
    st.error(f"VE = {VE:.3e}  (≥ 1e-3, critério V3 do plano ✗ — equilíbrio não confiável)")

st.subheader("Escalares")
B_pol_max_gauss = np.max(np.sqrt(Br**2 + Bth**2))
B_tor_max_gauss = np.max(np.abs(Bphi))
frac_torus = tor.closed_torus_volume_fraction(u, rho, r, theta, u_c) if u_c is not None else 0.0

scalars_display = {
    "M/M☉": M / units.M_SUN,
    "R_eq (km)": units.cm_to_km(R_eq),
    "R_pol (km)": units.cm_to_km(R_pol),
    "R_pol/R_eq": R_pol / R_eq if R_eq > 0 else float("nan"),
    "ρc confirmado (g/cm³)": rho[0, 0],
    "ρ média (g/cm³)": rho_mean,
    "W (erg)": W,
    "E_int = ∫P dV (erg)": Pi,
    "E_mag (erg)": E_mag,
    "E_pol (erg)": E_pol,
    "E_tor (erg)": E_tor,
    "E_mag/|W|": E_mag / abs(W) if W != 0 else float("nan"),
    "B_pol,max (G)": B_pol_max_gauss,
    # Br,Bth em r=0 sao zerados por construcao (poloidal_field, singularidade
    # de coordenada); usa o primeiro ponto de grade com r>0 como proxy do centro
    "B_central (G)": np.sqrt(Br[1, 0] ** 2 + Bth[1, 0] ** 2),
    "B_tor,max (G)": B_tor_max_gauss,
    "fração de volume do toro": frac_torus,
    "VE": VE,
}


def _format_scalar(key, value):
    """Regra de formatacao (R4): campo em gauss -> notacao cientifica; raios
    em km -> casas decimais fixas; resto -> notacao cientifica generica.
    Valores em scalars_display ja estao nas unidades de exibicao (km, G) —
    aqui so' decide a STRING, via units.py (ponto unico de verdade)."""
    if not isinstance(value, (int, float, np.floating)):
        return value
    if key.endswith("(G)"):
        return units.format_gauss(value)
    if key.endswith("(km)"):
        return units.format_km_value(value)
    return f"{value:.4e}"


st.table({"quantidade": list(scalars_display.keys()),
          "valor": [_format_scalar(k, v) for k, v in scalars_display.items()]})

st.subheader("Bt/Bp — duas definições")
ratio_energy, ratio_amp = tor.bt_bp_ratios(Br, Bth, Bphi, r, theta)
c1, c2 = st.columns(2)
c1.metric("Bt/Bp (energia) = E_tor/E_pol", f"{ratio_energy:.4f}")
c2.metric("Bt/Bp (amplitude) = max|Bφ|/max|Bpol|", f"{ratio_amp:.4f}")
st.caption(
    "As duas diferem por ordens de grandeza porque o toroidal fica confinado "
    "a um volume pequeno (D6) — a literatura é descuidada sobre qual usa."
)

st.subheader("Figuras (plano meridional)")
f1, f2, f3 = st.columns(3)
with f1:
    st.pyplot(plots.plot_density(rho, r, theta, H=H))
with f2:
    st.pyplot(plots.plot_flux_contours(u, r, theta, u_c=u_c))
with f3:
    if np.any(Bphi != 0):
        st.pyplot(plots.plot_toroidal(Bphi, r, theta))
    else:
        st.info("sem campo toroidal (Bt/Bp alvo = 0)")

# ---------------- persistencia (R2) ----------------
st.divider()
if st.button("salvar esta corrida"):
    scalars_json = {k: float(v) for k, v in scalars_display.items()}
    fields = {"rho": rho, "Phi": Phi, "u": u, "H": H, "Bphi": Bphi, "r": r, "theta": theta}
    h = store.save_run(params, scalars_json, fields)
    st.success(f"corrida salva: {h}")
