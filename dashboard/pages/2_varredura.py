"""Aba 2 — Varredura: grade de parametros (rho_c, k0), execucao paralela,
cache por hash, diagrama M-R. R1: fisica so' via scf.* (aqui, sweep_worker)."""

import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import streamlit as st

_DASHBOARD_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DASHBOARD_DIR))
sys.path.insert(0, str(_DASHBOARD_DIR.parent / "scf"))

import store
import seed
import sweep_worker

st.set_page_config(page_title="Varredura — wd-magnetizada", layout="wide")
st.title("Aba 2 — Varredura")

M_SUN = 1.989e33
REF_DIR = _DASHBOARD_DIR / "data" / "referencias"

st.sidebar.header("Grade")
st.sidebar.markdown("**ρc (g/cm³)**")
rho_c_min = st.sidebar.number_input("ρc mín", value=1e8, format="%.2e")
rho_c_max = st.sidebar.number_input("ρc máx", value=1e12, format="%.2e")
n_rho = st.sidebar.slider("N (ρc)", 2, 20, 8)

st.sidebar.markdown("**k0**")
incluir_k0_zero = st.sidebar.checkbox("incluir k0=0", value=True)
k0_min = st.sidebar.number_input("k0 mín (>0)", value=1e-16, format="%.2e")
k0_max = st.sidebar.number_input("k0 máx", value=1e-12, format="%.2e")
n_k0 = st.sidebar.slider("N (k0>0)", 1, 20, 6)

mu_e = st.sidebar.number_input("μₑ", min_value=1.0, max_value=2.5, value=2.0, step=0.1)
Nr = st.sidebar.select_slider("Nr", options=[65, 129, 161], value=129)
Ntheta = st.sidebar.select_slider("Ntheta", options=[33, 65, 129], value=65)
lmax = st.sidebar.slider("l_max", 4, 32, 16)
n_workers = st.sidebar.slider("processos paralelos", 1, 8, 4)

st.caption(
    "Varre (ρc, k0) — zeta (toroidal) fica de fora da grade por enquanto "
    "(cada ponto aqui é uma configuração *poloidal*; o toroidal é imposto "
    "depois, por cima, na Aba 3)."
)

rho_c_grid = np.geomspace(rho_c_min, rho_c_max, n_rho)
k0_grid = list(np.geomspace(k0_min, k0_max, n_k0))
if incluir_k0_zero:
    k0_grid = [0.0] + k0_grid

total_points = len(rho_c_grid) * len(k0_grid)
st.write(f"Grade: {len(rho_c_grid)} × {len(k0_grid)} = {total_points} pontos")

if st.button("rodar varredura", type="primary"):
    param_list = []
    for rho_c in rho_c_grid:
        for k0 in k0_grid:
            param_list.append({
                "rho_c": float(rho_c), "k0": float(k0), "mu_e": mu_e,
                "R_guess": seed.r_guess(rho_c), "Nr": Nr, "Ntheta": Ntheta,
                "lmax": lmax, "tol": 1e-6, "max_iter": 200,
            })

    cache_params_key = lambda p: {k: v for k, v in p.items() if k != "R_guess"}

    to_run, cached_results = [], []
    for p in param_list:
        h = store.run_exists(cache_params_key(p))
        if h:
            cached_results.append((p, h))
        else:
            to_run.append(p)

    st.write(f"{len(cached_results)} já em cache, {len(to_run)} a calcular")

    progress = st.progress(0.0)
    status = st.empty()
    n_done, n_ok, n_fail = 0, 0, 0
    rows = []
    fail_rows = []

    for p, h in cached_results:
        run = store.load_run(h)
        rows.append({"hash": h, **run["params"], **run["scalars"]})
        n_ok += 1
        n_done += 1

    if to_run:
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            futures = {ex.submit(sweep_worker.run_one, p): p for p in to_run}
            for fut in as_completed(futures):
                p = futures[fut]
                try:
                    res = fut.result()
                except Exception as e:
                    res = {"converged": False, "rho_c": p["rho_c"], "k0": p["k0"], "error": str(e)}

                n_done += 1
                if res["converged"]:
                    h = store.save_run(cache_params_key(p), res["scalars"], res["fields"])
                    rows.append({"hash": h, **cache_params_key(p), **res["scalars"]})
                    n_ok += 1
                else:
                    fail_rows.append({"rho_c": p["rho_c"], "k0": p["k0"]})
                    n_fail += 1

                progress.progress(n_done / total_points)
                status.text(f"{n_done}/{total_points}  (convergiu: {n_ok}, falhou: {n_fail})")

    st.session_state["sweep_rows"] = rows
    st.session_state["sweep_fail_rows"] = fail_rows
    st.success(f"varredura concluída: {n_ok}/{total_points} convergiram "
               f"({n_ok / total_points:.0%})")

rows = st.session_state.get("sweep_rows", [])
fail_rows = st.session_state.get("sweep_fail_rows", [])

if not rows:
    st.info("Rode a varredura para ver os gráficos.")
    st.stop()

df = pd.DataFrame(rows)

if fail_rows:
    with st.expander(f"{len(fail_rows)} ponto(s) que não convergiram"):
        st.dataframe(pd.DataFrame(fail_rows))

st.subheader("Diagrama M-R")
fig_mr = px.scatter(df, x="R_eq (km)", y="M/M☉", color="k0",
                     color_continuous_scale="Viridis",
                     hover_data=["hash", "rho_c", "VE"])
fig_mr.add_hline(y=1.44, line_dash="dash", line_color="gray",
                  annotation_text="limite de Chandrasekhar (μₑ=2, sem campo)")
ref_file = REF_DIR / "bera_bhattacharya_2014.csv"
if ref_file.exists():
    ref_df = pd.read_csv(ref_file)
    fig_mr.add_trace(go.Scatter(x=ref_df["R_km"], y=ref_df["M_Msun"], mode="markers",
                                 marker=dict(symbol="x", color="red", size=8),
                                 name="Bera & Bhattacharya 2014"))
else:
    st.caption(f"(sem overlay de literatura — arquivo não encontrado: {ref_file})")
event_mr = st.plotly_chart(fig_mr, on_select="rerun", key="mr_chart")

st.subheader("M vs ρc, colorido por E_mag/|W|")
fig_m_rhoc = px.scatter(df, x="rho_c", y="M/M☉", color="E_mag/|W|",
                         color_continuous_scale="Inferno", log_x=True,
                         hover_data=["hash", "k0", "VE"])
st.plotly_chart(fig_m_rhoc, key="m_rhoc_chart")

st.subheader("Mapa de calor de VE sobre a grade")
if len(df["k0"].unique()) > 1 and len(df["rho_c"].unique()) > 1:
    pivot = df.pivot_table(index="k0", columns="rho_c", values="VE", aggfunc="mean")
    fig_ve = go.Figure(data=go.Heatmap(
        z=np.log10(pivot.values), x=pivot.columns, y=pivot.index,
        colorscale="RdBu_r", colorbar=dict(title="log10(VE)")))
    fig_ve.update_xaxes(type="log", title="ρc")
    fig_ve.update_yaxes(type="log", title="k0")
    st.plotly_chart(fig_ve, key="ve_heatmap")
else:
    st.caption("grade pequena demais para mapa de calor 2D")

st.subheader("Carregar um ponto na Aba 1")
selected_hash = st.selectbox("ou escolha pelo hash", df["hash"].tolist())
if event_mr and event_mr.get("selection", {}).get("points"):
    pt = event_mr["selection"]["points"][0]
    idx = pt.get("point_index")
    if idx is not None and idx < len(df):
        selected_hash = df.iloc[idx]["hash"]
        st.write(f"selecionado no gráfico: {selected_hash}")

if st.button("carregar equilíbrio selecionado na Aba 1"):
    run = store.load_run(selected_hash)
    st.session_state["reload_run_params"] = {
        "rho_c": run["params"]["rho_c"], "k0": run["params"]["k0"],
        "mu_e": run["params"].get("mu_e", 2.0), "zeta_target_ratio": 0.0, "m_tor": 1,
        "Nr": run["params"].get("Nr", 129), "Ntheta": run["params"].get("Ntheta", 129),
        "lmax": run["params"].get("lmax", 16), "tol": run["params"].get("tol", 1e-6),
        "max_iter": run["params"].get("max_iter", 200),
    }
    st.switch_page("pages/1_equilibrio.py")
