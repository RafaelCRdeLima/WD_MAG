"""Aba 4 — Registro de corridas: tabela, comparacao, recarregar, marcar referencia."""

import sys
from pathlib import Path

import streamlit as st

_DASHBOARD_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DASHBOARD_DIR))

import store
import units

st.set_page_config(page_title="Registro — wd-magnetizada", layout="wide")
st.title("Aba 4 — Registro de corridas")

idx = store.load_index()
if idx.empty:
    st.info("Nenhuma corrida salva ainda. Vá para a Aba 1, rode um equilíbrio e clique "
            "em \"salvar esta corrida\".")
    st.stop()

st.subheader(f"{len(idx)} corrida(s)")
# R4: campo em gauss (notacao cientifica), raios em km (2 casas) — as colunas
# ja guardam os valores nessas unidades (Aba 1/2); aqui so' formata a exibicao.
_column_config = {"reference": st.column_config.CheckboxColumn("referência")}
for col in idx.columns:
    if col.endswith("(G)"):
        _column_config[col] = st.column_config.NumberColumn(col, format="%.3e")
    elif col.endswith("(km)"):
        _column_config[col] = st.column_config.NumberColumn(col, format="%.2f")

edited = st.data_editor(
    idx, width="stretch", num_rows="fixed",
    column_config=_column_config,
    disabled=[c for c in idx.columns if c != "reference"],
    key="runs_editor",
)

if not edited["reference"].equals(idx["reference"]):
    changed = edited[edited["reference"] != idx["reference"]]
    for _, row in changed.iterrows():
        store.mark_reference(row["hash"], bool(row["reference"]))
    st.rerun()

st.divider()
st.subheader("Comparação lado a lado")
hashes = idx["hash"].tolist()
c1, c2 = st.columns(2)
h1 = c1.selectbox("corrida A", hashes, index=0)
h2 = c2.selectbox("corrida B", hashes, index=min(1, len(hashes) - 1))

run1 = store.load_run(h1)
run2 = store.load_run(h2)


def _format_scalars(scalars: dict) -> dict:
    """Mesma regra R4 da Aba 1: gauss em notacao cientifica, km com 2 casas."""
    out = {}
    for k, v in scalars.items():
        if not isinstance(v, (int, float)):
            out[k] = v
        elif k.endswith("(G)"):
            out[k] = units.format_gauss(v)
        elif k.endswith("(km)"):
            out[k] = units.format_km_value(v)
        else:
            out[k] = v
    return out


comp_col1, comp_col2 = st.columns(2)
with comp_col1:
    st.markdown(f"**{h1}**")
    st.json(run1["params"])
    st.json(_format_scalars(run1["scalars"]))
    if st.button("recarregar A na Aba 1", key="reload_a"):
        st.session_state["reload_run_params"] = run1["params"]
        st.switch_page("pages/1_equilibrio.py")
with comp_col2:
    st.markdown(f"**{h2}**")
    st.json(run2["params"])
    st.json(_format_scalars(run2["scalars"]))
    if st.button("recarregar B na Aba 1", key="reload_b"):
        st.session_state["reload_run_params"] = run2["params"]
        st.switch_page("pages/1_equilibrio.py")
