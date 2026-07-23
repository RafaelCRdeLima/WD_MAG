"""Aba 4 — Registro de corridas: tabela, comparacao, recarregar, marcar referencia."""

import sys
from pathlib import Path

import streamlit as st

_DASHBOARD_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DASHBOARD_DIR))

import store

st.set_page_config(page_title="Registro — wd-magnetizada", layout="wide")
st.title("Aba 4 — Registro de corridas")

idx = store.load_index()
if idx.empty:
    st.info("Nenhuma corrida salva ainda. Vá para a Aba 1, rode um equilíbrio e clique "
            "em \"salvar esta corrida\".")
    st.stop()

st.subheader(f"{len(idx)} corrida(s)")
edited = st.data_editor(
    idx, width="stretch", num_rows="fixed",
    column_config={"reference": st.column_config.CheckboxColumn("referência")},
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

comp_col1, comp_col2 = st.columns(2)
with comp_col1:
    st.markdown(f"**{h1}**")
    st.json(run1["params"])
    st.json(run1["scalars"])
    if st.button("recarregar A na Aba 1", key="reload_a"):
        st.session_state["reload_run_params"] = run1["params"]
        st.switch_page("pages/1_equilibrio.py")
with comp_col2:
    st.markdown(f"**{h2}**")
    st.json(run2["params"])
    st.json(run2["scalars"])
    if st.button("recarregar B na Aba 1", key="reload_b"):
        st.session_state["reload_run_params"] = run2["params"]
        st.switch_page("pages/1_equilibrio.py")
