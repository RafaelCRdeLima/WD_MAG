"""Dashboard de controle do SCF magnetizado. Uso: streamlit run dashboard/app.py"""

import streamlit as st

st.set_page_config(page_title="wd-magnetizada — SCF", layout="wide")

st.title("Anã branca magnetizada — painel de controle do SCF")
st.markdown(
    """
Use o menu lateral para navegar entre as abas:

- **Equilíbrio** — execução única do SCF, inspeção de um ponto (ρc, k0)
- **Varredura** — grade de parâmetros, diagrama M-R
- **Exportação** — gera dado inicial (HDF5) e `inputs` do Castro
- **Registro** — histórico de corridas, comparação, referências

Toda a física vem de `scf.*` (`eos`, `poisson`, `gradshafranov`, `scf`,
`diagnostics`, `toroidal`) — este dashboard só explora, persiste e exporta.
Ver `plano_wd_magnetizada.md` para o plano completo e as decisões de projeto.
"""
)

st.info(
    "Nota de projeto: a receita original do plano (fixar a constante de "
    "Bernoulli via H=0 na superfície) mostrou-se instável para esta EOS "
    "perto do limite de Chandrasekhar. O SCF aqui usa a parametrização "
    "(ρc, k0) — densidade central e amplitude do campo, como duas entradas "
    "independentes — que é estável e validada. Ver scf/scf.py."
)
