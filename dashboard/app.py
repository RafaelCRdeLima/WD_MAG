"""Control dashboard for the magnetized SCF. Usage: streamlit run dashboard/app.py"""

import streamlit as st

st.set_page_config(page_title="wd-magnetizada — SCF", layout="wide")

st.title("Magnetized white dwarf — SCF control panel")
st.markdown(
    """
Use the sidebar menu to navigate between tabs:

- **Equilibrium** — single SCF run, inspect one point (rho_c, k0)
- **Sweep** — parameter grid, M-R diagram
- **Export** — generates initial data (HDF5) and Castro `inputs`
- **Runs** — run history, comparison, references
- **Braithwaite** — stability via dynamical relaxation of a random field
  (Castro, 3D) — orchestration only; most of it is disabled pending the
  Castro build (Phase 0)

All physics comes from `scf.*` (`eos`, `poisson`, `gradshafranov`, `scf`,
`diagnostics`, `toroidal`) — this dashboard only explores, persists, and
exports. See `plano_wd_magnetizada.md` for the full plan and project
decisions.
"""
)

st.info(
    "Project note: the plan's original recipe (fixing the Bernoulli "
    "constant via H=0 at the surface) turned out to be unstable for this "
    "EOS near the Chandrasekhar limit. The SCF here uses the (rho_c, k0) "
    "parametrization — central density and field amplitude, as two "
    "independent inputs — which is stable and validated. See scf/scf.py."
)
