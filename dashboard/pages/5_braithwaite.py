"""Tab 5 — Braithwaite: stability of the magnetized equilibrium via
dynamical relaxation of a random field (Braithwaite & Nordlund 2006),
not via imposing a chosen Bt/Bp a priori. Orchestrates the Castro side of
the project the way Tab 2 orchestrates the SCF sweep -- this tab computes
nothing itself (R1); it selects/records parameters and, once Phase 0
unblocks it, launches/monitors Castro and displays what Castro produced.

Honest skeleton (per the closeout prompt that created this tab): every
section that depends on Castro is disabled and labeled with what's
missing, not mocked. If there is no data, the section says so -- it does
not fabricate a plausible-looking placeholder number.

Plan (for context; Steps 2-5 are NOT implemented here, only Step 1 is
live):
  Step 0 -- Castro builds, USE_MHD=TRUE, validated by Orszag-Tang.
            The current blocker. See docs/investigations/ or the
            installation checklist for where that stands.
  Step 1 -- field-free background star from the SCF, exported via Tab 3.
            Live today -- this tab's section 1.
  Step 2 -- seed a RANDOM vector potential in the interior (multi-scale,
            confined to the star, amplitude->0 before the surface,
            edge-centered for machine-precision div B=0, E_mag/|W|~0.1-0.2,
            no imposed symmetry). New Castro problem-setup code, written
            once Phase 0 unblocks it. This tab's section 2 records the
            configuration for that future code to consume -- it does not
            generate anything yet.
  Step 3 -- evolve in Castro, watch the random field relax over a few
            Alfven times.
  Step 4 -- measure: surviving Bt/Bp (the central result), retained
            magnetic energy fraction, resulting surface dipole (the
            collaborator deliverable -- a relaxed field has a poloidal
            component and hence an exterior dipole, unlike pure toroidal),
            final mass.
  Step 5 -- mandatory resolution study, 128^3 vs 256^3 -- in ideal MHD,
            reconnection is controlled by numerical resistivity, so the
            result only holds if it's stable under refinement.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

_DASHBOARD_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DASHBOARD_DIR))

import store
import units

st.set_page_config(page_title="Braithwaite — wd-magnetizada", layout="wide")
st.title("Tab 5 — Braithwaite (stability via dynamical relaxation)")

st.info(
    "Why this method, not imposing Bt/Bp: instead of picking a Bt/Bp ratio "
    "and testing whether it's stable, a random field is released and the "
    "3D dynamics picks the configuration that survives -- Braithwaite & "
    "Nordlund (2006)'s approach for Ap stars. It answers *which* mixed "
    "field is stable, not *whether* a chosen one is. This runs on Castro "
    "(3D dynamical evolution) -- there is no equivalent in the 2D "
    "equilibrium SCF. Every section below that needs Castro is disabled "
    "until Phase 0 (Castro build, validated by Orszag-Tang) unblocks it."
)

# ---------------------------------------------------------------------
# Step 1 — background star (live: this is just an equilibrium, already works)
# ---------------------------------------------------------------------
st.header("1. Background star")
st.caption(
    "A field-free equilibrium from the SCF (Tab 1/2), exported via Tab 3. "
    "This step needs nothing from Castro -- it works today."
)

idx = store.load_index()
if idx.empty:
    st.warning("No saved runs yet. Go to Tab 1 or 2, converge a field-free "
               "equilibrium (field=none, rotation=none), and save it.")
    field_free_idx = pd.DataFrame()
else:
    def _col(df, name, default=0.0):
        return df[name] if name in df.columns else pd.Series(default, index=df.index)

    k0_col = _col(idx, "param_k0").fillna(0.0)
    K_col = _col(idx, "param_K_tor").fillna(0.0)
    Om_col = _col(idx, "param_Omega_c").fillna(0.0)
    field_free_mask = (k0_col == 0.0) & (K_col == 0.0) & (Om_col == 0.0)
    field_free_idx = idx[field_free_mask]

if field_free_idx.empty:
    st.warning(
        "No field-free, non-rotating run found in the registry (checked "
        "param_k0, param_K_tor, param_Omega_c all == 0). Braithwaite needs "
        "a plain background star to seed the random field into -- go make "
        "one in Tab 1 (field=none, rotation=none) and save it."
    )
    background_hash = None
else:
    background_hash = st.selectbox(
        "background star (field-free equilibrium)",
        field_free_idx["hash"].tolist(),
        format_func=lambda h: (
            f"{h}  (rho_c={field_free_idx.set_index('hash').loc[h, 'param_rho_c']:.3e} g/cm³, "
            f"M={field_free_idx.set_index('hash').loc[h, 'M/M_sun']:.4f} M_sun)"
        ) if "param_rho_c" in field_free_idx.columns and "M/M_sun" in field_free_idx.columns else h,
    )
    run = store.load_run(background_hash)
    R_eq_km = run["scalars"].get("R_eq (km)", float("nan"))
    R_pol_km = run["scalars"].get("R_pol (km)", float("nan"))
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("M/M_sun", f"{run['scalars'].get('M/M_sun', float('nan')):.4f}")
    c2.metric("rho_c (g/cm³)", f"{run['params'].get('rho_c', float('nan')):.3e}")
    c3.metric("R_eq (km)", f"{R_eq_km:.2f}" if np.isfinite(R_eq_km) else "n/a")
    c4.metric("VE", f"{run['scalars'].get('VE', float('nan')):.3e}")
    if run["scalars"].get("VE", 1.0) >= 1e-3:
        st.error("This run's VE >= 1e-3 -- not a reliable equilibrium. Pick a "
                  "different one (same R5 rule as Tab 3's export gate).")

st.divider()

# ---------------------------------------------------------------------
# Step 2 — random field seed CONFIGURATION (not the generator)
# ---------------------------------------------------------------------
st.header("2. Random field seed (Step 2)")
st.caption(
    "The generator (multi-scale vector potential, edge-centered for "
    "machine-precision div B=0, confined to the star) is new Castro "
    "problem-setup code that does not exist yet -- it is the first "
    "focused prompt once Phase 0 unblocks. This section only records the "
    "PARAMETERS that generator will eventually consume, with provenance "
    "(git commit, RNG seed), so nothing has to be redesigned when it's "
    "written."
)

sc1, sc2 = st.columns(2)
E_mag_target = sc1.slider(
    "E_mag/|W| target", 0.05, 0.30, 0.15, 0.01,
    help="Plan's stated range: ~0.1-0.2.")
n_modes = sc2.slider(
    "number of modes (spectrum breadth)", 1, 50, 10,
    help="How many independent scales the random vector potential mixes. "
         "Not yet meaningful without the generator -- placeholder for its "
         "eventual spectrum control.")
dominant_scale_over_R = sc1.slider(
    "dominant scale / R_star", 0.05, 1.0, 0.3, 0.05,
    help="Characteristic scale of the field pattern, as a fraction of the "
         "stellar radius.")
rng_seed = sc2.number_input(
    "RNG seed (required)", min_value=0, value=42, step=1,
    help="Recorded with the configuration -- an unrecorded seed makes the "
         "eventual run unreproducible.")

st.button(
    "generate seed (requires Castro)", disabled=True,
    help="Disabled: the generator lives in Castro's problem setup and does "
         "not exist until Phase 0 (Castro build) unblocks it. This button "
         "will call that generator once it's written -- not reimplemented "
         "here (R1)."
)

if background_hash is not None:
    if st.button("save this seed configuration"):
        cfg = {
            "background_run_hash": background_hash,
            "E_mag_over_W_target": E_mag_target,
            "n_modes": n_modes,
            "dominant_scale_over_R": dominant_scale_over_R,
            "rng_seed": int(rng_seed),
        }
        h = store.save_braithwaite_seed_config(cfg)
        st.success(f"Configuration saved: `{h}` (not yet consumed -- no generator to run it through)")
else:
    st.caption("Pick a background star above before saving a seed configuration.")

saved_configs = store.load_braithwaite_seed_configs()
if not saved_configs.empty:
    st.markdown("**Saved seed configurations**")
    st.dataframe(saved_configs, width="stretch")
else:
    st.caption("No seed configurations saved yet.")

st.divider()

# ---------------------------------------------------------------------
# Steps 3-4 — evolution and diagnostics (all pending Castro)
# ---------------------------------------------------------------------
st.header("3. Evolution and diagnostics (Steps 3-4)")
st.caption(
    "Everything below needs a working Castro build (Phase 0) plus the "
    "Step 2 generator. Placeholders keep each result's future position "
    "fixed rather than having this tab redesigned around the first "
    "result."
)

st.subheader("Launch / monitor")
lc1, lc2 = st.columns(2)
lc1.button("launch Castro evolution run", disabled=True,
           help="Pending Castro build (Phase 0) and a generated seed (Step 2).")
lc2.button("refresh run status", disabled=True,
           help="Pending an actual running Castro job to monitor.")
st.caption("aguardando build do Castro (Phase 0)")

st.subheader("Outcome diagnostics (Step 4)")
d1, d2, d3, d4 = st.columns(4)
for col, label in zip((d1, d2, d3, d4),
                       ("Bt/Bp (survival)", "E_mag retained (%)",
                        "surface dipole (G)", "final M/M_sun")):
    with col:
        st.metric(label, "—")
        st.caption("pending Castro")
st.caption(
    "Surface dipole is the collaborator deliverable this step produces: a "
    "relaxed mixed field has a poloidal component and hence an exterior "
    "dipole for magnetic braking to act on, which the pure self-consistent "
    "toroidal branch (docs/teoria.md Sec 6.2) does not have."
)

st.divider()

# ---------------------------------------------------------------------
# Step 5 — resolution study (pending Castro, and pending Steps 3-4 output)
# ---------------------------------------------------------------------
st.header("4. Resolution study (Step 5)")
st.warning(
    "Mandatory, not optional: in ideal MHD, reconnection during the "
    "relaxation is controlled by numerical resistivity, not physical "
    "resistivity -- so the surviving Bt/Bp, retained energy, and surface "
    "dipole above are only trustworthy if they are stable under mesh "
    "refinement. A single-resolution result is not a result."
)
rc1, rc2 = st.columns(2)
rc1.metric("128³ result", "—")
rc2.metric("256³ result", "—")
st.caption("aguardando build do Castro (Phase 0) e uma primeira run das Etapas 3-4")
