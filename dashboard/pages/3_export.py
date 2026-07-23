"""Tab 3 — Export: toroidal (D6), derived scales, torus resolution check,
generates HDF5 + Castro inputs + run_manifest.json.

R3: only writes files, never launches Castro runs.
R5: export blocked if VE > 1e-3, no option to force it.
"""

import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import streamlit as st

_DASHBOARD_DIR = Path(__file__).resolve().parent.parent
_REPO_ROOT = _DASHBOARD_DIR.parent
sys.path.insert(0, str(_REPO_ROOT / "scf"))
sys.path.insert(0, str(_DASHBOARD_DIR))

import diagnostics as diag
import toroidal as tor
import units
import store
from eos import sound_speed, x_of_enthalpy

st.set_page_config(page_title="Export — wd-magnetizada", layout="wide")
st.title("Tab 3 — Export to Castro")

idx = store.load_index()
if idx.empty:
    st.info("No saved runs. Go to Tab 1 or 2 and save an equilibrium first.")
    st.stop()

run_hash = st.selectbox("equilibrium to export", idx["hash"].tolist())
run = store.load_run(run_hash)
rho, Phi, u, H = run["fields"]["rho"], run["fields"]["Phi"], run["fields"]["u"], run["fields"]["H"]
r, theta = run["fields"]["r"], run["fields"]["theta"]
mu_e = run["params"].get("mu_e", 2.0)
rho_c = run["params"]["rho_c"]

Br, Bth = diag.poloidal_field(u, r, theta)
M = diag.volume_integral(rho, r, theta)
R_eq, R_pol = diag.equatorial_polar_radii(rho, r, theta)

st.subheader("Toroidal field (D6)")
c1, c2 = st.columns(2)
zeta_target = c1.slider("Bt/Bp target (energy)", 0.0, 2.0, 0.3, 0.05, key="export_zeta_target")
m_tor = c2.slider("m_tor", 1, 4, 1, key="export_mtor")

Bphi = np.zeros_like(rho)
u_c = None
if zeta_target > 0:
    try:
        Bphi, zeta_used, u_c = tor.solve_zeta_for_energy_ratio(u, rho, r, theta, zeta_target, m_tor=m_tor)
        ratio_energy, ratio_amp = tor.bt_bp_ratios(Br, Bth, Bphi, r, theta)
        cc1, cc2, cc3 = st.columns(3)
        cc1.metric("Bt/Bp (energy)", f"{ratio_energy:.4f}")
        cc2.metric("Bt/Bp (amplitude)", f"{ratio_amp:.4f}")
        cc3.metric("zeta", f"{zeta_used:.3e}")
        st.caption(
            f"Continuity of beta*beta' at the torus boundary: guaranteed "
            f"analytically for m_tor>=1 (B_phi proportional to "
            f"(u-u_c)^{{m_tor+1}}, and its derivative goes to 0 at u=u_c). "
            f"m_tor={m_tor} ✓"
        )
    except ValueError as e:
        st.warning(str(e))
else:
    st.caption("Bt/Bp target = 0 → no toroidal field.")

_, _, E_mag = diag.magnetic_energies(Br, Bth, Bphi, r, theta)
W = diag.gravitational_energy(rho, Phi, r, theta)
Pi = diag.pressure_integral(H, r, theta, mu_e)
VE, _, _, _ = diag.virial_error(rho, Phi, H, Br, Bth, Bphi, r, theta, mu_e)

st.subheader("Derived scales")
rho_mean_est = M / (4.0 / 3.0 * np.pi * R_eq**2 * R_pol) if R_eq > 0 else float("nan")
B_mean = np.sqrt(diag.volume_integral(Br**2 + Bth**2 + Bphi**2, r, theta) / (4.0 / 3.0 * np.pi * R_eq**2 * R_pol)) if R_eq > 0 else 0.0
t_dyn = units.dynamical_time(M, R_eq)
v_A = units.alfven_speed(B_mean, rho_mean_est)
t_alf = units.alfven_time(R_eq, v_A) if v_A > 0 else float("inf")
alf_over_dyn = t_alf / t_dyn if t_dyn > 0 else float("nan")

d0, d1, d2, d3, d4 = st.columns(5)
d0.metric("mean B", units.format_gauss(B_mean))
d1.metric("t_dyn (s)", f"{t_dyn:.3e}")
d2.metric("v_A (cm/s)", f"{v_A:.3e}")
d3.metric("t_Alfven (s)", f"{t_alf:.3e}" if np.isfinite(t_alf) else "inf")
d4.metric("t_Alf/t_dyn", f"{alf_over_dyn:.3f}" if np.isfinite(alf_over_dyn) else "inf")
if 0.3 < alf_over_dyn < 3:
    st.success("t_Alf/t_dyn ~ 1 → strong-field regime (D4).")

st.subheader("Torus resolution check")
n_cell = st.select_slider("cells per side (Castro box)", options=[128, 256, 384], value=128)
box_size_Req = st.number_input("box size (in R_eq)", value=4.0, step=0.5)
box_half = box_size_Req * R_eq / 2.0
cell_size = (2 * box_half) / n_cell

if u_c is not None:
    r_in, r_out, thickness = tor.torus_radial_extent(u, rho, r, theta, u_c)
    n_cells_torus = thickness / cell_size if cell_size > 0 else 0
    st.write(f"Radial thickness of the torus (at the equator): {units.format_km(thickness)} "
             f"(Castro cell: {units.format_km(cell_size)}) "
             f"→ **{n_cells_torus:.1f} cells** at {n_cell}³ resolution")
    if n_cells_torus < 10:
        st.error(f"< 10 cells across the torus — {n_cell}³ resolution is "
                 f"probably insufficient for this configuration.")
    else:
        st.success(f">= 10 cells across the torus ({n_cells_torus:.1f}).")
else:
    st.caption("no toroidal field — torus resolution check does not apply.")

st.subheader("Box parameters")
small_dens = rho_c * 1e-8
sponge_radius = 1.5 * R_eq
t_sponge = 5 * t_dyn
x_c = x_of_enthalpy(H[0, 0], mu_e)
c_s_center = sound_speed(x_c, mu_e)
perturb_amp = 1e-4 * c_s_center
seed_rng = st.number_input("RNG seed (required)", min_value=0, value=42, step=1)
if np.isfinite(t_alf):
    stop_time_alf = st.number_input("stop_time (in t_Alfven)", value=5.0, step=1.0)
    stop_time_s = stop_time_alf * t_alf
else:
    st.caption("no field (t_Alfven = inf) — stop_time set in t_dyn (Phase 4, hydro only).")
    stop_time_alf = st.number_input("stop_time (in t_dyn)", value=10.0, step=1.0)
    stop_time_s = stop_time_alf * t_dyn

box_params_display = {
    "castro.small_dens (g/cm³)": f"{small_dens:.4e}",
    "sponge start radius (km)": units.format_km(sponge_radius),
    "damping duration (s, ~5 t_dyn)": f"{t_sponge:.4e}",
    "perturbation amplitude (cm/s, 1e-4 c_s)": f"{perturb_amp:.4e}",
    "RNG seed": str(int(seed_rng)),
    "stop_time (s)": f"{stop_time_s:.4e}",
}
st.table({"parameter": list(box_params_display.keys()),
          "value": list(box_params_display.values())})

# ---------------- export (R5: blocked if VE > 1e-3) ----------------
st.divider()
st.subheader("Export")
if VE >= 1e-3:
    st.error(f"Export disabled: VE = {VE:.3e} >= 1e-3 (plan's V3). "
             f"This equilibrium is not reliable enough to become Castro "
             f"initial data. There is no option to force it.")
else:
    st.success(f"VE = {VE:.3e} < 1e-3 — export enabled.")
    if st.button("generate HDF5 + inputs + manifest", type="primary"):
        # source_run_* preserves the parameters of the loaded poloidal
        # equilibrium (see source_run_hash in the manifest); export_* are
        # this tab's own controls (toroidal imposed here, Castro box) —
        # distinct names on purpose, avoids colliding with run["params"]["k0"] etc.
        export_params = {
            f"source_run_{k}": v for k, v in run["params"].items()
        }
        export_params.update({
            "export_zeta_target_ratio": zeta_target, "export_m_tor": m_tor,
            "n_cell": n_cell, "box_size_Req": box_size_Req,
            "seed": int(seed_rng), "stop_time_in_units": stop_time_alf,
        })
        export_hash = store.params_hash(export_params)
        out_dir = _DASHBOARD_DIR / "exports" / export_hash
        out_dir.mkdir(parents=True, exist_ok=True)

        A_phi = np.divide(u, r[:, None] * np.sin(theta)[None, :],
                           out=np.zeros_like(u), where=(r[:, None] * np.sin(theta)[None, :]) > 0)

        h5_path = out_dir / "initial_data.h5"
        with h5py.File(h5_path, "w") as f:
            f.create_dataset("r", data=r)
            f.create_dataset("theta", data=theta)
            f.create_dataset("rho", data=rho)
            f.create_dataset("A_phi", data=A_phi)
            f.create_dataset("B_phi", data=Bphi)
            f.attrs["M_Msun"] = M / units.M_SUN
            f.attrs["R_eq_cm"] = R_eq
            f.attrs["R_pol_cm"] = R_pol
            f.attrs["rho_c_gcm3"] = rho_c
            f.attrs["E_mag_erg"] = E_mag
            f.attrs["VE"] = VE
            f.attrs["units"] = "Gaussian CGS; B in gauss (not B'=B/sqrt(4pi) as used by Castro)"

        inputs_template = (_DASHBOARD_DIR / "templates" / "inputs.template").read_text()
        inputs_text = inputs_template.format(
            stop_time=stop_time_s, box_half=box_half, n_cell=n_cell,
            small_dens=small_dens, sponge_radius=sponge_radius, sponge_timescale=t_sponge,
            seed=int(seed_rng), perturb_amplitude=perturb_amp,
            initial_data_file=str(h5_path),
        )
        (out_dir / "inputs").write_text(inputs_text)

        def _git_hash(path):
            try:
                out = subprocess.run(["git", "-C", str(path), "rev-parse", "HEAD"],
                                      capture_output=True, text=True, timeout=5)
                return out.stdout.strip() if out.returncode == 0 else "no-git"
            except Exception:
                return "no-git"

        manifest = {
            "export_hash": export_hash, "source_run_hash": run_hash,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "git_commit_scf_dashboard": _git_hash(_REPO_ROOT),
            "params": export_params,
            "derived": {"t_dyn_s": t_dyn, "v_A_cms": v_A, "t_Alfven_s": t_alf,
                        "t_Alf_over_t_dyn": alf_over_dyn, "VE": VE},
            "paths": {"hdf5": str(h5_path), "inputs": str(out_dir / "inputs")},
        }
        (out_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2, default=float))

        st.success(f"Exported to `{out_dir}`: `initial_data.h5`, `inputs`, `run_manifest.json`")
