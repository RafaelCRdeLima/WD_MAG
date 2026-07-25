"""Tab 2 — Sweep: parameter grid over up to two selectable axes (rho_c,
k0, K, Omega_c), parallel execution, cache by hash. R1: physics only via
scf.* (here, sweep_worker).

Mode selectors (rotation / field) mirror Tab 1 exactly -- same strings,
same mutual-exclusion rules (poloidal XOR toroidal self-consistent;
rotation none/rigid/differential). For differential rotation, A/R_eq is
ALWAYS a fixed value, never a sweep axis -- opening it as an axis would
explode the grid (3 physical axes x mesh resolution is already a lot of
SCF solves)."""

import io
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import streamlit as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatterSciNotation, ScalarFormatter, NullFormatter

_DASHBOARD_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DASHBOARD_DIR))
sys.path.insert(0, str(_DASHBOARD_DIR.parent / "scf"))

import store
import seed
import sweep_worker
import eos

st.set_page_config(page_title="Sweep — wd-magnetizada", layout="wide")
st.title("Tab 2 — Sweep")

REF_DIR = _DASHBOARD_DIR / "data" / "references"

AXIS_LABELS = {"rho_c": "rho_c (g/cm³)", "k0": "k0", "K": "K (toroidal)", "Omega_c": "Omega_c (rad/s)"}


def _log_rho_axis(fig, axis="x"):
    """Genuine log10 axis for rho_c: major ticks at every power of ten
    (dtick=1 on a log axis), labels in scientific notation (10^8, 10^9,
    ...). Plotly's default exponentformat ("B", SI-style byte/currency
    suffixes) turns these into "100M", "1B", "10B" on the tick labels and
    "0.2T"/"0.4T"/"0.6T" on heatmap gridlines -- reads as data size or
    money, not a log-scaled physical quantity, and was never set to
    "power" explicitly on any chart in this file."""
    kw = dict(type="log", dtick=1, exponentformat="power", title="rho_c (g/cm³)")
    (fig.update_xaxes if axis == "x" else fig.update_yaxes)(**kw)


def _msun_axis(fig, axis="y"):
    """M/M_sun tick labels: plain notation, explicit so this axis is
    never subject to the same SI-suffix default as the rho_c axes on
    the same figures (M/M_sun values are O(1) here and unlikely to
    trigger it in practice, but left implicit it's one parameter change
    away from doing so silently)."""
    kw = dict(exponentformat="none", title="M/M_sun")
    (fig.update_xaxes if axis == "x" else fig.update_yaxes)(**kw)


def _publication_figure_pdf(df, rho_c_neutronization, mu_e):
    """Vector-PDF, matplotlib version of the priority plot -- for the
    collaborator and the paper, not the interactive dashboard. Grouped by
    E_tor/|W| (the portable form, docs/teoria.md Sec 6.2e), not raw K
    (means nothing outside this project's specific
    B_phi=K*rho^m*varpi^(2m-1) ansatz). Astronomy convention: serif
    fonts, single-column-readable size, no embedded title (the figure
    caption belongs in the paper's text, not burned into the image)."""
    rc = {
        "font.family": "serif",
        "font.size": 9,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7,
        "mathtext.fontset": "dejavuserif",
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,  # embed as real text, not paths -- editable/searchable in the PDF
    }
    with plt.rc_context(rc):
        fig, ax = plt.subplots(figsize=(3.5, 3.1))

        df = df.copy()
        df["_ve_group"] = df["E_mag/|W|"].round(2)
        groups = sorted(df.groupby("_ve_group"), key=lambda kv: kv[0])
        cmap = plt.get_cmap("viridis")
        n_nonzero = max(sum(1 for g, _ in groups if g > 0), 1)
        color_i = 0
        for group_val, sub in groups:
            sub = sub.sort_values("rho_c")
            is_zero = group_val <= 0
            if is_zero:
                color, lw, ms, zorder = "black", 1.6, 4.0, 3
                label = "no field (K=0)"
            else:
                color = cmap(0.15 + 0.75 * color_i / max(n_nonzero - 1, 1))
                lw, ms, zorder = 1.1, 3.0, 2
                label = rf"$E_{{\rm tor}}/|W| \approx {group_val:.2f}$"
                color_i += 1
            ax.plot(sub["rho_c"], sub["M/M_sun"], marker="o", markersize=ms,
                     linewidth=lw, color=color, zorder=zorder, label=label)

        ax.axvline(rho_c_neutronization, color="firebrick", linestyle="--",
                    linewidth=1.0, zorder=1)
        ax.text(rho_c_neutronization * 1.15, 0.96,
                 rf"inverse $\beta$-decay threshold ($\mu_e={mu_e:g}$)",
                 transform=ax.get_xaxis_transform(), rotation=90,
                 va="top", ha="left", fontsize=7, color="firebrick")

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$\rho_c\ (\mathrm{g\,cm^{-3}})$")
        ax.set_ylabel(r"$M/M_\odot$")
        # x spans several decades -- genuine power-of-ten scientific
        # notation (10^10, 10^11, ...) is the right call here.
        ax.xaxis.set_major_locator(LogLocator(base=10.0))
        ax.xaxis.set_major_formatter(LogFormatterSciNotation(base=10.0))
        # y (M/M_sun) is O(1) and stays within one decade -- log-scaled
        # SPACING as requested, but "1.4x10^0" reads as over-formatting
        # for a number that small; plain decimal ticks (still
        # log-spaced) is what published mass-density figures actually
        # use for an axis in this range.
        ax.yaxis.set_major_locator(LogLocator(base=10.0, subs=np.arange(1.0, 10.0)))
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(which="both", direction="in", top=True, right=True)
        ax.legend(frameon=False, loc="best", handlelength=1.6)
        fig.tight_layout(pad=0.4)

        buf = io.BytesIO()
        fig.savefig(buf, format="pdf")
        plt.close(fig)
        buf.seek(0)
        return buf


def _axis_grid_ui(axis):
    """Min/max/N controls for `axis` used as a SWEPT dimension. Returns a
    list of values (may include 0.0 for k0/K, matching Tab 1's "field
    off" convention)."""
    if axis == "rho_c":
        lo = st.sidebar.number_input("rho_c min", value=1e8, format="%.2e", key="ax_rhoc_min")
        hi = st.sidebar.number_input("rho_c max", value=1e12, format="%.2e", key="ax_rhoc_max")
        n = st.sidebar.slider("N (rho_c)", 2, 20, 8, key="ax_rhoc_n")
        return list(np.geomspace(lo, hi, n))
    if axis == "k0":
        include_zero = st.sidebar.checkbox("include k0=0", value=True, key="ax_k0_zero")
        lo = st.sidebar.number_input("k0 min (>0)", value=1e-16, format="%.2e", key="ax_k0_min")
        hi = st.sidebar.number_input("k0 max", value=1e-12, format="%.2e", key="ax_k0_max")
        n = st.sidebar.slider("N (k0>0)", 1, 20, 6, key="ax_k0_n")
        grid = list(np.geomspace(lo, hi, n))
        return ([0.0] + grid) if include_zero else grid
    if axis == "K":
        include_zero = st.sidebar.checkbox("include K=0", value=True, key="ax_K_zero")
        lo = st.sidebar.slider("log10(K) min", -8.0, -2.0, -6.0, 0.1, key="ax_K_lomin")
        hi = st.sidebar.slider("log10(K) max", -8.0, -2.0, -3.0, 0.1, key="ax_K_lomax")
        n = st.sidebar.slider("N (K>0)", 1, 12, 5, key="ax_K_n")
        grid = list(10.0 ** np.linspace(min(lo, hi), max(lo, hi), n))
        return ([0.0] + grid) if include_zero else grid
    if axis == "Omega_c":
        lo = st.sidebar.number_input("Omega_c min", value=0.0, step=0.5, key="ax_om_min")
        hi = st.sidebar.number_input("Omega_c max", value=20.0, step=0.5, key="ax_om_max")
        n = st.sidebar.slider("N (Omega_c)", 2, 20, 6, key="ax_om_n")
        return list(np.linspace(lo, hi, n))
    raise ValueError(axis)


def _fixed_value_ui(axis):
    """Single-value control for `axis` when it is NOT swept this run."""
    if axis == "rho_c":
        return st.sidebar.number_input("rho_c (fixed, g/cm³)", value=1e12, format="%.2e", key="fx_rhoc")
    if axis == "k0":
        log_k0 = st.sidebar.slider("log10(|k0|) (fixed)", -20.0, -8.0, -14.0, 0.1, key="fx_k0")
        return 10.0 ** log_k0
    if axis == "K":
        log_K = st.sidebar.slider("log10(K) (fixed)", -8.0, -2.0, -3.0, 0.1, key="fx_K")
        return 10.0 ** log_K
    if axis == "Omega_c":
        return st.sidebar.slider("Omega_c (fixed, rad/s)", 0.0, 45.0, 1.0, 0.1, key="fx_om")
    raise ValueError(axis)


st.sidebar.header("Modes (same as Tab 1)")
rotation_mode = st.sidebar.radio("rotation", ["none", "rigid", "differential"], horizontal=True)
A_over_Req = 0.0
if rotation_mode == "differential":
    A_over_Req = st.sidebar.slider(
        "A / R_eq (fixed law -- Omega_c is the only rotation axis available)",
        0.05, 3.0, 0.3, 0.05)

field_mode = st.sidebar.radio("field", ["none", "poloidal", "toroidal (self-consistent)"])
m_tor_sc = 1.0
lmax = 16
if field_mode == "poloidal":
    lmax = st.sidebar.slider("l_max", 4, 32, 16)
elif field_mode == "toroidal (self-consistent)":
    m_tor_sc = st.sidebar.slider("m (toroidal power law)", 1.0, 3.0, 1.0, 0.5)

st.sidebar.header("Sweep axes")
available_axes = ["rho_c"]
if field_mode == "poloidal":
    available_axes.append("k0")
elif field_mode == "toroidal (self-consistent)":
    available_axes.append("K")
if rotation_mode != "none":
    available_axes.append("Omega_c")

default_axes = available_axes[:2]
axes = st.sidebar.multiselect(
    "axes to sweep (pick 1 or 2)", options=available_axes,
    default=default_axes,
    help="A is never an axis (fixed law for differential rotation, see above). "
         "Only axes compatible with the current mode selection are offered."
)
if len(axes) == 0:
    st.sidebar.error("pick at least one axis")
    st.stop()
if len(axes) > 2:
    st.sidebar.error("pick at most two axes -- more explodes the grid")
    st.stop()

axis_grids = {}
fixed_values = {}
for axis in axes:
    st.sidebar.markdown(f"**{AXIS_LABELS[axis]} (swept)**")
    axis_grids[axis] = _axis_grid_ui(axis)
for axis in available_axes:
    if axis not in axes:
        st.sidebar.markdown(f"**{AXIS_LABELS[axis]} (fixed)**")
        fixed_values[axis] = _fixed_value_ui(axis)
# axes never offered because the relevant mode is off -- fixed at "no effect"
fixed_values.setdefault("k0", 0.0)
fixed_values.setdefault("K", 0.0)
fixed_values.setdefault("Omega_c", 0.0)
fixed_values.setdefault("rho_c", 1e12)

mu_e = st.sidebar.number_input("mu_e", min_value=1.0, max_value=2.5, value=2.0, step=0.1)
Nr = st.sidebar.select_slider("Nr", options=[65, 129, 161], value=129)
Ntheta = st.sidebar.select_slider("Ntheta", options=[33, 65, 129], value=65)
n_workers = st.sidebar.slider("parallel processes", 1, 8, 4)

rho_c_neutronization = eos.neutronization_threshold_rho_c(mu_e)
st.caption(
    f"Points with rho_c >= {rho_c_neutronization:.2e} g/cm³ (mu_e={mu_e}) are above the "
    "inverse beta-decay threshold (Boshkayev et al. 2013) -- shown, not filtered, "
    "marked distinctly on every plot below (see docs/teoria.md §1.12)."
)


def _grid_for(axis):
    return axis_grids[axis] if axis in axis_grids else [fixed_values[axis]]


rho_c_grid = _grid_for("rho_c")
k0_grid = _grid_for("k0")
K_grid = _grid_for("K")
Omega_c_grid = _grid_for("Omega_c")

total_points = len(rho_c_grid) * len(k0_grid) * len(K_grid) * len(Omega_c_grid)
st.write(f"Grid: {len(rho_c_grid)} (rho_c) x {len(k0_grid)} (k0) x {len(K_grid)} (K) "
         f"x {len(Omega_c_grid)} (Omega_c) = {total_points} points")

if st.button("run sweep", type="primary"):
    param_list = []
    for rho_c in rho_c_grid:
        for k0 in k0_grid:
            for K_tor in K_grid:
                for Omega_c in Omega_c_grid:
                    param_list.append({
                        "rho_c": float(rho_c), "k0": float(k0), "K_tor": float(K_tor),
                        "m_tor_sc": m_tor_sc, "Omega_c": float(Omega_c),
                        "A_over_Req": A_over_Req if rotation_mode == "differential" else 0.0,
                        "field_mode": field_mode, "rotation_mode": rotation_mode,
                        "mu_e": mu_e, "R_guess": seed.r_guess(rho_c), "Nr": Nr, "Ntheta": Ntheta,
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

    st.write(f"{len(cached_results)} already cached, {len(to_run)} to compute")

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
                    res = {"converged": False, "rho_c": p["rho_c"], "k0": p["k0"],
                           "K_tor": p["K_tor"], "Omega_c": p["Omega_c"], "error": str(e)}

                n_done += 1
                if res["converged"]:
                    h = store.save_run(cache_params_key(p), res["scalars"], res["fields"])
                    rows.append({"hash": h, **cache_params_key(p), **res["scalars"]})
                    n_ok += 1
                else:
                    fail_rows.append({"rho_c": p["rho_c"], "k0": p["k0"],
                                      "K_tor": p["K_tor"], "Omega_c": p["Omega_c"]})
                    n_fail += 1

                progress.progress(n_done / total_points)
                status.text(f"{n_done}/{total_points}  (converged: {n_ok}, failed: {n_fail})")

    st.session_state["sweep_rows"] = rows
    st.session_state["sweep_fail_rows"] = fail_rows
    st.success(f"sweep complete: {n_ok}/{total_points} converged "
               f"({n_ok / total_points:.0%})")

rows = st.session_state.get("sweep_rows", [])
fail_rows = st.session_state.get("sweep_fail_rows", [])

if not rows:
    st.info("Run the sweep to see the plots.")
    st.stop()

df = pd.DataFrame(rows)
if "rho_c_valid" not in df.columns:
    # cached rows from a pre-item-3 schema version would lack this column,
    # but SCHEMA_VERSION was bumped precisely to force a cache miss in that
    # case -- this branch is a defensive fallback, not the expected path
    df["rho_c_valid"] = True

if fail_rows:
    with st.expander(f"{len(fail_rows)} point(s) that did not converge"):
        st.dataframe(pd.DataFrame(fail_rows))

# ---------------- priority plot: M vs rho_c, one curve per K -------------
# This is the figure that goes to the collaborator (see closeout prompt,
# item 4). Built first, before any other plot, on purpose.
st.subheader("M vs rho_c, one curve per K (priority plot)")
# K is internal to this project's B_phi=K*rho^m*varpi^(2m-1) ansatz and
# means nothing outside it -- E_tor/|W| is the portable form (see
# docs/teoria.md Sec 6.2e: E_tor/|W| varies <1% across 2.5 decades in
# rho_c at fixed K, so it's the quantity that actually transfers to a
# collaborator, not K). No inverse solve exists yet for this branch
# (unlike the D6 poloidal branch's solve_zeta_for_energy_ratio()), so
# this can't be a real INPUT sweep axis -- E_mag/|W| (==E_tor/|W| here,
# since this branch has no poloidal component) is already computed per
# point, so grouping/labeling the already-swept K points by it is exact
# and does not need that missing piece.
_group_by_ve = (field_mode == "toroidal (self-consistent)" and "K_tor" in df.columns
                 and st.checkbox("group curves by E_tor/|W| instead of raw K", value=False))
if "K_tor" in df.columns and df["rho_c"].nunique() > 1:
    fig_priority = go.Figure()
    if _group_by_ve:
        df["_group_key"] = df["E_mag/|W|"].round(2)
        group_col = "_group_key"
    else:
        group_col = "K_tor"

    # K by itself doesn't say anything about the field -- neither does
    # E_tor/|W| alone, since that's a fraction of |W|, not a field
    # strength in gauss. B_tor,max is what answers "how strong a field,
    # in what units", so it has to be readable directly off this figure,
    # not just in a separate scalars table. Markers are colored by
    # B_tor,max (continuous); the K=0 (field-free) curve is excluded from
    # that color scale's range -- B_tor,max=0 there exactly, and including
    # it would compress the color range used by every field-carrying
    # curve for no reason.
    _has_b = "B_tor,max (G)" in df.columns
    _nonzero_b = df.loc[df["K_tor"] > 0, "B_tor,max (G)"] if _has_b else pd.Series(dtype=float)
    _b_cmin, _b_cmax = (float(_nonzero_b.min()), float(_nonzero_b.max())) if len(_nonzero_b) else (0.0, 1.0)
    _line_palette = px.colors.qualitative.Dark24
    _shown_colorbar = False

    for gi, (group_val, sub) in enumerate(df.sort_values("rho_c").groupby(group_col)):
        sub = sub.sort_values("rho_c")
        sub_valid = sub[sub["rho_c_valid"]]
        sub_invalid = sub[~sub["rho_c_valid"]]
        K_val = sub["K_tor"].iloc[0]
        e_tor_w_mean = sub["E_mag/|W|"].mean()
        is_zero = K_val <= 0

        b_range_str = ""
        if _has_b and not is_zero:
            b_lo, b_hi = sub["B_tor,max (G)"].min(), sub["B_tor,max (G)"].max()
            b_range_str = (f", B_tor,max≈{b_lo:.1e} G" if b_lo == b_hi
                            else f", B_tor,max: {b_lo:.1e}–{b_hi:.1e} G")

        if _group_by_ve:
            label = (f"E_tor/|W|≈{group_val:.2f}{b_range_str}" if group_val > 0
                      else "no toroidal field (K=0)")
        else:
            label = (f"K={K_val:.2e} (E_tor/|W|≈{e_tor_w_mean:.3f}{b_range_str})" if K_val > 0
                      else "K=0 (no toroidal field)")

        line_color = "black" if is_zero else _line_palette[gi % len(_line_palette)]
        if is_zero or not _has_b:
            marker = dict(size=8, color=line_color)
        else:
            marker = dict(size=9, color=sub["B_tor,max (G)"], colorscale="Plasma",
                           cmin=_b_cmin, cmax=_b_cmax, showscale=not _shown_colorbar,
                           colorbar=(dict(title="B_tor,max (G)", tickformat=".0e",
                                          exponentformat="power")
                                     if not _shown_colorbar else None),
                           line=dict(width=1.2, color=line_color))
            _shown_colorbar = True

        fig_priority.add_trace(go.Scatter(
            x=sub["rho_c"], y=sub["M/M_sun"], mode="lines+markers",
            name=label, line=dict(width=2, color=line_color), marker=marker))
        if len(sub_invalid) > 0:
            fig_priority.add_trace(go.Scatter(
                x=sub_invalid["rho_c"], y=sub_invalid["M/M_sun"], mode="markers",
                marker=dict(symbol="x", size=10, color="red"),
                name=f"{label} (rho_c >= neutronization threshold)",
                showlegend=bool(group_val == df[group_col].iloc[0])))

        # The central question this whole branch exists to answer: what
        # order of B_tor, at what rho_c, gives 2 Msun. Answering it
        # requires reading two numbers off two different curves at two
        # different rho_c today -- mark the M=2 crossing directly, with
        # both numbers attached, so it doesn't.
        if _has_b and len(sub) >= 2:
            rho_arr = sub["rho_c"].to_numpy()
            M_arr = sub["M/M_sun"].to_numpy()
            B_arr = sub["B_tor,max (G)"].to_numpy()
            sign_change = np.where(np.diff(np.sign(M_arr - 2.0)) != 0)[0]
            for ci in sign_change:
                m0, m1 = M_arr[ci], M_arr[ci + 1]
                if m1 == m0:
                    continue
                frac = (2.0 - m0) / (m1 - m0)
                logr0, logr1 = np.log10(rho_arr[ci]), np.log10(rho_arr[ci + 1])
                rho_cross = 10 ** (logr0 + frac * (logr1 - logr0))
                b0, b1 = B_arr[ci], B_arr[ci + 1]
                if b0 > 0 and b1 > 0:
                    B_cross = 10 ** (np.log10(b0) + frac * (np.log10(b1) - np.log10(b0)))
                else:
                    B_cross = b0 + frac * (b1 - b0)
                fig_priority.add_trace(go.Scatter(
                    x=[rho_cross], y=[2.0], mode="markers", showlegend=False,
                    marker=dict(symbol="star", size=15, color=line_color,
                                line=dict(width=1.5, color="black")),
                    hovertemplate=(f"M=2 Msun here<br>rho_c≈{rho_cross:.3e} g/cm³"
                                    f"<br>B_tor,max≈{B_cross:.3e} G<extra></extra>")))
                fig_priority.add_annotation(
                    x=rho_cross, y=2.0, xref="x", yref="y",
                    text=f"ρc≈{rho_cross:.2e}<br>B_tor,max≈{B_cross:.2e} G",
                    showarrow=True, arrowhead=2, arrowcolor=line_color,
                    ax=30, ay=-35 - 25 * (gi % 3), font=dict(size=10, color=line_color),
                    bgcolor="white", opacity=0.9)

    fig_priority.add_vline(x=rho_c_neutronization, line_dash="dash", line_color="red",
                            annotation_text="inverse beta-decay threshold "
                                            f"(mu_e={mu_e}, Boshkayev+2013)")
    fig_priority.add_hline(y=2.0, line_dash="dot", line_color="gray", line_width=1,
                            annotation_text="M = 2 M_sun", annotation_position="bottom right")
    _log_rho_axis(fig_priority, "x")
    _msun_axis(fig_priority, "y")
    st.plotly_chart(fig_priority, key="priority_chart")
    st.caption(
        "Red x markers / dashed vertical line: points at or past the inverse "
        "beta-decay (neutronization) threshold for this mu_e -- shown for "
        "context, not physically trustworthy past that line (item 3, "
        "docs/teoria.md §1.12). Not filtered out of the sweep. Marker color "
        "is B_tor,max (G) at that point -- not constant along a curve, "
        "since the same K gives a different field strength at different "
        "rho_c (docs/teoria.md §6.2e). Star markers + annotations: where "
        "each curve crosses M = 2 M_sun, linearly interpolated between the "
        "two bracketing grid points (in log10(rho_c) for rho_c, in "
        "log10(B) for B_tor,max)."
    )
    _pub_pdf = _publication_figure_pdf(df, rho_c_neutronization, mu_e)
    st.download_button(
        "export figure for publication (PDF)", data=_pub_pdf,
        file_name="M_vs_rhoc_toroidal.pdf", mime="application/pdf",
        help="Vector PDF, matplotlib (not this interactive Plotly chart): "
             "log-log, one curve per E_tor/|W| (portable form, not raw K), "
             "K=0 curve highlighted, labeled neutronization threshold, no "
             "embedded title -- the figure for the collaborator/paper."
    )
else:
    st.caption(
        "Needs rho_c swept with more than one value to be meaningful "
        "(K is shown per-curve even when only K=0 was run)."
    )

st.subheader("M-R diagram")
_color_col = "k0" if field_mode == "poloidal" else ("K_tor" if field_mode == "toroidal (self-consistent)" else "Omega_c")
df["validity"] = df["rho_c_valid"].map({True: "valid", False: "past neutronization threshold"})
fig_mr = px.scatter(df, x="R_eq (km)", y="M/M_sun", color=_color_col,
                     symbol="validity",
                     color_continuous_scale="Viridis",
                     hover_data={"hash": True, "rho_c": ":.3e", "VE": ":.3e",
                                 "B_pol,max (G)": ":.3e"})
fig_mr.update_xaxes(title="R_eq (km)")
_msun_axis(fig_mr, "y")
fig_mr.add_hline(y=1.44, line_dash="dash", line_color="gray",
                  annotation_text="Chandrasekhar limit (mu_e=2, no field)")
ref_file = REF_DIR / "bera_bhattacharya_2014.csv"
if ref_file.exists():
    ref_df = pd.read_csv(ref_file)
    fig_mr.add_trace(go.Scatter(x=ref_df["R_km"], y=ref_df["M_Msun"], mode="markers",
                                 marker=dict(symbol="cross", color="red", size=8),
                                 name="Bera & Bhattacharya 2014"))
else:
    st.caption(f"(no literature overlay — file not found: {ref_file})")
event_mr = st.plotly_chart(fig_mr, on_select="rerun", key="mr_chart")

st.subheader("M vs rho_c, colored by E_mag/|W|")
fig_m_rhoc = px.scatter(df, x="rho_c", y="M/M_sun", color="E_mag/|W|",
                         color_continuous_scale="Inferno", log_x=True,
                         hover_data={"hash": True, "k0": ":.3e", "VE": ":.3e",
                                     "B_pol,max (G)": ":.3e"})
_log_rho_axis(fig_m_rhoc, "x")
_msun_axis(fig_m_rhoc, "y")
st.plotly_chart(fig_m_rhoc, key="m_rhoc_chart")

st.subheader("M vs rho_c, colored by B_pol,max (G)")
fig_m_b = px.scatter(df, x="rho_c", y="M/M_sun", color="B_pol,max (G)",
                      color_continuous_scale="Plasma", log_x=True,
                      hover_data={"hash": True, "k0": ":.3e", "VE": ":.3e"})
fig_m_b.update_layout(coloraxis_colorbar=dict(tickformat=".2e", title="B (G)"))
_log_rho_axis(fig_m_b, "x")
_msun_axis(fig_m_b, "y")
st.plotly_chart(fig_m_b, key="m_b_chart")

if field_mode == "poloidal":
    st.subheader("Dipolarity (B_polo/B_eq) vs k0")
    st.caption(
        "B_polo/B_eq at the stellar surface — exactly 2 for a pure dipole; "
        "departure measures multipole contamination as the poloidal field "
        "strengthens (see diagnostics.surface_dipolarity, Tab 1). k0=0 points "
        "are omitted (log x-axis, and the field-free case has no field to "
        "measure a ratio of)."
    )
    df_dip = df[df["k0"] > 0]
    if len(df_dip) > 0 and "dipolarity" in df_dip.columns:
        fig_dip = px.scatter(df_dip, x="k0", y="dipolarity", color="rho_c",
                              color_continuous_scale="Viridis", log_x=True,
                              hover_data={"hash": True, "rho_c": ":.3e", "VE": ":.3e",
                                          "B_polo (G)": ":.3e", "B_eq (G)": ":.3e"})
        fig_dip.update_layout(coloraxis_colorbar=dict(tickformat=".0e", title="rho_c (g/cm³)",
                                                        exponentformat="power"))
        fig_dip.add_hline(y=2.0, line_dash="dash", line_color="gray",
                           annotation_text="pure dipole (=2)")
        st.plotly_chart(fig_dip, key="dipolarity_chart")
    else:
        st.caption("no k0>0 points with dipolarity data yet — run the sweep.")

if rotation_mode != "none":
    st.subheader("T/|W| and mass-loss ratio vs Omega_c")
    fig_rot = go.Figure()
    fig_rot.add_trace(go.Scatter(x=df["Omega_c (rad/s)"], y=df["T/|W|"], mode="markers",
                                  name="T/|W|", marker=dict(color="royalblue")))
    fig_rot.add_trace(go.Scatter(x=df["Omega_c (rad/s)"], y=df["equatorial mass-loss ratio"],
                                  mode="markers", name="mass-loss ratio",
                                  marker=dict(color="firebrick"), yaxis="y2"))
    fig_rot.update_layout(
        xaxis_title="Omega_c (rad/s)", yaxis_title="T/|W|",
        yaxis2=dict(title="mass-loss ratio", overlaying="y", side="right"))
    fig_rot.add_hline(y=0.14, line_dash="dash", line_color="royalblue",
                       annotation_text="secular instability (0.14)")
    st.plotly_chart(fig_rot, key="rot_chart")

st.subheader("VE heat map over the grid")
_axis_column = {"rho_c": "rho_c", "k0": "k0", "K": "K_tor", "Omega_c": "Omega_c"}
if len(axes) == 2:
    ax1, ax2 = axes
    col1, col2 = _axis_column[ax1], _axis_column[ax2]
    if df[col1].nunique() > 1 and df[col2].nunique() > 1:
        pivot = df.pivot_table(index=col2, columns=col1, values="VE", aggfunc="mean")
        fig_ve = go.Figure(data=go.Heatmap(
            z=np.log10(pivot.values), x=pivot.columns, y=pivot.index,
            colorscale="RdBu_r", colorbar=dict(title="log10(VE)")))
        fig_ve.update_xaxes(title=AXIS_LABELS[ax1])
        fig_ve.update_yaxes(title=AXIS_LABELS[ax2])
        # this heatmap's x/y are the raw pivot coordinate values, plotted
        # on a LINEAR axis by default -- for rho_c (spanning decades)
        # that's the same currency/data-size-looking mislabeling as the
        # scatter plots above (e.g. gridlines at "0.2T/0.4T/0.6T"), not
        # just a tick-format issue: type="log" was never set here at all.
        if ax1 == "rho_c":
            _log_rho_axis(fig_ve, "x")
        if ax2 == "rho_c":
            _log_rho_axis(fig_ve, "y")
        st.plotly_chart(fig_ve, key="ve_heatmap")
    else:
        st.caption("grid too small for a 2D heat map")
else:
    st.caption("only one axis swept — no 2D heat map to show")

st.subheader("Load a point in Tab 1")
selected_hash = st.selectbox("or pick by hash", df["hash"].tolist())
if event_mr and event_mr.get("selection", {}).get("points"):
    pt = event_mr["selection"]["points"][0]
    idx = pt.get("point_index")
    if idx is not None and idx < len(df):
        selected_hash = df.iloc[idx]["hash"]
        st.write(f"selected on the plot: {selected_hash}")

if st.button("load selected equilibrium in Tab 1"):
    run = store.load_run(selected_hash)
    p = run["params"]
    st.session_state["reload_run_params"] = {
        "rho_c": p["rho_c"], "k0": p.get("k0", 0.0),
        "mu_e": p.get("mu_e", 2.0), "zeta_target_ratio": 0.0, "m_tor": 1,
        "K_tor": p.get("K_tor", 0.0), "m_tor_sc": p.get("m_tor_sc", 1.0),
        "Omega_c": p.get("Omega_c", 0.0), "A_over_Req": p.get("A_over_Req", 0.0),
        "Nr": p.get("Nr", 129), "Ntheta": p.get("Ntheta", 129),
        "lmax": p.get("lmax", 16), "tol": p.get("tol", 1e-6),
        "max_iter": p.get("max_iter", 200),
    }
    st.switch_page("pages/1_equilibrium.py")
