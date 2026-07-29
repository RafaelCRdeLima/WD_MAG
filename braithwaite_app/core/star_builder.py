"""The three functions behind Passo 1.5 (background star construction).
Each is standalone and tested against real numbers from this session
before any UI wiring -- see tests/test_star_builder.py.

This module does NOT reimplement any physics -- it only reproduces
calculations that were done by hand (ad hoc Python one-liners, or manual
grep/run checks) during this session, as real, tested functions.
"""

import re
import subprocess
from pathlib import Path


class GravityPatchMissing(RuntimeError):
    """Raised by verify_gravity_patch when the r=0 division-by-zero patch
    (Source/gravity/Gravity.cpp) is not present in source or not proven
    present in the built executable's behavior. Never caught silently --
    a half-shift launch without this patch crashes with a division by
    zero at the exact cell where the star's density peak sits.
    """


# ---------------------------------------------------------------------
# (a) half-shift domain geometry
# ---------------------------------------------------------------------
def half_shift_domain(n_cell: int, box_half_width_cm: float) -> tuple[float, float]:
    """Domain bounds shifted by half a cell so a cell CENTER (not a
    vertex) lands exactly on the star's center at the origin -- the fix
    that reduced the peak-density sampling deficit from -4.78% to
    -1.16% at 64^3 (see docs/teoria.md Sec 6.7/6.8).

    dx = 2*box_half_width_cm / n_cell
    prob_lo = -(n_cell+1)/2 * dx
    prob_hi = prob_lo + n_cell*dx

    Requires n_cell even -- AMReX/Castro hard-require this at the base
    level ("defBaseLevel: must have even number of cells", a real error
    hit this session, not a guess). Odd n_cell is not a valid alternative
    to reach the same cell-centered geometry.
    """
    if n_cell % 2 != 0:
        raise ValueError(
            f"n_cell={n_cell} is odd -- Castro requires an even n_cell at the "
            "base level (confirmed via a real 'defBaseLevel: must have even "
            "number of cells' error, not assumed). Odd n_cell is not a valid "
            "way to center a cell on the origin."
        )
    dx = 2.0 * box_half_width_cm / n_cell
    prob_lo = -(n_cell + 1) / 2.0 * dx
    prob_hi = prob_lo + n_cell * dx
    return prob_lo, prob_hi


# ---------------------------------------------------------------------
# (b) measurement window = [t_field_relax, X], never [0, X]
# ---------------------------------------------------------------------
def _last_leaving_crossing(series: list[tuple[float, float]], threshold: float) -> float:
    """series: sorted [(t_ttdyn, pct_deviation), ...]. Returns the t/t_dyn
    of the FIRST sample outside +-threshold (the "leaving band" row
    itself) -- reproduces exactly the ad hoc crossing-point script used
    repeatedly this session (e.g. the 1%/2% boundaries for seed 42:
    0.573/1.128 t_dyn, verified below against real log data; that script
    printed the crossing row's own t/t_dyn, not the prior in-band sample).
    """
    was_in = True
    for t, pct in series:
        is_in = abs(pct) <= threshold
        if was_in and not is_in:
            return t
        was_in = is_in
    return series[-1][0]  # never left the band across the whole series


def find_measurement_window(
    rho_c_series: list[tuple[float, float]],
    rho_c_ic: float,
    t_field_relax_ttdyn: float = 0.4,
    rho_c_thresholds: tuple[float, float] = (0.01, 0.02),
) -> dict:
    """rho_c_series: [(t_ttdyn, rho_c), ...], sorted by t. rho_c_ic: the
    true t=0 value (measured, not the requested target).

    t_field_relax_ttdyn default (0.4) is the value validated THIS SESSION
    specifically for E_mag/|W|=0.15 (E_tor/E_mag reached its quasi-steady
    band by t/t_dyn~0.4-1.1 for that field strength) -- not a universal
    constant. Callers studying a different field energy should pass a
    recalibrated value once that timescale has actually been measured for
    it; this function does not silently assume 0.4 generalizes.

    Returns:
      {"valid": False, "reason": str}  -- star drifts before the field
        would have relaxed; there is no honest measurement window.
      {"valid": True, "window": (t_field_relax_ttdyn, X), "X_1pct": ...,
        "X_2pct": ...}
    """
    if not rho_c_series:
        raise ValueError("empty rho_c_series")

    series_sorted = sorted(rho_c_series, key=lambda r: r[0])
    pct_series = [(t, 100.0 * (rho - rho_c_ic) / rho_c_ic) for t, rho in series_sorted]

    t1, t2 = sorted(rho_c_thresholds)
    x_1pct = _last_leaving_crossing(pct_series, 100.0 * t1)
    x_2pct = _last_leaving_crossing(pct_series, 100.0 * t2)

    # X is the wider (2%) boundary -- the operational validity limit used
    # throughout this session's later analysis.
    X = x_2pct

    if X <= t_field_relax_ttdyn:
        return {
            "valid": False,
            "reason": (
                f"estrela deriva antes do campo relaxar (X={X:.3f} t/t_dyn, "
                f"t_field_relax={t_field_relax_ttdyn:.3f} t/t_dyn) -- sem "
                "janela de medição válida"
            ),
            "X_1pct": x_1pct,
            "X_2pct": x_2pct,
        }

    return {
        "valid": True,
        "window": (t_field_relax_ttdyn, X),
        "X_1pct": x_1pct,
        "X_2pct": x_2pct,
    }


def parse_rho_c_log(log_path: Path, t_dyn_s: float) -> list[tuple[float, float]]:
    """Reads Castro's own 'TIME= ... MAXIMUM DENSITY = ...' lines (the raw,
    whole-domain max -- immune to the Div_B-style ghost-cell bug, so no
    interior masking needed for density) and returns [(t/t_dyn, rho_c), ...].
    """
    rows = []
    with open(log_path) as f:
        for line in f:
            if line.startswith("TIME=") and "MAXIMUM DENSITY" in line:
                parts = line.split()
                t = float(parts[1])
                rho = float(parts[-1])
                rows.append((t / t_dyn_s, rho))
    return rows


# ---------------------------------------------------------------------
# (c) gravity r=0 patch verification -- mandatory, two layers
# ---------------------------------------------------------------------
_PATCH_SIGNATURE = "if (r > 0.0_rt)"


def _check_patch_static(castro_source_dir: Path) -> None:
    gravity_cpp = Path(castro_source_dir) / "Source" / "gravity" / "Gravity.cpp"
    if not gravity_cpp.exists():
        raise GravityPatchMissing(f"{gravity_cpp} not found -- is castro_source_dir correct?")
    text = gravity_cpp.read_text()
    if _PATCH_SIGNATURE not in text:
        raise GravityPatchMissing(
            f"Gravity r=0 patch not found in {gravity_cpp} (missing guard "
            f"'{_PATCH_SIGNATURE}' in interpolate_monopole_grav). The Castro "
            "checkout was reset/recloned without reapplying the patch "
            "documented in castro_problems/wd_braithwaite/CASTRO_CORE_PATCHES.md. "
            "Reapply it before using half-shift geometry -- otherwise the "
            "build will crash with a division by zero at the star's center."
        )


def verify_gravity_patch(
    castro_source_dir: Path,
    executable_path: Path,
    finterior_path: Path,
    fline_path: Path,
    static_ic_inputs: Path,
    scratch_dir: Path,
) -> None:
    """Two layers, both mandatory:
      1. Static: grep Gravity.cpp for the patch guard (catches "checkout
         lost the patch"). Fast (milliseconds), no build/run needed.
      2. Functional: run the same static IC check (max_step=0, half-shift
         geometry) used to validate the patch this session, and confirm
         grav_x(r=0) == 0.0 exactly via `fline`. Catches "source has the
         patch but this binary is stale" -- text in a file does not prove
         what is compiled. This is the authoritative check; call sites
         should cache its result by executable hash and not repeat it
         per seed.

    Raises GravityPatchMissing with a clear message on either failure.
    Does not launch or gate anything itself beyond raising.
    """
    _check_patch_static(castro_source_dir)

    scratch_dir = Path(scratch_dir)
    scratch_dir.mkdir(parents=True, exist_ok=True)
    plotfile = scratch_dir / "plt_gravity_patch_check00000"
    subprocess.run(
        [
            "rm", "-rf", str(plotfile),
        ],
        check=False,
    )
    result = subprocess.run(
        [
            str(executable_path), str(static_ic_inputs),
            "max_step=0", "stop_time=0.0",
            f"amr.plot_file={scratch_dir / 'plt_gravity_patch_check'}",
            "amr.plot_int=1", "amr.check_int=-1",
            "amr.derive_plot_vars=density grav_x grav_y grav_z",
        ],
        cwd=str(Path(executable_path).parent),
        capture_output=True, text=True, timeout=120,
    )
    if result.returncode != 0 or not plotfile.exists():
        raise GravityPatchMissing(
            "Functional gravity-patch check failed to run cleanly (static IC "
            f"check aborted or produced no plotfile). stderr tail:\n"
            f"{result.stderr[-2000:]}"
        )

    fline_out = subprocess.run(
        [str(fline_path), "-v", "grav_x", "-j", "32", "-k", "32",
         "-ilo", "32", "-ihi", "32", str(plotfile)],
        capture_output=True, text=True, timeout=30,
    )
    m = re.search(r"grav_x=([\-0-9.eE+]+)", fline_out.stdout)
    if not m:
        raise GravityPatchMissing(
            f"Could not read grav_x from fline output: {fline_out.stdout!r} "
            f"{fline_out.stderr!r}"
        )
    grav_x_at_center = float(m.group(1))
    if grav_x_at_center != 0.0:
        raise GravityPatchMissing(
            f"Functional check failed: grav_x(r=0) = {grav_x_at_center!r}, "
            "expected exactly 0.0. The executable does not behave as if the "
            "r=0 gravity patch is compiled in, even though the source has it "
            "-- rebuild before using half-shift geometry."
        )
