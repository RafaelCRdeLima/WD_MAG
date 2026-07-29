"""Extraction wrappers -- call the existing command-line tools
(`finterior`), never reimplement the physics. Every value returned here
must match what running the tool directly on the same plotfile gives;
if it doesn't, that's the "second calculation path diverging" bug class
this project has already paid for more than once -- see
tests/test_extraction.py, which checks exactly this.
"""

import re
import subprocess
from pathlib import Path

FINTERIOR = Path(
    "/home/rafael/wd-magnetizada/castro/external/amrex/Tools/Plotfile/finterior.gnu.ex"
)

_TIME_RE = re.compile(r"time = ([\d.eE+-]+)")
_INTERIOR_RE = re.compile(
    r"INTERIOR.*?min = ([\d.eE+-]+)\s+max = ([\d.eE+-]+)(?:\s+volume_sum = ([\d.eE+-]+))?"
)


def run_finterior(var: str, plotfile: Path, margin: int = 2, finterior_path: Path = FINTERIOR) -> dict:
    """Thin wrapper over the finterior CLI tool. Returns exactly what the
    tool reports (time, interior min/max/volume_sum) -- no recomputation.
    """
    result = subprocess.run(
        [str(finterior_path), "-v", var, "-m", str(margin), str(plotfile)],
        capture_output=True, text=True, timeout=60,
    )
    out = result.stdout
    t_match = _TIME_RE.search(out)
    interior_match = _INTERIOR_RE.search(out)
    if t_match is None or interior_match is None:
        raise RuntimeError(f"finterior output not parseable for {var} on {plotfile}:\n{out}\n{result.stderr}")
    return {
        "time": float(t_match.group(1)),
        "min": float(interior_match.group(1)),
        "max": float(interior_match.group(2)),
        "volume_sum": float(interior_match.group(3)) if interior_match.group(3) else None,
    }


def extract_full_series(run_dir: Path, run_id: str, W_abs_erg: float, t_dyn_s: float) -> list[dict]:
    """Every plotfile for a run, in time order -- for the per-seed
    inspection view (reads live plotfiles, not the persisted summary
    row; the summary is one point, this is the full trajectory needed
    for the four time-series plots).
    """
    plotfiles = sorted(Path(run_dir).glob(f"plt_{run_id}0*"))
    series = []
    for plotfile in plotfiles:
        try:
            series.append(extract_field_measurement(plotfile, W_abs_erg, t_dyn_s))
        except RuntimeError:
            continue  # partially-written plotfile (e.g. still being flushed) -- skip
    return series


def extract_field_measurement(plotfile: Path, W_abs_erg: float, t_dyn_s: float) -> dict:
    """E_mag/|W| and E_tor/E_mag at one plotfile, via finterior's
    interior-masked volume_sum (never the raw/edge value -- the
    ca_derdivb ghost-cell bug this session found applies to every
    derived variable computed the same way, not just Div_B).
    """
    emag = run_finterior("emag_density", plotfile)
    etor = run_finterior("etor_density", plotfile)
    divb = run_finterior("Div_B", plotfile)

    t_ttdyn = emag["time"] / t_dyn_s
    e_mag_over_w = emag["volume_sum"] / W_abs_erg if emag["volume_sum"] is not None else float("nan")
    e_tor_over_emag = (
        etor["volume_sum"] / emag["volume_sum"]
        if emag["volume_sum"] else float("nan")
    )

    return {
        "t_ttdyn": t_ttdyn,
        "E_mag_over_W": e_mag_over_w,
        "E_tor_over_Emag": e_tor_over_emag,
        "divB_interior_max": max(abs(divb["min"]), abs(divb["max"])),
    }
