"""Verifies the extraction wrapper against a direct, independent call to
the same CLI tool -- catches the "second calculation path diverging"
bug class this project has paid for before. Uses a real plotfile
already on disk from this session's runs.
"""

import re
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core.extraction import FINTERIOR, run_finterior

WD_BRAITHWAITE = Path("/home/rafael/wd-magnetizada/castro/Exec/science/wd_braithwaite")


def _find_a_real_plotfile():
    candidates = sorted(WD_BRAITHWAITE.glob("plt_*0000"))
    return candidates[0] if candidates else None


def test_run_finterior_matches_direct_cli_call():
    plotfile = _find_a_real_plotfile()
    if plotfile is None:
        pytest.skip("no plotfile on disk to test extraction against")

    wrapped = run_finterior("density", plotfile)

    raw = subprocess.run(
        [str(FINTERIOR), "-v", "density", "-m", "2", str(plotfile)],
        capture_output=True, text=True, timeout=60,
    ).stdout
    t = float(re.search(r"time = ([\d.eE+-]+)", raw).group(1))
    m = re.search(r"INTERIOR.*?min = ([\d.eE+-]+)\s+max = ([\d.eE+-]+)", raw)
    direct_min, direct_max = float(m.group(1)), float(m.group(2))

    assert wrapped["time"] == pytest.approx(t)
    assert wrapped["min"] == pytest.approx(direct_min)
    assert wrapped["max"] == pytest.approx(direct_max)
