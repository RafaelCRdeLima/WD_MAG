"""Guards against docs/teoria.md and docs/teoria.tex drifting apart.
teoria.tex is hand-derived from teoria.md, not auto-generated (see the
sync note at the top of both files) -- this project was burned by exactly
this drift once: teoria.tex/teoria.pdf carried the wrong sign for the
magnetic virial identity for a full commit while teoria.md already had
the correction, and PDF is the format that circulates to the
collaborator. A PDF with an authoritative look but a wrong sign is worse
than no PDF.

Compares git commit timestamps: if teoria.md was committed more recently
than teoria.tex, teoria.tex is presumed stale.

LIMITATION: this only catches the drift once both files are committed
(git log reads committed history, not the working tree) -- it is not a
pre-commit hook, just a regression test that fails the next time the test
suite runs after a bad commit. Good enough given this project has no CI;
run `pytest` before considering a doc change finished.
"""

import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs"


def _last_commit_time(path):
    out = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "log", "-1", "--format=%ct", "--", str(path)],
        capture_output=True, text=True, timeout=5,
    )
    if out.returncode != 0 or not out.stdout.strip():
        return None
    return int(out.stdout.strip())


def test_teoria_tex_not_older_than_teoria_md():
    md_time = _last_commit_time(DOCS_DIR / "teoria.md")
    tex_time = _last_commit_time(DOCS_DIR / "teoria.tex")
    if md_time is None or tex_time is None:
        pytest.skip("teoria.md/teoria.tex have no commit history yet")
    assert tex_time >= md_time, (
        f"docs/teoria.md was committed more recently (unix {md_time}) than "
        f"docs/teoria.tex (unix {tex_time}) -- teoria.tex is hand-derived "
        "from teoria.md and must be updated (content AND a recompiled "
        "teoria.pdf) in the same commit as any teoria.md content change, "
        "or in a later one. See the sync note at the top of both files."
    )


if __name__ == "__main__":
    test_teoria_tex_not_older_than_teoria_md()
    print("OK")
