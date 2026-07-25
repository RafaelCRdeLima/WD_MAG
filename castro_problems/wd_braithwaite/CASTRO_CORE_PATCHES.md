# Local patches to Castro's own `Source/` (not upstream, not covered by `scripts/sync_wd_braithwaite.sh`)

`scripts/sync_wd_braithwaite.sh` only mirrors files inside
`Exec/science/wd_braithwaite/`. The patch below lives in Castro's own
`Source/driver/`, so it is not captured by that script and would be
lost on a fresh `git clone`/`git submodule update` of the Castro
checkout. Re-apply manually if that ever happens.

## `Source/driver/Castro_io.cpp` — missing `#include <extern_parameters.H>`

**Symptom:** switching `EOS_DIR` to `ztwd` (see this problem's
`GNUmakefile`) breaks the build with `'network_rp' has not been
declared` inside the auto-generated `extern_job_info_tests.H`.

**Root cause:** `writeJobInfo()` in `Castro_io.cpp` includes
`extern_job_info_tests.H`, which references the `network_rp` namespace
declared in the auto-generated `extern_parameters.H`. Every EOS this
project had used until now (`gamma_law`) happens to `#include
<extern_parameters.H>` itself (for its own `eos_rp` parameters), so
this was always transitively available in `Castro_io.cpp`. `ztwd` has
zero runtime parameters (`EOS/ztwd/_parameters` is effectively empty)
and so never includes it, and nothing else in that translation unit
does either -- a real, EOS-module-specific gap in Castro itself, not
something in our problem code. No other Exec problem in this checkout
uses `ztwd`, which is presumably why this was never hit before.

**Fix applied** (insert after the existing includes near the top of
`Source/driver/Castro_io.cpp`, before line ~15):

```cpp
#include <extern_parameters.H>
```

See `docs/teoria.md` Sec 6.5 for the full context (the EOS mismatch
investigation this patch was needed for).

## `Div_B` derive registration -- ghost-cell gap at box boundaries

**Symptom:** `Div_B` (Castro's own built-in derived plotfile variable)
shows astronomical, physically nonsensical extrema (up to ~1e17 in one
run, ~1e76 in another under `-np 4`) that appear, grow, and sometimes
partially subside across a run -- looked exactly like a broken CT
(constrained transport) scheme, enough to kill a run that (per the
investigation below) was actually fine.

**Root cause, confirmed with `tools/finterior.cpp`** (mirrored here,
build with `cd external/amrex/Tools/Plotfile && make programs=finterior`,
not covered by `sync_wd_braithwaite.sh` since it lives outside
`Exec/science/wd_braithwaite/`): `derive_lst.add("Div_B", ...,
ca_derdivb, the_same_box)` in `Source/driver/Castro_setup.cpp` requests
no extra ghost cells for the derive computation, but `ca_derdivb`
(`Source/driver/Derive.cpp`) reads `dat(i+1,j,k)` / `(i,j+1,k)` /
`(i,j,k+1)` -- one cell beyond the box's valid region. At the true
domain boundary and at every internal box-to-box boundary (this
project's grids: `max_grid_size=32` on a 64^3 domain -> 8 boxes), that
read can land in unfilled ghost data, producing an extremum that is
essentially reading uninitialized/stale memory as a `Real`.

Confirmed directly: masking out the 2 cells nearest any box boundary
(`finterior.cpp`, margin=2) on the SAME plotfiles that showed ~1e17
gives max|Div_B| ~1e-10 (absolute, essentially float64 noise) in the
interior, at every single time sampled across a full run -- the actual
CT-evolved field is divergence-free to machine precision throughout;
only the box-boundary cells in the derived diagnostic are wrong. Not
fixed upstream (would require either growing the derive's requested box
by 1, `grow_box_by_one`, matching how "magvort" is registered a few
lines below it in the same file, or ensuring the source MultiFab's
ghost cells are filled before the derive runs) -- for now, always read
`Div_B` through `finterior.cpp` (or equivalent box-boundary masking),
never raw `fextrema`, when using it as a correctness check.

See `docs/teoria.md` Sec 6.6 for the full investigation (including how
this was distinguished from the global-damping and AMR-regrid
hypotheses that were tested and ruled out along the way).
