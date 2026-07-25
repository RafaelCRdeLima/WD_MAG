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
