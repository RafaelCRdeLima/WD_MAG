# Patches to FLASH 4.8 needed to build this problem

`flash_crosscheck/FLASH4.8/` is gitignored (420 MB, extracted from the
`flash48:latest` Docker image at `/opt/FLASH4.8`), so anything changed
inside it has to be recorded here or it is lost. Same reasoning as
`castro_problems/wd_braithwaite/CASTRO_CORE_PATCHES.md`.

## 1. `use iso_c_binding` without `intrinsic` breaks the build

**Symptom.** Every compile of a file that uses `iso_c_binding` dies with

```
f951: Fatal Error: Reading module 'iso_c_binding.mod' at line 1 column 1: Unexpected EOF
```

on a *clean* object directory, so it is not stale artifacts.

**Cause.** FLASH's setup scans `use` statements to build the dependency
list, and treats `iso_c_binding` as a module it has to produce. The
generated Makefile then satisfies that prerequisite with

```
touch iso_c_binding.mod
```

which leaves a zero-byte `iso_c_binding.mod` in the object directory.
gfortran finds that file ahead of its own intrinsic module and fails
reading it. The scanner only makes this mistake for the bare form: where
the source already says `use, intrinsic :: iso_c_binding`, the dependency
is correctly skipped. `source/Cpp/Cpp_interface.F90` has both forms, on
lines 5 and 14, which is why the failure looks arbitrary.

**Patch.** Rewrite the bare form as the intrinsic form -- semantically
identical Fortran, and it removes the bogus dependency:

```
sed -i 's/^\(\s*\)use iso_c_binding\s*$/\1use, intrinsic :: iso_c_binding/' \
    source/Cpp/Cpp_interface.F90 \
    source/Cpp/CppMain/Cpp_strings.F90 \
    source/Cpp/Cpp_strings.F90 \
    source/RuntimeParameters/CppAPI/RuntimeParameters_cAPI.F90
```

Other files under `source/Simulation/SimulationMain/python/` and
`unitTest/` have the same bare form but are not pulled into this setup, so
they are left alone.

## 2. `sites/rafael` points at host paths that do not exist in the container

`sites/rafael/Makefile.h` sets the compiler to
`/home/rafael/micromamba/envs/flash_gcc12/bin/mpif90`. Inside the image
the environment is at `/opt/micromamba/...` and there is no `/home/rafael`,
so every compile fails with `Error 127` (command not found).

**Fix.** A separate site, so the existing one is left untouched for
whatever native build it belongs to:

```
cp -r sites/rafael sites/wd_docker
sed -i 's|/home/rafael/micromamba|/opt/micromamba|g' sites/wd_docker/Makefile.h
```

## 3. The container runs as root by default

Without `--user`, `setup` and `make` leave root-owned files in the
mounted tree that the host user then cannot delete. Always run with

```
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -v "$PWD/FLASH4.8:/opt/FLASH4.8" flash48:latest ...
```

## 4. Runtime, not build: two EOS settings this problem needs

Neither is a patch to FLASH, but both are non-obvious and the run aborts
without them.

- **`helm_table.dat` must be in the run directory.** The object directory
  gets a symlink to `source/physics/Eos/EosMain/Helmholtz/helm_table.dat`;
  copying only `flash4` to a separate run directory leaves it behind and
  `Eos_init` aborts.
- **`eos_coulombMult = 0.0`, `eos_coulombAbort = .false.`** The Coulomb
  correction drives the total pressure negative in the cold ambient and
  aborts the run. Disabling it also brings the evolved EOS *closer* to the
  ztwd model handed in, which has no Coulomb term.
- **The ambient cannot be a hard vacuum.** At `rho = 1e-6` the Helmholtz
  Newton-Raphson fails to converge (50 iterations, then abort). The
  ambient was raised to `1.0 g/cm^3` with `smlrho = 1e-2` -- nine orders
  below `rho_c`, so negligible in mass, but enough for the EOS inversion.

## 5. The half-shift geometry is blocked by the multipole centre snapping

**Not a patch -- an unsolved incompatibility, recorded so the next attempt
does not start from scratch.**

The half-shift geometry (docs/teoria.md Sec 6.6) is what makes the
measurement window exist: it takes the initial-condition deficit from
`-4.78%` to `-1.16%` and the initial transient from `+2.76%` to `+0.81%`,
and without it neither code yields a valid window at all. FLASH reproduces
the IC number (`-1.23%`), so the geometry itself ports. What does not port
is the gravity.

**Diagnosed, not guessed**, by the same check Castro's Sec 6.7 patch used
-- dump the potential along a line through the origin and look:

| | `gpot(centre)` | `gpot(neighbour)` | ratio | `g_x` at the neighbour |
|---|---|---|---|---|
| `Multipole` (legacy) | `-5.02e22` | `-2.41e18` | `2.1e4` | `-1.64e15` cm/s^2 |
| `Multipole_new` | `-2.38e18` | `-2.32e18` | `1.025` | `-5.36e9` cm/s^2 |

The legacy solver puts a four-order-of-magnitude spike in the potential at
`r = 0`, which is why `dt` was throttled `100x` with the limiter pinned at
the origin from step 1. **Switching to `Multipole_new` (setup shortcut
`+newMpole`) removes it** and restores `dt` to `9.5e-3`. That much is a
real fix and is kept.

What remains is subtler. With `Multipole_new` the potential is smooth but
symmetric about `j - 1/2`, half a cell to the *left* of the star, giving
`g_x(0) = -1.87e9` cm/s^2 where symmetry demands zero. The cause is in
`gr_mpoleCen3Dcartesian.F90:336`: the solver computes its expansion centre
and then **snaps it to the nearest cell face**. For a vertex-centred
(symmetric) domain that is a no-op, since the star's centre already sits
on a face. For the half-shift domain the star's centre sits on a cell
*centre*, so the snap moves the expansion centre `dx/2` away and the star
feels a spurious central acceleration.

So the two requirements are in direct conflict: good central-density
sampling wants the star's centre at a cell centre, and FLASH's multipole
wants it at a cell face.

**The obvious patch does not work.** Redirecting the snap to cell centres
(`Grid_getCellCoords(..., CENTER, ...)` instead of `FACES`) builds cleanly
and then segfaults at initialization -- with the expansion centre on a
cell centre, at least one cell sits at radius exactly zero, which the
inner-zone radial binning does not survive. That is the same class of
problem Castro solved with a `r == 0` guard, but here it lives inside the
solver's binning machinery rather than in a single acceleration formula.
The patch was reverted; the tree builds and runs.

Making FLASH correct for this problem therefore means adapting
`Multipole_new`'s inner-zone binning to tolerate a zero-radius cell, which
is a real piece of work inside the solver, not a two-line change.

## Build and run, end to end

```
scf/.venv/bin/python3 flash_crosscheck/make_wd_model.py
cp -r flash_crosscheck/WDHydrostatic FLASH4.8/source/Simulation/SimulationMain/

docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -v "$PWD/FLASH4.8:/opt/FLASH4.8" flash48:latest sh -lc \
  'cd /opt/FLASH4.8 && ./setup WDHydrostatic -auto -3d -site=wd_docker \
   -objdir=objdir_wd && cd objdir_wd && make -j 6'

# run directory needs flash4, flash.par, wd_model.dat, helm_table.dat
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -v "$PWD/run:/work" -w /work flash48:latest sh -lc 'mpirun -np 8 ./flash4'
```
