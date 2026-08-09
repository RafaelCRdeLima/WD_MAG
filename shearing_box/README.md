# Shearing box — measuring what the global runs cannot resolve

## Why this exists

The global runs do not resolve the MRI. Measured from the volume-typical
vertical field, the quality factor $Q = \lambda_{\rm MRI}/\Delta x$ peaks at
$6.5$ around $t = 4.5$ s and sits at $0.3$–$0.5$ over $t = 40$–$78$ s — the
window in which the differential-rotation steepening is measured
(`investigations/mri_wavelength.py`, DIARIO §6.18). So "the differential
rotation survives" and "we have no way to erase it" are not distinguished by
those runs.

Refining does not fix it: cost scales as $N^4$. A local box does, because in
the incompressible approximation the acoustic timestep disappears — and that
approximation is not marginal here but ideal, with $v/c_s \sim 4\times10^{-4}$
and $v_A/c_s \sim 10^{-5}$. Measured below: **80 s per orbit at $64^3$ on 8
threads, so 100 orbits is about 2.2 hours.** The compressible equivalent was
months (DIARIO §6.11c, §6.12).

What the box **cannot** do: the Tayler instability. That is driven by the
curvature of toroidal field lines about the rotation axis, and a local
Cartesian box has no curvature, so the $m=1$ kink cannot grow here by
construction. Since the kink is what destroys the field in the global runs,
nothing in this directory speaks to that result.

## Provenance — settled

**We have the official v6.0, and the mirror turned out to be equivalent.**

Getting it needed a human: the author's page,
`https://ipag.osug.fr/~lesurg/snoopy.html`, sits behind Anubis proof-of-work
anti-bot protection, which exists to stop automated fetching and was not
circumvented. Every path on that host returns the challenge page, Software
Heritage has no copy, and the Grenoble GitLab was unreachable. Rafael
downloaded `snoopy_v6.0.tgz` through a browser, where Anubis clears in a
second. It is GPL-3.0 — the licence was never the obstacle, only the bot
protection.

`src/snoopy-v6.0-official/` is that tarball: 2.4 MB, dated 30 March 2011,
shipping an embedded `.git` whose log runs to *"minor bug in SGS models"* and
*"Added explicit viscosity/resistivity/th diffusion"*.

`src/snoopy-rmoleary-mirror/` is
[github.com/rmoleary/snoopy](https://github.com/rmoleary/snoopy), used while
the official was unreachable. Diffed against it:

| file | difference |
|---|---|
| `src/timestep.c` | **none** |
| `src/shear.c` | **none** |
| `src/gfft.c` | **none** |
| `src/problem/mri/gvars.h` | **none** |
| `src/problem/mri/snoopy.cfg` | **none** |
| `src/mainloop.c`, `particles.c`, `snoopy.c` | particle handling only |

The fork's entire change is a `fld` → `fldi` fix in the particle velocity loop
plus an `init_particles` call under `#ifdef WITH_PARTICLES` — **and
`WITH_PARTICLES` is not defined in the MRI problem, so none of it is
compiled** in our configuration. Both binaries run the smoke test to a
`timevar` that is identical byte for byte.

One earlier claim here was wrong and is corrected: the mirror was dated to
v5.0 from the string in `output_vtk.c`. That string is stale and present in the
official tarball too. The real version is in `src/snoopy.c:189`, and both trees
are v6.0.

`build.sh` now builds the official tree. The mirror is kept only as the record
of this comparison.

## Build

    ./build.sh            # FFTW3 + SNOOPY, problem 'mri'
    ./build.sh shear_dynamo

Nothing leaves this directory and nothing needs root. Two non-obvious points,
both encoded in the script:

- **The include path must go in `CFLAGS`, not `CPPFLAGS`.** `configure` probes
  for `fftw3.h` using `CPPFLAGS` and reports success, but `Makefile.in` never
  substitutes `CPPFLAGS`, so a `CPPFLAGS`-only invocation configures cleanly
  and then fails every single compile.
- **`-malign-double` is dropped** from the stock flags. It is an i386 option
  that on x86-64 changes struct layout against the SysV ABI, and the FFTW built
  here does not use it.

FFTW is built from source because Ubuntu 24.04 ships `libfftw3-double3` but not
`libfftw3-dev` — runtime libraries with no header and no `.so` symlink.

## Running

**Create `data/` first or the run segfaults.** `output_vtk()` does
`fopen("data/vNNNN.vtk", "w")` and never checks the result, so a missing
directory becomes a null `FILE*` and a crash inside `fwrite` at the first
snapshot — after `t = 0` has already been written to `timevar`, which makes it
look like a physics failure rather than a missing directory.

    mkdir -p runs/<name>/data
    cd runs/<name>
    cp ../../src/snoopy-v6.0-official/src/problem/mri/snoopy.cfg .
    OMP_NUM_THREADS=8 ../../src/snoopy-v6.0-official/snoopy

## Parameters, and how they map to our star

SNOOPY is non-dimensional: $\Omega = 1$, box measured in code units. **In
incompressible MHD there is no $\beta$** — the thermal pressure is a Lagrange
multiplier enforcing $\nabla\cdot v = 0$ and the sound speed is infinite. The
$\beta \sim 10^9$ that made a compressible box unaffordable is not a parameter
of the physics at all. What is left:

| SNOOPY key | meaning | our value |
|---|---|---|
| `shear` | $q = -d\ln\Omega/d\ln\varpi$ | **0 to 2** (Komatsu $j$-constant), stock cfg has 1.5 = Keplerian |
| `omega` | $\Omega$, sets the time unit | 1 by convention; one orbit is $2\pi$ |
| `boxsize` | in units where $\lambda_{\rm MRI}$ follows from $B$ | to set from $\lambda_{\rm MRI} = 2.6\times10^6$ cm |
| `reynolds` | $1/\nu$ | — |
| `reynolds_magnetic` | $1/\eta$ | — |

Note the inversion: $\mathrm{Pm} = \nu/\eta =$ `reynolds_magnetic`/`reynolds`.
The stock MRI config has both at 1000, i.e. $\mathrm{Pm} = 1$. For
$\mathrm{Pm} = 0.6$ at $\mathrm{Re} = 1000$, set `reynolds_magnetic = 600`.

**$\mathrm{Pm} \approx 750$, computed for our conditions**
(`investigations/magnetic_prandtl.py`, DIARIO §6.19). An earlier note here
said $\mathrm{Pm} \sim 0.6$, taken from the literature on the convective zone
of a *cool crystallising* CO white dwarf. That is a different object: at our
mean density and an assumed remnant temperature of $10^8$ K, electron
transport gives $\nu = 2.4$ and $\eta = 3.3\times10^{-3}$ cm$^2$/s. The result
is robust — $\mathrm{Pm}$ stays above $30$ even for a Coulomb logarithm of
$5$, and above $190$ across the whole $\rho = 5\times10^7$–$10^9$,
$T = 10^7$–$10^9$ range.

Two consequences, both favourable:

- **The $\mathrm{Pm}_{\rm crit} \sim 2$–$4$ worry is gone.** We are two to
  three orders above it, not below.
- **We are in the same high-$\mathrm{Pm}$ regime as the protoneutron-star
  boxes** where the MInIT parasitic coefficients were calibrated, so they have
  a much better chance of transferring than feared.

**But the box still cannot reach it.** $\mathrm{Pm} = 750$ needs
$\mathrm{Rm} = 750\,\mathrm{Re}$, and a DNS reaching $\mathrm{Re} \sim 10^3$
would need $\mathrm{Rm} \sim 10^6$. So the box brackets rather than matches:
scan $\mathrm{Pm} = 1, 2, 4, 8, 16$ at fixed $\mathrm{Re}$, measure the trend,
and extrapolate with the extrapolation stated. **Transport rises with
$\mathrm{Pm}$ in the disc literature, so a box at $\mathrm{Pm} = 16$ gives a
LOWER BOUND on the transport at $750$** — and a lower bound is enough if it
already erases the differential rotation.

The box therefore does not remove extrapolation, it moves it from $\beta$ to
$(\mathrm{Rm}, \mathrm{Pm})$, in a direction we can sign.

## What comes out

`timevar` carries `vxvy` and `bxby` — the Reynolds and Maxwell stresses, which
are what the closure needs. `em` and `ev` are the magnetic and kinetic
energies.

## Smoke test, 2026-08-08 (official v6.0 build)

Stock MRI problem, $64^3$, box $(4,4,1)$, $q = 1.5$,
$\mathrm{Re} = \mathrm{Rm} = 1000$, run to $t = 12.6$ (two orbits), 8 threads:

| | |
|---|---|
| wall time | 161 s → **80 s/orbit**, ~2.2 h for 100 orbits |
| $E_{\rm mag}$ | $1.04\times10^{-2} \to 0.60$, a factor of 58 |
| Maxwell stress `bxby` | $-0.309$ |
| Reynolds stress `vxvy` | $2.46\times10^{-3}$ |

The MRI grows, the stresses are produced, and the throughput matches the
estimate that made this route viable. This is a smoke test and nothing more:
it is the stock Keplerian disc configuration, not our star, and it has not been
validated against a published run.

One trap for whoever edits `build.sh`: **do not probe the binary for a version
string.** SNOOPY parses no arguments and ignores unknown ones, so
`snoopy --version` does not print and exit — it silently starts a full
simulation in the current directory. That wedged a build until the process was
killed by PID; `pkill -f "snoopy --version"` does not work either, because the
pattern matches the killing shell's own command line and it kills itself.

## Next

1. ~~Recompute $\nu$, $\eta$ and $\mathrm{Pm}$.~~ Done: $\mathrm{Pm} \approx 750$.
2. ~~Download the official v6.0 and diff.~~ Done; see Provenance.
3. Reproduce a published $q = 1.5$ MRI run — validating our use, not the code.
4. Scan $q$, $\mathrm{Pm}$, and the net-toroidal-to-vertical flux ratio, which
   for us runs 2 to 9.5 in amplitude.
5. Extract $\alpha^{\rm PI}$, $\beta^{\rm PI}$ and compare with the published
   $-1.4$ and $-0.8$.
