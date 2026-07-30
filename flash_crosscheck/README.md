# Castro vs FLASH: does an Eulerian Godunov scheme hold a hydrostatic star?

## The question this exists to answer

The Braithwaite measurements in `papers/wd-toroidal-poloidal/` are bounded
by a validity window of `~0.23` Alfvén crossing times, and that window is
set by the background star's own drift: `rho_c` leaves its `±2%` band by
`t/t_dyn = 1.128` and never settles (paper Fig. 1). The open question is
whether that drift is a Castro-specific defect or generic to Eulerian
Godunov schemes evolving a self-gravitating hydrostatic star.

It decides something concrete: whether migrating to another code is worth
it. Castro's MHD is single-level only (`Docs/source/mhd.rst`), so AMR is
unavailable, which is a real reason to look at FLASH — but only if FLASH
holds the star better. If it drifts the same way, the window is a property
of the method and the migration buys nothing.

## Why a gamma-law polytrope, and not the real star

FLASH ships `Gamma`, `Helmholtz`, `Multigamma` and `Tabulated` EOS units.
There is no `ztwd`. Giving Castro a ztwd star and FLASH a Helmholtz star
would mean any difference in `rho_c(t)` could be the equation of state
rather than the scheme — the same confound that caused the Step 3 collapse
(`docs/teoria.md` Sec 6.4), diagnosed there only after it had been blamed
on the magnetic field.

A gamma-law polytrope exists in *both* codes (Castro's Microphysics
`gamma_law`, FLASH's `Eos/Gamma`), so the EOS drops out and only the
scheme is left.

`n = 3/2` (`gamma = 5/3`), not `n = 3`: an `n = 3` polytrope is marginally
stable, so a drifting `rho_c` could be genuine near-neutral physics rather
than a numerical artifact. `n = 3/2` is unambiguously dynamically stable,
so any drift measured is the scheme.

## The star, and the one caveat that matters

`make_polytrope.py` solves Lane-Emden and writes `model_polytrope.dat` in
the single-header-line format Castro's `model_parser` already reads (the
same format `scf/castro_model_writer.py` emits), so the Castro side needs
no new reader and the FLASH side reads the same file.

Verified against the textbook values for `n = 3/2`:
`xi1 = 3.653754`, `-xi1^2 theta'(xi1) = 2.714055`.

Scaled to the same central density and radius as the wd_braithwaite
background star, so that **the number of cells across the star is
identical to the production run** (30.8 at `64^3` on the same domain) --
resolution being the parameter most likely to drive the drift:

| | polytrope | wd_braithwaite (ztwd) |
|---|---|---|
| `rho_c` (g/cm^3) | `9.884e8` | `9.884e8` |
| `R` (cm) | `2.36e8` | `2.36e8` |
| `M` | `4.57 Msun` | `1.35 Msun` |
| `t_dyn` (s) | `0.1472` | `0.2758` |
| `rho_c/rho_mean` | **`6.0`** | **`~170`** |

**The caveat: this polytrope is an easier problem than the real star.**
Its central condensation is `6` against roughly `170` for the ztwd star,
and central condensation is exactly what makes hydrostatic balance hard to
hold on a grid. This is not a tuning choice that could be fixed -- for a
single-index polytrope `rho_c/rho_mean` is a function of `n` alone
(`5.99` at `n = 3/2`, `54` at `n = 3`), and reaching `170` requires
`n > 3`, which is the dynamically unstable regime this test is
specifically constructed to avoid. At `rho_c ~ 1e9` no single-gamma
polytrope resembles a white dwarf: the physically-motivated
non-relativistic degenerate polytrope (`K` fixed by `mu_e = 2`) comes out
at `~16 Msun`, because non-relativistic degeneracy is badly invalid at
this density.

**So the test is asymmetric in what it can conclude, and that is stated up
front rather than discovered later:**

- **Both codes drift alike → decisive.** The window is generic to the
  method class, and migrating codes does not fix it.
- **FLASH holds the star → inconclusive.** It may only mean the polytrope
  was too easy. The harder test (matching central condensation, hence
  matching the EOS) would then be required before believing it.

The cheap experiment is worth running because the decisive outcome is also
the one that saves the most work.

## What is measured

Field-free, self-gravity only, no MHD, no seed field. `rho_c(t)` from each
code, in units of that star's own `t_dyn`, over `4 t_dyn`. Same domain
(`+-4.90e8` cm), same `64^3`, same gamma, same `abar`/`zbar`, same 1D
model file.

## Environment

FLASH 4.8 is not extracted on the host: it lives in the Docker image
`flash48:latest`, tree at `/opt/FLASH4.8`, with `sites/rafael` configured
and toolchain from `/opt/micromamba/envs/flash_gcc12` (`mpif90`, `h5pfc`).
The image's `object/` is an existing CCSN build and is left alone -- work
uses a separate `-objdir` on a mounted volume, or nothing survives the
container.

## Superseded: the polytrope route

`make_polytrope.py` and the gamma-law design above were built first and are
kept because the reasoning still applies if the Helmholtz route stalls. It
was superseded on instruction: use the EOS most appropriate to a white
dwarf that FLASH actually has, which is **Helmholtz**. See
`make_wd_model.py` and `WDHydrostatic/`.

The EOS-mismatch problem that motivated the polytrope was solved a
different way: the star is integrated in HSE with the project's own ztwd
EOS at a temperature low enough that Helmholtz's extra terms are
negligible, and `make_wd_model.py` MEASURES that
(`P_ion/P_deg = 9.8e-4`, `P_rad/P_deg = 1.2e-12` at `T = 1e7` K) instead of
assuming it. That also removes the central-condensation caveat: the star
is the real ztwd structure, validated to `0.026%` in mass against the
production star.

## Result with equivalent damping

`WDHydrostatic/Simulation_adjustEvolution.F90` reproduces Castro's
`problem_source.H`: global velocity damping on every cell,
`rate = 1/(0.2 t_dyn)`, cosine ramp-off, kinetic-energy-consistent energy
update, with the reference run's window (`[0, 20 t_dyn]`, ramp from 18) so
damping is fully on throughout. Two deliberate differences, both in the
file's header: it is operator-split (FLASH's driver has no in-hydro source
hook, so it goes through `Simulation_adjustEvolution`, which the driver
calls every step for exactly this) and the per-step decay is `exp(-r dt)`
rather than Castro's explicit `1 - r dt`.

Central density deviation, both damped, same law, same parameters:

| `t/t_dyn` | FLASH | Castro |
|---|---|---|
| 0.418 | `-1.02%` | `-0.49%` |
| 0.583 | `-2.64%` | `-1.03%` |
| 0.767 | `-4.63%` | `-1.43%` |
| 0.863 | `-5.63%` | `-1.61%` |

**FLASH drifts about 3.5x more than Castro at the same `t/t_dyn`, and then
fails at `t/t_dyn = 0.876`** -- timestep collapse to `1e-12` with the EOS
driven off-table. The failing cell sits at `r = 2.7e8` cm, just outside
`R = 2.46e8`: the star/vacuum interface, not the interior. Data in
`flash_rhoc_damped.csv` (undamped first attempt in `flash_rhoc.csv`).

**This does not support migrating, but it is not a fair fight either.**
The Castro setup carries real tuning this FLASH setup does not: the
discretization chain of Sec 6.6 (half-shift geometry, volume averaging),
the gravity `r=0` patch of Sec 6.7, and floors tuned over many runs. This
FLASH problem is a first attempt whose choices were made to get it to run
at all -- ambient at `1 g/cm^3` (four orders above Castro's, forced by
Helmholtz refusing to invert a harder vacuum), Coulomb correction off,
`mpole_lmax = 0`, no exterior sponge -- and whose initial mapping already
lands `-4.7%` off the 1D central density, the same class of error Sec 6.6
fixed on the Castro side.

So the honest reading: on the question that motivated this, *does another
code hold the star better*, the first answer is no. Before that becomes a
conclusion, the surface/vacuum treatment has to be fixed, because that is
what is killing the run and it is also the most likely source of the extra
drift.

## Status


- [x] FLASH 4.8 extracted from the Docker image, patched, and built with
      Helmholtz + Poisson/Multipole + unsplit hydro (`FLASH_CORE_PATCHES.md`)
- [x] `WDHydrostatic` simulation unit; 1D ztwd model validated against the
      production star (M to 0.026%)
- [x] First run: 39 steps to `t/t_dyn = 0.773`, then a timestep collapse
      and an out-of-table EOS state -- `flash_rhoc.csv`
- [x] Equivalent damping added and the comparison made like-for-like:
      FLASH drifts ~3.5x more, and fails at `t/t_dyn = 0.876`
- [ ] Fix the star/vacuum interface, which is where the run dies and the
      likeliest source of the extra drift: raise/soften the ambient,
      `eos_forceConstantInput`, an exterior sponge, or `use_hybridRiemann`
- [ ] Port the Sec 6.6 mapping fix (the IC starts `-4.7%` off `rho_c`)
- [ ] Only then treat the drift ratio as a statement about the schemes
