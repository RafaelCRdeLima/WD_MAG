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

## Status

- [x] Polytrope generator, verified against Lane-Emden tabulated values
- [x] `model_polytrope.dat` written
- [ ] FLASH Simulation unit (reads the model, maps to 3D, self-gravity)
- [ ] FLASH build in a separate objdir on a mounted volume
- [ ] Castro side rebuilt with `EOS_DIR := gamma_law` on the same model
- [ ] `rho_c(t)` comparison
