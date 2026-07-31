# Magnetized ultramassive white dwarfs: the SN Ia / AIC investigation

Working plan for the collaboration with Jorge Rueda, Laura (Chile), Bani and
Zenia (India). Written in English because it is collaboration-facing.

Status: **Phase 1 starting.** Everything in Sects. 3 and 4 is done and
measured; Sects. 5–7 are the plan.

---

## 1. The question

Ultramassive white dwarfs (1.4 to slightly above 2 M☉) carrying strong
magnetic fields — toroidal-dominated interior, poloidal component giving a
surface dipole — evolve under magnetic braking. As support is removed the
star contracts and the central density rises until it meets one of two
fates:

- **carbon ignition** → thermonuclear runaway → overluminous SN Ia
- **electron capture** → loss of pressure support → collapse to a neutron star

The collaboration has a first set of results using a *phenomenological*
interior field, some prescribed function of the matter density. Our
contribution is to replace that with a field obtained from Maxwell's
equations — the Grad–Shafranov / self-consistent-field machinery in this
repository — and to say what that substitution changes.

## 2. Architecture: this is secular, not hydrodynamic

Magnetic braking acts on Myr–Gyr timescales. That cannot be integrated in a
hydrodynamics code. The evolution must therefore be a **quasi-static
sequence of equilibria** parametrized by decreasing angular momentum, with a
*local* thermal criterion evaluated at each step:

```
    family of equilibria (M, J, magnetic flux)        <- our contribution
              |
              v
    secular walk: J decreases by braking torque       <- sequence, not simulation
              |
              v
    at each step: eps_nuc(C+C, screened) vs eps_nu    <- FLASH microphysics as a library
              |
              v
    runaway condition d eps_nuc/dT > d eps_nu/dT      <- the branch point
```

A hydrodynamics code (FLASH, Castro) earns its place in exactly two places:
**testing whether a configuration is dynamically stable at all** (Phase 1),
and **verifying the critical configuration** once it is identified (Phase 5).
Not in the evolution itself.

This should be confirmed with Jorge explicitly. It is the kind of assumption
that is cheap to align on now and expensive to discover later.

---

## 3. What we already have

### 3.1 Structure solver (`scf/`)

| module | what it does |
|---|---|
| `scf.py` | Hachisu self-consistent-field iteration, pluggable force terms |
| `poisson.py` | Poisson solve, spherical harmonics, `lmax` truncation |
| `gradshafranov.py` | Δ\* inversion on the P_l^1 angular basis, radial Green's function |
| `eos.py` | zero-temperature degenerate EOS, `x(H)` inversion, neutronization threshold |
| `diagnostics.py` | virial error, magnetic energies, surface dipolarity, equatorial/polar radii |
| `terms/toroidal_sc.py` | self-consistent toroidal branch, B_φ = K ρ^m ϖ^(2m−1) |
| `terms/poloidal.py`, `terms/rotation.py` | imposed poloidal, rotation laws |
| `castro_model_writer.py` | **1D radial exporter — field-free and spherically symmetric only** |

Note on the last row: it is not a matter of "opening a gate". The writer
asserts `check_field_free_non_rotating` *and* `check_spherical_symmetry`,
then averages over θ to produce a 1D profile. It discards precisely what we
need to export. Phase 1 requires a new 2D axisymmetric writer.

`gradshafranov.py` carries a corrected radial recursion (ratio form, not
absolute powers of r) with a regression test suite in
`scf/tests/test_gradshafranov_overflow.py`. The original form lost the outer
term entirely at `lmax = 16` for sources that do not vanish at the origin.

### 3.2 Investigations (`investigations/`), all with committed data

| script | what it establishes |
|---|---|
| `mixed_2msun.py` | the certified 2 M☉ configuration and the poloidal trade-off |
| `vector_potential_export.py` | the A → curl → B construction gives ∇·B = 1.7e−16 |
| `barotropic_ceiling.py` | the barotropic twisted-torus ceiling, measured here |
| `confinement_cost.py` | what confining the toroidal field costs in peak field |
| `braking_sequence.py` | ρ_c vs magnetic support at fixed M = 2 M☉ |
| `confinement_gain_sequence.py` | the confinement gain per point along that sequence |
| `nonbarotropic_equilibrium.py` | the non-barotropic solver, validated and incomplete |

### 3.3 Papers (`papers/`)

- `wd-toroidal-poloidal/` — the certified 2 M☉ mixed configuration
- `wd-toroidal/` — the toroidal-only theory, including a rotating result
- `wd-braithwaite-relaxation/` — dynamical relaxation of a random seed field
- `wd-nonbarotropic/` — beyond barotropy (started this week)

### 3.4 Hydrodynamics

**Castro** (`castro/`, `amrex/`, `microphysics/`, `castro_problems/wd_braithwaite`).
MHD via constrained transport, single-level only. Known constraints already
mapped: half-shift geometry (a cell *centre* at r = 0), the gravity r = 0
patch, the well-balancing gap, and `small_dens = 1.0e4`.

**FLASH 4.8** (`flash_crosscheck/FLASH4.8`, plus Docker images
`flash48:latest`, `flash48-ccsn:latest`, `flash48_ccsn_neutrinos:latest`).
Verified present in the source tree — this is the decisive asset for this
project, because the physics the collaboration wants to evolve is already
there:

| unit | what it provides |
|---|---|
| `Burn/BurnMain/nuclearBurn/{Iso7,Aprox13,Aprox19}` | networks including ¹²C+¹²C |
| `bn_sneutx.F90` | neutrino losses, Itoh et al. (1996): pair, plasma, photo, bremsstrahlung, recombination |
| `bn_screen4.F90` | screening, Graboske et al. (1973) weak + Alastuey & Jancovici (1978) strong |
| `bn_ecapnuc.F90`, `bn_mazurek.F90` | electron capture — the collapse channel |
| `Eos/EosMain/Helmholtz` | finite-temperature EOS |
| `sourceTerms/{Cool,Deleptonize}` | auxiliary sinks |

Our own setup `flash_crosscheck/WDHydrostatic/` (Config, Simulation_data,
Simulation_init, Simulation_initBlock, Simulation_initSpecies,
Simulation_adjustEvolution, Makefile, flash.par) plus `make_wd_model.py`,
`compare_codes.py`, and `FLASH_CORE_PATCHES.md` recording the core patches
that were required.

Obstacles already hit and documented, so they are not rediscovered:
`iso_c_binding` needs `use, intrinsic ::`; objdir ownership under Docker
needs `--user`; `helm_table.dat` must be present; Coulomb corrections can
drive the pressure negative; and the near-vacuum ambient is too hard for the
Helmholtz inversion to invert.

---

## 4. What we have measured

Every number below is from a committed script with committed output.

### 4.1 The certified 2 M☉ configuration

| quantity | value |
|---|---|
| ρ_c, μ_e | 10⁹ g cm⁻³, 2 (below the 1.94×10¹⁰ neutronization threshold) |
| M | 2.0072 M☉ |
| virial error | 9.9×10⁻⁵ (gate: < 10⁻³) |
| E_tor/\|W\| | 0.203 |
| max\|B_φ\| | 4.28×10¹³ G = 0.969 B_c |
| exterior dipole | 5.9×10¹⁰ G (at k₀ = 10⁻¹³) |
| B_t/B_p | 884 in energy, 26 in amplitude |
| R_eq, R_pol | 5.207×10⁸, 5.736×10⁸ cm (prolate, ratio 1.10) |

Mesh convergence: mass changes 0.15%, 0.04%, 0.01% over four halvings of Δr,
converging to 2.007 M☉.

### 4.2 What barotropy forbids

The azimuthal Lorentz force must vanish on its own, which forces
ϖB_φ = β(u) and confines the toroidal field to the closed-line region
u > u_s. **That survives any equation of state.** What barotropy adds is
that f_L/ρ must be a gradient, pinning the poloidal source to ρ M′(u).

Measured consequence, on a field-free 1.346 M☉ star:

| ζ | max converged E_tor/E_mag |
|---|---|
| 1.0 | **0.040** |
| 1.1 | 0.029 |
| 2.0 | 0.003 |

At that ceiling E_tor/\|W\| = 10⁻⁵, against the 0.203 that 2 M☉ requires —
**short by four orders of magnitude**. A self-consistent barotropic mixed
equilibrium cannot support an ultramassive white dwarf. The closed region is
37.7% of the stellar volume, so this is not a lack of space; it is the form
of the source.

### 4.3 Confinement is cheaper than the volume argument suggests

Reproducing E_tor/\|W\| = 0.203 inside ~37% of the volume needs a *lower*
peak field than the space-filling law, because peak field at fixed energy is
set by how peaked the profile is, not only by how much volume it has. The
gain, recomputed per point along the sequence, runs 2.69 → 2.93.

### 4.4 The braking sequence, at fixed M = 2 M☉

Eleven points, ρ_c from 8×10⁸ to 2×10¹⁰ g cm⁻³.

**The sequence is soft.** K_tor moves 17% and E_mag/\|W\| 21% (0.208 → 0.164)
while ρ_c moves by a factor of 25 — the marginal polytrope showing through,
Γ ≈ 4/3, mass nearly independent of ρ_c. For the scenario: once braking
starts removing support, the star crosses the whole window quickly. It does
not linger.

**Where the lines fall:**

| line | ρ_c (g cm⁻³) | E_mag/\|W\| | space-filling | confined |
|---|---|---|---|---|
| C+C pycnonuclear | 3.0×10⁹ | 0.181 | 1.98 B_c | **0.72 B_c** |
| ²⁰Ne capture | 9.6×10⁹ | 0.168 | 4.25 B_c | 1.49 B_c |
| ¹⁶O capture | 1.94×10¹⁰ | 0.164 | 6.73 B_c | 2.30 B_c |

with B_c = 4.414×10¹³ G. The confined branch crosses B_c at ρ_c = 5.1×10⁹.

**This splits the two branches, and the split is at the composition.**

- **C/O core** → meets C+C first, at 3×10⁹, and never reaches the capture
  densities. Confined, the ignition point sits at 0.72 B_c, *below* the
  Landau threshold, so a field-independent EOS is defensible. This branch is
  calculable with Skye or PC (finite temperature + Coulomb).
- **O/Ne core** → no carbon to burn, goes straight to ²⁰Ne capture at
  9.6×10⁹, which is 1.49 B_c even confined. This branch needs a *magnetized*
  EOS with Landau quantization and P_⊥ ≠ P_∥.

C/O at these masses implies a **merger origin**; single-star evolution
produces O/Ne. This is a premise of the SN Ia scenario, not a detail.

### 4.5 The export path is de-risked

For an axisymmetric field the vector potential is analytic:
A_φ = u/ϖ and A_z = −∫₀^ϖ B_φ dϖ′. Taking the discrete curl on the staggered
mesh makes ∇·B = 0 an *identity of the discretization* rather than something
to clean up. Measured on a 64³ mesh: **∇·B = 1.7×10⁻¹⁶** (normalized), with
83.7% of the peak field retained by the Cartesian sampling.

---

## 5. Open tensions and risks

**R1 — Tayler instability. This is the largest risk and it is quantitative.**
The scenario wants a surface dipole of ~10⁹ G together with a toroidal field
strong enough to hold 2 M☉. Since E_pol ∝ k₀² at fixed E_tor, and we measure
E_t/E_p = 884 at B_pole = 5.9×10¹⁰ G, scaling to B_pole = 10⁹ G gives

> **E_tor/E_pol ≈ 3×10⁶**

which is the configuration most exposed to the Tayler m = 1 instability. And
the braking torque goes as B_pole², so a weak dipole *also* makes braking
slow, giving the instability more time. Both effects push the same way. If
the configuration rearranges on an Alfvén time (seconds to minutes), it never
sits still long enough to brake for Myr and the scenario needs rethinking.

Worth asking Jorge whether 10⁹ G is an observational constraint or a
conservative choice: confinement lets us reach B_t/B_p ≈ 1 in amplitude with
*both* peaks at ~0.36 B_c, and a 10¹¹–10¹² G dipole would be both more stable
and faster-braking.

**R2 — the relaxation measurement window.** Our Braithwaite study reached only
~0.23 Alfvén crossing times before the window closed. Any stability statement
from that study is bounded by this, and Phase 1 must do better or say so.

**R3 — the collapse branch is above B_c even confined.** Section 4.4. If part
of the population is O/Ne, that fraction requires the expensive EOS.

**R4 — the non-barotropic equilibrium is underdetermined.** The free function
(P on the polar axis, one value per equipotential) is genuine physical
freedom — a real star's stratification comes from its thermal and
compositional history. For one prescribed field there is a one-function
family of equilibria with different masses. Closing it needs a physical
stratification condition (convective stability is the candidate), not a
better solver.

**R5 — EOS/screening inconsistency in the default FLASH combination.** FLASH
uses Helmholtz (no Coulomb) for the EOS and Alastuey–Jancovici for screening;
they are not derived from the same free energy, and the classical AJ form
overestimates the enhancement in the strongly quantum regime where an
ultramassive crystallized core lives. Since the ignition density depends
exponentially on the enhancement factor, this feeds straight into the branch
point. Skye and PC both derive EOS and screening from one free energy.

---

## 6. The phases

### Phase 0 — settle the premises *(cheap, email)*

- Composition: C/O or O/Ne? Decides the branch **and** the EOS requirement (4.4).
- Is the 10⁹ G dipole a hard constraint? See R1.
- Confirm the secular architecture of Sect. 2.
- Agree the mass range and the field ranges.

### Phase 1 — stability. **This is go/no-go.** *(starting now)*

If the configuration that supports 2 M☉ does not survive on an Alfvén time,
everything downstream is moot. This is prior to the EOS work, and cheaper.

1. **A 2D axisymmetric model exporter.** New code: write (ρ, A_φ, A_z) on the
   meridional grid, plus provenance. The existing 1D writer cannot be
   extended — it asserts spherical symmetry and averages over θ.
2. **Problem-side reader and curl.** The consuming problem interpolates A onto
   cell edges and takes the discrete curl onto faces, so ∇·B = 0 by
   construction. Already proven in Python (4.5); this transcribes it.
3. **3D MHD run**, looking for the m = 1 Tayler mode. Diagnostics: the
   azimuthal mode decomposition of the field, E_tor/E_mag against time, and
   the elapsed time in Alfvén units — which must be reported, given R2.
4. **Verdict**: does the toroidal-dominated configuration survive, and at what
   B_t/B_p does it become stable if it does not?

Exit criterion: a measured survival time in Alfvén units for at least two
values of B_t/B_p, one of them the value the 10⁹ G dipole implies.

### Phase 2 — the equilibrium family

The braking track is at **fixed mass and fixed flux, varying rotation** — not
what Sect. 4.4 computed, which varies the magnetic support and therefore maps
the locus rather than the track. Rotation is already in the solver and a
rotating configuration is already certified in `papers/wd-toroidal/`.

Deliverable: a family of confined mixed configurations over (M, J), exported
in the Phase 1 format.

### Phase 3 — the ignition criterion

FLASH's routines used **as a library**, not as a simulation: `bn_burner` with
`bn_screen4` for C+C, `bn_sneutx` for Itoh neutrino losses, `bn_ecapnuc` for
the capture channel. Evaluate ε_nuc against ε_ν and the runaway condition at
each point of the Phase 2 family.

**This is the only phase that requires a finite-temperature EOS.** The
structure of a degenerate ultramassive white dwarf is weakly sensitive to T,
so Phases 1 and 2 can proceed on the current cold EOS and the two work
fronts stay decoupled. Recommendation for this phase: **Skye** (differentiable,
self-consistent crystallization, Coulomb and screening from one free energy),
or **PC** if a self-contained port without MESA is preferred. Not HELM alone,
per R5.

### Phase 4 — the branch outcome

Walk the Phase 2 family under a braking torque, apply the Phase 3 criterion,
record which line is crossed first and at what time. This is the scientific
result.

### Phase 5 — hydrodynamic verification

Take the critical configuration and verify it in FLASH: that it is
dynamically stable up to the branch point, and that ignition proceeds as the
local criterion predicts. Here the known obstacles of Sect. 3.4 return.

---

## 7. What we need from the collaboration

- The answers in Phase 0.
- Their current phenomenological B(ρ) prescription, so the substitution can be
  compared like for like rather than in isolation.
- Their braking torque prescription and the assumed initial rotation.
- Whether the ultramassive configurations are to be built by us (SCF) or by
  them (FLASH), and in which direction the model files should flow.
