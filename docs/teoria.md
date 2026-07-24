# Theory behind the `wd-magnetizada` dashboard

> This is the plain-text version (source of truth, good for `git diff`
> and terminal reading). For a PDF with book-quality typesetting (real
> equations, not verbatim text), see `teoria.tex` — compile with
> `xelatex teoria.tex` (twice, for the table of contents) — or open the
> already-compiled `teoria.pdf`. The content is the same in all three;
> `teoria.tex` was written by hand from this file, not auto-generated, so
> content changes must be replicated manually in both.

This document explains what the dashboard computes: the equations behind
each number and each figure, and the code that implements them. It is not
an introduction to stellar physics or MHD — it assumes the reader already
knows that. What this document does is bridge the theory (section 1) and
the real code (`scf/`, `dashboard/`).

See `plano_wd_magnetizada.md` for the full project plan (decisions
D1–D6, architecture, phases). This document is about the *implemented
physics*, not the work plan.

---

## Table of symbols and units

| Symbol | Meaning | Unit (CGS) | Name in code |
|---|---|---|---|
| `r`, `θ` | spherical coordinates | cm, rad | `r`, `theta` |
| `ϖ` | cylindrical radius, `ϖ = r sinθ` | cm | `omega` (see note below) |
| `z` | height, `z = r cosθ` | cm | `z` (only in `plots.py`) |
| `ρ` | mass density | g cm⁻³ | `rho` |
| `x` | normalized Fermi momentum, `x = p_F/(m_e c)` | dimensionless | `x` |
| `P` | pressure | dyn cm⁻² = erg cm⁻³ | `pressure(x)` |
| `H` | specific enthalpy | erg g⁻¹ = cm² s⁻² | `H` |
| `Φ` | gravitational potential | erg g⁻¹ = cm² s⁻² | `Phi` |
| `A_φ` | φ component of the vector potential | G cm | `A_phi` |
| `u` | flux function, `u = ϖ A_φ` | G cm³ | `u` |
| `B_r, B_θ, B_φ` | magnetic field components | G (gauss) | `Br`, `Bth`/`Btheta`, `Bphi` |
| `f(u)` | poloidal current function, `f(u) = k₀` | g^(-1/2) cm^(1/2) s^-1 | `k0` |
| `M(u)` | poloidal Bernoulli potential, `M(u) = k₀ u` | erg g⁻¹ | `M_u` (local variable in `scf.py`) |
| `β(u)` | toroidal current function, `β = ϖ B_φ` | G cm | not computed separately — `B_φ` already comes out of `impose_toroidal()` |
| `ζ` | amplitude of the imposed toroidal field | (units of `u^{-m}` × G) | `zeta` |
| `m` | exponent of the imposed toroidal field | dimensionless, integer ≥ 1 | `m_tor` |
| `u_c` | value of `u` on the last closed line | G cm³ | `u_c` |
| `C` | Bernoulli constant | erg g⁻¹ | `C` |
| `ρc` | central density (input parameter) | g cm⁻³ | `rho_c` |
| `k₀` | poloidal field amplitude (input parameter) | g^(-1/2) cm^(1/2) s^-1 | `k0` |
| `μₑ` | mean molecular weight per electron, `Y_e = 1/μₑ` | dimensionless | `mu_e` |
| `M` | total stellar mass | g (displayed in M_sun) | `M` (⚠ collides with `M(u)` in the source notation — see note) |
| `R_eq`, `R_pol` | equatorial and polar radii (where `H=0`) | cm (displayed in km) | `R_eq`, `R_pol` |
| `W` | gravitational energy | erg | `W` |
| `Π` | internal energy, `∫P dV` | erg | `Pi` |
| `E_pol`, `E_tor`, `ℳ` | poloidal, toroidal, total magnetic energies | erg | `E_pol`, `E_tor`, `E_mag` |
| `VE` | virial error | dimensionless | `VE` |
| `G` | gravitational constant | 6.674×10⁻⁸ cm³ g⁻¹ s⁻² | `G_CONST` |
| `A`, `B` | EOS constants | see §1.1 | `A_CONST`, `B_of_mu_e(mu_e)` |
| `l` | degree of the Legendre expansion | integer ≥ 0 (Poisson) or ≥ 1 (GS) | `l` |
| `t_dyn`, `v_A`, `t_Alfven` | derived time/velocity scales | s, cm/s, s | `t_dyn`, `v_A`, `t_alf` |

**Notation notes:**
- The code uses the variable `omega` for `ϖ` (cylindrical radius), **not**
  for angular velocity — there is no rotation in this project (D3, static
  star). It is a naming choice specific to this code; whenever you read
  `omega` in `scf.py`/`gradshafranov.py`, it always means `ϖ = r sinθ`.
- `M` is used in the source theory for two different objects: the total
  stellar mass (§3.7 of the original prompt) and the Bernoulli potential
  `M(u)` (§3.4–3.5). The code resolves the collision by naming the mass
  `M` (returned by `scf.total_mass()`) and the potential `M_u` (a local
  variable inside `hachisu_scf()`, never exposed outside the function). In
  this document, whenever there is risk of ambiguity, the text spells out
  "mass" or "Bernoulli potential" in full.
- `Bth` appears as `Btheta` in `diagnostics.py` (function parameter) and
  `Bth` on the dashboard pages (local variable) — same object, two names,
  because of `flake8`/readability in the two contexts.

---

## 1. Common theoretical core

### 1.1 Equation of state

Fully degenerate electron gas at T = 0:

```
P = A [ x(2x² − 3)(x² + 1)^{1/2} + 3 sinh⁻¹ x ]
ρ = B x³
x ≡ p_F / (m_e c)

A = 6.01 × 10²²  dyn cm⁻²
B = 9.82 × 10⁵ / Y_e  g cm⁻³     (Y_e = 0.5 → B ≈ 1.96 × 10⁶)
```

`A` is a fixed physical constant (a combination of `m_e`, `c`, `ℏ`); `B`
depends on composition only through `Y_e` (electrons per baryon). The
dashboard parametrizes by `μₑ = 1/Y_e` instead of `Y_e` directly.

→ `scf/eos.py :: pressure()` (P(x)), `density()` (ρ(x)), `B_of_mu_e()` (B(μₑ))

Specific enthalpy, analytic — this is what makes the SCF possible:

```
H = ∫ dP/ρ = (8A/B) [ √(1 + x²) − 1 ]
```

→ `scf/eos.py :: enthalpy()`

Inverse, applied at every SCF iteration:

```
x = √[ (1 + HB/8A)² − 1 ]        (H ≤ 0 => ρ = 0)
```

→ `scf/eos.py :: x_of_enthalpy()`, `density_of_enthalpy()`

The `−1` term normalizes H = 0 at the surface (the ρ=0 boundary of the star).

**Effective polytropic index.** Near the surface `x ≪ 1` gives `ρ ∝ H^{3/2}`,
i.e. `n = 3/2`. In the relativistic core `x ≫ 1` gives `ρ ∝ H³`, i.e.
`n = 3`. The EOS slides between the two, and `n = 3` is the marginal case
where the mass becomes independent of the central density — the origin
of the Chandrasekhar limit. This fact is the direct cause of the
parametrization decision in §1.4/§1.5.

**Sound speed.** Added after the original theoretical spine (it was not in
the plan), because it is needed for the amplitude of the perturbation
exported in Tab 3:

```
c_s = √(dP/dρ) ,   dP/dx = 8A x⁴/√(1+x²)   (analytic derivative of P(x))
```

→ `scf/eos.py :: sound_speed()`, validated against a finite difference in
`scf/tests/test_sound_speed.py`

---

### 1.2 Self-gravity

```
∇²Φ = 4πGρ
```

Solved by expansion in Legendre polynomials. With
`ρ(r,θ) = Σ_l ρ_l(r) P_l(cosθ)`:

```
Φ_l(r) = − (4πG / (2l+1)) [ r^{−(l+1)} ∫₀^r ρ_l(r′) r′^{l+2} dr′
                          + r^{l}     ∫_r^∞ ρ_l(r′) r′^{1−l} dr′ ]
```

The weights `r′^{l+2}` and `r′^{1−l}` carry the `r′²` of the volume element.

→ `scf/poisson.py :: solve_poisson()` (exact implementation of this
formula, D_l/E_l computed by cumulative trapezoidal sums); `legendre_matrix()`
computes `P_l(cosθ)`.

Validated against the closed-form solution of the uniform sphere and
Newton's shell theorem in `scf/tests/test_poisson.py`.

`G = 6.674 × 10⁻⁸ cm³ g⁻¹ s⁻¹` — a fixed physical constant (rounded CODATA
reference value), not a project parameter. → `scf/poisson.py ::
G_CONST`, repeated in `dashboard/units.py :: G_CONST` (same value, two
modules, because `poisson.py` cannot depend on the dashboard — R1).

---

### 1.3 Magnetic field and the flux function

With `ϖ = r sinθ`, the flux function is defined as

```
u = ϖ A_φ
```

and the axisymmetric field splits into poloidal plus toroidal:

```
B = (1/ϖ) [ ∇u × ê_φ  +  β(u) ê_φ ]
```

In spherical components:

```
B_r = (1 / (r² sinθ)) ∂u/∂θ
B_θ = − (1 / (r sinθ)) ∂u/∂r
B_φ = β(u) / ϖ
```

→ `scf/diagnostics.py :: poloidal_field()` implements `B_r`, `B_θ` (via
finite differences of `u`, protected against the coordinate singularity
at `ϖ→0`). `B_φ` is not computed from an intermediate `β(u)`; it comes
directly out of `scf/toroidal.py :: impose_toroidal()`, which already
implements the chosen functional form for `β` (see §1.9) substituted into
the formula above.

The general form `B = (1/ϖ)[∇u×ê_φ + β ê_φ]` does not appear anywhere in
the code as a vector identity — the code goes straight to the three
scalar components.

**Geometric interpretation: the contours of `u` are the poloidal field
lines.** That's why Tab 1 plots them (§2).

---

### 1.4 Grad-Shafranov

```
Δ* u = − 4π ϖ² ρ f(u) − β β′(u)

Δ* = ∂²/∂ϖ² − (1/ϖ) ∂/∂ϖ + ∂²/∂z²
```

`f(u) = dM/du` is the current function generating the poloidal field;
`β(u)` generates the toroidal field. The simplest choice, and the one used
here, is `f(u) = k₀` constant (Lander & Jones 2009). Since the toroidal
field is imposed after the poloidal field has converged (§1.9, D6) and
does not enter the GS equation actually solved by the SCF, the term
`−ββ′(u)` **does not appear in the source the code actually assembles** —
the implemented source is only `−4πϖ²ρk₀`.

→ `scf/gradshafranov.py :: solve_gradshafranov()` — source passed by
`scf/scf.py :: hachisu_scf()` as `source = -4*np.pi*omega2*rho*k0`.

The current density corresponding to this source is

```
J_φ = c ρ ϖ f(u)
```

This formula **does not correspond to any function in the code** — `J_φ`
is never computed as a named quantity. It was used only analytically,
during the debugging that found the bug described in the box below, to
assemble the right-hand side of the integral Ampère's-law test (§1.4,
derivative-free tests).

**Solution by associated-Legendre expansion.** The operator `Δ*` separates
into `P_l¹(cosθ)` functions, **not** `P_l(cosθ)`. The radial structure of
the Green's function is analogous to that of Poisson, but the integral
weights are **not the same**:

```
weights in Poisson (for Φ):         r′^{l+2}   and   r′^{1−l}
weights in Grad-Shafranov (for u):  r′^{l+1}   and   r′^{−l}
```

→ `scf/gradshafranov.py :: solve_gradshafranov()`, lines where `D_l`/`E_l`
are assembled with `r ** (l + 1)` and `r ** (-l)`.

> **Why it's like this — the Green's-function exponents (G4).**
> The difference of one power comes from the `ϖ` weight of the GS source
> compared to the `r′²` of the volume element in Poisson, and depends on
> whether the equation is solved for `u` or for `A_φ`. **This was a real
> bug in this project:** an earlier version of `gradshafranov.py` used
> `r′^{l+2}` and `r′^{1−l}` (mistakenly copied from the structure of
> `poisson.py`, which has one extra power of `r` because of the
> substitution `χ = rΦ_l` used in reducing the scalar Laplacian — `Δ*`
> does not need that substitution). The bug inflated the poloidal field by
> a factor of ~7000× and the magnetic energy by ~5×10⁷.
>
> The bug survived piece-by-piece validation — including a "by-hand"
> closed form — because that closed form was built by reusing the same
> indicial equation (`r^{l+1}`, `r^{−l}` for the homogeneous solutions,
> which **were correct**) and therefore inherited the same normalization
> of the internal integrals, which was wrong. Only two tests that use
> neither the Green's function nor a second derivative — consistency via
> Ampère's law — revealed the fixed factor. See the full note (with the
> reasoning chain) at the top of `scf/gradshafranov.py`.
>
> **Operational lesson:** "solve with the formula, check against the same
> formula" is not independent validation, even across different
> files/functions, if both inherit the same normalization from a shared
> derivation.

**Derivative-free normalization tests** — the defense against this class
of error:

```
Flux:     u(ϖ,z) = ∫₀^ϖ B_z(ϖ′,z) ϖ′ dϖ′
Ampère:   oint_C B·dl = 4π k₀ ∫_S ρ ϖ dA
```

- **Ampère**: implemented and in the repository → `scf/tests/test_gradshafranov.py
  :: test_ampere_law()`. Uses a synthetic case with a pure `l=1`-mode
  source (not the full EOS/SCF), computes `oint B·dl` on a rectangular
  loop in `(r,θ)` and compares it against `−∫∫ [source/sinθ] dr dθ` (an
  equivalent form derived from Ampère's law for this geometry). Closes to
  ~2%.
- **Flux**: *proposed* during debugging (see the conversation history that
  fixed the bug) but **not implemented as a separate test** — the Ampère
  test alone already revealed and confirmed the bug fix, so the flux
  consistency test was never actually written. Recorded as a gap in §7.

---

### 1.5 Bernoulli and the parametrization

```
H + Φ − M(u) = C ,        M(u) = ∫ f(u) du
```

With `f = k₀` constant, `M(u) = k₀ u`.

The integration constant is fixed **at the center**:

```
C = H(ρc) + Φc − M(u_c)
```

where `Φc`, `u_c` here denote the values at the center (`r=0`) —
**note**: this use of `u_c` (value of `u` at the center) is different
from the `u_c` used in §1.9 (value of `u` on the last closed line, at the
surface). They are the same symbol for two different evaluation points;
the code never needs both at once, but the reader should note the
collision. Since `ϖ = r sinθ` vanishes at the center for any `θ`,
`u_c(center) = 0` always and the term `M(u_c)` disappears — it stays in
the formula for generality, not because it contributes.

→ `scf/scf.py :: hachisu_scf()` — line `C = H_c + Phi[0, 0] - M_u[0, 0]`.

> **Why it's like this — the (ρc, k₀) parametrization (G4).**
> The classical Hachisu recipe fixes `C` by imposing `H = 0` at two
> surface points (pole and equator). This is ill-posed for this EOS:
> since `n → 3` in the core (§1.1), the mass becomes nearly independent
> of the central density and the radius plunges, so solving for `C` from
> the radius inverts a nearly-singular map. Tested and confirmed in this
> project: the Picard iteration under that recipe is linearly unstable
> across the whole range `ρc = 10⁶`–`10¹⁰` g/cm³, with or without
> sub-relaxation.
>
> The adopted parametrization is by **(ρc, k₀)**, two independent inputs,
> with no surface condition — `C` is fixed by a local condition (`H` at
> the center), not a global one (an integral over the whole mass). This
> change alone solved the instability: geometric convergence to machine
> precision in ~12 iterations, versus exponential divergence of the old
> recipe.
>
> **Caveat:** in the strong-field regime, when the density peak migrates
> away from the center, this anchor also loses the physical meaning that
> makes it stable — see §6.

---

### 1.6 The SCF loop

```
1.  ρ ← initial guess (spherical n = 3 polytrope)
2.  A_φ ← 0
3.  repeat:
4.      Φ   ← Poisson(ρ)
5.      A_φ ← GradShafranov(ρ, f)
6.      u   ← ϖ A_φ ;  M ← ∫ f du
7.      C   ← H(ρc) + Φc − M(u_c)
8.      H   ← C − Φ + M(u)
9.      ρ_new ← EOS⁻¹(H)          [H ≤ 0 => ρ = 0]
10.     ρ   ← (1−ω) ρ + ω ρ_new        ω ≈ 0.3
11. until max|Δρ|/ρc < tol
```

→ `scf/scf.py :: hachisu_scf()` implements steps 1, 3–9 and 11 literally
(step 2 is implicit: `u` starts at zero before the first iteration).

**Discrepancy with the real code, step 10:** the implementation **does
not do sub-relaxation**. The corresponding line in `hachisu_scf()` is
`rho = rho_new` — a direct replacement, equivalent to `ω = 1`, not
`ω ≈ 0.3`. There is no `ω` parameter in the function signature. This is
not an oversight: it is a direct consequence of the change described in
the §1.5 box. The old recipe (two surface conditions) required
sub-relaxation to try to stabilize an iteration that was unstable
regardless; the `(ρc, k₀)` parametrization converges by direct
replacement because the instability that sub-relaxation was trying to
mask was removed at the root. See the project note at the top of
`scf/scf.py` for the full history (including tests showing that
sub-relaxing the old recipe did not fix the problem).

The convergence criterion (`tol`) is a numerical, not physical, parameter
— the dashboard exposes values from `1e-4` to `1e-8` (Tab 1), with `1e-6`
as the default.

---

### 1.7 Virial and energy diagnostics

```
W = ½ ∫ ρ Φ dV                      (gravitational, negative)
Π = ∫ P dV                          (internal)
E_pol = ∫ (B_r² + B_θ²) / 8π  dV
E_tor = ∫ B_φ² / 8π  dV
ℳ = E_pol + E_tor
```

→ `scf/diagnostics.py :: gravitational_energy()` (W), `pressure_integral()`
(Π), `magnetic_energies()` (E_pol, E_tor, ℳ — called `E_mag` in the
code). All use `volume_integral()` as their base (`dV = r² sinθ dr dθ dφ`,
`φ` integrated analytically to `2π`).

Scalar virial theorem for a static configuration:

```
W + 3Π + ℳ = 0
```

Virial error, used as a quality gate:

```
VE = | W + 3Π + ℳ | / |W|          acceptance criterion: VE < 10⁻³
```

→ `scf/diagnostics.py :: virial_error()`. The `10⁻³` threshold is a
project convention (V3 of the plan), not a physical constant — it appears
hard-coded as a comparison in `dashboard/pages/1_equilibrium.py` and
`3_export.py`, not in `units.py` nor in `diagnostics.py` (no physics
module defines this threshold as a named constant — it is a UI/gate
decision repeated in two pages).

**Magnetic virial identity.** Exact, and it is the test that links the
magnetic sector to the gravitational one:

```
∫ ρ ∇M(u) · r  dV  =  ∫ B²/8π  dV
```

This identity **does not correspond to any function in the code**. It was
verified numerically in an ad hoc way during the debugging of the §1.4 bug
(comparing `k0 * ∫ρ r ∂u/∂r dV` with `∫B²/8π dV` for a converged
solution), but there is no `scf/diagnostics.py` function that computes it
nor a committed test that exercises it. Recorded in §7 as a gap.

> **Mandatory conceptual note.** `M(u)` **is not** the magnetic energy. It
> is the specific potential of the Lorentz force — only the part that
> enters the Bernoulli equation. There is no *local* identity between
> `M(u)` and `B²/8π`. The *global* identity above is what forces the two
> quantities to agree in order of magnitude. Measured in this project, in
> the linear regime (perturbative field), the ratio
> `(E_mag/|W|) / (M_u/H_c)` is ≈ 0.5 and is constant in `k₀` (confirmed by
> doubling `k₀`: both ratios quadruple, the ratio between them does not
> change); it drifts to ~0.38 as the field stops being a small
> perturbation (`k₀ ≳ 10⁻¹²` in the `ρc=10⁹`, `R≈3×10⁸` cm
> configuration). This was the observation that, together with the
> Ampère test, confirmed that the bug in the §1.4 box was real and not a
> units artifact.

**Physical limit.** `ℳ < |W|` is a rigid constraint: `E_mag/|W| ≥ 1` is
an impossible configuration (it would violate the virial with `Π ≥ 0`).
Use it as a sanity anchor when reading any result — if the dashboard
shows `E_mag/|W|` above a few tenths, be suspicious before believing it.

---

### 1.8 The two Bt/Bp ratios

```
Bt/Bp (energy)     = E_tor / E_pol
Bt/Bp (amplitude)  = max|B_φ| / max|B_pol|
```

→ `scf/toroidal.py :: bt_bp_ratios()`

They differ by orders of magnitude because the toroidal field is confined
to a small volume (§1.9). The literature frequently does not say which
one it uses. **The dashboard always shows both, labeled** (Tab 1 and
Tab 3).

---

### 1.9 The toroidal field and the twisted torus

In the barotropic formulation, `β = β(u)`: the toroidal function depends
only on the flux function. It is this condition that cancels the φ
component of the Lorentz force and allows the Bernoulli integral to exist
(§1.5).

Geometric consequence: outside the star there is no current, so `B_φ = 0`
there; since `β = β(u)`, it must vanish on every flux line that escapes
the surface. **The toroidal field is automatically confined to the region
of closed poloidal lines** — the twisted torus is not imposed by
geometric decree, it falls out of the consistency of the equation.

Adopted functional form:

```
β(u) = ζ (u − u_c)^{m+1} Θ(u − u_c) ,     m ≥ 1
```

with `u_c` the value of `u` on the last closed line (here, `u_c` = the
value of `u` at the surface — see the notation distinction in §1.5). The
exponent `≥ 1` keeps `ββ′` continuous at the edge of the torus.

→ `scf/toroidal.py :: impose_toroidal()` implements this form **already
substituted into `B_φ = β/ϖ`**: `B_φ = ζ(u−u_c)^{m+1}/ϖ` for `u > u_c`,
`0` outside. `find_uc()` implements the search for `u_c` as the
**maximum of `u` along the entire stellar surface** (`H=0`, sweeping `θ`
from pole to equator) — it is this specific choice of "last closed line"
that the code uses; contours with `u` larger than that maximum, by
definition, do not touch the surface anywhere and remain entirely
interior.

`ζ` is not a direct user input parameter — the dashboard asks for the
target `Bt/Bp` (energy) ratio and solves for `ζ` to reach it, exploiting
the fact that `B_φ` is linear in `ζ` (so `E_tor` is quadratic):

→ `scf/toroidal.py :: solve_zeta_for_energy_ratio()`

> **Why the toroidal field is imposed rather than solved for (G4).**
> The condition `β = β(u)` confines the toroidal field to a small-volume
> torus, and the resulting energy ratio — if `β` were extracted from a
> general barotropic closure condition instead of chosen freely — comes
> out to a few percent. **It is not possible to reach Bt/Bp ~ 1/2 that
> way.** That's why the project solves the barotropic SCF only for the
> poloidal field (§1.4–1.6, with `f(u)=k₀`) and imposes the toroidal
> field on top, at the desired ratio (`ζ` solved for the target), loading
> it into Castro out of exact equilibrium and relaxing it with damping
> (sponge, Tab 3). For a dynamical study, exact equilibrium is not
> required — it just needs to be close enough that the transient doesn't
> destroy the topology.

**Torus volume fraction** (how much of the star has `u > u_c`):
→ `scf/toroidal.py :: closed_torus_volume_fraction()`

**Torus radial extent** (used to check mesh resolution, Tab 3):
→ `scf/toroidal.py :: torus_radial_extent()`, measured along the equator
by default.

---

### 1.10 Time scales and units

```
t_dyn    = √( R_eq³ / GM )
v_A      = ⟨B⟩ / √(4π ρ̄)
t_Alfven = R_eq / v_A
```

→ `dashboard/units.py :: dynamical_time()`, `alfven_speed()`,
`alfven_time()`. `⟨B⟩` and `ρ̄` (volume averages) are computed in
`dashboard/pages/3_export.py` directly with `diagnostics.volume_integral()`
— there is no dedicated `mean_field()` function in `scf/`.

The ratio `t_Alfven / t_dyn` measures the cost of the simulation: in the
weak-field regime it reaches 10³–10⁴, which is prohibitive; in the
strong-field regime of this project, with `E_mag/|W| ~ 0.1`, it stays at
2–3 (D4 of the plan).

**Natural field unit.** From `B²/8π ~ Gρ²R²`:

```
B_unit = R ρ √(8πG)
```

For `ρc ~ 10⁹` and `R ~ 10⁸` cm this gives `~10¹⁴` G, which is also the
order of a white dwarf's virial field. **This formula is not implemented
anywhere in the code** — it was used only as an order-of-magnitude
estimate during the debugging of the §1.4 bug, to check whether the
reported field made physical sense before looking for the numerical
error. **Sanity anchor:** with `E_mag/|W| ~ 0.1`, the field shown by the
dashboard should be in the `10¹³`–`10¹⁴` G range; if it comes out very
different from that, be suspicious before believing it (it was exactly
this suspicion that found the bug).

**Castro convention:** the field is loaded as `B′ = B / √(4π)`. The
dashboard **always displays gauss** and converts only at export time
(Tab 3).

→ `dashboard/units.py :: gauss_to_castro()`, `castro_to_gauss()`. Checked
in this work cycle: the factor matches the plan (`√(4π) ≈ 3.5449`).
**Status note:** `gauss_to_castro()` exists and is correct, but **it is
not called anywhere in the current export pipeline**
(`dashboard/pages/3_export.py` writes `B_phi` to the HDF5 file directly
in gauss, with an explicit `units` attribute documenting this). Converting
to the Castro `B′` convention is the responsibility of Castro's
`problem_initialize.H` (not yet written — Phase 0 of the plan is
pending), which will read the HDF5 file in gauss and apply
`gauss_to_castro()` (or the C++ equivalent) when assembling the internal
state. See §7.

All display-unit conversions (gauss, km) — both the number and the string
formatting — live in `dashboard/units.py`, the single source of truth
(dashboard rule R4): `cm_to_km()`, `g_to_msun()`, `format_gauss()`,
`format_km()`, `format_km_value()`.

---

### 1.11 Rotation and the self-consistent toroidal branch

Extension added after the original poloidal-only spine, motivated by the
collaboration's interest in ~2 M☉ white dwarfs as SN Ia progenitors:
differential rotation plus a toroidal field, because the two deformations
oppose each other (rotation flattens, toroidal elongates), which pushes
back the mass-shedding limit.

**Architecture.** `scf.hachisu_scf()` takes three optional *terms*
(`scf/terms/`) instead of one function per physics combination — each
term contributes an additive piece to the same Bernoulli equation:

```
H + Φ − C_rot(ϖ) − M_pol(u) + M_tor(ρϖ²) = C
```

rearranged for the loop: `H = C − Φ + C_rot(ϖ) + M_pol(u) − M_tor(ρϖ²)`.
`rotation`/`poloidal`/`toroidal` = `None` turns a term off (zero
contribution, zero energy). `poloidal` and `toroidal` are mutually
exclusive (out of scope: self-consistent mixed field — same reason as D6,
barotropy does not deliver a free Bt/Bp). With all three `None`, the loop
reproduces the pre-extension code **bit for bit** (checked against a
frozen copy of the old `scf.py`, `scf/tests/test_regression_v0.py`).

**Rotation** (`scf/terms/rotation.py`): j-constant law,
`Ω(ϖ) = Ω_c A²/(A²+ϖ²)`, with rigid rotation as the `A→∞` limit, handled
analytically (not a large-but-finite `A` substitute). Its Bernoulli
potential has a closed form, `C_rot(ϖ) = (Ω_c²A²/2) ϖ²/(A²+ϖ²)`, reducing
to `½Ω_c²ϖ²` as `A→∞`.

The virial theorem gains the rotational kinetic term:

```
2T + W + 3Π + ℳ = 0 ,      T = ½ ∫ ρ Ω²(ϖ) ϖ² dV
VE = |2T + W + 3Π + ℳ| / |W|
```

→ `diagnostics.virial_error()` (now takes `T=0.0`, default reduces
exactly to the old formula) and `diagnostics.virial_error_terms()` (the
term-aware wrapper — pulls `T`/`Br`/`Bth`/`Bphi` from whichever terms are
active and calls `virial_error()`, so the residual formula itself lives
in exactly one place).

**Stability/termination diagnostics:**
- `T/|W|`: ≥0.14 secular non-axisymmetric instability threshold
  (Ostriker & Bodenheimer 1973), ≥0.27 dynamical bar-mode threshold. Same
  traffic-light convention as VE; blocks export (Tab 3) at ≥0.14, same as
  VE's R5.
- Equatorial mass-loss ratio: `Ω²(R_eq)R_eq / (dΦ/dr at R_eq, equator)` —
  → 1 signals the effective gravity vanishing at the equator (Keplerian
  breakup). → `diagnostics.equatorial_mass_loss_ratio()`.

**Self-consistent purely toroidal field** (`scf/terms/toroidal_sc.py`):
`B_φ = K ρ^m ϖ^{2m-1}`, `m≥1`. Derived from the Lorentz force (Maxwell
stress tensor, verified symbolically with SymPy for `m=1,2,3/2` —
residual exactly 0; `scf/tests/test_toroidal_sc.py`):

```
M_tor(s) = [m K² / (4π(2m−1))] s^{2m−1} ,    s = ρϖ²
```

Unlike `M_pol(u)`, `M_tor` depends on `ρ` — the *unknown* being solved for
at each grid point — so step 9 of the SCF loop (§1.6) stops being a direct
EOS inversion and becomes a per-point root find (`H(ρ)+M_tor(ρϖ²)=$RHS$`,
monotonic increasing in `ρ`, safe to bracket with `scipy.optimize.brentq`
— → `scf.py :: _solve_rho_implicit()`). This branch does **not** use
Grad-Shafranov at all — no flux function, no `Δ*`, no Green's function;
`B_φ` is an algebraic function of the local `ρ`, not a PDE solution.

Validated: turning on the field at fixed `ρc` increases the mass (sign
check); pure toroidal field produces **prolate** deformation
(`R_pol > R_eq`, opposite of the poloidal/oblate case), VE closes
(`~4×10⁻⁴`).

> **Correction — the sign of the magnetic virial identity, and why it
> is not "M_tor depends on ρ" (G4/G1).**
> The magnetic virial identity for this branch is
> ```
> ∫ ρ ∇M_tor · r  dV  =  − ∫ B_φ²/8π  dV
> ```
> a **minus** sign — confirmed by an independent derivation (Maxwell
> stress tensor divergence: `∫r·f_L dV = +∫B²/8π dV` in general, for any
> purely toroidal field, regardless of any Bernoulli convention) and by
> clean numerical convergence with grid refinement (2.27%→1.00%→0.36%→
> 0.14% at nr=65→257).
>
> The first version of this note attributed the sign to "`M_tor` depends
> on `ρ`, `M_pol` does not." **That attribution was wrong**, caught on
> review. Here is the actual mechanism: writing the master Bernoulli as
> `H + Φ − C_rot − M_pol + M_tor = C` (note the signs: `M_pol` enters with
> a minus, `M_tor` with a plus) and comparing its gradient against the
> momentum equation gives `F_L,tor/ρ = −∇M_tor` — the minus comes directly
> from `M_tor` having been written with a **plus** sign in the master
> equation, full stop. Substituting into the sign-fixed, convention-free
> Maxwell identity `∫r·F_L dV = +∫B²/8π dV` is what produces the minus
> sign in the virial identity. Had the master equation instead been
> written with `−M_tor`, the virial identity would come out `+`, with
> `M_tor` depending on `ρ` in exactly the same way either time. **`ρ`
> entering `M_tor` is why the SCF inversion becomes implicit (a fact about
> the algorithm's structure) — it has nothing to do with which sign the
> virial identity carries (a fact about which convention was chosen when
> the term was first written into the Bernoulli equation).** Two different
> questions, easy to conflate, and conflating them is exactly the kind of
> mistake this project has already been burned by once (§1.4's Green's-
> function bug survived because independent-looking checks secretly
> shared one derivation). Get the causal story right here or a future
> reader "fixing" the (correct) code to match the (wrong) causal
> explanation is a real, specific risk — see `scf/terms/toroidal_sc.py`
> for the full derivation.

**V-R1/V-R2 status (rigid vs. differential rotation, no field).** V-R2
(differential, j-constant, `A=0.3 R_guess`) reaches `2.19 M☉` against the
`~2.2 M☉` target (Yoon & Langer 2005) — **0.40% error**, reached with
`mass_loss_ratio~0.13`, comfortably away from breakup.
`scf/tests/test_differential_rotation.py`.

V-R1 (rigid) **does not validate** against the `~1.5 M☉` target (Hachisu
et al. 2012; Boshkayev et al. 2013). At `ρc=10¹²`, stepping `Ω_c` upward
with continuation, the sequence terminates *numerically* at
`R_pol/R_eq≈0.932`, `Ω_c≈26.47`, with **mass_loss_ratio only 0.135** —
far from the `→1` that a genuine mass-shedding termination requires, and
far from the Roche-model expectation for a centrally-condensed (`n≈3`)
configuration, `R_pol/R_eq=2/3≈0.667` at real breakup. The cause is not
resolved. Two suspects, not yet investigated (not on the project's
critical path — Jorge wants differential rotation, which already works):
outer-layer radial resolution, and — most likely — `Ω_c` itself being a
bad control parameter this close to the sequence's end (a documented
phenomenon in the rotating-SCF literature; Hachisu's own method sidesteps
it by parametrizing on axis ratio and solving for `Ω²`, not the reverse).
A direct test confirms this: bisecting for a target axis ratio instead of
imposing `Ω_c` directly saturates at the same point — no smaller axis
ratio is reachable by continuation in `Ω_c`, for any target. The
mechanism itself is not in doubt: `(M−M₀)/M₀` tracks `T/|W|` with a
stable proportionality coefficient (~3.0) across the whole tested range,
and V-R6 (below) independently confirms `T`. See
`scf/tests/test_rotation.py` for the full numbers and reasoning chain.

> **Bug found and fixed during this investigation (unrelated to the V-R1
> question above, but real and previously unnoticed).**
> `diagnostics.surface_radius()` (and everything built on it —
> `equatorial_polar_radii()`, `surface_field()`, `surface_dipolarity()`,
> `toroidal.find_uc()`, `impose_toroidal()`) used to take `ρ` and
> linear-interpolate where it crosses zero. But `ρ = EOS⁻¹(H)` is clipped
> to *exactly* `0.0` for `H≤0` (`eos.density_of_enthalpy`) — so the grid
> point just past the surface always has `ρ=0.0` exactly, which makes the
> interpolation fraction collapse to exactly `1.0` and the function always
> return a raw grid point, never a true sub-grid position. Radii were
> silently grid-quantized project-wide. Fixed by interpolating on `H`
> instead (continuous, crosses zero smoothly between grid points, never
> clipped) — verified to shift `R_eq`/`R_pol` by tens of km on a modest
> mesh (`nr=65`), not a cosmetic correction. Found by noticing that a
> rotation sequence's `R_pol/R_eq` was landing on suspiciously exact small
> integer ratios (`37/39`, `38/39`, ...) — the classic signature of grid
> quantization, not physics.

**T (Fase 1, V-R6).** Computed two independent ways — once via
`rotation.energy()`, once via a from-scratch volume integral written in
the test, not calling any shared code beyond `Ω(ϖ)` — and the two agree
to machine precision (`rel_diff=0.00e+00`).

---

## 2. Tab 1 — Equilibrium

→ `dashboard/pages/1_equilibrium.py`

Runs a single SCF pass (§1.6) and shows the result. Uses the entire
theoretical core (§1.1–1.10).

**Input parameters** (sidebar): `ρc`, `μₑ`, `k₀` (optional, toggles the
poloidal field on/off), target `Bt/Bp` ratio and `m` (toroidal, optional),
plus the numerical mesh parameters (`Nr`, `Ntheta`, `l_max`, `tol`,
`max_iter`). `θ` always covers the full `[0, π]` — no equatorial symmetry,
a plan decision (D3) so as not to mask asymmetric modes (`m=1`).

**`k₀` range.** Not known a priori (it depends on `ρc`, `R`, `μₑ` in a
non-trivial way — see the §1.4 bug, which for a long time made the
apparent range look wrong by 3–4 orders of magnitude). The dashboard
probes empirically by raising `k₀` geometrically until `VE > 10⁻³` or the
SCF stops converging, on a coarse mesh, and stores the result in a cache
at `dashboard/k0_range_cache.json`, indexed by `(ρc, μₑ)`.
→ `dashboard/pages/1_equilibrium.py :: _estimate_k0_max()`

**Displayed scalars** ("Scalars" table):

| Scalar on screen | Definition | Function |
|---|---|---|
| `M/M_sun` | total mass / `M_SUN` | `scf.total_mass()`, `units.M_SUN` |
| `R_eq`, `R_pol` (km) | radius where `ρ→0`, at the equator and the pole | `diagnostics.equatorial_polar_radii()` |
| `R_pol/R_eq` | flattening | direct ratio |
| `rho_c confirmed` | `ρ[r=0]` post-convergence | should match the input `ρc` |
| `mean rho` | `M / volume of the ellipsoid (R_eq, R_eq, R_pol)` | geometric approximation, not an exact integral |
| `W`, `E_int=Π`, `E_mag`, `E_pol`, `E_tor` | §1.7 | `diagnostics.virial_error()`, `magnetic_energies()` |
| `E_mag/|W|` | dimensionless field strength | direct ratio (see the sanity anchor in §1.7/1.10) |
| `B_pol,max`, `B_central`, `B_tor,max` | field values, gauss | maxima/pointwise values of `Br`,`Bth`,`Bphi` |
| `torus volume fraction` | §1.9 | `toroidal.closed_torus_volume_fraction()` |
| `VE` | §1.7 | `diagnostics.virial_error()` |

**Convergence**: two plots (`max|Δρ|/ρc` and `VE`, both vs. iteration,
log scale) → `plots.plot_convergence()`, `plots.plot_virial_history()`.
The `VE` history only exists if `hachisu_scf(..., track_virial=True)` —
it costs extra integrals every iteration, so it is optional (see §1.6).

**Bt/Bp**: the two ratios of §1.8, side by side, always labeled.

**Meridional-plane figures** (§1.3):
- **density**: color map of `ρ`, with the `H=0` boundary (the stellar
  surface) overlaid in cyan → `plots.plot_density()`
- **poloidal field lines**: contours of `u` — physically, each contour
  *is* a poloidal field line (§1.3), because `B` is tangent to the curves
  of constant `u` by construction (`B·∇u = 0` follows directly from the
  formulas for `B_r`,`B_θ` in terms of derivatives of `u`). The last
  closed line (`u_c`, §1.9) is highlighted in red when a toroidal field is
  imposed → `plots.plot_flux_contours()`
- **toroidal field**: color map of `B_φ`, shows the confined torus →
  `plots.plot_toroidal()`

All radial axes in km, field always in gauss with a colorbar in
scientific notation (rule R4, see §1.10).

---

## 3. Tab 2 — Sweep

→ `dashboard/pages/2_sweep.py`, `dashboard/sweep_worker.py`

Runs the SCF (§1.6) on a `(ρc, k₀)` grid in parallel
(`ProcessPoolExecutor`), with caching by parameter hash
(`store.run_exists()`/`store.save_run()`, §5). Each grid point is an
independent call to `scf.hachisu_scf()` followed by the same diagnostics
as Tab 1 — no new physics, just orchestration.

→ `dashboard/sweep_worker.py :: run_one()` — the picklable function each
grid worker process runs.

**What an equilibrium sequence is.** Fixing `ρc` and varying `k₀` (or vice
versa) produces a sequence of equilibrium configurations. **The sequence
ends** when the SCF stops converging or when `VE` exceeds the acceptance
threshold (§1.7) and does not improve with resolution — see the measured
numbers in §6. A sequence ending is, in itself, a physical result (not an
error to be fixed): it signals the limit of validity of that family of
equilibria.

**M-R diagram.** `R_eq` (km) on the x-axis, `M/M_sun` on the y-axis,
colored by `k₀`. The horizontal line at `1.44 M_sun` marks the
Chandrasekhar limit **without a field** (`μₑ=2`) — it is not the physical
limit of the magnetized sequence, it's just a reading reference. Optional
overlay of literature data points from
`dashboard/data/references/bera_bhattacharya_2014.csv`, if the file
exists (it does not currently exist in this repository — no data has been
digitized; the dashboard works without it and warns that it's missing).

**Why the relevant mass maximum is over the whole plane, not a slice.**
This is the point that most often causes an incorrect comparison with the
literature. `M_max(k₀=0)` already depends on `ρc`: it only approaches the
Chandrasekhar limit asymptotically, for high `ρc` (`~10¹¹`–`10¹²` g/cm³ —
see `scf/tests/test_scf_v1.py`, validated to 0.78% in that range). A slice
at low `ρc` (for example `ρc=10⁹`, used during this project's debugging)
gives `M(k₀=0) = 1.39 M_sun` — **already below** the field-free
Chandrasekhar limit, and any mass measured along that slice (with or
without a field) is not comparable to the literature reference number
(`M_max ~ 1.9 M_sun` in Bera & Bhattacharya 2014), which is the maximum
taken over the **entire** `(ρc, k₀)` plane, not over a line with `k₀`
varying at fixed `ρc`. This tab's sweep is exactly what allows that
comparison to be made correctly — by mapping the plane, not a slice.

**Grid convergence.** Points that fail to converge are recorded (not
silently dropped) and shown in a separate expander; the grid's convergence
rate is, itself, information about where the family of solutions ends.

**VE heat map.** Over the `(ρc, k₀)` grid, in `log₁₀(VE)` — reveals where
the method is at the edge of validity without having to read numbers one
by one.

---

## 4. Tab 3 — Export

→ `dashboard/pages/3_export.py`

Takes an already-saved equilibrium (Tab 1 or 2), imposes the toroidal
field (§1.9) and generates the three output artifacts: an HDF5 initial
data file, Castro's `inputs`, and `run_manifest.json`.

**Toroidal imposition.** Same function as §1.9
(`toroidal.solve_zeta_for_energy_ratio()`), with this tab's own controls
(the target ratio may differ from the one used when the equilibrium was
saved). Continuity of `ββ′` at the torus edge is **guaranteed
analytically** for `m ≥ 1` (it is not recomputed numerically here — it is
a property of the functional form of §1.9, `(u-u_c)^{m+1}` and its
derivative go to zero at `u=u_c` for any `m≥1`); the page only confirms
that `m_tor≥1` was respected.

**Why the field is initialized from the vector potential, not from B at
cell centers.** This is a Castro-side decision (D3 of the plan, section
5), not something the SCF resolves — but the dashboard already prepares
the data in that format: the exported HDF5 file contains `A_φ` (computed
as `u/ϖ` from the converged flux function), **not** `B_r`, `B_θ`
directly. Castro's constrained-transport algorithm keeps `B` on the mesh
faces; initializing via `∇×A` interpolated onto the edges guarantees
`∇·B = 0` to machine precision by construction. Initializing with `B`
values at cell centers does not give that guarantee. **What's missing:**
interpolating `A_φ` onto the Cartesian mesh edges and computing `∇×A`
there is the responsibility of Castro's `problem_initialize_mhd_data.H`,
which has not been written yet (Phase 0 of the plan is pending) — the
dashboard delivers `A_φ` on a spherical `(r,θ)` mesh; interpolation onto
Castro's Cartesian mesh happens on the other side.

**B′ = B/√(4π) convention.** See §1.10. The HDF5 file exported by this
tab stores `B_φ` in pure gauss, with a `units` attribute in the header
saying so explicitly — the conversion to the Castro convention does not
yet happen in this pipeline (see the status note in §1.10 and the gap in
§7).

**Derived scales and simulation cost.** `t_dyn`, `v_A`, `t_Alfven` and the
`t_Alfven/t_dyn` ratio from §1.10, computed from the `⟨B⟩` and `ρ̄` of the
loaded equilibrium. A success badge appears when the ratio is between
0.3 and 3 — the strong-field regime that makes the simulation cheap (D4
of the plan).

**Torus resolution check.** Given the number of cells per side of
Castro's box (128/256/384) and the box size (a multiple of `R_eq`),
computes how many cells cross the torus's radial extent
(`toroidal.torus_radial_extent()`, measured at the equator). Fewer than
10 cells triggers a warning — the torus is small compared to the star and
it is the torus that needs to be resolved (plan, section 6).

**Box parameters.** `castro.small_dens` is chosen as `ρc × 10⁻⁸` — this is
a **convention of this dashboard's design**, not a value derived from
physics nor cited in the plan; there is no documented theoretical
justification for the `10⁻⁸` exponent beyond it being a density floor low
enough not to perturb the star and high enough not to generate an absurd
`v_A` in the numerical vacuum (the "fluff problem", plan section 5). The
sponge start radius (`1.5 R_eq`), the damping duration (`5 t_dyn`) and the
perturbation amplitude (`10⁻⁴ c_s`, using `eos.sound_speed()` at the
center) come directly from the values cited in the plan (section 5), and
are not recalculated/optimized by this dashboard.

**R5 — export blocked if `VE ≥ 10⁻³`**, with no override option.
Implemented as a simple `if`/`else` around the export button — there is no
override mechanism at any layer.

---

## 5. Tab 4 — Runs

→ `dashboard/pages/4_runs.py`, `dashboard/store.py`

No physics — this module only persists what the other tabs have already
computed (dashboard rule R1/R3).

**What gets saved**, per run, in `dashboard/runs/<hash>/`:

| File | Content | Function |
|---|---|---|
| `params.json` | full input parameters | `store.save_run()` |
| `scalars.json` | derived scalars (same as the Tab 1 table) | same |
| `fields.npz` | `rho, Phi, u, H, Bphi, r, theta` on the mesh | same |
| `manifest.json` | hash, timestamp, git, dependencies | same |

`dashboard/runs/index.csv` aggregates hash + parameters + scalars for all
runs, so Tab 4 can load quickly without opening every directory
(`store.load_index()`).

**Why the git hash matters.** `manifest.json` records
`git_commit_hash(REPO_ROOT)` — the commit of `scf/` **and** `dashboard/`
(the same git repository covers both, initialized specifically to give
this dashboard provenance; `amrex/`, `castro/`, `microphysics/` are
excluded via `.gitignore`, being their own repositories). A scalar without
its associated commit cannot be reproduced with confidence — if the code
changes (for example, the §1.4 bug being fixed), results saved before and
after **are not comparable**, even with the same input parameters.
`git_dirty()` is also recorded — it flags whether there were uncommitted
changes at the time of the run, which makes exact reproduction impossible
even knowing the hash.

→ `dashboard/store.py :: git_commit_hash()`, `git_dirty()`,
`dependency_versions()` (versions of `numpy`, `scipy`, `streamlit`,
`plotly`, `h5py`, `python`).

**Cache.** The hash is `sha256(json(params, sort_keys=True))[:12]`
(`store.params_hash()`) — deterministic in the parameters, used both to
name the run's directory and for Tab 2 to skip already-computed points.

**Schema version:** **there is no** `schema_version` field (or
equivalent) in `manifest.json` or `scalars.json` currently. When the set
of scalars changed during this project (for example, when
`B_pol,max (G)` was added to Tab 2's schema), runs saved before the
change ended up with missing columns in `index.csv`, breaking charts that
expected the new column — this actually happened during development and
was worked around by manually deleting the old cache, not by a
versioning mechanism. Recorded in §7 as a gap.

**Features:** filterable/sortable table (`st.data_editor`, formatted per
column per rule R4 — gauss in scientific notation, km with 2 decimal
places), side-by-side comparison of two runs, reload a run into Tab 1
(via `st.session_state["reload_run_params"]` + `st.switch_page`), mark as
reference (`store.mark_reference()` — today it only records the flag; no
other tab reads `reference` yet, see §7).

---

## 6. Known limitations

Measured numbers, not vague descriptions — all in the
`ρc = 10⁹` g/cm³, `R ≈ 3×10⁸` cm, `nr=161`, `ntheta=65`, `l_max=16`
configuration unless otherwise noted:

- **The `VE < 10⁻³` criterion fails above `k₀ ≈ 2×10⁻¹²`** in this
  configuration. The last valid point is `k₀ ≈ 1.6×10⁻¹²`,
  `M ≈ 1.50 M_sun`, `VE ≈ 6.6×10⁻⁴`. At `k₀ ≈ 2.3×10⁻¹²`
  (`M ≈ 2.02 M_sun`), `VE ≈ 1.57×10⁻³` — above the threshold.
- **At that same point (`k₀ ≈ 2.3×10⁻¹²`) the density peak migrates away
  from the center** — from `r_idx=1` (essentially the origin) to
  `r_idx≈25` (a real radius, not a grid artifact) — and `R_pol/R_eq`
  drops to `~0.61`: genuine equatorial evacuation. The anchor at
  `ρc(r=0)` (§1.5) loses physical meaning right there, because the center
  is no longer the point of maximum density.
- **Whether the `VE` failure above is insufficient resolution or a
  genuine sequence termination is an open question — but there is strong
  evidence for the second hypothesis**: a convergence study done in this
  project (`l_max` from 16 to 48, mesh from `129²` to `385²`, nearly 9×
  more points) kept `VE` at `1.2`–`1.6×10⁻³`, **with no downward trend**
  — a plateau above the threshold, not a curve converging to zero. This
  is the expected signature of a physical termination, not of
  under-resolution.
- **All the sequence results (`k₀` varying) documented above are from a
  single `ρc = 10⁹` slice**, where `M(k₀=0) = 1.39 M_sun` — below the
  field-free Chandrasekhar limit itself (`1.44 M_sun`). Comparisons with
  literature mass maxima (`M_max ~ 1.9 M_sun`, Bera & Bhattacharya 2014)
  require sweeping the entire `(ρc, k₀)` plane (Tab 2), not a slice — see
  §3.
- **Castro's `B′ = B/√(4π)` conversion is not applied anywhere in the
  current export pipeline** (§1.10, §4) — the exported HDF5 file stores
  `B` in pure gauss, documented via an attribute. This is left as the
  responsibility of Castro's `problem_initialize.H` (not yet written).
- **Phase 0 of the plan (building Castro with `USE_MHD=TRUE`) is
  pending** — system dependencies (`gfortran`, `libhdf5-openmpi-dev`,
  `libopenmpi-dev`) have not yet been installed in this environment.
  Nothing in Tab 3 has been tested against a real Castro build.
- **Rigid rotation (V-R1) does not validate against the literature
  (~1.5 M☉)** — the `Ω_c`-controlled sequence terminates numerically at
  `R_pol/R_eq≈0.93` with mass-loss ratio only `0.135` (should be `→1` at
  real breakup; Roche model predicts `R_pol/R_eq=2/3` for this EOS's
  central concentration). Not on the project's critical path (differential
  rotation, V-R2, already reaches the `~2.2 M☉` target at 0.40% error, far
  from breakup) — see §1.11.

---

## 7. Open questions

Gaps identified while writing this document — not filled in on our own
initiative (rule G1):

1. **Flux consistency test** (`u(ϖ,z) = ∫₀^ϖ B_z ϖ′ dϖ′`, §1.4) was
   proposed during debugging of the Green's-function bug but never
   implemented as a separate test. The Ampère test alone was enough to
   find and confirm the bug; the flux test would serve as a second,
   independent line of defense.
2. **Magnetic virial identity** (`∫ρ∇M(u)·r dV = ∫B²/8π dV`, §1.7) has no
   dedicated function in `diagnostics.py` and no committed test — it was
   only verified numerically in an ad hoc way during one debugging
   session.
3. **The formula `J_φ = cρϖf(u)`** (§1.4) does not correspond to any
   function — it is never computed as a named quantity in the code.
4. **`B_unit = Rρ√(8πG)`** (§1.10, natural field unit) is not implemented
   — it was used only as an order-of-magnitude estimate during
   debugging.
5. ~~**No `schema_version`** in `manifest.json`/`scalars.json` (§5)~~ —
   **fixed**: `store.SCHEMA_VERSION` now guards `run_exists()`, so a
   schema mismatch is a cache miss instead of a silent hit with missing
   columns (added when rotation/toroidal_sc extended the scalar set).
6. **`store.mark_reference()`** records the `reference` flag in the
   index, but no other tab (in particular Tab 2, which per the original
   prompt should show reference runs on its charts) reads that flag yet.
7. **`castro.small_dens = ρc × 10⁻⁸`** (§4) is a convention with no
   documented theoretical justification — it works as a starting point,
   not as a calibrated value.
8. **The original question that motivated the whole §1.4 bug
   investigation** — whether the ratio `(E_mag/|W|)/(M_u/H_c) ≈ 0.5`
   (linear regime) has a closed-form derivation, or is just what is
   observed numerically for this EOS and this parameter range — was not
   resolved analytically, only confirmed as self-consistent (constant
   under `k₀ → 2k₀`).
9. **`gauss_to_castro()`/`castro_to_gauss()`** exist and are correct
   (checked in this work cycle), but are not called by any export code
   yet — see the limitation in §6.
10. **Rigid-rotation sequence termination (§1.11, V-R1)** — not resolved.
    The standard remedy (parametrize the near-terminal sequence by axis
    ratio instead of `Ω_c`, solving `Ω_c` by an outer root-find — the
    technique Hachisu's own method uses) was identified and spot-checked
    (a direct bisection on axis ratio saturates at the same point found by
    stepping `Ω_c`, confirming `Ω_c` degenerates as a control parameter
    there) but not implemented as production code — out of scope for now
    since differential rotation (the project's actual target) does not
    hit this problem in the tested range.

---

## 8. References

**SCF method and magnetized equilibria**
- Hachisu, I. 1986, ApJS 61, 479 — the SCF method
- Tomimura, Y. & Eriguchi, Y. 2005, MNRAS — twisted torus, canonical reference
- Lander, S. K. & Jones, D. I. 2009, MNRAS — mixed fields, free functions
- Lander, S. K. & Jones, D. I. 2012 — stability of mixed fields

**Magnetized white dwarfs**
- Das, U. & Mukhopadhyay, B. 2014, MNRAS 445, 3951 — SCF for magnetized WDs
- Bera, P. & Bhattacharya, D. 2014, MNRAS — self-consistent M-R with Lorentz force
- Bera, P. & Bhattacharya, D. 2016, MNRAS 456, 3375 — field geometry
- Bera, P. & Bhattacharya, D. 2017, MNRAS 465, 4026 — perturbation study
- Nityananda, R. & Konar, S. 2014 — critique of super-Chandrasekhar models
- Coelho, J. G. et al. 2014 — virial limits

**Stability**
- Markey, P. & Tayler, R. J. 1973 — the m = 1 instability
- Braithwaite, J. & Spruit, H. C. 2004 — relaxation to a stable configuration
- Braithwaite, J. & Nordlund, Å. 2006

**Codes**
- Castro: https://github.com/AMReX-Astro/Castro
- MHD documentation: https://amrex-astro.github.io/Castro/docs/mhd.html
- XNS: https://www.arcetri.inaf.it/science/ahead/XNS/
