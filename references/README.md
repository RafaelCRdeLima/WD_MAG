# References

Downloaded literature on super-Chandrasekhar magnetized white dwarfs. All four
PDFs are the arXiv preprints; the two MNRAS and one ApJ links are behind
publisher bot protection, so the accepted versions were fetched from arXiv
instead. Journal references are given below and should be the ones cited.

**Banibrata Mukhopadhyay is an author on all four.** They are one research
programme seen at four stages, and reading them together is more informative
than reading any one: the group's case for super-Chandrasekhar white dwarfs
rests on magnetic support, and each paper adds a piece — rotation and GR,
observability, matter anisotropy, formation channel.

---

## 1. Subramanian & Mukhopadhyay (2015)

**GRMHD formulation of highly super-Chandrasekhar rotating magnetised white
dwarfs: stable configurations of non-spherical white dwarfs**
MNRAS 454, 752–765 · [arXiv:1507.01606](https://arxiv.org/abs/1507.01606) · 15 pp.

Axisymmetric stationary equilibria of differentially rotating, magnetized
polytropic compact stars in full general relativity, ideal MHD, computed with a
modified version of the open-source **XNS** code (originally a neutron-star
code; the modification is to make it work for white dwarfs). Sequences are
built by varying one physical quantity at a time, separately for purely
toroidal and purely poloidal field geometries, with uniform and differential
rotation.

Main results, for central density $\simeq 2\times10^{10}$ g cm$^{-3}$ and
fields up to $\sim5\times10^{14}$ G:

- Differential rotation alone gives $1.4$–$1.8\,M_\odot$, equatorial radii
  $1200$–$1450$ km, oblate with $R_{\rm pol}/R_{\rm eq}$ down to $0.6$.
- Toroidal field alone gives $1.4$–$2.3\,M_\odot$. The *surface* stays nearly
  spherical while the *interior* becomes extremely prolate, and the mean
  density drops by nearly a factor of ten — bloated low-density outer regions.
- Toroidal field **plus** differential rotation reaches $1.4$–$3.1\,M_\odot$,
  radii to $\sim3300$ km, peak field $3.6\times10^{14}$ G, surface equatorial
  period $\sim10$ s.
- Mass and radius grow with field and with central $\Omega$ apparently without
  bound, but the ratio of kinetic to gravitational energy and the surface
  angular velocity both peak — which the authors read as the likely onset of
  an instability.
- **Polar hollows**: when the central angular velocity exceeds that of the
  corresponding uniformly rotating star at mass shedding, the density
  isocontours turn concave near the poles. Flatter rotation profiles deform the
  star more but produce smaller hollows; the hollows vanish as the rotation
  approaches uniform.

**Why this one matters most here.** This is the paper our configuration comes
out of. Same rotation law, $j(\Omega) = A^2(\Omega_c - \Omega)$, same regime:
$2$–$3\,M_\odot$, interior field $\sim10^{14}$ G, surface period $1$–$10$ s.
Our star sits inside their family — $2.005\,M_\odot$, $3.2\times10^{13}$ G
toroidal, surface period $4.3$ s — and the $T/|W| = 0.14$–$0.16$ threshold
quoted in our report is theirs.

And this is exactly where our work goes beyond it. **These are equilibria.**
Stability is inferred from energy ratios and from the peaks in
$\Omega_{\rm surf}(B_{\max})$, never tested dynamically. A configuration can
satisfy every equilibrium and energy criterion in this paper and still be torn
apart by a non-axisymmetric magnetic instability in a few Alfvén times — which
is what our 3D runs find happening to the field, though not to the star.

---

## 2. Gupta, Mukhopadhyay & Tout (2020)

**Suppression of luminosity and mass–radius relation of highly magnetized white
dwarfs**
MNRAS 496, 894–902 · [arXiv:2006.02449](https://arxiv.org/abs/2006.02449) · 11 pp.

Newtonian, non-rotating. Solves magnetostatic equilibrium together with photon
diffusion and mass conservation, so that the interface between the degenerate
core and the ideal-gas envelope is found **self-consistently** rather than
assumed — which the authors note is what essentially all earlier work did.
Magnetic opacity replaces Kramers' opacity where appropriate.

- Without a field and at fixed temperature across the interface, the
  Chandrasekhar limit survives up to $L \sim 10^{-2}L_\odot$; photon diffusion
  raises the mass of large-radius white dwarfs.
- With central field $\sim10^{14}$ G, super-Chandrasekhar masses $\sim1.9\,M_\odot$
  appear even at constant interface temperature.
- The headline: **small-radius magnetized white dwarfs stay super-Chandrasekhar
  down to $L \sim 10^{-20}L_\odot$**, while their large-radius counterparts on
  the same mass–radius relation converge back to Chandrasekhar's result as $L$
  falls. So these objects could exist and be undetectable simply by being too
  faint. Gravitational waves (LISA) are offered as the way to find them.

**Relevance to us.** Indirect, and it is the odd one out — this is about
observability, not structure or dynamics. But it is the one paper of the four
that takes the thermal structure seriously, and that touches our most awkward
caveat: our EOS is barotropic (`ztwd`, zero-temperature), so the
$6\times10^{49}$ erg released by the field have no physical destination. This
paper is a reminder that a degenerate core plus an ideal-gas envelope behaves
differently from a single barotrope, and it is a natural reference if we ever
justify moving to a non-barotropic EOS.

---

## 3. Deb, Mukhopadhyay & Weber (2022)

**Anisotropic magnetized white dwarfs: unifying under- and over-luminous
peculiar and standard type Ia supernovae**
ApJ 926, 1 · [arXiv:2112.03938](https://arxiv.org/abs/2112.03938) · 11 pp.

General relativistic, applying to white dwarfs the framework the same authors
built for neutron and quark stars in a companion paper. Two anisotropies act
together: local anisotropy of dense matter (parameter $\kappa$) and the
anisotropy induced by the magnetic field itself, whose **orientation** is
treated as a free choice.

- Static equilibrium is not achieved unless *both* anisotropies are included —
  the magnetic one alone does not suffice.
- Field orientation changes the sign of the effect: a **transverse** field
  raises the stellar mass with increasing field strength, a **radial** field
  lowers it.
- The resulting mass–radius relations span from sub- to super-Chandrasekhar,
  with maximum masses reaching $\sim2.4$–$2.8\,M_\odot$ for
  $B_0 \simeq 3.8\times10^{14}$ G depending on $\kappa$.
- The claimed pay-off is a single framework covering over-luminous *and*
  under-luminous SNe Ia, and a warning against treating white dwarfs associated
  with SNe Ia as standard candles.

**Relevance to us.** Structure only, no dynamics, and the anisotropy
parameterisation is a modelling choice rather than something our ideal-MHD run
contains. The transferable point is that **field geometry is not a detail**: it
decides the sign of the mass shift here, just as it decides the deformation
(prolate vs oblate) in paper 1 and the stability in our simulations. Our field
is toroidal-dominated by a factor $10^7$ in energy, which is the most extreme
corner of that space.

---

## 4. Zuraiq, Kumar, Hackett, Bhattarai, Tout & Mukhopadhyay (2024)

**Simulating super-Chandrasekhar white dwarfs**
Astrophysics and Space Science Proceedings (Springer), conference volume ·
[arXiv:2411.18692](https://arxiv.org/abs/2411.18692) · 12 pp.

A proceedings contribution, so shorter and more of a progress report than the
other three. One-dimensional stellar evolution with the Cambridge **STARS**
code, modified to carry a toroidal magnetic field profile in the hydrostatic
balance plus a cooling treatment. The field is taken to matter only after the
star collapses and the field grows, not during the main sequence, and it is
introduced in stages so the model can relax.

- Formation route followed: a $3\,M_\odot$ ZAMS star evolved through the AGB
  with mass loss enabled, leaving a $0.4$–$0.5\,M_\odot$ carbon–oxygen white
  dwarf, which then **accretes** to grow.
- Because the field modifies hydrostatic balance, there is no single
  Chandrasekhar limit but a **series of mass limits** set by the field, with a
  possible ultimate limit once all instabilities are excluded.
- Accretion onto a magnetized versus a non-magnetized white dwarf ends at
  different limiting masses.

**Relevance to us.** This is the formation-channel paper — the one that asks
where such a star comes from, which the equilibrium papers do not. Worth noting
that it reaches the opposite conclusion in spirit from ours on the point that
matters: it is 1D and spherical, so a toroidal field enters only as a pressure
term and *cannot* go unstable. Its "satisfying underlying stability" is
stability against radial/thermal criteria, not against the $m=1$ modes that
destroy the field in our 3D runs. Useful as background on the accretion
channel, and as a contrast: our configuration is merger-remnant-like, spun far
faster than anything this route produces.

---

## What the four have in common, and where our work sits

All four compute **structure**: equilibrium configurations, mass–radius
relations, evolutionary tracks. None of them evolves the magnetic field in
three dimensions, and none can. Stability, where it is discussed at all, is
assessed through energy ratios, radial criteria, or the appearance of turning
points along a sequence.

That is precisely the gap our simulations address. The question these papers
leave open is not whether a $2\,M_\odot$ magnetized white dwarf can be
constructed — it can, repeatedly, by several methods — but whether it survives
once you let it evolve in 3D. Our answer so far, at two resolutions:

- **The star survives.** Mass conserved to $0.02\%$ over $60$ s, and the
  differential rotation that supports it is never erased.
- **The ordered field does not.** It is disrupted by an $m=1$ instability
  within a few Alfvén times, and the magnetic energy falls by three orders of
  magnitude.

If that holds up, it bears directly on every one of these papers: the
configurations are constructible and dynamically viable as *stars*, but the
$10^{14}$ G ordered interior field that does the supporting in the models is
not something a 3D star will hold onto. What our runs cannot yet say is how
much of the decay is physical rather than numerical resistivity — see the
convergence section of `reports/report_rot192_rot256.pdf`.
