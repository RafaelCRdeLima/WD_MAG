"""Campaign HZ against its own control: the heating is not the field's.

Run:  scf/.venv/bin/python3 investigations/plot_hz_control.py

Reads thermal_hz192.csv (192^3 under helmholtz, WITH the field, to t = 1.27 s),
thermal_hz192ctl.csv (the same star with problem.field_scale = 0, to t = 12 s)
and field_hz192.csv (E_mag over the with-field run).

WHAT THIS SETTLES

The 256^3 helmholtz window found a shell at 3.4e9 K and DIARIO 10.1 left two
readings of it open: (a) the field settling out of magnetohydrostatic
equilibrium, which ztwd discarded and helmholtz keeps, or (b) numerics. The
control differs from the run in ONE line -- field_scale = 0 -- so the star it
evolves is in equilibrium by construction: gas pressure, gravity and rotation
is exactly what the SCF balanced.

It heats identically. Over the 1.27 s the two share, E_ion agrees to 1-2% at
every output, and where it differs the CONTROL IS HOTTER. So the answer is
(b), and campaign HZ has answered its own question in the negative.

Two independent reasons, either of which is sufficient:

  1. ENERGY. E_ion rises by 5.09e50 erg with the field and 5.15e50 without.
     The field's entire budget is 6.06e49 erg, so the rise is 8.4x more energy
     than there is field to pay for it. E_ion is a lower bound on the internal
     energy -- at 4e9 K the radiation term is comparable -- so the true excess
     is larger.

  2. THE CONTROL. There is no field in it at all.

The panels:

  (a) rho_max. Both fall by a factor of 140 in 1.27 s and the control, followed
      to t = 12 s, reaches 1.07e6 -- a factor of 2800. This is not a star
      settling, it is a star being destroyed, and it happens with no field
      present.

  (b) E_ion against the field's budget, drawn as a horizontal line. The
      measurement is above the line by t = 0.1 s and above 8x it by t = 0.56 s.

  (c) Mass fraction above 2e9 K. DIARIO 10.1 argued the heat sat in a shell of
      about 5% of the mass, which the field could afford. It does not: 94% of
      the mass passes 2e9 K. The shell reading came from one snapshot at
      t = 0.36 s, before the front had crossed the star.

  (d) E_mag in the with-field run. It falls from 6.05e49 to 1.35e49 -- the
      field really does lose 4.7e49 erg over the window -- and that release is
      irrelevant, because the control heats the same without it.

WHAT IT DOES NOT SETTLE

Why the star is destroyed. The energy is either real, in which case it is
gravitational and the ztwd profile is not an equilibrium under helmholtz, or
the temperature field is corrupt and the star is being driven by a bad
pressure. Both are (b) and neither is the field. The cheap discriminator is
not a run: it is comparing P_helmholtz(rho, 1e7) against P_ztwd(rho) along the
model profile, which is arithmetic and needs no cluster.

The first output already points at the handoff. At t = 0.043 s the hottest
cell inside the star sits at rho = 1.0e5, right at the density cut, and the
AMBIENT is at 5.1e8 K -- hotter than anywhere in the interior. The disturbance
starts at the surface, where the pressure scale height is smallest and a
percent-level mismatch in the equation of state does the most damage, and only
reaches rho ~ 4e8 by t = 0.22 s.
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

# The magnetic energy the field has to spend, from the initial model. The same
# number the reports quote for "the energy with no thermal destination".
E_MAG_BUDGET = 6.06e49

# Column indices in the fthermal output, which is documented in its own header
# block in tools/fthermal.cpp.
T, RHO_MAX, T_CORE, T_MEAN, T_MAX, T_AMB, RHO_AT_TMAX = 0, 1, 2, 3, 4, 5, 6
F_2E9, E_ION, MASS = 13, 14, 15


def read(name):
    """Rows of a fthermal/fbtbp CSV, comments dropped, as floats."""
    path = os.path.join(HERE, name)
    with open(path) as fh:
        rows = [r for r in csv.reader(fh) if r and not r[0].lstrip().startswith("#")]
    return [[float(x) for x in r] for r in rows]


def col(rows, i):
    return [r[i] for r in rows]


hz = read("thermal_hz192.csv")
ct = read("thermal_hz192ctl.csv")
fd = read("field_hz192.csv")

t_end = hz[-1][T]          # where the with-field run stopped, 1.27 s

fig, ax = plt.subplots(2, 2, figsize=(11.5, 8.4))
fig.suptitle(
    "Campanha HZ contra o próprio controle: o aquecimento não é do campo",
    fontsize=13,
)

# The whole result is that these two curves lie on top of each other, so they
# are drawn to make that readable rather than to hide one under the other: the
# control is a wide pale band, the with-field run a thin dark line on top of
# it. Two lines of equal weight would just paint the second over the first and
# the figure would look like a single measurement.
kw_ct = dict(color="steelblue", lw=4.5, alpha=0.45, solid_capstyle="round",
             label=r"controle, $B=0$ (hz192ctl), até 12 s")
kw_hz = dict(color="crimson", lw=1.4, ls="--", marker="o", ms=3.2,
             label="com campo (hz192), até 1.27 s")

# ---------------------------------------------------------------- (a) rho_max
a = ax[0, 0]
a.semilogy(col(ct, T), col(ct, RHO_MAX), **kw_ct)
a.semilogy(col(hz, T), col(hz, RHO_MAX), **kw_hz)
a.axhline(3.0e9, color="0.5", ls=":", lw=1)
a.text(2.6, 3.8e9, r"$\rho_{\max}(0)=3.0\times10^{9}$", fontsize=8, color="0.35")
a.annotate(r"fator $2800$", xy=(4.6, 1.9e5), xytext=(6.5, 8.0e5), fontsize=8,
           arrowprops=dict(arrowstyle="->", color="0.3", lw=0.9))
a.set_ylim(1.0e5, 1.2e10)
a.set_xlabel("t (s)")
a.set_ylabel(r"$\rho_{\max}$ (g cm$^{-3}$)")
a.set_title("(a) a estrela é destruída, com campo ou sem", fontsize=10)
a.legend(fontsize=8, loc="lower left")
a.grid(alpha=0.3)

# ------------------------------------------------------------------ (b) E_ion
b = ax[0, 1]
b.semilogy(col(ct, T), col(ct, E_ION), **kw_ct)
b.semilogy(col(hz, T), col(hz, E_ION), **kw_hz)
b.axhline(E_MAG_BUDGET, color="k", ls="--", lw=1.2)
b.text(3.4, E_MAG_BUDGET * 1.3,
       r"todo o campo: $6.06\times10^{49}$ erg", fontsize=8)
b.annotate(
    r"$8.4\times$ o que o campo tem",
    xy=(0.62, 5.11e50), xytext=(2.9, 3.2e50), fontsize=8,
    arrowprops=dict(arrowstyle="->", color="0.3", lw=0.9),
)
b.set_xlabel("t (s)")
b.set_ylabel(r"$E_{\rm ion}$ (erg)")
b.set_title("(b) a energia térmica excede o que o campo tem", fontsize=10)
b.grid(alpha=0.3)

# ------------------------------------------------------- (c) mass above 2e9 K
c = ax[1, 0]
c.plot(col(ct, T), col(ct, F_2E9), **kw_ct)
c.plot(col(hz, T), col(hz, F_2E9), **kw_hz)
c.axhline(0.05, color="0.5", ls=":", lw=1)
c.text(1.02, 0.10, "a casca de ~5% que a §10.1 supôs",
       fontsize=8, color="0.35")
c.set_xlim(0, 3.0)
c.set_ylim(-0.03, 1.03)
c.set_xlabel("t (s)")
c.set_ylabel(r"fração de massa com $T > 2\times10^{9}$ K")
c.set_title("(c) não é uma casca: é 94% da massa", fontsize=10)
c.grid(alpha=0.3)

# ------------------------------------------------------------------ (d) E_mag
d = ax[1, 1]
e_mag = [r[1] + r[2] for r in fd]          # E_tor + E_pol
d.semilogy(col(fd, T), e_mag, color="crimson", lw=1.8, marker="o", ms=3.5,
           label=r"$E_{\rm mag}$, run com campo")
d.set_xlabel("t (s)")
d.set_ylabel(r"$E_{\rm mag}$ (erg)")
d.set_title(r"(d) o campo perde $4.7\times10^{49}$ erg, e não é o que aquece",
            fontsize=10)
d.legend(fontsize=8)
d.grid(alpha=0.3)

fig.tight_layout(rect=(0, 0, 1, 0.95))
for ext in ("pdf", "png"):
    out = os.path.join(HERE, f"hz_control.{ext}")
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")

# The numbers the text quotes, printed so they cannot drift away from the
# figure they are supposed to describe.
de_hz = max(col(hz, E_ION)) - hz[0][E_ION]
de_ct = max(col(ct, E_ION)) - ct[0][E_ION]
print(f"\ndE_ion, com campo   = {de_hz:.3e} erg  ({de_hz / E_MAG_BUDGET:.2f}x o campo)")
print(f"dE_ion, controle    = {de_ct:.3e} erg  ({de_ct / E_MAG_BUDGET:.2f}x o campo)")
print(f"dE_mag medido       = {e_mag[0] - e_mag[-1]:.3e} erg")
print(f"rho_max, controle   : {ct[0][RHO_MAX]:.3e} -> {ct[-1][RHO_MAX]:.3e} "
      f"(fator {ct[0][RHO_MAX] / ct[-1][RHO_MAX]:.0f}) em {ct[-1][T]:.1f} s")
print(f"M(rho>1e5), controle: {ct[0][MASS]:.3f} -> {ct[-1][MASS]:.3f} Msun")
print(f"max f(T>2e9)        : {max(col(hz, F_2E9)):.3f} com campo, "
      f"{max(col(ct, F_2E9)):.3f} no controle")
