#!/usr/bin/env python3
"""Magnetic Prandtl number in the interior of our 2 Msun white dwarf.

Why
---
Pm = nu/eta is the parameter a shearing box has to match, now that the
incompressible formulation has removed beta from the problem entirely
(DIARIO 6.14). It also decides whether the published MInIT closure
coefficients -- calibrated for protoneutron stars -- transfer to a white
dwarf, and whether we sit above or below the Pm ~ 2-4 threshold at which the
disc literature finds MRI turbulence dying.

The value in circulation, Pm ~ 0.58, is quoted for the convective region of a
COOL CRYSTALLISING CO white dwarf.  Ours is a hot merger remnant at 30-1000x
that density.  This recomputes it.

Physics
-------
Electron transport dominates in degenerate matter.  Shternin (2008),
arXiv:0803.3893, Eq. (1)-(2):

    eta_visc = n_e v_F p_F / (5 nu_e),      nu_e = nu_ee + nu_ei

with p_F = hbar (3 pi^2 n_e)^(1/3), m* = E_F/c^2, v_F = p_F/m*.

The key asymmetry, and the reason Pm is not simply 1: electron-electron
collisions conserve total electron momentum, so they do NOT degrade a current
and do not enter the electrical conductivity -- but they DO degrade shear
stress and so enter the viscosity.  Hence

    sigma  = n_e e^2 / (m* nu_ei)                 (ei only)
    eta_mag = c^2 / (4 pi sigma)
    nu_kin = eta_visc / rho
    Pm     = nu_kin / eta_mag

nu_ee is taken from Shternin's Eq. (9) with the asymptotic I_eta of his
Eqs. (14)-(15); nu_ei from the standard degenerate form
nu_ei = 4 pi Z e^4 n_i Lambda / (p_F^2 v_F).

Both nu and eta scale as 1/nu_coll, so Pm ~ 1/nu_coll^2.  Because nu_ei is
proportional to n_i while nu_kin and sigma each carry a compensating n_e, the
ion density cancels and Pm depends on temperature only through nu_ee and the
Coulomb logarithm -- it is far less T-sensitive than one might guess.

Honest limitations
------------------
1. Lambda_ei is the weak point.  Pm goes as 1/Lambda^2 and Lambda is of order
   unity in dense matter but not accurately known without conductivity tables
   (Potekhin's condegin, Itoh et al.).  It is a free parameter here and the
   sensitivity is printed.
2. Ion-ion correlations in a strongly coupled liquid suppress nu_ei below the
   uncorrelated estimate, which would push Pm UP.  Gamma is printed so the
   regime is visible.
3. Ion viscosity and radiative viscosity are ignored.  In the degenerate
   interior electrons dominate; in a cool convective envelope they may not,
   which is the likeliest reason the literature value differs.
4. Our EOS is barotropic and carries no temperature at all, so T is an
   assumption about a merger remnant, not a simulation output.  Hence a scan.
"""

import math

# CGS
HBAR = 1.054572e-27
C = 2.997925e10
ME = 9.109384e-28
E = 4.803205e-10
KB = 1.380649e-16
MU = 1.660539e-24

MU_E = 2.0     # carbon/oxygen
Z_ION = 6.0    # carbon; oxygen (Z=8) shifts nu_ei, printed as sensitivity
A_ION = 12.0


def electron_state(rho):
    n_e = rho / (MU_E * MU)
    p_F = HBAR * (3.0 * math.pi**2 * n_e) ** (1.0 / 3.0)
    x = p_F / (ME * C)
    E_F = ME * C**2 * math.sqrt(1.0 + x * x)
    m_star = E_F / C**2
    v_F = p_F / m_star
    return n_e, p_F, x, E_F, m_star, v_F


def plasma_temp(n_e, m_star):
    omega_pe = math.sqrt(4.0 * math.pi * E**2 * n_e / m_star)
    return HBAR * omega_pe / KB


def I_eta(u, theta):
    """Shternin 2008 Eqs. (14) cold / (15) warm asymptotics.

    theta = sqrt(3) T_pe/T, so theta >> 1 is the COLD plasma (T << T_pe),
    regimes II/IV; theta <~ 1 is warm, regimes I/III.  Only the asymptotes are
    implemented, so values near theta ~ 1 are interpolation-free and rough.
    """
    if theta >= 1.0:                                    # cold, Eq. (15)
        I_l = math.pi**3 / (12.0 * theta)
        xi = 1.813
        I_t = xi * u ** (10.0 / 3.0) / theta ** (2.0 / 3.0)
        I_lt = math.pi**3 * u**2 / (6.0 * theta)
    else:                                               # warm, Eq. (14)
        L = math.log(1.0 / theta)
        I_l = (2.0 / 3.0) * (L + 1.919)
        I_t = (1.0 / 3.0) * (L + 2.742)
        I_lt = (2.0 / 3.0) * (L + 2.052)
    return I_l + I_t + I_lt


def transport(rho, T, lam=1.0, Z=Z_ION, A=A_ION):
    n_e, p_F, x, E_F, m_star, v_F = electron_state(rho)
    u = v_F / C
    T_pe = plasma_temp(n_e, m_star)
    theta = math.sqrt(3.0) * T_pe / T

    alpha = E**2 / (HBAR * C)
    nu_ee = (12.0 * alpha**2 / (math.pi * HBAR)) * KB * T * (C / v_F) ** 2 \
        * I_eta(u, theta)

    n_i = rho / (A * MU)
    nu_ei = 4.0 * math.pi * Z * E**4 * n_i * lam / (p_F**2 * v_F)

    eta_visc = n_e * v_F * p_F / (5.0 * (nu_ee + nu_ei))
    nu_kin = eta_visc / rho

    sigma = n_e * E**2 / (m_star * nu_ei)          # ee does not degrade current
    eta_mag = C**2 / (4.0 * math.pi * sigma)

    a_i = (3.0 / (4.0 * math.pi * n_i)) ** (1.0 / 3.0)
    gamma = (Z * E) ** 2 / (a_i * KB * T)

    return dict(x=x, u=u, T_pe=T_pe, theta=theta, gamma=gamma,
                nu_ee=nu_ee, nu_ei=nu_ei, nu_kin=nu_kin, eta_mag=eta_mag,
                sigma=sigma, Pm=nu_kin / eta_mag)


def main():
    print("Electron transport in degenerate C/O matter, mu_e = 2, Z = 6\n")

    print(f"{'rho':>9} {'T':>8} {'x=pF/mc':>8} {'Gamma':>7} {'theta':>7} "
          f"{'nu_ee/nu_ei':>12} {'nu(cm2/s)':>10} {'eta(cm2/s)':>11} {'Pm':>9}")
    for rho in (5.0e7, 1.0e8, 5.0e8, 1.0e9):
        for T in (1.0e7, 1.0e8, 5.0e8, 1.0e9):
            r = transport(rho, T)
            print(f"{rho:>9.1e} {T:>8.1e} {r['x']:>8.2f} {r['gamma']:>7.1f} "
                  f"{r['theta']:>7.2f} {r['nu_ee']/r['nu_ei']:>12.3f} "
                  f"{r['nu_kin']:>10.3g} {r['eta_mag']:>11.3g} {r['Pm']:>9.4g}")
        print()

    # our star: mean density from the 256^3 run, assumed remnant temperature
    RHO, T = 4.8e7, 1.0e8
    r = transport(RHO, T)
    print(f"Our star, rho_mean = {RHO:.1e} g/cm^3, T = {T:.0e} K assumed:")
    print(f"  relativity x = p_F/m_e c   {r['x']:.2f}  (relativistic)")
    print(f"  ion coupling Gamma          {r['gamma']:.1f}  "
          f"({'liquid' if r['gamma'] < 175 else 'CRYSTAL'}, melting at ~175)")
    print(f"  nu_ee/nu_ei                 {r['nu_ee']/r['nu_ei']:.3f}  "
          f"(ee is a correction, not the limiter)")
    print(f"  sigma                       {r['sigma']:.3g} s^-1")
    print(f"  nu                          {r['nu_kin']:.3g} cm^2/s")
    print(f"  eta                         {r['eta_mag']:.3g} cm^2/s")
    print(f"  Pm                          {r['Pm']:.4g}")

    print("\nSensitivity to the Coulomb logarithm (Pm ~ 1/Lambda^2):")
    for lam in (0.5, 1.0, 2.0, 3.0, 5.0):
        r = transport(RHO, T, lam=lam)
        print(f"  Lambda = {lam:>4.1f}   nu = {r['nu_kin']:>9.3g}   "
              f"eta = {r['eta_mag']:>9.3g}   Pm = {r['Pm']:>9.4g}")

    print("\nSensitivity to composition:")
    for Z, A, name in ((6.0, 12.0, "carbon"), (8.0, 16.0, "oxygen")):
        r = transport(RHO, T, Z=Z, A=A)
        print(f"  {name:>7} (Z={Z:.0f})   Pm = {r['Pm']:.4g}")

    # Reynolds numbers on the stellar scale, for the record
    r = transport(RHO, T)
    L, V = 3.2e8, 2.2e9      # R_eq and the rotation velocity there
    print(f"\nOn the stellar scale L = {L:.1e} cm, V = {V:.1e} cm/s:")
    print(f"  Re = {L*V/r['nu_kin']:.3g}")
    print(f"  Rm = {L*V/r['eta_mag']:.3g}")


if __name__ == "__main__":
    main()
