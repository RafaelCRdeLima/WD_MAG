!!****if* source/Simulation/SimulationMain/WDHydrostatic/Simulation_adjustEvolution
!!
!! NAME
!!   Simulation_adjustEvolution
!!
!! DESCRIPTION
!!   Global velocity damping, the equivalent of the Castro problem's
!!   problem_source.H (docs/teoria.md Sec 6.4). Without it the comparison
!!   against the Castro reference curve is meaningless, because that curve
!!   was produced with damping active over [0, 20 t_dyn] -- the run log
!!   says so: "global damping active [0, 5.516124195 s], ramping off from
!!   4.964511776 s", which is 20 and 18 t_dyn.
!!
!!   Castro's form, reproduced here:
!!
!!     rate    = 1 / (damping_timescale_in_tdyn * t_dyn)
!!     ramp(t) = 1                              for t <= ramp_start
!!               0.5*(1 + cos(pi*x))            for ramp_start < t < end
!!               (x = (t - ramp_start)/(end - ramp_start))
!!     d(rho u)/dt = -rate * ramp * rho u
!!     d(rho E)/dt = u . d(rho u)/dt
!!
!!   Applied to EVERY cell, not just the exterior -- that is the point of
!!   it, and what distinguishes it from a sponge.
!!
!!   Two deliberate differences from Castro, both recorded because they are
!!   the reasons this is "equivalent" and not "identical":
!!
!!   1. Castro adds this as a source term inside the hydro update
!!      (castro.add_ext_src=1). FLASH has no equivalent hook in this
!!      driver, so it is applied operator-split, once per step, through
!!      Simulation_adjustEvolution -- which the driver calls every step
!!      for exactly this kind of thing. First-order equivalent.
!!
!!   2. The decay over a step is taken as exp(-rate*ramp*dt) rather than
!!      Castro's explicit (1 - rate*ramp*dt). It is the exact solution of
!!      the same ODE over the step and is unconditionally stable, whereas
!!      the explicit form breaks down once rate*dt > 1. At the timestep
!!      this run actually takes (rate*dt ~ 0.16) the two differ by ~1.5%
!!      in the per-step factor.
!!
!!   The energy update removes kinetic energy only: EINT is untouched and
!!   ENER is rebuilt from it, which is the same kinetic-energy-consistent
!!   treatment as Castro's SrE = u . Sr. DENS and EINT unchanged means
!!   PRES and TEMP are unchanged too, so no Eos call is needed.
!!
!!   SECOND TERM: the exterior sponge, ported from Castro_sponge.cpp.
!!   The global damping above is not what keeps Castro's ambient in place --
!!   the sponge is, and leaving it out is why the first FLASH runs died at
!!   the star/vacuum interface. The ambient has no pressure support worth
!!   speaking of, so it free-falls onto the star, arrives at near free-fall
!!   speed, and dumps enough specific energy into near-empty cells to drive
!!   the EOS off its table. Castro damps that with a density-gated sponge on
!!   a timescale of 1e-4 s -- roughly 600 times more aggressive than the
!!   global damping, and applied only where the density is low:
!!
!!     alpha  = dt / sponge_timescale
!!     f      = 0                                   for rho > upper_density
!!              0.5*(1 - cos(pi*(rho-upper)/(lower-upper)))
!!                                                  for lower <= rho <= upper
!!              1                                   for rho < lower_density
!!     (rho v) -> (rho v) / (1 + alpha*f)
!!
!!   The update is implicit in exactly Castro's form, which is what makes it
!!   stable at alpha ~ 90 (which is where dt/1e-4 lands here). Same
!!   kinetic-energy-consistent energy treatment.
!!
!!***

subroutine Simulation_adjustEvolution(blkcnt, blklst, nstep, dt, stime)

  use Simulation_data, ONLY : sim_dampTimescaleTdyn, sim_dampEndTdyn, &
                              sim_dampRampStartTdyn, sim_tDyn, sim_meshMe, &
                              sim_spongeTimescale, sim_spongeUpperDens, &
                              sim_spongeLowerDens
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
                             Grid_releaseBlkPtr

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: blkcnt
  integer, intent(in) :: blklst(blkcnt)
  integer, intent(in) :: nstep
  real, intent(in) :: dt
  real, intent(in) :: stime

  integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
  real, pointer, dimension(:,:,:,:) :: solnData
  integer :: b, i, j, k
  real :: tEnd, tRamp, rate, ramp, xr, fac, ekin
  real :: alpha, deltaRho, rho, sf, spFac
  logical :: doDamp, doSponge

  tEnd  = sim_dampEndTdyn       * sim_tDyn
  tRamp = sim_dampRampStartTdyn * sim_tDyn

  doDamp = (sim_dampTimescaleTdyn > 0.0) .and. (stime < tEnd)
  doSponge = (sim_spongeTimescale > 0.0) .and. &
             (sim_spongeUpperDens > 0.0) .and. (sim_spongeLowerDens > 0.0)

  if (.not. (doDamp .or. doSponge)) return

  fac = 1.0
  ramp = 0.0
  if (doDamp) then
     rate = 1.0 / (sim_dampTimescaleTdyn * sim_tDyn)
     ramp = 1.0
     if (stime > tRamp) then
        xr = (stime - tRamp) / (tEnd - tRamp)
        ramp = 0.5 * (1.0 + cos(PI * xr))
     end if
     fac = exp(-rate * ramp * dt)
  end if

  alpha = 0.0
  deltaRho = 1.0
  if (doSponge) then
     alpha = dt / sim_spongeTimescale
     deltaRho = sim_spongeLowerDens - sim_spongeUpperDens
  end if

  do b = 1, blkcnt
     call Grid_getBlkIndexLimits(blklst(b), blkLimits, blkLimitsGC)
     call Grid_getBlkPtr(blklst(b), solnData, CENTER)

     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)

              ! global damping, then the exterior sponge on top of it
              spFac = 1.0
              if (doSponge) then
                 rho = solnData(DENS_VAR, i, j, k)
                 if (rho > sim_spongeUpperDens) then
                    sf = 0.0
                 else if (rho >= sim_spongeLowerDens) then
                    sf = 0.5 * (1.0 - cos(PI * (rho - sim_spongeUpperDens) / deltaRho))
                 else
                    sf = 1.0
                 end if
                 spFac = 1.0 / (1.0 + alpha * sf)
              end if

              solnData(VELX_VAR, i, j, k) = fac * spFac * solnData(VELX_VAR, i, j, k)
              solnData(VELY_VAR, i, j, k) = fac * spFac * solnData(VELY_VAR, i, j, k)
              solnData(VELZ_VAR, i, j, k) = fac * spFac * solnData(VELZ_VAR, i, j, k)

              ekin = 0.5 * (solnData(VELX_VAR, i, j, k)**2 &
                          + solnData(VELY_VAR, i, j, k)**2 &
                          + solnData(VELZ_VAR, i, j, k)**2)
              solnData(ENER_VAR, i, j, k) = solnData(EINT_VAR, i, j, k) + ekin

           end do
        end do
     end do

     call Grid_releaseBlkPtr(blklst(b), solnData, CENTER)
  end do

  if (sim_meshMe == MASTER_PE .and. mod(nstep, 20) == 1) then
     print *, '[WDHydrostatic] t/t_dyn =', stime / sim_tDyn, &
              ' damp factor =', fac, ' sponge alpha =', alpha
  end if

end subroutine Simulation_adjustEvolution
