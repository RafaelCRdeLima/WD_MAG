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
!!***

subroutine Simulation_adjustEvolution(blkcnt, blklst, nstep, dt, stime)

  use Simulation_data, ONLY : sim_dampTimescaleTdyn, sim_dampEndTdyn, &
                              sim_dampRampStartTdyn, sim_tDyn, sim_meshMe
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

  if (sim_dampTimescaleTdyn <= 0.0) return

  tEnd  = sim_dampEndTdyn       * sim_tDyn
  tRamp = sim_dampRampStartTdyn * sim_tDyn

  if (stime >= tEnd) return

  rate = 1.0 / (sim_dampTimescaleTdyn * sim_tDyn)

  ramp = 1.0
  if (stime > tRamp) then
     xr = (stime - tRamp) / (tEnd - tRamp)
     ramp = 0.5 * (1.0 + cos(PI * xr))
  end if

  fac = exp(-rate * ramp * dt)

  do b = 1, blkcnt
     call Grid_getBlkIndexLimits(blklst(b), blkLimits, blkLimitsGC)
     call Grid_getBlkPtr(blklst(b), solnData, CENTER)

     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)

              solnData(VELX_VAR, i, j, k) = fac * solnData(VELX_VAR, i, j, k)
              solnData(VELY_VAR, i, j, k) = fac * solnData(VELY_VAR, i, j, k)
              solnData(VELZ_VAR, i, j, k) = fac * solnData(VELZ_VAR, i, j, k)

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
     print *, '[WDHydrostatic] damping: t/t_dyn =', stime / sim_tDyn, &
              ' ramp =', ramp, ' factor =', fac
  end if

end subroutine Simulation_adjustEvolution
