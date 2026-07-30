!!****if* source/Simulation/SimulationMain/WDHydrostatic/Simulation_initBlock
!!
!! NAME
!!   Simulation_initBlock
!!
!! DESCRIPTION
!!   Maps the 1D white-dwarf model onto a 3D block: for each cell, the
!!   density and temperature are averaged over sim_nSubZones^3 sub-points,
!!   which matters because point-sampling the centre of a cell in a
!!   steeply-stratified core misses the cell's mean by a measurable amount.
!!   Castro's own version of this problem was found to reproduce the target
!!   central density to -4.78% under point sampling
!!   (docs/teoline.md Sec 6.6), so volume averaging is not a refinement
!!   here, it is the difference between starting from the intended star and
!!   starting from a different one.
!!
!!   Velocities are set to zero: the star is meant to be in equilibrium,
!!   and any motion that appears is the answer to the question being asked.
!!
!!   Pressure and internal energy are NOT computed here. DENS, TEMP and the
!!   composition are set and Eos_wrapped is called in MODE_DENS_TEMP, so
!!   the thermodynamic state is whatever Helmholtz says it is -- computing
!!   P by hand would reintroduce exactly the EOS mismatch this setup exists
!!   to avoid.
!!
!!***

subroutine Simulation_initBlock(blockID)

  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getCellCoords, &
                             Grid_getBlkPtr, Grid_releaseBlkPtr
  use Eos_interface, ONLY : Eos_wrapped

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  integer, intent(in) :: blockID

  integer, dimension(2, MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX, sizeY, sizeZ
  real, allocatable, dimension(:) :: xCen, yCen, zCen
  real, allocatable, dimension(:) :: xL, yL, zL, xR, yR, zR
  real, pointer, dimension(:,:,:,:) :: solnData

  integer :: i, j, k, ii, jj, kk, n
  real    :: dxSub, dySub, dzSub, xx, yy, zz, rad
  real    :: rhoSum, tempSum, rhoI, tempI
  logical :: gcell = .true.

  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

  sizeX = blkLimitsGC(HIGH, IAXIS) - blkLimitsGC(LOW, IAXIS) + 1
  sizeY = blkLimitsGC(HIGH, JAXIS) - blkLimitsGC(LOW, JAXIS) + 1
  sizeZ = blkLimitsGC(HIGH, KAXIS) - blkLimitsGC(LOW, KAXIS) + 1

  allocate(xCen(sizeX), yCen(sizeY), zCen(sizeZ))
  allocate(xL(sizeX), yL(sizeY), zL(sizeZ))
  allocate(xR(sizeX), yR(sizeY), zR(sizeZ))

  call Grid_getCellCoords(IAXIS, blockID, CENTER,     gcell, xCen, sizeX)
  call Grid_getCellCoords(IAXIS, blockID, LEFT_EDGE,  gcell, xL,   sizeX)
  call Grid_getCellCoords(IAXIS, blockID, RIGHT_EDGE, gcell, xR,   sizeX)
  call Grid_getCellCoords(JAXIS, blockID, CENTER,     gcell, yCen, sizeY)
  call Grid_getCellCoords(JAXIS, blockID, LEFT_EDGE,  gcell, yL,   sizeY)
  call Grid_getCellCoords(JAXIS, blockID, RIGHT_EDGE, gcell, yR,   sizeY)
  call Grid_getCellCoords(KAXIS, blockID, CENTER,     gcell, zCen, sizeZ)
  call Grid_getCellCoords(KAXIS, blockID, LEFT_EDGE,  gcell, zL,   sizeZ)
  call Grid_getCellCoords(KAXIS, blockID, RIGHT_EDGE, gcell, zR,   sizeZ)

  call Grid_getBlkPtr(blockID, solnData, CENTER)

  n = sim_nSubZones

  do k = blkLimitsGC(LOW, KAXIS), blkLimitsGC(HIGH, KAXIS)
     dzSub = (zR(k - blkLimitsGC(LOW,KAXIS) + 1) - zL(k - blkLimitsGC(LOW,KAXIS) + 1)) / real(n)
     do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        dySub = (yR(j - blkLimitsGC(LOW,JAXIS) + 1) - yL(j - blkLimitsGC(LOW,JAXIS) + 1)) / real(n)
        do i = blkLimitsGC(LOW, IAXIS), blkLimitsGC(HIGH, IAXIS)
           dxSub = (xR(i - blkLimitsGC(LOW,IAXIS) + 1) - xL(i - blkLimitsGC(LOW,IAXIS) + 1)) / real(n)

           rhoSum  = 0.0
           tempSum = 0.0
           do kk = 1, n
              zz = zL(k - blkLimitsGC(LOW,KAXIS) + 1) + (real(kk) - 0.5) * dzSub - sim_zctr
              do jj = 1, n
                 yy = yL(j - blkLimitsGC(LOW,JAXIS) + 1) + (real(jj) - 0.5) * dySub - sim_yctr
                 do ii = 1, n
                    xx = xL(i - blkLimitsGC(LOW,IAXIS) + 1) + (real(ii) - 0.5) * dxSub - sim_xctr
                    rad = sqrt(xx*xx + yy*yy + zz*zz)
                    call sim_interp1d(rad, rhoI, tempI)
                    rhoSum  = rhoSum  + rhoI
                    tempSum = tempSum + tempI
                 end do
              end do
           end do

           solnData(DENS_VAR, i, j, k) = max(rhoSum / real(n*n*n), sim_smallRho)
           solnData(TEMP_VAR, i, j, k) = max(tempSum / real(n*n*n), sim_smallT)
           solnData(VELX_VAR, i, j, k) = 0.0
           solnData(VELY_VAR, i, j, k) = 0.0
           solnData(VELZ_VAR, i, j, k) = 0.0
           solnData(SPECIES_BEGIN, i, j, k) = 1.0

        end do
     end do
  end do

  call Grid_releaseBlkPtr(blockID, solnData, CENTER)

  ! Let Helmholtz set PRES/EINT/ENER/GAMC/GAME from (rho, T, composition).
  call Eos_wrapped(MODE_DENS_TEMP, blkLimitsGC, blockID)

  deallocate(xCen, yCen, zCen, xL, yL, zL, xR, yR, zR)

contains

  !! Linear interpolation of the 1D table, by bisection. Outside the table
  !! the ambient values are used, which is the vacuum floor, not an
  !! extrapolation of the stellar profile.
  subroutine sim_interp1d(rad, rhoOut, tempOut)
    real, intent(in)  :: rad
    real, intent(out) :: rhoOut, tempOut
    integer :: lo, hi, mid
    real    :: w

    if (rad <= sim_rTab(1)) then
       rhoOut  = sim_rhoTab(1)
       tempOut = sim_tempTab(1)
       return
    end if
    if (rad >= sim_rTab(sim_nPts)) then
       rhoOut  = sim_rhoAmbient
       tempOut = sim_tempAmbient
       return
    end if

    lo = 1
    hi = sim_nPts
    do while (hi - lo > 1)
       mid = (lo + hi) / 2
       if (sim_rTab(mid) > rad) then
          hi = mid
       else
          lo = mid
       end if
    end do

    w = (rad - sim_rTab(lo)) / (sim_rTab(hi) - sim_rTab(lo))
    rhoOut  = (1.0 - w) * sim_rhoTab(lo)  + w * sim_rhoTab(hi)
    tempOut = (1.0 - w) * sim_tempTab(lo) + w * sim_tempTab(hi)
  end subroutine sim_interp1d

end subroutine Simulation_initBlock
