!!****if* source/Simulation/SimulationMain/WDHydrostatic/Simulation_init
!!
!! NAME
!!   Simulation_init
!!
!! DESCRIPTION
!!   Reads the runtime parameters and the 1D white-dwarf model produced by
!!   flash_crosscheck/make_wd_model.py. The file format is deliberately
!!   minimal -- comment lines beginning with '#', then rows of
!!   "r density temperature" -- because both the writer and this reader are
!!   ours, and the FLASH-native multi-header format buys nothing here.
!!
!!   Every rank reads the file. It is a couple of thousand lines, so the
!!   simplicity is worth more than avoiding the redundant reads.
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype

  implicit none

#include "constants.h"
#include "Flash.h"

  character(len=256) :: line
  integer :: lun, ios, i
  real    :: rr, dd, tt

  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get('sim_modelFile',   sim_modelFile)
  call RuntimeParameters_get('sim_nSubZones',   sim_nSubZones)
  call RuntimeParameters_get('sim_rhoAmbient',  sim_rhoAmbient)
  call RuntimeParameters_get('sim_tempAmbient', sim_tempAmbient)
  call RuntimeParameters_get('sim_abar',        sim_abar)
  call RuntimeParameters_get('sim_zbar',        sim_zbar)
  call RuntimeParameters_get('sim_xctr',        sim_xctr)
  call RuntimeParameters_get('sim_yctr',        sim_yctr)
  call RuntimeParameters_get('sim_zctr',        sim_zctr)
  call RuntimeParameters_get('smlrho',          sim_smallRho)
  call RuntimeParameters_get('smallt',          sim_smallT)

  lun = 33
  open(unit=lun, file=trim(sim_modelFile), status='old', iostat=ios)
  if (ios /= 0) then
     call Driver_abortFlash('[Simulation_init] cannot open sim_modelFile: ' &
                            // trim(sim_modelFile))
  end if

  sim_nPts = 0
  do
     read(lun, '(a)', iostat=ios) line
     if (ios /= 0) exit
     line = adjustl(line)
     if (len_trim(line) == 0) cycle
     if (line(1:1) == '#') cycle

     read(line, *, iostat=ios) rr, dd, tt
     if (ios /= 0) cycle

     if (sim_nPts >= SIM_NPTS_MAX) then
        call Driver_abortFlash('[Simulation_init] model file longer than SIM_NPTS_MAX')
     end if
     sim_nPts = sim_nPts + 1
     sim_rTab(sim_nPts)    = rr
     sim_rhoTab(sim_nPts)  = max(dd, sim_smallRho)
     sim_tempTab(sim_nPts) = max(tt, sim_smallT)
  end do
  close(lun)

  if (sim_nPts < 2) then
     call Driver_abortFlash('[Simulation_init] model file has fewer than 2 usable rows')
  end if

  ! Monotonic radius is assumed by the bisection in Simulation_initBlock;
  ! check it rather than trust the generator.
  do i = 2, sim_nPts
     if (sim_rTab(i) <= sim_rTab(i-1)) then
        call Driver_abortFlash('[Simulation_init] model radii are not strictly increasing')
     end if
  end do

  if (sim_meshMe == MASTER_PE) then
     print *, '[WDHydrostatic] model file: ', trim(sim_modelFile)
     print *, '[WDHydrostatic] points read = ', sim_nPts
     print *, '[WDHydrostatic] r range     = ', sim_rTab(1), sim_rTab(sim_nPts)
     print *, '[WDHydrostatic] rho_c       = ', sim_rhoTab(1)
     print *, '[WDHydrostatic] T (uniform) = ', sim_tempTab(1)
     print *, '[WDHydrostatic] abar, zbar  = ', sim_abar, sim_zbar
     print *, '[WDHydrostatic] nSubZones   = ', sim_nSubZones
  end if

end subroutine Simulation_init
