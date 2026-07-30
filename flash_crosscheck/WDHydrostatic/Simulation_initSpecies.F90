!!****if* source/Simulation/SimulationMain/WDHydrostatic/Simulation_initSpecies
!!
!! NAME
!!   Simulation_initSpecies
!!
!! DESCRIPTION
!!   One species, He4-like: A = 4, Z = 2, so mu_e = A/Z = 2. That is the
!!   same composition the Castro problem uses (mu2.net), which is what
!!   makes the two runs describe the same star rather than two stars that
!!   happen to share a central density.
!!
!!   A and Z come from runtime parameters rather than being hardwired, so
!!   the composition cannot silently disagree with the model file that
!!   make_wd_model.py built under the same assumption.
!!
!!***

subroutine Simulation_initSpecies()

  use Multispecies_interface, ONLY : Multispecies_setProperty
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "Multispecies.h"
#include "Flash.h"

  real :: abar, zbar

  call RuntimeParameters_get('sim_abar', abar)
  call RuntimeParameters_get('sim_zbar', zbar)

  call Multispecies_setProperty(HE4_SPEC, A, abar)
  call Multispecies_setProperty(HE4_SPEC, Z, zbar)
  call Multispecies_setProperty(HE4_SPEC, GAMMA, 5.0e0/3.0e0)
  call Multispecies_setProperty(HE4_SPEC, E, zbar)

end subroutine Simulation_initSpecies
