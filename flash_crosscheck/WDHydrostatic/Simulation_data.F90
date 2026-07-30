!!****if* source/Simulation/SimulationMain/WDHydrostatic/Simulation_data
!!
!! NAME
!!   Simulation_data
!!
!! DESCRIPTION
!!   Runtime parameters and the 1D model table for the hydrostatic white
!!   dwarf problem. The table is read once in Simulation_init and
!!   interpolated per cell in Simulation_initBlock.
!!
!!***

module Simulation_data

  implicit none

#include "constants.h"

  integer, parameter :: SIM_NPTS_MAX = 8192

  ! runtime parameters
  character(len=MAX_STRING_LENGTH), save :: sim_modelFile
  integer, save :: sim_nSubZones
  real, save    :: sim_rhoAmbient, sim_tempAmbient
  real, save    :: sim_abar, sim_zbar
  real, save    :: sim_xctr, sim_yctr, sim_zctr
  real, save    :: sim_smallRho, sim_smallT

  ! global velocity damping, mirroring Castro's problem_source.H
  real, save    :: sim_dampTimescaleTdyn, sim_dampEndTdyn
  real, save    :: sim_dampRampStartTdyn, sim_tDyn

  ! the 1D model, radius-ordered
  integer, save :: sim_nPts
  real, save    :: sim_rTab(SIM_NPTS_MAX)
  real, save    :: sim_rhoTab(SIM_NPTS_MAX)
  real, save    :: sim_tempTab(SIM_NPTS_MAX)

  integer, save :: sim_meshMe

end module Simulation_data
