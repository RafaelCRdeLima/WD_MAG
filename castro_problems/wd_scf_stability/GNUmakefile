PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE

USE_MHD    = TRUE
# Step 3 (evolution): the star must be held together by gravity during
# relaxation, or it just flies apart from unbalanced pressure -- Step
# 2's USE_GRAV=FALSE was fine for a static (max_step=0) IC check, but
# is wrong for any real time evolution. MonopoleGrav (see inputs) is
# the isolated-star gravity solve, same choice as
# Exec/gravity_tests/StarGrav.
USE_GRAV   = TRUE

USE_MODEL_PARSER = TRUE
NUM_MODELS = 1

# define the location of the CASTRO top directory
CASTRO_HOME  ?= ../../..

# This sets the EOS directory in $(MICROPHYSICS_HOME)/EOS -- ztwd
# (zero-temperature white dwarf: Chandrasekhar 1935 degenerate-electron
# P(rho), Hachisu 1986 enthalpy), the SAME formula scf/eos.py uses to
# build the background star. Originally this was gamma_law (a
# placeholder, since Step 2 only needed a density profile to confine
# the seed field into, not accurate thermodynamics) -- that choice
# caused a real bug found in Step 3: the star was built in hydrostatic
# equilibrium under ztwd but evolved under gamma_law, so it was never
# in equilibrium under the EOS actually governing its evolution and
# collapsed from t=0 (confirmed with a field-free E_mag=0 control run,
# which collapsed identically -- see docs/teoria.md Sec 6.4/6.5).
# ztwd ignores temperature entirely for eos_input_rt (P depends only on
# rho and mu_e), which also makes model.dat's placeholder temperature
# column physically inert rather than merely approximate.
EOS_DIR     := ztwd

# This sets the network directory in $(MICROPHYSICS_HOME)/networks --
# mu2.net (in this directory, not Microphysics/networks/general_null/)
# is a single species "X" with aion=4, zion=2, giving mu_e=2.0 to match
# scf/eos.py's default. Species name kept as "X" so model.dat's column
# header doesn't need to change, only its composition.
NETWORK_DIR := general_null
NETWORK_INPUTS = mu2.net

PROBLEM_DIR ?= ./

Bpack   := $(PROBLEM_DIR)/Make.package
Blocs   := $(PROBLEM_DIR)

include $(CASTRO_HOME)/Exec/Make.Castro
