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

# This sets the EOS directory in $(MICROPHYSICS_HOME)/EOS -- gamma_law,
# not helmholtz: Step 2 only needs a background density profile to
# confine the seed field into (see problem_initialize.H), not a
# physically accurate degenerate-matter EOS inside Castro. See
# scf/castro_model_writer.py's module docstring for why the model
# file's temperature/pressure columns are placeholders.
EOS_DIR     := gamma_law

# This sets the network directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := general_null
NETWORK_INPUTS = gammalaw.net

PROBLEM_DIR ?= ./

Bpack   := $(PROBLEM_DIR)/Make.package
Blocs   := $(PROBLEM_DIR)

include $(CASTRO_HOME)/Exec/Make.Castro
