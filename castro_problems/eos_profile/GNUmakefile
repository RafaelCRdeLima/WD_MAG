# The DIARIO 11 check, which has been open since August: does the ztwd profile
# the SCF built remain an equilibrium under helmholtz?
#
# Section 11 closed the HZ campaign negatively -- the field-free control heated
# identically, so the heating was not magnetic -- and pointed at the EOS
# handoff as the suspect: the star is balanced against P_ztwd(rho) and then
# evolved against P_helmholtz(rho, T). If those differ by a percent, the star
# is out of equilibrium everywhere from t = 0, and the surface, where the scale
# height is smallest, is where it shows first. That is exactly where the heat
# appeared (DIARIO 11).
#
# The comparison has to use the SAME EOS the runs used, not a reimplementation,
# so this builds against Microphysics twice:
#     make EOS_DIR=ztwd       -> main3d.gnu.ztwd.ex
#     make EOS_DIR=helmholtz  -> main3d.gnu.helmholtz.ex
# Both read the same density list and print P, so the difference is the answer.
#
# mu2.net, the same composition the runs carry: A = 4, Z = 2, hence mu_e = 2,
# which is what ztwd assumes.
PRECISION  = DOUBLE
PROFILE    = FALSE
DEBUG      = FALSE
DIM        = 3
COMP       = gnu
USE_MPI    = FALSE
USE_OMP    = FALSE
USE_REACT  = FALSE

EBASE = main
BL_NO_FORT = TRUE

MICROPHYSICS_HOME ?= ../../microphysics

EOS_DIR     ?= ztwd
NETWORK_DIR := general_null
NETWORK_INPUTS = mu2.net

EXTERN_SEARCH += . $(MICROPHYSICS_HOME)/unit_test

DEFINES += -DEOS_NAME=\"$(EOS_DIR)\"

Bpack := ./Make.package
Blocs := .

include $(MICROPHYSICS_HOME)/unit_test/Make.unit_test
