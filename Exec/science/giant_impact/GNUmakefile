PRECISION ?= DOUBLE
PROFILE ?= FALSE

DEBUG ?= FALSE
BACKTRACE ?= FALSE
TEST ?= FALSE

USE_MPI ?= TRUE
USE_OMP ?= FALSE

COMP ?= gnu

DIM ?= 3

USE_GRAV ?= TRUE

USE_CXX_MODEL_PARSER = TRUE
NUM_MODELS = 2
MAX_NPTS_MODEL ?= 100000

GPU_COMPATIBLE_PROBLEM = TRUE

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos

EOS_DIR ?= tillotson

# This sets the network directory in $(MICROPHYSICS_HOME)/networks

NETWORK_DIR ?= general_null
NETWORK_INPUTS = gammalaw.net

# Location where the source files are coming from.

PROBLEM_DIR ?= $(CASTRO_HOME)/Exec/science/giant_impact

Blocs = $(PROBLEM_DIR)
Bpack = $(PROBLEM_DIR)/Make.package

include $(CASTRO_HOME)/Exec/Make.Castro
