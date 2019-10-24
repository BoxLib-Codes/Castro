PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE
OPTIMIZE   = TRUE 

DIM        = 2

COMP	   = pgi

USE_MPI    = TRUE
USE_OMP    = FALSE
USE_GPU    = FALSE
USE_CUDA   = FALSE
USE_ACC    = FALSE

USE_THORNADO_GPU = TRUE
USE_THORNADO_ACC = TRUE
ifeq ($(USE_THORNADO_GPU),TRUE)
  INCLUDE_LOCATIONS += $(OLCF_CUDA_ROOT)/include
  LIBRARY_LOCATIONS += $(OLCF_CUDA_ROOT)/lib64
  LIBRARIES += -lcusolver -lcublas -lcudart -lcuda
endif
ifeq ($(USE_THORNADO_ACC),TRUE)
  GENERIC_PGI_FLAGS += -acc -Minfo=accel -mcmodel=medium -ta=tesla:cc70,ptxinfo,cuda10.1
endif

USE_HDF5   = TRUE
HDF5_HOME  ?= /usr/local/Cellar/hdf5/1.10.4
ifdef OLCF_HDF5_ROOT
  HDF5_HOME := $(OLCF_HDF5_ROOT)
endif
ifeq ($(HOST),manda)
  HDF5_HOME := /home/almgren/bin/miniconda2
endif

LAPACK_HOME ?= /Users/pmoesta/codes/ecp/lapack-release
ifdef OLCF_NETLIB_LAPACK_ROOT
  LAPACK_HOME := $(OLCF_NETLIB_LAPACK_ROOT)/lib64
endif
ifdef OLCF_ESSL_ROOT
  INCLUDE_LOCATIONS += $(OLCF_ESSL_ROOT)/include
  LIBRARY_LOCATIONS += $(OLCF_ESSL_ROOT)/lib64
  LIBRARIES += -lessl
endif

USE_THORNADO = TRUE
THORNADO_HOME ?= $(CASTRO_HOME)/thornado
MOMENT_CLOSURE ?= MAXIMUM_ENTROPY_CB
NEUTRINO_MATTER_SOLVER ?= FIXED_POINT_COUPLED

USE_WEAKLIB = TRUE
WEAKLIB_HOME ?= $(CASTRO_HOME)/weaklib

# define the location of the CASTRO top directory
CASTRO_HOME ?= ../../..
MICROPHYSICS_HOME ?= ../../../../Microphysics

ifdef MICROPHYSICS_HOME

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_DIR     := weaklib

# This sets the network directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := nuclear_matter

else

$(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)

endif

# This sets the network directory in Castro/Networks

Bpack   := ./Make.package 

include $(CASTRO_HOME)/Exec/Make.Castro
