# Makefile for BUBBLYSPH
# 

# OPTIONS:
# dims        3D (3) or 2D (2) simulation                             default: 2
# kernel      Quintic spline (1) or Wendland C2 (2)                   default: 2
# precision   Double (2) or single (1) precision                      default: 1
# bubbles     Include bubbles (1) or don't (0)                        default: 0
# zwalls      Wall BCs in 3rd dimension (1) or periodic (0)           default: 0
# slipwalls   Wall BCs are slip (1) or no-slip (0)                    default: 0
# restart     Restart from DUMP file (1) or start from scratch (0)    default: 0
# srs         Use turbulence closure model (1) or don't (0)           default: 1
# isotcomp    Isothermal compressibility (1) or incompressible (0)    default: 1
# wavemaker   A wavemaker at left and damping at right (1) or not (0) default: 0
# in3d        Take 3D input files (1) or don't (0)                    default: 0
# schwaiger   Use Schwaiger operator (1) or don't (0)                 default: 0
# acoustics   Calculate acoustic signal (1) or don't (0)              default: 0

# EXAMPLE USAGE:
# make dims=2 bubbles=0 slipwalls=1 srs=0 wavemaker=1 isotcomp=1


FC := gfortran
LD := gfortran
CFLAGS := -Wall -O3 -g -m64
FFLAGS :=-fopenmp -fbounds-check -ffpe-trap=zero -O3 -Wall -g -J./obj -I./obj -m64
## dimensions: default is 2
ifeq ($(dims), 3)
FFLAGS += -Ddim3
endif

## Schwaiger operator?
ifneq ($(schwaiger),0)
FFLAGS += -Dschwaiger
endif

## Kernel: default is 2 (Wendland)
ifeq ($(kernel),1)
FFLAGS += -Dkernel=1
else
FFLAGS += -Dkernel=2
endif
# precision level: default is 1, single (set precision=2 for double)
ifneq ($(precision), 2)
FFLAGS += -Dsingle
endif
## include bubbles?: default is no, set bubbles=1 for bubbles
ifeq ($(bubbles),1)
FFLAGS += -Dbubbles
endif
## 3rd dimension walls: default is periodic. set thirdwalls=1 for walls in Z-direction
ifeq ($(zwalls),1)
FFLAGS += -Dzwall
endif
## slipwalls: default is no
ifeq ($(slipwalls),1)
FFLAGS += -Dslip
endif
## Restart
ifeq ($(restart),1)
FFLAGS += -Drestart
endif
## Turbulence closure model
ifneq ($(srs),0)
FFLAGS += -Dsrs
endif
## Isothermal compressibility
ifneq ($(isotcomp),0)
FFLAGS += -Disotcomp
endif
## Wavemaker
ifeq ($(wavemaker),1)
FFLAGS += -Dwavemaker
endif
## Take 3D input files
ifeq ($(in3d),1)
FFLAGS += -Din3d
endif
## Acoustics
ifeq ($(acoustics),1)
FFLAGS += -Dacoustics
endif

LDFLAGS := -fopenmp -m64

SUB_DIRS := para base
SRC_DIR  := $(addprefix source/,$(SUB_DIRS))

#common_parameter needs to be first as mod file is depended upon by nearly everything.
OBJ_FILES := obj/kind_parameters.o obj/common_parameter.o obj/common_2d.o obj/cell_griddin.o 
OBJ_FILES += obj/bubble_acoustics.o
OBJ_FILES += obj/sphtools.o obj/input.o 
OBJ_FILES += obj/freesurf_mod.o obj/mirror_boundaries_mod.o 
OBJ_FILES += obj/bubble_boundaries.o
OBJ_FILES += obj/gradient_calculation.o obj/neighbour_finding_mod.o obj/shifting.o obj/linear_solver.o obj/bubble_evolution.o
OBJ_FILES += obj/svd_lib.o
OBJ_FILES += obj/predictor_mod.o obj/output.o
OBJ_FILES += obj/evolve.o
OBJ_FILES += $(foreach sdir,$(SRC_DIR),$(patsubst $(sdir)/%.F90,obj/%.o,$(wildcard $(sdir)/*.F90)))

HDEPS := $(OBJ_FILES:.o=.d)

vpath %.F90 $(SRC_DIR)

#-------


default: bubblysph
bubblysph: $(OBJ_FILES)
	$(LD) -o $@ $^ $(LDFLAGS)

obj/%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $<

-include $(HDEPS)


clean:
	rm -vf ./obj/*.o
	rm -vf ./obj/*.mod
	rm -vf ./bubblysph
	rm -vf ./DUMP
	rm -vf ./data_directory/statistics/*
	rm -vf ./data_directory/PART*
	rm -vf ./data_directory/BUBBLE*
	rm -vf ./data_directory/SLICE*
	rm -vf ./data_directory/WGRID*
	rm -vf ./data_directory/UGRID*	
	rm -vf ./data_directory/WSURF*	
	rm -vf ./data_directory/WELEV*
	rm -vf ./paraview_files/*
	rm -vf ./fort.*
	rm -vf ./mean*
	rm -vf ./u_r*	


