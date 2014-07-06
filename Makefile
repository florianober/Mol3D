################################################################
#                Makefile for mol3d based on mc3dv 4.x
#						  J. Sauter 2010 & F. Ober 2012
################################################################
#
# Usage examples:
#
# 	Build from scratch
#       > make new
#
#	Build for debugging
#	> make CO=debug
#
#	Build fast code
#	> make CO=fast
#
#	Build with Intel(R) Fortran Compiler [default = intel fortran]
#	> make FC=ifort
#
#	Build parallel on n cores:
#	> make -j n
#
#	And, of course, combinations thereof:
#	> make CO=fast FC=gfortran new -j 2
#
#	Get rid of everything build
#	> make clobber
#	
#	Also the directory of source files can be adjusted
#	> make SRC_DIR=othersrc/ new
#
#	Redo the last make 
#	> make likelast
#
################################################################
# 
# Defaults
#

#~ FC = gfortran
FC = ifort
RM = rm
MAKE = make
AWK = awk
DEPGEN = makedeps 
DEP_FILE = dependencies.dep
CO = normal
#
MAIN_EXE = mol3d
SRC_DIR = ./src/
BUILD_DIR = ./build/
LOG=$(addprefix $(BUILD_DIR),last_build)
#
#
# compiler dependet flags
#
ifeq ($(FC),g95)  # just for testing
  CFLAGS = -march=nocona -ffast-math -funroll-loops -O3
  DEPFLAGS = -I$(BUILD_DIR) -I$(BUILD_DIR)
  OFLAGS = 
endif

ifeq ($(FC),gfortran)
  CFLAGS = -Wall -O2
  DEPFLAGS = -J$(BUILD_DIR) -I$(BUILD_DIR) -lcfitsio
  OFLAGS = -x f95-cpp-input
  #OFLAGS = 
  ifeq ($(CO),debug)
    CFLAGS = -pg -O3 -fbounds-check -pedantic -Wall
  endif
  ifeq ($(CO),paradebug)
    CFLAGS = -march=native -ffast-math -funroll-loops -O3 -fopenmp -Wall -fbounds-check -pedantic
  endif
  ifeq ($(CO),para)
    CFLAGS = -march=native -ffast-math -funroll-loops -O3 -fopenmp -Wall
  endif
  ifeq ($(CO),fast)
    CFLAGS = -march=native -mtune=native -ffast-math -funroll-loops -O3
  endif
endif
#
ifeq ($(FC),ifort)
  CFLAGS = -O3 #-warn all -warn errors
  DEPFLAGS = -module $(BUILD_DIR) -I$(BUILD_DIR) -lcfitsio
  OFLAGS = -fpp
  ifeq ($(CO),debug)
#    CFLAGS = -pg -check bounds -check uninit -std -warn all -warn errors -WB -zero -traceback
#    CFLAGS = -pg -check bounds -check uninit -std -warn all -WB -zero -traceback
    CFLAGS = -pg -O3 -check bounds -check uninit -std -warn all -warn nodeclarations -WB -zero -traceback
  endif
  ifeq ($(CO),fast)
    #CFLAGS = -fast
    CFLAGS = -ip -ipo -O3 -axAVX
    #~ CFLAGS = -O3 -axAVX
  endif
  ifeq ($(CO),para)
    CFLAGS = -ip -ipo -O3 -axAVX -openmp
  endif
  ifeq ($(CO),paradebug)
    CFLAGS = -openmp -pg -check bounds -check uninit -std -warn all -warn nodeclarations -WB -zero -traceback
  endif
endif
#
#
# The generic target that makes everything
#
all: $(MAIN_EXE)
#
#
# A target to delete all object files
#
clean:
	@echo "Cleaning building directory"
	@$(RM) -rf $(BUILD_DIR)/*mod &> /dev/null
	@$(RM) -rf $(BUILD_DIR)/*o &> /dev/null
	@$(RM) -rf $(BUILD_DIR)/*dep &> /dev/null
#
#
# A target to delete everything built
#
clobber: clean
	@echo "Removing executable"
	@$(RM) -rf $(MAIN_EXE) &> /dev/null
	@$(RM) -rf $(BUILD_DIR)/last_build &> /dev/null
#
#
# A target to enforce building everything from scratch
#
new:	clobber
	$(MAKE)
#
#
# A target to build exaclty like the last time
# (relies on the file last.build in the BUILD_DIR)
#
likelast:
	@ a1=$(shell grep "Type" $(LOG) | awk '{print $$3}');\
	  a2=$(shell grep "Compiler" $(LOG) | awk '{print $$2}');\
	  a3=$(shell grep "Source" $(LOG) | awk '{print $$3}');\
	  make CO=$$a1 FC=$$a2 SRC_DIR=$$a3;
#
#
# A target the helps with the Makefile
#
help:
	@head -n 35 Makefile
#
#
##################################################################################
#
#
#  Definition of source files (their order does not matter)
#
MAIN_SRC = mol3d.f90
#
MODULE_SRC = 	basic_type.f90 \
		datatype.f90 \
		math_mod.f90 \
		var_globalnew.f90 \
		common_type.f90 \
		randgen_type.f90 \
		fluxes_type.f90 \
		model_type.f90 \
		grid_type.f90 \
		gas_type.f90 \
		dust_type.f90 \
		photon_type.f90 \
		interact_mod.f90 \
		immediate_mod.f90 \
		scatter_mod.f90 \
		initiate.f90 \
		tools_mod.f90 \
		string_mod.f90 \
		model_mod.f90 \
		source_type.f90 \
		grd_mod.f90 \
		fileio.f90 \
		simulation_mod.f90 \
		transfer_mod.f90 \
		lvlpop_mod.f90 \
		temp_mod.f90 \
		error_mod.f90 \
		start_mod.f90 \
		linkedlist_mod.f90 \
		parser_mod.f90
#
# The names of object files is derived from the source files
#
MODULE_OBJ = $(MODULE_SRC:f90=o)
MAIN_OBJ = $(MAIN_SRC:f90=o)
#
#
# Target for the main executable file
#
$(MAIN_EXE): $(addprefix $(BUILD_DIR),$(MODULE_OBJ)) $(addprefix $(BUILD_DIR),$(MAIN_OBJ)) 
ifeq ($(CO),perf) 
		@scalasca -instrument $(FC) $(CFLAGS) $(DEPFLAGS) -o $@  $^
else
		@$(FC) $(CFLAGS) $(DEPFLAGS) -o $@  $^
endif
	@date > $(LOG)
	@echo "Compiler:  " $(FC) >> $(LOG)
	@echo "Build Type:" $(CO) >> $(LOG)
	@echo "CFLAGS:    " $(CFLAGS) >> $(LOG)
	@echo "DEPFLAGS:  " $(DEPFLAGS) >> $(LOG)
	@echo "Source dir:" $(SRC_DIR) >> $(LOG)
#
#
# How to build an object file from an f90 source file
#
$(addprefix $(BUILD_DIR),%.o): $(addprefix $(SRC_DIR),%.f90) Makefile
ifeq ($(CO),paradebug) 
		@scalasca -instrument $(FC) $(CFLAGS) $(DEPFLAGS) $(OFLAGS) -c $< -o $@ 
else 
		@$(FC) $(CFLAGS) $(DEPFLAGS) $(OFLAGS) -c $< -o $@
endif
#
#
# A target that determines all dependencies
#
$(addprefix $(BUILD_DIR),$(DEP_FILE)): $(addprefix $(SRC_DIR),$(MODULE_SRC)) $(addprefix $(SRC_DIR),$(MAIN_SRC))
	@$(AWK) -f $(DEPGEN) builddir=$(BUILD_DIR) $^ > $@
$(addprefix $(SRC_DIR),$(MODULE_SRC)):

#
#
# Include all dependencies. GNU Make always checks if the included dependency file is up-to-date
# before anything else is done. If necessary, the dependency file is updated (by the target above)
# And the Makefile is rerun.
#
-include $(addprefix $(BUILD_DIR),$(DEP_FILE))
