# Makefile

# CF90: Fortran compiler
# Available options:
#   * GFORT (serial version only)
#   * IFORT (serial & parallel & debug version)
#   * MPIF90 (parallel version only)
CF90 = GFORT

# MPI_PARALLEL: Variable for parallel processing (only for some compilers)
# Available options:
#   * YES
#   * NO
MPI_PARALLEL = NO

# OMP_THREADS: Use OpenMP for threading
# Available options:
#   * DYN: wave function parallelization (threads)
#   * YES: use threads for FFT
#   * NO
OMP_THREADS = NO

# OpenMP debugging, just for me... To be deleted after testing...
OMP_DEBUG = NO

# TYPE_FFT: FFT solvers
# Available options:
#   * NETLIB
#   * FFTW
#   * MKL
TYPE_FFT = FFTW

# MKL_THREADS: Enable OpenMP threading in MKL (only used if TYPE_FFT=MKL)
# Available options:
#   * YES
#   * NO
MKL_THREADS = YES

# DEBUG: enable debugging
# Available options:
#   * YES
#   * NO
DEBUG = NO

# LINK_STATIC: select static linkage of the binary
# Available options:
#   * YES
#   * NO
LINK_STATIC = NO

### SET MAIN LIBRARY LOCATIONS
# Set the location to the root of the FFTW subroutine package
FFTW_PATH = $(HOME)/Developer/fftw-3.3.8-gnu

# Set the location to the root of Intel's MKL distribution
MKL_PATH = $(MKLROOT)

# Set the location to the MKL FFTW wrappers in case MKL is used
MKLFFTWW_PATH = $(HOME)/Developer/fftw3xf

#####################################################################
#               Code selection										#
#####################################################################

# Code version
CODE_VERSION = 74
# grid: FFT or finite difference (FFT options must be NO for finite differences)
FINDIFF = NO
NUMEROV = NO
COUDOUB = YES
COUDOUB3D = NO
# Switch to extended model with polarizable raregas
RAREGAS = NO
# Include private sector
EXTENDED = NO
FSIC = NO

CODE_OPTIONS = '-DIVERSION=$(CODE_VERSION)'
ifeq ($(FINDIFF),YES)
	CODE_OPTIONS += -Dfindiff
endif
ifeq ($(NUMEROV),YES)
    CODE_OPTIONS += -Dnumerov
endif
ifeq ($(COUDOUB),YES)
	CODE_OPTIONS += -Dcoudoub
endif
ifeq ($(COUDOUB3D),YES)
    CODE_OPTIONS += -Dcoudoub3D
endif
ifeq ($(RAREGAS),YES)
    CODE_OPTIONS += -Draregas
endif
ifeq ($(EXTENDED),YES)
    CODE_OPTIONS += -Dextended
endif
ifeq ($(FSIC),YES)
    CODE_OPTIONS += -Dfsic
endif

#####################################################################
#                Process parallelization options                    #
#####################################################################

# MPI
ifeq ($(MPI_PARALLEL), YES)
	USE_MPI = YES
else
	USE_MPI = NO
endif

# OpenMP
ifneq ($(OMP_THREADS), NO)
	OMP = YES
else
	OMP = NO
endif

#####################################################################
#                     Parallelisation activation 	                #
#####################################################################

OMPADDL =

ifeq "$(CF90)" "MPIF90"
	COMPILER = mpif90
endif

ifeq "$(CF90)" "GFORT"
	COMPILER = gfortran
	ifeq "$(OMP)" "YES"
		OMPADDL = -fopenmp
	endif
endif

ifeq "$(CF90)" "IFORT"
	COMPILER = ifort
	ifeq "$(USE_MPI)" "YES"
		LDLIBS = -lmpi
  	endif
 	ifeq "$(OMP)" "YES"
		OMPADDL = -qopenmp
  	endif
endif

#####################################################################
#            Compiler- and subroutine-dependent options             #
#####################################################################

#OPT2: setting for the FFT package, needs forced double precision
#OPT3: setting for critical soubroutines which do not stand optimization

ifeq "$(CF90)" "GFORT"
	OPT1 = $(OMPADDL) -w -O3 -mfpmath=sse -fdefault-real-8 -fdefault-double-8
	OPT2 = $(OMPADDL) -w -O3 -mfpmath=sse -fdefault-real-8 -fdefault-double-8
	OPT3 = $(OMPADDL) -w -g -fdefault-real-8 -fdefault-double-8
	ifeq "$(DEBUG)" "YES"
		OPT1 = $(OMPADDL) -pg -w -g -fbacktrace -fdefault-real-8 -fdefault-double-8 -fcheck=bounds
		OPT2 = $(OMPADDL) -pg -w -g -fbacktrace -fdefault-real-8 -fdefault-double-8 -fcheck=bounds
		OPT3 = $(OMPADDL) -pg -w -g -fbacktrace -fdefault-real-8 -fdefault-double-8 -fcheck=bounds
  	endif
endif

ifeq "$(CF90)" "IFORT"
	STATIC = -static
	OPT1 = $(OMPADDL) -fpp -w -xW -O3 -ip -no-prec-div -align all -autodouble
	OPT2 = $(OMPADDL) -fpp -w -xW -O3 -ip -no-prec-div -align all -autodouble
	OPT3 = $(OMPADDL) -fpp -w -g -align all -autodouble
	ifeq "$(DEBUG)" "YES"
		OPT1 = $(OMPADDL) -pg -fpp -w -g -CB -traceback -align all -autodouble
		OPT2 =$(OMPADDL) -pg -fpp -w -g -CB -traceback -align all -autodouble
		OPT3 =$(OMPADDL) -pg -fpp -w -g -align all -autodouble
	endif
endif

ifeq "$(CF90)" "MPIF90"
	STATIC = -static
	OPT1 = -fpp -w -O2 -no-prec-div -align all -autodouble
	ifeq "$(FCOMP)" "gfortran"
		OPT1 = -w -O3 -mfpmath=sse -fdefault-real-8 -fdefault-double-8
	endif
	OPT2 = $(OPT1)
	OPT3 = -fpp -w -g -align all -autodouble
	ifeq "$(DEBUG)" "YES"
		OPT1 = -pg -fpp -w -g -CB -traceback -align all -autodouble
		OPT2 = $(OPT1) -autodouble
		OPT3 = -pg -fpp -w -g -align all -autodouble
	endif
endif

#####################################################################
#  Final pre-processor and compiler flags for MPI/OMP/FFT and MKL   #
#####################################################################

# (MKL is used to provide FFTW-support through its wrappers)
ifeq ($(TYPE_FFT), FFTW)
	calls_FFTW = YES
else ifeq ($(TYPE_FFT), MKL)
	calls_FFTW = YES
else
	calls_FFTW = NO
endif

# Switches used to preprocess files.
ifeq ($(calls_FFTW), YES)
	FFTWADD = -Dfftw_cpu
	fftw_cpu_value := 1
	ifeq ($(TYPE_FFT), FFTW)
		FFTWADD += -Dfftwnomkl
		fftwnomkl_value := 1
	endif
else ifeq ($(TYPE_FFT), NETLIB)
	FFTWADD = -Dnetlib_fft
	netlib_fft_value := 1
endif

ifeq ($(MPI_PARALLEL), YES)
	MPIADD := -Dparayes
	parayes_value := 1
else
	MPIADD = -Dparano
	parano_value := 1
endif

OMPADD =
ifeq ($(OMP), YES)
	OMPADD = -Dparopenmp
	paropenmp_value := 1
	ifeq ($(OMP_THREADS), DYN)
		OMPADD += -Ddynopenmp
		dynopenmp_value :=1
	endif
endif

ifeq ($(OMP_DEBUG), YES)
	OMPADD += -Domp_debug
endif

# Final compiler flags
COMPILERFLAGS1 = $(LDLIBS) $(OPT1) $(strip $(MPIADD) $(OMPADD) $(FFTWADD)) $(CODE_OPTIONS)
COMPILERFLAGS2 = $(LDLIBS) $(OPT2) $(strip $(MPIADD) $(OMPADD) $(FFTWADD)) $(CODE_OPTIONS)
COMPILERFLAGS3 = $(LDLIBS) $(OPT3) $(strip $(MPIADD) $(OMPADD) $(FFTWADD)) $(CODE_OPTIONS)

#####################################################################
#                Linker configuration for FFT and MKL               #
#####################################################################

LINKER = $(COMPILER)

### DERIVED FROM MAIN LIBRARY LOCATIONS AT TOP OF MAKEFILE
FFTW_LIBDIR = $(FFTW_PATH)/lib
FFTW_INCLUDE = $(FFTW_PATH)/include
MKL_LIBDIR = $(MKL_PATH)/lib
MKL_INCLUDE = $(MKL_PATH)/include/fftw

# Flags for linking using the MKL libraries.
ifeq ($(TYPE_FFT), MKL)
			
	MKLFLAGS =
	FFTWFLAGS =

	ifneq "$(MKLFFTWW_PATH)" ""
		MKLFLAGS += -L$(MKLFFTWW_PATH)
  	endif
	
	ifeq ($(CF90), IFORT)
		MKLFLAGS += -lfftw3xf_intel
	else
		MKLFLAGS += -lfftw3xf_gnu
	endif

 	ifneq "$(MKL_LIBDIR)" ""
		MKLFLAGS += -L$(MKL_LIBDIR)
	endif
  
	ifneq "$(MKL_INCLUDE)" ""
		MKLFLAGS += -I$(MKL_INCLUDE)
	endif

	ifeq ($(LINK_STATIC), YES) # Static MKL => group flags must be passed to the linker
		MKLFLAGS += -Wl,--start-group
  	endif

	MKLFLAGS += -lmkl_intel_lp64

	ifeq ($(MKL_THREADS), YES) # OpenMP threading support in MKL: Threading layer
		MKLFLAGS += -lmkl_intel_thread
  	else
		MKLFLAGS += -lmkl_sequential
	endif

	MKLFLAGS += -lmkl_core

	ifeq ($(LINK_STATIC), YES) # Static MKL
		MKLFLAGS += -Wl,--end-group
	endif

	ifeq ($(MKL_THREADS), YES) # OpenMP threading support in MKL: Run-time lib
		MKLFLAGS += -liomp5
	endif

	FFTWFLAGS = $(MKLFLAGS)

# Flags for directly linking the original FFTW libraries.
else ifeq ($(TYPE_FFT), FFTW)
	FFTWFLAGS += -L$(FFTW_LIBDIR) -lfftw3
	ifeq ($(OMP), YES)
		FFTWFLAGS += -lfftw3_omp
	endif
endif

LINKERFLAGS = $(LDLIBS) $(FFTWFLAGS) $(OPT1)

ifeq ($(LINK_STATIC), YES)
	LINKERFLAGS += $(STATIC)
endif

#####################################################################
#                          Executable name                          #
#####################################################################

EXEC = qdd

#####################################################################
#                 Phony- and default targets                        #
#####################################################################

.PHONY: all clean distclean cleanall cleanobj cleanmod # NO TRUE FILE TARGET PREREQUISITS MAY
											 # APPEAR HERE, UNLESS YOU WANT THEM TO BE
											 # REBUILT EVERY TIME!

.DEFAULT_GOAL := all # In this case the default target is already pointing to 'all'
					 # because it is setup to be the first target. However, if 'all'
					 # weren't to be the first target, this statement will enforce it to
					 # be still the default target.

all: $(EXEC)

clean: cleanall

cleanall: cleanobj cleanmod

distclean: cleanall
	@rm -vf ../../bin/$(EXEC)

cleanobj:
	@rm -vf *.o

cleanmod:
	@rm -vf *.mod

#####################################################################
#                             Checks                                #
#####################################################################

# List of valid compilers.
valid_CF90 := GFORT IFORT MPIF90

# List of valid MPI_PARALLEL:
valid_MPI_PARALLEL := YES NO

# List of valid OMP_THREADS:
valid_OMP_THREADS := DYN YES NO

# List of valid TYPE_FFT values.
valid_FFT := FFTW NETLIB MKL

#initial_checks: Makefile
initial_checks:
	@echo "##############################################################################"
	@echo "Lists of known compilation parameters."
	@echo "##############################################################################"
	@echo "Parameters defined in Makefile (Solvers and extensions):"
	@echo "------------------------------------------------------------------------------"
	@echo "findiff       $(FINDIFF)"
	@echo "numerov       $(NUMEROV)"
	@echo "coudoub       $(COUDOUB)"
	@echo "coudoub3D     $(COUDOUB3D)"
	@echo "raregas       $(RAREGAS)"
	@echo "##############################################################################"
	@echo "Parameters defined in Makefile (Parallelization, includes and libs):"
	@echo "------------------------------------------------------------------------------"
	@echo "CF90          $(CF90)"
	@echo "MPI_PARALLEL  $(MPI_PARALLEL)"
	@echo "OMP_THREADS   $(OMP_THREADS)"
	@echo "TYPE_FFT      $(TYPE_FFT)"
	@echo "DEBUG         $(DEBUG)"
	@echo "LINK_STATIC   $(LINK_STATIC)"
ifneq ($(TYPE_FFT), NETLIB)
ifeq ($(TYPE_FFT), MKL)
	@echo "MKL_THREADS   $(MKL_THREADS)"
	@echo "MKL_PATH      $(MKL_PATH)"
	@echo "MKLFFTWW_PATH $(MKLFFTWW_PATH)"
	@echo "MKL_LIBDIR    $(MKL_LIBDIR)"
	@echo "MKL_INCLUDE   $(MKL_INCLUDE)"
else
	@echo "FFTW_PATH     $(FFTW_PATH)"
	@echo "FFTW_LIBDIR   $(FFTW_LIBDIR)"
	@echo "FFTW_INCLUDE  $(FFTW_INCLUDE)"
endif
endif
	@echo "##############################################################################"
	@echo "-D flags enabled (1) by Makefile rules:"
	@echo "------------------------------------------------------------------------------"
	@echo "parayes       $(parayes_value)"
	@echo "parano        $(parano_value)"
	@echo "netlib_fft    $(netlib_fft_value)"
	@echo "fftw_cpu      $(fftw_cpu_value)"
	@echo "fftwnomkl     $(fftwnomkl_value)"
	@echo "paropenmp     $(paropenmp_value)"
	@echo "dynopenmp     $(dynopenmp_value)"
	@echo "##############################################################################"

# Check that the fortran compiler is defined.
ifneq ($(CF90), $(filter $(CF90), $(valid_CF90)))
	$(error ERROR: Unknown fortran compiler in Makefile (CF90=$(CF90)))
endif

# Check that NETLIB and OpenMP are not both active.
ifeq ($(TYPE_FFT), NETLIB)
  ifeq ($(OMP), YES)
    $(error ERROR: NETLIB does not support OpenMP threading. Please use either FFTW or MKL, or disable OpenMP threading)
  endif
endif

# Check that the TYPE_FFT is valid:
ifneq ($(TYPE_FFT), $(filter $(TYPE_FFT), $(valid_FFT)))
	$(error ERROR: Unknown TYPE_FFT in Makefile ($(TYPE_FFT)))
endif

# Check that the parallelization and threading options are valid:
ifneq ($(MPI_PARALLEL), $(filter $(MPI_PARALLEL), $(valid_MPI_PARALLEL)))
	$(error ERROR: Unknown MPI_PARALLEL value in Makefile ($(MPI_PARALLEL)))
endif

ifneq ($(OMP_THREADS), $(filter $(OMP_THREADS), $(valid_OMP_THREADS)))
	$(error ERROR: Unknown OMP_THREADS value in Makefile ($(OMP_THREADS)))
endif

# Check that variables which should be either YES or NO have proper values:
ifneq ($(DEBUG),       $(filter $(DEBUG),       YES NO))
	$(error ERROR: Unknown $(strip DEBUG       value) in Makefile ($(DEBUG)))
endif

ifneq ($(LINK_STATIC), $(filter $(LINK_STATIC), YES NO))
	$(error ERROR: Unknown $(strip LINK_STATIC value) in Makefile ($(LINK_STATIC)))
endif

ifneq ($(MKL),         $(filter $(MKL)        , YES NO))
	$(error ERROR: Unknown $(strip MKL         value) in Makefile ($(MKL)))
endif

# coudoub3D is only enabled when using FFTW calls.
ifneq ($(COUDOUB3D),NO)
  ifneq ($(calls_FFTW), YES)
	$(error ERROR: coudoub3D requires the use of FFTW)
  endif
endif

# Check collision of MPI and OpenMP:
ifeq ($(USE_MPI), YES)
  ifeq ($(OMP), YES)
	$(error ERROR: code not ready for simultaneous MPI and OpenMP parallelizations)
  endif
endif

	@echo "Done with the initial checks, starting compilation..."
	@echo "#######################################################"

#####################################################################
#                  Targets and common dependencies                  #
#####################################################################

OBJINT = params.o main.o kinetic.o restart.o restartc.o init.o\
       static.o dynamic.o lda.o util.o abso_bc.o\
       pseudosoft.o pseudogoed.o ionmd.o forces.o\
       carlo.o localize.o localizer.o\
       sicnew.o sicnewc.o rho.o rhoc.o nonloc.o nonlocc.o\
       schmid.o zeroforce.o loc_mfield.o givens.o subgrids.o\
       parallele.o rta.o\
       HEeigensystem.o mini.o zdiag.o orthmat.o

ifeq ($(EXTENDED),YES)
	OBJINT += util_extended.o attachement.o
endif

ifeq ($(FSIC),YES)
	OBJINT += 2stUT.o 2st_util.o 2stUTc.o
endif       

ifeq ($(RAREGAS), YES)
	OBJINT +=  functions.o pot_substrate.o forces_substrate.o md_substrate.o\
		short.o image.o lattice.o util_substrate.o
endif

# List of modules:
MODLIST = params.mod kinetic.mod
OBJINT += coulsolv.o
MODLIST += coulsolv.mod

# Object file for FFT:
ifeq ($(TYPE_FFT), NETLIB)
	FFT_OBJ = fftpack.o
endif

ifeq ($(FSIC),YES)
	MODLIST += 2st_util.mod twostr.mod twost.mod
endif

OBJS = $(OBJINT) $(FFT_OBJ)

# List of objects excluding params.o:
OBJS_NPAR = $(filter-out params.o, $(OBJS))

# Objects that depend on kinetic.mod and/or coulsolv.mod,
# with no other special requirements:
OBJS_KIN  = dynamic.o ionmd.o util.o zeroforce.o
OBJS_COUL = image.o lda.o loc_mfield.o
OBJS_KINCOUL = init.o main.o pseudosoft.o
DEPS_KINCOUL = kinetic.mod coulsolv.mod

### I AM NOT SURE WHAT THIS IS FOR. IT JUST SEEM THAT EVERYTHING GETS REBUILT WHENEVER THE
### THE MAKEFILE IS TOUCHED
# Dependencies on Makefile and params.mod
# Dependencies are inherited, so there's no need to repeat them
# in the targets below.
# $(MODLIST):   Makefile 
# params.o:     Makefile 
# $(OBJS_NPAR): Makefile params.mod

# Dependencies for objects that depend on kinetic.o and/or coulsolv.o,
# with no other special requirements:
$(OBJS_KIN):     kinetic.mod
$(OBJS_COUL):    coulsolv.mod
$(OBJS_KINCOUL): $(DEPS_KINCOUL)

# Set FFTW-dependent modules:
ifeq ($(calls_FFTW), YES)
coulsolv.mod: fftw.mod
kinetic.mod: fftw.mod
endif

### DEBUG INFO
# $(info $$FFT_OBJ is [${FFT_OBJ}])
# $(info )
# $(info $$OBJINT is [${OBJINT}])
# $(info )
# $(info $$OBJS is [${OBJS}])
# $(info )
# $(info $$OBJS_NPAR is [${OBJS_NPAR}])
# $(info )
# $(info $$OBJS_KIN is [${OBJS_KIN}])
# $(info )
# $(info $$OBJS_COUL is [${OBJS_COUL}])
# $(info )
# $(info $$OBJS_KINCOUL is [${OBJS_KINCOUL}])
# $(info )
# $(info $$DEPS_KINCOUL is [${DEPS_KINCOUL}])
# $(info )

#####################################################################
#                   Compilation and linkage rules                   #
#####################################################################

$(EXEC): initial_checks $(OBJS)
	@echo Linking executable $@
	$(LINKER) -o $@ $(strip $(OBJS) $(LINKERFLAGS))
	mv -fv $(EXEC) ../../bin/

# Implicit rule for objects and modules together:
%.o %.mod: %.F90
#	@echo " "
#	@echo "IMPLICIT compilation selected for object or module:" $<
	$(COMPILER) $(COMPILERFLAGS1) -c $<
#	@echo " "

# Additions for twostsic:
ifeq ($(FSIC),YES)
FSICDEPEND = 2st_util.o 2stUTc.o 2stUT.o

2st_util.o 2st_util.mod: fsic/2st_util.F90 kinetic.mod orthmat.o
# 	@echo " "
# 	@echo "EXPLICIT compilation for file..."
	$(COMPILER) $(COMPILERFLAGS1) -DREALSWITCH -o 2st_util.o -c $<
# 	@echo " "

2stUTc.o twost.mod: fsic/2stUT.F90 kinetic.mod twostr.mod
# 	@echo " "
# 	@echo "EXPLICIT compilation for file..."
	$(COMPILER) $(COMPILERFLAGS1) -DCOMPLEXSWITCH -o 2stUTc.o -c $<
# 	@echo " "

2stUT.o twostr.mod: fsic/2stUT.F90 kinetic.mod 2st_util.mod
# 	@echo "EXPLICIT compilation for file..."
	$(COMPILER) $(COMPILERFLAGS1) -DREALSWITCH -o 2stUT.o -c $<
endif

# Explicit rules
# (complex dependencies, REAL/COMPLEX switch, non-default COMPILERFLAGS).
kinetic.o : kinetic.F90 fft.F90 findiff/findiff.F90 $(FFT_OBJ)

coulsolv.o: coulsolv.F90 falr.o coulex.o kinetic.mod $(FFT_OBJ)

static.o: static.F90 pseudosoft.F90 $(DEPS_KINCOUL) $(FSICDEPEND)

dyanamic.o: dynamic.F90 pseudosoft.o twostr.o twost.o $(DEPS_KINCOUL) $(FSICDEPEND)

rho.o: rho.F90
	$(COMPILER) $(COMPILERFLAGS1) -DREALSWITCH    -o $@ -c $<

rhoc.o: rho.F90 kinetic.mod
	$(COMPILER) $(COMPILERFLAGS1) -DCOMPLEXSWITCH -o $@ -c $<

localizer.o: localize.F90 kinetic.mod
	$(COMPILER) $(COMPILERFLAGS1) -DREALSWITCH    -o $@ -c $<

localize.o: localize.F90 kinetic.mod
	$(COMPILER) $(COMPILERFLAGS1) -DCOMPLEXSWITCH -o $@ -c $<

sicnew.o: sicnew.F90 $(DEPS_KINCOUL)
	$(COMPILER) $(COMPILERFLAGS1) -DREALSWITCH    -o $@ -c $<

sicnewc.o: sicnew.F90 $(DEPS_KINCOUL)
	$(COMPILER) $(COMPILERFLAGS1) -DCOMPLEXSWITCH -o $@ -c $<

nonloc.o: nonloc.F90 kinetic.mod
	$(COMPILER) $(COMPILERFLAGS1) -DREALSWITCH    -o $@ -c $<

nonlocc.o: nonloc.F90 kinetic.mod
	$(COMPILER) $(COMPILERFLAGS1) -DCOMPLEXSWITCH -o $@ -c $<

fftw.o fftw.mod: fftw.F90
ifeq ($(TYPE_FFT), MKL)
	$(COMPILER) $(COMPILERFLAGS1) -I$(MKL_INCLUDE) -c $<
else ifeq ($(TYPE_FFT), FFTW)
	$(COMPILER) $(COMPILERFLAGS1) -I$(FFTW_INCLUDE) -c $<
endif

ifeq ($(TYPE_FFT), NETLIB)
fftpack.o: fftpack.F90 fftpack2.F90
	$(COMPILER) $(COMPILERFLAGS2) -c $<
endif

main.o: main.F90 util.o orthmat.o $(FSICDEPEND)
	$(COMPILER) $(COMPILERFLAGS3) -c $<

givens.o: givens.F90
	$(COMPILER) $(COMPILERFLAGS3) -c $<

restart.o: restart.F90 kinetic.mod $(FSICDEPEND)
	$(COMPILER) $(COMPILERFLAGS3) -DREALSWITCH    -o $@ -c $<

restartc.o: restart.F90 kinetic.mod $(FSICDEPEND)
	$(COMPILER) $(COMPILERFLAGS3) -DCOMPLEXSWITCH -o $@ -c $<

parallele.o: parallele.F90 kinetic.mod
	$(COMPILER) $(COMPILERFLAGS3) -c $<

rta.o: rta.F90 kinetic.mod params.mod
	$(COMPILER) $(COMPILERFLAGS3) -c $<

orthmat.o: orthmat.F90
	$(COMPILER) $(COMPILERFLAGS3) -c -cpp $<

zdiag.o: zdiag.f
	$(COMPILER) $(COMPILERFLAGS3) -c $<

mini.o: mini.f
	$(COMPILER) $(COMPILERFLAGS3) -c $<

HEeigensystem.o: HEeigensystem.f
	$(COMPILER) $(COMPILERFLAGS3) -c -cpp $<

ifeq ($(EXTENDED),YES)
attachement.o: extended/attachement.F90
	$(COMPILER) $(COMPILERFLAGS3) -c $<

util_extended.o: extended/util_extended.F90 params.o
	$(COMPILER) $(COMPILERFLAGS3) -c $<
endif

# QM/MM routines
ifeq ($(RAREGAS), YES)
forces_substrate.o: QMMM/forces_substrate.F90 params.mod
	$(COMPILER) $(COMPILERFLAGS1) -c $<

util_substrate.o: QMMM/util_substrate.F90 params.mod
	$(COMPILER) $(COMPILERFLAGS1) -c $<

functions.o: QMMM/functions.F90 params.mod 
	$(COMPILER) $(COMPILERFLAGS1) -c $<

pot_substrate.o: QMMM/pot_substrate.F90 params.mod
	$(COMPILER) $(COMPILERFLAGS1) -c $<

md_substrate.o: QMMM/md_substrate.F90 params.mod
	$(COMPILER) $(COMPILERFLAGS1) -c $<

short.o: QMMM/short.F90 params.mod
	$(COMPILER) $(COMPILERFLAGS1) -c $<

image.o: QMMM/image.F90 params.mod
	$(COMPILER) $(COMPILERFLAGS1) -c $<

lattice.o: QMMM/lattice.F90 params.mod
	$(COMPILER) $(COMPILERFLAGS1) -c $<

subgrids.o: subgrids.F90 params.mod
	$(COMPILER) $(COMPILERFLAGS1) -c $<
endif
