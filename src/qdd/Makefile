
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

# OMP_THREADS: Invoke OpenMP
# Available options:
#   * DYN: wave function parallelization (threads)
#   * YES: use threads for FFT
#   * NO
OMP_THREADS = YES

# TYPE_FFT: FFT solver
# Available options:
#   * NETLIB
#   * FFTW
#   * MKL
TYPE_FFT = FFTW

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

# MKL_THREADS: Enable MKL threading (only used for TYPE_FFT=FFTW)
# Available options:
#   * YES
#   * NO
MKL_THREADS = NO

# MKL_WRAPPERS: Path to the MKL FFTW wrappers (only used for TYPE_FFT=FFTW)
# Default value:
#MKL_WRAPPERS = /usr/local/FFTW/FFTW3/mkl/interfaces/fftw3xf
#for CCC
MKL_WRAPPERS = $(HOME)/Developer/fftw3xf

FFTW3_LIBDIR = /Users/fmgjcoppens/fftw3-icc/lib

# MKLPATH: Path to the MKL libraries (only used for TYPE_FFT=FFTW)
# Default value:
#MKLPATH = /usr/local/FFTW/FFTW3/mkl/lib/em64t
#for CCC
MKLPATH = $(MKLROOT)/lib
# Path for hyperion:
# MKLPATH = /divers/intel/Compiler/12.1/mkl/lib/intel64

# MKLINCLUDE: Include path for the MKL (only used for TYPE_FFT=FFTW )
# Default value:
#MKLINCLUDE = /usr/local/FFTW/FFTW3/mkl/include/fftw
#for CCC
MKLINCLUDE = $(MKLROOT)/include/fftw


# For large boxe sizes, try option
# '-i-dynamic -mcmodel=medium' or '-shared-intel -mcmodel=medium'

#####################################################################
#               Code selection					                    #
#####################################################################

# Code version
CODE_VERSION = 74
# grid: FFT or finite difference (FFT options must be NO for finite differences)
FINDIFF = NO
NUMEROV = NO
COUDOUB=YES
COUDOUB3D=NO
# Switch to extended model with polarizable raregas
RAREGAS = NO
# Include private sector
EXTENDED = YES
FSIC=NO

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
ifeq ($(MPI_PARALLEL), NO)
	USE_MPI = NO
else
	USE_MPI = YES
endif

# OpenMP
ifeq ($(OMP_THREADS), NO)
	OMP = NO
else
	OMP = YES
endif

#####################################################################
#                             Compiler                              #
#####################################################################

# Insert the run command for invoking your compiler
# ifeq "$(CF90)" "OWN"
# 	COMPILER  = ...
# 	LDLIBS    =
#   ifeq "$(OMP)" "YES"
# 	OMPADDL = ...
#   endif
# endif

OMPADDL =

ifeq "$(CF90)" "MPIF90"
#	COMPILER  = /opt/mpich2-intel/bin/mpif90
#	LDLIBS    =
	COMPILER  = mpif90
	LDLIBS    =
endif


ifeq "$(CF90)" "GFORT"
	COMPILER  = gfortran
	LDLIBS    =
  ifeq "$(OMP)" "YES"
	OMPADDL = -fopenmp
  endif
endif


ifeq "$(CF90)" "GFORT1"
	COMPILER  = gfortran
	LDLIBS    =
  ifeq "$(OMP)" "YES"
	OMPADDL = -fopenmp
  endif
endif


ifeq "$(CF90)" "IFORT"
	COMPILER  = ifort
	LDLIBS    =
  ifeq "$(USE_MPI)" "YES"
	LDLIBS    = -lmpi
  endif
  ifeq "$(OMP)" "YES"
	OMPADDL =  -qopenmp
  endif
endif


ifeq "$(CF90)" "IFORT1"
	COMPILER  = ifort
	LDLIBS    =
  ifeq "$(USE_MPI)" "YES"
	LDLIBS    = -lmpi
  endif
  ifeq "$(OMP)" "YES"
	OMPADDL = -qopenmp
  endif
endif


ifeq "$(CF90)" "XLF_R"
	COMPILER  = xlf90_r
	LDLIBS    =
  ifeq "$(USE_MPI)" "YES"
	COMPILER  = mpxlf90_r
	LDLIBS    = -lpesslsmp -lesslsmp
  endif
endif


#####################################################################
#                    Compiler-dependent options                     #
#####################################################################

#OPT2: setting for the FFT package, needs forced double precision
#OPT3: setting for critical soubroutines which do not stand optimization

# Insert your actual compiler her with its options.
# If wanted fill also the debugging options.
# ifeq "$(OWN)" "YES"
# 	OPT1  = ...
# 	OPT2  = ...
# 	OPT3  = ...
# ifeq "$(DEBUG)" "YES"
# 	OPT1  = ...
# 	OPT2  = ...
# 	OPT3  = ...
# endif
# endif

ifeq "$(CF90)" "GFORT"
	OPT1  =  $(OMPADDL) -w -O3 -mfpmath=sse -fdefault-real-8 -fdefault-double-8 -mcmodel=large
	OPT2  =  -w -O3 -mfpmath=sse -fdefault-real-8 -fdefault-double-8
	OPT3  =  -w -g  -fdefault-real-8 -fdefault-double-8
#	OPT1  = -cpp -w -O3 -march=k8-sse3 -mfpmath=sse
#	OPT2  = -cpp -w -O3 -march=k8-sse3 -mfpmath=sse -fdefault-real-8
#	OPT3  = -cpp -w -g
  ifeq "$(DEBUG)" "YES"
	OPT1  =  -pg -w -g -fbacktrace -fdefault-real-8 -fdefault-double-8 -fcheck=bounds
	OPT2  =  -pg -w -g -fbacktrace -fdefault-real-8 -fdefault-double-8 -fcheck=bounds
	OPT3  =  -pg -w -g -fbacktrace -fdefault-real-8 -fdefault-double-8 -fcheck=bounds
  endif
endif


ifeq "$(CF90)" "GFORT1"
	OPT1  = $(OMPADDL) -w -O3 -msse4.2 -mfpmath=sse -ffast-math -fdefault-real-8 -fdefault-double-8 -finline-functions -funroll-loops
	OPT2  = $(OPT1) -fdefault-real-8 -fdefault-double-8
	OPT3  =  -w -g  -fdefault-real-8 -fdefault-double-8
  ifeq "$(DEBUG)" "YES"
	OPT1  =  -pg -cpp -w -g -fbacktrace -fdefault-real-8 -fdefault-double-8
	OPT2  =  -pg -cpp -w -g -fbacktrace -fdefault-real-8 -fdefault-double-8
	OPT3  =  -pg -cpp -w -g
  endif
endif


ifeq "$(CF90)" "IFORT"
	STATIC = -static
#	OPT1  = -fpp -w -xW -O3 -ip -no-prec-div -align all -autodouble -static
#	OPT2  = -fpp -w -xW -O3 -ip -no-prec-div -align all -autodouble -static
#	OPT3  = -fpp -w  -pg      -g               -align all -autodouble -static
	OPT1  = $(OMPADDL) -fpp -w -xW -O3 -ip -no-prec-div -align all -autodouble
	OPT2  =  -fpp -w -xW -O3 -ip -no-prec-div -align all -autodouble
	OPT3  =  -fpp -w       -g                 -align all -autodouble
  ifeq "$(DEBUG)" "YES"
	OPT1  =  -pg -fpp -w -g -CB -traceback -align all -autodouble
	OPT2  =  -pg -fpp -w -g -CB -traceback -align all -autodouble
	OPT3  =  -pg -fpp -w -g                -align all -autodouble
  endif
endif


ifeq "$(CF90)" "IFORT1"
	STATIC = -static
	OPT1  = -fpp -w -axAVX -msse4.2 -O3 -ip -no-prec-div -align all -autodouble
	OPT2  = -fpp -w -axAVX -msse4.2 -O3 -ip -no-prec-div -align all -autodouble
	OPT3  = -fpp -w        -g                            -align all -autodouble
  ifeq "$(DEBUG)" "YES"
	OPT1  =  -fpp -w -g -CB -traceback -align all -autodouble
	OPT2  =   $(OPT1)
	OPT3  =  -fpp -w -g -align all -autodouble
  endif
endif


ifeq "$(CF90)" "MPIF90"
        STATIC = -static
        OPT1  =  -fpp -w -O2 -no-prec-div -align all -autodouble
	ifeq "$(FCOMP)" "gfortran"
		OPT1  =  -w -O3 -mfpmath=sse -fdefault-real-8 -fdefault-double-8
	endif
        OPT2  =  $(OPT1)
        OPT3  =  -fpp -w -g               -align all -autodouble
  ifeq "$(DEBUG)" "YES"
        OPT1  =  -pg -fpp -w -g -CB -traceback -align all -autodouble
        OPT2  =  $(OPT1) -autodouble
        OPT3  =  -pg -fpp -w -g                -align all -autodouble
  endif
endif


ifeq "$(CF90)" "XLF_R"
	OPT1  = -d -w -qport=mod -q64 -qarch=pwr6 -qtune=pwr6 -O3 -qstrict
	OPT2  = -d -w -qport=mod -q64 -qarch=pwr6 -qtune=pwr6 -O3 -qstrict -qdpc -qautodbl=dbl
	OPT3  = -d -w -qport=mod -q64 -qarch=pwr6 -qtune=pwr6
endif


#####################################################################
#                        Compilation rules                          #
#####################################################################

# MKL only used to provide FFTW support through wrappers.
ifeq ($(TYPE_FFT), FFTW)
	calls_FFTW = YES
else
  ifeq ($(TYPE_FFT), MKL)
	calls_FFTW = YES
  else
	calls_FFTW = NO
  endif
endif

# Switches used to preprocess files.
ifeq ($(calls_FFTW), YES)
	FFTWADD = -Dfftw_cpu
	fftw_cpu_value := 1
  ifeq ($(TYPE_FFT), FFTW)
	FFTWADD += -Dfftwnomkl
	fftwnomkl_value := 1
  endif
else
  ifeq ($(TYPE_FFT), NETLIB)
	FFTWADD = -Dnetlib_fft
	netlib_fft_value := 1
  endif
endif

ifeq ($(MPI_PARALLEL), YES)
	MPIADD := -Dparayes
	parayes_value := 1
else
	MPIADD = -Dparano
	parano_value := 1
endif

ifeq ($(OMP), YES)
	OMPADD = -Dparopenmp
	paropenmp_value := 1
  ifeq ($(OMP_THREADS), DYN)
	OMPADD += -Ddynopenmp
	dynopenmp_value :=1
  endif
else
	OMPADD :=
endif

# Compiler flags
COMPILERFLAGS1 = $(LDLIBS) $(OPT1) $(strip $(MPIADD) $(OMPADD) $(FFTWADD)) $(CODE_OPTIONS)
COMPILERFLAGS2 = $(LDLIBS) $(OPT2) $(strip $(MPIADD) $(OMPADD) $(FFTWADD)) $(CODE_OPTIONS)
COMPILERFLAGS3 = $(LDLIBS) $(OPT3) $(strip $(MPIADD) $(OMPADD) $(FFTWADD)) $(CODE_OPTIONS)

# This will also be used in the compilation (with -DREALSWITCH / -DCOMPLEXSWITCH)

#####################################################################
#                             Linker                                #
#####################################################################

LINKER = $(COMPILER)

# Flags for linking using the MKL libraries.
ifeq ($(TYPE_FFT), MKL)

	MKLFLAGS =
	FFTWFLAGS =


  ifneq "$(MKL_WRAPPERS)" ""
	MKLFLAGS += -L$(MKL_WRAPPERS)
  endif
	MKLFLAGS += -lfftw3xf_intel

  ifneq "$(MKLPATH)" ""
	MKLFLAGS += -L$(MKLPATH)
  endif
  ifneq "$(MKLINCLUDE)" ""
	MKLFLAGS += -I$(MKLINCLUDE)
  endif

  ifeq ($(LINK_STATIC), YES) # Static MKL => group flags must be passed to the linker
	MKLFLAGS += -Wl,--start-group
  endif

	MKLFLAGS += -lmkl_intel_lp64

  ifeq ($(MKL_THREADS), YES)
	MKLFLAGS += -lmkl_intel_thread
  else
	MKLFLAGS += -lmkl_sequential
  endif

	MKLFLAGS += -lmkl_core

  ifeq ($(LINK_STATIC), YES) # Static MKL
	MKLFLAGS += -Wl,--end-group
  endif

  ifeq ($(MKL_THREADS), YES)
	MKLFLAGS += -liomp5
  endif

	MKLFLAGS += -lpthread
	FFTWFLAGS = $(MKLFLAGS) -L$(FFTW3_LIBDIR)

# Flags for directly linking the original FFTW libraries.
else
  ifeq ($(TYPE_FFT), FFTW)
    ifeq ($(OMP), YES)
	FFTWFLAGS += -lfftw3_omp
    endif
	#FFTWFLAGS += -lfftw3 -lm -L/usr/local/FFTW/fftw-3.3.3/lib -L$(FFTW3_LIBDIR)
	#FFTWFLAGS += -lfftw3 -lm -L$(FFTW3_LIBDIR)
	FFTWFLAGS += -lfftw3 -lm 
	#FFTWFLAGS += -L/apps/fftw3/3.3.5-intel17.0/lib -lfftw3 -lm

  endif
endif


LINKERFLAGS = $(LDLIBS) $(FFTWFLAGS) $(OPT1) $(OMPADDL)

ifeq ($(LINK_STATIC), YES)
	LINKERFLAGS += $(STATIC)
endif


#####################################################################
#                          Executable name                          #
#####################################################################

# Set the date
DATE =$(shell date +%d.%m.%y)

# Set type of parallelization.
ifeq ($(MPI_PARALLEL), YES)
	MPI_NAME = MPI
else
  ifeq ($(OMP), YES)
	MPI_NAME = OMP
  else
	MPI_NAME = Mono
  endif
endif

EXEC = qdd

#####################################################################
#                 Phony targets and empty recipes                   #
#####################################################################

.PHONY: all initial_checks clean $(EXEC)

.DEFAULT_GOAL: all

all: $(EXEC)

clean:
	@rm -rfv *.o *.mod *.il ../../bin/$(EXEC)

# Prevent make from trying to update the source files:
makefile: ;
#define.mk: ;
#define.h: ;
*.F90: ;
#*.cu: ;

#####################################################################
#                             Checks                                #
#####################################################################

# List of valid compilers.
valid_CF90 := GFORT GFORT1 IFORT IFORT1 MPIF90 XLF_R

# List of valid MPI_PARALLEL:
valid_MPI_PARALLEL := YES NO

# List of valid OMP_THREADS:
valid_OMP_THREADS := DYN YES NO

# List of valid TYPE_FFT values.
valid_FFT := FFTW NETLIB MKL

#findiff_value  :=$(shell echo -e "\\n\#if(findiff)  \\nfindiff  \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
#numerov_value  :=$(shell echo -e "\\n\#if(numerov)  \\nnumerov  \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
#coufou_value   :=$(shell echo -e "\\n\#if(coufou)   \\ncoufou   \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
#coudoub_value  :=$(shell echo -e "\\n\#if(coudoub)  \\ncoudoub  \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
#coudoub3D_value:=$(shell echo -e "\\n\#if(coudoub3D)\\ncoudoub3D\\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
#raregas_value  :=$(shell echo -e "\\n\#if(raregas)  \\nraregas  \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)

initial_checks: makefile
	@echo "#######################################################"
	@echo  Lists of known compilation parameters.
	@echo "#######################################################"
	@echo  Parameters defined in define.h:
	@echo "--------------------------------------------"
	@echo "findiff       $(FINDIFF)"
	@echo "numerov       $(NUMEROV)"
	@echo "coufou        $(COUFOU)"
	@echo "coudoub       $(COUDOUB)"
	@echo "coudoub3D     $(COUDOUB3D)"
	@echo "raregas       $(RAREGAS)"
	@echo "#######################################################"
	@echo  Parameters defined in makefile / define.mk:
	@echo "--------------------------------------------"
	@echo "CF90          $(CF90)"
	@echo "MPI_PARALLEL  $(MPI_PARALLEL)"
	@echo "OMP_THREADS   $(OMP_THREADS)"
	@echo "TYPE_FFT      $(TYPE_FFT)"
	@echo "DEBUG         $(DEBUG)"
	@echo "LINK_STATIC   $(LINK_STATIC)"
	@echo "MKL_THREADS   $(MKL_THREADS)"
	@echo "MKL_WRAPPERS  $(MKL_WRAPPERS)"
	@echo "MKLPATH       $(MKLPATH)"
	@echo "MKLINCLUDE    $(MKLINCLUDE)"
	@echo "#######################################################"
	@echo "-D flags enabled (1) by makefile rules:"
	@echo "--------------------------------------------"
	@echo "parayes       $(parayes_value)"
	@echo "parano        $(parano_value)"
	@echo "netlib_fft    $(netlib_fft_value)"
	@echo "fftw_cpu      $(fftw_cpu_value)"
	@echo "fftwnomkl     $(fftwnomkl_value)"
	@echo "paropenmp     $(paropenmp_value)"
	@echo "dynopenmp     $(dynopenmp_value)"
	@echo "#######################################################"

# Check that the fortran compiler is defined.
ifneq ($(CF90), $(filter $(CF90), $(valid_CF90)))
	$(error ERROR: Unknown fortran compiler in makefile / define.mk (CF90=$(CF90)))
endif

# Check that the TYPE_FFT is valid:
ifneq ($(TYPE_FFT), $(filter $(TYPE_FFT), $(valid_FFT)))
	$(error ERROR: Unknown TYPE_FFT in makefile / define.mk ($(TYPE_FFT)))
endif

# Check that the parallelization and threading options are valid:
ifneq ($(MPI_PARALLEL), $(filter $(MPI_PARALLEL), $(valid_MPI_PARALLEL)))
	$(error ERROR: Unknown MPI_PARALLEL value in makefile / define.mk ($(MPI_PARALLEL)))
endif
ifneq ($(OMP_THREADS), $(filter $(OMP_THREADS), $(valid_OMP_THREADS)))
	$(error ERROR: Unknown OMP_THREADS value in makefile / define.mk ($(OMP_THREADS)))
endif

# Check that variables which should be either YES or NO have proper values:
ifneq ($(DEBUG),       $(filter $(DEBUG),       YES NO))
	$(error ERROR: Unknown $(strip DEBUG       value) in makefile / define.mk ($(DEBUG)))
endif
ifneq ($(LINK_STATIC), $(filter $(LINK_STATIC), YES NO))
	$(error ERROR: Unknown $(strip LINK_STATIC value) in makefile / define.mk ($(LINK_STATIC)))
endif
ifneq ($(MKL),         $(filter $(MKL)        , YES NO))
	$(error ERROR: Unknown $(strip MKL         value) in makefile / define.mk ($(MKL)))
endif

# coufou, coudoub: only one can be enabled.
#ifneq ($(COUFOU),NO)
#  ifneq ($(COUDOUB),NO)
#	$(error ERROR: coufou and coudoub cannot be enabled simultaneously in define.h)
#  endif
#endif

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

OBJINT = main.o params.o kinetic.o restart.o restartc.o init.o\
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
#   ifneq ($(TYPE_FFT), cuFFT)
OBJINT += coulsolv.o
MODLIST += coulsolv.mod
#   endif
# Object file for FFT:
ifeq ($(TYPE_FFT), NETLIB)
  FFT_OBJ = fftpack.o
endif
ifeq ($(calls_FFTW), YES)
  FFT_OBJ = fftw.o
  MODLIST += fftw.mod
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
#  ifneq ($(TYPE_FFT), cuFFT)
DEPS_KINCOUL = kinetic.mod coulsolv.mod
#  else
#  DEPS_KINCOUL = kinetic.mod
#  endif

# Dependencies on makefile, define.mk, define.h and params.mod
# Dependencies are inherited, so there's no need to repeat them
# in the targets below.
#$(MODLIST):   makefile define.mk define.h
#params.o:     makefile define.mk define.h
#$(OBJS_NPAR): makefile define.mk define.h params.mod
$(MODLIST):   makefile 
params.o:     makefile 
$(OBJS_NPAR): makefile params.mod

# Dependencies for objects that depend on kinetic.o and/or coulsolv.o,
# with no other special requirements:
$(OBJS_KIN):     kinetic.mod
#  ifneq ($(TYPE_FFT), cuFFT)
$(OBJS_COUL):    coulsolv.mod
#  endif
$(OBJS_KINCOUL): $(DEPS_KINCOUL)

# Set FFTW-dependent modules:
ifeq ($(calls_FFTW), YES)
coulsolv.mod: fftw.mod
kinetic.mod: fftw.mod
endif

#####################################################################
#                   Compilation and linkage rules                   #
#####################################################################

.SUFFIXES:
.SUFFIXES: .F90 .f .o .mod

$(EXEC): initial_checks $(OBJS)
	@echo Linking executable $@
	$(LINKER) -o $@ $(strip $(OBJS) $(LINKERFLAGS))
	mv -fv $(EXEC) ../../bin/

%.o: %.mod

# Implicit rule for objects:
%.o: %.F90
	$(COMPILER) $(COMPILERFLAGS1) -c $<

# Implicit rule for modules:
%.o %.mod: %.F90
	$(COMPILER) $(COMPILERFLAGS1) -c $<

# Additions for twostsic:
ifeq ($(FSIC),YES)

FSICDEPEND = 2st_util.o 2stUTc.o 2stUT.o
#OBJS_2ST_UTIL = main.o static.o
#OBJS_TWOSTSIC     = dynamic.F90 main.o
#OBJS_TWOSTSIC_R   = dynamic.o main.o restart.o static.o

#$(OBJS_2ST_UTIL): 2st_util.mod
#$(OBJS_TWOSTSIC):     twost.mod
#$(OBJS_TWOSTSIC_R):   twostr.mod

2st_util.o 2st_util.mod: fsic/2st_util.F90 kinetic.mod orthmat.o 
	$(COMPILER) $(COMPILERFLAGS1) -DREALSWITCH    -o 2st_util.o -c $<

2stUTc.o twost.mod: fsic/2stUT.F90 kinetic.mod twostr.mod 
	$(COMPILER) $(COMPILERFLAGS1) -DCOMPLEXSWITCH -o 2stUTc.o -c $<

2stUT.o twostr.mod: fsic/2stUT.F90 kinetic.mod 2st_util.mod
	$(COMPILER) $(COMPILERFLAGS1) -DREALSWITCH    -o 2stUT.o -c $<

endif

# Particular rules
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

ifeq ($(TYPE_FFT), NETLIB)
fftpack.o: fftpack.F90 fftpack2.F90
	$(COMPILER) $(COMPILERFLAGS2) -c $<
endif


main.o: main.F90 util.o orthmat.o $(FSICDEPEND)
	$(COMPILER) $(COMPILERFLAGS3) -c $<

givens.o: givens.F90
	$(COMPILER) $(COMPILERFLAGS3) -c $<

#multi_exp.o: multi_exp.F90 params.mod
#	$(COMPILER) $(COMPILERFLAGS3) -c $<

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

ifeq ($(RAREGAS), YES)

# QM/MM routines

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

