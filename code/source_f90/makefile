# Makefile

#####################################################################
#                             Defaults                              #
#####################################################################
# PLEASE, DO NOT MODIFY THESE PARAMETERS HERE, EDIT THE INCLUDE'D   #
# FILE (define.mk) IF YOU NEED TO CHANGE ANY OF THEM                #
#####################################################################

# CF90: Fortran compiler
# Available options:
#   * GFORT (serial version only)
#   * IFORT (serial & parallel & debug version)
#   * XLF_R (serial & parallel version)
#   * MPIF90 (parallel version only)
CF90 = IFORT

# MPI_PARALLEL: Variable for parallel processing (only for some compilers)
# Available options:
#   * YES
#   * SIM: simulate parallelization (simultaneous mono jobs)
#   * NO
MPI_PARALLEL = NO

# OMP_THREADS: Invoke OpenMP
# Available options:
#   * DYN: wave function parallelization (threads)
#   * YES: use threads for FFT
#   * NO
OMP_THREADS = NO

# TYPE_FFT: FFT solver
# Available options:
#   * NETLIB
#   * FFTW
#   * MKL
#   * cuFFT
TYPE_FFT = cuFFT

# DEBUG: enable debugging
# Available options:
#   * YES
#   * NO
DEBUG  = NO

# LINK_STATIC: select static linkage of the binary
# Available options:
#   * YES
#   * NO
LINK_STATIC = YES

# MKL_THREADS: Enable MKL threading (only used for TYPE_FFT = FFTW)
# Available options:
#   * YES
#   * NO
MKL_THREADS = YES

# MKL_WRAPPERS: Path to the MKL FFTW wrappers (only used for TYPE_FFT = FFTW)
# Default value:
MKL_WRAPPERS =

# MKLPATH: Path to the MKL libraries (only used for TYPE_FFT = FFTW)
# Default value:
MKLPATH =
# Path for hyperion:
# MKLPATH = /divers/intel/Compiler/12.1/mkl/lib/intel64

# MKLINCLUDE: Include path for the MKL (only used for TYPE_FFT = FFTW )
# Default value:
MKLINCLUDE =

# For large boxe sizes, try option
# '-i-dynamic -mcmodel=medium' or '-shared-intel -mcmodel=medium'

#####################################################################
#                 Include customized makefile                       #
#####################################################################

# Include file with customized compilation parameters
include define.mk

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
	OMP= NO
else
	OMP= YES
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
	OMPADDL =  -openmp
  endif
endif


ifeq "$(CF90)" "IFORT1"
	COMPILER  = ifort
	LDLIBS    =
  ifeq "$(USE_MPI)" "YES"
	LDLIBS    = -lmpi
  endif
  ifeq "$(OMP)" "YES"
	OMPADDL =  -openmp
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

ifeq "$(TYPE_FFT)" "cuFFT"
	GPUCOMPILER  = nvcc
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
	OPT1  =  -w -O3 -mfpmath=sse -fdefault-real-8 -fdefault-double-8
	OPT2  =  -w -O3 -mfpmath=sse -fdefault-real-8 -fdefault-double-8
	OPT3  =  -w -g  -fdefault-real-8 -fdefault-double-8
#	OPT1  = -cpp -w -O3 -march=k8-sse3 -mfpmath=sse
#	OPT2  = -cpp -w -O3 -march=k8-sse3 -mfpmath=sse -fdefault-real-8
#	OPT3  = -cpp -w -g
  ifeq "$(DEBUG)" "YES"
	OPT1  =  -pg -w -g -fbacktrace -fdefault-real-8 -fdefault-double-8
	OPT2  =  -pg -w -g -fbacktrace -fdefault-real-8 -fdefault-double-8
	OPT3  =  -pg -w -g -fbacktrace -fdefault-real-8 -fdefault-double-8
  endif
endif


ifeq "$(CF90)" "GFORT1"
	OPT1  =  -w -O3 -msse4.2 -mfpmath=sse -ffast-math -fdefault-real-8 -fdefault-double-8 -finline-functions -funroll-loops
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
	OPT1  =  -fpp -w -xW -O3 -ip -no-prec-div -align all -autodouble
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

ifeq "$(TYPE_FFT)" "cuFFT"
	#GPUCOMPILERFLAGS = -O3 --gpu-architecture sm_13 
	#GPUCOMPILERFLAGS = -O3 --gpu-architecture sm_20 
	#GPUCOMPILERFLAGS = -O3 --gpu-architecture sm_30 
	OPTGPU = -O3 --gpu-architecture sm_35 
  ifeq "$(DEBUG)" "YES"
       	#GPUCOMPILERFLAGS = -g -G --gpu-architecture sm_13 
	#GPUCOMPILERFLAGS = -g -G --gpu-architecture sm_20 
	#GPUCOMPILERFLAGS = -g -G --gpu-architecture sm_30 
	OPTGPU = -g -G --gpu-architecture sm_35 
  endif
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
  ifeq ($(TYPE_FFT), cuFFT)
	FFTWADD = -Dfftw_gpu
	fftw_gpu_value := 1
  endif
endif

ifeq ($(MPI_PARALLEL), YES)
	MPIADD := -Dparayes
	parayes_value := 1
else
	MPIADD = -Dparano
	parano_value := 1
  ifeq ($(MPI_PARALLEL), SIM)
	MPIADD += -Dsimpara
	simpara_value := 1
  endif
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
COMPILERFLAGS1 = $(LDLIBS) $(OPT1) $(strip $(MPIADD) $(OMPADD) $(FFTWADD))
COMPILERFLAGS2 = $(LDLIBS) $(OPT2) $(strip $(MPIADD) $(OMPADD) $(FFTWADD))
COMPILERFLAGS3 = $(LDLIBS) $(OPT3) $(strip $(MPIADD) $(OMPADD) $(FFTWADD))
GPUCOMPILERFLAGS = $(OPTGPU) $(MPIADD) $(OMPADD)
# This will also be used in the compilation (with -DREALSWITCH / -DCOMPLEXSWITCH)
IDRIS =
ifeq "$(CF90)" "XLF_R"
	IDRIS    = -WF,
endif

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
	FFTWFLAGS = $(MKLFLAGS)

# Flags for directly linking the original FFTW libraries.
else
  ifeq ($(TYPE_FFT), FFTW)
    ifeq ($(OMP), YES)
	FFTWFLAGS += -lfftw3_omp
    endif
	FFTWFLAGS += -lfftw3 -lm
  endif
endif

#Flags for directly linking the cuda libraries
ifeq ($(TYPE_FFT), cuFFT)
        CUDAFLAGS += -lcufft -lcudart
endif

LINKERFLAGS = $(LDLIBS) $(FFTWFLAGS) $(OPT1) $(OMPADDL) $(CUDAFLAGS)

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
        EXEC = essai.par
else
  ifeq ($(MPI_PARALLEL), SIM)
	MPI_NAME = SIM
        EXEC = essai.sim
  else
    ifeq ($(OMP), YES)
	MPI_NAME = OMP
    else
	MPI_NAME = Mono
    endif
    EXEC = essai.seq
  endif
endif

#####################################################################
#                 Phony targets and empty recipes                   #
#####################################################################

.PHONY: all initial_checks clean

.DEFAULT_GOAL: all

all: $(EXEC)

clean:
	@rm -rfv *.o *.mod *.il $(shell ls *.i | grep -v makefile.i) *.f

# Prevent make from trying to update the source files:
makefile: ;
define.mk: ;
define.h: ;
*.F90: ;
*.cu: ;
define_cuda.h: ;
params_gpu.h: ;

#####################################################################
#                             Checks                                #
#####################################################################

# List of valid compilers.
valid_CF90 := GFORT GFORT1 IFORT IFORT1 MPIF90 XLF_R

# List of valid MPI_PARALLEL:
valid_MPI_PARALLEL := YES SIM NO

# List of valid OMP_THREADS:
valid_OMP_THREADS := DYN YES NO

# List of valid TYPE_FFT values.
valid_FFT := FFTW NETLIB MKL cuFFT

gridfft_value  :=$(shell echo -e "\\n\#if(gridfft)  \\ngridfft  \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
findiff_value  :=$(shell echo -e "\\n\#if(findiff)  \\nfindiff  \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
numerov_value  :=$(shell echo -e "\\n\#if(numerov)  \\nnumerov  \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
coufou_value   :=$(shell echo -e "\\n\#if(coufou)   \\ncoufou   \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
coudoub_value  :=$(shell echo -e "\\n\#if(coudoub)  \\ncoudoub  \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
coudoub3D_value:=$(shell echo -e "\\n\#if(coudoub3D)\\ncoudoub3D\\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
twostsic_value :=$(shell echo -e "\\n\#if(twostsic) \\ntwostsic \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
raregas_value  :=$(shell echo -e "\\n\#if(raregas)  \\nraregas  \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
lda_gpu_value  :=$(shell echo -e "\\n\#if(lda_gpu)  \\nlda_gpu  \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)

initial_checks: makefile define.mk define.h
	@echo "#######################################################"
	@echo  Lists of known compilation parameters.
	@echo "#######################################################"
	@echo  Parameters defined in define.h:
	@echo "--------------------------------------------"
	@echo "gridfft       $(gridfft_value)"
	@echo "findiff       $(findiff_value)"
	@echo "numerov       $(numerov_value)"
	@echo "coufou        $(coufou_value)"
	@echo "coudoub       $(coudoub_value)"
	@echo "coudoub3D     $(coudoub3D_value)"
	@echo "twostsic      $(twostsic_value)"
	@echo "raregas       $(raregas_value)"
	@echo "lda_gpu       $(lda_gpu_value)"
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
	@echo "simpara       $(simpara_value)"
	@echo "netlib_fft    $(netlib_fft_value)"
	@echo "fftw_cpu      $(fftw_cpu_value)"
	@echo "fftw_gpu      $(fftw_gpu_value)"
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

# gridfft, findiff, numerov: only one can be enabled.
ifneq ($(gridfft_value), 0)
  ifneq ($(findiff_value), 0)
	$(error ERROR: gridfft and findiff cannot be enabled simultaneously in define.h)
  endif
  ifneq ($(numerov_value), 0)
	$(error ERROR: gridfft and numerov cannot be enabled simultaneously in define.h)
  endif
else
	$(warning WARNING: gridfft disabled, so TYPE_FFT=$(TYPE_FFT) will be disregarded)
  ifeq ($(findiff_value), 0)
    ifeq ($(numerov_value), 0)
	$(error ERROR: gridfft, findiff and numerov cannot be simultaneously disabled in define.h)
    endif
  endif
#  ifneq ($(findiff_value), 0)
#    ifneq ($(numerov_value), 0)
#	$(error ERROR: findiff and numerov cannot be enabled simultaneously in define.h)
#    endif
#  endif
endif

# coufou, coudoub: only one can be enabled.
ifneq ($(coufou_value), 0)
  ifneq ($(coudoub_value), 0)
	$(error ERROR: coufou and coudoub cannot be enabled simultaneously in define.h)
  endif
endif

# coudoub3D is only enabled when using FFTW calls.
ifneq ($(coudoub3D_value), 0)
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

# lda_gpu is only enabled when using cuFFT calls.
ifneq ($(lda_gpu_value),0)
  ifneq ($(TYPE_FFT),cuFFT)
       $(error ERROR: lda_gpu requires the use of GPU)
  endif
endif
	@echo "Done with the initial checks, starting compilation..."
	@echo "#######################################################"

#####################################################################
#                  Targets and common dependencies                  #
#####################################################################

OBJINT = main.o params.o kinetic.o restart.o restartc.o init.o\
       generlcgo.o\
       static.o dynamic.o lda.o util.o abso_bc.o\
       pseudosoft.o pseudogoed.o ionmd.o forces.o\
       carlo.o localize.o localizer.o subgrids.o analyse.o\
       sicnew.o sicnewc.o rho.o rhoc.o nonloc.o nonlocc.o\
       schmid.o zeroforce.o functions.o\
       loc_mfield.o givens.o\
       parallele.o expevol.o 2stUT.o \
       pot_substrate.o forces_substrate.o md_substrate.o\
       short.o image.o lattice.o attachement.o

# List of modules:
MODLIST = params.mod kinetic.mod
ifneq ($(TYPE_FFT), cuFFT)
OBJINT += coulsolv.o
MODLIST += coulsolv.mod
endif
# Object file for FFT:
ifeq ($(TYPE_FFT), NETLIB)
  FFT_OBJ = fftpack.o
endif
ifeq ($(calls_FFTW), YES)
  FFT_OBJ = fftw.o
  MODLIST += fftw.mod
endif
ifeq ($(TYPE_FFT), cuFFT)
  FFT_OBJ += cuda.o cuda_alloc.o
  MODLIST += cuda_alloc.mod
endif
# Sic-conditional objects:
ifneq ($(twostsic_value), 0)
  OBJINT  += localize_rad.o 2stUTc.o
  MODLIST += localize_rad.mod twostr.mod twost.mod
endif

OBJS = $(OBJINT) $(FFT_OBJ)

# List of objects excluding params.o:
OBJS_NPAR = $(filter-out params.o, $(OBJS))

# Objects that depend on kinetic.mod and/or coulsolv.mod,
# with no other special requirements:
OBJS_KIN  = dynamic.o expevol.o ionmd.o util.o zeroforce.o
OBJS_COUL = image.o lda.o loc_mfield.o
OBJS_KINCOUL = init.o main.o pseudosoft.o
ifneq ($(TYPE_FFT), cuFFT)
DEPS_KINCOUL = kinetic.mod coulsolv.mod
else
DEPS_KINCOUL = kinetic.mod
endif
ifeq ($(TYPE_FFT), cuFFT)
  DEPS_KINCOUL += cuda_alloc.mod
endif

# Dependencies on makefile, define.mk, define.h and params.mod
# Dependencies are inherited, so there's no need to repeat them
# in the targets below.
$(MODLIST):   makefile define.mk define.h
params.o:     makefile define.mk define.h
$(OBJS_NPAR): makefile define.mk define.h params.mod

# Dependencies for objects that depend on kinetic.o and/or coulsolv.o,
# with no other special requirements:
$(OBJS_KIN):     kinetic.mod
ifneq ($(TYPE_FFT), cuFFT)
$(OBJS_COUL):    coulsolv.mod
endif
$(OBJS_KINCOUL): $(DEPS_KINCOUL)

# Set FFTW-dependent modules:
ifeq ($(calls_FFTW), YES)
coulsolv.mod: fftw.mod
kinetic.mod: fftw.mod
endif
ifeq "$(TYPE_FFT)" "cuFFT"
kinetic.mod: cuda_alloc.mod
endif

#####################################################################
#                   Compilation and linkage rules                   #
#####################################################################

.SUFFIXES:
.SUFFIXES: .F90 .o .mod

$(EXEC): initial_checks $(OBJS)
	@echo Linking executable $@
	$(LINKER) -o $@ $(strip $(OBJS) $(LINKERFLAGS))
	mv -fv $(EXEC) ../

%.o: %.mod

# Implicit rule for objects:
%.o: %.F90
	$(COMPILER) $(COMPILERFLAGS1) -c $<

# Implicit rule for modules:
%.o %.mod: %.F90
	$(COMPILER) $(COMPILERFLAGS1) -c $<

# Additions for twostsic:
ifneq ($(twostsic_value), 0)

OBJS_LOCALIZE_RAD = main.o static.o
OBJS_TWOSTSIC     = dynamic.F90 expevol.F90 main.o
OBJS_TWOSTSIC_R   = dynamic.o main.o restart.o static.o

$(OBJS_LOCALIZE_RAD): localize_rad.mod
$(OBJS_TWOSTSIC):     twost.mod
$(OBJS_TWOSTSIC_R):   twostr.mod

localize_rad.o localize_rad.mod: localize_rad.F90 kinetic.mod
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DREALSWITCH    -o localize_rad.o -c $<

2stUTc.o twost.mod: 2stUT.F90 kinetic.mod twostr.mod
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DCOMPLEXSWITCH -o 2stUTc.o -c $<

2stUT.o twostr.mod: 2stUT.F90 kinetic.mod localize_rad.mod
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DREALSWITCH    -o 2stUT.o -c $<
else
2stUT.o: 2stUT.F90 kinetic.mod
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DREALSWITCH    -o $@ -c $<
endif

# Particular rules
# (complex dependencies, REAL/COMPLEX switch, non-default COMPILERFLAGS).

kinetic.o : kinetic.F90 fft.F90 findiff.F90 $(FFT_OBJ)

ifneq "$(TYPE_FFT)" "cuFFT"
coulsolv.o: coulsolv.F90 falr.F90 coulex.F90 findiff-sinft.F90 kinetic.mod $(FFT_OBJ)
endif
static.o: pseudosoft.F90 $(DEPS_KINCOUL)

rho.o: rho.F90
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DREALSWITCH    -o $@ -c $<

rhoc.o: rho.F90 kinetic.mod
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DCOMPLEXSWITCH -o $@ -c $<

localizer.o: localize.F90 kinetic.mod
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DREALSWITCH    -o $@ -c $<

localize.o: localize.F90 kinetic.mod
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DCOMPLEXSWITCH -o $@ -c $<

sicnew.o: sicnew.F90 $(DEPS_KINCOUL)
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DREALSWITCH    -o $@ -c $<

sicnewc.o: sicnew.F90 $(DEPS_KINCOUL)
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DCOMPLEXSWITCH -o $@ -c $<

nonloc.o: nonloc.F90 kinetic.mod
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DREALSWITCH    -o $@ -c $<

nonlocc.o: nonloc.F90 kinetic.mod
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DCOMPLEXSWITCH -o $@ -c $<

ifeq ($(TYPE_FFT), NETLIB)
fftpack.o: fftpack.F90 fftpack2.F90
	$(COMPILER) $(COMPILERFLAGS2) -c $<
endif

givens.o: givens.F90
	$(COMPILER) $(COMPILERFLAGS3) -c $<

restart.o: restart.F90 kinetic.mod
	$(COMPILER) $(COMPILERFLAGS3) $(IDRIS)-DREALSWITCH    -o $@ -c $<

restartc.o: restart.F90 kinetic.mod
	$(COMPILER) $(COMPILERFLAGS3) $(IDRIS)-DCOMPLEXSWITCH -o $@ -c $<

parallele.o: parallele.F90 kinetic.mod
	$(COMPILER) $(COMPILERFLAGS3) -c $<

ifeq "$(TYPE_FFT)" "cuFFT"
cuda.o: cuda.cu define_cuda.h params_gpu.h static_gpu.cu makefile
	$(GPUCOMPILER) $(GPUCOMPILERFLAGS) -c $<
endif

