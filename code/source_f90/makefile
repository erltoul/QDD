# Makefile v2.4 (2010.05.05)
# CF90  = GFORT (serial version only)
#         IFORT (serial & parallel & debug version)
#         XLF_R (serial & parallel version)
#         MPIF90 (parallel version only)

# the compiler is set here
CF90    = IFORT
# the option for parallel processing (only for some compilers)
USE_MPI = NO
# Which FFT solver
TYPE_FFT = FFTW
# debugging option
DEBUG   = NO
# Link static or not.
LINK_STATIC = YES
# invoke openmp
OMP = NO
# Use MKL libraries
MKL = YES
# MKL Threading
MKL_PARALLEL = YES
# Path to MKL libraries
MKLPATH =/divers/intel/Compiler/12.1/mkl/lib/intel64
MKLPATH =
# For large boxe sizes, try option
# '-i-dynamic -mcmodel=medium' or '-shared-intel -mcmodel=medium'

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
	COMPILER  = /opt/mpich2-intel/bin/mpif90
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

#####################################################################
#                        Compilation rules                          #
#####################################################################

# Switches used to preprocess files.
ifeq "$(TYPE_FFT)" "FFTW"
	FFTWADD = -Dfftw_cpu
endif
ifeq "$(TYPE_FFT)" "NETLIB"
	FFTWADD = -Dnetlib_fft
endif

OMPADD =
ifeq "$(OMP)" "YES"
	OMPADD = $(OMPADDL) -Dparopenmp
endif

# Compiler flags
COMPILERFLAGS1 =  $(LDLIBS) $(OPT1) $(OMPADD) $(FFTWADD)
COMPILERFLAGS2 =  $(LDLIBS) $(OPT2) $(OMPADD) $(FFTWADD)
COMPILERFLAGS3 =  $(LDLIBS) $(OPT3) $(OMPADD) $(FFTWADD)

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
ifeq ($(MKL), YES)

	MKLFLAGS =

  ifeq ($(TYPE_FFT), FFTW) # MKL FFTW wrappers
    ifneq "$(MKL_WRAPPERS)" ""
	MKLFLAGS += -L$(MKL_WRAPPERS)
    endif
	MKLFLAGS += -lfftw3xf_intel
  endif

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

  ifeq ($(MKL_PARALLEL), YES)
	MKLFLAGS += -lmkl_intel_thread
  else
	MKLFLAGS += -lmkl_sequential
  endif

	MKLFLAGS += -lmkl_core

  ifeq ($(LINK_STATIC), YES) # Static MKL
	MKLFLAGS += -Wl,--end-group
  endif

  ifeq ($(MKL_PARALLEL), YES)
	MKLFLAGS += -liomp5
  endif

	MKLFLAGS += -lpthread

# Flags for directly linking the original FFTW libraries.
else
  ifeq ($(TYPE_FFT), FFTW)
    ifeq ($(OMP), YES)
	LINKERFLAGS += -lfftw3_omp
    endif
	LINKERFLAGS += -lfftw3 -lm
  endif
endif

LINKERFLAGS = $(LDLIBS) $(MKLFLAGS) $(OPT1) $(OMPADDL)

ifeq ($(LINK_STATIC), YES)
	LINKERFLAGS += $(STATIC)
endif


#####################################################################
#                          Executable name                          #
#####################################################################

# Set the date
DATE =$(shell date +%d.%m.%y)

# Set type of parallelization.
ifeq "$(USE_MPI)" "YES"
	MPI_NAME = Mpi
else
  ifeq "$(OMP)" "YES"
	MPI_NAME = OMP
  else
	MPI_NAME = Mono
  endif
endif

EXEC = tdks_v$(DATE)_$(MPI_NAME)_$(TYPE_FFT).bin

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
define.h: ;
*.F90: ;

#####################################################################
#                             Checks                                #
#####################################################################

# List of valid compilers.
valid_CF90 := GFORT GFORT1 IFORT IFORT1 MPIF90 XLF_R

# List of valid TYPE_FFT values.
valid_FFT := FFTW NETLIB

gridfft_value  :=$(shell echo -e "\\n\#if(gridfft)  \\ngridfft  \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
findiff_value  :=$(shell echo -e "\\n\#if(findiff)  \\nfindiff  \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
numerov_value  :=$(shell echo -e "\\n\#if(numerov)  \\nnumerov  \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
coufou_value   :=$(shell echo -e "\\n\#if(coufou)   \\ncoufou   \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
coudoub_value  :=$(shell echo -e "\\n\#if(coudoub)  \\ncoudoub  \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
coudoub3D_value:=$(shell echo -e "\\n\#if(coudoub3D)\\ncoudoub3D\\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
parayes_value  :=$(shell echo -e "\\n\#if(parayes)  \\nparayes  \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
parano_value   :=$(shell echo -e "\\n\#if(parano)   \\nparano   \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
simpara_value  :=$(shell echo -e "\\n\#if(simpara)  \\nsimpara  \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
#	dyn1:=`grep -i '^[[:space:]]*\#[[:space:]]*define[[:space:]]\+dynopenmp[[:space:]]\+[[:digit:]]' define.h | tail -n1 | grep -ic '^[[:space:]]*\#[[:space:]]*define[[:space:]]\+dynopenmp[[:space:]]\+1'`
dynopenmp_value:=$(shell echo -e "\\n\#if(dynopenmp)\\ndynopenmp\\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
fullsic_value  :=$(shell echo -e "\\n\#if(fullsic)  \\nfullsic  \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
symmcon_value  :=$(shell echo -e "\\n\#if(symmcon)  \\nsymmcon  \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
twostsic_value :=$(shell echo -e "\\n\#if(twostsic) \\ntwostsic \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
raregas_value  :=$(shell echo -e "\\n\#if(raregas)  \\nraregas  \\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)
fftwnomkl_value:=$(shell echo -e "\\n\#if(fftwnomkl)\\nfftwnomkl\\n\#else\\n0\\n\#endif" | cat define.h - | gcc -E -undef - 2>/dev/null | tail -n1)

initial_checks: makefile define.h
	@echo "#######################################################"
	@echo  Lists of known compilation parameters.
	@echo "#######################################################"
	@echo  Parameters defined in define.h:
	@echo "-------------------------------------"
	@echo "gridfft       $(gridfft_value)"
	@echo "findiff       $(findiff_value)"
	@echo "numerov       $(numerov_value)"
	@echo "coufou        $(coufou_value)"
	@echo "coudoub       $(coudoub_value)"
	@echo "coudoub3D     $(coudoub3D_value)"
	@echo "parayes       $(parayes_value)"
	@echo "parano        $(parano_value)"
	@echo "simpara       $(simpara_value)"
	@echo "dynopenmp     $(dynopenmp_value)"
	@echo "fullsic       $(fullsic_value)"
	@echo "symmcond      $(symmcon_value)"
	@echo "twostsic      $(twostsic_value)"
	@echo "raregas       $(raregas_value)"
	@echo "fftwnomkl     $(fftwnomkl_value)"
	@echo "#######################################################"
	@echo  Parameters defined in makefile:
	@echo "-------------------------------------"
	@echo "CF90          $(CF90)"
	@echo "USE_MPI       $(USE_MPI)"
	@echo "TYPE_FFT      $(TYPE_FFT)"
	@echo "DEBUG         $(DEBUG)"
	@echo "LINK_STATIC   $(LINK_STATIC)"
	@echo "OMP           $(OMP)"
	@echo "MKL           $(MKL)"
	@echo "#######################################################"

# Check that the fortran compiler is defined.
ifneq ($(CF90), $(filter $(CF90), $(valid_CF90)))
	$(error ERROR: Unknown fortran compiler in makefile (CF90=$(CF90)))
endif

# Check that the TYPE_FFT is valid:
ifneq ($(TYPE_FFT), $(filter $(TYPE_FFT), $(valid_FFT)))
	$(error ERROR: Unknown TYPE_FFT in makefile ($(TYPE_FFT)))
endif

# Check that variables which should be either YES or NO have proper values:
ifneq ($(USE_MPI),     $(filter $(USE_MPI),     YES NO))
	$(error ERROR: Unknown $(strip USE_MPI     value) in makefile ($(USE_MPI)))
endif
ifneq ($(DEBUG),       $(filter $(DEBUG),       YES NO))
	$(error ERROR: Unknown $(strip DEBUG       value) in makefile ($(DEBUG)))
endif
ifneq ($(LINK_STATIC), $(filter $(LINK_STATIC), YES NO))
	$(error ERROR: Unknown $(strip LINK_STATIC value) in makefile ($(LINK_STATIC)))
endif
ifneq ($(OMP),         $(filter $(OMP),         YES NO))
	$(error ERROR: Unknown $(strip OMP         value) in makefile ($(OMP)))
endif
ifneq ($(MKL),         $(filter $(MKL)        , YES NO))
	$(error ERROR: Unknown $(strip MKL         value) in makefile ($(MKL)))
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

# coudoub3D is only enabled when using FFTW.
ifneq ($(coudoub3D_value), 0)
  ifneq ($(TYPE_FFT), FFTW)
	$(error ERROR: coudoub3D requires the use of FFTW)
  endif
endif

# If parayes is enabled in define.h,
ifneq ($(parayes_value), 0)
# parano must be disabled,
  ifneq ($(parano_value), 0)
	$(error ERROR: parayes and parano cannot be enabled simultaneously in define.h)
  endif
# simpara must be disabled,
  ifneq ($(simpara_value), 0)
	$(error ERROR: parayes and simpara cannot be enabled simultaneously in define.h)
  endif
# and MPI compilation must be enabled;
  ifneq ($(USE_MPI), YES)
	$(error ERROR: parayes requires USE_MPI=YES in makefile)
  endif
# if parayes is disabled,
else
# then parano must be enabled.
  ifeq ($(parano_value), 0)
	$(error ERROR: either parayes or parano should be enabled in define.h)
  endif
endif

# If simpara is enabled in define.h,
ifneq ($(simpara_value), 0)
# parano must also be enabled
  ifeq ($(parano_value), 0)
	$(error ERROR: in order to enable simpara, you also need to enable parano in define.h)
  endif
# and MPI compilation must be enabled.
  ifneq ($(USE_MPI), YES)
	$(error ERROR: simpara requires USE_MPI=YES in makefile)
  endif
endif

# If parano is enabled in define.h, and simpara is disabled, don't use MPI:
ifneq ($(parano_value), 0)
  ifeq ($(simpara_value), 0)
    ifneq ($(USE_MPI), NO)
	$(error ERROR: if parano is enabled and simpara is disabled, USE_MPI should be NO in the makefile)
    endif
  endif
endif

# Check collision of MPI and OpenMP:
ifeq ($(USE_MPI), YES)
  ifeq ($(OMP), YES)
	$(error ERROR: code not ready for simultaneous MPI and OpenMP parallelizations)
  endif
endif

# dynopenmp requires OpenMP:
ifneq ($(dynopenmp_value), 0)
  ifneq ($(OMP), YES)
	$(error ERROR: You are trying to compile the code with dynopenmp enabled in define.h, but OMP is disabled in the makefile)
  endif
endif

# Checks connecting FFTW, MKL, and fftwnomkl:
ifeq ($(FFT_TYPE), FFTW)
  ifeq ($(MKL), YES)
    ifneq ($(fftwnomkl_value), 0)
	$(error ERROR: If you compile with FFTW and MKL, you must disable fftwnomkl in define.h)
    endif
  else
    ifeq ($(fftwnomkl_value), 0)
	$(error ERROR: If you compile with FFTW and MKL, you must enable fftwnomkl in define.h)
    endif
  endif
endif


	@echo "Done with the initial checks, starting compilation..."
	@echo "#######################################################"

#####################################################################
#                       Targets and dependencies                    #
#####################################################################

OBJINT = main.o params.o kinetic.o restart.o restartc.o init.o\
       generlcgo.o coulsolv.o\
       static.o dynamic.o lda.o util.o abso_bc.o\
       pseudosoft.o pseudogoed.o ionmd.o forces.o\
       carlo.o localize.o localizer.o subgrids.o analyse.o\
       sicnew.o sicnewc.o rho.o rhoc.o nonloc.o nonlocc.o\
       schmid.o zeroforce.o functions.o\
       localize_rad.o symmcond_step.o loc_mfield.o givens.o\
       parallele.o expevol.o 2stUT.o 2stUTc.o\
       pot_substrate.o forces_substrate.o md_substrate.o\
       short.o image.o lattice.o

# Object file for FFT:
ifeq "$(TYPE_FFT)" "NETLIB"
FFT_OBJ = fftpack.o
endif
ifeq "$(TYPE_FFT)" "FFTW"
FFT_OBJ = fftw.o
endif

OBJS      = $(OBJINT) $(FFT_OBJ)

# List of objects excluding params.o:
OBJS_NPAR = $(filter-out params.o, $(OBJS))

# Objects that depend on kinetic.o and coulsolv.o, with no other special requirements:
OBJS_COUL = loc_mfield.o pseudosoft.o init.o

DEPS_COUL = kinetic.o coulsolv.o

OBJH = makefile define.h params.F90 kinetic.F90 coulsolv.F90

#####################################################################
#                   Compilation and linkage rules                   #
#####################################################################

.SUFFIXES:
.SUFFIXES: .F90 .o

$(EXEC): initial_checks $(OBJS)
	@echo Linking executable $@
	$(LINKER) -o $@ $(strip $(OBJS) $(LINKERFLAGS))
	mv -fv $(EXEC) ../


# Dependencies on makefile, define.h and params.o
# Dependencies are inherited, so there's no need to repeat them
# in the targets below.
params.o: makefile define.h

$(OBJS_NPAR): makefile define.h params.o

# Rule for objects that depend on kinetic.o and coulsolv.o, with no other special requirements:
$(OBJS_COUL): $(DEPS_COUL)

%.o: %.mod

# Implicit rule
.F90.o:
	$(COMPILER) $(COMPILERFLAGS1) -c $<

# Particular rules
# (complex dependencies, REAL/COMPLEX switch, non-default COMPILERFLAGS).

main.o: main.F90 $(DEPS_COUL) 2stUT.o 2stUTc.o localize_rad.o

kinetic.o : kinetic.F90 fft.F90 findiff.F90 $(FFT_OBJ)

coulsolv.o: coulsolv.F90 falr.F90 coulex.F90 findiff-sinft.F90 kinetic.o $(FFT_OBJ)

static.o: pseudosoft.F90 kinetic.o localize_rad.o 2stUT.o

localize_rad.o: localize_rad.F90 kinetic.o
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DREALSWITCH    -o $@ -c $<

rho.o: rho.F90
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DREALSWITCH    -o $@ -c $<

rhoc.o: rho.F90
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DCOMPLEXSWITCH -o $@ -c $<

localizer.o: localize.F90
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DREALSWITCH    -o $@ -c $<

localize.o: localize.F90
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DCOMPLEXSWITCH -o $@ -c $<

sicnew.o: sicnew.F90 $(DEPS_COUL) symmcond_step.o
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DREALSWITCH    -o $@ -c $<

sicnewc.o: sicnew.F90 $(DEPS_COUL) symmcond_step.o
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DCOMPLEXSWITCH -o $@ -c $<

nonloc.o: nonloc.F90 kinetic.o
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DREALSWITCH    -o $@ -c $<

nonlocc.o: nonloc.F90 kinetic.o
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DCOMPLEXSWITCH -o $@ -c $<

2stUT.o: 2stUT.F90 $(DEPS_COUL) symmcond_step.o
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DREALSWITCH    -o $@ -c $<

2stUTc.o: 2stUT.F90 $(DEPS_COUL) symmcond_step.o
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DCOMPLEXSWITCH -o $@ -c $<

ifeq "$(TYPE_FFT)" "NETLIB"
fftpack.o: fftpack.F90 fftpack2.F90
	$(COMPILER) $(COMPILERFLAGS2) -c $<
endif

givens.o: givens.F90
	$(COMPILER) $(COMPILERFLAGS3) -c $<

restart.o: restart.F90 kinetic.o
	$(COMPILER) $(COMPILERFLAGS3) $(IDRIS)-DREALSWITCH    -o $@ -c $<

restartc.o: restart.F90 kinetic.o
	$(COMPILER) $(COMPILERFLAGS3) $(IDRIS)-DCOMPLEXSWITCH -o $@ -c $<

parallele.o: parallele.F90 kinetic.o
	$(COMPILER) $(COMPILERFLAGS3) -c $<

