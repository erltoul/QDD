# Makefile v2.4 (2010.05.05)
# CF90  = GFORT (serial version only)
#         IFORT (serial & parallel & debug version)
#         XLF_R (serial & parallel version)
#         MPIF90 (parallel version only)

# Insert the date
DATE    = 030112
# the compiler is set here 
CF90    = IFORT1
# the option for parallel processing (only for some compilers)
USE_MPI = NO
# Which FFT solver
TYPE_FFT = FFTW
# debugging option
DEBUG   = NO

#####################################################################
#                             Options                               #
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
	OPT1  =  -pg -w -O3 -mfpmath=sse -fdefault-real-8 -fdefault-double-8
	OPT2  =  -pg -w -O3 -mfpmath=sse -fdefault-real-8 -fdefault-double-8
	OPT3  =  -pg -w -g  -fdefault-real-8 -fdefault-double-8
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
	OPT1  =  -pg -w -O3 -msse4.2 -mfpmath=sse -fdefault-real-8 -fdefault-double-8 -ffast-math -finline-functions -funroll-loops
	OPT2  =  -pg -w -O3 -msse4.2 -mfpmath=sse -ffast-math -fdefault-real-8 -fdefault-double-8 -finline-functions -funroll-loops -fdefault-real-8 -fdefault-double-8
	OPT3  =  -pg -w -g  -fdefault-real-8 -fdefault-double-8
ifeq "$(DEBUG)" "YES"
	OPT1  =  -pg -cpp -w -g -fbacktrace                 
	OPT2  =  -pg -cpp -w -g -fbacktrace -fdefault-real-8
	OPT3  =  -pg -cpp -w -g                             
endif
endif


ifeq "$(CF90)" "IFORT"
#	OPT1  = -fpp -w -xW -O3 -ip -no-prec-div -align all -autodouble -static
#	OPT2  = -fpp -w -xW -O3 -ip -no-prec-div -align all -autodouble -static
#	OPT3  = -fpp -w  -pg      -g               -align all -autodouble -static
	OPT1  =  -pg -fpp -w -xW -O3 -ip -no-prec-div -align all -autodouble
	OPT2  =  -pg -fpp -w -xW -O3 -ip -no-prec-div -align all -autodouble
	OPT3  =  -pg -fpp -w        -g               -align all -autodouble
ifeq "$(DEBUG)" "YES"
	OPT1  =  -pg -fpp -w -g -CB -traceback -align all -autodouble
	OPT2  =  -pg -fpp -w -g -CB -traceback -align all -autodouble
	OPT3  =  -pg -fpp -w -g                -align all -autodouble
endif
endif


ifeq "$(CF90)" "IFORT1"
	OPT1  = -fpp -w -axaxv -msse4.2 -O3 -ip -no-prec-div -align all -autodouble -static
	OPT2  = -fpp -w -axaxv -msse4.2 -O3 -ip -no-prec-div -align all -autodouble -static
	OPT3  = -fpp -w        -g               -align all -autodouble -static
#	OPT1  =  -pg -fpp -w -axTP -msse2 -O3 -ip -no-prec-div -align all -autodouble -static
#	OPT2  =  -pg -fpp -w -axTP -msse2 -O3 -ip -no-prec-div -align all -autodouble -static
#	OPT3  =  -pg -fpp -w        -g               -align all -autodouble -static
ifeq "$(DEBUG)" "YES"
	OPT1  =  -pg -fpp -w -g -CB -traceback -align all -autodouble
	OPT2  =  -pg -fpp -w -g -CB -traceback -align all -autodouble
	OPT3  =  -pg -fpp -w -g                -align all -autodouble
endif
endif


ifeq "$(CF90)" "MPIF90"
        OPT1  =  -pg -fpp -w -O2 -no-prec-div -align all -autodouble -static
        OPT2  =  -pg -fpp -w -O2 -no-prec-div -align all -autodouble -static
        OPT3  =  -pg -fpp -w -g               -align all -autodouble -static
ifeq "$(DEBUG)" "YES"
        OPT1  =  -pg -fpp -w -g -CB -traceback -align all -autodouble
        OPT2  =  -pg -fpp -w -g -CB -traceback -align all -autodouble
        OPT3  =  -pg -fpp -w -g                -align all -autodouble
endif
endif


ifeq "$(CF90)" "XLF_R"
	OPT1  = -d -w -qport=mod -q64 -qarch=pwr6 -qtune=pwr6 -O3 -qstrict                    
	OPT2  = -d -w -qport=mod -q64 -qarch=pwr6 -qtune=pwr6 -O3 -qstrict -qdpc -qautodbl=dbl
	OPT3  = -d -w -qport=mod -q64 -qarch=pwr6 -qtune=pwr6                                 
endif


#####################################################################
#                             Compiler                              #
#####################################################################

# Insert the run command for invokung your compiler
# ifeq "$(CF90)" "OWN"
# 	COMPILER  = ...
# 	LDLIBS    = 
# endif


ifeq "$(CF90)" "MPIF90"
	COMPILER  = mpif90
	LDLIBS    = 
endif


ifeq "$(CF90)" "GFORT"
	COMPILER  = gfortran
	LDLIBS    =  
endif


ifeq "$(CF90)" "GFORT1"
	COMPILER  = gfortran
	LDLIBS    = 
endif

ifeq "$(CF90)" "IFORT"
	COMPILER  = ifort
	LDLIBS    = 
ifeq "$(USE_MPI)" "YES"
	LDLIBS    = -lmpi 
endif
endif

ifeq "$(CF90)" "IFORT1"
	COMPILER  = ifort
	LDLIBS    =  
ifeq "$(USE_MPI)" "YES"
	LDLIBS    = -lmpi 
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

COMPILERFLAGS1 =  $(LDLIBS)$(OPT1)
COMPILERFLAGS2 =  $(LDLIBS)$(OPT2) 
COMPILERFLAGS3 =  $(LDLIBS)$(OPT3) 

LINKER         = $(COMPILER)
LINKERFLAGS    = $(COMPILERFLAGS1) 


#####################################################################
#                              Rules                                #
#####################################################################

ifeq "$(USE_MPI)" "NO"
	MPI_NAME = Mono
endif
ifeq "$(USE_MPI)" "YES"
	MPI_NAME = Mpi
endif

IDRIS =
ifeq "$(CF90)" "XLF_R"
	IDRIS    = -WF,
endif

.PHONY: all clean

.SUFFIXES: .F90 .o

#EXEC = essai.seq
EXEC = tdks_v$(DATE)_$(MPI_NAME)_$(TYPE_FFT).bin

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

ifeq "$(TYPE_FFT)" "NETLIB"
OBJS = $(OBJINT) fftpack.o
endif

ifeq "$(TYPE_FFT)" "FFTW"
OBJS = $(OBJINT) fftw.o
endif

OBJH = define.h params.F90 kinetic.F90 coulsolv.F90 makefile

all:$(EXEC)
$(EXEC): $(OBJS)
	@echo Linking executable $@

ifeq "$(TYPE_FFT)" "FFTW"
	$(LINKER) $(LINKERFLAGS) -o $@ $(OBJS) -lfftw3 -lm
endif
ifeq "$(TYPE_FFT)" "NETLIB"
	$(LINKER) $(LINKERFLAGS) -o $@ $(OBJS)
endif

	mv $(EXEC) ../

$(OBJS):  define.h makefile params.F90

%.o: %.mod

.F90.o:
	$(COMPILER) $(COMPILERFLAGS1) -c $< 
main.o: main.F90 define.h makefile params.o kinetic.o coulsolv.o 2stUT.o 2stUTc.o localize_rad.o

ifeq "$(TYPE_FFT)" "NETLIB"
kinetic.o : kinetic.F90 fft.F90 findiff.F90 define.h makefile params.o fftpack.o
endif

ifeq "$(TYPE_FFT)" "FFTW"
#ifeq "$(DIM)" "1d"
#kinetic.o : kinetic.F90 fft.fftw1d.F90 findiff.F90 define.h makefile params.o fftw.o
#else
kinetic.o : kinetic.F90 fft.F90 findiff.F90 define.h makefile params.o fftw.o
#endif
endif

ifeq "$(TYPE_FFT)" "FFTW"
coulsolv.o: coulsolv.F90 falr.F90 coulex.fftw.F90 findiff-sinft.F90 define.h makefile params.o kinetic.o
endif

ifeq "$(TYPE_FFT)" "NETLIB"
coulsolv.o: coulsolv.F90 falr.F90 coulex.F90 coulex.fftw.F90 findiff-sinft.F90 define.h makefile params.o kinetic.o fftpack.o
endif

loc_mfield.o: loc_mfield.F90 define.h makefile params.o kinetic.o coulsolv.o

pseudosoft.o: pseudosoft.F90 define.h makefile params.o kinetic.o coulsolv.o

static.o: pseudosoft.F90 define.h makefile params.o kinetic.o localize_rad.o 2stUT.o

localize_rad.o: localize_rad.F90 define.h makefile params.o kinetic.o 
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DREALSWITCH    -o $@ -c $< 
rho.o: rho.F90 define.h makefile params.o
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DREALSWITCH    -o $@ -c $< 
rhoc.o: rho.F90 define.h makefile params.o
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DCOMPLEXSWITCH -o $@ -c $< 

localizer.o: localize.F90 define.h makefile params.o
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DREALSWITCH    -o $@ -c $< 

localize.o: localize.F90 define.h makefile params.o
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DCOMPLEXSWITCH -o $@ -c $< 

sicnew.o: sicnew.F90 define.h makefile params.o kinetic.o coulsolv.o symmcond_step.o
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DREALSWITCH    -o $@ -c $< 

sicnewc.o: sicnew.F90 define.h makefile params.o kinetic.o coulsolv.o symmcond_step.o
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DCOMPLEXSWITCH -o $@ -c $< 

nonloc.o: nonloc.F90 define.h makefile params.o kinetic.o
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DREALSWITCH    -o $@ -c $< 

nonlocc.o: nonloc.F90 define.h makefile params.o kinetic.o
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DCOMPLEXSWITCH -o $@ -c $< 

2stUT.o: 2stUT.F90 define.h makefile params.o kinetic.o coulsolv.o symmcond_step.o
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DREALSWITCH    -o $@ -c $< 
 
2stUTc.o: 2stUT.F90 define.h makefile params.o kinetic.o coulsolv.o symmcond_step.o
	$(COMPILER) $(COMPILERFLAGS1) $(IDRIS)-DCOMPLEXSWITCH -o $@ -c $< 


ifeq "$(TYPE_FFT)" "NETLIB"
fftpack.o: fftpack.F90 fftpack2.F90 define.h makefile params.o
	$(COMPILER) $(COMPILERFLAGS2) -c $< 
endif

givens.o: givens.F90 define.h makefile params.o
	$(COMPILER) $(COMPILERFLAGS3) -c $< 

init.o: init.F90 define.h makefile params.o kinetic.o coulsolv.o
	$(COMPILER) $(COMPILERFLAGS3) -c $< 

restart.o: restart.F90 define.h makefile params.o kinetic.o
	$(COMPILER) $(COMPILERFLAGS3) $(IDRIS)-DREALSWITCH    -o $@ -c $< 

restartc.o: restart.F90 define.h makefile params.o kinetic.o
	$(COMPILER) $(COMPILERFLAGS3) $(IDRIS)-DCOMPLEXSWITCH -o $@ -c $< 

parallele.o: parallele.F90 define.h makefile params.o kinetic.o
	$(COMPILER) $(COMPILERFLAGS3) -c $< 

clean:
	@rm -rf *.o *.mod *.il *.i *.f
