# for Linux
COMPILER       = ifort
#COMPILERFLAGS  =   -g -CB  -align all -w -autodouble -fpp -traceback
COMPILERFLAGS  =   -O3 -xW -p -align all -no-prec-div -w -fpp 
#COMPILERFLAGS  =   -O3 -xW -qp -traceback -no-prec-div -w -autodouble -fpp 
#COMPILERFLAGS  =   -O3 -CB -traceback -no-prec-div -w -autodouble -fpp 
#COMPILERFLAGS  =   -O1  -autodouble 
#COMPILERFLAGS  =   -O3 -xW  -w -autodouble -fpp -pg
#COMPILERFLAGS = -O2 -w -autodouble -xW -save
#COMPILERFLAGS  =   -O2 -w -autodouble -ipo -fpp
#COMPILERFLAGS2  =  -O3 -xW -traceback -ip  -w -autodouble  -fpp
#COMPILERFLAGS2  =  -O3 -xW -traceback -ip -align all -w -autodouble  -fpp
COMPILERFLAGS2  =  -O3 -xW -ip -p -w -fpp
#COMPILERFLAGS2  =  -g -CB  -w -autodouble  -fpp -traceback
#COMPILERFLAGS2  =  -O1  -w -autodouble  -fpp

#  setting for critical soubroutines which do not stand optimization
COMPILERFLAGS3  =  -g -w -p -align all -fpp 
#COMPILERFLAGS3  =  -g -w -autodouble -align all -fpp -traceback

#  setting for the FFT package, needs forced double precision
COMPILERFLAGS4 = $(COMPILERFLAGS) -autodouble


LINKER         = $(COMPILER)
LINKERFLAGS    = $(COMPILERFLAGS)




SOURCE_TRAILER = .F90
EXEC           = essai-f90
OBJSPERS	   = main.o params.o kinetic.o restart.o init.o\
                     generlcgo.o fftpack.o coulsolv.o\
                     static.o dynamic.o lda.o util.o abso_bc.o\
                     pseudosoft.o pseudogoed.o ionmd.o forces.o\
                     carlo.o localize.o subgrids.o analyse.o\
                     sic.o sicc.o rho.o rhoc.o nonloc.o nonlocc.o\
                     schmid.o zeroforce.o functions.o\
                     localize_rad.o symmcond_step.o loc_mfield.o givens.o\
                     parallele.o expevol.o 2stUT.o 2stUTc.o\
                     pot_substrate.o forces_substrate.o md_substrate.o\
                     short.o image.o lattice.o
OBJS		   = $(OBJSPERS) $(OBJSGEN) 
SHELL = /bin/sh

#---------------------------------------------------------	
all:$(EXEC) 
$(EXEC):$(OBJS) 	
	$(LINKER) $(LINKERFLAGS) -o $@  $(OBJS) $(STDLIB)
	mv $(EXEC) ../.
#
$(OBJSPERS) : define.h  params.F90 kinetic.F90 coulsolv.F90 Makefile
.SUFFIXES:  .o .F90

%.o: %.mod

.f90.o:                                                     
	$(COMPILER) $(COMPILERFLAGS) -c $<    

.F90.o:                                                     
	$(COMPILER) $(COMPILERFLAGS) $(PREFLAGS) -c $<    

main.o: main.F90 define.h  params.o kinetic.o coulsolv.o 2stUT.o


localize_rad.o: localize_rad.F90 define.h  params.o kinetic.o Makefile
	$(COMPILER) $(COMPILERFLAGS) -DREALSWITCH  -c -o localize_rad.o localize_rad.F90

sic.o: sicnew.F90 define.h  params.o  kinetic.o coulsolv.o symmcond_step.o Makefile
	$(COMPILER) $(COMPILERFLAGS) -DREALSWITCH  -c -o sic.o sicnew.F90
	$(COMPILER) $(COMPILERFLAGS) -DCOMPLEXSWITCH  -c -o sicc.o sicnew.F90

rho.o: rho.F90 define.h  params.o Makefile
	$(COMPILER) $(COMPILERFLAGS) -DREALSWITCH  -c rho.F90
	$(COMPILER) $(COMPILERFLAGS) -DCOMPLEXSWITCH  -c -o rhoc.o rho.F90

nonloc.o: nonloc.F90 define.h  params.o kinetic.o Makefile
	$(COMPILER) $(COMPILERFLAGS) -DREALSWITCH  -c nonloc.F90
	$(COMPILER) $(COMPILERFLAGS) -DCOMPLEXSWITCH  -c -o nonlocc.o nonloc.F90

2stUT.o: 2stUT.F90 define.h  params.o kinetic.o coulsolv.o symmcond_step.o Makefile
	$(COMPILER) $(COMPILERFLAGS) -DREALSWITCH  -c 2stUT.F90
	$(COMPILER) $(COMPILERFLAGS) -DCOMPLEXSWITCH  -c -o 2stUTc.o 2stUT.F90


givens.o: givens.F90 params.o Makefile define.h
	$(COMPILER) $(COMPILERFLAGS3) -c -o givens.o givens.F90

init.o: init.F90 params.o kinetic.o coulsolv.o Makefile define.h
	$(COMPILER) $(COMPILERFLAGS3) -c -o init.o init.F90

pseudosoft.o: pseudosoft.F90 params.o kinetic.o coulsolv.o Makefile define.h

loc_mfield.o: loc_mfield.F90 params.o kinetic.o coulsolv.o Makefile define.h

coulsolv.o: coulsolv.F90 falr.F90 coulex.F90 findiff-sinft.F90 params.o kinetic.o fftpack.o Makefile define.h

kinetic.o : kinetic.F90 fft.F90 findiff.F90 params.o  fftpack.o Makefile define.h

fftpack.o: fftpack.F90 fftpack2.F90 params.o Makefile define.h
	$(COMPILER) $(COMPILERFLAGS4) -c -o fftpack.o fftpack.F90
