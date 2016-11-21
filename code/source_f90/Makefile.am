# Compiler-dependent options
if GFORT
if USE_DEBUG
OPT1  =  -pg -w -g -fbacktrace -fdefault-real-8 -fdefault-double-8
OPT2  =  $(OPT1)
OPT3  =  $(OPT2)
else
OPT1=-w -O3 -mfpmath=sse -fdefault-real-8 -fdefault-double-8
OPT2=-w -O3 -mfpmath=sse -fdefault-real-8 -fdefault-double-8
OPT3=-w -g  -fdefault-real-8 -fdefault-double-8
endif
endif

if IFORT
if USE_DEBUG
OPT1  =  -pg -fpp -w -g -CB -traceback -align all -autodouble
OPT2  =  -pg -fpp -w -g -CB -traceback -align all -autodouble
OPT3  =  -pg -fpp -w -g                -align all -autodouble
else
OPT1  = -fpp -w -xW -O3 -ip -no-prec-div -align all -autodouble
OPT2  = -fpp -w -xW -O3 -ip -no-prec-div -align all -autodouble
OPT3  = -fpp -w       -g                 -align all -autodouble
#~ OPT1  =-fpp -w -axsse4.2 -msse4.2 -O3 -ip -no-prec-div -align all -autodouble
#~ OPT2  =-fpp -w -axsse4.2 -msse4.2 -O3 -ip -no-prec-div -align all -autodouble
#~ OPT3  =-fpp -w        -g                            -align all -autodouble
endif
endif


if XLF_R
OPT1  =-d -w -qport=mod -q64 -qarch=pwr6 -qtune=pwr6 -O3 -qstrict
OPT2  =-d -w -qport=mod -q64 -qarch=pwr6 -qtune=pwr6 -O3 -qstrict -qdpc -qautodbl=dbl
OPT3  =-d -w -qport=mod -q64 -qarch=pwr6 -qtune=pwr6
endif

# underlying compiler for mpif90
if MPIF90
if MPI_GFORT
if USE_DEBUG
OPT1  =  -pg -w -g -fbacktrace -fdefault-real-8 -fdefault-double-8
OPT2  =  $(OPT1)
OPT3  =  $(OPT2)
else
OPT1=-w -O3 -mfpmath=sse -fdefault-real-8 -fdefault-double-8
OPT2=-w -O3 -mfpmath=sse -fdefault-real-8 -fdefault-double-8
OPT3=-w -g  -fdefault-real-8 -fdefault-double-8
endif
endif

if MPI_IFORT
if USE_DEBUG
OPT1  =  -pg -fpp -w -g -CB -traceback -align all -autodouble
OPT2  =  -pg -fpp -w -g -CB -traceback -align all -autodouble
OPT3  =  -pg -fpp -w -g                -align all -autodouble
else
OPT1  = -fpp -w -xW -O3 -ip -no-prec-div -align all -autodouble
OPT2  = -fpp -w -xW -O3 -ip -no-prec-div -align all -autodouble
OPT3  = -fpp -w       -g                 -align all -autodouble
endif
endif

if MPI_XLF_R
if USE_DEBUG
OPT1  =-d -w -qport=mod -q64 -qarch=pwr6 -qtune=pwr6 -O3 -qstrict
OPT2  =-d -w -qport=mod -q64 -qarch=pwr6 -qtune=pwr6 -O3 -qstrict -qdpc -qautodbl=dbl
OPT3  =-d -w -qport=mod -q64 -qarch=pwr6 -qtune=pwr6
else
OPT1  =-d -w -qport=mod -q64 -qarch=pwr6 -qtune=pwr6 -O3 -qstrict
OPT2  =-d -w -qport=mod -q64 -qarch=pwr6 -qtune=pwr6 -O3 -qstrict -qdpc -qautodbl=dbl
OPT3  =-d -w -qport=mod -q64 -qarch=pwr6 -qtune=pwr6
endif
endif
endif

if USE_CUDA
if USE_DEBUG
#OPTGPU = -g -G --gpu-architecture sm_13 
OPTGPU = -g -G --gpu-architecture sm_20 
#OPTGPU = -g -G --gpu-architecture sm_30 
#OPTGPU = -g -G --gpu-architecture sm_35 
else
#OPTGPU = -O3 --gpu-architecture sm_13 
#OPTGPU = -O3 --gpu-architecture sm_20 
#OPTGPU = -O3 --gpu-architecture sm_30 
OPTGPU = -O3 --gpu-architecture sm_35  -maxrregcount=72 --use_fast_math --ptxas-options=-v 
endif
endif
#for CCC with older GPU cards 
# OPTGPU = -O3 --gpu-architecture sm_20  -maxrregcount=20 --use_fast_math --ptxas-options=-v
#endif


if !XLF_R
if USE_STATIC
STATIC=-static
endif
endif
#FFT type dependent options
#~ LDLIBS=-L/usr/local/FFTW/fftw-3.3.3/lib

#Compiler options
COMPILERFLAGS1= $(STATIC) $(LD_MPI) $(OPENMP_FCFLAGS) $(OPT1) $(DMPI) $(DOMP) $(DFFT) 
COMPILERFLAGS2= $(STATIC) $(LD_MPI) $(OPENMP_FCFLAGS) $(OPT2) $(DMPI) $(DOMP) $(DFFT) 
COMPILERFLAGS3= $(STATIC) $(LD_MPI) $(OPENMP_FCFLAGS) $(OPT3) $(DMPI) $(DOMP) $(DFFT) 
if USE_CUDA
GPUCOMPILERFLAGS =  $(OPTGPU) $(DMPI) $(DOMP)
endif
#default compilation rules 
AM_FCFLAGS=$(COMPILERFLAGS1)


#default linker flags
LDLIBS =
LINKERFLAGS=-I. 

#Add Static linking
if USE_STATIC
LINKERFLAGS+=$(STATIC)
endif

#Add MPI linkerflags
if USE_MPI
LDLIBS+=$(LD_MPI)
endif

#Add FFTW Linkerflags
if CALLS_FFTW
LDLIBS+=$(LD_FFTW)
endif
#Add CUDA linkerflags
if USE_CUDA
LDLIBS+=$(LD_CUDA)
endif

DEPENDENCIES = define.h
########################################################################
#						 libraries list
########################################################################
# FFT libraries
if USE_NETLIB
OBJFFTPACK =fftpack.$(OBJEXT)
endif
if CALLS_FFTW
OBJFFTW =fftw.$(OBJEXT)
endif
if USE_CUDA
OBJCUDA= cuda.$(OBJEXT) cuda_alloc.$(OBJEXT)
GPU_COMPILER=nvcc
else
OBJCOULSOLV =coulsolv.$(OBJEXT)
endif
########################################################################
# 		        	Program name and complation rules
########################################################################

OBJKINPARA =  params.$(OBJEXT) kinetic.$(OBJEXT)

OBJCOMMON = $(OBJKINPARA)\
						restart.$(OBJEXT) restartc.$(OBJEXT)\
						static.$(OBJEXT) dynamic.$(OBJEXT) lda.$(OBJEXT) util.$(OBJEXT) abso_bc.$(OBJEXT) \
						pseudosoft.$(OBJEXT) pseudogoed.$(OBJEXT) ionmd.$(OBJEXT) forces.$(OBJEXT)\
 		        carlo.$(OBJEXT) localize.$(OBJEXT) localizer.$(OBJEXT)  subgrids.$(OBJEXT)  analyse.$(OBJEXT)\
						rho.$(OBJEXT) rhoc.$(OBJEXT) nonloc.$(OBJEXT) nonlocc.$(OBJEXT)\
						schmid.$(OBJEXT) zeroforce.$(OBJEXT) functions.$(OBJEXT)\
						loc_mfield.$(OBJEXT) givens.$(OBJEXT)\
						parallele.$(OBJEXT) expevol.$(OBJEXT) \
						pot_substrate.$(OBJEXT) forces_substrate.$(OBJEXT) md_substrate.$(OBJEXT)\
						short.$(OBJEXT) image.$(OBJEXT) lattice.$(OBJEXT) attachement.$(OBJEXT)\
						sicnew.$(OBJEXT) sicnewc.$(OBJEXT)\
						init.$(OBJEXT) generlcgo.$(OBJEXT)\
						2stUT.$(OBJEXT) orthmat.$(OBJEXT)
#sic-conditional objects:
if USE_TWOSTSIC
OBJCOMMON+= 2stUTc.$(OBJEXT) localize_rad.$(OBJEXT)
endif
# Programs
if USE_NETLIB
bin_PROGRAMS= pwteleman
pwteleman_SOURCES = main.F90
pwteleman_FCFLAGS = $(COMPILERFLAGS1) 
pwteleman_LDADD = $(OBJFFTPACK) $(OBJCOULSOLV) $(OBJCOMMON) $(LDLIBS)
pwteleman_LDFLAGS = $(LINKERFLAGS)
pwteleman-main.$(OBJEXT):  $(OBJFFTPACK) $(OBJCOULSOLV) $(OBJCOMMON)
endif
if USE_FFTW
bin_PROGRAMS = fftw_pwteleman
fftw_pwteleman_SOURCES = main.F90
fftw_pwteleman_FCFLAGS = $(COMPILERFLAGS1) 
fftw_pwteleman_LDADD = $(OBJFFTW) $(OBJCOULSOLV)  $(OBJCOMMON) $(LDLIBS)
fftw_pwteleman_LDFLAGS = $(LINKERFLAGS)
fftw_pwteleman-main.$(OBJEXT): $(OBJFFTW) $(OBJCOULSOLV)  $(OBJCOMMON)
endif
if USE_MKL
bin_PROGRAMS = mkl_pwteleman
mkl_pwteleman_SOURCES = main.F90
mkl_pwteleman_FCFLAGS = $(COMPILERFLAGS1)
mkl_pwteleman_LDADD = $(OBJFFTW) $(OBJCOULSOLV)  $(OBJCOMMON) $(LDLIBS)
mkl_pwteleman_LDFLAGS= $(LINKERFLAGS)
mkl_pwteleman-main.$(OBJEXT): $(OBJFFTW) $(OBJCOULSOLV)  $(OBJCOMMON)
endif
if USE_CUDA
bin_PROGRAMS= gpu_pwteleman
gpu_pwteleman_SOURCES = main.F90
gpu_pwteleman_FCFLAGS = $(COMPILERFLAGS1) 
gpu_pwteleman_LDADD = $(OBJCUDA) $(OBJCOMMON) $(LDLIBS)
gpu_pwteleman_LDFLAGS = $(LINKERFLAGS)
gpu_pwteleman-main.$(OBJEXT): $(OBJCUDA) $(OBJCOMMON) $(LDLIBS)
endif

########################################################################
#               Compilation rules for objects
########################################################################       
# twostr and twost modules are contained in 2stUT.o and 2stUTc.o, respectively. 
# Dependencies concerning objects are needed,  due to numerous #include statement 

# params.o
params.$(OBJEXT): params.F90 define.h surf.F90 pseudo.F90 Makefile
	$(FC) $(COMPILERFLAGS1) -c -o $@ $< 

# fftpack.o
fftpack.$(OBJEXT): fftpack.F90 define.h fftpack2.F90 params.$(OBJEXT) Makefile
	$(FC) $(COMPILERFLAGS2) -c -o $@ $<

# fftw.o
fftw.$(OBJEXT): fftw.F90 define.h fftw3.f03 params.$(OBJEXT) Makefile
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<

cuda_alloc.$(OBJEXT): cuda_alloc.F90 Makefile
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<

#kinetic.o
if USE_NETLIB
kinetic.$(OBJEXT): kinetic.F90 define.h fft.F90 findiff.F90 $(OBJFFTPACK) params.$(OBJEXT) Makefile
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<
endif
if CALLS_FFTW
kinetic.$(OBJEXT): kinetic.F90 define.h fft.F90 findiff.F90 $(OBJFFTW) params.$(OBJEXT)  Makefile
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<
endif
if USE_CUDA
kinetic.$(OBJEXT): kinetic.F90 define.h fft.F90 findiff.F90  $(OBJCUDA) params.$(OBJEXT) Makefile
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<
endif

# coulsolv.o
if USE_NETLIB
coulsolv.$(OBJEXT): coulsolv.F90 define.h falr.F90 coulex.F90 findiff-sinft.F90 $(OBJFFTPACK) params.$(OBJEXT) Makefile
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<
endif
if CALLS_FFTW
coulsolv.$(OBJEXT): coulsolv.F90 define.h falr.F90 coulex.F90 findiff-sinft.F90 $(OBJFFTW) $(OBJKINPARA) Makefile
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<
endif
# coulsolv does not work with cuda (yet ?)
#~ if USE_CUDA
#~ coulsolv.$(OBJEXT): coulsolv.F90 define.h falr.F90 coulex.F90 findiff-sinft.F90 $(OBJCUDA) params.o Makefile
#~ 	$(FC) $(COMPILERFLAGS1) -c -o $@ $< $(AM_LDFLAGS)
#~ endif

# sicnew.o
sicnew.$(OBJEXT): sicnew.F90 define.h $(OBJCOULSOLV) $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS1) $(IDRIS) -DREALSWITCH -c -o $@ $<

# sicnewc.o
sicnewc.$(OBJEXT): sicnew.F90 define.h $(OBJCOULSOLV) $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS1) $(IDRIS) -DCOMPLEXSWITCH -c -o $@ $<

# restart.o
if USE_TWOSTSIC
restart.$(OBJEXT): restart.F90 define.h 2stUTc.$(OBJEXT) 2stUT.$(OBJEXT) $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS3) $(IDRIS) -DREALSWITCH -c -o $@ $<
else
restart.$(OBJEXT): restart.F90 define.h $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS3) $(IDRIS) -DREALSWITCH -c -o $@ $<
endif

# restartc.o
restartc.$(OBJEXT): restart.F90 define.h $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS3) $(IDRIS) -DCOMPLEXSWITCH -c -o $@ $<

# init.o
init.$(OBJEXT): init.F90 define.h $(OBJCOULSOLV) $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<

# generlcgo.o
generlcgo.$(OBJEXT): generlcgo.F90 define.h params.$(OBJEXT)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<

# localize.o
localize.$(OBJEXT): localize.F90 define.h $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS1) $(IDRIS) -DCOMPLEXSWITCH -c -o $@ $<

# localizer.o
localizer.$(OBJEXT): localize.F90 define.h $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS1) $(IDRIS) -DREALSWITCH -c -o $@ $<

# localize_rad.o
localize_rad.$(OBJEXT): localize_rad.F90 define.h orthmat.$(OBJEXT) $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS1) $(IDRIS) -DREALSWITCH -c -o $@ $<

# static.o
if USE_CUDA
if USE_TWOSTSIC
static.$(OBJEXT): static.F90 define.h cuda_alloc.$(OBJEXT) 2stUT.$(OBJEXT) localize_rad.$(OBJEXT) pseudosoft.$(OBJEXT) $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS1) $(IDRIS) -DREALSWITCH -c -o $@ $<
else
static.$(OBJEXT): static.F90 define.h cuda_alloc.$(OBJEXT) pseudosoft.$(OBJEXT) $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<
endif
else
if USE_TWOSTSIC
static.$(OBJEXT): static.F90 define.h 2stUT.$(OBJEXT) localize_rad.$(OBJEXT) pseudosoft.$(OBJEXT) $(OBJKINPARA) $(OBJCOULSOLV)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<
else
static.$(OBJEXT): static.F90 define.h  pseudosoft.$(OBJEXT) $(OBJKINPARA) $(OBJCOULSOLV)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<
endif
endif

# orthmat.o
orthmat.$(OBJEXT): orthmat.F90 define.h params.$(OBJEXT)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<

# 2stUT.o
if USE_TWOSTSIC
2stUT.$(OBJEXT): 2stUT.F90 define.h localize_rad.$(OBJEXT)  orthmat.$(OBJEXT) $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS1) $(IDRIS) -DREALSWITCH -c -o $@ $<
else
2stUT.$(OBJEXT): 2stUT.F90 define.h  orthmat.$(OBJEXT) $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS1) $(IDRIS) -DREALSWITCH -c -o $@ $<
endif

# 2stUTc.o 
if USE_TWOSTSIC
2stUTc.$(OBJEXT): 2stUT.F90 define.h 2stUT.$(OBJEXT) orthmat.$(OBJEXT) $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS1) $(IDRIS) -DCOMPLEXSWITCH -c -o $@ $<
endif

# dynamic.o
if USE_TWOSTSIC
dynamic.$(OBJEXT): dynamic.F90 define.h  2stUTc.$(OBJEXT) 2stUT.$(OBJEXT) orthmat.$(OBJEXT) $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<
else
dynamic.$(OBJEXT): dynamic.F90 define.h $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<
endif
# lda.o
lda.$(OBJEXT): lda.F90 define.h $(OBJCOULSOLV) params.$(OBJEXT)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<

# util.o
util.$(OBJEXT): util.F90 define.h params.$(OBJEXT)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<

# abso_bc.o

abso_bc.$(OBJEXT): abso_bc.F90 define.h params.$(OBJEXT)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<

# pseudosoft.o
pseudosoft.$(OBJEXT): pseudosoft.F90 define.h $(OBJCOULSOLV) $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<

# pseudogoed.o
pseudogoed.$(OBJEXT): pseudogoed.F90 params.$(OBJEXT) define.h
	$(FC) $(COMPILERFLAGS1) -c -o $@ $< 

# ionmd.o
ionmd.$(OBJEXT): ionmd.F90 define.h $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<

# forces.o
forces.$(OBJEXT): forces.F90 define.h params.$(OBJEXT)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<

# carlo.o  #Beware, peraps implicitly dependent to parallele.F90 ?
carlo.$(OBJEXT): carlo.F90 define.h params.$(OBJEXT)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<

# subgrids.o
subgrids.$(OBJEXT): subgrids.F90 define.h params.$(OBJEXT)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<

# analyse.o
analyse.$(OBJEXT): analyse.F90 define.h params.$(OBJEXT)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<

# rho.o
rho.$(OBJEXT): rho.F90 define.h $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS1) $(IDRIS) -DREALSWITCH  -c -o $@ $<

# rhoc.o
rhoc.$(OBJEXT): rho.F90 define.h $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS1) $(IDRIS) -DCOMPLEXSWITCH -c -o $@ $<

# nonloc.o 
nonloc.$(OBJEXT): nonloc.F90 define.h params.$(OBJEXT)
	$(FC) $(COMPILERFLAGS1) $(IDRIS) -DREALSWITCH  -c -o $@ $<

# nonlocc.o 
nonlocc.$(OBJEXT): nonloc.F90 define.h params.$(OBJEXT)
	$(FC) $(COMPILERFLAGS1) $(IDRIS) -DCOMPLEXSWITCH -c -o $@ $<
	
# schmid.o
schmid.$(OBJEXT): schmid.F90 define.h params.$(OBJEXT)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<

# zeroforce.o
zeroforce.$(OBJEXT): zeroforce.F90 define.h $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<

# functions.o
functions.$(OBJEXT): functions.F90 define.h params.$(OBJEXT)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<
	
# loc_mfield.o
loc_mfield.$(OBJEXT): loc_mfield.F90 define.h $(OBJCOULSOLV) $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<

# givens.o 
givens.$(OBJEXT): givens.F90 define.h params.$(OBJEXT)
	$(FC) $(COMPILERFLAGS3) -c -o $@ $<

# parallele.o
parallele.$(OBJEXT): parallele.F90 define.h $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS3) -c -o $@ $<

# expevol.o
if USE_TWOSTSIC
expevol.$(OBJEXT): expevol.F90 define.h 2stUTc.$(OBJEXT) $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<
else
expevol.$(OBJEXT): expevol.F90 define.h  $(OBJKINPARA)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<
endif

# pot_substrate.o
pot_substrate.$(OBJEXT): pot_substrate.F90 define.h params.$(OBJEXT)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<

# forces_substrate.o
forces_substrate.$(OBJEXT): forces_substrate.F90 define.h params.$(OBJEXT)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<
	
# md_substrate.o
md_substrate.$(OBJEXT): md_substrate.F90 define.h params.$(OBJEXT)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<
	
# short.o
short.$(OBJEXT): short.F90 define.h params.$(OBJEXT)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<

# image.o
image.$(OBJEXT): image.F90 define.h $(OBJCOULSOLV) params.$(OBJEXT)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<
	
# lattice.o
lattice.$(OBJEXT): lattice.F90 define.h params.$(OBJEXT)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<

# attachement.o
attachement.$(OBJEXT): attachement.F90 define.h params.$(OBJEXT)
	$(FC) $(COMPILERFLAGS1) -c -o $@ $<

#cuda.o
cuda.$(OBJEXT): cuda.cu define_cuda.h values_gpu.h params_gpu.h static_gpu.cu
	$(GPU_COMPILER) $(GPUCOMPILERFLAGS) -c -o $@ $<



#~ initial_checks: makefile define.mk define.h
#~ 	@echo "#######################################################"
#~ 	@echo  Lists of known compilation parameters.
#~ 	@echo "#######################################################"
#~ 	@echo  Parameters defined in define.h:
#~ 	@echo "--------------------------------------------"
#~ 	@echo "gridfft       $(gridfft_value)"
#~ 	@echo "findiff       $(findiff_value)"
#~ 	@echo "numerov       $(numerov_value)"
#~ 	@echo "coufou        $(coufou_value)"
#~ 	@echo "coudoub       $(coudoub_value)"
#~ 	@echo "coudoub3D     $(coudoub3D_value)"
#~ 	@echo "twostsic      $(twostsic_value)"
#~ 	@echo "raregas       $(raregas_value)"
#~ 	@echo "lda_gpu       $(lda_gpu_value)"
#~ 	@echo "#######################################################"
#~ 	@echo  Parameters defined in makefile / define.mk:
#~ 	@echo "--------------------------------------------"
#~ 	@echo "CF90          $(CF90)"
#~ 	@echo "MPI_PARALLEL  $(MPI_PARALLEL)"
#~ 	@echo "OMP_THREADS   $(OMP_THREADS)"
#~ 	@echo "TYPE_FFT      $(TYPE_FFT)"
#~ 	@echo "DEBUG         $(DEBUG)"
#~ 	@echo "LINK_STATIC   $(LINK_STATIC)"
#~ 	@echo "MKL_THREADS   $(MKL_THREADS)"
#~ 	@echo "MKL_WRAPPERS  $(MKL_WRAPPERS)"
#~ 	@echo "MKLPATH       $(MKLPATH)"
#~ 	@echo "MKLINCLUDE    $(MKLINCLUDE)"
#~ 	@echo "#######################################################"
#~ 	@echo "-D flags enabled (1) by makefile rules:"
#~ 	@echo "--------------------------------------------"
#~ 	@echo "parayes       $(parayes_value)"
#~ 	@echo "parano        $(parano_value)"
#~ 	@echo "simpara       $(simpara_value)"
#~ 	@echo "netlib_fft    $(netlib_fft_value)"
#~ 	@echo "fftw_cpu      $(fftw_cpu_value)"
#~ 	@echo "fftw_gpu      $(fftw_gpu_value)"
#~ 	@echo "fftwnomkl     $(fftwnomkl_value)"
#~ 	@echo "paropenmp     $(paropenmp_value)"
#~ 	@echo "dynopenmp     $(dynopenmp_value)"
#~ 	@echo "#######################################################"

CLEANFILES=*.mod *.MOD
