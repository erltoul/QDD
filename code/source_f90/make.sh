#!/bin/bash

#export GXX_ROOT=/usr/lib64/gcc/x86_64-suse-linux/4.3/ #Compilation on Hyperion

if [ -z $1 ] ; then
    nb_procs=1
else
    nb_procs=$1
fi

if [ -z $2 ] ; then
    type_check=netlib
else
    type_check=$2
fi

#if [ -z $3 ] ; then
#    dim_check=1d
#else
#    dim_check=$3
#fi

if [ ! -f define.h ] || [ ! -f params.F90 ] || [ ! -f makefile ] ; then
    echo "some files are missing..."
    exit 0
fi

sed -i -e 's/netlib_fft.*/netlib_fft 0/' define.h
sed -i -e 's/fftw_cpu.*/fftw_cpu 0/' define.h
#sed -i -e 's/fftw3d_cpu.*/fftw3d_cpu 0/' define.h
sed -i -e 's/fftw_gpu.*/fftw_gpu 0/' define.h
#sed -i -e 's/fftw3d_gpu.*/fftw3d_gpu 0/' define.h

if [ $nb_procs = 1 ] ; then
    sed -i -e 's/parano.*/parano 1/' define.h
    sed -i -e 's/parayes.*/parayes 0/' define.h
    sed -i -e 's/hamdiag.*/hamdiag 1/' define.h
    sed -i -e 's/mdshort.*/mdshort 1/' define.h
    sed -i -e 's/madelungonly.*/madelungonly 1/' define.h
#    sed -i -e 's/knode=.*/knode=1/' params.F90
    sed -i -e 's/EXEC = essai\..*/EXEC = essai\.seq/' makefile
    sed -i -e 's/CF90    = .*/CF90    = IFORT/' makefile
    sed -i -e 's/USE_MPI = .*/USE_MPI = NO/' makefile
    if [ $type_check = netlib ] ; then
        sed -i -e 's/netlib_fft.*/netlib_fft 1/' define.h
        sed -i -e 's/TYPE_FFT = .*/TYPE_FFT = NETLIB/' makefile
#        sed -i -e 's/DIM = .*/DIM = 1d/' makefile
    fi
    if [ $type_check = fftw ] ; then
#    	if [ $dim_check = 1d ] ; then
#        	sed -i -e 's/DIM = .*/DIM = 1d/' makefile
        sed -i -e 's/fftw_cpu.*/fftw_cpu 1/' define.h
#    	fi
#    	if [ $dim_check = 3d ] ; then
#        	sed -i -e 's/DIM = .*/DIM = 3d/' makefile
#        	sed -i -e 's/fftw3d_cpu.*/fftw3d_cpu 1/' define.h
#    	fi
        sed -i -e 's/TYPE_FFT = .*/TYPE_FFT = FFTW/' makefile
    fi
    if [ $type_check = gpu ] ; then
#    	if [ $dim_check = 1d ] ; then
#        	sed -i -e 's/DIM = .*/DIM = 1d/' makefile
        sed -i -e 's/fftw_gpu.*/fftw_gpu 1/' define.h
#    	fi
#    	if [ $dim_check = 3d ] ; then
#        	sed -i -e 's/DIM = .*/DIM = 3d/' makefile
#        	sed -i -e 's/fftw3d_gpu.*/fftw3d_gpu 1/' define.h
#	fi
        sed -i -e 's/TYPE_FFT = .*/TYPE_FFT = cuFFT/' makefile
    fi
else
    sed -i -e 's/parano.*/parano 0/' define.h
    sed -i -e 's/parayes.*/parayes 1/' define.h
    sed -i -e 's/hamdiag.*/hamdiag 0/' define.h
    sed -i -e 's/mdshort.*/mdshort 0/' define.h
    sed -i -e 's/madelungonly.*/madelungonly 0/' define.h
    sed -i -e 's/knode=.*/knode='$nb_procs'/' params.F90
    sed -i -e 's/EXEC = essai\..*/EXEC = essai\.par/' makefile
    sed -i -e 's/CF90    = .*/CF90    = GFORT/' makefile
    sed -i -e 's/USE_MPI = .*/USE_MPI = YES/' makefile
fi

make clean
make

if [ $nb_procs = 1 ] ; then
    echo "the executable for sequential code is named 'essai.seq'"
else
    echo "the executable for parallel code is named 'essai.par'"
fi
