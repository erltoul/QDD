#!/bin/bash

# $1 parallele  0 1 2
# $2 type_check:  gpu netlib fftw mkl
# $3 debug: none debug
# $4 mkl: none mkl


NAME=$TELEMANNAME
if [ -z $NAME ]  ; then
	NAME=essai
	echo "no name given - prefix changed to $NAME"
fi


para=0
if [ $1 ] ; then
 para=$1
fi

if [ -z $2 ] ; then
    type_check=none
else
    type_check=$2
fi

debug=0
openmp=0
mkl=0
if [ -z $3 ] ; then
    debug=none
else
	if [ $3 = "debug" ] ; then
    		debugebug=$3
    	        echo "DEBUG is on"
	fi
if [ $3 = "openmp" ] ; then
	openmp=$3
        echo "OPENMP is on - set OMP_NUM_THREADS"
fi
fi
if [ -z $4 ] ; then
	mkl=none 
else
if [ $4 = "mkl" ] ; then
	mkl=$4
        echo "MKL_OPENMP is on- set MKL_NUM_THREADS"
fi
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

export PCF=GFORT
if which ifort ; then
        export PCF=IFORT
fi


# For parallel setups: if mpif90 is available, use it; else, use ifort.
if [ $para = 1 ] || [ $para = 2 ] ; then
    if which mpif90 ; then
        export PCF=MPIF90
    else
        export PCF=GFORT
    fi
fi

echo 
# This should detect the Message Passing Toolkit (MPT) of SGI.  MPT
# provides SGI's own implementation of MPI, which cannot be linked
# with -static.
case "$LD_LIBRARY_PATH" in
    *sgi?mpt*)
        if [ $para = 1 ] || [ $para = 2 ] ; then
            echo '*** SGI cluster detected, LINK_STATIC set to NO ***'
            sed -i -e 's/LINK_STATIC = .*/LINK_STATIC = NO/' makefile
        elif [ $para = 0 ] ; then
            sed -i -e 's/LINK_STATIC = .*/LINK_STATIC = YES/' makefile
        fi
        ;;
    *)
#        sed -i -e 's/LINK_STATIC = .*/LINK_STATIC = YES/' makefile
        if [ $type_check = gpu ] ; then
            sed -i -e 's/LINK_STATIC = .*/LINK_STATIC = NO/' makefile
        fi
        ;;
esac

if [ $debug = debug ] ; then
    sed -i -e 's/DEBUG = .*/DEBUG = YES/' makefile
else
    sed -i -e 's/DEBUG = .*/DEBUG = NO/' makefile
fi
if [ $openmp = openmp ] ; then
    sed -i -e 's/OMP_THREADS = .*/OMP_THREADS = DYN/' makefile
else
    sed -i -e 's/OMP_THREADS = .*/OMP_THREADS = NO/' makefile
fi
if [ $mkl = mkl ] ; then
    sed -i -e 's/MKL_THREADS = .*/MKL_THREADS = YES/' makefile
else
    sed -i -e 's/MKL_THREADS = .*/MKL_THREADS = NO/' makefile
fi





if [ $para = 0 ] ; then
    echo '*** serial compilation, simpara = no ***'
    echo "*** the executable for serial code is named '$NAME.seq' ***"
#    sed -i -e 's/parano[[:space:]]*[0-9]\+/parano 1/' define.h
#    sed -i -e 's/parayes.*/parayes 0/' define.h
#    sed -i -e 's/simpara.*/simpara 0/' define.h
    sed -i -e "/^ *#/!s/\(EXEC\) *=.*/\\1 = $NAME.seq/g" makefile
#    sed -i -e 's/CF90 = .*/CF90 = GFORT/' makefile
    sed -i -e 's/CF90 = .*/CF90 = '$PCF'/' makefile
    sed -i -e 's/MPI_PARALLEL = .*/MPI_PARALLEL = NO/' makefile
elif [ $para = 1 ] ; then
    echo '*** parallel compilation, simpara = no ***'
    echo "*** the executable for parallel code is name '$NAME.par' ***"
#    sed -i -e 's/parano[[:space:]]*[0-9]\+/parano 0/' define.h
#    sed -i -e 's/parayes.*/parayes 1/' define.h
#    sed -i -e 's/simpara.*/simpara 0/' define.h
    sed -i -e "/^ *#/!s/\(EXEC\) *=.*/\\1 = $NAME.par/g" makefile
    sed -i -e 's/CF90 = .*/CF90 = '$PCF'/' makefile
    sed -i -e 's/MPI_PARALLEL = .*/MPI_PARALLEL = YES/' makefile
elif [ $para = 2 ] ; then
    echo '*** parallel compilation, simpara = yes ***'
    echo "*** the executable for sim. parallel code is named '$NAME.sim' ***"
#    sed -i -e 's/parano[[:space:]]*[0-9]\+/parano 1/' define.h
#    sed -i -e 's/parayes.*/parayes 0/' define.h
#    sed -i -e 's/simpara.*/simpara 1/' define.h
    sed -i -e "/^ *#/!s/\(EXEC\) *=.*/\\1 = $NAME.sim/g" makefile
    sed -i -e 's/CF90 = .*/CF90 = '$PCF'/' makefile
    sed -i -e 's/MPI_PARALLEL = .*/MPI_PARALLEL = SIM/' makefile
else
    echo "(0) serial; (1) parallel; (2) simpara"
    exit 0 
fi

if [ $type_check = netlib ] ; then
    sed -i -e 's/netlib_fft.*/netlib_fft 1/' define.h
    sed -i -e 's/TYPE_FFT = .*/TYPE_FFT = NETLIB/' makefile
#        sed -i -e 's/DIM = .*/DIM = 1d/' makefile
elif [ $type_check = fftw ] ; then
#    	if [ $dim_check = 1d ] ; then
#        	sed -i -e 's/DIM = .*/DIM = 1d/' makefile
    sed -i -e 's/fftw_cpu.*/fftw_cpu 1/' define.h
#    	fi
#    	if [ $dim_check = 3d ] ; then
#        	sed -i -e 's/DIM = .*/DIM = 3d/' makefile
#        	sed -i -e 's/fftw3d_cpu.*/fftw3d_cpu 1/' define.h
#    	fi
    sed -i -e 's/TYPE_FFT = .*/TYPE_FFT = FFTW/' makefile
elif [ $type_check = mkl ] ; then
    sed -i -e 's/fftw_cpu.*/fftw_cpu 1/' define.h
    sed -i -e 's/TYPE_FFT = .*/TYPE_FFT = MKL/' makefile
elif [ $type_check = gpu ] ; then
    sed -i -e 's/fftw_cpu.*/fftw_gpu 1/' define.h
    sed -i -e 's/TYPE_FFT = .*/TYPE_FFT = cuFFT/' makefile
    cp define.h define_cuda.h
    sed -i -e 's/!/\/\//g' define_cuda.h #define_cuda.h is just define.h turned into C++
else
    echo "Empty second argument, preserving FFTW lib options."
    exit 0 
fi


#make clean
make
#make clean

