#!/bin/bash
NAME=essai

para=0
if [ $1 ] ; then
 para=$1
fi

if [ -z $2 ] ; then
    type_check=""
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

# For parallel setups: if mpif90 is available, use it; else, use ifort.
if [ $para = 1 ] || [ $para = 2 ] ; then
    if which mpif90 ; then
        PCF=MPIF90
    else
        PCF=IFORT
    fi
fi

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
        sed -i -e 's/LINK_STATIC = .*/LINK_STATIC = YES/' makefile
        ;;
esac

if [ $para = 0 ] ; then
    echo '*** serial compilation, simpara = no ***'
    echo "*** the executable for serial code is named '$NAME.seq' ***"
    sed -i -e 's/parano.*/parano 1/' define.h
    sed -i -e 's/parayes.*/parayes 0/' define.h
    sed -i -e 's/simpara.*/simpara 0/' define.h
    sed -i -e "/^ *#/!s/\(EXEC\) *=.*/\\1 = $NAME.seq/g" makefile
    sed -i -e 's/CF90    = .*/CF90    = IFORT/' makefile
    sed -i -e 's/USE_MPI = .*/USE_MPI = NO/' makefile
elif [ $para = 1 ] ; then
    echo '*** parallel compilation, simpara = no ***'
    echo "*** the executable for parallel code is name '$NAME.par' ***"
    sed -i -e 's/parano.*/parano 0/' define.h
    sed -i -e 's/parayes.*/parayes 1/' define.h
    sed -i -e 's/simpara.*/simpara 0/' define.h
    sed -i -e "/^ *#/!s/\(EXEC\) *=.*/\\1 = $NAME.par/g" makefile
    sed -i -e 's/CF90    = .*/CF90    = '$PCF'/' makefile
    sed -i -e 's/USE_MPI = .*/USE_MPI = YES/' makefile
elif [ $para = 2 ] ; then
    echo '*** parallel compilation, simpara = yes ***'
    echo "*** the executable for sim. parallel code is named '$NAME.sim' ***"
    sed -i -e 's/parano.*/parano 1/' define.h
    sed -i -e 's/parayes.*/parayes 0/' define.h
    sed -i -e 's/simpara.*/simpara 1/' define.h
    sed -i -e "/^ *#/!s/\(EXEC\) *=.*/\\1 = $NAME.sim/g" makefile
    sed -i -e 's/CF90    = .*/CF90    = '$PCF'/' makefile
    sed -i -e 's/USE_MPI = .*/USE_MPI = YES/' makefile
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
else
    echo "Empty second argument, preserving FFTW lib options."
    exit 0 
fi


make clean
make
make clean

