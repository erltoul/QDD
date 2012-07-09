#!/bin/bash
NAME=essai

para=0
if [ $1 ] ; then
 para=$1
fi

if [ ! -f define.h ] || [ ! -f params.F90 ] || [ ! -f makefile ] ; then
    echo "some files are missing..."
    exit 0
fi

if [ $para = 1 ] || [ $para = 2 ] ; then
    if which mpif90 ; then
        PCF=MPIF90
    else
        PCF=IFORT
    fi
fi

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
    echo '* parallel compilation, simpara = no ***'
    echo "* the executable for parallel code is name '$NAME.par' ***"
    sed -i -e 's/parano.*/parano 0/' define.h
    sed -i -e 's/parayes.*/parayes 1/' define.h
    sed -i -e 's/simpara.*/simpara 0/' define.h
    sed -i -e "/^ *#/!s/\(EXEC\) *=.*/\\1 = $NAME.par/g" makefile
    sed -i -e 's/CF90    = .*/CF90    = '$PCF'/' makefile
    sed -i -e 's/USE_MPI = .*/USE_MPI = YES/' makefile
elif [ $para = 2 ] ; then
    echo '* parallel compilation, simpara = yes ***'
    echo "* the executable for sim. parallel code is named '$NAME.sim' ***"
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

make clean
make
make clean
