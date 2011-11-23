#!/bin/bash

if [ -z $1 ] ; then
    nb_procs=1
else
    nb_procs=$1
fi

if [ ! -f define.h ] || [ ! -f params.F90 ] || [ ! -f makefile ] ; then
    echo "some files are missing..."
    exit 0
fi

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
else
    sed -i -e 's/parano.*/parano 0/' define.h
    sed -i -e 's/parayes.*/parayes 1/' define.h
    sed -i -e 's/hamdiag.*/hamdiag 0/' define.h
    sed -i -e 's/mdshort.*/mdshort 0/' define.h
    sed -i -e 's/madelungonly.*/madelungonly 0/' define.h
    sed -i -e 's/knode=.*/knode='$nb_procs'/' params.F90
    sed -i -e 's/EXEC = essai\..*/EXEC = essai\.par/' makefile
    sed -i -e 's/CF90    = .*/CF90    = MPIF90/' makefile
    sed -i -e 's/USE_MPI = .*/USE_MPI = YES/' makefile
fi

make clean
make

if [ $nb_procs = 1 ] ; then
    echo "the executable for sequential code is named 'essai.seq'"
else
    echo "the executable for parallel code is named 'essai.par'"
fi