#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.69])
AC_INIT([pw-teleman], [73], [francois.labbe@univ-lemans.fr])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_MACRO_DIR([m4])
# Checks for programs.
#~ AC_LANG(Fortran)
#~ AC_FC_SRCEXT(F90)
#~ AC_FC_FUNC(cos)
########################################################################
#                         Check compiler
########################################################################
AC_ARG_WITH([compiler],
  [AS_HELP_STRING([--with-compiler],
  [Use specified compiler, choose among gfortran, ifort, mpif90 and xlf_r])],
  [])
#Choice of compiler
AC_CHECK_PROGS( [best_compiler], [ifort gfortran xlf_r])
if test -z $with_compiler ; then 
  #search for ifort as default compiler. If not found, search for gfortran or xlfr .
  echo "No compiler specified. Will use $best_compiler."
  FC="$best_compiler"
else 
  if test "x$with_compiler" = xyes; then
    # --with-compiler specified without a compiler name.
    AC_MSG_FAILURE([Please specify a compiler when using option --with-compiler. 
    Usage :  --with-compiler=<compilo> . 
    <compilo> can be gfortran, ifort, xlf_r or mpif90])
  else
    # check if required compiler exists on the system
    echo "user required use of $with_compiler"
    AC_CHECK_PROG( [found_compiler], [$with_compiler], [yes], [no])
  fi
  if test "x$found_compiler" = xno ; then
    AC_MSG_FAILURE([$with_compiler was not found on your system. Make sure $with_compiler is installed, or choose an other fortran compiler.])
  else 
    FC="$with_compiler"
  fi
fi

echo "will use $FC"
AC_PROG_FC([$FC])
#~ AC_PROG_FC([$with_compiler])
    #~ AC_MSG_ERROR([ /!\ need to check how MPI it makes its choice (and how to know what it choose)])


AM_CONDITIONAL([IFORT], [test "x$FC" = xifort])
AM_CONDITIONAL([GFORT], [test "x$FC" = xgfort])
AM_CONDITIONAL([XLF_R], [test "x$FC" = xxlf_r])
AM_CONDITIONAL([MPIF90], [test "x$FC" = xmpif90])
########################################################################
#                Special IDRIS flag for compiler xlf_r
########################################################################
if test "x$found_compiler" = xxlf_r ; then
  IDRIS=-WF
else 
  IDRIS=
fi
AC_SUBST([IDRIS]) 
########################################################################
#                   Check bitness and OS (according to compiler)
########################################################################
FL_CHECK_OS
########################################################################
#                   Check options and libraries
########################################################################
#=======================================================================
#                        Static compilation
#=======================================================================
#  Warning: may need -Wl,etc... if ill located library are used (cuda not in /usr/lib ...)
AC_ARG_WITH([static],
  [AS_HELP_STRING([--with-static],
    [--with-static : Use static compilation : All libraries are contained in the final executable. This makes a big executable.])],
    [],
    [with_static=no])
echo "Use static linking: $with_static"
AS_IF([test "x$with_static" = xyes],
        [],
        [test "x$with_static" = xno],
        [],
        [AC_MSG_FAILURE([--with-static=$with_static is not a valid option. Please use --with-static , --with-static=yes or --with-static=no.])]
      )
AM_CONDITIONAL([USE_STATIC],[test "x$with_static" = xyes])
#=======================================================================
# Choice of parallelism  (using --with-para)
#=======================================================================
AC_ARG_WITH([para],
            [AS_HELP_STRING([--with-para],
            [--with-para : use threads for FFT parallelism. Use --with-para=dyn for simultaneous mono jobs. Only for some compilers.])],
            [],
            [with_para=no])
echo "Use parallelism : $with_para"
AS_IF([test "x$with_para" != xno],
      [
        if test  "x$with_para" != xyes   &&  test "x$with_para" != xsim; then 
          AC_MSG_FAILURE([ $with_para is not a valid parallelism option. Use --with-para or --with-para=sim if you want some parallelism.])
        fi
        #nsert here flags for use of MPI with gfortran, irfot and xlf_r
      ])
#=======================================================================
# Check if OpenMP was asked using --with-openmp
#=======================================================================
AC_ARG_WITH([openmp],
  [AS_HELP_STRING([--with-openmp],
  [--with-openmp : Invoke OpenMP. --with-openmp : use threads for FFT.  --with-openmp=dyn : wave function parallelization.])],
  [],
  [with_openmp=no])
echo "Use OpenMp : $with_openmp"
#Flags for use of openmp
DOMP=
OPENMP_FCFLAGS=
AS_IF([test "x$with_openmp" != xno],
      [
        if test  "x$with_openmp" != xyes   &&  test "x$with_openmp" != xdyn; then 
          AC_MSG_FAILURE([ $with_openmp is not a valid parallelism option. Use --with-openmp or --with-openmp=dyn if you want wave function parallelization.])
        fi
        # Preprocessor flags ( -D<flag> )
        DOMP=-Dparopenmp
        if test "x$with_openmp" = xdyn; then 
          DOMP="$DOMP -Ddynopenmp"
        fi
        echo "DOMP: $DOMP"
        # OpenMP flags (-fopenmp, -openmp...)  in accordance with compiler in variable OPENMP_FCFLAGS :
        AC_OPENMP
      ])
#substitue variables with values in Makefile.in
AC_SUBST(DOMP)
AC_SUBST(OPENMP_FCFLAGS)
#~ AM_CONDITIONAL([USE_OMP],[test "x$with_fft" = xNETLIB])

#=======================================================================
# Check option --with-fft, sets fftw flags and -D flags for fft
#=======================================================================
# Read options 
AC_ARG_WITH([fft],
  [AS_HELP_STRING([--with-fft],
  [Choice of library for FFT: --with-fft=<choice> , with choice in : NETLIB, FFTW, MKL, CUFFT. If not set, default library is NETLIB.])],
  [],
  [with_fft=NETLIB])
echo "Use FFT : $with_fft"
AC_ARG_WITH([mklpath],
  [AS_HELP_STRING([--with-mklpath],
  [If MKL is the required FFT library, --with-mklpath=<path> specifies the path to MKL, if not set in MKLROOT.])],
  [mklpath="$with_mklpath"],
  [mklpath=no])
AC_ARG_WITH([wrappers],
  [AS_HELP_STRING([--with-wrappers],
  [If MKL is the required FFT library, --with-wrappers=<path> option specifies the path to intel fftw3 wrapppers library libfftw3xf_intel.a.  If path to libfftw3xf_intel.a is not given, wrappers are expected to be located in $MKLROOT/interfaces/fftw3xf. If option --with-wrapper is not present, mkl will use default wrappers for fftW3 anyway.]
  )],
  [],
  [with_wrappers=no])
AC_ARG_WITH([mkl-threads],
  [AS_HELP_STRING([--with-mklthreads],
  [If --with-fft=MKL, uses MKL threads. No effect with other FFT. ])],
  [mklthreads=yes],
  [mklthreads=no])
  echo "with wrappers = $with_wrappers"
  echo "with mkl threads = $mklthreads"
# Error messages if MKL threads and MKL wrappers are callded without MKL
AS_IF([test "x$with_fft" != "xMKL"], 
      [AS_IF([test "x$with_wrappers" != xno],
            AC_MSG_FAILURE([Option --with-wrappers is only available with option --with-fft=MKL]))
      AS_IF([test "x$mklthreads" != xno], 
            AC_MSG_FAILURE([Option --with-wrappers is only available with option --with-fft=MKL]))]
      )
# Creation of flags
MKL_IS_HERE=unknown
FFTW_IS_HERE=unknown
CUDA_IS_HERE=unknown
AS_IF(
      #.......................NETLIB....................................
      [test "x$with_fft" = xNETLIB],
      [      
        DFFT="-Dnetlib_fft"
      ],
      #........................MKL......................................
      [test "x$with_fft" = xMKL],
      [
      # PATH to MKL
      AS_IF([test "x$mklpath" != xno],
              [mklroot=${mklpath}],
              [mklroot=${MKLROOT}]
            )
      AS_IF([test "x$mklroot" = x],[AC_MSG_FAILURE([Path to MKL is empty. Load intel module or give MKL location using --with-mklpath=<path/to/mkl>])])
      # Linker flags for wrappers, includes and library
      AS_IF([test "x$with_wrappers" = xyes],
            [MKL_WRAPPERS="$mklroot/interfaces/fftw3xf"],
            [MKL_WRAPPERS="$with_wrappers"])
      MKL_INCLUDE="$mklroot/include"
      AS_IF([test "x$flag64" = xyes],
            [MKL_LIBPATH="$mklroot/lib/intel64"],
            [MKL_LIBPATH="$mklroot/lib/ia32"] 
              )
      #~ AS_IF([test "x$with_wrappers" != xno],
            #~ [
      CFLAGS="$CFLAGS "
      LDFLAGS="$LDFLAGS -L$MKL_WRAPPERS -L$MKL_LIBPATH -I$MKL_INCLUDE"
            #~ ],
            #~ [LDFLAGS="$LDFLAGS -L$MKL_LIBPATH -I$MKL_INCLUDE"])
             
      #Check common libraries
      AC_CHECK_LIB([dl], [main],[],[MKL_IS_HERE=no],[])
      AC_CHECK_LIB([m], [main],[],[MKL_IS_HERE=no])
      AC_CHECK_LIB([pthread], [main],[],[MKL_IS_HERE=no])
      AC_CHECK_LIB([iomp5], [main], [], [MKL_IS_HERE=no], [])      
      #check mkl libraries
      MKLFLAGS=
      AC_CHECK_LIB([mkl_core], [main], 
                    [MKLFLAGS="-lmkl_core $MKLFLAGS"
                      LIBS="-lmkl_core $LIBS"],
                    [MKL_IS_HERE=no], [-lmkl_sequential])
      AC_CHECK_LIB([mkl_intel_lp64], [main],
                    [MKLFLAGS="-lmkl_intel_lp64 $MKLFLAGS"
                      LIBS="-lmkl_intel_lp64 $LIBS"],
                    [MKL_IS_HERE=no], [-lmkl_sequential])
      echo "Use MKL threads: $mklthread"
      AS_IF([test "x$mklthreads" = xno],
        [AC_CHECK_LIB([mkl_sequential], [main], 
                      [MKLFLAGS="-lmkl_sequential  $MKLFLAGS"
                      LIBS="-lmkl_sequential $LIBS"],
                       [MKL_IS_HERE=no], [-lmkl_core])],
        [AC_CHECK_LIB([mkl_intel_thread], [main], [MKLFLAGS="-lmkl_intel_thread $MKLFLAGS"], [MKL_IS_HERE=no], [-lmkl_intel_lp64 -lmkl_core])]
        )
      AS_IF([test "$MKL_IS_HERE" = xno], 
            [AC_MSG_FAILURE([One or more library of MKL was not found or conflicted with an other library.])],
            [MKL_IS_HERE=yes]
            )
      # Static libraries : Transmit group flags to the linker
      AS_IF([test "x$with_static" = xyes],[MKLFLAGS="-Wl,--start-group $MKLFLAGS -Wl,--end-group "],[])
      # Check wrapper library
      #~ AS_IF([test "x$with_wrappers" != xno],
            #~ [
            AC_CHECK_LIB([fftw3xf_intel], [main], [],
                      [AC_MSG_FAILURE([Could not find -lfftw3xf_intel library. Make sure the library is located in $MKLROOT, or give path to wrappers using --with-wrappers=<path/to/libfftw3xf_intel.a> .  Instructions on how to build your own wrappers library can be found on Intel Website : https://software.intel.com/en-us/node/522277], [])]
                      )
            echo "$MKLFLAGS"
             LD_FFTW="-lfftw3xf_intel $MKLFLAGS -liomp5 -lpthread -lm -ldl"
             #~ ],
            #~ [LD_FFTW="$MKLFLAGS -liomp5 -lpthread -lm -ldl"])
      AC_SUBST([LD_FFTW])
      DFFT="-Dfftw_cpu"
      ],
      #........................FFTW.....................................
      [test "x$with_fft" = xFFTW],
      [    
        # Switches used to preprocess files
        DFFT="-Dfftw_cpu -Dfftwnomkl"
        # Flags for FFTW
        if test "x$with_openmp" = xyes || test "x$with_openmp" = xdyn; then
          # with openmp
          fftw_flags="-lfftw3_omp"
          AC_CHECK_LIB([fftw3_omp], [main], [FFTW_IS_HERE=yes], 
              [AC_MSG_FAILURE([FFTW was asked, but test for $fftw_flags failed.
                    Please run configure with option LDFLAGS=-L<path to fftw>, or install fftw in default library path.])]
          )
        else
          # without openmp
          fftw_flags="-lfftw3 -lm"
          AC_CHECK_LIB([fftw3], [main],
              [AC_CHECK_LIB([m], [main], [FFTW_IS_HERE=yes],
                  [AC_MSG_FAILURE([FFTW was asked, but test for -lm failed.])
                  ])
              ],
              [AC_MSG_FAILURE([FFTW was asked, but test for -lfftw3 failed. Please run configure with option LDFLAGS=-L<path to fftw>, or install fftw in default library path.])]
          )
        fi
        LD_FFTW="$fftw_flags"
        AC_SUBST([LD_FFTW])
        echo "FFTW flags : $LD_FFTW" 
      ],
      #........................CUDA.....................................
      [test "x$with_fft" = xCUFFT],
      [
        # Ecrire macro de localisation de cuda.
        FL_LOCATE_CUDA
        AC_CHECK_LIB([stdc++], [main], [], [AC_MSG_FAILURE([Standart C++ library -lstdc++ library was not found.])])
        AC_SUBST([LD_CUDA], ["-lcudart -lcufft -lstdc++"])
        DFFT="-Dfftw_gpu"
        #insert tests for cudart, etc...
      ],
      #......................INVALID....................................
      [
        AC_MSG_FAILURE([ Unvalid FFT option. --with-fft=<choice> .  <choice> can can be : NETLIB, FFTW, MKL or CUFFT. If not set, default library is NETLIB. ])
      ]
      )

LIBS=
AC_SUBST(DFFT)
AM_CONDITIONAL([USE_NETLIB],[test "x$with_fft" = xNETLIB])
AM_CONDITIONAL([USE_FFTW],[test "x$FFTW_IS_HERE" = xyes ])
AM_CONDITIONAL([USE_MKL], [test "x$MKL_IS_HERE" = xyes ])
AM_CONDITIONAL([USE_CUDA],[test "x$CUDA_IS_HERE" = xyes ])
AM_CONDITIONAL([CALLS_FFTW],[test "x$FFTW_IS_HERE" = xyes || test "x$MKL_IS_HERE" = xyes ])

# par defaut :
# a modifier avec des conditions selon bibliotheque utilisée
#~ AC_SUBST(DFFTW, "-Dfftw_cpu")
#=======================================================================
AC_ARG_WITH([mpi],
  [AS_HELP_STRING([--with-mpi],
  [mpi not ready yet.])],
  [],
  [with_mpi=no])
echo "Use MPI : $with_mpi"
AM_CONDITIONAL([USE_MPI],[test "x$with_mpi" = xyes ])

#~ LD_MPI=-lmpi
LD_MPI=
#~ DMPI=
DMPI=-Dparano
AC_SUBST(DMPI)
AC_SUBST(LD_MPI)

AC_SUBST(DIDRIS,[])
# inserrer condition netlib, fftw, cuda, etc..."
#~ AC_SUBST(DFFTWADD)
# Checks for header files.
AC_CONFIG_HEADERS([define.h])
#~ AC_CONFIG_HEADERS([config.h])
# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_PROG_RANLIB
AM_PROG_AR
AC_CONFIG_FILES([Makefile
                 code/source_f90/Makefile])
AC_OUTPUT
