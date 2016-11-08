#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.69])
AC_INIT([pw-teleman], [73], [francois.labbe@univ-lemans.fr])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_MACRO_DIR([m4])
# Remark : Unable to properly use Fortran options of Autotools, due to preprocessor flags in the code.
# Preprocessor use is NOT natural in fortran.
#~ AC_LANG(Fortran)
#~ AC_FC_SRCEXT(F90)
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
  echo "No compiler specified by user."
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
    #~ AC_MSG_ERROR([ /!\ need to check how MPI makes its choice (and how to know what it choose)])


AM_CONDITIONAL([IFORT], [test "x$FC" = xifort])
AM_CONDITIONAL([GFORT], [test "x$FC" = xgfortran])
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
#                   Read #define in header
########################################################################
# Check for header.
FL_READ_HEADER([code/source_f90/define.h])
echo "#################################################################"
echo "Parameters defined in define.h:"
echo "-------------------------------------"
echo "gridfft:      $gridfft"
echo "findiff:      $findiff"
echo "numerov:      $numerov"
echo "coufou:       $coufou"
echo "coudoub:      $coudoub"
echo "coudoub3D:    $coudoub3D"
echo "twostsic:     $twostsic"
echo "cmplxsic:     $cmplxsic"
echo "raregas:      $raregas"
echo "lda_gpu:      $lda_gpu"
echo "asynclaunch:  $asynclaunch"
echo "#################################################################"
# Check  coherence of options
# gridfft, findfft, numerov : only one can be enalble
AS_IF([test "x$gridfft" = x1], 
      [
      AS_IF([test "x$findiff" = x1 ], [AC_MSG_FAILURE([gridfft and findiff cannot be enabled simulatneously in define.h])])
      AS_IF([test "x$numerov" = x1 ], [AC_MSG_FAILURE([gridfft and numerov cannot be enabled simulatneously in define.h])])
      ],
      [
      AS_IF([test "x$findiff" = x0], 
            [
            AS_IF([test "x$numerov" = x0 ], [AC_MSG_FAILURE([gridfft, findiff and numerov cannot be simulatneously disabled in define.h])])
            ]
            )
      ])
# coufou, coudoub: only one can be enabled
AS_IF([test "x$coufou" = x1], 
      [
      AS_IF([test "x$coudoub" = x1 ], [AC_MSG_FAILURE([coufou and coudoub cannot be enabled simulatneously in define.h])])
      ])
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
#                               MPI
#=======================================================================
#Choice of parallelism
AC_ARG_WITH([para],
            [AS_HELP_STRING([--with-para],
            [--with-para : Use MPI for FFT parallelism. Use --with-para=dyn for simultaneous mono jobs. Only for some compilers.])],
            [],
            [with_para=no])
#check and flags
AS_IF([test "x$with_para" != xno],
      [
      AS_IF([test  "x$with_para" != xyes],
            [AS_IF([test  "x$with_para" != xsim],
                  [AC_MSG_FAILURE([ $with_para is not a valid parallelism option. Use --with-para or --with-para=sim if you want some parallelism.])]
                  )
            ])
      LD_MPI=
      AS_IF([test "x$FC" = xifort], [LD_MPI="-lmpi"])
      AS_IF([test "x$FC" = xxlf_r],
            [
            FC=mpxlf90_r
            LD_MPI="-lpesslsmp -lesslsmp"
            ])
      AC_SUBST(LD_MPI)
      ])
AS_IF([test "x$with_para" = xno],
      [DMPI="-Dparano"],
      [test "x$with_para" = xyes],
      [DMPI="-Dparayes"],
      [test "x$with_para" = xsim],
      [DMPI="-Dparano -Dsimpara"],
      [DMPI=""])
AC_MSG_NOTICE([Using parallelism : $with_para])
AC_MSG_NOTICE([Compiler taking MPI into account: $FC])
AM_CONDITIONAL([USE_MPI],[test "x$with_para" != xno])
AC_SUBST(DMPI)
#=======================================================================
# Choice of OpenMP
#=======================================================================
AC_ARG_WITH([openmp],
  [AS_HELP_STRING([--with-openmp],
  [--with-openmp : Invoke OpenMP. --with-openmp : use threads for FFT.  --with-openmp=dyn : wave function parallelization.])],
  [],
  [with_openmp=no])
echo "Use OpenMp : $with_openmp"
#check collision of MPI and OMP
AS_IF([test "x$with_para" != xno],
      [AS_IF([test "x$with_openmp" != xno], 
        [
        AC_MSG_FAILURE([Code not ready for simultaneous MPI and OpenMP parallelizations])
        ])
      ])
#Flags for use of openmp
DOMP=
OPENMP_FCFLAGS=
AS_IF([test "x$with_openmp" != xno],
      [
      AS_IF([test "x$with_openmp" != xyes], 
            [
            AS_IF([test "x$with_openmp" != xdyn],
                  [
                  AC_MSG_FAILURE([$with_openmp is not a valid parallelism option. Use --with-openmp or --with-openmp=dyn if you want wave function parallelization.])
                  ])
            ])
# Preprocessor flags ( -D<flag> )
      DOMP="-Dparopenmp"
      AS_IF([test "x$with_openmp" = xdyn], [DOMP="$DOMP -Ddynopenmp"])
# OpenMP flags (-fopenmp, -openmp...)  in accordance with compiler in variable OPENMP_FCFLAGS :
      ])
#substitue variables with values in Makefile.in
AC_LANG(Fortran)
AC_OPENMP
# Because Autotools can't properlly handle library tests :
AC_LANG(C)
AC_SUBST(DOMP)
AC_SUBST(OPENMP_FCFLAGS)
echo "OPENMP_FCFLAGS : $OPENMP_FCFLAGS"
#=======================================================================
# Check option --with-fft, sets fftw flags and -D flags for fft
#=======================================================================
# Read options 
#FFT option:
AC_ARG_WITH([fft],
  [AS_HELP_STRING([--with-fft],
    [Choice of library for FFT: --with-fft=<choice> , with choice in : NETLIB, FFTW, MKL, CUFFT. If not set, default library is NETLIB.])
  ],
  [
    AS_IF([test "x$gridfft" = x0], [AC_MSG_FAILURE([Use of option --with-fft is incompatible with gridfft=0 set in define.h])])
    AS_IF([test "x$with_fft" = xyes],
      [
      with_fft=NETLIB
      AC_MSG_NOTICE([Ignoring empty option --with-fft. Using default.])
      ])
  ],
  [
    AS_IF([test "x$gridfft" = x1], 
          [with_fft=NETLIB],
          [with_fft=no])
  ])
AS_IF([test "x$with_fft" != xno], [AC_MSG_NOTICE([Use FFT : $with_fft]) ])
# coudoub3D only enabled when using FFTW cals:
AS_IF([test "x$coudoub3D" = x1 ],
      [AS_IF([test "x$with_fft" != xFFTW],
              [AS_IF([test "x$with_fft" != xMKL],
                      [
                      AC_MSG_FAILURE([coudoub3D requires the use of FFTW (or MKL with fftw wrappers)])
                      ])
              ])
      ])
# lda_gpu is only enabled when using  CUFFT calls
AS_IF([test "x$lda_gpu" = x1 ],
      [AS_IF([test "x$with_fft" != xCUFFT],
              [
              AC_MSG_FAILURE([lda_gpu requires the use of GPU])
              ])
      ])
#MKL path option:
AC_ARG_WITH([mklpath],
  [AS_HELP_STRING([--with-mklpath],
  [If MKL is the required FFT library, --with-mklpath=<path> specifies the path to MKL, if not set in MKLROOT.])],
  [
  AS_IF([test "x$with_fft" != xMKL], 
        [
        AC_MSG_NOTICE([ Call for option --with-mklpath is irrelevant.] ) 
        AC_MSG_NOTICE([ --with-mklpath shall be used to specify mkl path only when mkl is used as FFT library] )
        AC_MSG_NOTICE([ Ignoring option --with-mklpath] ) 
        ])
  mklpath="$with_mklpath"
  ],
  [mklpath=no])
# MKL wrappers option
AC_ARG_WITH([wrappers],
  [AS_HELP_STRING([--with-wrappers],
  [If MKL is the required FFT library, --with-wrappers=<path> option specifies the path to intel fftw3 wrapppers library libfftw3xf_intel.a.  If path to libfftw3xf_intel.a is not given, wrappers are expected to be located in $MKLROOT/interfaces/fftw3xf. If option --with-wrapper is not present, mkl will use default wrappers for fftW3 anyway.])],
  [
  AS_IF([test "x$with_fft" != xMKL], 
      [
      AC_MSG_NOTICE([ Call for option --with-wrappers is irrelevant.] ) 
      AC_MSG_NOTICE([ --with-wrappers shall be used to specify mkl wrappers path only when mkl is used as FFT library] )
      AC_MSG_NOTICE([ Ignoring option --with-wrappers] ) 
      ])
  ],
  [with_wrappers=no])
# MKL threads option
AC_ARG_WITH([mklthreads],
  [AS_HELP_STRING([--with-mklthreads],
  [If --with-fft=MKL, uses MKL threads. No effect with other FFT. ])],
  [
  mklthreads=yes
  AS_IF([test "x$with_fft" != xMKL], 
    [
    AC_MSG_NOTICE([ Call for --with-mkl-threads is irrelevant.] ) 
    AC_MSG_NOTICE([ --with-mkl-threads shall be used only when mkl is used as FFT library] )
    AC_MSG_NOTICE([ Ignoring option --with-mkl-threads] ) 
    ])
  ],
  [mklthreads=no])
  echo "with wrappers = $with_wrappers"
  echo "with mklthreads = $mklthreads"
#~ # Error messages if MKL threads and MKL wrappers are callded without MKL
#~ AS_IF([test "x$with_fft" != "xMKL"], 
      #~ [AS_IF([test "x$with_wrappers" != xno],
            #~ AC_MSG_FAILURE([Option --with-wrappers is only available with option --with-fft=MKL]))
      #~ AS_IF([test "x$mklthreads" != xno], 
            #~ AC_MSG_FAILURE([Option --with-wrappers is only available with option --with-fft=MKL]))]
      #~ )
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
      AS_IF([test "x$mklthreads" = xyes],
            [AC_CHECK_LIB([iomp5], [main], [], [MKL_IS_HERE=no], [])
            ])
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
      AS_IF([test "x$mklthreads" = xno],
            [AC_CHECK_LIB([mkl_sequential], [main], 
                          [MKLFLAGS="-lmkl_sequential  $MKLFLAGS"
                          LIBS="-lmkl_sequential $LIBS"],
                           [MKL_IS_HERE=no], [-lmkl_core])],
            [AC_CHECK_LIB([mkl_intel_thread], [main], 
                          [
                          MKLFLAGS="-lmkl_intel_thread $MKLFLAGS"
                          LIBS="-lmkl_intel_thread $LIBS"
                          ], [MKL_IS_HERE=no], [-lmkl_intel_lp64 -lmkl_core])]
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
                      [AC_MSG_FAILURE([Could not find -lfftw3xf_intel library. Make sure the library is located in $MKLROOT/interfaces/fftw3xf, or give path to wrappers using --with-wrappers=<path/to/libfftw3xf_intel.a> .  Instructions on how to build your own wrappers library can be found on Intel Website : https://software.intel.com/en-us/node/522277], [])]
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
      #......................NoFFT......................................
      [test "x$with_fft" = xno],
      [
        AC_MSG_NOTICE([No FFT needed])
      ],
      #......................INVALID....................................
      [
        AC_MSG_FAILURE([Unvalid FFT option. --with-fft=<choice> .  <choice> can can be : NETLIB, FFTW, MKL or CUFFT. If not set, default library is NETLIB. ])
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


AC_SUBST(DIDRIS,[])
AC_CONFIG_HEADERS([config.h])
AC_PROG_RANLIB
AM_PROG_AR
AC_CONFIG_FILES([Makefile
                 code/source_f90/Makefile])
AC_OUTPUT
