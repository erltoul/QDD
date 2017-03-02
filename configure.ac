#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.69])
AC_INIT([pw-teleman], [022017], [francois.labbe@univ-lemans.fr])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_MACRO_DIR([m4])
# Remark : Unable to properly use Fortran options of Autotools, due to preprocessor flags in the code.
# Preprocessor use is NOT natural in fortran.
#~ AC_LANG(Fortran)
#~ AC_FC_SRCEXT(F90)
########################################################################
# Prevent configure tu set -g as default FCFLAGS
########################################################################
AS_IF([test "x$FCFLAGS" = x], [FCFLAGS=])
########################################################################
#                         Check compiler
########################################################################
AC_ARG_WITH([compiler],
  [AS_HELP_STRING([--with-compiler],
  [Use specified compiler, choose among gfortran, ifort, mpifort (mpif90) and xlf_r])],
  [])
#Choice of compiler
#search for ifort as default compiler. If not found, search for gfortran or xlfr .
AS_IF([test -z $with_compiler],
      [
      AC_MSG_NOTICE([No compiler specified by user.])
      AC_CHECK_PROGS( [best_compiler], [ifort gfortran xlf_r])
      FC="$best_compiler"
      ],
      [
  AS_IF([test "x$with_compiler" = xyes],
        [AC_MSG_FAILURE([Please specify a compiler when using option --with-compiler. Usage :  --with-compiler=<compilo> . <compilo> can be gfortran, ifort, xlf_r or mpifort (mpif90)])],
        [
        AC_MSG_CHECKING([for $with_compiler, as User asked])
        AC_CHECK_PROG( [found_compiler], [$with_compiler], [yes], [no])
        ]
        )
  AS_IF([test "x$found_compiler" = xno],
        [AC_MSG_FAILURE([$with_compiler was not found on your system. Make sure $with_compiler is installed, or choose an other fortran compiler.])],
        [FC="$with_compiler"]
        )
  ])
AC_PROG_FC([$FC])
AM_CONDITIONAL([IFORT], [test "x$FC" = xifort])
AM_CONDITIONAL([GFORT], [test "x$FC" = xgfortran])
AM_CONDITIONAL([XLF_R], [test "x$FC" = xxlf_r])
AM_CONDITIONAL([MPIFORT], [test "x$FC" = xmpifort || test "x$FC" = xmpif90])

#Find the underlying compiler used by mpifort:
FCOMP=unknown
AS_IF([test "x$FC" = xmpifort || test "x$FC" = xmpif90  ],
      [
      $FC 2>&1 | grep -i "ifort" > /dev/null 2>&1 && FCOMP=ifort
      $FC 2>&1 | grep -i "gfortran" > /dev/null 2>&1 && FCOMP=gfortran
      $FC 2>&1 | grep -i "xlf_r" > /dev/null 2>&1 && FCOMP=xlf_r
      AS_IF([test "x$FCOMP" = xunknown],
            [AC_MSG_FAILURE([Could not identify compiler $true_comp used by mpifort. ])],
            [AC_MSG_NOTICE([$FC uses $FCOMP as compiler])]
            )
      ])
AM_CONDITIONAL([MPI_IFORT], [test "x$FCOMP" = xifort])
AM_CONDITIONAL([MPI_GFORT], [test "x$FCOMP" = xgfortran])
AM_CONDITIONAL([MPI_XLF_R], [test "x$FCOMP" = xxlf_r])

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
# Check  coherence of header options
# gridfft, findfft, numerov : only one can be enabled
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
#                        Debug option
#=======================================================================
AC_ARG_WITH([debug],
    [AS_HELP_STRING([--with-debug],[--with-debug : Use compiler specfic debug option during compilation])
    ],
  [],
  [with_debug=no])
AC_MSG_NOTICE([Use Debug option at compile time: $with_debug])
AM_CONDITIONAL([USE_DEBUG],[test "x$with_debug" = xyes])
#=======================================================================
#                 Deprecated optimization option
#=======================================================================
AC_ARG_WITH([xW],
    [AS_HELP_STRING([--with-xW],[--with-xW : use deprecated optimization option -xW instead of -axSSE4.2 -msse4.2])
    ],
  [],
  [with_xW=no])
AC_MSG_NOTICE([Use deprecated optimization -xW at compile time: $with_xW])
AM_CONDITIONAL([USE_XW],[test "x$with_xW" = xyes])
#=======================================================================
#                        Static compilation
#=======================================================================
AC_ARG_WITH([static],
  [AS_HELP_STRING([--with-static],
    [--with-static : Use static compilation : All libraries are contained in the final executable. This makes a big executable.])],
  [],
  [with_static=no])
AC_MSG_NOTICE([Use static linking: $with_static])
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
            [--with-para : Use MPI for FFT parallelism. Use --with-para=sim for simultaneous mono jobs.  --withpara=world for Crank Nicolson parallel run. Only for some compilers.])],
            [],
            [with_para=no])
#check and flags
AS_IF([test "x$with_para" != xno],
      [
      AS_IF([test  "x$with_para" != xyes],
          [AS_IF([test  "x$with_para" != xworld],
                  [AS_IF([test  "x$with_para" != xsim],
                        [AC_MSG_FAILURE([$with_para is not a valid parallelism option. Use --with-para or --with-para=sim if you want some parallelism.])]
                        )
                ])
          ])
      LD_MPI=
      AS_IF([test "x$FC" = xgfortran],
            [
            AC_MSG_FAILURE([Sorry, parallel option using MPI is not available with gfortran yet. Try using mpifort or mpif90 (configured with gfortran)])
            LD_MPI=" -lmpi"
            #~ FL_CHECK_MPIF_HEADER
            ]
            )
      AS_IF([test "x$FC" = xifort],
            [
            #~ AC_CHECK_LIB([mpi],[main],
                    #~ [AC_CHECK_LIB([mpi_mpifh],[main],[LD_MPI="-lmpi_mpifh -lmpi"],[AC_MSG_FAILURE([-lmpi_mpif flag failed])])]
                    #~ ,[AC_MSG_FAILURE([-lmpi flag failed])])
            LD_MPI="-lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi"
            FL_CHECK_MPIF_HEADER
            ])
      AS_IF([test "x$FC" = xxlf_r],
            [
            FC=mpxlf90_r
            LD_MPI="-lpesslsmp -lesslsmp"
            ])
      AC_SUBST(LD_MPI)
      ])


AS_IF([test "x$with_para" = xno],
      [DMPI="-Dparano"
       EXEEXT=".seq"],
      [test "x$with_para" = xyes],
      [DMPI="-Dparayes"
       EXEEXT=".par"],
      [test "x$with_para" = xsim],
      [DMPI="-Dparano -Dsimpara"
       EXEEXT=".sim"],
      [test "x$with_para" = xworld],
      [DMPI="-Dparano -Dparaworld"
       EXEEXT=".world"]
      )
      
AC_MSG_NOTICE([Using parallelism: $with_para])
AC_MSG_NOTICE([Binary extension will be : "$EXEEXT"])
AS_IF([test  "x$with_para" != xno],
      [AC_MSG_NOTICE([Compiler taking MPI into account: $FC])]
      )
AS_IF([test "x$LD_MPI" != x],
      [AC_MSG_NOTICE([Additional flag for parallelism: LD_MPI=$LD_MPI])]
      )

AM_CONDITIONAL([USE_MPI],[test "x$with_para" != xno])
AC_SUBST(DMPI)
AC_SUBST(EXEEXT)
#=======================================================================
# Choice of OpenMP
#=======================================================================
AC_ARG_WITH([openmp],
  [AS_HELP_STRING([--with-openmp],
  [--with-openmp : Invoke OpenMP. --with-openmp : use threads for FFT.  --with-openmp=dyn : wave function parallelization.])],
  [],
  [with_openmp=no])
AC_MSG_NOTICE([Use OpenMp : $with_openmp])
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
      AC_LANG_PUSH([Fortran])
      AC_OPENMP
      # Because Autotools can't properlly handle library tests :
      AC_LANG_POP([Fortran])
      ])
AC_SUBST(DOMP)
AC_SUBST(OPENMP_FCFLAGS)
#=======================================================================
# Check and set FFT (options and -D flags)
#=======================================================================
# FFT related  option:
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
# CUFFT only availaible when using IFORT (name of variables in cuda.cu : module_mp_variable)
#~ AS_IF([test "x$with_fft" = xCUFFT],
      #~ [AS_IF([test "x$FC" != xifort], 
             #~ [AC_MSG_FAILURE([FFT using CUDA is only available with compiler ifort.]
                              #~ [Sorry.])
              #~ ])
      #~ ])
          


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
      AC_MSG_NOTICE([ --with-wrappers shall be used to specify mkl wrappers path only when MKL is used as FFT library] )
      AC_MSG_NOTICE([ Ignoring option --with-wrappers] ) 
      ])
  ],
  [with_wrappers=yes])
  
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
AC_MSG_NOTICE([wrappers = $with_wrappers])
AC_MSG_NOTICE([mklthreads = $mklthreads])
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
        FL_CHECK_MKL_LIBS
        AC_SUBST([LD_FFTW])
        DFFT="-Dfftw_cpu"
      ],
      #........................FFTW.....................................
      [test "x$with_fft" = xFFTW],
      [    
        FL_CHECK_FFTW_LIBS
        AC_SUBST([LD_FFTW])
        DFFT="-Dfftw_cpu -Dfftwnomkl"
      ],
      #........................CUDA.....................................
      [test "x$with_fft" = xCUFFT],
      [
        FL_LOCATE_CUDA
        AC_MSG_NOTICE([cuda linker flags: LD_CUDA=$LD_CUDA])
        DFFT="-Dfftw_gpu"
        cp code/source_f90/define.h code/source_f90/define_cuda.h
        sed -i -e 's/!/\/\//g' code/source_f90/define_cuda.h #define_cuda.h is just define.h turned into C++
      ],
      #......................NoFFT......................................
      [test "x$with_fft" = xno],
      [
        AC_MSG_NOTICE([No FFT needed])
      ],
      #......................INVALID....................................
      [
        AC_MSG_FAILURE([Invalid FFT option. --with-fft=<choice> .  <choice> can be : NETLIB, FFTW, MKL or CUFFT. If not set, default library is NETLIB. ])
      ]
      )
LIBS=
AC_SUBST(DFFT)
AM_CONDITIONAL([USE_NETLIB],[test "x$with_fft" = xNETLIB])
AM_CONDITIONAL([USE_FFTW],[test "x$FFTW_IS_HERE" = xyes ])
AM_CONDITIONAL([USE_MKL], [test "x$MKL_IS_HERE" = xyes ])
AM_CONDITIONAL([USE_CUDA],[test "x$CUDA_IS_HERE" = xyes ])
AM_CONDITIONAL([CALLS_FFTW],[test "x$FFTW_IS_HERE" = xyes || test "x$MKL_IS_HERE" = xyes ])

AS_IF([test "x$FFTW_IS_HERE" = xyes || test "x$MKL_IS_HERE" = xyes ],
      [AC_MSG_NOTICE([FFTW linker flags: LD_FFTW=$LD_FFTW])])

AC_MSG_NOTICE([preprocessor flags: CPPFLAGS=$CPPFLAGS])
AC_MSG_NOTICE([FFTW preprocessor flags : $DFFT])
AC_MSG_NOTICE([Compiler flags:    FCFLAGS=$FCFLAGS])
AC_MSG_NOTICE([linker flags:      LDFLAGS=$LDFLAGS])



#~ echo "FCFLAGS: $FCFLAGS"
#~ echo "LDFLAGS: $LDFLAGS"
#~ echo "LD_MPI: $LD_MPI" 
#~ echo "LD_FFTW: $LD_FFTW" 
AC_CONFIG_HEADERS([define.h])
AC_CONFIG_FILES([Makefile
                 code/source_f90/Makefile])
AC_OUTPUT
