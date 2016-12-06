#Make all checks for MKL library
AC_DEFUN([FL_CHECK_MKL_LIBS],
[
LDLIBS=
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
  CPPFLAGS="$CPPFLAGS  -I$MKL_INCLUDE"
  CFLAGS="$CFLAGS "
  LDFLAGS="$LDFLAGS -L$MKL_WRAPPERS -L$MKL_LIBPATH -Wl,-rpath,$MKL_WRAPPERS -Wl,-rpath,$MKL_LIBPATH"
  
  #Check common libraries
  AC_CHECK_LIB([dl], [main],[],[MKL_IS_HERE=no],[])
  AC_CHECK_LIB([m], [main],[],[MKL_IS_HERE=no])
  AC_CHECK_LIB([pthread], [main],[],[MKL_IS_HERE=no])
  # Add iomp5 library if mkl threads are used
  AS_IF([test "x$mklthreads" = xyes],
        [AS_IF([test "x$flag64" = xyes],
              [OMP5_LIBPATH="$mklroot/../compiler/lib/intel64"],
              [OMP5_LIBPATH="$mklroot/../compiler/lib/ia32"])
        LDFLAGS="-L$OMP5_LIBPATH  $LDFLAGS -Wl,-rpath,$OMP5_LIBPATH"
        AC_CHECK_LIB([iomp5], [main], 
                      [MKLFLAGS="-liomp5"
                      LIBS="-liomp5 $LIBS"],
                      [AC_MSG_FAILURE([Library libiomp5 not found])]
                    )
        ],
        [MKLFLAGS=])
  #check mkl libraries
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
  AS_IF([test "x$MKL_IS_HERE" = xno], 
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
  LD_FFTW="-lfftw3xf_intel $MKLFLAGS -lpthread -lm"
         #~ ],
        #~ [LD_FFTW="$MKLFLAGS -liomp5 -lpthread -lm -ldl"])
])
