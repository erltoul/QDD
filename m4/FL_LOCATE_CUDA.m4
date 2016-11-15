# Autodetection of cuda libraries path.
# If libraries work at first try, do nothing.
# If libraries do not xork at first try, tries and locate -lcudart, and adds to LDFLAGS.
# Works only if nvcc is installed in <path>/bin/nvcc and libs in <path>/lib/ or path <path>/lib64/

AC_DEFUN([FL_LOCATE_CUDA],[
  AC_CHECK_LIB([cudart], [main], [],
                [
                AC_PATH_PROG( [NVCC_PATH], [nvcc ], [missing])
                AS_IF([test "x$NVCC_PATH" = xmissing],
                      [AC_MSG_FAILURE([Could not find nvcc. Please check cuda install.])
                      ])
                #Attention, prevoir une autre solution pour MacOS ?
                NVCC_PATH="$(which nvcc)"
                NVCC_PATH="$(dirname $NVCC_PATH)" 
                NVCC_PATH="$(readlink -f $NVCC_PATH/..)"
                AC_MSG_NOTICE([nvcc has been located in $NVCC_PATH])
                #ajouter choix lib ou lib64
                AS_IF([test "x$BITS_SIZE" = x64],
                      [LDFLAGS="$LDFLAGS -L$NVCC_PATH/lib64"],
                      [LDFLAGS="$LDFLAGS -L$NVCC_PATH/lib"]
                      )
                AC_MSG_NOTICE([Will add $LDFLAGS to linker search path])
                unset ac_cv_lib_cudart_main
                AC_CHECK_LIB([cudart], [main], [],[AC_MSG_FAILURE([could not use -lcudart. Please check cuda install, or add LDFLAGS=-L<path/to/cuda> to ./configure options.])])
                ])
  AC_CHECK_LIB([cufft], [main], [CUDA_IS_HERE=yes],
               [AC_MSG_FAILURE([CUFFT was asked, but test for -lcufft failed, probably because library libcufft.so was not found.])
               ])
])
