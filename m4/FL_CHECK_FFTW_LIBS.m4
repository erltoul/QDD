#Make all checks for FFTW library
AC_DEFUN([FL_CHECK_FFTW_LIBS],
[
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
])
