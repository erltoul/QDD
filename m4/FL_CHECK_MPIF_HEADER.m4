#Check include of mpif.h
AC_DEFUN([FL_CHECK_MPIF_HEADER],
[
LDLIBS="$LD_MPI"
echo "Checking presence of mpif.h for compiler $FC"
AC_LANG_PUSH([Fortran])
AC_MSG_NOTICE([[Testing COMPILATION including 'mpif.h']])
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([],  
                  [[
                  IMPLICIT NONE
                  INCLUDE 'mpif.h']]
                  )
  ],
  [echo "success"],
  [AC_MSG_FAILURE([
    [Could not COMPILE simple fortran program including 'mpif.h'. ]
    [If you really want to use $FC, you should give includes path:]
    [    ./configure  --with-para --with-compiler=$FC  CPPFLAGS=-I/path/to/mpi/include]
    []
    [You should use mpifort instead of $FC.]
    [     ./configure  --with-para --with-compiler=mpifort]
    [Refer to mpifort documentation if you want to change subjacent compiler used by mpifort.]
    []
    [mpif90 and mpiF77 are obsolete versions of mpifort.]
                    ])]

  )
AC_MSG_NOTICE([[Testing LINKING using $LD_MPI]])
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([],  
                  [[
                  IMPLICIT NONE
                  INCLUDE 'mpif.h']]
                  )
  ],
  [echo "success"],
  [AC_MSG_FAILURE([
    [Could not LINK simple fortran program including 'mpif.h'. ]
    [If you really want to use $FC, you should give libraries path:]
    [    ./configure  --with-para --with-compiler=$FC  LDFLAGS=-L/path/to/lib]
    []
    [You should use mpifort instead of $FC.]
    [     ./configure  --with-para --with-compiler=mpifort]
    [Refer to mpifort documentation if you want to change subjacent compiler used by mpifort.]
    []
    [mpif90 and mpiF77 are obsolete versions of mpifort.]
                    ])]
  )
AC_LANG_POP([Fortran])
LDLIBS=
])
