#read variable defined in header
# FL_READ_DEF("header", variable)
# header must be given with cotes " ", and with full path from folder containing configure.ac.
AC_DEFUN([FL_READ_DEF],
[
AC_EGREP_CPP(string_to_match,
  [
#include $1
#if($2)
string_to_match
#endif
  ],
  [value=1],
  [value=0])
  AM_CONDITIONAL([USE_$2], [test "x$value" = x1])
  $2=$value
  AC_SUBST($2)
])

