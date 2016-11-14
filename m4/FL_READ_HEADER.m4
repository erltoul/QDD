#Reads "#define"s from header and set variables accordingly.
AC_DEFUN([FL_READ_HEADER],
[
  FL_READ_DEF(["$1"] , [gridfft])
  FL_READ_DEF(["$1"] , [findiff])
  FL_READ_DEF(["$1"] , [numerov])
  FL_READ_DEF(["$1"] , [coufou])
  FL_READ_DEF(["$1"] , [coudoub])
  FL_READ_DEF(["$1"] , [coudoub3D])
  FL_READ_DEF(["$1"] , [twostsic])
  FL_READ_DEF(["$1"] , [cmplxsic])
  FL_READ_DEF(["$1"] , [raregas])
  FL_READ_DEF(["$1"] , [lda_gpu])
  FL_READ_DEF(["$1"] , [asynclaunch])
  AM_CONDITIONAL( [USE_GRIDFFT],  [test "x$gridfft" = x1 ])
  AM_CONDITIONAL( [USE_FINDIFF],  [test "x$findiff" = x1 ])
  AM_CONDITIONAL( [USE_NUMEROV],  [test "x$numerov" = x1 ])
  AM_CONDITIONAL( [USE_COUFOU],  [test "x$gcoufou" = x1 ])
  AM_CONDITIONAL( [USE_COUDOUB],  [test "x$coudoub" = x1 ])
  AM_CONDITIONAL( [USE_COUDOUB3D],  [test "x$coudoub3D" = x1 ])
  AM_CONDITIONAL( [USE_TWOSTSIC],  [test "x$twostsic" = x1 ])
  AM_CONDITIONAL( [USE_CMPLXSIC],  [test "x$cmplxsic" = x1 ])
  AM_CONDITIONAL( [USE_RAREGAS],  [test "x$graregas" = x1 ])
  AM_CONDITIONAL( [USE_LDA_GPU],  [test "x$lda_gpu" = x1 ])
  AM_CONDITIONAL( [USE_ASYNCLAUNCH],  [test "x$asynclaunch" = x1 ])
]
)

