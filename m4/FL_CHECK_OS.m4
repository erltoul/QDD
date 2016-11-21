# Check type of OS (Linux/MacOS), and bitnes (32/64bits)
#Detect 32 bits of 64 bit config, as well as 
AC_DEFUN([FL_CHECK_OS],
[
AS_IF([test "x$FC" = xmpif90], 
      [COMPILER="$(mpif90 --showme:command)"],
      [COMPILER="$FC"]
      )

AS_IF([test "x$COMPILER" = xifort],
      [
      # check if INTEL compiler uses x86_64 architecture
      ifort -V 2>&1 | grep -i "Intel(R) 64" > /dev/null 2>&1 && flag64=yes
      # check if the platform is OSX
      ifort -dumpmachine 2>&1 | grep -i "darwin" > /dev/null 2>&1 && flagosx=yes
      ],
      [test "x$COMPILER" = xgfortran],
      [
      gfortran -v 2>&1 | grep -i "build=x86_64" >/dev/null 2>&1 && flag64=yes
      # check if the platform is OSX
      gfortran -dumpmachine 2>&1 | grep -i "darwin" > /dev/null 2>&1 && flagosx=yes
      ],
      [test "x$COMPILER" = xxlf_r],
      [
      AC_MSG_ERROR([Not ready for this compiler yet !]) 
      ]    
      )
AS_IF([test "x$flag64" = xyes], [BITS_SIZE=64], [BITS_SIZE=32])
AS_IF([test "x$flagosx" = xyes], [OS_TYPE=MacOSX ])
AC_MSG_NOTICE([Detected $BITS_SIZE bits config using info on compiler $FC])
])
