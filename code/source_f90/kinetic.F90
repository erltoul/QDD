#include"define.h"

#if(gridfft)

#if(fftw_cpu)
INCLUDE "fft.fftw3d.F90"
#endif

!#if(fftw_cpu)
!INCLUDE "fft.fftw1d.F90"
!#endif

!#if(fftw3d_cpu)
!INCLUDE "fft.fftw3d.F90"
!#endif

#if(netlib_fft)
INCLUDE "fft.F90"
#endif

#endif !gridfft

#if(findiff|numerov)
INCLUDE "findiff.F90"
#endif

