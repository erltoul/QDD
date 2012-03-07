#include"define.h"
 
#if(gridfft)
#if(coufou)

#if(fftw_cpu)
INCLUDE 'falr.fftw.F90'
#endif

!#if(fftw3d_cpu)
!INCLUDE 'falr.fftw.F90'
!#endif

#if(fftw_gpu)
INCLUDE 'falr.fftw.F90'
#endif

!#if(fftw3d_gpu)
!INCLUDE 'falr.fftw.F90'
!#endif

#if(netlib_fft)
INCLUDE 'falr.F90'
#endif
#endif !coufou

#if(coudoub)
#if(fftw_cpu)
INCLUDE 'coulex.fftw.F90'
#endif

!#if(fftw3d_cpu)
!INCLUDE 'coulex.fftw3d.F90'
!#endif

#if(fftw_gpu)
INCLUDE 'coulex.fftwgpu3d.F90'
#endif

!#if(fftw_gpu)
!INCLUDE 'coulex.fftwgpu.F90'
!#endif

!#if(fftw3d_gpu)
!INCLUDE 'coulex.fftwgpu3d.F90'
!#endif

#if(netlib_fft)
INCLUDE 'coulex.F90'
#endif

#endif !coudoub
#endif !gridfft

#if(findiff|numerov)
!INCLUDE 'gridcoul.F90'
INCLUDE "findiff-sinft.F90"
#endif
