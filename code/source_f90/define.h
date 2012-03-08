!
!       code version:
!
#define IVERSION 73
!
! grid: FFT or finite difference
!
#define gridfft 1
#define findiff 0
#define numerov 0
!
! coulomb solver  (set all to 0 for 'findiff' or 'numerov')
!
#define coufou 0
#define coudoub 1
!
! parallel or serial:
!
#define parayes 0
#define parano 1
#define simpara 0
!
! full SIC and variants (outdated)
!
#define fullsic 0
#define symmcond 0
#define twostsic 0
!
!  
! switch to extended model with polarizable raregas
!
#define raregas 0
!
!choose fft solver (only one !)
!
#define netlib_fft 0
#define fftw_cpu 1
!#define fftw3d_cpu 0 !Used only to determine wich 1d or 3d FFTS are the fastest

