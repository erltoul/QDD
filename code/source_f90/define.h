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
!choose fft solver (choose only one)
!For the moment, don't use iffastpropag with GPU
!
#define netlib_fft 0
#define fftw_cpu 0
!#define fftw3d_cpu 0
#define fftw_gpu 1
!#define fftw3d_gpu 0

!For the moment compute lda on GPU seems to be unstable and can freeze the computer
#define lda_gpu 0
