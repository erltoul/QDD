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
!   switch to 3D FFTW in Coulomb solver (preliminary option)
#define coudoub3D 0
!
! full SIC and localized SIC
!
#define twostsic 0
#define cmplxsic 0
!
!  
! switch to extended model with polarizable raregas
!
#define raregas 0
!
! switch to old version of 'kinprop' with interlaced 1D FFT
!#define oldkinprop 0

