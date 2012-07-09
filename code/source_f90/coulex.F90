MODULE coulsolv
#if(fftw_cpu|fftw_gpu)
USE, intrinsic :: iso_c_binding
#endif
USE params, ONLY: DP
#if(fftw_gpu)
USE cuda_alloc
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

SAVE
INTEGER,PRIVATE :: kxmax,kymax,kzmax,ksmax
! kxmax must be the largest
INTEGER,PRIVATE :: kdfull
INTEGER,PRIVATE :: kdred
INTEGER,PRIVATE :: kfft2
!INTEGER,PARAMETER,PRIVATE :: kddoub=kdfull
INTEGER,PRIVATE :: kfft,kfftx,kffty,kfftz
!INTEGER,PARAMETER,PRIVATE :: kdcorf=(kxmax/2+1)*(kymax/2+1)*(kzmax/2+1)
! include block: xkgrid
REAL(DP),ALLOCATABLE,PRIVATE :: xval(:),yval(:),zval(:)
REAL(DP),ALLOCATABLE,PRIVATE :: xt2(:),yt2(:),zt2(:)
REAL(DP),PRIVATE :: dx,dy,dz,dxsp,grnorm,fnorm
INTEGER,PRIVATE :: nx,ny,nz,nx1,ny1,nz1, nxr,nxi,nyr,nyi,nzr,nzi,nxy1,nxyz
INTEGER,PRIVATE :: nxhigh,nxlow,nyhigh,nylow,nzhigh,nzlow
!COMMON /xgrid/ xval,yval,zval,xt2,yt2,zt2,dx,dy,dz,dxsp,grnorm,fnorm,  &
!    nx,ny,nz,nx1,ny1,nz1, nxr,nxi,nyr,nyi,nzr,nzi,nxy1,nxyz,  &
!    nxhigh,nxlow,nyhigh,nylow,nzhigh,nzlow

INTEGER,ALLOCATABLE,PRIVATE :: ikm(:,:)
REAL(DP),PRIVATE :: dkx,dky,dkz,akmax,dksp,ecut
INTEGER,PRIVATE :: nxk,nxklo,nxkhi,nksp,nkxyz
!COMMON /kgrid/ akv2r,akv2i, dkx,dky,dkz,akmax,dksp,  &
!    ikm,nxk,nxklo,nxkhi,nksp,nkxyz,ecut

!COMMON /fftini/wrkx,wrky,wrkz,wsavex,wsavey,wsavez,ifacx,ifacy,ifacz

! include block: option
REAL(DP),PARAMETER,PRIVATE :: zero=0D0
REAL(DP),PARAMETER,PRIVATE :: pi=3.141592653589793D0

!COMMON /fftcom/fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax)

#if(netlib_fft)
REAL(DP),ALLOCATABLE,PRIVATE :: akv2r(:),akv2i(:)
REAL(DP),ALLOCATABLE,PRIVATE :: wrkx(:),wrky(:),wrkz(:)
REAL(DP),ALLOCATABLE,PRIVATE :: wsavex(:),wsavey(:),wsavez(:)
INTEGER,ALLOCATABLE,PRIVATE :: ifacx(:),ifacy(:),ifacz(:)
COMPLEX(DP),ALLOCATABLE,PRIVATE :: fftax(:),fftay(:),fftb(:,:)
#endif

#if(fftw_cpu)
REAL(DP),ALLOCATABLE,PRIVATE :: akv2r(:),akv2i(:)
type(C_PTR), PRIVATE :: pforwx,pforwy,pforwz,pbackx,pbacky,pbackz
INTEGER(C_INT), PRIVATE :: wisdomtest
COMPLEX(C_DOUBLE_COMPLEX),ALLOCATABLE,PRIVATE :: fftax(:),fftay(:),fftb(:,:)
#endif

#if(fftw_gpu)
REAL(C_DOUBLE),ALLOCATABLE,PRIVATE :: akv2r(:),akv2i(:)
INTEGER*8, PRIVATE :: pfft
INTEGER,PARAMETER,PRIVATE :: batch=1
COMPLEX(C_DOUBLE_COMPLEX),PRIVATE, POINTER :: fftac(:,:,:),akvc(:,:,:)
TYPE(C_PTR), PRIVATE :: c_p_fftac,c_p_akvc
COMPLEX(C_DOUBLE_COMPLEX),PRIVATE,POINTER :: gpu_fftac(:,:,:),gpu_akvc(:,:,:)
TYPE(C_PTR),PRIVATE :: c_gpu_fftac,c_gpu_akvc
INTEGER, PRIVATE :: res
COMPLEX(C_DOUBLE_COMPLEX), PARAMETER :: size_cmplx=CMPLX(0D0,0D0)
#endif

CONTAINS

SUBROUTINE init_coul(dx0,dy0,dz0,nx0,ny0,nz0)

!-----------------------------------------------------------------------


!     read grid parameters from file or simply initialize them
!     note that the Coulomb solver doubles the grid internally
kxmax=2*nx0;kymax=2*ny0;kzmax=2*nz0;ksmax=kxmax
kdfull=nx0*ny0*nz0
kdred=kxmax*kymax*kzmax
kfft=ksmax;kfftx=kxmax;kffty=kymax;kfftz=kzmax
kfft2=kfft*2

nx=nx0  !/2
ny=ny0  !/2
nz=nz0  !/2
dx=dx0
dy=dy0
dz=dz0

ALLOCATE(xval(kxmax),yval(kymax),zval(kzmax))
ALLOCATE(xt2(kxmax),yt2(kymax),zt2(kzmax))
ALLOCATE(akv2r(kdred),akv2i(kdred))
ALLOCATE(ikm(kxmax,kymax))
#if(netlib_fft)
ALLOCATE(fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax))
ALLOCATE(wrkx(kfft2),wrky(kfft2),wrkz(kfft2))
ALLOCATE(wsavex(kfft2),wsavey(kfft2),wsavez(kfft2))
ALLOCATE(ifacx(kfft2),ifacy(kfft2),ifacz(kfft2))
#endif
#if(fftw_cpu)
ALLOCATE(fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax))
#endif
#if(fftw_gpu)
CALL my_cuda_allocate_coulex() ! Pinned memory allocation to make CPU>GPU and GPU>CPU transfers faster
#endif

!     call input routine fftinp, which initializes the grid and fft tabl

CALL fftinp

RETURN
END SUBROUTINE init_coul

!-----fftinp------------------------------------------------------------

SUBROUTINE fftinp
IMPLICIT REAL(DP) (A-H,O-Z)

!     initializes work tables for FFT

!     grid parameters nx,ny,nz,dx,dy,dz,ecut must have been read or
!     initialized before !

!-----------------------------------------------------------------------

!     initialize grid in coordinate space

nx1=nx+1
ny1=ny+1
nz1=nz+1
nxr=nx+nx
nyr=ny+ny
nzr=nz+nz
nxi=nxr
nyi=nyr
nzi=nzr
nxy1=nxi*nyi
nxyz=nxi*nyi*nzi
nkxyz=nxi*nyi*nzi

!     grid lengths must match with parameters in incs

IF(kxmax < nxi) THEN
  WRITE(6,'(a)') ' ERROR: parameter   kxmax   too small'
  STOP ' error in parameter: KXMAX in COULEX too small'
ELSE IF(kymax < nyi) THEN
  WRITE(6,'(a)') ' ERROR: parameter   kymax   too small'
  STOP ' error in parameter: KYMAX in COULEX too small'
ELSE IF(kzmax < nzi) THEN
  WRITE(6,'(a)') ' ERROR: parameter   kzmax   too small'
  STOP ' error in parameter: KZMAX in COULEX too small'
END IF

!     initialize grid in fourier space

dkx=pi/(dx*REAL(nx))
dky=pi/(dy*REAL(ny))
dkz=pi/(dz*REAL(nz))

dxsp=dx*dy*dz
dksp=dkx*dky*dkz

grnorm=SQRT(dxsp/dksp)
fnorm=1.0/SQRT(REAL(nx*ny*nz))

!test      akmax=sqrt(3*(nx*nx)*dx*dx)+2.0
!test      nxk=int(akmax/dkx)+1
!test      if(nxk.gt.nx1) nxk=nx1
nxk=nx1

!     built Greens function in Fourier space
!     by Fourier transformation from real space

ikzero = nxy1*(nz-1)+nxi*(ny-1)+nx
write(*,*) ' nzi,nyi,nxi,nx,ny,nz,ikzero=',nzi,nyi,nxi,nx,ny,nz,ikzero
ii=0
xz1=-nz*dz
DO i3=1,nzi
  xz1=xz1+dz
  xz2=xz1*xz1
  xy1=-ny*dy
  DO i2=1,nyi
    xy1=xy1+dy
    xy2=xy1*xy1
    xx1=-nx*dx
    DO i1=1,nxi
      xx1=xx1+dx
      xx2=xx1*xx1
      ak2=xx2+xy2+xz2
      ii=ii+1
!        write(*,*) ' i1,i2,i3,ii=',i1,i2,i3,ii
      IF(ii /= ikzero) THEN
        akv2r(ii) =  1D0/SQRT(ak2)
      ELSE
!              akv2r(ii) = (6D0*pi/(dx*dy*dz))**(1D0/3D0)  ! spherical approx
!              akv2r(ii) = 1.19003868*(dx*dy*dz)**(-1D0/3D0)
        akv2r(ii) = 2.34D0*1.19003868D0*(dx*dy*dz)**(-1D0/3D0)  ! empirical
      END IF
      akv2i(ii) = 0D0
    END DO
  END DO
END DO
nksp=ii

#if(netlib_fft|fftw_cpu)
CALL fourf(akv2r(1),akv2i(1))

!  CALL prifld(akv2r,'k**2 for Co')
!  CALL prifld(akv2i,'k**2 for Co')
!nxyf=nx2*ny2
!nyf=nx2
ind=0 !nz*nxyf+ny*nyf
WRITE(6,*) 'k**2 along x'
DO i1=1,nxi
  WRITE(6,'(f8.2,2(1pg13.5))') (i1-nx)*dx,akv2r(i1+ind),akv2i(i1+ind)
END DO
#endif

#if(fftw_gpu)
CALL fourfakv(akv2r(1),akv2i(1))

res = cudaFreeHost(c_p_akvc)
DEALLOCATE(akv2r,akv2i)       !Only gpu_akvc will be used from now
#endif

RETURN
END SUBROUTINE fftinp

!-------------------------------------------------------------------

SUBROUTINE falr(rhoinp,chpfalr,nxdum,nydum,nzdum,kdum)

IMPLICIT REAL(DP) (A-H,O-Z)


REAL(DP), INTENT(IN)                     :: rhoinp(kdfull)
REAL(DP), INTENT(OUT)                    :: chpfalr(kdfull)
INTEGER, INTENT(IN)                  :: nxdum
INTEGER, INTENT(IN)                  :: nydum
INTEGER, INTENT(IN)                  :: nzdum
INTEGER, INTENT(IN)                  :: kdum

REAL(DP),ALLOCATABLE :: rhokr(:),rhoki(:)

ALLOCATE(rhokr(kdred),rhoki(kdred))

!     call a routine written by you which writes your density field
!     on the array rho.
!     remember not to send your original density array to the fcs.
!     in this case we have a homogeneously charged sphere .


CALL rhofld(rhoinp,rhokr,rhoki)


!     call coufou, which contains the fcs procedure.

CALL coufou2(rhokr,rhoki)

!     call a routine written by you which outputs the results of the fcs
!     and maybe some other things to an output file or the screen.

CALL result(chpfalr,rhokr,rhoki)

DEALLOCATE(rhokr,rhoki)



END SUBROUTINE falr


!-----rhofld------------------------------------------------------------

SUBROUTINE rhofld(rhoinp,rhokr,rhoki)
IMPLICIT REAL(DP) (A-H,O-Z)

!     copy density on complex array of double extnesion in x,y,z


REAL(DP), INTENT(IN)                         :: rhoinp(kdfull)
REAL(DP), INTENT(OUT)                        :: rhokr(kdred)
REAL(DP), INTENT(OUT)                        :: rhoki(kdred)


rhokr=0D0
rhoki=0D0
i0=0
DO i3=1,nz
  DO i2=1,ny
    ii = (i3-1)*nxi*nyi+(i2-1)*nxi
    DO i1=1,nx
      ii=ii+1
      i0 = i0+1
      rhokr(ii)=rhoinp(i0)
    END DO
  END DO
END DO


RETURN
END SUBROUTINE rhofld


!-----result------------------------------------------------------------

SUBROUTINE result(chpfalr,rhokr,rhoki)

IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(OUT)                        :: chpfalr(kdfull)
REAL(DP), INTENT(IN)                         :: rhokr(kdred)
REAL(DP), INTENT(IN OUT)                     :: rhoki(kdred)



!     copy Coulomb field back to standard grid

ii=0
i0=0
DO i3=1,nzi
  DO i2=1,nyi
    DO i1=1,nxi
      ii=ii+1
      IF(i3 <= nz .AND. i2 <= ny .AND. i1 <= nx) THEN
        i0 = i0+1
        chpfalr(i0) = 2D0*rhokr(ii)
      END IF
    END DO
  END DO
END DO

RETURN
END SUBROUTINE result


!-----cofows------------------------------------------------------------

SUBROUTINE coufou2(rhokr,rhoki)

IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: rhokr(kdred)
REAL(DP), INTENT(IN OUT)                     :: rhoki(kdred)


LOGICAL,PARAMETER :: tprint=.false.
LOGICAL,PARAMETER :: rqplot=.false.

!------------------------------------------------------------------------------

!     fourier transformation of the density

CALL fourf(rhokr,rhoki)

!     calculation of the coulomb field (writing on the density field)
#if(netlib_fft|fftw_cpu)
DO ik=1,kdred
  SAVE2     = akv2r(ik)*rhokr(ik)+akv2i(ik)*rhoki(ik)
  rhoki(ik) = akv2r(ik)*rhoki(ik)+akv2i(ik)*rhokr(ik)
  rhokr(ik) = SAVE2
END DO
#endif
#if(fftw_gpu)
CALL multiply_ak(gpu_fftac,gpu_akvc,kdred)
#endif
!     fourier back transformation

CALL fourb(rhokr,rhoki)

RETURN

END SUBROUTINE coufou2

#if(fftw_gpu)
!-----fourfakv---------------------------------------------------------

SUBROUTINE fourfakv(pskr,pski)

IMPLICIT REAL(DP) (A-H,O-Z)


REAL(DP), INTENT(OUT)                        :: pskr(kdred)
REAL(DP), INTENT(OUT)                        :: pski(kdred)
INTEGER :: typefft=3

!     fourier forward transformation
!     I/O: pskr   real part of the wave-function
!          pski   imaginary part of the wave-function

DATA  mini/0/              ! flag for initialization
!----------------------------------------------------------------------

!     check initialization

IF(mini == 0) THEN
  CALL cuda_plan_3d(pfft,kxmax,kymax,kzmax)
  mini  = kxmax*kymax*kzmax
  WRITE(7,'(a)') ' fft initialized '
ELSE IF(mini /= kxmax*kymax*kzmax) THEN
  STOP ' n2 in four3d not as initialized!'
END IF

nzzh=(nzi-1)*nxy1
nyyh=(nyi-1)*nxi
!test      sqh=sqrt(0.5)

tnorm=grnorm*fnorm

CALL copyr1dto3d(pskr,pski,akvc,kfftx,kffty,kfftz)

CALL copy_on_gpu(akvc,gpu_akvc,kdred)

CALL run_fft_for3d(pfft,gpu_akvc,typefft)

CALL multiply_gpu(gpu_akvc,kdred,tnorm)

!test2      write(6,'(a)') ' TNORM=',tnor

CALL cpu_time(time_fin)
!WRITE(*,*) 'FOURF end: time=',time_fin-time_ini

RETURN
END SUBROUTINE fourfakv
#endif

!-----fourf-------------------------------------------------------fourf

SUBROUTINE fourf(pskr,pski)
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT REAL(DP) (A-H,O-Z)


REAL(DP), INTENT(OUT)                        :: pskr(kdred)
REAL(DP), INTENT(OUT)                        :: pski(kdred)


!     fourier forward transformation
!     I/O: pskr   real part of the wave-function
!          pski   imaginary part of the wave-function

#if(netlib_fft|fftw_cpu)
DATA  mxini,myini,mzini/0,0,0/              ! flag for initialization
#endif
#if(fftw_gpu)
DATA  mini/0/              ! flag for initialization
#endif
!----------------------------------------------------------------------


!     check initialization
#if(netlib_fft)
IF(mxini == 0) THEN
  CALL dcfti1 (kfftx,wsavex,ifacx)
  mxini  = kfftx
  WRITE(7,'(a)') ' x-fft initialized '
ELSE IF(mxini /= kfftx) THEN
  STOP ' nx2 in four3d not as initialized!'
END IF
IF(myini == 0) THEN
  CALL dcfti1 (kffty,wsavey,ifacy)
  myini  = kffty
  WRITE(7,'(a)') ' y-fft initialized '
ELSE IF(myini /= kffty) THEN
  STOP ' ny2 in four3d not as initialized!'
END IF
IF(mzini == 0) THEN
  CALL dcfti1 (kfftz,wsavez,ifacz)              !! here was
  mzini  = kfftz                              !! the bug
  WRITE(7,'(a)') ' z-fft initialized '
ELSE IF(mzini /= kfftz) THEN
  STOP ' nz2 in four3d not as initialized!'
END IF
#endif

#if(fftw_cpu)
IF(mxini == 0) THEN
  wisdomtest=fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
  IF (wisdomtest == 0) THEN
    WRITE(6,*) 'wisdom_fftw.dat not found, creating it'
    WRITE(7,*) 'wisdom_fftw.dat not found, creating it'
  END IF
  pforwx=fftw_plan_dft_1d(kfftx,fftax,fftax,FFTW_FORWARD,FFTW_EXHAUSTIVE)
  pbackx=fftw_plan_dft_1d(kfftx,fftax,fftax,FFTW_BACKWARD,FFTW_EXHAUSTIVE)
  mxini  = kfftx
  WRITE(7,'(a)') ' x-fft initialized '
ELSE IF(mxini /= kfftx) THEN
  STOP ' nx2 in four3d not as initialized!'
END IF
IF(myini == 0) THEN
  pforwy=fftw_plan_dft_1d(kffty,fftay,fftay,FFTW_FORWARD,FFTW_EXHAUSTIVE)
  pbacky=fftw_plan_dft_1d(kffty,fftay,fftay,FFTW_BACKWARD,FFTW_EXHAUSTIVE)
  myini  = kffty
  WRITE(7,'(a)') ' y-fft initialized '
ELSE IF(myini /= kffty) THEN
  STOP ' ny2 in four3d not as initialized!'
END IF
IF(mzini == 0) THEN
  pforwz=fftw_plan_dft_1d(kfftz,fftb,fftb,FFTW_FORWARD,FFTW_EXHAUSTIVE)
  pbackz=fftw_plan_dft_1d(kfftz,fftb,fftb,FFTW_BACKWARD,FFTW_EXHAUSTIVE)
  mzini  = kfftz
  wisdomtest=fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
  IF (wisdomtest == 0) THEN
    WRITE(6,*) 'Error exporting wisdom to file wisdom_fftw.dat'
    WRITE(7,*) 'Error exporting wisdom to file wisdom_fftw.dat'
  END IF
  WRITE(7,'(a)') ' z-fft initialized '
ELSE IF(mzini /= kfftz) THEN
  STOP ' nz2 in four3d not as initialized!'
END IF
#endif

nzzh=(nzi-1)*nxy1
nyyh=(nyi-1)*nxi
!test      sqh=sqrt(0.5)

#if(netlib_fft|fftw_cpu)
tnorm=grnorm*fnorm

CALL fftx(pskr,pski)
CALL ffty(pskr,pski)
CALL fftz(pskr,pski)

DO i1=1,nkxyz
  pskr(i1)=tnorm*pskr(i1)
  pski(i1)=tnorm*pski(i1)
END DO
#endif
#if(fftw_gpu)
CALL fft(pskr,pski)
#endif

!test2      write(6,'(a)') ' TNORM=',tnor

CALL cpu_time(time_fin)
!WRITE(*,*) 'FOURF end: time=',time_fin-time_ini

RETURN
END SUBROUTINE fourf

!-----fourb------------------------------------------------------------

SUBROUTINE fourb(pskr,pski)
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(OUT)                        :: pskr(kdred)
REAL(DP), INTENT(OUT)                        :: pski(kdred)


!     fourier backward transformation
!     I/O:  pskr   real part of the wave-function
!           pski   imaginary part of the wave-function
!----------------------------------------------------------------------

nzzh=(nzi-1)*nxy1
nyyh=(nyi-1)*nx1
!      sq2=sqrt(2.0)

#if(netlib_fft|fftw_cpu)
tnorm=fnorm/(8D0*grnorm)*pi**1.5D0
!      tnorm=fnorm/(16D0*grnorm)

DO i=1,nkxyz
  pskr(i)=tnorm*pskr(i)
  pski(i)=tnorm*pski(i)
END DO

CALL ffbz(pskr,pski)
CALL ffby(pskr,pski)
CALL ffbx(pskr,pski)
#endif
#if(fftw_gpu)
CALL ffb(pskr,pski)
#endif

RETURN
END SUBROUTINE fourb

#if(netlib_fft|fftw_cpu)
!-----fftx-------------------------------------------------------------

SUBROUTINE fftx(psxr,psxi)
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: psxr(kdred)
REAL(DP), INTENT(IN OUT)                        :: psxi(kdred)

!     performs the fourier-transformation in x-direction
!     the input-wave-function (psxr,psxi) (i.e. real and imaginary part)
!     is overwritten by the fourier-transformed wave-function

!----------------------------------------------------------------------

nx11=nx1+1

i30=-nxy1

DO i3=1,nzi
  i30=i30+nxy1
  i0=i30-nxi
  
  DO i2=1,nyi
    i0=i0+nxi
    
!         composition of the wave-function
!         positive space
    ii=i0+nx-1
    DO i1=1,nx1
      ii=ii+1
      fftax(i1) = psxr(ii)
    END DO
    
!         negative space
    ii=i0
    DO i1=nx11,nxi
      ii=ii+1
      fftax(i1) = psxr(ii)
    END DO
    
!         execution of the fourier-transformation

#if(netlib_fft)
    CALL dcftf1 (kfftx,fftax,wrkx,wsavex,ifacx)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pforwx,fftax,fftax)
#endif
!         decomposition of the wave-function
!        positive space
    
    ii=i0+nx-1
    DO i1=1,nx1
      ii=ii+1
      psxr(ii) = REAL(fftax(i1))
      psxi(ii) = AIMAG(fftax(i1))
    END DO
!        negative space
    ii=i0
    DO  i1=nx11,nxi
      ii=ii+1
      psxr(ii) = REAL(fftax(i1))
      psxi(ii) = AIMAG(fftax(i1))
    END DO
    
  END DO
END DO


RETURN
END SUBROUTINE fftx


!-----ffty--------------------------------------------------------------

SUBROUTINE ffty(psxr,psxi)
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                        :: psxr(kdred)
REAL(DP), INTENT(IN OUT)                     :: psxi(kdred)

!     performs the fourier-transformation in y-direction
!     the input-wave-function (psxr,psxi) (i.e. real and imaginary part)
!     is overwritten by the fourier-transformed wave-function
!     yp is the parity in y-direction (input!)

!----------------------------------------------------------------------

ny11=ny1+1

IF(nxk < nx1) THEN
  nxklo=nx-nxk+1
ELSE
  nxklo=1
END IF
nxkhi=nx+nxk-1

i30=-nxy1
DO i3=1,nzi
  i30=i30+nxy1
!test        do 10 i1=nxklo,nxkhi
  DO i1=nx,nxi ! 1,nx1
    i0=i1+i30
    
!         composition of the wave-function
!         positive space
    ii=i0+nxi*(ny-2)
    DO i2=1,ny1
      ii=ii+nxi
      fftay(i2)=CMPLX(psxr(ii),psxi(ii),DP)
    END DO
!         negative space
    ii=i0-nxi
    DO i2=ny11,nyi
      ii=ii+nxi
      fftay(i2)=CMPLX(psxr(ii),psxi(ii),DP)
    END DO
    
!         execution of the fourier-transformation

#if(netlib_fft)
    CALL dcftf1 (kffty,fftay,wrky,wsavey,ifacy)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pforwy,fftay,fftay)
#endif

!         decomposition of the wave-function
!         positive space
    ii=i0+nxi*(ny-2)
    DO i2=1,ny1
      ii=ii+nxi
      psxr(ii)=REAL(fftay(i2))
      psxi(ii)=AIMAG(fftay(i2))
    END DO
!         negative space
    ii=i0-nxi
    DO i2=ny11,nyi
      ii=ii+nxi
      psxr(ii)=REAL(fftay(i2))
      psxi(ii)=AIMAG(fftay(i2))
    END DO
    
    
  END DO
END DO


RETURN
END SUBROUTINE ffty
!-----fftz-------------------------------------------------------------

SUBROUTINE fftz(psxr,psxi)
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                        :: psxr(kdred)
REAL(DP), INTENT(IN OUT)                     :: psxi(kdred)
INTEGER :: nxyf 
INTEGER :: nyf  
INTEGER :: nzh  

!     performs the fourier-transformation in z-direction
!     the input-wave-function (psxr,psxi) (i.e. real and imaginary part)
!     is overwritten by the fourier-transformed wave-function

!----------------------------------------------------------------------

nxyf = kfftx*kffty
nyf  = kfftx
nzh  = kfftz/2

DO i2=1,kffty
  DO i3=1,kfftz
    i3m   = MOD(i3+nzh,kfftz)+1
!test          do i1=1,kfftx
    DO i1=nx,nxi ! 1,nx1
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftb(i3m,i1) = CMPLX(psxr(ind),psxi(ind),DP)
    END DO
  END DO
!test        do i1=1,kfftx
  DO i1=nx,nxi ! 1,nx1
#if(netlib_fft)
    CALL dcftf1 (kfftz,fftb(1,i1),wrkz,wsavez,ifacz)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pforwz,fftb(1,i1),fftb(1,i1))
#endif
  END DO
  DO i3=1,kfftz
    i3m   = MOD(i3+nzh,kfftz)+1
!test          do i1=1,kfftx
    DO i1=nx,nxi ! 1,nx1
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      psxr(ind)= REAL(fftb(i3m,i1))
      psxi(ind)= AIMAG(fftb(i3m,i1))
    END DO
  END DO
!        i1 = kfftx/2
!        if(i2.eq.kffty/2)
!     &   write(*,'(a,50(/1x,2(2(1pg12.4,2x))))')
!     &         'PSXR:',(psxr((i3-1)*nxyf+(i2-1)*nyf+i1),i3=1,kfftz)
END DO

RETURN
END SUBROUTINE fftz
!-----ffbz-------------------------------------------------------------

SUBROUTINE ffbz(psxr,psxi)
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                        :: psxr(kdred)
REAL(DP), INTENT(IN OUT)                     :: psxi(kdred)
INTEGER :: nxyf 
INTEGER :: nyf  
INTEGER :: nzh  

!----------------------------------------------------------------------

nxyf = kfftx*kffty
nyf  = kfftx
nzh  = kfftz/2

DO i2=1,kffty
  DO i3=1,kfftz
    i3m   = MOD(i3+nzh,kfftz)+1
    DO i1=nx,nxi ! 1,nx1
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftb(i3m,i1) = CMPLX(psxr(ind),psxi(ind),DP)
    END DO
  END DO
  DO i1=nx,nxi ! 1,nx1
#if(netlib_fft)
    CALL dcftb1 (kfftz,fftb(1,i1),wrkz,wsavez,ifacz)    ! basic fft
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pbackz,fftb(1,i1),fftb(1,i1))
#endif
  END DO
  DO i3=1,kfftz                  ! copy back
    i3m   = MOD(i3+nzh,kfftz)+1
    DO i1=nx,nxi ! 1,nx1
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      psxr(ind) = REAL(fftb(i3m,i1))
      psxi(ind) = AIMAG(fftb(i3m,i1))
    END DO
  END DO
END DO

RETURN
END SUBROUTINE ffbz
!-----ffby-------------------------------------------------------------

SUBROUTINE ffby(psxr,psxi)
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                        :: psxr(kdred)
REAL(DP), INTENT(IN OUT)                     :: psxi(kdred)

!----------------------------------------------------------------------

ny11=ny1+1

IF(nxk < nx1) THEN
  nxklo=nx-nxk+1
ELSE
  nxklo=1
END IF
nxkhi=nx+nxk-1

i30=-nxy1

!DO i3=1,nzi
DO i3=1,nz1
  i30=i30+nxy1
!test        do 10 i1=nxklo,nxkhi
  DO i1=nx,nxi ! 1,nx1
    i0=i1+i30
    
!         composition of the complex wave-function
!         positive space
    ii=i0+nxi*(ny-2)
    DO i2=1,ny1
      ii=ii+nxi
      fftay(i2)=CMPLX(psxr(ii),psxi(ii),DP)
    END DO
!         negative space
    ii=i0-nxi
    DO i2=ny11,nyi
      ii=ii+nxi
      fftay(i2)=CMPLX(psxr(ii),psxi(ii),DP)
    END DO
    
!         execution
#if(netlib_fft)
    CALL dcftb1 (kffty,fftay,wrky,wsavey,ifacy)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pbacky,fftay,fftay)
#endif
!         decomposition of the inverse transformed wave-function
!         positive space
    ii=i0+nxi*(ny-2)
    DO i2=1,ny1
      ii=ii+nxi
      psxr(ii)=REAL(fftay(i2))
      psxi(ii)=AIMAG(fftay(i2))
    END DO
!         negative space
    ii=i0-nxi
    DO i2=ny11,nyi
      ii=ii+nxi
      psxr(ii)=REAL(fftay(i2))
      psxi(ii)=AIMAG(fftay(i2))
    END DO 
  END DO
END DO


RETURN
END SUBROUTINE ffby
!-----ffbx-------------------------------------------------------------

SUBROUTINE ffbx(psxr,psxi)
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                        :: psxr(kdred)
REAL(DP), INTENT(IN OUT)                     :: psxi(kdred)

!----------------------------------------------------------------------

nx11=nx1+1

i30=-nxy1

!DO i3=1,nzi
DO i3=1,nz1
  i30=i30+nxy1
  i0=i30-nxi
!  DO i2=1,nyi
  DO i2=1,ny1
    i0=i0+nxi
    
!         composition of the complex wave-function
!        positive space
    ii=i0+nx-1
    DO i1=1,nx1
      ii=ii+1
      fftax(i1)=CMPLX(psxr(ii),psxi(ii),DP)
    END DO
!        negative space
    DO i1=nx11,nxi
      fftax(i1)=CONJG(fftax(nxi-i1+2))
    END DO
    
!        execution
#if(netlib_fft)
    CALL dcftb1 (kfftx,fftax,wrkx,wsavex,ifacx)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pbackx,fftax,fftax)
#endif
!         decomposition of the inverse transformed wave-function
!         positive space
    ii=i0+nx-1
    DO i1=1,nx1
      ii=ii+1
      psxr(ii)=REAL(fftax(i1))
      psxi(ii)=AIMAG(fftax(i1))
    END DO
!         negative space
    ii=i0
    DO i1=nx11,nxi
      ii=ii+1
      psxr(ii)=REAL(fftax(i1))
      psxi(ii)=AIMAG(fftax(i1))
    END DO
  END DO  
END DO
  
RETURN
END SUBROUTINE ffbx
#endif

#if(fftw_gpu)
!-----fft--------------------------------------------------------------

SUBROUTINE fft(psxr,psxi)

USE, intrinsic :: iso_c_binding
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: psxr(kdred)
REAL(DP), INTENT(IN OUT)                     :: psxi(kdred)
INTEGER :: typefft=3

tnorm=grnorm*fnorm

CALL copyr1dto3d(psxr,psxi,fftac,kfftx,kffty,kfftz)

CALL copy_on_gpu(fftac,gpu_fftac,kdred)

CALL run_fft_for3d(pfft,gpu_fftac,typefft)

CALL multiply_gpu(gpu_fftac,kdred,tnorm)

RETURN
END SUBROUTINE fft

!-----ffb--------------------------------------------------------------

SUBROUTINE ffb(psxr,psxi)

USE, intrinsic :: iso_c_binding
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: psxr(kdred)
REAL(DP), INTENT(IN OUT)                     :: psxi(kdred)
INTEGER :: typefft=3

!----------------------------------------------------------------------

tnorm=fnorm/(8D0*grnorm)*pi**1.5D0

CALL run_fft_back3d(pfft,gpu_fftac,typefft)

CALL multiply_gpu(gpu_fftac,kdred,tnorm)

CALL copy_from_gpu(fftac,gpu_fftac,kdred)

CALL copyr3dto1d(fftac,psxr,psxi,kfftx,kffty,kfftz)

RETURN
END SUBROUTINE ffb

! ******************************

SUBROUTINE copyr1dto3d(rpart,impart,vec3d,nbx2,nby2,nbz2)

! ******************************

USE params
USE, intrinsic :: iso_c_binding

REAL(DP), INTENT(IN)                      :: rpart(kdred),impart(kdred)
COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT)    :: vec3d(nbx2,nby2,nbz2)
INTEGER :: nxyfn,nyfn,nnx,nny,nnz

nxyfn = nbx2*nbz2
nyfn  = kfftx
nnx2=nx+nx
nny2=ny+ny
nnz2=nz+nz

DO i3=1,nbz2
  DO i2=1,nby2
    DO i1=1,nbx2
      ind=(i3-1)*nxyfn+(i2-1)*nyfn+i1
      vec3d(MOD(i1+nnx2,nbx2)+1,MOD(i2+nny2,nby2)+1,MOD(i3+nnz2,nbz2)+1)=CMPLX(rpart(ind),impart(ind),DP)
    END DO
  END DO
END DO

RETURN
END SUBROUTINE copyr1dto3d

! ******************************

SUBROUTINE copyr3dto1d(vec3d,rpart,impart,nbx2,nby2,nbz2)

! ******************************

USE params
USE, intrinsic :: iso_c_binding

COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN)     :: vec3d(nbx2,nby2,nbz2)
REAL(DP), INTENT(OUT)                     :: rpart(kdred),impart(kdred)
INTEGER :: nxyfn,nyfn,nnx,nny,nnz

nxyfn = nbx2*nbz2
nyfn  = kfftx
nnx2=nx+nx
nny2=ny+ny
nnz2=nz+nz

DO i3=1,nbz2
  DO i2=1,nby2
    DO i1=1,nbx2
      ind=(i3-1)*nxyfn+(i2-1)*nyfn+i1
      rpart(ind)=REAL(vec3d(MOD(i1+nnx2,nbx2)+1,MOD(i2+nny2,nby2)+1,MOD(i3+nnz2,nbz2)+1))
      impart(ind)=AIMAG(vec3d(MOD(i1+nnx2,nbx2)+1,MOD(i2+nyy2,nby2)+1,MOD(i3+nnz2,nbz2)+1))
    END DO
  END DO
END DO

RETURN
END SUBROUTINE copyr3dto1d

! ******************************

SUBROUTINE my_cuda_allocate_coulex()

! ******************************

res = cudaMallocHost (c_p_fftac,kdred*sizeof(size_cmplx))
CALL c_f_pointer(c_p_fftac,fftac,[kxmax,kymax,kzmax])
res = cudaMallocHost (c_p_akvc,kdred*sizeof(size_cmplx))
CALL c_f_pointer(c_p_akvc,akvc,[kxmax,kymax,kzmax])

res = cudaMalloc(c_gpu_fftac,kdred*sizeof(size_cmplx))
CALL c_f_pointer(c_gpu_fftac,gpu_fftac,[kdred])
res = cudaMalloc(c_gpu_akvc,kdred*sizeof(size_cmplx))
CALL c_f_pointer(c_gpu_akvc,gpu_akvc,[kxmax,kymax,kzmax])

RETURN
END SUBROUTINE my_cuda_allocate_coulex

SUBROUTINE coulsolv_end()

res = cudaFreeHost(c_p_fftac)
res = cudaFree(c_gpu_fftac)
res = cudaFree(c_gpu_akvc)
CALL kill_plan(pfft)

DEALLOCATE(xval,yval,zval)
DEALLOCATE(xt2,yt2,zt2)
DEALLOCATE(ikm)

RETURN
END SUBROUTINE coulsolv_end
#endif

#if(fftw_cpu)
SUBROUTINE coulsolv_end()

USE FFTW

CALL fftw_destroy_plan(pforwx)
CALL fftw_destroy_plan(pforwy)
CALL fftw_destroy_plan(pforwz)
CALL fftw_destroy_plan(pbackx)
CALL fftw_destroy_plan(pbacky)
CALL fftw_destroy_plan(pbackz)

DEALLOCATE(xval,yval,zval)
DEALLOCATE(xt2,yt2,zt2)
DEALLOCATE(ikm)

RETURN
END SUBROUTINE coulsolv_end
#endif

END MODULE coulsolv
