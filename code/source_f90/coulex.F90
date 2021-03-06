MODULE coulsolv
USE params, ONLY: DP
IMPLICIT REAL(DP) (A-H,O-Z)

SAVE
INTEGER,PRIVATE :: kxmax,kymax,kzmax, ksmax
! kxamx must be the largest
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

REAL(DP),ALLOCATABLE,PRIVATE :: akv2r(:),akv2i(:)
INTEGER,ALLOCATABLE,PRIVATE :: ikm(:,:)
REAL(DP),PRIVATE :: dkx,dky,dkz,akmax,dksp,ecut
INTEGER,PRIVATE :: nxk,nxklo,nxkhi,nksp,nkxyz
!COMMON /kgrid/ akv2r,akv2i, dkx,dky,dkz,akmax,dksp,  &
!    ikm,nxk,nxklo,nxkhi,nksp,nkxyz,ecut

REAL(DP),ALLOCATABLE,PRIVATE :: wrkx(:),wrky(:),wrkz(:)
REAL(DP),ALLOCATABLE,PRIVATE :: wsavex(:),wsavey(:),wsavez(:)
INTEGER,ALLOCATABLE,PRIVATE :: ifacx(:),ifacy(:),ifacz(:)
!COMMON /fftini/wrkx,wrky,wrkz,wsavex,wsavey,wsavez,ifacx,ifacy,ifacz

! include block: option
REAL(DP),PARAMETER,PRIVATE :: zero=0D0
REAL(DP),PARAMETER,PRIVATE :: pi=3.141592653589793D0

COMPLEX(DP),ALLOCATABLE,PRIVATE :: fftax(:),fftay(:),fftb(:,:)
!COMMON /fftcom/fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax)

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
ALLOCATE(wrkx(kfft2),wrky(kfft2),wrkz(kfft2))
ALLOCATE(wsavex(kfft2),wsavey(kfft2),wsavez(kfft2))
ALLOCATE(ifacx(kfft2),ifacy(kfft2),ifacz(kfft2))
ALLOCATE(fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax))
ALLOCATE(akv2r(kdred),akv2i(kdred))
ALLOCATE(ikm(kxmax,kymax))


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
WRITE(*,*) ' dkx,dky,dkz,dksp=',dkx,dky,dkz,dksp

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



RETURN
END SUBROUTINE fftinp

!-------------------------------------------------------------------

SUBROUTINE falr(rhoinp,chpfalr,nxdum,nydum,nzdum,kdum)

IMPLICIT REAL(DP) (A-H,O-Z)


REAL(DP), INTENT(IN)                     :: rhoinp(kdfull)
REAL(DP), INTENT(OUT)                     :: chpfalr(kdfull)
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

DO ik=1,kdred
  SAVE2     = akv2r(ik)*rhokr(ik)+akv2i(ik)*rhoki(ik)
  rhoki(ik) = akv2r(ik)*rhoki(ik)+akv2i(ik)*rhokr(ik)
  rhokr(ik) = SAVE2
END DO

!     fourier back transformation

CALL fourb(rhokr,rhoki)

RETURN

END SUBROUTINE coufou2




!-----fourf-------------------------------------------------------fourf

SUBROUTINE fourf(pskr,pski)
IMPLICIT REAL(DP) (A-H,O-Z)


REAL(DP), INTENT(OUT)                        :: pskr(kdred)
REAL(DP), INTENT(OUT)                        :: pski(kdred)


!     fourier forward transformation
!     I/O: pskr   real part of the wave-function
!          pski   imaginary part of the wave-function

DATA  mxini,myini,mzini/0,0,0/              ! flag for initialization

!----------------------------------------------------------------------


!     check initialization

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
nzzh=(nzi-1)*nxy1
nyyh=(nyi-1)*nxi
!test      sqh=sqrt(0.5)
tnorm=grnorm*fnorm

CALL fftx(pskr,pski)
CALL ffty(pskr,pski)
CALL fftz(pskr,pski)


!test2      write(6,'(a)') ' TNORM=',tnorm
DO i1=1,nkxyz
  pskr(i1)=tnorm*pskr(i1)
  pski(i1)=tnorm*pski(i1)
END DO
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
tnorm=fnorm/(8D0*grnorm)*pi**1.5D0
!      tnorm=fnorm/(16D0*grnorm)

DO i=1,nkxyz
  pskr(i)=tnorm*pskr(i)
  pski(i)=tnorm*pski(i)
END DO


CALL ffbz(pskr,pski)
CALL ffby(pskr,pski)
CALL ffbx(pskr,pski)


RETURN
END SUBROUTINE fourb


!-----fftx-------------------------------------------------------------

SUBROUTINE fftx(psxr,psxi)
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
    
    CALL dcftf1 (kfftx,fftax,wrkx,wsavex,ifacx)
    
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
    CALL dcftf1 (kffty,fftay,wrky,wsavey,ifacy)
    
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
    CALL dcftf1 (kfftz,fftb(1,i1),wrkz,wsavez,ifacz)
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
    CALL dcftb1 (kfftz,fftb(1,i1),wrkz,wsavez,ifacz)    ! basic fft
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
    CALL dcftb1 (kffty,fftay,wrky,wsavey,ifacy)
    
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
    CALL dcftb1 (kfftx,fftax,wrkx,wsavex,ifacx)
    
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

END MODULE coulsolv
