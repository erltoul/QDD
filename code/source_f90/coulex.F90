MODULE coulsolv
#if(fftw_cpu)
USE, intrinsic :: iso_c_binding
USE FFTW
USE kinetic, ONLY: FFTW_planflag
#endif
USE params, ONLY: DP,numthr
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

#if(coudoub3D && fftw_cpu)
COMPLEX(DP),ALLOCATABLE,PRIVATE :: akv2(:,:,:)
#else
REAL(DP),ALLOCATABLE,PRIVATE :: akv2r(:),akv2i(:)
#endif
INTEGER,ALLOCATABLE,PRIVATE :: ikm(:,:)
REAL(DP),PRIVATE :: dkx,dky,dkz,akmax,dksp,ecut
INTEGER,PRIVATE :: nxk,nxklo,nxkhi,nksp,nkxyz
!COMMON /kgrid/ akv2r,akv2i, dkx,dky,dkz,akmax,dksp,  &
!    ikm,nxk,nxklo,nxkhi,nksp,nkxyz,ecut

#if(netlib_fft)
REAL(DP),ALLOCATABLE,PRIVATE :: wrkx(:),wrky(:),wrkz(:)
REAL(DP),ALLOCATABLE,PRIVATE :: wsavex(:),wsavey(:),wsavez(:)
INTEGER,ALLOCATABLE,PRIVATE :: ifacx(:),ifacy(:),ifacz(:)
COMPLEX(DP),ALLOCATABLE,PRIVATE :: fftax(:),fftay(:),fftb(:,:)
#endif

#if(fftw_cpu)
INTEGER,PRIVATE,SAVE :: nini=0
INTEGER(C_INT), PRIVATE :: wisdomtest
type(C_PTR), PRIVATE :: pforwx,pforwy,pforwz,pbackx,pbacky,pbackz
COMPLEX(C_DOUBLE_COMPLEX),ALLOCATABLE,PRIVATE :: fftax(:),fftay(:),fftb(:,:)
type(C_PTR), PRIVATE :: pforw,pback,pforwc
COMPLEX(C_DOUBLE_COMPLEX), PRIVATE, ALLOCATABLE :: ffta(:,:,:)
REAL(DP), PRIVATE, ALLOCATABLE :: rffta(:,:,:)
#endif


! include block: option
REAL(DP),PARAMETER,PRIVATE :: zero=0D0
REAL(DP),PARAMETER,PRIVATE :: pi=3.141592653589793D0

!COMMON /fftcom/fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax)

CONTAINS

SUBROUTINE init_coul(dx0,dy0,dz0,nx0,ny0,nz0)

!-----------------------------------------------------------------------

LOGICAL,PARAMETER :: tcoultest=.false.
REAL(DP),ALLOCATABLE :: rhotest(:),ctest(:)

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
ALLOCATE(fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax))
#if(coudoub3D && fftw_cpu)
ALLOCATE(ffta(kxmax,kymax,kzmax),akv2(kxmax,kymax,kzmax))
ALLOCATE(rffta(kxmax,kymax,kzmax))
#else
ALLOCATE(akv2r(kdred),akv2i(kdred))
#endif
ALLOCATE(ikm(kxmax,kymax))
#if(netlib_fft)
ALLOCATE(wrkx(kfft2),wrky(kfft2),wrkz(kfft2))
ALLOCATE(wsavex(kfft2),wsavey(kfft2),wsavez(kfft2))
ALLOCATE(ifacx(kfft2),ifacy(kfft2),ifacz(kfft2))
#endif

#if(coudoub3D && fftw_cpu)
#if(paropenmp)
  call dfftw_init_threads(iret)
!  numthr = 4
  call dfftw_plan_with_nthreads(numthr)
  WRITE(*,*) ' init Coul FFTW threads: iret=',iret,', nr. of threads=',numthr
#endif
IF (nini==0) THEN
#if(fftwnomkl)
  wisdomtest=fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw_coul.dat'//C_NULL_CHAR)
#endif
  IF (wisdomtest == 0) THEN
    wisdomtest = fftw_import_system_wisdom()
    IF(wisdomtest == 0 ) THEN
      WRITE(6,*) 'wisdom_fftw.dat not found, creating it'
      WRITE(7,*) 'wisdom_fftw.dat not found, creating it'
    ELSE
      WRITE(*,*) 'Coulex: wisdom from system'  
    END IF
  END IF
!  pforw=fftw_plan_dft_3d(nz2,ny2,nx2,ffta,ffta,FFTW_FORWARD,FFTW_planflag+FFTW_UNALIGNED)
!  pforw=fftw_plan_dft_3d(kxmax,kymax,kzmax,ffta,ffta,FFTW_FORWARD,FFTW_planflag)
!  pback=fftw_plan_dft_3d(kxmax,kymax,kzmax,ffta,ffta,FFTW_BACKWARD,FFTW_planflag)
  pforw=fftw_plan_dft_r2c_3d(kxmax,kymax,kzmax,rffta,ffta,FFTW_planflag)
  pback=fftw_plan_dft_c2r_3d(kxmax,kymax,kzmax,ffta,rffta,FFTW_planflag)
  pforwc=fftw_plan_dft_3d(kxmax,kymax,kzmax,ffta,ffta,FFTW_FORWARD,FFTW_planflag)
  nini  = kxmax*kymax*kzmax
  WRITE(*,*) ' Coul-Solv initialized nini=',nini,kxmax,kymax,kzmax
ELSE IF(nini /= kxmax*kymax*kzmax) THEN
  WRITE(*,*) ' nini,nx2,ny2,nz2=',nini,nx2,ny2,nz2
  STOP ' nx2, ny2 or/and nz2 in four3d not as initialized!'
END IF
#if(fftwnomkl)
wisdomtest=fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw_coul.dat'//C_NULL_CHAR)
#endif
IF(wisdomtest==0) THEN
  WRITE(*,*) ' export wisdom_fftw_coul.dat failed'
ELSE
  WRITE(*,*) ' export wisdom_fftw_coul.dat successfull'
END IF
CALL fftw_forget_wisdom
#endif

!     call input routine fftinp, which initializes the grid and fft tabl

CALL fftinp

! test section
IF(tcoultest) THEN
  ALLOCATE(rhotest(nx*ny*nz),ctest(nx*ny*nz))
  rhotest = 0D0
  ii = 0
  DO i3=1,nz;  DO i2=1,ny;  DO i1=1,nx
    ii=ii+1
    IF(i3== nz/2 .AND. i2==ny/2 .AND. i1==nx/2) rhotest(ii) =1D0/(dx*dy*dz)
  END DO; END DO; END DO
  charge = SUM(rhotest)*(dx*dy*dz)
  WRITE(*,*) '# test Coulomb for point charge:',charge
  CALL falr(rhotest,ctest,nxdum,nydum,nzdum,kdum)
  ii = 0
  DO i3=1,nz;  DO i2=1,ny;  DO i1=1,nx
    ii=ii+1
    IF(i3== nz/2 .AND. i2==ny/2)  WRITE(*,*) (i1-nx/2)*dx,rhotest(ii), &
      ctest(ii) ! *(i1-nx/2)*dx/2D0
  END DO; END DO; END DO
  STOP
  DEALLOCATE(rhotest,ctest)  
END IF

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

!     initialize grid in Fourier space

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
#if(coudoub3D && fftw_cpu)
factor = 2D0*(dx*dy*dz)/(nx*ny*nz*8D0)  ! 2D0 for e**2
ii=0
DO i3=1,nzi
  IF(i3<=nz) THEN
    xz1=(i3-1)*dz
  ELSE
    xz1=(i3-nzi-1)*dz
  END IF
  xz2=xz1*xz1
  DO i2=1,nyi
    IF(i2<=ny) THEN
      xy1=(i2-1)*dy
    ELSE
      xy1=(i2-nyi-1)*dy
    END IF
    xy2=xy1*xy1
    DO i1=1,nxi
      IF(i1<=nx) THEN
        xx1=(i1-1)*dx
      ELSE
        xx1=(i1-nxi-1)*dx
      END IF
      xx2=xx1*xx1
      ak2=xx2+xy2+xz2
      ii=ii+1
!        write(*,*) ' i1,i2,i3,ii=',i1,i2,i3,ii
      IF(ak2>dx**2/10D0) THEN
        akv2(i1,i2,i3) =  factor/SQRT(ak2) 
        rffta(i1,i2,i3) =  factor/SQRT(ak2) 
      ELSE
        akv2(i1,i2,i3) =  factor*2.34D0*1.19003868D0*(dx*dy*dz)**(-1D0/3D0) 
        rffta(i1,i2,i3) =  factor*2.34D0*1.19003868D0*(dx*dy*dz)**(-1D0/3D0) 
      END IF
    END DO
  END DO
END DO
nksp=ii
#else
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
#endif

#if(coudoub3D && fftw_cpu)
!ffta=akv2
CALL fftw_execute_dft_r2c(pforw,rffta,ffta)
akv2=ffta
#else
CALL fourf(akv2r(1),akv2i(1))
#endif
!  CALL prifld(akv2r,'k**2 for Co')
!  CALL prifld(akv2i,'k**2 for Co')
!nxyf=nx2*ny2
!nyf=nx2
ind=0 !nz*nxyf+ny*nyf
#if(coudoub3D && fftw_cpu)
WRITE(6,*) '1/k**2 along x'
DO i1=1,nxi
  WRITE(6,'(f8.2,2(1pg13.5))') (i1-nx)*dx,akv2(i1,1,1)
END DO
WRITE(6,*) '1/k**2 along z'
DO i3=1,nzi
  WRITE(6,'(f8.2,2(1pg13.5))') (i3-nz)*dz,akv2(1,1,i3)
END DO
#else
WRITE(6,*) '1/k**2 along x'
DO i1=1,nxi
  WRITE(6,'(f8.2,2(1pg13.5))') (i1-nx)*dx,akv2r(i1+ind),akv2i(i1+ind)
END DO
#endif

RETURN
END SUBROUTINE fftinp

#if(coudoub3D && fftw_cpu)
!-------------------------------------------------------------------

SUBROUTINE falr(rhoinp,chpfalr,nxdum,nydum,nzdum,kdum)
IMPLICIT REAL(DP) (A-H,O-Z)

! Coulomb solver using FFTW


REAL(DP), INTENT(IN)                     :: rhoinp(kdfull)
REAL(DP), INTENT(OUT)                     :: chpfalr(kdfull)
INTEGER, INTENT(IN)                  :: nxdum
INTEGER, INTENT(IN)                  :: nydum
INTEGER, INTENT(IN)                  :: nzdum
INTEGER, INTENT(IN)                  :: kdum

!WRITE(*,*) ' running new Coulomb'

! map density to 2**3 larger grid, immediately on FFTW3 work space

rffta = 0D0
i0=0
DO i3=1,nz
  DO i2=1,ny
    DO i1=1,nx
      i0 = i0+1
      rffta(i1,i2,i3)=rhoinp(i0)
    END DO
  END DO
END DO

! FFT forward, multiply with Green's function, FFT backward

CALL fftw_execute_dft_r2c(pforw,rffta,ffta)

ffta = akv2*ffta

CALL fftw_execute_dft_c2r(pback,ffta,rffta)

! map back to standard grid, augment with normalization factor


!facnr =SQRT(8D0*pi*pi*pi)/SQRT(REAL(kxmax*kymax*kzmax,DP)) * 2D0
i0=0
DO i3=1,nz
  DO i2=1,ny
    DO i1=1,nx
      i0 = i0+1
      chpfalr(i0) = rffta(i1,i2,i3)   ! *facnr   ! all scaling in 'akv2'
    END DO
  END DO
END DO


END SUBROUTINE falr
#else
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

!WRITE(*,*) ' running old Coulomb'

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

!     copy density on complex array of double extension in x,y,z


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

!     Fourier transformation of the density

CALL fourf(rhokr,rhoki)

!     calculation of the Coulomb field (writing on the density field)

DO ik=1,kdred
  SAVE2     = akv2r(ik)*rhokr(ik)+akv2i(ik)*rhoki(ik)
  rhoki(ik) = akv2r(ik)*rhoki(ik)+akv2i(ik)*rhokr(ik)
  rhokr(ik) = SAVE2
END DO

!     Fourier back transformation

CALL fourb(rhokr,rhoki)

RETURN

END SUBROUTINE coufou2

#endif



!-----fourf-------------------------------------------------------fourf

SUBROUTINE fourf(pskr,pski)
#if(fftw_cpu)
USE FFTW
#endif
#if(parayes)
USE params, only : myn
#endif
IMPLICIT REAL(DP) (A-H,O-Z)
#if(parayes)
INCLUDE 'mpif.h'
REAL(DP) :: is(mpi_status_size)
#endif

REAL(DP), INTENT(OUT)                        :: pskr(kdred)
REAL(DP), INTENT(OUT)                        :: pski(kdred)


!     Fourier forward transformation
!     I/O: pskr   real part of the wave-function
!          pski   imaginary part of the wave-function

DATA  mxini,myini,mzini/0,0,0/              ! flag for initialization

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
#if(parano)
IF(mxini == 0) THEN
#if(fftwnomkl)
  wisdomtest=fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
#endif
  IF (wisdomtest == 0) THEN
    WRITE(6,*) 'wisdom_fftw.dat not found, creating it'
    WRITE(7,*) 'wisdom_fftw.dat not found, creating it'
  END IF
  pforwx=fftw_plan_dft_1d(kfftx,fftax,fftax,FFTW_FORWARD,FFTW_planflag)
  pbackx=fftw_plan_dft_1d(kfftx,fftax,fftax,FFTW_BACKWARD,FFTW_planflag)
  mxini  = kfftx
  WRITE(7,'(a)') ' x-fft initialized '
ELSE IF(mxini /= kfftx) THEN
  STOP ' nx2 in four3d not as initialized!'
END IF
IF(myini == 0) THEN
  pforwy=fftw_plan_dft_1d(kffty,fftay,fftay,FFTW_FORWARD,FFTW_planflag)
  pbacky=fftw_plan_dft_1d(kffty,fftay,fftay,FFTW_BACKWARD,FFTW_planflag)
  myini  = kffty
  WRITE(7,'(a)') ' y-fft initialized '
ELSE IF(myini /= kffty) THEN
  STOP ' ny2 in four3d not as initialized!'
END IF
IF(mzini == 0) THEN
  pforwz=fftw_plan_dft_1d(kfftz,fftb,fftb,FFTW_FORWARD,FFTW_planflag)
  pbackz=fftw_plan_dft_1d(kfftz,fftb,fftb,FFTW_BACKWARD,FFTW_planflag)
  mzini  = kfftz
  IF (wisdomtest == 0) THEN
    WRITE(6,*) 'Error exporting wisdom to file wisdom_fftw.dat'
    WRITE(7,*) 'Error exporting wisdom to file wisdom_fftw.dat'
  END IF
  WRITE(7,'(a)') ' z-fft initialized '
ELSE IF(mzini /= kfftz) THEN
  STOP ' nz2 in four3d not as initialized!'
END IF
#if(fftwnomkl)
wisdomtest=fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw_coul.dat'//C_NULL_CHAR)
#endif
IF(wisdomtest==0) THEN
  WRITE(*,*) ' export wisdom_fftw_coul.dat failed'
ELSE
  WRITE(*,*) ' export wisdom_fftw_coul.dat successfull'
END IF
CALL fftw_forget_wisdom
#endif

#if(parayes)
IF(mxini == 0) THEN
  IF(myn == 0) THEN
    !Master node creates wisdom if necessary...
#if(fftwnomkl)
    wisdomtest=fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
#endif
    IF (wisdomtest == 0) THEN
      WRITE(6,*) 'wisdom_fftw.dat not found, creating it'
      WRITE(7,*) 'wisdom_fftw.dat not found, creating it'
    END IF
    pforwx=fftw_plan_dft_1d(kfftx,fftax,fftax,FFTW_FORWARD,FFTW_planflag)
    pbackx=fftw_plan_dft_1d(kfftx,fftax,fftax,FFTW_BACKWARD,FFTW_planflag)
    mxini  = kfftx
#if(fftwnomkl)
    wisdomtest=fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
#endif
    IF (wisdomtest == 0) THEN
      WRITE(6,*) 'Error exporting wisdom to file wisdom_fftw.dat'
      WRITE(7,*) 'Error exporting wisdom to file wisdom_fftw.dat'
    END IF
  ENDIF
  CALL mpi_barrier(mpi_comm_world,mpi_ierror)
  !... then other nodes use it
  IF(myn /=0 ) THEN
#if(fftwnomkl)
    wisdomtest=fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
#endif
    pforwx=fftw_plan_dft_1d(kfftx,fftax,fftax,FFTW_FORWARD,FFTW_planflag)
    pbackx=fftw_plan_dft_1d(kfftx,fftax,fftax,FFTW_BACKWARD,FFTW_planflag)
    mxini  = kfftx
  ENDIF
ELSE IF(mxini /= kfftx) THEN
  STOP ' nx2 in four3d not as initialized!'
END IF
IF(myini == 0) THEN
  IF(myn == 0) THEN
    pforwy=fftw_plan_dft_1d(kffty,fftay,fftay,FFTW_FORWARD,FFTW_planflag)
    pbacky=fftw_plan_dft_1d(kffty,fftay,fftay,FFTW_BACKWARD,FFTW_planflag)
    myini  = kffty
#if(fftwnomkl)
    wisdomtest=fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
#endif
    IF (wisdomtest == 0) THEN
      WRITE(6,*) 'Error exporting wisdom to file wisdom_fftw.dat'
      WRITE(7,*) 'Error exporting wisdom to file wisdom_fftw.dat'
    END IF
  ENDIF
  CALL mpi_barrier(mpi_comm_world,mpi_ierror)
  IF(myn /= 0) THEN
#if(fftwnomkl)
    wisdomtest=fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
#endif
    pforwy=fftw_plan_dft_1d(kffty,fftay,fftay,FFTW_FORWARD,FFTW_planflag)
    pbacky=fftw_plan_dft_1d(kffty,fftay,fftay,FFTW_BACKWARD,FFTW_planflag)
    myini  = kffty
  ENDIF
ELSE IF(myini /= kffty) THEN
  STOP ' ny2 in four3d not as initialized!'
END IF
IF(mzini == 0) THEN
  IF(myn == 0) THEN
    pforwz=fftw_plan_dft_1d(kfftz,fftb,fftb,FFTW_FORWARD,FFTW_planflag)
    pbackz=fftw_plan_dft_1d(kfftz,fftb,fftb,FFTW_BACKWARD,FFTW_planflag)
    mzini  = kfftz
#if(fftwnomkl)
    wisdomtest=fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
#endif
    IF (wisdomtest == 0) THEN
      WRITE(6,*) 'Error exporting wisdom to file wisdom_fftw.dat'
      WRITE(7,*) 'Error exporting wisdom to file wisdom_fftw.dat'
    END IF
  ENDIF
  CALL mpi_barrier(mpi_comm_world,mpi_ierror)
  IF(myn /= 0) THEN
#if(fftwnomkl)
    wisdomtest=fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
#endif
    pforwz=fftw_plan_dft_1d(kfftz,fftb,fftb,FFTW_FORWARD,FFTW_planflag)
    pbackz=fftw_plan_dft_1d(kfftz,fftb,fftb,FFTW_BACKWARD,FFTW_planflag)
    mzini  = kfftz
  ENDIF
ELSE IF(mzini /= kfftz) THEN
  STOP ' nz2 in four3d not as initialized!'
END IF
#endif
#endif


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


!     Fourier backward transformation
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
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: psxr(kdred)
REAL(DP), INTENT(IN OUT)                        :: psxi(kdred)

!     performs the Fourier-transformation in x-direction
!     the input-wave-function (psxr,psxi) (i.e. real and imaginary part)
!     is overwritten by the Fourier-transformed wave-function

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
    
!         execution of the Fourier-transformation

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

!     performs the Fourier-transformation in y-direction
!     the input-wave-function (psxr,psxi) (i.e. real and imaginary part)
!     is overwritten by the Fourier-transformed wave-function
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
    
!         execution of the Fourier-transformation

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

!     performs the Fourier-transformation in z-direction
!     the input-wave-function (psxr,psxi) (i.e. real and imaginary part)
!     is overwritten by the Fourier-transformed wave-function

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
DEALLOCATE(fftax,fftay,fftb)
DEALLOCATE(ikm)

RETURN
END SUBROUTINE coulsolv_end
#endif

END MODULE coulsolv
