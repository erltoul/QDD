MODULE kinetic
#if(fftw_cpu|fftw_gpu)
USE, intrinsic :: iso_c_binding
#endif
USE params, ONLY: DP
#if(fftw_gpu)
USE cuda_alloc
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

SAVE
!     arrays for kinetic energy and electronic propagation
!       kinetic energy coefficients in strange ordered fourier space
!     akv  = fouier-field for 0.5*k^2
!     ak   = fourier-field for exp(i*dt*(h^2/2m)*k^2)

COMPLEX(DP),ALLOCATABLE :: akpropx(:),akpropy(:),akpropz(:)
INTEGER,PRIVATE,ALLOCATABLE :: modx(:),mody(:),modz(:)
COMPLEX(DP),PARAMETER,PRIVATE :: eye=(0D0,1D0)
REAL(DP),PARAMETER,PRIVATE :: PI=3.141592653589793D0


INTEGER, PRIVATE :: kfft,kfft2,kdfull2

#if(netlib_fft)
COMPLEX(DP),ALLOCATABLE :: ak(:)
REAL(DP),ALLOCATABLE :: akv(:)
COMPLEX(DP), ALLOCATABLE :: akx(:),aky(:),akz(:)
COMPLEX(DP), PRIVATE, ALLOCATABLE :: fftax(:),fftay(:),fftb(:,:)
REAL(DP), PRIVATE, ALLOCATABLE :: wrkx(:),wrky(:),wrkz(:)
REAL(DP), PRIVATE, ALLOCATABLE :: wsavex(:),wsavey(:),wsavez(:)
INTEGER, PRIVATE, ALLOCATABLE :: ifacx(:),ifacy(:),ifacz(:)
#endif
#if(fftw_cpu)
COMPLEX(DP),ALLOCATABLE :: ak(:)
REAL(DP),ALLOCATABLE :: akv(:)
COMPLEX(DP), ALLOCATABLE :: akx(:),aky(:),akz(:)
COMPLEX(C_DOUBLE_COMPLEX), PRIVATE, ALLOCATABLE :: fftax(:),fftay(:),fftaz(:),fftb(:,:),ffta(:,:,:)
type(C_PTR), PRIVATE :: pforwx,pforwy,pforwz,pforwz1,pbackx,pbacky,pbackz,pbackz1
type(C_PTR), PRIVATE :: pforw,pback
INTEGER(C_INT), PRIVATE :: wisdomtest
#endif
#if(fftw_gpu)
COMPLEX(C_DOUBLE_COMPLEX), POINTER :: fftax(:),fftay(:),fftaz(:),fftb(:,:),ffta(:,:,:),ffta2(:,:,:)
TYPE(C_PTR), PRIVATE :: c_p_fftax,c_p_fftay,c_p_fftaz,c_p_fftb,c_p_ffta,c_p_ffta2
COMPLEX(C_DOUBLE_COMPLEX), POINTER :: gpu_ffta(:,:,:),gpu_ffta2(:,:,:),gpu_ffta_int(:,:,:)
TYPE(C_PTR),PRIVATE :: c_gpu_ffta,c_gpu_ffta2,c_gpu_ffta_int

INTEGER*8, PRIVATE :: px,py,pz
INTEGER*8,PRIVATE :: p

INTEGER, PARAMETER, PRIVATE :: batch=1
INTEGER, PRIVATE :: res
COMPLEX(DP), PARAMETER :: size_cmplx=CMPLX(0D0,0D0,DP)
REAL(DP), PARAMETER :: size_double=(0D0,DP)

COMPLEX(C_DOUBLE_COMPLEX), POINTER :: gpu_akfft(:,:,:)
TYPE(C_PTR),PRIVATE :: c_gpu_akfft

REAL(C_DOUBLE),POINTER :: gpu_akvfft(:,:,:)
TYPE(C_PTR),PRIVATE :: c_gpu_akvfft

COMPLEX(C_DOUBLE_COMPLEX),POINTER :: gpu_akxfft(:),gpu_akyfft(:),gpu_akzfft(:)
TYPE(C_PTR),PRIVATE :: c_gpu_akxfft,c_gpu_akyfft,c_gpu_akzfft
#endif

CONTAINS
!-----init_grid_fft-----------------------------------------------------

SUBROUTINE init_grid_fft(dx0,dy0,dz0,nx0,ny0,nz0,dt1,h2m)

!     initialize details for FFT

#if(fftw_cpu)
USE FFTW
#endif

REAL(DP) :: dt1,h2m
#if(fftw_cpu)
INTEGER, SAVE ::  nxini=0,nyini=0,nzini=0,nini=0 ! flag for initialization
#endif

nx2=nx0;ny2=ny0;nz2=nz0
kxmax=nx0;kymax=ny0;kzmax=nz0
nx=nx2/2;ny=ny2/2;nz=nz2/2

kfft=2*kxmax
kfft2=kfft*2+1
kdfull2=kxmax*kymax*kzmax

dkx=pi/(dx0*nx)
dky=pi/(dy0*ny)
dkz=pi/(dz0*nz)

ALLOCATE(akpropx(kxmax),akpropy(kymax),akpropz(kzmax))
ALLOCATE(modx(kxmax),mody(kymax),modz(kzmax))

#if(netlib_fft)
ALLOCATE(ak(kdfull2),akv(kdfull2))
ALLOCATE(akx(kdfull2),aky(kdfull2),akz(kdfull2))!,rakx(kdfull2),raky(kdfull2),rakz(kdfull2))
ALLOCATE(fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax))
ALLOCATE(wrkx(kfft2),wrky(kfft2),wrkz(kfft2))
ALLOCATE(wsavex(kfft2),wsavey(kfft2),wsavez(kfft2))
ALLOCATE(ifacx(kfft2),ifacy(kfft2),ifacz(kfft2))
#endif
#if(fftw_cpu)
ALLOCATE(ak(kdfull2),akv(kdfull2))
ALLOCATE(akx(kdfull2),aky(kdfull2),akz(kdfull2))!,rakx(kdfull2),raky(kdfull2),rakz(kdfull2))
ALLOCATE(fftax(kxmax),fftay(kymax),fftaz(kzmax),fftb(kzmax,kxmax),ffta(kxmax,kymax,kzmax))
#endif
#if(fftw_gpu)
CALL my_cuda_allocate(nx2,ny2,nz2) !Pinned memory allocation to make CPU>GPU and GPU>CPU transfers faster
#endif

WRITE(7,*) 'h bar squared over two m electron',h2m
WRITE(7,*) ' testprint EYE=',eye
WRITE(7,*) ' dkx,dky,dkz=',dkx,dky,dkz
WRITE(7,*) ' testprint: nx2,ny2,nz2=',nx2,ny2,nz2
WRITE(7,*) ' testprint: kdfull2,kfft2=',kdfull2,kfft2

!     prepare k**2 and kinetic propagation factor in 3D momentum space
#if(netlib_fft|fftw_cpu)
ind=0
DO i3=1,nz2
  IF(i3 >= (nz+1)) THEN
    zkz=(i3-nz2-1)*dkz
  ELSE
    zkz=(i3-1)*dkz
  END IF
  DO i2=1,ny2
    IF(i2 >= (ny+1)) THEN
      zky=(i2-ny2-1)*dky
    ELSE
      zky=(i2-1)*dky
    END IF
    DO i1=1,nx2
      IF(i1 >= (nx+1)) THEN
        zkx=(i1-nx2-1)*dkx
      ELSE
        zkx=(i1-1)*dkx
      END IF
      ind=ind+1
      ak(ind)=EXP(-eye*dt1*(zkx**2+zky**2+zkz**2)*h2m)
      akv(ind)=(zkx**2+zky**2+zkz**2)*h2m
      akx(ind)=-zkx*eye
      aky(ind)=-zky*eye
      akz(ind)=-zkz*eye
    END DO
  END DO
END DO
#endif

#if(fftw_gpu)
CALL build_kgpu(nx2,ny2,nz2,h2m,dt1,dkx,dky,dkz,gpu_akfft,gpu_akvfft,gpu_akxfft,gpu_akyfft,gpu_akzfft)
#endif
!     prepare kinetic propagation factors in 1D momentum spaces

DO i3=1,nz2
  IF(i3 >= (nz+1)) THEN
    zkz=(i3-nz2-1)*dkz
  ELSE
    zkz=(i3-1)*dkz
  END IF
  akpropz(i3)=EXP(-eye*dt1*zkz**2*h2m)
END DO

DO i2=1,ny2
  IF(i2 >= (ny+1)) THEN
    zky=(i2-ny2-1)*dky
  ELSE
    zky=(i2-1)*dky
  END IF
  akpropy(i2)=EXP(-eye*dt1*zky**2*h2m)
END DO

DO i1=1,nx2
  IF(i1 >= (nx+1)) THEN
    zkx=(i1-nx2-1)*dkx
  ELSE
    zkx=(i1-1)*dkx
  END IF
  akpropx(i1)=EXP(-eye*dt1*zkx**2*h2m)
END DO

#if(fftw_cpu)
IF (nini==0) THEN
  wisdomtest=fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
  IF (wisdomtest == 0) THEN
    WRITE(6,*) 'wisdom_fftw.dat not found, creating it'
    WRITE(7,*) 'wisdom_fftw.dat not found, creating it'
  END IF
  pforw=fftw_plan_dft_3d(nz2,ny2,nx2,ffta,ffta,FFTW_FORWARD,FFTW_EXHAUSTIVE)
  pback=fftw_plan_dft_3d(nz2,ny2,nx2,ffta,ffta,FFTW_BACKWARD,FFTW_EXHAUSTIVE)
  nini  = nx2*ny2*nz2
ELSE IF(nini /= nx2*ny2*nz2) THEN
  STOP ' nx2, ny2 or/and nz2 in four3d not as initialized!'
END IF
IF(nxini == 0) THEN
  pforwx=fftw_plan_dft_1d(nx2,fftax,fftax,FFTW_FORWARD,FFTW_EXHAUSTIVE)
  pbackx=fftw_plan_dft_1d(nx2,fftax,fftax,FFTW_BACKWARD,FFTW_EXHAUSTIVE)
  nxini  = nx2
!       write(6,'(a)') ' x-fft initialized '
ELSE IF(nxini /= nx2) THEN
  STOP ' nx2 in four3d not as initialized!'
END IF
IF(nyini == 0) THEN
  pforwy=fftw_plan_dft_1d(ny2,fftay,fftay,FFTW_FORWARD,FFTW_EXHAUSTIVE)
  pbacky=fftw_plan_dft_1d(ny2,fftay,fftay,FFTW_BACKWARD,FFTW_EXHAUSTIVE)
  nyini  = ny2
!       write(6,'(a)') ' y-fft initialized '
ELSE IF(nyini /= ny2) THEN
  STOP ' ny2 in four3d not as initialized!'
END IF
IF(nzini == 0) THEN
  pforwz=fftw_plan_dft_1d(nz2,fftb(1,i1),fftb(1,i1),FFTW_FORWARD,FFTW_EXHAUSTIVE)
  pbackz=fftw_plan_dft_1d(nz2,fftb(1,i1),fftb(1,i1),FFTW_BACKWARD,FFTW_EXHAUSTIVE)
  pforwz1=fftw_plan_dft_1d(nz2,fftaz,fftaz,FFTW_FORWARD,FFTW_EXHAUSTIVE)
  pbackz1=fftw_plan_dft_1d(nz2,fftaz,fftaz,FFTW_BACKWARD,FFTW_EXHAUSTIVE)
  nzini  = nz2
  wisdomtest=fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
  IF (wisdomtest == 0) THEN
    WRITE(6,*) 'Error exporting wisdom to file wisdom_fftw.dat'
    WRITE(7,*) 'Error exporting wisdom to file wisdom_fftw.dat'
  END IF
!       write(6,'(a)') ' z-fft initialized '
ELSE IF(nzini /= nz2) THEN
  STOP ' nz2 in four3d not as initialized!'
END IF
#endif

#if(fftw_gpu)
IF (nini==0) THEN
  CALL cuda_plan_3d(p,nx2,ny2,nz2)
  nini  = nx2*ny2*nz2
ELSE IF(nini /= nx2*ny2*nz2) THEN
  STOP ' nx2, ny2 or/and nz2 in four3d not as initialized!'
END IF
IF(nxini == 0) THEN
  CALL cuda_plan_1d(px,nx2,batch)
  nxini  = nx2
!       write(6,'(a)') ' x-fft initialized '
ELSE IF(nxini /= nx2) THEN
  STOP ' nx2 in four3d not as initialized!'
END IF
IF(nyini == 0) THEN 
  CALL cuda_plan_1d(py,ny2,batch)
  nyini  = ny2
!       write(6,'(a)') ' y-fft initialized '
ELSE IF(nyini /= ny2) THEN
  STOP ' ny2 in four3d not as initialized!'
END IF
IF(nzini == 0) THEN
  CALL cuda_plan_1d(pz,nz2,batch)
  nzini  = nz2
!       write(6,'(a)') ' z-fft initialized '
ELSE IF(nzini /= nz2) THEN
  STOP ' nz2 in four3d not as initialized!'
END IF
#endif

DO i1=1,nx2
  modx(i1)=MOD(i1+nx,nx2)+1
END DO

DO i2=1,ny2
  mody(i2)=MOD(i2+ny,ny2)+1
END DO

DO i3=1,nz2
  modz(i3)=MOD(i3+nz,nz2)+1
END DO

END SUBROUTINE init_grid_fft

! ******************************

SUBROUTINE  kinprop(q1,q2)

! ******************************

!       propagation with exp(-i*dt*e_kin)

USE params
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN OUT)                  :: q1(kdfull2)
COMPLEX(DP), INTENT(OUT)                     :: q2(kdfull2)

#if(netlib_fft)
DATA  nxini,nyini,nzini/0,0,0/ ! flag for initialization
#endif

tnorm=1D0/SQRT(8D0*pi*pi*pi*REAL(nx2*ny2*nz2,DP))

#if(fftw_gpu)
STOP 'fastpropag not yet implemented with GPU'
#endif
!       check initialization
#if(netlib_fft)
IF(nxini == 0) THEN
  CALL dcfti1 (nx2,wsavex,ifacx)
  nxini  = nx2
!       write(6,'(a)') ' x-fft initialized '
ELSE IF(nxini /= nx2) THEN
  STOP ' nx2 in four3d not as initialized!'
END IF
IF(nyini == 0) THEN
  CALL dcfti1 (ny2,wsavey,ifacy)
  nyini  = ny2
!       write(6,'(a)') ' y-fft initialized '
ELSE IF(nyini /= ny2) THEN
  STOP ' ny2 in four3d not as initialized!'
END IF
IF(nzini == 0) THEN
  CALL dcfti1 (nz2,wsavez,ifacz)
  nzini  = nz2
!       write(6,'(a)') ' z-fft initialized '
ELSE IF(nzini /= nz2) THEN
  STOP ' nz2 in four3d not as initialized!'
END IF
#endif

!       propagation in x-direction

xfnorm = 1D0/nx2
DO i3=1,nz2
  DO i2=1,ny2
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftax(modx(i1))=q1(ind) ! copy to workspace
    END DO
#if(netlib_fft)
    CALL dcftf1 (nx2,fftax,wrkx,wsavex,ifacx) ! basic fft
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pforwx,fftax,fftax)
#endif
#if(fftw_gpu)
    CALL run_fft_for(px,fftax,fftax,nx2)
#endif
    DO i1=1,nx2
      fftax(i1) = akpropx(i1)*fftax(i1)
    END DO
#if(netlib_fft)
    CALL dcftb1 (nx2,fftax,wrkx,wsavex,ifacx) ! basic back fft
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pbackx,fftax,fftax)
#endif
#if(fftw_gpu)
    CALL run_fft_back(px,fftax,fftax,nx2)
#endif
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q2(ind)= fftax(modx(i1))*xfnorm
    END DO
  END DO
END DO

!      transformation in y-direction

yfnorm = 1D0/ny2
DO i3=1,nz2
  DO i1=1,nx2
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftay(mody(i2)) = q2(ind)
    END DO
#if(netlib_fft)
    CALL dcftf1 (ny2,fftay,wrky,wsavey,ifacy)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pforwy,fftay,fftay)
#endif
#if(fftw_gpu)
    CALL run_fft_for(py,fftay,fftay,ny2)
#endif
    DO i2=1,ny2
      fftay(i2) = akpropy(i2)*fftay(i2)
    END DO
#if(netlib_fft)
    CALL dcftb1 (ny2,fftay,wrky,wsavey,ifacy)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pbacky,fftay,fftay)
#endif
#if(fftw_gpu)
    CALL run_fft_back(py,fftay,fftay,ny2)
#endif
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q2(ind)= fftay(mody(i2))*yfnorm
    END DO
  END DO
END DO

!       propagation in z-direction

zfnorm = 1D0/nz2
DO i2=1,ny2
  DO i3=1,nz2
    i3m = modz(i3)
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftb(i3m,i1) = q2(ind)
    END DO
  END DO
  DO i1=1,nx2
#if(netlib_fft)
    CALL dcftf1 (nz2,fftb(1,i1),wrkz,wsavez,ifacz)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pforwz,fftb(1,i1),fftb(1,i1))
#endif
#if(fftw_gpu)
    CALL run_fft_for(pz,fftb(1,i1),fftb(1,i1),nz2)
#endif
    DO i3=1,nz2
      fftb(i3,i1) = akpropz(i3)*fftb(i3,i1)
    END DO
#if(netlib_fft)
    CALL dcftb1 (nz2,fftb(1,i1),wrkz,wsavez,ifacz)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pbackz,fftb(1,i1),fftb(1,i1))
#endif
#if(fftw_gpu)
    CALL run_fft_back(pz,fftb(1,i1),fftb(1,i1),nz2)
#endif
  END DO
  DO i3=1,nz2
    i3m = modz(i3)
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q1(ind)= fftb(i3m,i1)*zfnorm
    END DO
  END DO
END DO

RETURN
END SUBROUTINE  kinprop


! ******************************

SUBROUTINE  gradient(fin,gradfout,idirec)

! ******************************


!  The gradient of complex field 'fin' in direction
!  'idirec' (x =1, y=2, z=3).
!  The fields are given in Fourier space and the gradient is
!   applied as product with 'kx', 'ky', 'kz'.
!

USE params
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN OUT)                      :: fin(kdfull2)
COMPLEX(DP), INTENT(IN OUT)                     :: gradfout(kdfull2)
INTEGER, INTENT(IN)                          :: idirec


! ************************************************************

#if(netlib_fft|fftw_cpu)
  IF(idirec == 1) THEN
!       x-derivative

!    dkx=pi/(dx*REAL(nx,DP))
!    ind=0
!    DO i3=1,nz2
!      DO i2=1,ny2
!        DO i1=1,nx2
!          IF(i1 >= (nx+1)) THEN
!            zkx=(i1-nx2-1)*dkx
!          ELSE
!            zkx=(i1-1)*dkx
!          END IF
!          ind=ind+1
!          gradfout(ind) = eye*zkx*fin(ind)
          gradfout = -akx*fin
!        END DO
!      END DO
!    END DO

  ELSEIF(idirec == 2) THEN
!       y-derivative

!    dky=pi/(dy*REAL(ny,DP))
!    ind=0
!    DO i3=1,nz2
!      DO i2=1,ny2
!        IF(i2 >= (ny+1)) THEN
!          zky=(i2-ny2-1)*dky
!        ELSE
!          zky=(i2-1)*dky
!        END IF
!        DO i1=1,nx2
!          ind=ind+1
!          gradfout(ind) = eye*zky*fin(ind)
          gradfout = -aky*fin
!        END DO
!      END DO
!    END DO

  ELSEIF(idirec == 3) THEN
!       z-derivative

!    ind=0
!    dkz=pi/(dz*REAL(nz,DP))
!    DO i3=1,nz2
!      IF(i3 >= (nz+1)) THEN
!        zkz=(i3-nz2-1)*dkz
!      ELSE
!        zkz=(i3-1)*dkz
!      END IF
!      DO i2=1,ny2
!        DO i1=1,nx2
!          ind=ind+1
!          gradfout(ind) = eye*zkz*fin(ind)
          gradfout = -akz*fin
!        END DO
!      END DO
!    END DO

  ELSE
    STOP ' RGRADIENT called with invalid IDIREC'
  ENDIF
#endif

#if(fftw_gpu)
!gradient computed with the "multiply_ak2" cuda function
#endif

RETURN
END SUBROUTINE  gradient

! ******************************

SUBROUTINE  xgradient_rspace(fin,gradfout)

! ******************************

!  The gradient of the complex field 'fin' in x-direction.
!  The fields are given in coordinate space. The gradient is
!  evaluated in k_x space.
!

USE params
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN OUT)                      :: fin(kdfull2)
COMPLEX(DP), INTENT(IN OUT)                     :: gradfout(kdfull2)

! ************************************************************

dkx=pi/(dx*nx)
DO i3=1,nz2
  DO i2=1,ny2
!                 forward transform along x
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftax(modx(i1))=fin(ind)  ! copy to workspace
    END DO
#if(netlib_fft)
    CALL dcftf1 (nx2,fftax,wrkx,wsavex,ifacx)    ! basic fft
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pforwx,fftax,fftax)    ! basic fft
#endif
#if(fftw_gpu)
    CALL run_fft_for(px,fftax,fftax,nx2)    ! basic fft
#endif
!                 multiply by k_x factor in k_x space
    DO i1=1,nx2
      IF(i1 >= (nx+1)) THEN
        zkx=(i1-nx2-1)*dkx
      ELSE
        zkx=(i1-1)*dkx
      END IF
      fftax(i1) = fftax(i1)*eye*zkx
    END DO
!                 backward transform along x
#if(netlib_fft)
    CALL dcftb1 (nx2,fftax,wrkx,wsavex,ifacx)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pbackx,fftax,fftax)
#endif
#if(fftw_gpu)
    CALL run_fft_back(px,fftax,fftax,nx2)
#endif
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      gradfout(ind)= fftax(modx(i1))/nx2
    END DO
  END DO
END DO

RETURN

END SUBROUTINE  xgradient_rspace


! ******************************

SUBROUTINE  ygradient_rspace(fin,gradfout)

! ******************************


!  The gradient of the complex field 'fin' in y-direction.
!  The fields are given in coordinate space. The gradient is
!  evaluated in k_y space.
!

USE params
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN OUT)                      :: fin(kdfull2)
COMPLEX(DP), INTENT(IN OUT)                     :: gradfout(kdfull2)

! ************************************************************

dky=pi/(dy*ny)
DO i3=1,nz2
  DO i1=1,nx2
!                 forward transform along y
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftay(mody(i2)) = fin(ind)
    END DO
#if(netlib_fft)
    CALL dcftf1 (ny2,fftay,wrky,wsavey,ifacy)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pforwy,fftay,fftay)
#endif
#if(fftw_gpu)
    CALL run_fft_for(py,fftay,fftay,ny2)
#endif
!                 multiply by k_y factor in k_y space
    DO i2=1,ny2
      IF(i2 >= (ny+1)) THEN
        zky=(i2-ny2-1)*dky
      ELSE
        zky=(i2-1)*dky
      END IF
      fftay(i2) = fftay(i2)*eye*zky
    END DO
!                 backward transform along y
#if(netlib_fft)
    CALL dcftb1 (ny2,fftay,wrky,wsavey,ifacy)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pbacky,fftay,fftay)
#endif
#if(fftw_gpu)
    CALL run_fft_back(py,fftay,fftay,ny2)
#endif
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      gradfout(ind)= fftay(mody(i2))/ny2
    END DO
  END DO
END DO

RETURN

END SUBROUTINE  ygradient_rspace


! ******************************

SUBROUTINE  zgradient_rspace(fin,gradfout)

! ******************************


!  The gradient of the complex field 'fin' in z-direction.
!  The fields are given in coordinate space. The gradient is
!  evaluated in k_z space.
!

USE params
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN OUT)                      :: fin(kdfull2)
COMPLEX(DP), INTENT(IN OUT)                     :: gradfout(kdfull2)
#if(netlib_fft)
COMPLEX(DP), ALLOCATABLE :: fftaz(:)
#endif



! ************************************************************

#if(netlib_fft)
ALLOCATE(fftaz(nz2))
#endif

dkz=pi/(dz*nz)
DO i2=1,ny2
  DO i1=1,nx2
!                 forward transform along z
    DO i3=1,nz2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      i3m = modz(i3)
      fftaz(i3m) = fin(ind)
    END DO
#if(netlib_fft)
    CALL dcftf1 (nz2,fftaz,wrkz,wsavez,ifacz)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pforwz,fftaz,fftaz)
#endif
#if(fftw_gpu)
    CALL run_fft_for(pz,fftb(1,i1),fftb(1,i1),nz2)
#endif
!                 multiply by k_z factor in k_z space
    DO i3=1,nz2
      IF(i3 >= (nz+1)) THEN
        zkz=(i3-nz2-1)*dkz
      ELSE
        zkz=(i3-1)*dkz
      END IF
      fftaz(i3) = fftaz(i3)*eye*zkz
    END DO
!                 backward transform along z
#if(netlib_fft)
    CALL dcftb1 (nz2,fftaz,wrkz,wsavez,ifacz)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pbackz,fftaz,fftaz)
#endif
#if(fftw_gpu)
    CALL run_fft_back(pz,fftb(1,i1),fftb(1,i1),nz2)
#endif
    DO i3=1,nz2                  ! copy back
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      gradfout(ind)= fftaz(modz(i3))/nz2
    END DO
  END DO
END DO

#if(netlib_fft)
DEALLOCATE(fftaz)
#endif

RETURN

END SUBROUTINE  zgradient_rspace

! ******************************
#if(netlib_fft|fftw_cpu)
SUBROUTINE  fftf(q1,q2)
#endif
#if(fftw_gpu)
SUBROUTINE  fftf(q1,q2,fft,gpu_fft,copyback)
#endif

! ******************************

USE params
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN)                      :: q1(kdfull2)
COMPLEX(DP), INTENT(OUT)                     :: q2(kdfull2)
#if(fftw_gpu)
COMPLEX(C_DOUBLE_COMPLEX)                    :: fft(nx2,ny2,nz2)
COMPLEX(C_DOUBLE_COMPLEX)                    :: gpu_fft(kdfull2)
LOGICAL, INTENT(IN)                          :: copyback
INTEGER :: typefft=1
#endif
!INTEGER, PARAMETER :: kfft=2*kxmax
!INTEGER, PARAMETER :: kfft2=kfft*2+1
!COMPLEX(DP) :: fftax,fftay,fftb
!COMMON /mfftini/fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax),  &
!    wrkx(kfft2),wrky(kfft2),wrkz(kfft2),  &
!    wsavex(kfft2),wsavey(kfft2),wsavez(kfft2),  &
!    ifacx(kfft2),ifacy(kfft2),ifacz(kfft2)

#if(netlib_fft)
INTEGER,SAVE :: nxini=0,nyini=0,nzini=0     ! flag for initialization
#endif

tnorm=1D0/SQRT(8D0*pi*pi*pi*REAL(nx2*ny2*nz2,DP))

#if(netlib_fft)
!     check initialization

! nx=nx2/2
! ny=ny2/2
! nz=nz2/2
IF(nxini == 0) THEN
  CALL dcfti1 (nx2,wsavex,ifacx)
  nxini  = nx2
!        write(6,'(a)') ' x-fft initialized '
ELSE IF(nxini /= nx2) THEN
  STOP ' nx2 in four3d not as initialized!'
END IF
IF(nyini == 0) THEN
  CALL dcfti1 (ny2,wsavey,ifacy)
  nyini  = ny2
!        write(6,'(a)') ' y-fft initialized '
ELSE IF(nyini /= ny2) THEN
  STOP ' ny2 in four3d not as initialized!'
END IF
IF(nzini == 0) THEN
  CALL dcfti1 (nz2,wsavez,ifacz)
  nzini  = nz2
!        write(6,'(a)') ' z-fft initialized '
ELSE IF(nzini /= nz2) THEN
  STOP ' nz2 in four3d not as initialized!'
END IF

!     transformation in x-direction

! nxyf=nx2*ny2
!       nyf=nx2
DO i3=1,nz2
  DO i2=1,ny2
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftax(modx(i1))=q1(ind)  ! copy to workspace
    END DO
    CALL dcftf1 (nx2,fftax,wrkx,wsavex,ifacx)    ! basic fft
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q2(ind)= fftax(i1)        !  copy back in strange order
    END DO
  END DO
END DO

!     transformation in y-direction

DO i3=1,nz2
  DO i1=1,nx2
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftay(mody(i2)) = q2(ind)
    END DO
    CALL dcftf1 (ny2,fftay,wrky,wsavey,ifacy)
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q2(ind)= fftay(i2)
    END DO
  END DO
END DO

!oldc
!oldc     transformation in z-direction
!oldc
!old      do i2=1,ny2
!old        do i1=1,nx2
!old          do i3=1,nz2
!old            ind=(i3-1)*nxyf+(i2-1)*nyf+i1
!old            ffta(mod(i3+nz,nz2)+1) = q2(ind)
!old          enddo
!old          call dcftf1 (nz2,ffta,wrkz,wsavez,ifacz)
!old          do i3=1,nz2
!old            ind=(i3-1)*nxyf+(i2-1)*nyf+i1
!old            q2(ind)= ffta(i3)*tnorm
!old          enddo
!old        enddo
!old      enddo

!     transformation in z-direction

DO i2=1,ny2
  DO i3=1,nz2
    i3m = modz(i3)
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftb(i3m,i1) = q2(ind)
    END DO
  END DO
  DO i1=1,nx2
    CALL dcftf1 (nz2,fftb(1,i1),wrkz,wsavez,ifacz)
  END DO
  DO i3=1,nz2
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q2(ind)= fftb(i3,i1)*tnorm
    END DO
  END DO
END DO
#endif

#if(fftw_cpu)
CALL copy1dto3d(q1,ffta,nx2,ny2,nz2)

CALL fftw_execute_dft(pforw,ffta,ffta)

CALL copy3dto1d(ffta,q2,tnorm,nx2,ny2,nz2)
#endif

#if(fftw_gpu)
CALL copy1dto3d(q1,fft,nx2,ny2,nz2)

CALL copy_on_gpu(fft,gpu_fft,kdfull2)

CALL run_fft_for3d(p,gpu_fft,typefft)

CALL multiply_gpu(gpu_fft,kdfull2,tnorm)

IF(copyback) THEN
  CALL copy_from_gpu(fft,gpu_fft,kdfull2)

  CALL copy3dto1d(fft,q2,nx2,ny2,nz2)
END IF
#endif

RETURN
END SUBROUTINE  fftf

! ******************************
#if(netlib_fft|fftw_cpu)
SUBROUTINE  fftback(q1,q2)
#endif
#if(fftw_gpu)
SUBROUTINE  fftback(q1,q2,fft,gpu_fft)
#endif
! ******************************

USE params
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN)                      :: q1(kdfull2)
COMPLEX(DP), INTENT(OUT)                     :: q2(kdfull2)
#if(fftw_gpu)
COMPLEX(C_DOUBLE_COMPLEX)                    :: fft(nx2,ny2,nz2)
COMPLEX(C_DOUBLE_COMPLEX)                    :: gpu_fft(kdfull2)
!LOGICAL, INTENT(IN)                          :: recopy
INTEGER                                      :: typefft=1
#endif
!INTEGER, PARAMETER :: kfft=2*kxmax
!INTEGER, PARAMETER :: kfft2=kfft*2+1
!COMPLEX(DP) :: fftax,fftay,fftb
!COMMON /mfftini/fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax),  &
!    wrkx(kfft2),wrky(kfft2),wrkz(kfft2),  &
!    wsavex(kfft2),wsavey(kfft2),wsavez(kfft2),  &
!    ifacx(kfft2),ifacy(kfft2),ifacz(kfft2)

! nxyf=nx2*ny2
!      nyf=nx2
!      nx=nx2/2
! ny=ny2/2
! nz=nz2/2

facnr =SQRT(8D0*pi*pi*pi)/SQRT(REAL(nx2*ny2*nz2,DP))

#if(netlib_fft)
!     transformation in z-direction

DO i2=1,ny2
  DO i3=1,nz2
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftb(i3,i1) = q1(ind)*facnr
    END DO
  END DO
  DO i1=1,nx2
    CALL dcftb1 (nz2,fftb(1,i1),wrkz,wsavez,ifacz)    ! basic fft
  END DO
  DO i3=1,nz2                  ! copy back
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q2(ind)= fftb(modz(i3),i1)
    END DO
  END DO
END DO

!     transformation in y-direction

DO i3=1,nz2
  DO i1=1,nx2
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftay(i2) = q2(ind)
    END DO
    CALL dcftb1 (ny2,fftay,wrky,wsavey,ifacy)
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q2(ind)= fftay(mody(i2))
    END DO
  END DO
END DO


!     transformation in x-direction

DO i3=1,nz2
  DO i2=1,ny2
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftax(i1)=q2(ind)
    END DO
    CALL dcftb1 (nx2,fftax,wrkx,wsavex,ifacx)
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q2(ind)= fftax(modx(i1))
    END DO
  END DO
END DO
#endif

#if(fftw_cpu)
CALL secopy1dto3d(q1,ffta,nx2,ny2,nz2)

CALL fftw_execute_dft(pback,ffta,ffta)

CALL secopy3dto1d(ffta,q2,facnr,nx2,ny2,nz2)
#endif

#if(fftw_gpu)
CALL run_fft_back3d(p,gpu_fft,typefft)

CALL multiply_gpu(gpu_fft,kdfull2,facnr)

CALL copy_from_gpu(fft,gpu_fft,kdfull2)

CALL secopy3dto1d(fft,q2,nx2,ny2,nz2)
#endif

RETURN
END SUBROUTINE  fftback


! 1 "fft.f"

! ******************************
#if(netlib_fft|fftw_cpu)
SUBROUTINE  rftf(q1,q2)
#endif
#if(fftw_gpu)
SUBROUTINE  rftf(q1,q2,fft,gpu_fft,copyback)
#endif
! ******************************

USE params
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: q1(kdfull2)
COMPLEX(DP), INTENT(OUT)                     :: q2(kdfull2)
!INTEGER, PARAMETER :: kfft=2*kxmax
!INTEGER, PARAMETER :: kfft2=kfft*2+1
!COMPLEX(DP) :: fftax,fftay,fftb
!COMMON /mfftini/fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax),  &
!    wrkx(kfft2),wrky(kfft2),wrkz(kfft2),  &
!    wsavex(kfft2),wsavey(kfft2),wsavez(kfft2),  &
!    ifacx(kfft2),ifacy(kfft2),ifacz(kfft2)

#if(netlib_fft)
INTEGER,SAVE :: nxini=0,nyini=0,nzini=0     ! flag for initialization
#endif
#if(fftw_gpu)
COMPLEX(C_DOUBLE_COMPLEX)                    :: fft(nx2,ny2,nz2)
COMPLEX(C_DOUBLE_COMPLEX)                    :: gpu_fft(kdfull2)
LOGICAL, INTENT(IN)                          :: copyback
INTEGER :: typefft=2
#endif

tnorm=1D0/SQRT(8D0*pi*pi*pi*REAL(nx2*ny2*nz2,DP))

#if(netlib_fft)
!     check initialization
! nx=nx2/2
! ny=ny2/2
! nz=nz2/2

IF(nxini == 0) THEN
  CALL dcfti1 (nx2,wsavex,ifacx)
  nxini  = nx2
!        write(6,'(a)') ' x-fft initialized '
ELSE IF(nxini /= nx2) THEN
  STOP ' nx2 in four3d not as initialized!'
END IF
IF(nyini == 0) THEN
  CALL dcfti1 (ny2,wsavey,ifacy)
  nyini  = ny2
!        write(6,'(a)') ' y-fft initialized '
ELSE IF(nyini /= ny2) THEN
  STOP ' ny2 in four3d not as initialized!'
END IF
IF(nzini == 0) THEN
  CALL dcfti1 (nz2,wsavez,ifacz)
  nzini  = nz2
!        write(6,'(a)') ' z-fft initialized '
ELSE IF(nzini /= nz2) THEN
  STOP ' nz2 in four3d not as initialized!'
END IF

!     transformation in x-direction

! nxyf=nx2*ny2
!      nyf=nx2
DO i3=1,nz2
  DO i2=1,ny2
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftax(modx(i1))=CMPLX(q1(ind),0D0,DP) ! copy to workspace
    END DO
    CALL dcftf1 (nx2,fftax,wrkx,wsavex,ifacx)    ! basic fft
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q2(ind)= fftax(i1)        !  copy back in strange order
    END DO
  END DO
END DO

!     transformation in y-direction

DO i3=1,nz2
  DO i1=1,nx2
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftay(mody(i2)) = q2(ind)
    END DO
    CALL dcftf1 (ny2,fftay,wrky,wsavey,ifacy)
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q2(ind)= fftay(i2)
    END DO
  END DO
END DO

!oldc
!oldc     transformation in z-direction
!oldc
!old      do i2=1,ny2
!old        do i3=1,nz2
!old          do i1=1,nx2
!old            ind=(i3-1)*nxyf+(i2-1)*nyf+i1
!old            fftb(mod(i3+nz,nz2)+1,i1) = q2(ind)
!old          enddo
!old        enddo
!old        do i1=1,nx2
!old          call dcftf1 (nz2,fftb(1,i1),wrkz,wsavez,ifacz)
!old        enddo
!old        do i3=1,nz2
!old          do i1=1,nx2
!old            ind=(i3-1)*nxyf+(i2-1)*nyf+i1
!old            q2(ind)= fftb(i3,i1)*tnorm
!old          enddo
!old        enddo
!old      enddo

!     transformation in z-direction

DO i2=1,ny2
  DO i3=1,nz2
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftb(modz(i3),i1) = q2(ind)
    END DO
  END DO
  DO i1=1,nx2
    CALL dcftf1 (nz2,fftb(1,i1),wrkz,wsavez,ifacz)
  END DO
  DO i3=1,nz2
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q2(ind)= fftb(i3,i1)*tnorm
    END DO
  END DO
END DO
#endif

#if(fftw_cpu)
CALL copyr1dto3d(q1,ffta,nx2,ny2,nz2)

CALL fftw_execute_dft(pforw,ffta,ffta)

CALL copy3dto1d(ffta,q2,tnorm,nx2,ny2,nz2)
#endif

#if(fftw_gpu)
CALL copyr1dto3d(q1,fft,nx2,ny2,nz2)

CALL copy_on_gpu(fft,gpu_fft,kdfull2)

CALL run_fft_for3d(p,gpu_fft,typefft)

CALL multiply_gpu(gpu_fft,kdfull2,tnorm)

IF(copyback) THEN
  CALL copy_from_gpu(fft,gpu_fft,kdfull2)

  CALL copy3dto1d(fft,q2,nx2,ny2,nz2)
ENDIF
#endif

RETURN
END SUBROUTINE  rftf

! ******************************

#if(netlib_fft|fftw_cpu)
SUBROUTINE  rfftback(q1,q3)
!SUBROUTINE  rfftback(q1,q2)
#endif
#if(fftw_gpu)
SUBROUTINE  rfftback(q1,q3,fft,gpu_fft,recopy)
#endif

! ******************************

USE params
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN)                      :: q1(kdfull2)
!REAL(DP), INTENT(OUT)                        :: q2(kdfull2)
REAL(DP), INTENT(OUT)                        :: q3(kdfull2)
!INTEGER, PARAMETER :: kfft=2*kxmax
!INTEGER, PARAMETER :: kfft2=kfft*2+1
!COMPLEX(DP) :: fftax,fftay,fftb
!COMMON /mfftini/fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax),  &
!    wrkx(kfft2),wrky(kfft2),wrkz(kfft2),  &
!    wsavex(kfft2),wsavey(kfft2),wsavez(kfft2),  &
!    ifacx(kfft2),ifacy(kfft2),ifacz(kfft2)

#if(netlib_fft)
COMPLEX(DP),ALLOCATABLE :: q2(:)
#endif
#if(fftw_gpu)
COMPLEX(C_DOUBLE_COMPLEX)                    :: fft(nx2,ny2,nz2)
COMPLEX(C_DOUBLE_COMPLEX)                    :: gpu_fft(kdfull2)
LOGICAL,INTENT(IN) :: recopy
INTEGER :: typefft=2
#endif
!      data  nxini,nyini,nzini/0,0,0/  ! flag for initialization
! nxyf=nx2*ny2
!      nyf=nx2
!      nx=nx2/2
! ny=ny2/2
! nz=nz2/2

facnr =SQRT(8D0*pi*pi*pi)/SQRT(REAL(nx2*ny2*nz2,DP))

#if(netlib_fft)
!oldc
!oldc     transformation in z-direction
!oldc
!old      do i2=1,ny2
!old        do i1=1,nx2
!old          do i3=1,nz2
!old            ind=(i3-1)*nxyf+(i2-1)*nyf+i1
!old            ffta(i3) = q1(ind)*facnr
!old          enddo
!old          call dcftb1 (nz2,ffta,wrkz,wsavez,ifacz)    ! basic fft
!old          do i3=1,nz2                  ! copy back
!old            ind=(i3-1)*nxyf+(i2-1)*nyf+i1
!old            q2(ind)= ffta(mod(i3+nz,nz2)+1)
!old          enddo
!old        enddo
!old      enddo

ALLOCATE(q2(kdfull2))

!     transformation in z-direction

DO i2=1,ny2
  DO i3=1,nz2
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftb(i3,i1) = q1(ind)*facnr
    END DO
  END DO
  DO i1=1,nx2
    CALL dcftb1 (nz2,fftb(1,i1),wrkz,wsavez,ifacz)    ! basic fft
  END DO
  DO i3=1,nz2                  ! copy back
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q2(ind)= fftb(modz(i3),i1)
    END DO
  END DO
END DO

!     transformation in y-direction

DO i3=1,nz2
  DO i1=1,nx2
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftay(i2) = q2(ind)
    END DO
    CALL dcftb1 (ny2,fftay,wrky,wsavey,ifacy)
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q2(ind)= fftay(mody(i2))
    END DO
  END DO
END DO


!     transformation in x-direction

DO i3=1,nz2
  DO i2=1,ny2
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftax(i1)=q2(ind)
    END DO
    CALL dcftb1 (nx2,fftax,wrkx,wsavex,ifacx)
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
!      q2(ind)= REAL(fftax(MOD(i1+nx,nx2)+1),DP)
      q3(ind)= REAL(fftax(modx(i1)),DP)
    END DO
  END DO
END DO

DEALLOCATE(q2)
#endif

#if(fftw_cpu)
CALL secopy1dto3d(q1,ffta,nx2,ny2,nz2)

CALL fftw_execute_dft(pback,ffta,ffta)

CALL copyr3dto1d(ffta,q3,facnr,nx2,ny2,nz2)
#endif

#if(fftw_gpu)
IF(recopy) THEN
  CALL secopy1dto3d(q1,fft,nx2,ny2,nz2)

  CALL copy_on_gpu(fft,gpu_fft,kdfull2)
ENDIF

CALL run_fft_back3d(p,gpu_fft,typefft)

CALL multiply_gpu(gpu_fft,kdfull2,facnr)

CALL copy_from_gpu(fft,gpu_fft,kdfull2)

CALL copyr3dto1d(fft,q3,nx2,ny2,nz2)
#endif

RETURN
END SUBROUTINE  rfftback

#if(fftw_cpu|fftw_gpu)

SUBROUTINE copy1dto3d(vec1d,vec3d,nbx2,nby2,nbz2)

! ******************************

USE params

COMPLEX(DP), INTENT(IN)                      :: vec1d(kdfull2)
COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT)       :: vec3d(nbx2,nby2,nbz2)

DO i3=1,nbz2
  DO i2=1,nby2
    DO i1=1,nbx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      vec3d(modx(i1),mody(i2),modz(i3))=vec1d(ind)
    END DO
  END DO
END DO

RETURN
END SUBROUTINE copy1dto3d

! ******************************

SUBROUTINE copyr1dto3d(vec1d,vec3d,nbx2,nby2,nbz2)

! ******************************

USE params

REAL(DP), INTENT(IN)                      :: vec1d(kdfull2)
COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT)    :: vec3d(nbx2,nby2,nbz2)
       
DO i3=1,nbz2
  DO i2=1,nby2
    DO i1=1,nbx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      vec3d(modx(i1),mody(i2),modz(i3))=CMPLX(vec1d(ind),0D0,DP) 
    END DO
  END DO
END DO

RETURN
END SUBROUTINE copyr1dto3d

! ******************************

SUBROUTINE secopy1dto3d(vec1d,vec3d,nbx2,nby2,nbz2)

! ******************************

USE params

COMPLEX(DP), INTENT(IN)                      :: vec1d(kdfull2)
COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT)       :: vec3d(nbx2,nby2,nbz2)

DO i3=1,nbz2
  DO i2=1,nby2
    DO i1=1,nbx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      vec3d(i1,i2,i3)=vec1d(ind)
    END DO
  END DO
END DO

RETURN
END SUBROUTINE secopy1dto3d

! ******************************
#if(fftw_gpu)
SUBROUTINE copyr1dtor3d(vec1d,vec3d,nbx2,nby2,nbz2)

! ******************************

USE params

REAL(DP), INTENT(IN)                      :: vec1d(kdfull2)
REAL(C_DOUBLE), INTENT(OUT)               :: vec3d(nbx2,nby2,nbz2)
       
DO i3=1,nbz2
  DO i2=1,nby2
    DO i1=1,nbx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      vec3d(modx(i1),mody(i2),modz(i3))=vec1d(ind)
    END DO
  END DO
END DO

RETURN
END SUBROUTINE copyr1dtor3d

! ******************************
#endif
#if(fftw_cpu)
SUBROUTINE copy3dto1d(vec3d,vec1d,coef,nbx2,nby2,nbz2)
#else
SUBROUTINE copy3dto1d(vec3d,vec1d,nbx2,nby2,nbz2)
#endif

! ******************************

USE params

COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN)        :: vec3d(nbx2,nby2,nbz2)
COMPLEX(DP), INTENT(OUT)                     :: vec1d(kdfull2)

DO i3=1,nbz2
  DO i2=1,nby2
    DO i1=1,nbx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
#if(fftw_cpu)
      vec1d(ind)=coef*ffta(i1,i2,i3)
#else
      vec1d(ind)=vec3d(i1,i2,i3)
#endif
    END DO
  END DO
END DO

RETURN
END SUBROUTINE copy3dto1d

! ******************************

#if(fftw_cpu)
SUBROUTINE copyr3dto1d(vec3d,vec1d,coef,nbx2,nby2,nbz2)
#else
SUBROUTINE copyr3dto1d(vec3d,vec1d,nbx2,nby2,nbz2)
#endif
! ******************************

USE params

COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN)     :: vec3d(nbx2,nby2,nbz2)
REAL(DP), INTENT(OUT)                     :: vec1d(kdfull2)

DO i3=1,nbz2
  DO i2=1,nby2
    DO i1=1,nbx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
#if(fftw_cpu)
      vec1d(ind)=REAL(vec3d(modx(i1),mody(i2),modz(i3)),DP)*coef
#else
      vec1d(ind)=REAL(vec3d(modx(i1),mody(i2),modz(i3)),DP)
#endif
    END DO
  END DO
END DO

RETURN
END SUBROUTINE copyr3dto1d

! ******************************

#if(fftw_cpu)
SUBROUTINE secopy3dto1d(vec3d,vec1d,coef,nbx2,nby2,nbz2)
#else
SUBROUTINE secopy3dto1d(vec3d,vec1d,nbx2,nby2,nbz2)
#endif

! ******************************

USE params

COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN)        :: vec3d(nbx2,nby2,nbz2)
COMPLEX(DP), INTENT(OUT)                     :: vec1d(kdfull2)

DO i3=1,nbz2
  DO i2=1,nby2
    DO i1=1,nbx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
#if(fftw_cpu)
      vec1d(ind)=vec3d(modx(i1),mody(i2),modz(i3))*coef
#else
      vec1d(ind)=vec3d(modx(i1),mody(i2),modz(i3))
#endif
    END DO
  END DO
END DO

RETURN
END SUBROUTINE secopy3dto1d

! ******************************

SUBROUTINE fft_end()

! ******************************

#if(fftw_cpu)
USE FFTW

CALL fftw_destroy_plan(pforw)
CALL fftw_destroy_plan(pback)
CALL fftw_destroy_plan(pforwx)
CALL fftw_destroy_plan(pforwy)
CALL fftw_destroy_plan(pforwz)
CALL fftw_destroy_plan(pforwz1)
CALL fftw_destroy_plan(pbackx)
CALL fftw_destroy_plan(pbacky)
CALL fftw_destroy_plan(pbackz)
CALL fftw_destroy_plan(pbackz1)

DEALLOCATE(fftax,fftay,fftaz,fftb,ffta)
DEALLOCATE(ak,akv)
#endif

#if(fftw_gpu)
res = cudaFreeHost(c_p_fftax)
res = cudaFreeHost(c_p_fftay)
res = cudaFreeHost(c_p_fftaz)
res = cudaFreeHost(c_p_fftb)
res = cudaFreeHost(c_p_ffta)
res = cudaFree(c_gpu_ffta)
res = cudaFree(c_gpu_akfft)
res = cudaFree(c_gpu_akvfft)
res = cudaFree(c_gpu_akxfft)
res = cudaFree(c_gpu_akyfft)
res = cudaFree(c_gpu_akzfft)

CALL kill_plan(px)
CALL kill_plan(py)
CALL kill_plan(pz)
CALL kill_plan(p)
#endif

DEALLOCATE(akpropx,akpropy,akpropz)

RETURN
END SUBROUTINE fft_end
#endif

#if(fftw_gpu)
! ******************************

SUBROUTINE my_cuda_allocate(nbx2,nby2,nbz2)

! ******************************

res = cudaMallocHost(c_p_fftax,nbx2*sizeof(size_cmplx))
CALL c_f_pointer(c_p_fftax,fftax,[nbx2])

res = cudaMallocHost(c_p_fftay,nby2*sizeof(size_cmplx))
CALL c_f_pointer(c_p_fftay,fftay,[nby2])

res = cudaMallocHost(c_p_fftaz,nbz2*sizeof(size_cmplx))
CALL c_f_pointer(c_p_fftaz,fftaz,[nbz2])

res = cudaMallocHost(c_p_fftb,nbz2*nbx2*sizeof(size_cmplx))
CALL c_f_pointer(c_p_fftb,fftb,[nbz2,nbx2])

res = cudaMallocHost(c_p_ffta,kdfull2*sizeof(size_cmplx))
CALL c_f_pointer(c_p_ffta,ffta,[nbx2,nby2,nbz2])

res = cudaMallocHost(c_p_ffta2,kdfull2*sizeof(size_cmplx))
CALL c_f_pointer(c_p_ffta2,ffta2,[nbx2,nby2,nbz2])

res = cudaMalloc(c_gpu_ffta,kdfull2*sizeof(size_cmplx))
CALL c_f_pointer(c_gpu_ffta,gpu_ffta,[kdfull2])

res = cudaMalloc(c_gpu_ffta2,kdfull2*sizeof(size_cmplx))
CALL c_f_pointer(c_gpu_ffta2,gpu_ffta2,[kdfull2])

res = cudaMalloc(c_gpu_ffta_int,kdfull2*sizeof(size_cmplx))
CALL c_f_pointer(c_gpu_ffta_int,gpu_ffta_int,[kdfull2])

res = cudaMalloc(c_gpu_akfft,kdfull2*sizeof(size_cmplx))
CALL c_f_pointer(c_gpu_akfft,gpu_akfft,[kdfull2])

res = cudaMalloc(c_gpu_akvfft,kdfull2*sizeof(size_double))
CALL c_f_pointer(c_gpu_akvfft,gpu_akvfft,[kdfull2])

res = cudaMalloc(c_gpu_akxfft,kdfull2*sizeof(size_cmplx))
CALL c_f_pointer(c_gpu_akxfft,gpu_akxfft,[kdfull2])

res = cudaMalloc(c_gpu_akyfft,kdfull2*sizeof(size_cmplx))
CALL c_f_pointer(c_gpu_akyfft,gpu_akyfft,[kdfull2])

res = cudaMalloc(c_gpu_akzfft,kdfull2*sizeof(size_cmplx))
CALL c_f_pointer(c_gpu_akzfft,gpu_akzfft,[kdfull2])

END SUBROUTINE my_cuda_allocate
#endif

END MODULE kinetic
