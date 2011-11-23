MODULE kinetic
USE params, ONLY: DP
IMPLICIT REAL(DP) (A-H,O-Z)

SAVE
!     arrays for kinetic energy and electronic propagation
!       kinetic energy coefficients in strange ordered fourier space
!     akv  = fouier-field for 0.5*k^2
!     ak   = fourier-field for exp(i*dt*(h^2/2m)*k^2)
COMPLEX(DP),ALLOCATABLE :: ak(:)
REAL(DP),ALLOCATABLE :: akv(:)
COMPLEX(DP),ALLOCATABLE :: akpropx(:),akpropy(:),akpropz(:)
COMPLEX(DP),PARAMETER,PRIVATE :: eye=(0D0,1D0)
REAL(DP),PARAMETER,PRIVATE :: PI=3.141592653589793D0


INTEGER, PRIVATE :: kfft,kfft2,kdfull2

COMPLEX(DP), PRIVATE, ALLOCATABLE :: fftax(:),fftay(:),fftb(:,:)
REAL(DP), PRIVATE, ALLOCATABLE :: wrkx(:),wrky(:),wrkz(:)
REAL(DP), PRIVATE, ALLOCATABLE :: wsavex(:),wsavey(:),wsavez(:)
INTEGER, PRIVATE, ALLOCATABLE :: ifacx(:),ifacy(:),ifacz(:)

CONTAINS
!-----init_grid_fft-----------------------------------------------------

SUBROUTINE init_grid_fft(dx0,dy0,dz0,nx0,ny0,nz0,dt1,h2m)

!     initialize details for FFT

REAL(DP) :: dt1,h2m

nx2=nx0;ny2=ny0;nz2=nz0
kxmax=nx0;kymax=ny0;kzmax=nz0
nx=nx2/2;ny=ny2/2;nz=nz2/2

kfft=2*kxmax
kfft2=kfft*2+1
kdfull2=kxmax*kymax*kzmax

dkx=pi/(dx0*nx)
dky=pi/(dy0*ny)
dkz=pi/(dz0*nz)


ALLOCATE(ak(kdfull2),akv(kdfull2))
ALLOCATE(akpropx(kxmax),akpropy(kymax),akpropz(kzmax))
ALLOCATE(fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax))
ALLOCATE(wrkx(kfft2),wrky(kfft2),wrkz(kfft2))
ALLOCATE(wsavex(kfft2),wsavey(kfft2),wsavez(kfft2))
ALLOCATE(ifacx(kfft2),ifacy(kfft2),ifacz(kfft2))

WRITE(7,*) 'h bar squared over two m electron',h2m
WRITE(7,*) ' testprint EYE=',eye
WRITE(7,*) ' dkx,dky,dkz=',dkx,dky,dkz
WRITE(7,*) ' testprint: nx2,ny2,nz2=',nx2,ny2,nz2
WRITE(7,*) ' testprint: kdfull2,kfft2=',kdfull2,kfft2

!     prepare k**2 and kinetic propagation factor in 3D momentum space

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
    END DO
  END DO
END DO

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

END SUBROUTINE init_grid_fft



! ******************************

SUBROUTINE  kinprop(q1,q2)

! ******************************

!       propagation with exp(-i*dt*e_kin)

USE params
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN OUT)                  :: q1(kdfull2)
COMPLEX(DP), INTENT(OUT)                     :: q2(kdfull2)

DATA  nxini,nyini,nzini/0,0,0/ ! flag for initialization

!       check initialization

tnorm=1D0/SQRT(8D0*pi*pi*pi*REAL(nx2*ny2*nz2,DP))
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

!       propagation in x-direction

xfnorm = 1D0/nx2
DO i3=1,nz2
  DO i2=1,ny2
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftax(MOD(i1+nx,nx2)+1)=q1(ind) ! copy to workspace
    END DO
    CALL dcftf1 (nx2,fftax,wrkx,wsavex,ifacx) ! basic fft
    DO i1=1,nx2
      fftax(i1) = akpropx(i1)*fftax(i1)
    END DO
    CALL dcftb1 (nx2,fftax,wrkx,wsavex,ifacx) ! basic back fft
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q2(ind)= fftax(MOD(i1+nx,nx2)+1)*xfnorm
    END DO
  END DO
END DO

!      transformation in y-direction

yfnorm = 1D0/ny2
DO i3=1,nz2
  DO i1=1,nx2
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftay(MOD(i2+ny,ny2)+1) = q2(ind)
    END DO
    CALL dcftf1 (ny2,fftay,wrky,wsavey,ifacy)
    DO i2=1,ny2
      fftay(i2) = akpropy(i2)*fftay(i2)
    END DO
    CALL dcftb1 (ny2,fftay,wrky,wsavey,ifacy)
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q2(ind)= fftay(MOD(i2+ny,ny2)+1)*yfnorm
    END DO
  END DO
END DO

!       propagation in z-direction

zfnorm = 1D0/nz2
DO i2=1,ny2
  DO i3=1,nz2
    i3m = MOD(i3+nz,nz2)+1
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftb(i3m,i1) = q2(ind)
    END DO
  END DO
  DO i1=1,nx2
    CALL dcftf1 (nz2,fftb(1,i1),wrkz,wsavez,ifacz)
    DO i3=1,nz2
      fftb(i3,i1) = akpropz(i3)*fftb(i3,i1)
    END DO
    CALL dcftb1 (nz2,fftb(1,i1),wrkz,wsavez,ifacz)
  END DO
  DO i3=1,nz2
    i3m = MOD(i3+nz,nz2)+1
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


  IF(idirec == 1) THEN
!       x-derivative
   
    dkx=pi/(dx*REAL(nx,DP))
    ind=0
    DO i3=1,nz2
      DO i2=1,ny2
        DO i1=1,nx2
          IF(i1 >= (nx+1)) THEN
            zkx=(i1-nx2-1)*dkx
          ELSE
            zkx=(i1-1)*dkx
          END IF
          ind=ind+1
          gradfout(ind) = eye*zkx*fin(ind)
        END DO
      END DO
    END DO

  ELSEIF(idirec == 2) THEN    
!       y-derivative
    
    dky=pi/(dy*REAL(ny,DP))
    ind=0
    DO i3=1,nz2
      DO i2=1,ny2
        IF(i2 >= (ny+1)) THEN
          zky=(i2-ny2-1)*dky
        ELSE
          zky=(i2-1)*dky
        END IF
        DO i1=1,nx2
          ind=ind+1
          gradfout(ind) = eye*zky*fin(ind)
        END DO
      END DO
    END DO

  ELSEIF(idirec == 3) THEN          
!       z-derivative
    
    ind=0
    dkz=pi/(dz*REAL(nz,DP))
    DO i3=1,nz2
      IF(i3 >= (nz+1)) THEN
        zkz=(i3-nz2-1)*dkz
      ELSE
        zkz=(i3-1)*dkz
      END IF
      DO i2=1,ny2
        DO i1=1,nx2
          ind=ind+1
          gradfout(ind) = eye*zkz*fin(ind)
        END DO
      END DO
    END DO
    
  ELSE
    STOP ' RGRADIENT called with invalid IDIREC'
  ENDIF

RETURN
END SUBROUTINE  gradient





! ******************************

SUBROUTINE  fftf(q1,q2)

! ******************************

USE params
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN)                      :: q1(kdfull2)
COMPLEX(DP), INTENT(OUT)                     :: q2(kdfull2)
!INTEGER, PARAMETER :: kfft=2*kxmax
!INTEGER, PARAMETER :: kfft2=kfft*2+1
!COMPLEX(DP) :: fftax,fftay,fftb
!COMMON /mfftini/fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax),  &
!    wrkx(kfft2),wrky(kfft2),wrkz(kfft2),  &
!    wsavex(kfft2),wsavey(kfft2),wsavez(kfft2),  &
!    ifacx(kfft2),ifacy(kfft2),ifacz(kfft2)

INTEGER,SAVE :: nxini=0,nyini=0,nzini=0     ! flag for initialization


!     check initialization

! nx=nx2/2
! ny=ny2/2
! nz=nz2/2
tnorm=1D0/SQRT(8D0*pi*pi*pi*REAL(nx2*ny2*nz2,DP))
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
      fftax(MOD(i1+nx,nx2)+1)=q1(ind)  ! copy to workspace
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
      fftay(MOD(i2+ny,ny2)+1) = q2(ind)
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
    i3m = MOD(i3+nz,nz2)+1
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

RETURN
END SUBROUTINE  fftf

! ******************************

SUBROUTINE  fftback(q1,q2)

! ******************************

USE params
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN)                      :: q1(kdfull2)
COMPLEX(DP), INTENT(OUT)                     :: q2(kdfull2)
!INTEGER, PARAMETER :: kfft=2*kxmax
!INTEGER, PARAMETER :: kfft2=kfft*2+1
!COMPLEX(DP) :: fftax,fftay,fftb
!COMMON /mfftini/fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax),  &
!    wrkx(kfft2),wrky(kfft2),wrkz(kfft2),  &
!    wsavex(kfft2),wsavey(kfft2),wsavez(kfft2),  &
!    ifacx(kfft2),ifacy(kfft2),ifacz(kfft2)

!      data  nxini,nyini,nzini/0,0,0/ ! flag for initialization
! nxyf=nx2*ny2
!      nyf=nx2
!      nx=nx2/2
! ny=ny2/2
! nz=nz2/2
facnr =SQRT(8D0*pi*pi*pi)/SQRT(REAL(nx2*ny2*nz2,DP))

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
      q2(ind)= fftb(MOD(i3+nz,nz2)+1,i1)
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
      q2(ind)= fftay(MOD(i2+ny,ny2)+1)
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
      q2(ind)= fftax(MOD(i1+nx,nx2)+1)
    END DO
  END DO
END DO
RETURN
END SUBROUTINE  fftback


! 1 "fft.f"

! ******************************

SUBROUTINE  rftf(q1,q2)

! ******************************

USE params
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

INTEGER,SAVE :: nxini=0,nyini=0,nzini=0     ! flag for initialization


!     check initialization

! nx=nx2/2
! ny=ny2/2
! nz=nz2/2
tnorm=1D0/SQRT(8D0*pi*pi*pi*REAL(nx2*ny2*nz2,DP))
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
      fftax(MOD(i1+nx,nx2)+1)=CMPLX(q1(ind),0D0,DP) ! copy to workspace
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
      fftay(MOD(i2+ny,ny2)+1) = q2(ind)
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
      fftb(MOD(i3+nz,nz2)+1,i1) = q2(ind)
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

RETURN
END SUBROUTINE  rftf

! ******************************

SUBROUTINE  rfftback(q1,q3)
!SUBROUTINE  rfftback(q1,q2)

! ******************************

USE params
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


COMPLEX(DP),ALLOCATABLE :: q2(:)



!      data  nxini,nyini,nzini/0,0,0/  ! flag for initialization
! nxyf=nx2*ny2
!      nyf=nx2
!      nx=nx2/2
! ny=ny2/2
! nz=nz2/2
facnr =SQRT(8D0*pi*pi*pi)/SQRT(REAL(nx2*ny2*nz2,DP))
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
      q2(ind)= fftb(MOD(i3+nz,nz2)+1,i1)
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
      q2(ind)= fftay(MOD(i2+ny,ny2)+1)
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
      q3(ind)= REAL(fftax(MOD(i1+nx,nx2)+1),DP)
    END DO
  END DO
END DO

DEALLOCATE(q2)


RETURN
END SUBROUTINE  rfftback
!INCLUDE 'fftpack.F90'
!INCLUDE 'fftpack2.F90'


END MODULE kinetic


