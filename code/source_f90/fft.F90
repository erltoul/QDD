MODULE kinetic
#if(fftw_cpu)
USE, intrinsic :: iso_c_binding
#endif
USE params, ONLY: DP,numthr
IMPLICIT REAL(DP) (A-H,O-Z)

SAVE
!     arrays for kinetic energy and electronic propagation
!       kinetic energy coefficients in strange ordered fourier space
!     akv  = fouier-field for 0.5*k^2
!     ak   = fourier-field for exp(i*dt*(h^2/2m)*k^2)
COMPLEX(DP),ALLOCATABLE :: ak(:)
REAL(DP),ALLOCATABLE :: akv(:)
COMPLEX(DP),ALLOCATABLE :: akpropx(:),akpropy(:),akpropz(:),akprop(:,:,:)
INTEGER,PRIVATE,ALLOCATABLE :: modx(:),mody(:),modz(:)
COMPLEX(DP),PARAMETER,PRIVATE :: eye=(0D0,1D0)
REAL(DP),PARAMETER,PRIVATE :: PI=3.141592653589793D0


INTEGER, PRIVATE :: kfft,kfft2,kdfull2
INTEGER, PRIVATE :: kxmax,kymax,kzmax
INTEGER, PRIVATE :: iret


#if(netlib_fft)
COMPLEX(DP), PRIVATE, ALLOCATABLE :: fftax(:),fftay(:),fftb(:,:)
REAL(DP), PRIVATE, ALLOCATABLE :: wrkx(:),wrky(:),wrkz(:)
REAL(DP), PRIVATE, ALLOCATABLE :: wsavex(:),wsavey(:),wsavez(:)
INTEGER, PRIVATE, ALLOCATABLE :: ifacx(:),ifacy(:),ifacz(:)
#endif
INTEGER,PUBLIC,SAVE :: FFTW_planflag
#if(fftw_cpu)
COMPLEX(C_DOUBLE_COMPLEX), PRIVATE, ALLOCATABLE :: fftax(:),fftay(:),fftaz(:),fftb(:,:)
COMPLEX(C_DOUBLE_COMPLEX), PRIVATE, ALLOCATABLE :: ffta(:,:,:,:)
type(C_PTR), PRIVATE :: pforwx,pforwy,pforwz,pforwz1,pbackx,pbacky,pbackz,pbackz1
type(C_PTR), PRIVATE,ALLOCATABLE :: pforw(:),pback(:)
INTEGER(C_INT), PRIVATE :: wisdomtest
#endif
#if(paropenmp)
INTEGER,PRIVATE,SAVE :: nacthr
#endif

CONTAINS
!-----init_grid_fft-----------------------------------------------------

SUBROUTINE init_grid_fft(dx0,dy0,dz0,nx0,ny0,nz0,dt1,h2m)

!     initialize details for FFT

#if(fftw_cpu)
USE FFTW
#if(parayes)
USE params, only : myn,numthr,nthr
INCLUDE 'mpif.h'
REAL(DP) :: is(mpi_status_size)
#endif
INTEGER, SAVE ::  nxini=0,nyini=0,nzini=0,nini=0 ! flag for initialization

#endif

REAL(DP) :: dt1,h2m
INTEGER :: omp_get_num_threads

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
ALLOCATE(akprop(kxmax,kymax,kzmax))
ALLOCATE(modx(kxmax),mody(kymax),modz(kzmax))
#if(netlib_fft)
ALLOCATE(fftax(kxmax),fftay(kymax),fftb(kzmax,kxmax))
ALLOCATE(wrkx(kfft2),wrky(kfft2),wrkz(kfft2))
ALLOCATE(wsavex(kfft2),wsavey(kfft2),wsavez(kfft2))
ALLOCATE(ifacx(kfft2),ifacy(kfft2),ifacz(kfft2))
#endif
#if(fftw_cpu)
WRITE(*,*) ' allocate with: kxmax,kymax,kzmax=',kxmax,kymax,kzmax
ALLOCATE(fftax(kxmax),fftay(kymax),fftaz(kzmax),fftb(kzmax,kxmax))
#if(paropenmp && dynopenmp)
  nacthr = numthr-1
#else
  nacthr = 0
#endif
ALLOCATE(ffta(kxmax,kymax,kzmax,0:nacthr),pforw(0:nacthr),pback(0:nacthr))
WRITE(*,*) ' FFTA allocated with NACTHR=',nacthr
!
!  central setting of FFTW planning expense --> edit here
!
!  FFTW_planflag = FFTW_MEASURE  
  FFTW_planflag = FFTW_EXHAUSTIVE
#endif

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
      akprop(i1,i2,i3)=ak(ind)
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


DO i3=1,nz2;  DO i2=1,ny2;    DO i1=1,nx2
  akprop(i1,i2,i3) = akpropx(i1)*akpropy(i2)*akpropz(i3)
END DO;  END DO;  END DO

#if(fftw_cpu)
#if(parano)
#if(paropenmp && !dynopenmp)
  call dfftw_init_threads(iret)
  WRITE(*,*) ' dfftw_init_threads: iret=',iret
!  numthr = 4
  call dfftw_plan_with_nthreads(numthr)
  WRITE(*,*) ' init fft FFTW threads: iret=',iret,', nr. of threads=',numthr,&
   omp_get_num_threads()
#endif
IF (nini==0) THEN
  wisdomtest=fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
  IF (wisdomtest == 0) THEN
    wisdomtest = fftw_import_system_wisdom()
    IF(wisdomtest == 0) THEN
      WRITE(6,*) 'wisdom_fftw.dat not found, creating it'
      WRITE(7,*) 'wisdom_fftw.dat not found, creating it'
    ELSE
      WRITE(*,*) 'wisdom from system'
    END IF
  END IF
!  pforw=fftw_plan_dft_3d(nz2,ny2,nx2,ffta,ffta,FFTW_FORWARD,FFTW_planflag+FFTW_UNALIGNED)
!      initialitze 3D plans, take care for independence in OpenMP
  DO i=0,nacthr
    pforw(i)=fftw_plan_dft_3d(nz2,ny2,nx2,ffta(:,:,:,i),ffta(:,:,:,i),FFTW_FORWARD,FFTW_planflag)
    pback(i)=fftw_plan_dft_3d(nz2,ny2,nx2,ffta(:,:,:,i),ffta(:,:,:,i),FFTW_BACKWARD,FFTW_planflag)
  END DO
  nini  = nx2*ny2*nz2
  WRITE(*,*) ' initialized nini=',nini,nx2,ny2,nz2
ELSE IF(nini /= nx2*ny2*nz2) THEN
  WRITE(*,*) ' nini,nx2,ny2,nz2=',nini,nx2,ny2,nz2
  STOP ' nx2, ny2 or/and nz2 in four3d not as initialized!'
END IF
IF(nxini == 0) THEN
  pforwx=fftw_plan_dft_1d(nx2,fftax,fftax,FFTW_FORWARD,FFTW_planflag)
  pbackx=fftw_plan_dft_1d(nx2,fftax,fftax,FFTW_BACKWARD,FFTW_planflag)
  nxini  = nx2
!       write(6,'(a)') ' x-fft initialized '
ELSE IF(nxini /= nx2) THEN
  STOP ' nx2 in four3d not as initialized!'
END IF
IF(nyini == 0) THEN
  pforwy=fftw_plan_dft_1d(ny2,fftay,fftay,FFTW_FORWARD,FFTW_planflag)
  pbacky=fftw_plan_dft_1d(ny2,fftay,fftay,FFTW_BACKWARD,FFTW_planflag)
  nyini  = ny2
!       write(6,'(a)') ' y-fft initialized '
ELSE IF(nyini /= ny2) THEN
  STOP ' ny2 in four3d not as initialized!'
END IF
IF(nzini == 0) THEN
  pforwz=fftw_plan_dft_1d(nz2,fftb(1,nx2),fftb(1,nx2),FFTW_FORWARD,FFTW_planflag)
  pbackz=fftw_plan_dft_1d(nz2,fftb(1,nx2),fftb(1,nx2),FFTW_BACKWARD,FFTW_planflag)
  pforwz1=fftw_plan_dft_1d(nz2,fftaz,fftaz,FFTW_FORWARD,FFTW_planflag)
  pbackz1=fftw_plan_dft_1d(nz2,fftaz,fftaz,FFTW_BACKWARD,FFTW_planflag)
  nzini  = nz2
!       write(6,'(a)') ' z-fft initialized '
ELSE IF(nzini /= nz2) THEN
  STOP ' nz2 in four3d not as initialized!'
END IF
wisdomtest=fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
IF (wisdomtest == 0) THEN
  WRITE(6,*) 'Error exporting wisdom to file wisdom_fftw.dat'
  WRITE(7,*) 'Error exporting wisdom to file wisdom_fftw.dat'
ELSE
  WRITE(*,*) ' suuccessfull export of wisdom to  wisdom_fftw.dat'
ENDIF
CALL fftw_forget_wisdom
#endif

#if(parayes)
IF (nini==0) THEN
  IF(myn == 0) THEN
    !Master node creates wisdom if necessary...
    wisdomtest=fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
    IF (wisdomtest == 0) THEN
      WRITE(6,*) 'wisdom_fftw.dat not found, creating it'
      WRITE(7,*) 'wisdom_fftw.dat not found, creating it'
    END IF
    pforw(0)=fftw_plan_dft_3d(nz2,ny2,nx2,ffta(:,:,:,0),ffta(:,:,:,0),FFTW_FORWARD,FFTW_planflag)
    pback(0)=fftw_plan_dft_3d(nz2,ny2,nx2,ffta(:,:,:,0),ffta(:,:,:,0),FFTW_BACKWARD,FFTW_planflag)
    nini  = nx2*ny2*nz2
    wisdomtest=fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
    IF (wisdomtest == 0) THEN
      WRITE(6,*) 'Error exporting wisdom to file wisdom_fftw.dat'
      WRITE(7,*) 'Error exporting wisdom to file wisdom_fftw.dat'
    END IF
  ENDIF
  CALL mpi_barrier(mpi_comm_world,mpi_ierror)
  !... then other nodes use it
  IF(myn/=0) THEN
    wisdomtest=fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
    pforw(0)=fftw_plan_dft_3d(nz2,ny2,nx2,ffta(:,:,:,0),ffta(:,:,:,0),FFTW_FORWARD,FFTW_planflag)
    pback(0)=fftw_plan_dft_3d(nz2,ny2,nx2,ffta(:,:,:,0),ffta(:,:,:,0),FFTW_BACKWARD,FFTW_planflag)
    nini  = nx2*ny2*nz2
  ENDIF
ELSE IF(nini /= nx2*ny2*nz2) THEN
  STOP ' nx2, ny2 or/and nz2 in four3d not as initialized!'
END IF
IF(nxini == 0) THEN
  IF(myn == 0) THEN
    pforwx=fftw_plan_dft_1d(nx2,fftax,fftax,FFTW_FORWARD,FFTW_planflag)
    pbackx=fftw_plan_dft_1d(nx2,fftax,fftax,FFTW_BACKWARD,FFTW_planflag)
    nxini  = nx2
    wisdomtest=fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
    IF (wisdomtest == 0) THEN
      WRITE(6,*) 'Error exporting wisdom to file wisdom_fftw.dat'
      WRITE(7,*) 'Error exporting wisdom to file wisdom_fftw.dat'
    END IF
  ENDIF
  CALL mpi_barrier(mpi_comm_world,mpi_ierror)
  IF(myn /= 0) THEN
    wisdomtest=fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
    pforwx=fftw_plan_dft_1d(nx2,fftax,fftax,FFTW_FORWARD,FFTW_planflag)
    pbackx=fftw_plan_dft_1d(nx2,fftax,fftax,FFTW_BACKWARD,FFTW_planflag)
    nxini  = nx2
  ENDIF
!       write(6,'(a)') ' x-fft initialized '
ELSE IF(nxini /= nx2) THEN
  STOP ' nx2 in four3d not as initialized!'
END IF
IF(nyini == 0) THEN
  IF(myn == 0) THEN
    pforwy=fftw_plan_dft_1d(ny2,fftay,fftay,FFTW_FORWARD,FFTW_planflag)
    pbacky=fftw_plan_dft_1d(ny2,fftay,fftay,FFTW_BACKWARD,FFTW_planflag)
    nyini  = ny2
    wisdomtest=fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
    IF (wisdomtest == 0) THEN
      WRITE(6,*) 'Error exporting wisdom to file wisdom_fftw.dat'
      WRITE(7,*) 'Error exporting wisdom to file wisdom_fftw.dat'
    END IF
  ENDIF
  CALL mpi_barrier(mpi_comm_world,mpi_ierror)
  IF(myn /= 0) THEN
    wisdomtest=fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
    pforwy=fftw_plan_dft_1d(ny2,fftay,fftay,FFTW_FORWARD,FFTW_planflag)
    pbacky=fftw_plan_dft_1d(ny2,fftay,fftay,FFTW_BACKWARD,FFTW_planflag)
    nyini  = ny2
  ENDIF
!       write(6,'(a)') ' y-fft initialized '
ELSE IF(nyini /= ny2) THEN
  STOP ' ny2 in four3d not as initialized!'
END IF
IF(nzini == 0) THEN
  IF(myn == 0) THEN
    pforwz=fftw_plan_dft_1d(nz2,fftb(1,nx2),fftb(1,nx2),FFTW_FORWARD,FFTW_planflag)
    pbackz=fftw_plan_dft_1d(nz2,fftb(1,nx2),fftb(1,nx2),FFTW_BACKWARD,FFTW_planflag)
    pforwz1=fftw_plan_dft_1d(nz2,fftaz,fftaz,FFTW_FORWARD,FFTW_planflag)
    pbackz1=fftw_plan_dft_1d(nz2,fftaz,fftaz,FFTW_BACKWARD,FFTW_planflag)
    nzini  = nz2
    wisdomtest=fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
    IF (wisdomtest == 0) THEN
      WRITE(6,*) 'Error exporting wisdom to file wisdom_fftw.dat'
      WRITE(7,*) 'Error exporting wisdom to file wisdom_fftw.dat'
    END IF
  ENDIF
  CALL mpi_barrier(mpi_comm_world,mpi_ierror)
  IF(myn /= 0) THEN
    wisdomtest=fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
    pforwz=fftw_plan_dft_1d(nz2,fftb(1,nx2),fftb(1,nx2),FFTW_FORWARD,FFTW_planflag)
    pbackz=fftw_plan_dft_1d(nz2,fftb(1,nx2),fftb(1,nx2),FFTW_BACKWARD,FFTW_planflag)
    pforwz1=fftw_plan_dft_1d(nz2,fftaz,fftaz,FFTW_FORWARD,FFTW_planflag)
    pbackz1=fftw_plan_dft_1d(nz2,fftaz,fftaz,FFTW_BACKWARD,FFTW_planflag)
    nzini  = nz2
  ENDIF
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
#endif

WRITE(*,*) ' end: fftay:',fftay

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
!COMPLEX(DP) ::ffftax(nx2,0:3),ffftay(ny2,0:3),ffftb(nz2,nx2,0:3)
INTEGER :: ithr

#if(netlib_fft)
DATA  nxini,nyini,nzini/0,0,0/ ! flag for initialization
#endif
COMPLEX(DP), ALLOCATABLE :: ffttax(:),ffttay(:),ffttaz(:),ffttb(:,:),fftta(:,:,:)

tnorm=1D0/SQRT(8D0*pi*pi*pi*REAL(nx2*ny2*nz2,DP))

!  here version using 3D FFTW
#if(fftw_cpu && !oldkinprop)
#if(paropenmp && dynopenmp)
  ithr = OMP_GET_THREAD_NUM()
  IF(ithr>nacthr) THEN
    WRITE(*,*) ' in kinprop: ithr,nacthr=',ithr,nacthr
    STOP "too large ITHR"
  END IF
!  WRITE(*,*) ' KINPROP: ithr=',ithr,nacthr
#else
  ithr = 0
#endif
facnr = 1D0/(nx2*ny2*nz2)
CALL copy1dto3d(q1,ffta(:,:,:,ithr),nx2,ny2,nz2)
!ffta(:,:,:,ithr)=reshape(q1,(/nx2,ny2,nz2/))
CALL fftw_execute_dft(pforw(ithr),ffta(:,:,:,ithr),ffta(:,:,:,ithr))
ffta(:,:,:,ithr) = akprop*ffta(:,:,:,ithr)
CALL fftw_execute_dft(pback(ithr),ffta(:,:,:,ithr),ffta(:,:,:,ithr))
!#if(paropenmp && dynopenmp)
CALL secopy3dto1d(ffta(:,:,:,ithr),q1,facnr,nx2,ny2,nz2)
!#else
!q1=reshape(ffta(:,:,:,ithr),(/kdfull2/))*facnr
!#endif
!WRITE(*,*) ' norms q1,q2=',SUM(CONJG(Q1)*Q1)*dvol,SUM(CONJG(Q2)*Q2)*dvol
#else

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
ALLOCATE(ffttax(nx2),ffttay(ny2),ffttaz(nz2),ffttb(nz2,nx2),fftta(nx2,ny2,nz2))

xfnorm = 1D0/nx2

DO i3=1,nz2
  DO i2=1,ny2
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
!      ffttax(MOD(i1+nx,nx2)+1)=q1(ind) ! copy to workspace
      fftax(modx(i1))=q1(ind) ! copy to workspace
    END DO
#if(netlib_fft)
    CALL dcftf1 (nx2,ffttax,wrkx,wsavex,ifacx) ! basic fft
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pforwx,ffttax,ffttax)
!    CALL fftw_execute_dft(pforwx,fftax(1),fftax(1))
#endif
    DO i1=1,nx2
      ffttax(i1) = akpropx(i1)*ffttax(i1)
    END DO
#if(netlib_fft)
    CALL dcftb1 (nx2,ffttax,wrkx,wsavex,ifacx) ! basic back fft
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pbackx,ffttax,ffttax)
!    CALL fftw_execute_dft(pbackx,fftax(1),fftax(1))
#endif
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
!      q2(ind)= ffttax(MOD(i1+nx,nx2)+1)*xfnorm
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
!      ffttay(MOD(i2+ny,ny2)+1) = q2(ind)
      fftay(mody(i2)) = q2(ind)
    END DO
#if(netlib_fft)
    CALL dcftf1 (ny2,ffttay,wrky,wsavey,ifacy)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pforwy,ffttay,ffttay)
!    CALL fftw_execute_dft(pforwy,fftay(1),fftay(1))
#endif
    DO i2=1,ny2
      ffttay(i2) = akpropy(i2)*ffttay(i2)
    END DO
#if(netlib_fft)
    CALL dcftb1 (ny2,ffttay,wrky,wsavey,ifacy)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pbacky,ffttay,ffttay)
!    CALL fftw_execute_dft(pbacky,fftay(1),fftay(1))
#endif
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
!      q2(ind)= ffttay(MOD(i2+ny,ny2)+1)*yfnorm
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
      ffttb(i3m,i1) = q2(ind)
    END DO
  END DO
  DO i1=1,nx2
#if(netlib_fft)
    CALL dcftf1 (nz2,ffttb(1,i1),wrkz,wsavez,ifacz)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pforwz,ffttb(1,i1),ffttb(1,i1))
#endif
    DO i3=1,nz2
      ffttb(i3,i1) = akpropz(i3)*ffttb(i3,i1)
    END DO
#if(netlib_fft)
    CALL dcftb1 (nz2,ffttb(1,i1),wrkz,wsavez,ifacz)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pbackz,ffttb(1,i1),ffttb(1,i1))
#endif
  END DO
  DO i3=1,nz2
    i3m = modz(i3)
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q1(ind)=ffttb(i3m,i1)*zfnorm
    END DO
  END DO
END DO

DEALLOCATE(ffttax,ffttay,ffttaz,ffttb,fftta)

#endif

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
    DO i3=1,nz2                  ! copy back
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      gradfout(ind)= fftaz(modz(i3))/nz2
    END DO
!
  END DO
END DO

#if(netlib_fft)
DEALLOCATE(fftaz)
#endif

RETURN

END SUBROUTINE  zgradient_rspace



! ******************************

SUBROUTINE  fftf(q1,q2)

! ******************************

USE params
#if(fftw_cpu)
USE FFTW
#endif
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
CALL copy1dto3d(q1,ffta(:,:,:,0),nx2,ny2,nz2)

CALL fftw_execute_dft(pforw(0),ffta,ffta)

CALL copy3dto1d(ffta(:,:,:,0),q2,tnorm,nx2,ny2,nz2)
#endif

RETURN
END SUBROUTINE  fftf

! ******************************

SUBROUTINE  fftback(q1,q2)

! ******************************

USE params
#if(fftw_cpu)
USE FFTW
#endif
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
CALL secopy1dto3d(q1,ffta(:,:,:,0),nx2,ny2,nz2)

CALL fftw_execute_dft(pback(0),ffta(:,:,:,0),ffta(:,:,:,0))

CALL secopy3dto1d(ffta(:,:,:,0),q2,facnr,nx2,ny2,nz2)
#endif

RETURN
END SUBROUTINE  fftback


! 1 "fft.f"

! ******************************

SUBROUTINE  rftf(q1,q2)

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
CALL copyr1dto3d(q1,ffta(:,:,:,0),nx2,ny2,nz2)

CALL fftw_execute_dft(pforw(0),ffta(:,:,:,0),ffta(:,:,:,0))

CALL copy3dto1d(ffta(:,:,:,0),q2,tnorm,nx2,ny2,nz2)
#endif

RETURN
END SUBROUTINE  rftf

! ******************************

SUBROUTINE  rfftback(q1,q3)
!SUBROUTINE  rfftback(q1,q2)

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
      q3(ind)= REAL(modx(i1)),DP)
    END DO
  END DO
END DO

DEALLOCATE(q2)
#endif

#if(fftw_cpu)
CALL secopy1dto3d(q1,ffta(:,:,:,0),nx2,ny2,nz2)

CALL fftw_execute_dft(pback(0),ffta(:,:,:,0),ffta(:,:,:,0))

CALL copyr3dto1d(ffta(:,:,:,0),q3,facnr,nx2,ny2,nz2)
#endif

RETURN
END SUBROUTINE  rfftback
!INCLUDE 'fftpack.F90'
!INCLUDE 'fftpack2.F90'

#if(fftw_cpu)

SUBROUTINE copy1dto3d(q1,ffta,nbx2,nby2,nbz2)

! ******************************

USE params

COMPLEX(DP), INTENT(IN)                      :: q1(kdfull2)
COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT)       :: ffta(nbx2,nby2,nbz2)
INTEGER,INTENT(IN) :: nbx2,nby2,nbz2

ind=0
DO i3=1,nbz2
  DO i2=1,nby2
    DO i1=1,nbx2
!      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
       ind=1+ind
!      ffta(MOD(i1+nx,nbx2)+1,MOD(i2+ny,nby2)+1,MOD(i3+nz,nbz2)+1)=q1(ind)
      ffta(i1,i2,i3)=q1(ind)
!DB      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
!DB      ffta(modx(i1),mody(i2),modz(i3))=q1(ind)
    END DO
  END DO
END DO

RETURN
END SUBROUTINE copy1dto3d

! ******************************

SUBROUTINE copy3dto1d(ffta,q2,coef,nbx2,nby2,nbz2)

! ******************************

USE params

COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN)        :: ffta(nbx2,nby2,nbz2)
COMPLEX(DP), INTENT(OUT)                     :: q2(kdfull2)

DO i3=1,nbz2
  DO i2=1,nby2
    DO i1=1,nbx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q2(ind)=coef*ffta(i1,i2,i3)
    END DO
  END DO
END DO

RETURN
END SUBROUTINE copy3dto1d

! ******************************

SUBROUTINE copyr1dto3d(q1,ffta,nbx2,nby2,nbz2)

! ******************************

USE params

REAL(DP), INTENT(IN)                      :: q1(kdfull2)
COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT)    :: ffta(nbx2,nby2,nbz2)
       
DO i3=1,nbz2
  DO i2=1,nby2
    DO i1=1,nbx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      ffta(modx(i1),mody(i2),modz(i3))=CMPLX(q1(ind),0D0,DP) 
    END DO
  END DO
END DO

RETURN
END SUBROUTINE copyr1dto3d

! ******************************

SUBROUTINE copyr3dto1d(ffta,q2,coef,nbx2,nby2,nbz2)

! ******************************

USE params

COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN)     :: ffta(nbx2,nby2,nbz2)
REAL(DP), INTENT(OUT)                     :: q2(kdfull2)

DO i3=1,nbz2
  DO i2=1,nby2
    DO i1=1,nbx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q2(ind)=REAL(ffta(modx(i1),mody(i2),modz(i3)),DP)*coef
    END DO
  END DO
END DO

RETURN
END SUBROUTINE copyr3dto1d

! ******************************

SUBROUTINE secopy1dto3d(q1,ffta,nbx2,nby2,nbz2)

! ******************************

USE params

COMPLEX(DP), INTENT(IN)                      :: q1(kdfull2)
COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT)       :: ffta(nbx2,nby2,nbz2)

DO i3=1,nbz2
  DO i2=1,nby2
    DO i1=1,nbx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      ffta(i1,i2,i3)=q1(ind)
    END DO
  END DO
END DO

RETURN
END SUBROUTINE secopy1dto3d

! ******************************

SUBROUTINE secopy3dto1d(ffta,q2,coef,nbx2,nby2,nbz2)

! ******************************

USE params

COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN)        :: ffta(nbx2,nby2,nbz2)
COMPLEX(DP), INTENT(OUT)                     :: q2(kdfull2)
REAL(DP),INTENT(IN) :: coef
INTEGER,INTENT(IN) :: nbx2,nby2,nbz2

ind=0
DO i3=1,nbz2
  DO i2=1,nby2
    DO i1=1,nbx2
!      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
!      q2(ind)=ffta(MOD(i1+nx,nbx2)+1,MOD(i2+ny,nby2)+1,MOD(i3+nz,nbz2)+1)*coef
      ind=1+ind
      q2(ind)=ffta(i1,i2,i3)*coef
!DB      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
!DB      q2(ind)=ffta(modx(i1),mody(i2),modz(i3))*coef
    END DO
  END DO
END DO

RETURN
END SUBROUTINE secopy3dto1d

! ******************************

SUBROUTINE fft_end()

! ******************************

USE FFTW

CALL fftw_destroy_plan(pforw(0))
CALL fftw_destroy_plan(pback(0))
CALL fftw_destroy_plan(pforwx)
CALL fftw_destroy_plan(pforwy)
CALL fftw_destroy_plan(pforwz)
CALL fftw_destroy_plan(pforwz1)
CALL fftw_destroy_plan(pbackx)
CALL fftw_destroy_plan(pbacky)
CALL fftw_destroy_plan(pbackz)
CALL fftw_destroy_plan(pbackz1)

DEALLOCATE(ak,akv)
DEALLOCATE(akpropx,akpropy,akpropz)
DEALLOCATE(fftax,fftay,fftaz,fftb,ffta)

RETURN
END SUBROUTINE fft_end

#endif

END MODULE kinetic
