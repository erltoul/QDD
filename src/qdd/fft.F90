!This file is a part of PW-TELEMAN project.
!PW-TELEMAN is a Time-Dependent Electronic Dynamics in Molecules And Nanosystems library.
!Copyright (C) 2011-2015  Paul-Gerhard Reinhard, Eric Suraud, Florent Calvayrac,
!Phuong Mai Dinh, David Brusson, Philipp Wopperer, José María Escartín Esteban.
!
!PW-Teleman is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!PW-Teleman is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with PW-Teleman.  If not, see <http://www.gnu.org/licenses/>.

MODULE kinetic

! This module provides basic variables and arrays for FFT definition
! of kinetic energand Coulomb solver. Valid for FFTw3 as well as NETLIB.

#if(fftw_cpu)
USE, intrinsic :: iso_c_binding
#endif
USE params, ONLY: DP,numthr,PI
IMPLICIT NONE

SAVE
!     Arrays for kinetic energy and electronic propagation.
!     Kinetic energy coefficients in strange ordered Fourier space
!       akv  = Fourier-field for 0.5*k^2
!       ak   = Fourier-field for exp(i*dt*(h^2/2m)*k^2)

INTEGER,PRIVATE,ALLOCATABLE :: modx(:),mody(:),modz(:)
COMPLEX(DP),PARAMETER,PRIVATE :: eye=(0D0,1D0)


INTEGER, PRIVATE :: kfft,kfft2,kdfull2
INTEGER, PRIVATE :: kxmax,kymax,kzmax
INTEGER, PRIVATE :: iret


REAL(DP),ALLOCATABLE :: akv(:)
COMPLEX(DP), ALLOCATABLE :: akx(:),aky(:),akz(:) !DB
COMPLEX(DP),ALLOCATABLE :: akpropx(:),akpropy(:),akpropz(:),akprop(:,:,:)
COMPLEX(DP),ALLOCATABLE :: ak(:)

#if(netlib_fft)
REAl(DP), PRIVATE, ALLOCATABLE :: fftax(:),fftay(:),fftb(:,:)   ! Complexes stored in real arrays for NETLIB FFT library
REAL(DP), PRIVATE, ALLOCATABLE :: wrkx(:),wrky(:),wrkz(:)
REAL(DP), PRIVATE, ALLOCATABLE :: wsavex(:),wsavey(:),wsavez(:)
INTEGER, PRIVATE, ALLOCATABLE :: ifacx(:),ifacy(:),ifacz(:)
#endif
#if(fftw_cpu)
INTEGER,PUBLIC,SAVE :: FFTW_planflag
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

!  Initialize details for FFT, plans for FFTW3.
!
!  Input:
!    dx0,dy0,dz0  =  grid spacing
!    nx0,ny0,nz0  =  nr. of grid points
!    dt1          =  electronic time step
!    h2m          =  electronic hbar*+2/2m

#if(fftw_cpu)
USE FFTW
#endif
#if(parayes)
USE params, only : myn,numthr,nthr,mpi_ierror
INCLUDE 'mpif.h'
REAL(DP) :: is(mpi_status_size)
#endif

REAL(DP),INTENT(IN):: dx0, dy0, dz0
INTEGER,INTENT(IN):: nx0, ny0, nz0
REAL(DP),INTENT(IN):: dt1, h2m

INTEGER:: nx, nx2, ny, ny2, nz, nz2
INTEGER:: i1, i2, i3
REAL(DP):: dkx, dky, dkz
REAL(DP):: zkx, zky, zkz

#if(fftw_cpu)
INTEGER, SAVE ::  nxini=0,nyini=0,nzini=0,nini=0 ! flag for initialization
#endif

#if(netlib_fft|fftw_cpu)
INTEGER:: ind 
#endif

#if(fftw_cpu)
INTEGER :: i  !, nacthr
#endif

#if(paropenmp && !dynopenmp)
INTEGER :: omp_get_num_threads
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

ALLOCATE(modx(kxmax),mody(kymax),modz(kzmax))
#if(netlib_fft)
ALLOCATE(ak(kdfull2),akv(kdfull2))
ALLOCATE(akx(kdfull2),aky(kdfull2),akz(kdfull2)) !DB
ALLOCATE(akpropx(kxmax),akpropy(kymax),akpropz(kzmax))
ALLOCATE(akprop(kxmax,kymax,kzmax))
ALLOCATE(fftax(2*kxmax),fftay(2*kymax),fftb(2*kzmax,kxmax)) !Takes 2 reals to store 1 complex
ALLOCATE(wrkx(kfft2),wrky(kfft2),wrkz(kfft2))
ALLOCATE(wsavex(kfft2),wsavey(kfft2),wsavez(kfft2))
ALLOCATE(ifacx(kfft2),ifacy(kfft2),ifacz(kfft2))
#endif
#if(fftw_cpu)
ALLOCATE(ak(kdfull2),akv(kdfull2))
ALLOCATE(akx(kdfull2),aky(kdfull2),akz(kdfull2)) !DB
ALLOCATE(akpropx(kxmax),akpropy(kymax),akpropz(kzmax))
ALLOCATE(akprop(kxmax,kymax,kzmax))
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
  FFTW_planflag = FFTW_MEASURE  
!  FFTW_planflag = FFTW_PATIENT
!  FFTW_planflag = FFTW_EXHAUSTIVE
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
#endif

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
#if(fftwnomkl)
  wisdomtest=fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
#endif
  IF (wisdomtest == 0) THEN
    wisdomtest = fftw_import_system_wisdom()
    IF(wisdomtest == 0) THEN
      WRITE(6,*) 'wisdom_fftw.dat not found, creating it'
      WRITE(7,*) 'wisdom_fftw.dat not found, creating it'
    ELSE
      WRITE(*,*) 'wisdom from system'
    END IF
  END IF

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
#if(fftwnomkl)
wisdomtest=fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
#endif
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
#if(fftwnomkl)
    wisdomtest=fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
#endif
    IF (wisdomtest == 0) THEN
      WRITE(6,*) 'wisdom_fftw.dat not found, creating it'
      WRITE(7,*) 'wisdom_fftw.dat not found, creating it'
    END IF
    pforw(0)=fftw_plan_dft_3d(nz2,ny2,nx2,ffta(:,:,:,0),ffta(:,:,:,0),FFTW_FORWARD,FFTW_planflag)
    pback(0)=fftw_plan_dft_3d(nz2,ny2,nx2,ffta(:,:,:,0),ffta(:,:,:,0),FFTW_BACKWARD,FFTW_planflag)
    nini  = nx2*ny2*nz2
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
  IF(myn/=0) THEN
#if(fftwnomkl)
    wisdomtest=fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
#endif
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

SUBROUTINE kinprop(q1,q2)

! ******************************

!  Propagation with kinetic energy exp(-i*dt*e_kin).
!
!  Input/Output:
!    q1     = coordinate-space array for propagated wavefunction 
!    q2     = complex wavefunction array as workspace

USE params
#if(fftw_cpu)
USE FFTW
#endif

IMPLICIT NONE

COMPLEX(DP), INTENT(IN OUT)                  :: q1(kdfull2)
COMPLEX(DP), INTENT(OUT)                     :: q2(kdfull2)

#if(netlib_fft)
INTEGER, SAVE ::  nxini=0,nyini=0,nzini=0 ! flag for initialization
REAL(DP), ALLOCATABLE :: ffttax(:),ffttay(:),ffttaz(:),ffttb(:,:) ! Complexes stored in real arrays for NETLIB FFT library
INTEGER:: i1, i2, i3, i1m, i2m, i3m, ind 
INTEGER :: ic,ir   ! Index for real and complex components when stored in ffttax, fftay...
COMPLEX(DP) :: cmplxfac
#endif
#if(fftw_cpu)
COMPLEX(DP), ALLOCATABLE :: ffttax(:),ffttay(:),ffttaz(:),ffttb(:,:)
REAL(DP)::facnr
INTEGER :: ithr
#endif
REAL(DP):: tnorm, xfnorm, yfnorm, zfnorm

!
! switch to fast(?)  propagation, using 3D FFT if possible
!
IF(iffastpropag == 1) THEN

tnorm=1D0/SQRT(8D0*pi*pi*pi*REAL(nx2*ny2*nz2,DP))

!
!  here version using 3D FFTW
!

#if(fftw_cpu)
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
#if(netlib_fft)
ALLOCATE(ffttax(2*nx2),ffttay(2*ny2),ffttaz(2*nz2),ffttb(2*nz2,nx2))
#endif
#if(fftw_cpu)
ALLOCATE(ffttax(nx2),ffttay(ny2),ffttaz(nz2),ffttb(nz2,nx2))
#endif
xfnorm = 1D0/nx2

DO i3=1,nz2
  DO i2=1,ny2
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      i1m=MOD(i1+nx,nx2)+1
#if(netlib_fft)
      ic= 2*i1m
      ir= ic-1
      ffttax(ir)=REAL(q1(ind),DP) ! copy to workspace
      ffttax(ic)=AIMAG(q1(ind)) ! copy to workspace
#endif
#if(fftw_cpu)
      ffttax(i1m)=q1(ind) ! copy to workspace
#endif
    END DO
#if(netlib_fft)
    CALL dcftf1 (nx2,ffttax,wrkx,wsavex,ifacx) ! basic fft
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pforwx,ffttax,ffttax)
#endif
    DO i1=1,nx2
#if(netlib_fft)
      ic=2*i1
      ir=ic-1
      cmplxfac = akpropx(i1)*CMPLX(ffttax(ir),ffttax(ic),DP)
      ffttax(ir) = REAL(cmplxfac,DP) 
      ffttax(ic) = AIMAG(cmplxfac)
#endif
#if(fftw_cpu)
      ffttax(i1) = akpropx(i1)*ffttax(i1)
#endif
    END DO
#if(netlib_fft)
    CALL dcftb1 (nx2,ffttax,wrkx,wsavex,ifacx) ! basic back fft
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pbackx,ffttax,ffttax)
#endif
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      i1m=MOD(i1+nx,nx2)+1
#if(netlib_fft)
      ic= 2*i1m
      ir= ic-1
      q2(ind)= xfnorm* CMPLX(ffttax(ir),ffttax(ic),DP)
#endif
#if(fftw_cpu)
      q2(ind)= xfnorm*ffttax(i1m)
#endif
    END DO
  END DO
END DO

!      transformation in y-direction

yfnorm = 1D0/ny2
DO i3=1,nz2
  DO i1=1,nx2
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      i2m=MOD(i2+ny,ny2)+1
#if(netlib_fft)
      ic=2*i2m
      ir=ic-1
      ffttay(ir) = REAL(q2(ind),DP)
      ffttay(ic) = AIMAG(q2(ind))
#endif
#if(fftw_cpu)
      ffttay(i2m) = q2(ind)
#endif
    END DO
#if(netlib_fft)
    CALL dcftf1 (ny2,ffttay,wrky,wsavey,ifacy)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pforwy,ffttay,ffttay)
#endif
    DO i2=1,ny2
#if(netlib_fft)
      ic=2*i2
      ir=ic-1
      cmplxfac = akpropy(i2)*CMPLX(ffttay(ir), ffttay(ic),DP)
      ffttay(ir) = REAL(cmplxfac,DP)
      ffttay(ic) = AIMAG(cmplxfac)
#endif
#if(fftw_cpu)
      ffttay(i2) = akpropy(i2)*ffttay(i2)
#endif
    END DO
#if(netlib_fft)
    CALL dcftb1 (ny2,ffttay,wrky,wsavey,ifacy)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pbacky,ffttay,ffttay)
#endif
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      i2m=MOD(i2+ny,ny2)+1
#if(netlib_fft)
      ic=2*i2m
      ir=ic-1
      q2(ind)= yfnorm*CMPLX(ffttay(ir),ffttay(ic),DP)
#endif
#if(fftw_cpu)
      q2(ind)= yfnorm*ffttay(i2m)
#endif
    END DO
  END DO
END DO

!       propagation in z-direction

zfnorm = 1D0/nz2
DO i2=1,ny2
  DO i3=1,nz2
!    i3m = modz(i3)
    i3m = MOD(i3+nz,nz2)+1
#if(netlib_fft)
    ic=2*i3m
    ir=ic-1
#endif
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
#if(netlib_fft)
      ffttb(ir,i1) = REAL(q2(ind),DP)
      ffttb(ic,i1) = AIMAG(q2(ind))
#endif
#if(fftw_cpu)
      ffttb(i3m,i1) = q2(ind)
#endif
    END DO
  END DO
  DO i1=1,nx2
#if(netlib_fft)
    CALL dcftf1 (nz2,ffttb(:,i1),wrkz,wsavez,ifacz)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pforwz,ffttb(1,i1),ffttb(1,i1))
#endif
    DO i3=1,nz2
#if(netlib_fft)
      ic=2*i3
      ir=ic-1
      cmplxfac= akpropz(i3)*CMPLX( ffttb(ir,i1), ffttb(ic,i1), DP)
      ffttb(ir,i1) = REAL(cmplxfac,DP)
      ffttb(ic,i1) = AIMAG(cmplxfac)
#endif
#if(fftw_cpu)
      ffttb(i3,i1) = akpropz(i3)*ffttb(i3,i1)
#endif
    END DO
#if(netlib_fft)
    CALL dcftb1 (nz2,ffttb(:,i1),wrkz,wsavez,ifacz)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pbackz,ffttb(1,i1),ffttb(1,i1))
#endif
  END DO
  DO i3=1,nz2
!    i3m = modz(i3)
    i3m = MOD(i3+nz,nz2)+1
#if(netlib_fft)
    ic=2*i3m
    ir=ic-1
#endif    
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
#if(netlib_fft)
      q1(ind)=zfnorm*CMPLX(ffttb(ir,i1),ffttb(ic,i1),DP)
#endif
#if(fftw_cpu)
      q1(ind)=zfnorm*ffttb(i3m,i1)
#endif
    END DO
  END DO
END DO

DEALLOCATE(ffttax,ffttay,ffttaz,ffttb)

#endif

!
! iffastpropag==0 switch
!
ELSE

#if(netlib_fft|fftw_cpu)
    CALL fftf(q1,q2)
    WRITE(*,*) ak(1),q2(1)
    q2(1:kdfull2) = ak*q2(1:kdfull2)
    CALL fftback(q2,q1)
#endif

END IF


RETURN
END SUBROUTINE  kinprop

SUBROUTINE calc_ekin(psin,ekinout)

! Calculates kinetic energy for single particle state with
!  complex wavefunction 'psin'.

USE params
IMPLICIT NONE


COMPLEX(DP), INTENT(IN)    :: psin(kdfull2)
REAL(DP), INTENT(OUT)      :: ekinout

REAL(DP) :: sum0
REAl(DP) :: sumk, sum0ex
INTEGER ::ii
REAL(DP) :: vol
COMPLEX(DP),DIMENSION(:),ALLOCATABLE :: psi2

!------------------------------------------------------------------

ALLOCATE(psi2(kdfull2))

CALL fftf(psin,psi2)
sum0 = 0D0
sumk = 0D0
DO ii=1,kdfull2
  vol   = REAL(psi2(ii),DP)*REAL(psi2(ii),DP) +AIMAG(psi2(ii))*AIMAG(psi2(ii))
  sum0  = vol + sum0
  sumk  = vol*akv(ii) + sumk
END DO
sum0ex = 1D0/((2D0*PI)**3*dx*dy*dz)
ekinout = sumk/sum0ex
!WRITE(6,*) ' sum0,sum0ex=',sum0,sum0ex


DEALLOCATE(psi2)

RETURN
END SUBROUTINE calc_ekin


! ******************************

SUBROUTINE  gradient(fin,gradfout,idirec)

!  The gradient of complex field 'fin' in direction 'idirec'
!    (x =1, y=2, z=3).
!  The fields are given in Fourier space and the gradient is
!   applied as product with 'kx', 'ky', 'kz'.
!

USE params
IMPLICIT NONE

COMPLEX(DP), INTENT(IN)    :: fin(kdfull2)
COMPLEX(DP), INTENT(OUT)   :: gradfout(kdfull2)
INTEGER, INTENT(IN)        :: idirec


! ************************************************************


#if(netlib_fft|fftw_cpu)
  IF(idirec == 1) THEN
!       x-derivative
          gradfout = -akx*fin
  ELSEIF(idirec == 2) THEN
!       y-derivative
          gradfout = -aky*fin
  ELSEIF(idirec == 3) THEN
!       z-derivative
          gradfout = -akz*fin
  ELSE
    STOP ' RGRADIENT called with invalid IDIREC'
  ENDIF
#endif


RETURN
END SUBROUTINE  gradient

! ******************************

SUBROUTINE  xgradient_rspace(fin,gradfout)

!  The gradient of the complex field 'fin' in x-direction.
!  The fields are given in coordinate space. The gradient is
!  evaluated in k_x space.
!

USE params
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT NONE

COMPLEX(DP), INTENT(IN)   :: fin(kdfull2)
COMPLEX(DP), INTENT(OUT)  :: gradfout(kdfull2)

INTEGER:: i1, i2, i3, ind 
#if(netlib_fft)
INTEGER:: ic, ir
#endif
REAL(DP):: dkx, zkx
! ************************************************************

dkx=pi/(dx*nx)
DO i3=1,nz2
  DO i2=1,ny2
!                 forward transform along x
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
#if(netlib_fft)
      ic=2*modx(i1)
      ir=ic-1
      fftax(ir)=REAL(fin(ind),DP)  ! copy to workspace
      fftax(ic)=AIMAG(fin(ind))  ! copy to workspace
#endif
#if(fftw_cpu)
      fftax(modx(i1))=fin(ind)  ! copy to workspace
#endif
    END DO
#if(netlib_fft)
    CALL dcftf1 (nx2,fftax,wrkx,wsavex,ifacx)    ! basic fft
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pforwx,fftax,fftax)    ! basic fft
#endif
!                 multiply by k_x factor in k_x space
    DO i1=1,nx2
      IF(i1 == (nx+1)) THEN
        zkx = 0D0
      ELSE IF(i1 > (nx+1)) THEN
        zkx=(i1-nx2-1)*dkx
      ELSE
        zkx=(i1-1)*dkx
      END IF
#if(netlib_fft)
      ic=2*i1
      ir=ic-1
      fftax(ir) = -fftax(ic)*zkx    ! *eye*zkx
      fftax(ic) =  fftax(ir)*zkx
#endif
#if(fftw_cpu)
      fftax(i1) = fftax(i1)*eye*zkx
#endif
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
#if(netlib_fft)
      ic=2*modx(i1)
      ir=ic-1
      gradfout(ind)= CMPLX(fftax(ir), fftax(ic),DP)/nx2
#endif
#if(fftw_cpu)
      gradfout(ind)= fftax(modx(i1))/nx2
#endif
    END DO
  END DO
END DO

RETURN

END SUBROUTINE  xgradient_rspace



SUBROUTINE  ygradient_rspace(fin,gradfout)

!  The gradient of the complex field 'fin' in y-direction.
!  The fields are given in coordinate space. The gradient is
!  evaluated in k_y space.
!

USE params
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT NONE

COMPLEX(DP), INTENT(IN)       :: fin(kdfull2)
COMPLEX(DP), INTENT(OUT)   :: gradfout(kdfull2)

#if(netlib_fft)
INTEGER:: ic, ir
#endif
INTEGER:: i1, i2, i3, ind 
REAL(DP):: dky,zky
! ************************************************************

dky=pi/(dy*ny)
DO i3=1,nz2
  DO i1=1,nx2
!                 forward transform along y
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
#if(netlib_fft)
      ic=mody(i2)
      ir=ic-1
      fftay(ir) = REAL(fin(ind),DP)
      fftay(ic) = AIMAG(fin(ind))
#endif
#if(fftw_cpu)
      fftay(mody(i2)) = fin(ind)
#endif
    END DO
#if(netlib_fft)
    CALL dcftf1 (ny2,fftay,wrky,wsavey,ifacy)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pforwy,fftay,fftay)
#endif
!                 multiply by k_y factor in k_y space
    DO i2=1,ny2
      IF(i2 == (ny+1)) THEN
        zky = 0D0
      ELSE IF(i2 > (ny+1)) THEN
        zky=(i2-ny2-1)*dky
      ELSE
        zky=(i2-1)*dky
      END IF
#if(netlib_fft)
      ic=2*i2
      ir=ic-1
      fftay(ir) = -fftay(ic)*zky  
      fftay(ic) =  fftay(ir)*zky
#endif
#if(fftw_cpu)
      fftay(i2) = fftay(i2)*eye*zky
#endif
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
#if(netlib_fft)
      ic=2*mody(i2)
      ir=ic-1
      gradfout(ind)= CMPLX(fftay(ir), fftay(ic), DP)/ny2
#endif
#if(fftw_cpu)
      gradfout(ind)= fftay(mody(i2))/ny2
#endif
    END DO
  END DO
END DO

RETURN

END SUBROUTINE  ygradient_rspace


SUBROUTINE  zgradient_rspace(fin,gradfout)

!  The gradient of the complex field 'fin' in z-direction.
!  The fields are given in coordinate space. The gradient is
!  evaluated in k_z space.
!

USE params
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT NONE

COMPLEX(DP), INTENT(IN)    :: fin(kdfull2)
COMPLEX(DP), INTENT(OUT)   :: gradfout(kdfull2)

#if(netlib_fft)
INTEGER:: ic, ir
#endif
#if(netlib_fft)
REAL(DP), ALLOCATABLE :: fftaz(:)
#endif
INTEGER:: i1, i2, i3, i3m, ind 
REAL(DP):: dkz,zkz


! ************************************************************

#if(netlib_fft)
ALLOCATE(fftaz(2*nz2))
#endif

dkz=pi/(dz*nz)
DO i2=1,ny2
  DO i1=1,nx2
!                 forward transform along z
    DO i3=1,nz2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      i3m = MOD(i3+nz,nz2)+1
!      i3m = modz(i3)
#if(netlib_fft)
      ic=2*i3m
      ir=ic-1
      fftaz(ir) = REAL(fin(ind),DP)
      fftaz(ic) = AIMAG(fin(ind))
#endif
#if(fftw_cpu)
      fftaz(i3m) = fin(ind)
#endif
    END DO
#if(netlib_fft)
    CALL dcftf1 (nz2,fftaz,wrkz,wsavez,ifacz)
#endif
#if(fftw_cpu)
    CALL fftw_execute_dft(pforwz,fftaz,fftaz)
#endif
!                 multiply by k_z factor in k_z space
    DO i3=1,nz2
      IF(i3 == (nz+1)) THEN
        zkz = 0D0
      ELSE IF(i3 > (nz+1)) THEN
        zkz=(i3-nz2-1)*dkz
      ELSE
        zkz=(i3-1)*dkz
      END IF
#if(netlib_fft)
      ic=2*i3
      ir=ic-1
      fftaz(ir) = -fftaz(ic)*zkz   !  eye*zkz
      fftaz(ic) =  fftaz(ir)*zkz
#endif
#if(fftw_cpu)
      fftaz(i3) = fftaz(i3)*eye*zkz
#endif
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
      i3m = MOD(i3+nz,nz2)+1
#if(fftw_cpu)
      gradfout(ind)= fftaz(modz(i3))/nz2
#endif
#if(netlib_fft)
      ic=2*i3m
      ir=ic-1
      gradfout(ind)= CMPLX(fftaz(ir),fftaz(ic),DP)/nz2
#endif
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

!#if(netlib_fft|fftw_cpu)
SUBROUTINE  fftf(q1,q2)
!#endif

! Forward 3D FFT from 'q1' in r-space to 'q2' in k-space.

USE params
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT NONE

COMPLEX(DP), INTENT(IN)      :: q1(kdfull2)
COMPLEX(DP), INTENT(OUT)     :: q2(kdfull2)

#if(netlib_fft)
INTEGER:: i1, i2, i3, i3m, ind 
INTEGER:: ic, ir
#endif
REAL(DP):: tnorm

#if(netlib_fft)
INTEGER,SAVE :: nxini=0,nyini=0,nzini=0     ! flag for initialization
#endif

tnorm=1D0/SQRT(8D0*pi*pi*pi*REAL(nx2*ny2*nz2,DP))

#if(netlib_fft)
!     check initialization

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

DO i3=1,nz2
  DO i2=1,ny2
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      ic=2*modx(i1)
      ir=ic-1
      fftax(ir)=REAL(q1(ind),DP)  ! copy to workspace
      fftax(ic)=AIMAG(q1(ind))
    END DO
    CALL dcftf1 (nx2,fftax,wrkx,wsavex,ifacx)    ! basic fft
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      ic=2*i1
      ir=ic-1
      q2(ind)= CMPLX(fftax(ir),fftax(ic),DP)        !  copy back in strange order
    END DO
  END DO
END DO

!     transformation in y-direction

DO i3=1,nz2
  DO i1=1,nx2
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      ic=2*mody(i2)
      ir=ic-1
      fftay(ir) = REAL(q2(ind),DP)
      fftay(ic) = AIMAG(q2(ind))
    END DO
    CALL dcftf1 (ny2,fftay,wrky,wsavey,ifacy)
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      ic=2*i2
      ir=ic-1
      q2(ind)= CMPLX(fftay(ir),fftay(ic),DP)
    END DO
  END DO
END DO


!     transformation in z-direction

DO i2=1,ny2
  DO i3=1,nz2
!    i3m = modz(i3)
    i3m = MOD(i3+nz,nz2)+1
    ic=2*i3m
    ir=ic-1
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftb(ir,i1) = REAL(q2(ind),DP)
      fftb(ic,i1) = AIMAG(q2(ind))
    END DO
  END DO
  DO i1=1,nx2
    CALL dcftf1 (nz2,fftb(:,i1),wrkz,wsavez,ifacz)
  END DO
  DO i3=1,nz2
    ic=2*i3
    ir=ic-1
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q2(ind)= tnorm*CMPLX(fftb(ir,i1),fftb(ic,i1),DP)
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


!#if(netlib_fft|fftw_cpu)
SUBROUTINE  fftback(q1,q2)
!#endif

! Backward 3D FFT from 'q1' in k-space to 'q2' in r-space.

USE params
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT NONE

COMPLEX(DP), INTENT(IN)                      :: q1(kdfull2)
COMPLEX(DP), INTENT(OUT)                     :: q2(kdfull2)
#if(netlib_fft)
INTEGER:: i1, i2, i3, i3m, ind 
INTEGER:: ic, ir
#endif
REAL(DP):: facnr

facnr =SQRT(8D0*pi*pi*pi)/SQRT(REAL(nx2*ny2*nz2,DP))

#if(netlib_fft)
!     transformation in z-direction
DO i2=1,ny2
  DO i3=1,nz2
    ic=2*i3
    ir=ic-1
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftb(ir,i1) = facnr*REAL(q1(ind),DP)
      fftb(ic,i1) = facnr*AIMAG(q1(ind))
    END DO
  END DO
  DO i1=1,nx2
    CALL dcftb1 (nz2,fftb(:,i1),wrkz,wsavez,ifacz)    ! basic fft
  END DO
  DO i3=1,nz2                  ! copy back
    i3m = MOD(i3+nz,nz2)+1
    ic=2*i3m
    ir=ic-1
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q2(ind)= CMPLX(fftb(ir,i1), fftb(ic,i1),DP)
    END DO
  END DO
END DO

!     transformation in y-direction

DO i3=1,nz2
  DO i1=1,nx2
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      ic=2*i2
      ir=ic-1
      fftay(ir) = REAL(q2(ind),DP)
      fftay(ic) = AIMAG(q2(ind))
    END DO
    CALL dcftb1 (ny2,fftay,wrky,wsavey,ifacy)
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      ic=2*mody(i2)
      ir=ic-1
      q2(ind)= CMPLX(fftay(ir),fftay(ic),DP)
    END DO
  END DO
END DO


!     transformation in x-direction

DO i3=1,nz2
  DO i2=1,ny2
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      ic=2*i1
      ir=ic-1
      fftax(ir)=REAL(q2(ind),DP)
      fftax(ic)=AIMAG(q2(ind))
    END DO
    CALL dcftb1 (nx2,fftax,wrkx,wsavex,ifacx)
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      ic=2*modx(i1)
      ir=ic-1
      q2(ind)= CMPLX(fftax(ir),fftax(ic),DP)
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




!#if(netlib_fft|fftw_cpu)
SUBROUTINE  rftf(q1,q2)
!#endif

! Forward 3D FFT from real r-space array 'q1' to complex k-space 'q2'.

USE params
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT NONE

REAL(DP), INTENT(IN)       :: q1(kdfull2)
COMPLEX(DP), INTENT(OUT)   :: q2(kdfull2)

#if(netlib_fft)
INTEGER,SAVE :: nxini=0,nyini=0,nzini=0     ! flag for initialization
INTEGER:: i1, i2, i3, ind 
INTEGER:: ic, ir
#endif
REAL(DP):: tnorm

tnorm=1D0/SQRT(8D0*pi*pi*pi*REAL(nx2*ny2*nz2,DP))

#if(netlib_fft)
!     check initialization
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

DO i3=1,nz2
  DO i2=1,ny2
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      ic=2*modx(i1)
      ir=ic-1
      fftax(ir)=q1(ind) ! copy to workspace
      fftax(ic)=0D0
    END DO
    CALL dcftf1 (nx2,fftax,wrkx,wsavex,ifacx)    ! basic fft
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      ic=2*i1
      ir=ic-1
      q2(ind)= CMPLX(fftax(ir),fftax(ic),DP)        !  copy back in strange order
    END DO
  END DO
END DO

!     transformation in y-direction

DO i3=1,nz2
  DO i1=1,nx2
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      ic=2*mody(i2)
      ir=ic-1
      fftay(ir) = REAL(q2(ind),DP)
      fftay(ic) = AIMAG(q2(ind))
    END DO
    CALL dcftf1 (ny2,fftay,wrky,wsavey,ifacy)
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      ic=2*i2
      ir=ic-1
      q2(ind)= CMPLX(fftay(ir),fftay(ic),DP)
    END DO
  END DO
END DO


!     transformation in z-direction

DO i2=1,ny2
  DO i3=1,nz2
    ic=2*modz(i3)
    ir=ic-1
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftb(ir,i1) = REAL(q2(ind),Dp)
      fftb(ic,i1) = AIMAG(q2(ind))
    END DO
  END DO
  DO i1=1,nx2
    CALL dcftf1 (nz2,fftb(1,i1),wrkz,wsavez,ifacz)
  END DO
  DO i3=1,nz2
    ic=2*i3
    ir=ic-1
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q2(ind)=tnorm*CMPLX(fftb(ir,i1),fftb(ic,i1),DP)
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


!#if(netlib_fft|fftw_cpu)
SUBROUTINE  rfftback(q1,q3)
!SUBROUTINE  rfftback(q1,q2)
!#endif

! Backward 3D FFT from complex k-space 'q1' to real r-space array 'q3'.

USE params
#if(fftw_cpu)
USE FFTW
#endif
IMPLICIT NONE

COMPLEX(DP), INTENT(IN)    :: q1(kdfull2)
REAL(DP), INTENT(OUT)      :: q3(kdfull2)

#if(netlib_fft)
COMPLEX(DP),ALLOCATABLE :: q2(:)
INTEGER :: i1, i2, i3, ind
INTEGER :: ir, ic
#endif
REAL(DP):: facnr

!      data  nxini,nyini,nzini/0,0,0/  ! flag for initialization
! nxyf=nx2*ny2
!      nyf=nx2
!      nx=nx2/2
! ny=ny2/2
! nz=nz2/2

facnr =SQRT(8D0*pi*pi*pi)/SQRT(REAL(nx2*ny2*nz2,DP))

#if(netlib_fft)

ALLOCATE(q2(kdfull2))

!     transformation in z-direction

DO i2=1,ny2
  DO i3=1,nz2
      ic=2*i3
      ir=ic-1
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      fftb(ir,i1) = facnr*REAL(q1(ind),DP)
      fftb(ic,i1) = facnr*AIMAG(q1(ind))
    END DO
  END DO
  DO i1=1,nx2
    CALL dcftb1 (nz2,fftb(:,i1),wrkz,wsavez,ifacz)    ! basic fft
  END DO
  DO i3=1,nz2                  ! copy back
    ic=2*modz(i3)
    ir=ic-1
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      q2(ind)= CMPLX(fftb(ir,i1), fftb(ic,i1),DP)
    END DO
  END DO
END DO

!     transformation in y-direction

DO i3=1,nz2
  DO i1=1,nx2
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      ic=2*i2
      ir=ic-1
      fftay(ir) = REAL(q2(ind),DP)
      fftay(ic) = AIMAG(q2(ind))
    END DO
    CALL dcftb1 (ny2,fftay,wrky,wsavey,ifacy)
    DO i2=1,ny2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      ic=2*mody(i2)
      ir=ic-1
      q2(ind) = CMPLX(fftay(ir),fftay(ic),DP)
    END DO
  END DO
END DO


!     transformation in x-direction

DO i3=1,nz2
  DO i2=1,ny2
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
      ic=2*i1
      ir=ic-1
      fftax(ir)=REAL(q2(ind),DP)
      fftax(ic)=AIMAG(q2(ind))
    END DO
    CALL dcftb1 (nx2,fftax,wrkx,wsavex,ifacx)
    DO i1=1,nx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
!      q2(ind)= REAL(fftax(MOD(i1+nx,nx2)+1),DP)
      ir=2*modx(i1)-1
      q3(ind)= fftax(ir)
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

#if(fftw_cpu)

SUBROUTINE copy1dto3d(vec1d,vec3d,nbx2,nby2,nbz2)

! Copies 3D array from linear storage to 3D storage (both complex).
! 
! Input:
!   vec1d           = array in linear storage
!   nbx2,nby2,nbz2  = dimensions of 3D array
! Output:
!   vec3d           =  array in 3D storage

USE params
IMPLICIT NONE

INTEGER,INTENT(IN)      :: nbx2
INTEGER,INTENT(IN)      :: nby2
INTEGER,INTENT(IN)      :: nbz2
COMPLEX(DP), INTENT(IN)                      :: vec1d(kdfull2)
COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT)       :: vec3d(nbx2,nby2,nbz2)

INTEGER:: i1, i2, i3, ind

ind=0
DO i3=1,nbz2
  DO i2=1,nby2
    DO i1=1,nbx2
      ind=ind+1
      vec3d(i1,i2,i3)=vec1d(ind)
    END DO
  END DO
END DO

RETURN
END SUBROUTINE copy1dto3d

! ******************************

SUBROUTINE copyr1dto3d(vec1d,vec3d,nbx2,nby2,nbz2)

! Copies 3D array from real & linear storage to complex & 3D storage.
! 
! Input:
!   vec1d           = array in linear storage
!   nbx2,nby2,nbz2  = dimensions of 3D array
! Output:
!   vec3d           =  array in 3D storage

USE params
IMPLICIT NONE

INTEGER,INTENT(IN)        :: nbx2
INTEGER,INTENT(IN)        :: nby2
INTEGER,INTENT(IN)        :: nbz2
REAL(DP), INTENT(IN)                      :: vec1d(kdfull2)
COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT)    :: vec3d(nbx2,nby2,nbz2)
INTEGER:: i1, i2, i3, ind

ind=0 
DO i3=1,nbz2
  DO i2=1,nby2
    DO i1=1,nbx2
      ind=ind+1
      vec3d(i1,i2,i3)=CMPLX(vec1d(ind),0D0,DP) 
    END DO
  END DO
END DO

RETURN
END SUBROUTINE copyr1dto3d

! ******************************

SUBROUTINE secopy1dto3d(vec1d,vec3d,nbx2,nby2,nbz2)

! Copies 3D array from linear storage to 3D storage (both complex).
! 
! Input:
!   vec1d           = array in linear storage
!   nbx2,nby2,nbz2  = dimensions of 3D array
! Output:
!   vec3d           =  array in 3D storage
!
! Does the same as 'copy1dto3d' but with faster (?) index computation.

USE params
IMPLICIT NONE

INTEGER,INTENT(IN)                           :: nbx2
INTEGER,INTENT(IN)                           :: nby2
INTEGER,INTENT(IN)                           :: nbz2
COMPLEX(DP), INTENT(IN)                      :: vec1d(kdfull2)
COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT)       :: vec3d(nbx2,nby2,nbz2)
INTEGER:: i1, i2, i3, ind

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
#if(fftw_cpu)
SUBROUTINE copy3dto1d(vec3d,vec1d,coef,nbx2,nby2,nbz2)
#else
SUBROUTINE copy3dto1d(vec3d,vec1d,nbx2,nby2,nbz2)
#endif

! Copies 3D array from 3D storage to linear storage (both complex).
! 
! Input:
!   vec3d           = array in 3D storage
!   nbx2,nby2,nbz2  = dimensions of 3D array
! Output:
!   vec1d           =  array in linear storage
!

USE params
IMPLICIT NONE

#if(fftw_cpu)
REAL(DP),INTENT(IN)                          ::coef
#endif
INTEGER,INTENT(IN)                           :: nbx2
INTEGER,INTENT(IN)                           :: nby2
INTEGER,INTENT(IN)                           :: nbz2
COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN)        :: vec3d(nbx2,nby2,nbz2)
COMPLEX(DP), INTENT(OUT)                     :: vec1d(kdfull2)
INTEGER:: i1, i2, i3, ind

DO i3=1,nbz2
  DO i2=1,nby2
    DO i1=1,nbx2
      ind=(i3-1)*nxyf+(i2-1)*nyf+i1
#if(fftw_cpu)
      vec1d(ind)=coef*vec3d(i1,i2,i3)
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

! Copies 3D array from 3D storage to linear storage.
! 
! Input:
!   vec3d           = complex array in 3D storage
!   nbx2,nby2,nbz2  = dimensions of 3D array
! Output:
!   vec1d           = real array in linear storage
!

USE params
IMPLICIT NONE

#if(fftw_cpu)
REAL(DP),INTENT(IN)                     ::coef
#endif
INTEGER,INTENT(IN)                      :: nbx2
INTEGER,INTENT(IN)                      :: nby2
INTEGER,INTENT(IN)                      :: nbz2
COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN)   :: vec3d(nbx2,nby2,nbz2)
REAL(DP), INTENT(OUT)                   :: vec1d(kdfull2)

INTEGER:: i1, i2, i3, ind

ind=0
DO i3=1,nbz2
  DO i2=1,nby2
    DO i1=1,nbx2
      ind=ind+1
      !ind=(i3-1)*nxyf+(i2-1)*nyf+i1
#if(fftw_cpu)
      !vec1d(ind)=REAL(vec3d(modx(i1),mody(i2),modz(i3)),DP)*coef
      vec1d(ind)=REAL(vec3d(i1,i2,i3),DP)*coef
#else
      !vec1d(ind)=REAL(vec3d(modx(i1),mody(i2),modz(i3)),DP)
      vec1d(ind)=REAL(vec3d(i1,i2,i3),DP)
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

! Copies 3D array from 3D storage to linear storage (both complex).
! 
! Input:
!   vec3d           = array in 3D storage
!   nbx2,nby2,nbz2  = dimensions of 3D array
! Output:
!   vec1d           =  array in linear storage
!
! Does the same as 'copy1dto3d' but with faster (?) index computation.


USE params
IMPLICIT NONE

#if(fftw_cpu)
REAL(DP),INTENT(IN)                     ::coef
#endif
INTEGER,INTENT(IN)                      :: nbx2
INTEGER,INTENT(IN)                      :: nby2
INTEGER,INTENT(IN)                      :: nbz2
COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN)   :: vec3d(nbx2,nby2,nbz2)
COMPLEX(DP), INTENT(OUT)                :: vec1d(kdfull2)

INTEGER:: i1, i2, i3, ind

ind=0 
DO i3=1,nbz2
  DO i2=1,nby2
    DO i1=1,nbx2
      ind=ind+1
      !ind=(i3-1)*nxyf+(i2-1)*nyf+i1
#if(fftw_cpu)
      vec1d(ind)=vec3d(i1,i2,i3)*coef
      !vec1d(ind)=vec3d(i1,i2,i3)*coef
#else
      vec1d(ind)=vec3d(i1,i2,i3)
      !vec1d(ind)=vec3d(modx(i1),mody(i2),modz(i3))
#endif
    END DO
  END DO
END DO

RETURN
END SUBROUTINE secopy3dto1d

! ******************************

SUBROUTINE fft_end()

! FFTW epilogue.


#if(fftw_cpu)
USE FFTW
IMPLICIT NONE

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
#endif


RETURN
END SUBROUTINE fft_end

#endif


END MODULE kinetic
