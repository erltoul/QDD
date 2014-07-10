!------------------------------------------------------------
!
PROGRAM analyze_MP
!
!     to analyze 'pMP' file for PES
!     note: one may have to set the stacksize to 'unlimited'
!           else the code crashes with a 'segmentation fault' 
!
INTEGER,PARAMETER :: DP=KIND(1D0)      ! precision  setting

INTEGER, PARAMETER :: komega=3000      ! nr. os mesh points in frequency space
REAL(DP), PARAMETER :: delomega=0.002  ! mesh size in frequency space

INTEGER, PARAMETER :: ktimes=100000     ! max. nr. of time points -- to be set
INTEGER :: kpoints=32                   ! nr. measuring points
INTEGER :: kstate=14                    ! nr. states
CHARACTER (LEN=13) ::  name = 'na8-jel' ! qualifier

COMPLEX(DP),ALLOCATABLE :: q0(:,:,:)
COMPLEX(DP),ALLOCATABLE :: accum(:,:)
REAL(DP),ALLOCATABLE :: str(:)
REAL(DP),ALLOCATABLE :: spstr(:)
COMPLEX(DP) :: cfac
REAL(DP) :: time,time2,omega,timfin,deltime
INTEGER :: iomega,i,j,nbe,ntimes,itimes
REAL(DP) :: timeinp(0:ktimes)


  WRITE(*,*) ' enter qualifier for the input file "pMP.<name>":'
  READ(5,'(a)') name

  OPEN(01,file='pMP.'//name)

  READ(01,'(1x)')
  READ(01,'(2x,2i6)') kpoints,kstate

  IF(kpoints > 50) STOP ' too many measuring points (limit is 50)'

  ALLOCATE(q0(0:ktimes,kpoints,kstate))
  ALLOCATE(accum(kpoints,kstate))
  ALLOCATE(str(kpoints))
  ALLOCATE(spstr(kstate))

  DO j=1,kpoints
    READ(01,'(1x)')
  END DO
  DO itimes=0,ktimes
    READ(01,'(1f14.5,100e18.8)',end=99) time,q0(itimes,:,1)
    timeinp(itimes) = time
    WRITE(7,'(f12.6)') time
    ntimes = itimes
    timfin = time
    DO j=2,kstate
      READ(01,'(1f14.5,100e18.8)',err=98,end=98) time2,q0(itimes,:,j)
      IF(time2.ne.time) goto 98
    END DO
  END DO
  STOP ' file count exhausted, enhance KTIMES'

 98   STOP ' inconsistent input file'
    
 99   CONTINUE

  ! check times
!  IF(timeinp(0) .NE. 0D0) STOP ' first time point must be 0D0'
  time2 = timeinp(1)-timeinp(0)
  DO itimes=1,ntimes
    IF(timeinp(itimes).NE.itimes*time2) WRITE(6,'(a,i10,2f14.5)') &
       ' wrong time:',itimes,timeinp(itimes),itimes*time2
  END DO

  CLOSE(01)

  OPEN(02,file='pPES.'//name)
  OPEN(03,file='pspPES.'//name)
  timfin = (timfin-timeinp(0))/0.048D0        ! time in 1/Ry
  deltim = timfin/ntimes
  DO iomega=1,komega
    omega = iomega*delomega

    accum = CMPLX(0D0,0D0,DP)
    DO itimes=0,ntimes
      cfac   = EXP(CMPLX(0D0,omega*itimes*deltim,DP))
      accum = accum + cfac*q0(itimes,:,:)
    END DO

    str = 0D0
    DO j=1,kstate
      str = str + accum(:,j)*conjg(accum(:,j))
      spstr(j) =  SUM(accum(:,j)*conjg(accum(:,j)))
    END DO    
!    WRITE(02,'(f8.3,20(1pg12.5))') omega,str
!    WRITE(6,'(f8.3,20(1pg12.5))') omega,str
    WRITE(02,'(f8.3,200(1pg12.5))') omega,SUM(str),str     
    WRITE(6,'(f8.3,50(1pg12.5))') omega,str
    WRITE(03,'(f8.3,50(1pg12.5))') omega,SUM(spstr(1:kstate-2)),spstr     


  END DO

  CLOSE(02)
  CLOSE(03)

!  DEALLOCATE(q0)
  DEALLOCATE(accum)
  DEALLOCATE(str)
  DEALLOCATE(spstr)

STOP
END
