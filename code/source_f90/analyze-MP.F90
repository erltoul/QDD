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
INTEGER :: kpoints=3                   ! nr. measuring points
INTEGER :: kstate=8                    ! nr. states
CHARACTER (LEN=13) ::  name = 'na8-jel' ! qualifier

COMPLEX(DP),ALLOCATABLE :: q0(:,:,:)
COMPLEX(DP),ALLOCATABLE :: accum(:,:)
REAL(DP),ALLOCATABLE :: str(:)
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

  DO j=1,kpoints
    READ(01,'(1x)')
  END DO
  DO itimes=0,ktimes
    READ(01,'(1f14.5,100e18.8)',end=99) time,q0(itimes,:,1)
    timeinp(itimes) = time
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
  IF(timeinp(0) .NE. 0D0) STOP ' first time point must be 0D0'
  time2 = timeinp(1)
  DO itimes=1,ntimes
    IF(timeinp(itimes).NE.itimes*time2) WRITE(6,'(a,i10,2f14.5)') &
       ' wrong time:',itimes,timeinp(itimes),itimes*time2
  END DO

  CLOSE(01)

  OPEN(02,file='pPES.'//name)
  timfin = timfin/0.048D0        ! time in 1/Ry
  deltim = timfin/ntimes
  DO iomega=1,komega
    omega = iomega*delomega

    accum = cmplx(0D0,0D0)
    DO itimes=0,ntimes
      cfac   = cexp(cmplx(0D0,omega*itimes*deltim))
      accum = accum + cfac*q0(itimes,:,:)
    END DO

    str = 0D0
    DO j=1,kstate
      str = str + accum(:,j)*conjg(accum(:,j))
    END DO    

!    WRITE(02,'(f8.3,20(1pg12.5))') omega,str
!    WRITE(6,'(f8.3,20(1pg12.5))') omega,str
    WRITE(02,'(f8.3,27(1pg12.5))') omega,str     
    WRITE(6,'(f8.3,27(1pg12.5))') omega,str


  END DO

  CLOSE(02)

!  DEALLOCATE(q0)
  DEALLOCATE(accum)
  DEALLOCATE(str)

STOP
END
