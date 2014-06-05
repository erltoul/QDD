#include"define.h"
 
SUBROUTINE genermowf(psiom,nmxst)

! generates orbital molecular wave function: psiom


USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
#if(parayes)
INCLUDE 'mpif.h'
#endif


REAL(DP), INTENT(OUT)                        :: psiom(kdfull2,kstate)
INTEGER, INTENT(OUT) :: nmxst(1:ng)       ! maximum nr. for each atom

INTEGER,ALLOCATABLE :: nactst(:)       ! keep track of actual atomic state
INTEGER :: nnodes(1:3,1:10)
DATA nnodes/0,0,0,   1,0,0,  0,1,0,   0,0,1,   2,0,0,  &
    1,1,0,   0,2,0,  1,0,1,   0,1,1,   0,0,2/

!---------------------------------------------------------------------

!      reset wavefunctions 'psiom'

ALLOCATE(nactst(1:ng))
nactst=0

WRITE(6,*) ' GENERLCGO entered. NION,NSTATE=',nion,nstate
WRITE(6,*) ' CH(ion):',(ch(np(ion)),ion=1,nion)
CALL flush(6)

DO nbr=1,nstate
  DO i=1,kdfull2
    psiom(i,nbr)=0D0
  END DO
END DO

!      book-keeping of number of states

nmaxval = 0
DO ion=1,nion
  nmxst(ion) = ch(np(ion))
  nmaxval = nmaxval + ch(np(ion))
END DO
ndiff = nmaxval - nclust
ncycle = ndiff/nion+1
IF(ndiff > 0) THEN
  nadd = +1
ELSE
  nadd = -1
END IF
!nadd=0     !  ?? check the initialization
IF(ndiff /= 0) THEN
  nadded = 0
  DO ncy=1,ncycle
    DO ion=nion,1,-1
      nmxst(ion) = nmxst(ion)-nadd
      nadded = nadded + nadd
      IF(nadded == ndiff) GO TO 19
    END DO
  END DO
  19     CONTINUE
END IF
WRITE(6,'(a,5i5)')  &
    ' ndiff,nmaxval,nstate,ncycle,nadd=',ndiff,nmaxval,nstate, ncycle,nadd
WRITE(6,'(a,100i3)') ' nmxst:',(nmxst(ion),ion=1,nion)
!WRITE(6,'(a,100i3)') '  ipol:',(ipol(ion),ion=1,nion)

!     loop through ions and fill electron states successively

!nmaxact = nstate/nion+1                ! max. states per atom

IF(numspin==2) THEN
  ind=0
  DO ion=1,nion
    DO iki=1,nmxst(ion)
        ind=ind+1
        ispin(ind) = 2-MOD(iki,2)
        IF (ipol(ion).eq.-1) ispin(ind) = 2-mod(ispin(ind)+1,2)
    END DO
  END DO
END IF


numstate = 0
DO ion=1,nion
  DO natlevel=1,nmxst(ion)
    IF(numstate == ksttot) GO TO 99
    numstate = 1 + numstate
    IF(numspin==2) THEN
      nactst(ion) = nactst(ion)+MOD(natlevel,2)
      write(6,*) 'ion,natlev,numst,ispin', &
                  ion,natlevel,numstate,ispin(numstate)
      write(7,*) 'ksttot,nstate,ion,natlev,numst,ispin', &
                  ksttot,nstate,ion,natlevel,numstate,ispin(numstate)
    ELSE
      nactst(ion) = nactst(ion)+1
    END IF
    IF(nactst(ion) > 10) STOP 'GENERMOWF: too high atomic level'
    IF(numstate > nstate) GO TO 99
!                                                 select nodes
    nstx = nnodes(initord(1,ion),nactst(ion))
    nsty = nnodes(initord(2,ion),nactst(ion))
    nstz = nnodes(initord(3,ion),nactst(ion))
!                                                 compute raw wf
    ind=0
    DO iz=1,nz2
      z1=((iz-nzsh)*dz-cz(ion))/radini(ion)
      DO iy=1,ny2
        y1=((iy-nysh)*dy-cy(ion))/radini(ion)
        DO ix=1,nx2
          x1=((ix-nxsh)*dx-cx(ion))/radini(ion)
          ind = ind+1
          psiom(ind,numstate) = (x1**nstx-0.5*MAX(0,nstx-1))  &
              *(y1**nsty-0.5*MAX(0,nsty-1)) *(z1**nstz-0.5*MAX(0,nstz-1))  &
              *EXP(-(x1*x1+y1*y1+z1*z1)*0.5)
        END DO
      END DO
    END DO
    WRITE(6,'(a,i3,a,a,5i3,1pg13.5)')  &
        ' electron state nr.',numstate,': ',  &
        ' ion,nstx,nsty,nstz=',ion,nstx,nsty,nstz,ispin(numstate), &
         SUM(psiom(:,numstate)**2)*dvol
  END DO
END DO
99   CONTINUE

DEALLOCATE(nactst)


!     Schmidt ortho-normalisation

!CALL schmidt(psiom)


DO nbe=1,nstate
  
  DO in=1,nbe
    IF((ispin(nbe) == ispin(in))) THEN
!      sum=0D0
!      DO i=1,nxyz
!        sum=sum+psiom(i,nbe)*psiom(i,in)
!      END DO
!      sum=sum*dvol
      cs = SUM(psiom(:,nbe)*psiom(:,in))*dvol      
      IF(in == nbe) THEN
        cs = 1D0/SQRT(cs)
!        DO i=1,nxyz
          psiom(:,nbe)=psiom(:,nbe)*cs
!        END DO
      ELSE
!        DO  i=1,nxyz
          psiom(:,nbe)=psiom(:,nbe)-cs*psiom(:,in)
!        END DO
      END IF
    END IF
  END DO
END DO

DO n=1,numstate
   WRITE(*,*) ' INITLCGO after Schmidt: state nr. norm=', &
   n,SUM(psiom(:,n)**2)*dvol
END DO
CALL flush(6)


RETURN
END SUBROUTINE genermowf

