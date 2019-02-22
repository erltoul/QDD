
!     **************************

REAL(DP) FUNCTION wfnorm_r(psi1,psi2)

!  Overlap of wavefunctions 'psi1' and 'psi2'.     REAL version

USE params
IMPLICIT NONE
COMPLEX(DP), INTENT(IN OUT)                  :: psi1(kdfull2)
COMPLEX(DP), INTENT(IN OUT)                  :: psi2(kdfull2)

INTEGER :: ii
REAL(DP) :: acc

acc = 0D0
DO ii=1,nxyz
  acc = REAL(psi1(ii),DP)*REAL(psi2(ii),DP) + AIMAG(psi1(ii))*AIMAG(psi2(ii)) + acc
END DO
wfnorm_r = acc*dvol
RETURN
END FUNCTION wfnorm_r
!     **************************

REAL(DP) FUNCTION wfnorm_i(psi1,psi2)

!  Overlap of wavefunctions 'psi1' and 'psi2'.     COMPLEX version

USE params
IMPLICIT NONE
COMPLEX(DP), INTENT(IN OUT)                  :: psi1(kdfull2)
COMPLEX(DP), INTENT(IN OUT)                  :: psi2(kdfull2)

INTEGER :: ii 
REAL(DP) :: acc
acc = 0D0
DO ii=1,nxyz
  acc = -REAL(psi1(ii),DP)*AIMAG(psi2(ii)) +AIMAG(psi1(ii))*REAL(psi2(ii),DP) + acc
END DO
wfnorm_i = acc*dvol
RETURN
END FUNCTION wfnorm_i





!**********************   from restart.F90  *******************

!     **************************

SUBROUTINE restherm()

!     **************************

USE params
USE kinetic
IMPLICIT NONE
INTEGER :: iact,ion

OPEN(UNIT=20,STATUS='unknown',FORM='unformatted',FILE='therm')
!     nxyz=nx2*ny2*nz2
READ(20) iact
IF(iact /= irest) THEN
  WRITE(7,*) 'iact=',iact
  STOP 'bad irest in for005dyn !'
END IF
DO ion=1,nion
  READ(20) cx(ion),cy(ion),cz(ion)
  READ(20) cpx(ion),cpy(ion),cpz(ion)
END DO
CLOSE(UNIT=20,STATUS='keep')
RETURN
END SUBROUTINE restherm






!     **************************

SUBROUTINE addcluster(psi,outna)

!     **************************

USE params
USE kinetic
IMPLICIT NONE

COMPLEX(DP), INTENT(IN OUT)              :: psi(kdfull2,kstate)
CHARACTER (LEN=13), INTENT(IN OUT)       :: outna

INTEGER :: i,iact,idum,ii,ion,k,nb,nclustt,niont,nspdwt,nstatet
OPEN(UNIT=ifile,STATUS='unknown',FORM='unformatted', FILE='save.'//outna)


!  read the iteration where the data has been saved last:


READ(ifile) iact,nstatet,nclustt,niont,nspdwt
DO i=1,nstatet
  READ(ifile) occup(i+nstate)
END DO


irest=iact


!  read wavefunctions:
IF(nclustt > 0)THEN
  DO nb=1,nstatet
    DO i=1,nxyz
      READ(ifile) psi(i,nb+nstate)
    END DO
  END DO
END IF

DO i=1,ksttot
  IF (i <= nstatet) THEN
    READ(ifile) ispin(i+nstate), nrel2abs(i+nstate),nabs2rel(i+nstate)
  ELSE
    READ(ifile) idum,idum,idum
  END IF
END DO




!  read protonic coordinates and momenta
WRITE(6,*) nion,niont
DO ion=1,niont
  READ(ifile) cx(ion+nion),cy(ion+nion),cz(ion+nion), np(ion+nion)
  READ(ifile) cpx(ion+nion),cpy(ion+nion),cpz(ion+nion), np(ion+nion)
END DO


#if(raregas)
IF (isurf /= 0) THEN

  DO i=1,nc
    READ(ifile) imobc(i)
    READ(ifile) xc(i),yc(i),zc(i)
    READ(ifile) pxc(i),pyc(i),pzc(i)
  END DO
  DO i=1,NE
    READ(ifile) imobe(i)
    READ(ifile) xe(i),ye(i),ze(i)
    READ(ifile) pxe(i),pye(i),pze(i)
  END DO
  DO i=1,nk
    READ(ifile) imobk(i)
    READ(ifile) xk(i),yk(i),zk(i)
    READ(ifile) pxk(i),pyk(i),pzk(i)
  END DO
  
  DO i=1,kdfull2
    READ(ifile) potfixedion(i)
  END DO
  
END IF
#endif


IF(nclustt > 0)THEN
  DO k=1,kmom
    READ(ifile) qe(k)
  END DO
  DO i=1,3
    READ(ifile) se(i)
  END DO
END IF


IF (nabsorb > 0) THEN
  DO ii=1,kdfull2
    READ(ifile) rhoabso(ii)
  END DO
END IF


CLOSE(UNIT=ifile,STATUS='keep')

nstate=nstate+nstatet
nclust=nclust+nclustt
nion=nion+niont
nspdw=nspdw+nspdwt
nion2=nion


RETURN
END SUBROUTINE addcluster


