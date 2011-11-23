#include "define.h"
 
#if(parano)


!     ******************************

SUBROUTINE schmidt(q0)

!     ******************************


!     serial version of Schmidt orthogonalization

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: q0(kdfull2,kstate)

REAL(DP) :: eord(kstate)
INTEGER :: isort(kstate)

LOGICAL, PARAMETER :: tord=.false.

!*********************************************************

!     sort the s.p. energies

DO n=1,nstate
  eord(n)   = amoy(n)
  isort(n) = n
END DO
!      write(6,'(a,10(/5g12.4))')
!     &  ' spe before:',(eord(n),n=1,nstate)
IF(tord) THEN
  DO n=1,nstate
    emin   = 1D32
    imin   = n
    DO i=n,nstate
      IF(eord(i) < emin) THEN
        emin   = eord(i)
        imin   = i
      END IF
    END DO
    isav         = isort(imin)
    eord(imin)   = eord(n)
    isort(imin)  = isort(n)
    eord(n)      = emin
    isort(n)     = isav
  END DO
!        write(6,'(a,10(/5g12.4))')
!     &    ' spe after:',(eord(n),n=1,nstate)
END IF

!     Schmidt ortho-normalisation

DO nbes=1,nstate
  nbe = isort(nbes)
  
  DO ncs=1,nstate
    ncc = isort(ncs)
    IF((ispin(nbe) == ispin(ncc)) .AND. ncc <= nbe) THEN
      cs=0.0
      DO i=1,nxyz
        cs=cs+q0(i,nbe)*q0(i,ncc)
      END DO
      cs=cs*dvol
      IF(ncc == nbe) THEN
        DO i=1,nxyz
          q0(i,nbe)=q0(i,nbe)/SQRT(cs)
        END DO
      ELSE
        DO  i=1,nxyz
          q0(i,nbe)=q0(i,nbe)-cs*q0(i,ncc)
        END DO
      END IF
    END IF
  END DO
END DO

RETURN
END SUBROUTINE schmidt

#endif


#if(parayes)


!     ******************************

SUBROUTINE schmidt(q0)

!     ******************************


!     parallel version of Schmidt orthogonalization

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
REAL(DP), INTENT(IN OUT)         :: q0(kdfull2,kstate)


INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)


REAL(DP), ALLOCATABLE :: q3(:)         ! workspace for foreign w.f.
INTEGER :: isp3

!     fields for some recursive algorithm, this is a dummy
!     operation which serves to inhibit loop unrolling

!old      integer ndumnode(-1:knodem)
!old      integer ndumstate(0:kstate)

LOGICAL :: tsync
LOGICAL, PARAMETER :: ttest=.false.

!*********************************************************

ALLOCATE(q3(kdfull2))

!     determine own node 'myn'

CALL  mpi_comm_rank(mpi_comm_world,myn,icode)
CALL mpi_barrier (mpi_comm_world, mpi_ierror)

!     loop over nodes and states: w.f. to be normalized

!     orthogonalize in each next node with wfs sent from previous nodes

IF(ttest) WRITE(*,*) 'SCHMID: myn=',myn
DO nod2=0,myn-1
  DO nbe2=1,nstate_node(nod2)
    IF(ttest) THEN
      WRITE(6,'(a,5i5)') 'SCHMID: wait myn,nbe2,nod2,tag=',  &
          myn,nbe2,nod2,2*nbe2+1
      WRITE(7,'(a,5i5)') 'SCHMID: wait myn,nbe2,nod2,tag=',  &
          myn,nbe2,nod2,2*nbe2+1
      CALL flush()
    END IF
    CALL mpi_recv(q3,kdfull2, mpi_double_precision,nod2,nbe2*2+1,  &
        mpi_comm_world,is,ic)
!test        nods = nod2
!test        call mpi_recv(isp3,1,mpi_integer,nod2,
!test     &                nbe2*2+2,mpi_comm_world,is,ic)
    isp3 = ispin_node(nbe2,nod2)
    IF(ttest) THEN
      WRITE(6,'(a,5i5)') 'SCHMID: received myn,nbe2,nod2,tag=',  &
          myn,nbe2,nod2,nbe2*2+1
      WRITE(7,*) 'SCHMID: received myn,nbe2,nod2=',myn,nbe2,nod2
      CALL flush()
    END IF
    DO nbe1=1,nstate_node(myn)
      IF(ispin(nrel2abs(nbe1)) == isp3) THEN
        CALL orthogonalize(q0(1,nbe1),q3)
        IF(ttest) THEN
          WRITE(7,'(a,4i5,3i3)')  &
              ' orthog: myn,nbe1,nod2,nbe2,spin2,spin1=',  &
              myn,nbe1,nod2,nbe2,isp3, ispin(nrel2abs(nbe1)),ispin_node(nbe1,myn)
          WRITE(6,'(a,4i5,3i3)')  &
              ' orthog: myn,nbe1,nod2,nbe2,spin2,spin1=',  &
              myn,nbe1,nod2,nbe2,isp3, ispin(nrel2abs(nbe1)),ispin_node(nbe1,myn)
        END IF
      END IF
    END DO
  END DO
END DO

!     ortho-normalize within node assuming orthogonality on
!     previous nodes (established above)

DO nbe1=1,nstate_node(myn)
  DO nbe2=1,nbe1-1
    IF(ispin(nrel2abs(nbe1)) == ispin(nrel2abs(nbe2))) THEN
      CALL orthogonalize(q0(1,nbe1),q0(1,nbe2))
      IF(ttest) WRITE(*,*) 'SCHMID orth. myn,nbe1,nbe2=', myn,nbe1,nbe2
      IF(ttest) WRITE(7,*) 'SCHMID orth. myn,nbe1,nbe2=', myn,nbe1,nbe2
    END IF
  END DO
  CALL normalize(q0(1,nbe1))
  IF(ttest) WRITE(*,*) 'SCHMID: norm. myn,nbe1=',myn,nbe1
  IF(ttest) WRITE(7,*) 'SCHMID: norm. myn,nbe1=',myn,nbe1
  
!     distribute to next nodes
  
  DO nod2=myn+1,knode-1
    IF(ttest) THEN
      WRITE(6,'(a,5i5)') 'SCHMID: before send myn,nbe1,nod2,tag=',  &
          myn,nbe1,nod2,nbe1*2+1
      WRITE(7,'(a,5i5)') 'SCHMID: before send myn,nbe1,nod2,tag=',  &
          myn,nbe1,nod2,nbe1*2+1
      CALL flush()
    END IF
    CALL mpi_ssend(q0(1,nbe1),kdfull2,mpi_double_precision,  &
        nod2,nbe1*2+1,mpi_comm_world,ic)
!test          nods = nod2
!test          call mpi_send(ispin(nrel2abs(nbe1)),1,mpi_integer,
!test     &                 nod2,nbe1*2+2,mpi_comm_world,ic)
    IF(ttest) THEN
      WRITE(6,'(a,5i5)') 'SCHMID: sent myn,nbe1,nod2,tag=',  &
          myn,nbe1,nod2,nbe1*2+1
      WRITE(7,*) 'SCHMID: sent myn,nbe1,nod2=',myn,nbe1,nod2
      CALL flush()
    END IF
  END DO
  
END DO

!     synchronize by sending and receiving termination signal
!     from last node

tsync = .false.
IF(tsync) THEN
  IF(myn == knode-1) THEN
    iend = 1
    DO nod2=0,knode-2
      CALL mpi_bsend(iend,1,mpi_integer, nod2,1,mpi_comm_world,ic)
      IF(ttest) WRITE(*,*) ' terminator sent myn=',myn
    END DO
  ELSE
    CALL mpi_recv(iend,1,mpi_integer,knode-1,  &
        mpi_any_tag,mpi_comm_world,is,ic)
    IF(ttest) WRITE(*,*) ' terminator received myn,iend=',myn,iend
!test        write(7,*) ' SCHMID: terminator received'
  END IF
END IF

CALL mpi_barrier (mpi_comm_world, mpi_ierror)
DEALLOCATE(q3)

RETURN
END SUBROUTINE schmidt



#endif




!     ******************************

SUBROUTINE normalize(qact)

!     ******************************


!     normalizes real wavefunction on 'qact'

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)


REAL(DP), INTENT(IN OUT)                     :: qact(kdfull2)


!*********************************************************

!     compute norm

cs=0D0
DO i=1,nxyz
  cs=cs+qact(i)*qact(i)
END DO
cs=cs*dvol

!     normalize

cs = 1D0/SQRT(cs)
DO i=1,nxyz
  qact(i)=qact(i)*cs
END DO

RETURN

END SUBROUTINE normalize

!     ******************************

SUBROUTINE orthogonalize(qact,qorth)

!     ******************************


!     orthogonalizes real wavefunction 'qact' on  'qorth'

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: qact(kdfull2)
REAL(DP), INTENT(IN)                         :: qorth(kdfull2)



!*********************************************************

!     compute overlap

cs=0D0
DO i=1,nxyz
  cs=cs+qact(i)*qorth(i)
END DO
cs=cs*dvol

!     orthogonalize

DO  i=1,nxyz
  qact(i)=qact(i)-cs*qorth(i)
END DO

RETURN

END SUBROUTINE orthogonalize
