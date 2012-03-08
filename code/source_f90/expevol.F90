#include "define.h"
 

!-----tstep_exp-------------------------------------------------------------

SUBROUTINE tstep_exp(q0,aloc,rho,it,qwork)

!     one electronic time step by TV splitting method.

!     'itsub' indicates the number of subiteration before
!     next analyzing step (routine 'info').
!     'itsub=1' is the first call in a series and
!     'itsub=ipasinf' is the last call.
!g     Now, 'ipasinf' gives the step of computation of the electronic force
!g     on the cluster ion. Here it is only important with option 'iffastpropag'.

!     For pure electronic propagation one has the option to
!     reduce the number of local unitary steps. The last half-step
!     is omitted (exept in case of the last call 'itsub=ipasinf')
!     and the first local step is doubled instead (except for the
!     first call 'itsub=1'). This option is switched on by the
!     compile time switch 'fastpropag'.
!g     Option 'iffastpropag' must be reserved for pure electronic dynamics,
!g     it is not adapted to the new electronic/ionic dynamics.

!g     For pure electronic propagation with 'iffastpropag', 'ipasinf' is then
!g     the number of iteration in which the number of local unitary steps is reduced

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN)                      :: q0(kdfull2,kstate)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)
!REAL(DP), INTENT(IN OUT)                     :: akv(kdfull2)
REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
INTEGER, INTENT(IN)                      :: it
COMPLEX(DP), INTENT(OUT)                     :: qwork(kdfull2,kstate)

COMPLEX(DP) :: q1(kdfull2)





!----------------------------------------------------------------------


#if(parayes)
STOP 'exponential evolution not yet parallelized'
#else
myn = 0
#endif

STOP

IF(nc+NE+nk > 0) STOP 'TSTEP_EXP not apppropriate for rare gas'

!      write(*,*) 'entering TSTEP_EXP'


!     one half time step to define new mean field
!     use exponential evolution to second order

IF(ifsicp==5) psisavex = q0

DO nb=1,nstate
  DO ind=1,nxyz
    qwork(ind,nb) = q0(ind,nb)
  END DO
  CALL exp_evol(qwork(1,nb),aloc,nb,4,dt1*0.5D0,q1)
END DO


!     compute mean field at half time step

CALL dyn_mfield(rho,aloc,qwork,dt1*0.5D0)


!     full time step to next wavefunctions
!     use exponential evolution to fourth order

itpri = MOD(it,ipasinf) + 1
DO nb=1,nstate
  CALL exp_evol(q0(1,nb),aloc,nb,4,dt1,q1)
END DO




!     possibly absorbing bounds    ?? really here ??


IF (ntref > 0 .AND. it > ntref) nabsorb = 0
!      tfs = tfs + (dt - dt1)*0.0484/(2.*ame)

!      if(itsub.ne.ipasinf) then
!         if(nabsorb.gt.0) then
!            iComeFromAbso=0
!            if (iangabso.ne.0) call calcrho(rho,q0)
!                                ! necessary for angular distribution
!            if (ispherAbso.eq.0) then
!               call abso(q0,it)
!            else
!               call spherAbso(q0,it)
!            endif
!         endif
!      endif



!     compute mean field at new time

CALL dyn_mfield(rho,aloc,q0,dt1)

RETURN
END SUBROUTINE tstep_exp

!-----exp_evol-------------------------------------------------------------

SUBROUTINE exp_evol(qact,aloc,nbe,norder,dtact,qwork)

!     Propagation of a wavefuntion by Taylor expanded exponential
!     evolution:
!       qact     = wavefunction on which H acts and resulting w.f.
!       aloc     = local potential for the actual spin component
!       ak       = kinetic energies in momentum space
!       nbe      = number of state
!       norder   = order of epxansion (4 recommended for full step))
!       dtact    = time step
!       qwork    = work space for complex wavefunction

!     Note: The propagation uses the action of the Hamiltonian
!           where the diagonal element (s.p.energy) is subtracted.
!           That diagonal element is evalutaed in the first order
!           call 'nterm=1'.

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN OUT)                  :: qact(kdfull2)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)
!REAL(DP), INTENT(IN OUT)                     :: akv(kdfull2)
INTEGER, INTENT(IN OUT)                  :: nbe
INTEGER, INTENT(IN)                      :: norder
REAL(DP), INTENT(IN OUT)                     :: dtact
COMPLEX(DP), INTENT(OUT)                     :: qwork(kdfull2)



COMPLEX(DP) :: dti,cfac
INTEGER :: ilocbas
!test      complex wfovlp,energexp

!----------------------------------------------------------------------

!      write(*,*) 'entering EXP_EVOL'

IF (ispin(nrel2abs(nbe)) == 1) THEN
  ilocbas = 1
ELSE IF (ispin(nrel2abs(nbe)) == 2) THEN
  ilocbas = nxyz+1
ELSE
  STOP " EXPEVOL: spin index must be 1 or 2"
END IF
dti = CMPLX(0D0,dtact,DP)
cfac = CMPLX(1D0,0D0,DP)
DO  i=1,nxyz
  qwork(i) = qact(i)
END DO
DO nterm=1,norder
  CALL hpsi(qwork,aloc(ilocbas),nbe,nterm)
  cfac = -dti/nterm*cfac
  DO  i=1,nxyz
    qact(i) = qact(i) + cfac*qwork(i)
  END DO
END DO
!      write(*,*) 'leave EXP_EVOL'

RETURN
END SUBROUTINE exp_evol


!-----hpsi  -------------------------------------------------------------

SUBROUTINE hpsi(qact,aloc,nbe,itpri)

!     Action of mean-field Hamiltonian on one s.p. wavefunction:
!       qact     = wavefunction on which H acts and resulting w.f.
!       aloc     = local potential for the actual spin component
!       ak       = kinetic energies in momentum space
!       nbe      = number of state
!       itpri    = switch for computing s.p. energies


USE params
USE kinetic
#if(twostsic)
USE twost
#endif
IMPLICIT REAL(DP) (A-H,O-Z)



COMPLEX(DP), INTENT(IN OUT)              :: qact(kdfull2)
REAL(DP), INTENT(IN)                     :: aloc(2*kdfull2)
!REAL(DP), INTENT(IN)                    :: akv(kdfull2)
INTEGER, INTENT(IN OUT)                  :: nbe
INTEGER, INTENT(IN OUT)                  :: itpri
COMPLEX(DP),ALLOCATABLE :: qex(:)

COMPLEX(DP) :: wfovlp
!                                   workspaces
COMPLEX(DP),ALLOCATABLE :: q1(:),q2(:)

LOGICAL :: tpri
#if(twostsic)
COMPLEX(DP) :: cf
#endif


!----------------------------------------------------------------------


!      write(*,*) 'entering HPSI'
tpri = itpri == 1


ALLOCATE(q1(kdfull2),q2(kdfull2))

!     action of kinetic energy

#if(gridfft)
CALL fftf(qact,q1)
DO  i=1,nxyz
  q1(i) = akv(i)*q1(i)
END DO
CALL fftback(q1,q2)
IF(tpri) ekinsp(nbe) = wfovlp(qact,q2)
#else
STOP ' HPSI not yet appropriate for finite differences'
#endif



!     action of potential and non-local PsP (optionally)

IF(ipsptyp == 1) THEN
  q1 = 0D0
  CALL nonlocalc(qact,q1,0)
  IF(tpri) enonlo(nbe) = wfovlp(qact,q1)
  DO  i=1,nxyz
    q1(i)=q1(i)+qact(i)*aloc(i)
  END DO
ELSE
  DO  i=1,nxyz
    q1(i)=qact(i)*aloc(i)
  END DO
END IF

IF(ifsicp==5) THEN
  ALLOCATE(qex(kdfull2))
  IF(tpri) epotbefore = wfovlp(qact,q1)
  CALL exchg(qact,qex,nbe)
  q1 = q1 + qex
!  IF(tpri) THEN
!    epotafter = wfovlp(qact,q1)
!    WRITE(6,'(a,i4,3(1pg13.5))') ' EXC: nbe,before,after,such=', &
!         nbe,epotbefore,epotafter,wfovlp(qact,qex)
!    CALL FLUSH(6)
!  END IF
  DEALLOCATE(qex)
END IF


!JM : subtract SIC potential for state NBE
#if(twostsic)
IF(ifsicp == 8) THEN
  
  is=ispin(nrel2abs(nbe))
  DO na=1,nstate
    IF(ispin(nrel2abs(na)) == is)THEN
      cf = wfovlp(psiut(1,na),qact)
      DO i=1,nxyz
        q1(i)=q1(i)-qnewut(i,na)*cf
      END DO
    END IF
  END DO
  
END IF
#endif
!JM


IF(tpri) THEN
  epotsp(nbe) = wfovlp(qact,q1)
  amoy(nbe) = ekinsp(nbe)+epotsp(nbe)
  q2 = q1+q2
  spvariance(nbe) = SQRT(REAL(wfovlp(q2,q2))-amoy(nbe)**2)
  qact = q2
ELSE
  qact = q1+q2
END IF


!espact =  amoy(nbe)


DEALLOCATE(q1,q2)



RETURN
END SUBROUTINE hpsi




SUBROUTINE hpsi_boostinv(qact,aloc,current,rho,nbe)

!     Action of boost-invariant Hamiltonian on one s.p. wavefunction:
!       qact     = wavefunction on which H acts and resulting w.f.
!       aloc     = local potential for the actual spin component
!       current  = average local momentum in x,y, and z-direction
!       nbe      = number of state
!     The routine requires that 'current' has been accumulated before.


USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)



COMPLEX(DP), INTENT(IN OUT)              :: qact(kdfull2)
REAL(DP), INTENT(IN)                     :: aloc(2*kdfull2)
REAL(DP), INTENT(IN)                     :: current(kdfull2,3)
REAL(DP), INTENT(IN)                     :: rho(2*kdfull2)
INTEGER, INTENT(IN OUT)                  :: nbe

COMPLEX(DP) :: wfovlp
!                                   workspaces
COMPLEX(DP),ALLOCATABLE :: q1(:),q2(:)

LOGICAL :: tpri
#if(twostsic)
COMPLEX(DP) :: cf
#endif


!----------------------------------------------------------------------

ALLOCATE(q1(kdfull2),q2(kdfull2))

q1=qact
CALL hpsi(qact,aloc,nbe,0)

CALL xgradient_rspace(q1,q2)
qact = qact - h2m*current(:,1)*q2 
q2 = current(:,1)*q1
CALL xgradient_rspace(q2,q2)
qact = qact - h2m*q2 


CALL ygradient_rspace(q1,q2)
qact = qact - h2m*current(:,2)*q2 
q2 = current(:,2)*q1
CALL ygradient_rspace(q2,q2)
qact = qact - h2m*q2 


CALL zgradient_rspace(q1,q2)
qact = qact - eye*h2m*current(:,3)*q2 
q2 = current(:,3)*q1
CALL zgradient_rspace(q2,q2)
qact = qact - eye* h2m*q2 

qact = qact + h2m* &
 (current(:,1)**2+current(:,2)**2+current(:,3)**2)

WRITE(*,*) ' HPSI_BOOSTINV: E_coll=',dvol*h2m*SUM(rho(:)* &
 (current(:,1)**2+current(:,2)**2+current(:,3)**2))

RETURN
END SUBROUTINE hpsi_boostinv