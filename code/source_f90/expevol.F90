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
!     is omitted (except in case of the last call 'itsub=ipasinf')
!     and the first local step is doubled instead (except for the
!     first call 'itsub=1'). This option is switched on by the
!     compile time switch 'fastpropag'.
!g     Option 'iffastpropag' must be reserved for pure electronic dynamics,
!g     it is not adapted to the new electronic/ionic dynamics.

!g     For pure electronic propagation with 'iffastpropag', 'ipasinf' is then
!g     the number of iteration in which the number of local unitary steps is reduced

USE params
USE twost, ONLY:tnearest
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN)                      :: q0(kdfull2,kstate)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)
!REAL(DP), INTENT(IN OUT)                     :: akv(kdfull2)
REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
INTEGER, INTENT(IN)                      :: it
COMPLEX(DP), INTENT(OUT)                     :: qwork(kdfull2,kstate)

COMPLEX(DP) :: q1(kdfull2)

! The parameter 'tnorotate' de-activates the subtraction of the
! Lagrangian matrix in the SIC step. The version of exponential
! evolution with subtraction of the Lagrangian matrix is found
! in 'exp_evolp'. The strategy needs yet testing. It was not
! really beneficial so far.
LOGICAL,PARAMETER :: tnorotate=.true.


!----------------------------------------------------------------------


#if(parayes)
STOP 'exponential evolution not yet parallelized'
#else
myn = 0
#endif

!STOP

#if(raregas)
IF(nc+NE+nk > 0) STOP 'TSTEP_EXP not appropriate for rare gas'
#endif

!      write(*,*) 'entering TSTEP_EXP'


!     one half time step to define new mean field
!     use exponential evolution to second order

IF(ifsicp==5) psisavex = q0

IF(tnorotate .OR. ifsicp .NE. 8) THEN
  DO nb=1,nstate
    qwork(:,nb) = q0(:,nb)
    CALL exp_evol(qwork(1,nb),aloc,nb,4,dt1*0.5D0,q1)
  END DO
ELSE
#if(twostsic)
  qwork = q0
  CALL exp_evolp(qwork,aloc,4,dt1/2D0,q1,q0)
#else
  STOP " IFSICP==8 reqires compilation with option twostsic"
#endif
END IF

#if(twostsic)
IF(tnearest) CALL eval_unitrot(qwork,q0)
#endif


!     compute mean field at half time step

CALL dyn_mfield(rho,aloc,qwork,dt1*0.5D0)


!     full time step to next wavefunctions
!     use exponential evolution to fourth order

itpri = MOD(it,ipasinf) + 1
IF(tnorotate .OR. ifsicp .NE. 8) THEN
  DO nb=1,nstate
    CALL exp_evol(q0(1,nb),aloc,nb,4,dt1,q1)
  END DO
ELSE
#if(twostsic)
  qwork = q0
  CALL exp_evolp(q0,aloc,4,dt1,q1,qwork)
#endif
END IF

#if(twostsic)
IF(tnearest) CALL eval_unitrot(q0,qwork)
#endif




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

!     Propagation of a wavefunction by Taylor expanded exponential
!     evolution:
!       qact     = wavefunction on which H acts and resulting w.f.
!       aloc     = local potential for the actual spin component
!       ak       = kinetic energies in momentum space
!       nbe      = number of state
!       norder   = order of epxansion (4 recommended for full step))
!       dtact    = time step

!     Note: The propagation uses the action of the Hamiltonian
!           where the diagonal element (s.p.energy) is subtracted.
!           That diagonal element is evaluated in the first order
!           call 'nterm=1'.

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN OUT)              :: qact(kdfull2)
REAL(DP), INTENT(IN OUT)                 :: aloc(2*kdfull2)
!REAL(DP), INTENT(IN OUT                 :: akv(kdfull2)
INTEGER, INTENT(IN)                      :: nbe
INTEGER, INTENT(IN)                      :: norder
REAL(DP), INTENT(IN OUT)                 :: dtact
COMPLEX(DP), INTENT(OUT)                 :: qwork(kdfull2)


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

#if(twostsic)
!-----exp_evolp-------------------------------------------------------------

SUBROUTINE exp_evolp(qact,aloc,norder,dtact,qwork,psi)

!     Propagation of a wavefuntion by Taylor expanded exponential
!     evolution:
!       qact     = wavefunction on which H acts and resulting w.f.
!       aloc     = local potential for the actual spin component
!       ak       = kinetic energies in momentum space
!       norder   = order of epxansion (4 recommended for full step))
!       dtact    = time step
!       qwork    = work space for complex wavefunction
!       psi      = set of wavefunctions before the step

!     Note: The propagation uses the action of the Hamiltonian
!           where the diagonal element (s.p.energy) is subtracted.
!           That diagonal element is evalutaed in the first order
!           call 'nterm=1'.

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN OUT)              :: qact(kdfull2,kstate)
REAL(DP), INTENT(IN OUT)                 :: aloc(2*kdfull2)
INTEGER, INTENT(IN)                      :: norder
REAL(DP), INTENT(IN)                     :: dtact
COMPLEX(DP), INTENT(OUT)                 :: qwork(kdfull2)
COMPLEX(DP), INTENT(IN)                  :: psi(kdfull2,kstate)

COMPLEX(DP),ALLOCATABLE :: chmatrix(:,:)

COMPLEX(DP) :: dti,cfac,cacc(kstate),wfovlp
INTEGER :: ilocbas
INTEGER :: nbe
!test      complex wfovlp,energexp

!----------------------------------------------------------------------

!write(*,*) 'entering EXP_EVOLP'

ALLOCATE(chmatrix(kstate,kstate))

dti = CMPLX(0D0,dtact,DP)

! compute H-matrix, store h*psi wavefunctions
DO nbe=1,nstate
  ilocbas = 1 + (ispin(nrel2abs(nbe))-1)*nxyz
  CALL hpsi(qact(1,nbe),aloc(ilocbas),nbe,1)
  DO nc=1,nstate
    IF(ispin(nrel2abs(nbe)) == ispin(nrel2abs(nc))) THEN
      chmatrix(nc,nbe) = wfovlp(psi(1,nc),qact(1,nbe))
    ELSE
      chmatrix(nc,nbe) = CMPLX(0D0,0D0)
    END IF
  END DO
END DO
! symmetrize H-matrix
DO nbe=1,nstate; DO nc=1,nbe-1
  chmatrix(nc,nbe) = (chmatrix(nc,nbe)+CONJG(chmatrix(nbe,nc)))/2D0
  chmatrix(nbe,nc) = CONJG(chmatrix(nc,nbe))
END DO; END DO

! now the Taylor expansion (recycle stored h*psi in first step)

!WRITE(*,*) 'CHMATRIX:'
!DO nbe=1,nstate
!  WRITE(*,*) chmatrix(:,nbe)
!END DO
DO nbe=1,nstate
  ilocbas = 1 + (ispin(nrel2abs(nbe))-1)*nxyz
  cfac = CMPLX(1D0,0D0,DP)
  DO nterm=1,norder
    ! H*psi and prepare subtraction factors
    IF(nterm==1) THEN
      qwork(:) = qact(:,nbe)
      qact(:,nbe) = psi(:,nbe)    ! restore original wavefunctions
      cacc(:) = chmatrix(:,nbe)
    ELSE
      DO nc=1,nstate
        IF(ispin(nrel2abs(nbe)) == ispin(nrel2abs(nc))) THEN
          cacc(nc) = CMPLX(0D0,0D0)
          DO na=1,nstate
            IF(ispin(nrel2abs(na)) == ispin(nrel2abs(nc))) THEN
              cacc(nc) = cacc(nc) + chmatrix(nc,na)*wfovlp(psi(1,na),qwork)
!              WRITE(*,*) ' NBE,NC,NA,ovlp:',nbe,nc,na,wfovlp(psi(1,na),qwork)
            END IF
          END DO
        END IF
      END DO
      CALL hpsi(qwork,aloc(ilocbas),nbe,2)
    END IF
    ! project H-matrix
!    WRITE(*,*) ' NBE,cacc:',nbe,cacc(1:nstate)
    DO nc=1,nstate
      IF(ispin(nrel2abs(nbe)) == ispin(nrel2abs(nc))) THEN
        qwork(:) = qwork(:) - psi(:,nc)*cacc(nc)
      END IF
    END DO
    ! accumulate to Taylor series
    cfac = -dti/nterm*cfac
    qact(:,nbe) = qact(:,nbe) + cfac*qwork(:)
  END DO
END DO

DEALLOCATE(chmatrix)

!write(*,*) 'leave EXP_EVOLP'

RETURN
END SUBROUTINE exp_evolp
!-----eval_unitrot-------------------------------------------------------------

SUBROUTINE eval_unitrot(qact,qold)

! Computes the piece of rotation amongst occupied states contained
! in the tansformation from 'qold' to 'qact'. The resulting unitary
! transformation is transferred via 'wfrotate'.

USE params
USE twost
IMPLICIT NONE

COMPLEX(DP), INTENT(IN)              :: qact(kdfull2,kstate)
COMPLEX(DP), INTENT(IN)              :: qold(kdfull2,kstate)

COMPLEX(DP) :: wfovlp,ovl
REAL(DP) :: rmo
INTEGER :: is,ni,na,nb,naeff,nbeff

!----------------------------------------------------------------------

wfrotate = 0D0
DO is=1,2
  ni = ndims(is)

! evaluate matrix from overlaps
  DO na=1,nstate
    IF(ispin(nrel2abs(na))==is) THEN
      naeff = na - (is-1)*ndims(1)
      DO nb=1,nstate
        IF(ispin(nrel2abs(nb)) == ispin(nrel2abs(na))) THEN
          nbeff = nb - (is-1)*ndims(1)
          wfrotate(naeff,nbeff,is) = wfovlp(qold(1,na),qact(1,nb))
        END IF
      END DO
    END IF
  END DO

! ortho-normalize
  DO nb=1,ni
    DO na=1,nb-1
      ovl=SUM(CONJG(wfrotate(1:ni,na,is))*wfrotate(1:ni,nb,is))
      wfrotate(:,nb,is) = wfrotate(:,nb,is)-wfrotate(:,na,is)*ovl
    END DO
    ovl=SUM(CONJG(wfrotate(1:ni,nb,is))*wfrotate(1:ni,nb,is))
    wfrotate(1:ni,nb,is) = wfrotate(1:ni,nb,is)*CMPLX(1D0/SQRT(REAL(ovl,DP)),0D0)
  END DO

! transpose for further use
  wfrotate(1:ni,1:ni,is) = CONJG(TRANSPOSE(wfrotate(1:ni,1:ni,is)))

  rmo = matdorth(wfrotate(1,1,is),kstate,ni)
  WRITE(*,*) 'is,unitarity wfrotate:',is,rmo
  IF(ABS(rmo)>1D-10) WRITE(*,*) ' WFROTATE:',wfrotate(1:ni,1:ni,is)

END DO

END SUBROUTINE eval_unitrot

#endif

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
