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

 

#ifdef REALSWITCH

MODULE twost_util

USE params
USE kinetic
USE orthmat
USE util, ONLY:wfovlp,project

INTEGER,SAVE,PRIVATE :: kdim

INTEGER,SAVE,PRIVATE :: ndim(2)


INTERFACE superpose_state
! automated choice between real, complex and complexsic versions of superpose_state
  MODULE PROCEDURE superpose_state_r, superpose_state_c, superpose_state_rc!, superpose_state_rc_dummy
END INTERFACE superpose_state

CONTAINS
!     ******************************

!----------------------------------------------------------------------
! REAL version
!----------------------------------------------------------------------
SUBROUTINE superpose_state_r(wfsup,coeff,q0,is)

!     Superposition to new state:
!       wfsup     new single-particle state
!       coeff     superposition coefficients
!       q0        set of s.p.states to be combined
!       is        spin of states

USE params
USE kinetic
IMPLICIT NONE

REAL(DP),INTENT(OUT)          :: wfsup(kdfull2)
REAL(DP),INTENT(IN)           :: coeff(kstate)
REAL(DP),INTENT(IN)           :: q0(kdfull2,kstate)
INTEGER,INTENT(IN)       ::is

INTEGER :: i, na, nas
!---------------------------------------------------------------------
wfsup(1:nxyz)=0.0_DP

DO nas=1,ndims(is)
  na = nas + (is-1)*ndims(1)
  DO i=1,nxyz
    wfsup(i) = wfsup(i) + coeff(nas)*q0(i,na) 
  END DO
END DO

RETURN
END SUBROUTINE superpose_state_r

!----------------------------------------------------------------------
! COMPLEX version
!----------------------------------------------------------------------
SUBROUTINE superpose_state_c(wfsup,coeff,q0,is)

!     Superposition to new state:
!       wfsup     new single-particle state
!       coeff     superposition coefficients
!       q0        set of s.p.states to be combined
!       is        spin of states

USE params
USE kinetic
IMPLICIT NONE

COMPLEX(DP),INTENT(OUT)       :: wfsup(kdfull2)
COMPLEX(DP),INTENT(IN)        :: coeff(kstate)
COMPLEX(DP),INTENT(IN)        :: q0(kdfull2,kstate)
INTEGER,INTENT(IN)       ::is

INTEGER :: i, na, nas
!---------------------------------------------------------------------

wfsup(1:nxyz)=(0.0_DP, 0.0_DP)

DO nas=1,ndims(is)
  na = nas + (is-1)*ndims(1)
  DO i=1,nxyz
    wfsup(i) = wfsup(i) + coeff(nas)*q0(i,na) 
  END DO
END DO

RETURN
END SUBROUTINE superpose_state_c

!----------------------------------------------------------------------
! cmplxsic version
!----------------------------------------------------------------------
SUBROUTINE superpose_state_rc(wfsup,coeff,q0,is)

!     Superposition to new state:
!       wfsup     new single-particle state
!       coeff     superposition coefficients
!       q0        set of s.p.states to be combined
!       is        spin of states

USE params
USE kinetic
IMPLICIT NONE

COMPLEX(DP),INTENT(OUT)       :: wfsup(kdfull2)
COMPLEX(DP),INTENT(IN)        :: coeff(kstate)
REAL(DP),INTENT(IN)           :: q0(kdfull2,kstate)
INTEGER,INTENT(IN)       ::is

INTEGER :: i, na, nas
!---------------------------------------------------------------------

wfsup(1:nxyz)=(0.0_DP, 0.0_DP)

DO nas=1,ndims(is)
  na = nas + (is-1)*ndims(1)
  DO i=1,nxyz
    wfsup(i) = wfsup(i) + coeff(nas)*q0(i,na) 
  END DO
END DO

RETURN
END SUBROUTINE superpose_state_rc


SUBROUTINE eval_unitrot(qact,qold)

! Computes the piece of rotation amongst occupied states contained
! in the tansformation from 'qold' to 'qact'. The resulting unitary
! transformation is transferred via 'wfrotate'.

USE params
USE util, ONLY:wfovlp
USE twost
USE orthmat, ONLY: matdorth
IMPLICIT NONE

COMPLEX(DP), INTENT(IN)   :: qact(kdfull2,kstate)
COMPLEX(DP), INTENT(IN)   :: qold(kdfull2,kstate)

COMPLEX(DP) :: ovl
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
          wfrotate(naeff,nbeff,is) = wfovlp(qold(:,na),qact(:,nb))
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
    wfrotate(1:ni,nb,is) = wfrotate(1:ni,nb,is)*CMPLX(1D0/SQRT(REAL(ovl,DP)),0D0,DP)
  END DO

! transpose for further use
  wfrotate(1:ni,1:ni,is) = CONJG(TRANSPOSE(wfrotate(1:ni,1:ni,is)))

  rmo = matdorth(wfrotate(:,:,is),kstate,ni)
  WRITE(*,*) 'is,unitarity wfrotate:',is,rmo
  IF(ABS(rmo)>1D-10) WRITE(*,*) ' WFROTATE:',wfrotate(1:ni,1:ni,is)

END DO

END SUBROUTINE eval_unitrot


END MODULE twost_util

#endif

