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
SUBROUTINE unitary_gradstepr(symm_mat,vect,ndima,norder)
#else
SUBROUTINE unitary_gradstep(symm_mat,vect,ndima,norder)
#endif
!
!  Unitary extrapolation of the gradient step.
!    symm_mat   = matrix of symmetry condition, rescaled by stepsize
!    vect       = matrix of state vectors to be propagated
!    ndima      = actual dimension of matrix
!    norder     = order of the Taylor expansion of the exponential
!
USE params
!USE twostr
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER,PARAMETER :: kdim=kstate

INTEGER,INTENT(IN) :: norder,ndima

#ifdef REALSWITCH
REAL(DP),INTENT(IN) :: symm_mat(kdim,kdim)
REAL(DP),INTENT(INOUT) :: vect(kdim,kdim)
REAL(DP) :: vectwork(kdim,kdim)
REAL(DP) :: symm_acc(kdim,kdim)
REAL(DP) :: sumv
REAL(DP) :: w(kdim)             ! workspace
#else
COMPLEX(DP),INTENT(IN) :: symm_mat(kdim,kdim)
COMPLEX(DP),INTENT(INOUT) :: vect(kdim,kdim)
COMPLEX(DP) :: vectwork(kdim,kdim)
COMPLEX(DP) :: symm_acc(kdim,kdim)
COMPLEX(DP) :: sumv
COMPLEX(DP) :: w(kdim)             ! workspace
#endif

INTEGER :: i,j,k,n      ! counters in DO loops
REAL(DP) :: facinv

DO i=1,ndima; DO j=1,ndima
  vectwork(i,j) = vect(i,j)
END DO; END DO

facinv = 1D0
DO n=1,norder

  facinv = facinv/n

  DO k=1,ndima
    w = 0D0
    DO i=1,ndima
      DO j=1,ndima
        w(i) = w(i) - symm_mat(i,j)*vectwork(k,j)
      END DO
    END DO
    DO i=1,ndima
      vectwork(k,i) = w(i)*facinv
    END DO
  END DO

  DO i=1,ndima; DO j=1,ndima
    vect(i,j) = vect(i,j)+vectwork(i,j)
  END DO; END DO

END DO

RETURN
END
