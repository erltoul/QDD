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

#include "define.h"
 
!-----calcrhot------------------------------------------------------

SUBROUTINE calcrhot(rho,q0)

!     density 'rho' for complex or real wavefunctions 'q0'

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

CLASS(*),TARGET, INTENT(IN) :: q0(kdfull2,kstate)
REAL(DP),TARGET, INTENT(OUT) :: rho(2*kdfull2)

!REAL(DP),DIMENSION(:),POINTER :: rh


!-----------------------------------------------------------------

!rh => rho

rho=0D0

DO nb=1,nstate
  ishift = (ispin(nrel2abs(nb))-1)*nxyz ! store spin=2 in upper block
  SELECT TYPE(q0)
   TYPE IS(REAL(DP)) 
    DO ind=1,nxyz
      rho(ind+ishift)=rho(ind+ishift)+ occup(nb)*q0(ind,nb)**2
    END DO
   TYPE IS(COMPLEX(DP)) 
    DO ind=1,nxyz
      rho(ind+ishift)=rho(ind+ishift)+ occup(nb)*&
       (REAL(q0(ind,nb))**2+AIMAG(q0(ind,nb))**2)
    END DO
  END SELECT 
END DO



RETURN
END SUBROUTINE calcrhot
