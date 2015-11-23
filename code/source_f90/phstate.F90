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
 
!-----phstate---------------------------------------------------------

SUBROUTINE phstate(psi)

!  Picks one particular 1ph state from the given wavefunctions 
!  'psi' and mixes this to a rotated ph configuration. The state
!  is rotated by 'phangle'.
!  The particle state is 'npstate' and the hole state 'nhstate'.
!  The old hole wavefunction is stored and returned on 'oldhole'.
!  The parameters 'phangle', 'npstate' and 'nhstate' are communicated
!  through module PARAMS.
!  'oldhole' is also communicated through PARAMS and allocated here.


USE params
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN OUT)   :: psi(kdfull2,kstate)

COMPLEX(DP) :: csave
INTEGER :: i

  IF(npstate > nstate) STOP ' PHSTATE: particle state out of range'
  IF(nhstate > nstate) STOP ' PHSTATE: hole state out of range'
  IF(occup(npstate) > 0.5D0) STOP 'PHSTATE: particle state already occupied'
  IF(occup(nhstate) < 0.5D0) STOP 'PHSTATE: hole state not occupied'

  ALLOCATE(oldhole(kdfull2))

  oldhole = psi(:,nhstate)
  ca = cos(phangle)
  sa = sin(phangle)
  psi(:,nhstate) = ca*psi(:,nhstate)+sa*psi(:,npstate)
  psi(:,npstate) = -sa*oldhole+ca*psi(:,npstate)


END SUBROUTINE phstate
