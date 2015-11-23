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

!#include "define.h"
 
! Code converted using TO_F90 by Alan Miller
! Date: 2010-04-16  Time: 12:39:07

PROGRAM postrungeometry
!USE params

outnam='Na6_MgO'
nion=6

OPEN(30,STATUS='old',FILE='pposion.'//outnam)
OPEN(40,STATUS='unknown',FILE='pgeomion.postrun')

DO ir=1,999999
  
  READ(30,*)
  
  DO i=1,nion
    READ(30,*,ERR=19,END=19) tfs,cx(i),cy(i),cz(i),dummy
  END DO
  
  WRITE(6,*) 'tfs = ',tfs
  
  CALL getclustergeometry
  
  WRITE(40,'(12e15.5)') tfs,comx,comy,comz,rmsion,dmdistion,  &
      qtion(1,1),qtion(2,2),qtion(3,3), qtion(1,2),qtion(1,3),qtion(2,3)
  
  
  
END DO

CLOSE(40)
CLOSE(30)



19   STOP
END PROGRAM postrungeometry
