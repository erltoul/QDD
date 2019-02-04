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

PARAMETER (pi=3.1415926)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2010-04-16  Time: 12:34:36
 
DO i=1,999
  READ(5,*,ERR=99,END=99) n,xx,yy,zz
  rms   = (xx+yy+zz)
  alpha = (zz+zz-xx-yy)*SQRT(5.0/(16.0*pi))/(rms)
  beta  = (xx-yy)*SQRT(15.0/(8.0*pi))/(rms)
  beta2 = SQRT(alpha*alpha+2.0*beta*beta)
  gamma = ATAN(1.4142136*beta/alpha)*180.0/pi
  IF(alpha < 0.0) gamma=180.0+gamma
  gamma = ABS(gamma)
  gamma = MOD(gamma,120.0)
  IF(gamma > 60.0 .AND. gamma <= 120.0) gamma = 120.0-gamma
  rms   = SQRT(rms)
  WRITE(6,'(1x,i5,f8.2,f8.3,f8.2)')  n,rms,beta2,gamma
END DO
99   CONTINUE
STOP
END SUBROUTINE pois

