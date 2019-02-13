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

!      integrate 1/r in rectangular box
 
IMPLICIT DOUBLE PRECISION (a-h,o-z)

DO iter=2,11
  ngrid = 2**iter
  
  acc = 0D0
  dx = 1D0/ngrid
  DO ix=1,ngrid
    x2 = (ix-0.5D0)*dx
    x2 = x2*x2
    DO iy=1,ngrid
      y2 = (iy-0.5D0)*dx
      y2 = y2*y2
      DO iz=1,ngrid
        z = (iz-0.5D0)*dx
        acc = acc + 1D0/SQRT(x2+y2+z*z)
      END DO
    END DO
  END DO
  WRITE(*,*) ' ngrid,integral=',ngrid,acc*dx**3
  
END DO
STOP
END
