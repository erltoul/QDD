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

!    Lennard-Jones for Ar
PARAMETER (v0_lj=4*7.61E-4)   /* 4 * 120 kelvin in units of [ry]    */
PARAMETER (sigma=6.54)        /* 0.34 nm  in units of [bohr]        */

DO i=20,140,2
  x = i*0.1
  pot = v0_lj*((sigma/x)**12.0-((sigma/x)**6.0))
  WRITE(6,'(1x,2f6.1,1pg13.5)')    x,x/2.0,pot
END DO
STOP
END
