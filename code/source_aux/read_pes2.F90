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

        PROGRAM test
          IMPLICIT REAL*8 (A-H,O-Z)
          INTEGER,PARAMETER :: nlines=100000000
          INTEGER,PARAMETER :: nheader=2 ! commentary lines
          INTEGER,PARAMETER :: ncomment=9 ! commentary lines 
          INTEGER,PARAMETER :: np=18 ! the spacing between 0 to 180 degrees
          INTEGER,DIMENSION(0:np) :: rc
          REAL*8,DIMENSION(0:np) :: pes,pes2
          REAL*8,DIMENSION(0:np) :: radius
          CHARACTER :: symbol
          
          OPEN(9,FILE="c60_e136_nz1001.ekinet") ! raw data
          OPEN(10,FILE="pes_e136_final.c60") ! with two head lines
          OPEN(11,FILE="pes_e136_rescaling.ek") ! new ouput for ARPES
          
          DO i0=1,ncomment
             READ(9,*)
          END DO

          DO i1=0,np
             READ(9,*) symbol,&
             xx,xx,xx,rc(i1),radius(i1),xx,xx
          END DO          

          DO i=1,nheader
            READ(10,*)
          END DO
          DO j=1,nlines
          
            READ(10,*,err=100,end=100) omega,pes 
            DO k=0,np
              pes2(k)=pes(k)/(rc(k)*8.0) ! 1/(rc(k)*8.0) is the scaling factor
              IF((k==0).OR.(k==np)) pes2(k)=pes(k) ! no scaling for polar points
            END DO
            WRITE(11,'(1pg14.6,18(1pg14.4))') omega,pes2
          ENDDO

100     continue
        END
