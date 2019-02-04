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

!!-----model-------------------------------------------------------------
 
! Code converted using TO_F90 by Alan Miller
! Date: 2010-04-16  Time: 12:37:44

SUBROUTINE model(a,nparmx,yfit,ndatmx)

REAL, INTENT(IN OUT)                     :: a(kparmx)
INTEGER, INTENT(IN)                      :: nparmx
REAL, INTENT(IN OUT)                     :: yfit(kdatmx)
INTEGER, INTENT(IN)                      :: ndatmx
IMPLICIT REAL*8(a-h,o-z)

!     compute chi**2 for Ar atom in SIC with local Goedecker

INTEGER, PARAMETER :: kparmx=5
INTEGER, PARAMETER :: kdatmx=10

CHARACTER (LEN=71) :: command

!-----------------------------------------------------------------------

IF(nparmx /= 3) STOP 'locgoe.f requires NPARMX=2'
IF(ndatmx /= 4) STOP 'locgoe.f requires NDATMX=4'
OPEN(13,STATUS='scratch',FORM='formatted')
!      write(6,'(a,4f6.1)')
!     &  'model: ',(A(i),i=1,nparmx)
WRITE(13,'(a,4f15.10)') 'fitpsp.bat ',0.8,(a(i),i=1,nparmx)
!      write(13,'(a,4f15.10)')
!     &  'fitpsp.bat ',0.6,A(1),0.0,A(2)
REWIND 13
READ(13,'(a)') command
CLOSE(13)
!      write(6,'(a)') command
CALL system(command)
OPEN(13,FILE='energies.res')
READ(13,'(1x)')
READ(13,'(1x)')
DO i=1,ndatmx
  READ(13,*) id1,id2,yfit(i)
END DO
CLOSE(13)
WRITE(6,'(a,4f7.2)') 'energies:',(yfit(i),i=1,ndatmx)
WRITE(7,'(a,4f7.2)') 'energies:',(yfit(i),i=1,ndatmx)
!      call system('tail -n 1 sumchi')

RETURN
END SUBROUTINE model

