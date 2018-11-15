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

 
!------------------------------------------------------------

SUBROUTINE makemadelung
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
!     get Madelung potential from data file of previous
!     run with iOnTheFly = 1

OPEN(131,STATUS='old',FILE='phim.par')
!      read(31,*) ndfull2

IF (ndfull2 /= kdfull2) THEN
  STOP 'grid for Madelung potential not compatible with settings'
END IF


DO i=1,kdfull2
  READ (131,*) dummy,dummy,dummy,phim(i),dummy,dummy
!       write(6,'(3f10.2,3e15.5)') phim(i)
!      write(6,*) i,phim(i)
END DO

CLOSE(131)

END SUBROUTINE makemadelung
!------------------------------------------------------------


!------------------------------------------------------------

SUBROUTINE maketestpot1
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

ddx = pi/a/5
ddy=ddx
ddz=ddy

ind=0

DO iz=minz,maxz
  z1=(iz-nz)*dz
  DO iy=miny,maxy
    y1=(iy-ny)*dy
    DO ix=minx,maxx
      x1=(ix-nx)*dx
      
      ind = ind + 1
      phim(ind) = - 0.7E0*SIN(ddx*x1)*SIN(ddy*y1)* SIN(ddz*z1)     +1
!                  phim(ind)=-(  ddx*x1**2+ddy*y1**2+ddz*z1**2)
      
      
      
    END DO
  END DO
END DO

!          call priphim


END SUBROUTINE maketestpot1
!------------------------------------------------------------



!------------------------------------------------------------

SUBROUTINE buildmadelung(iside,ilayers)
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

WRITE(6,*) 'Constructing Madelung potential.'

!         call  init


CALL buildlatt(iside,ilayers)

nions = nk + nc

WRITE(6,*) 'Calculating fixed Madelung potential'

CALL calcpot

CALL priphim


END SUBROUTINE buildmadelung
!------------------------------------------------------------


!------------------------------------------------------------

SUBROUTINE calcpot
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)


ind=0
DO iz=minz,maxz
  DO iy=miny,maxy
    DO ix=minx,maxx
      
      
      ind = ind + 1
      
      phim(ind) = 0.0D0
    END DO
  END DO
END DO

chganion = chgcore + chgval




DO i=1,nc
  
  xi = xc(i)
  yi = yc(i)
  zi = zc(i)
  
  
  IF (ABS(xi) > r2bx .OR. ABS(yi) > r2by .OR. ABS(zi) > r2bz) THEN
!       Ion außerhalb der Simulationsbox
    
    ind=0
    DO iz=minz,maxz
      z1=(iz-nz)*dz
      DO iy=miny,maxy
        y1=(iy-ny)*dy
        DO ix=minx,maxx
          x1=(ix-nx)*dx
          
          ind=ind+1
          r = SQRT((xi-x1)**2 + (yi-y1)**2 + (zi-z1)**2)
          IF (r <= 1E-5) THEN
            r=crcut
          END IF
          
          phim(ind) = phim(ind) + chganion/r
          phimv(ind) = phimv(ind) + chganion/r
          
        END DO
      END DO
    END DO
    
    
  ELSE
! Ion in der Simulationsbox
    
    ind=0
    DO iz=minz,maxz
      z1=(iz-nz)*dz
      DO iy=miny,maxy
        y1=(iy-ny)*dy
        DO ix=minx,maxx
          x1=(ix-nx)*dx
          
          ind=ind+1
          
          r = SQRT((xi-x1)**2 + (yi-y1)**2 + (zi-z1)**2)
          IF (r <= 1E-5) THEN
            r=crcut
          END IF
          
          phimv(ind) = phimv(ind) + chganion/r
          phimd(ind) = phimd(ind) + chganion/r
          
        END DO
      END DO
    END DO
    
    
    
    
    
  END IF
  
END DO








DO i=1,nk
  
  xi = xk(i)
  yi = yk(i)
  zi = zk(i)
  
  
  IF (ABS(xi) > r2bx .OR. ABS(yi) > r2by .OR. ABS(zi) > r2bz) THEN
    
    ind=0
    DO iz=minz,maxz
      z1=(iz-nz)*dz
      DO iy=miny,maxy
        y1=(iy-ny)*dy
        DO ix=minx,maxx
          x1=(ix-nx)*dx
          
          ind=ind+1
          
          
          r = SQRT((xi-x1)**2 + (yi-y1)**2 + (zi-z1)**2)
          IF (r <= 1E-5) THEN
            r=crcut
          END IF
          
          phim(ind) = phim(ind) + chgkat/r
          phimv(ind) = phimv(ind) + chgkat/r
          
        END DO
      END DO
    END DO
    
    
  ELSE
    
    ind=0
    DO iz=minz,maxz
      z1=(iz-nz)*dz
      DO iy=miny,maxy
        y1=(iy-ny)*dy
        DO ix=minx,maxx
          x1=(ix-nx)*dx
          
          ind=ind+1
          
          
          r = SQRT((xi-x1)**2 + (yi-y1)**2 + (zi-z1)**2)
          IF (r <= 1E-5) THEN
            r=crcut
          END IF
          
          phimv(ind) = phimv(ind) + chgkat/r
          phimd(ind) = phimd(ind) + chgkat/r
          
        END DO
      END DO
    END DO
    
    
    
    
  END IF
  
END DO


END SUBROUTINE calcpot
!------------------------------------------------------------

!------------------------------------------------------------
!      subroutine initMadelung
!------------------------------------------------------------
!#include "surf.inc"


!       ind=0
!       do iz=minz,maxz
!         z1=(iz-nz)*dz
!          do iy=miny,maxy
!           y1=(iy-ny)*dy
!           do ix=minx,maxx
!               x1=(ix-nx)*dx

!                  ind=ind+1

!           phim(ind)= 0.0

!            enddo
!           enddo
!          enddo


!          ndfull2 = ind

!          write (6,*) ndfull2, kdfull2

!          if (ndfull2.ne.kdfull2) then
!             stop 'unequal'
!          endif


!      end
!------------------------------------------------------------




!------------------------------------------------------------

SUBROUTINE priphim
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)


OPEN(131,STATUS='unknown',FILE='phim.par')
WRITE(22,*) ndfull2

ind=0
DO iz=minz,maxz
  z1=(iz-nz)*dz
  DO iy=miny,maxy
    y1=(iy-ny)*dy
    DO ix=minx,maxx
      x1=(ix-nx)*dx
      
      ind=ind+1
      
      WRITE(131,'(3f10.2,3e19.9)') x1,y1,z1,phim(ind) ,phimv(ind), phimd(ind)
      
      
    END DO
  END DO
END DO

CLOSE(131)



END SUBROUTINE priphim
!------------------------------------------------------------
