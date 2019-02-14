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

MODULE coulsolv
USE coulsolv_f
USE coulsolv_e 
LOGICAL,PUBLIC :: tcoulfalr=.FALSE.
!INTEGER,PARAMETER :: DP=KIND(1D0)  ! precision  setting

!INTERFACE
!  SUBROUTINE solv_poisson_i(rhoinp,chpfalr,kdum)
!  REAL(8), INTENT(IN)   :: rhoinp(*)
!  REAL(8), INTENT(OUT)  :: chpfalr(*)
!  INTEGER, INTENT(IN)    :: kdum   ! dummy variable
!  END SUBROUTINE solv_poisson_i
!END INTERFACE
!PROCEDURE(solv_poisson_i),POINTER :: solv_poisson => NULL()
!EXTERNAL :: solv_poisson_f,solv_poisson_e

CONTAINS

SUBROUTINE init_coul(dxi,dyi,dzi,nxi,nyi,nzi)

REAL(DP) :: dxi,dyi,dzi
INTEGER ::nxi,nyi,nzi

IF(tcoulfalr) THEN
  CALL init_coul_f(dxi,dyi,dzi,nxi,nyi,nzi)
!  solv_poisson => solv_poisson_f
ELSE
  CALL init_coul_e(dxi,dyi,dzi,nxi,nyi,nzi)
!  solv_poisson => solv_poisson_e
END IF

END SUBROUTINE init_coul


SUBROUTINE solv_poisson(rhoinp,chpfalr,kdum)
  REAL(8), INTENT(IN)   :: rhoinp(*)
  REAL(8), INTENT(OUT)  :: chpfalr(*)
  INTEGER, INTENT(IN)    :: kdum   ! dummy variable

  IF(tcoulfalr) THEN
    CALL solv_poisson_f(rhoinp,chpfalr,kdum)
  ELSE
    CALL solv_poisson_e(rhoinp,chpfalr,kdum)
  END IF

END SUBROUTINE solv_poisson

END MODULE coulsolv


#if(netlib_fft|fftw_cpu)

!#if(coufou)
#include "falr.F90"
!#endif

!#if(coudoub)
#include "coulex.F90"
!#endif

#endif

#if(findiff|numerov)
!INCLUDE "gridcoul.F90"
#include "solv_poisson.F90"
#endif
