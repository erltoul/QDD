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

MODULE orthmat
USE params
IMPLICIT NONE


INTERFACE matdorth
  MODULE PROCEDURE matdorth_r, matdorth_c
END INTERFACE matdorth

INTERFACE matunite
  MODULE PROCEDURE matunite_r, matunite_c
END INTERFACE matunite

INTERFACE matexp
  MODULE PROCEDURE matexp_r, matexp_c
END INTERFACE matexp

INTERFACE orthnorm
  MODULE PROCEDURE orthnorm_r, orthnorm_c
END INTERFACE orthnorm

INTERFACE vecovlp
  MODULE PROCEDURE vecovlp_r, vecovlp_c
END INTERFACE vecovlp

INTERFACE vecnorm
  MODULE PROCEDURE vecnorm_r, vecnorm_c
END INTERFACE vecnorm

CONTAINS

!_________________________________________Mat Orth__________________________________________________
!                                     sigma(Vi*cjg(Vj))
!-------------------------------------------------------
! REAL version
!-------------------------------------------------------
REAL(DP) FUNCTION matdorth_r(aa,n,ndim)

REAL(DP), INTENT(IN)     :: aa(n,n)
INTEGER, INTENT(IN OUT) :: n
INTEGER, INTENT(IN)     :: ndim
INTEGER                 :: i,j
matdorth_r=0D0
DO i=1,ndim
  DO j=i,ndim
    IF(i==j) THEN
      matdorth_r=matdorth_r+SUM(aa(i,1:ndim)*aa(j,1:ndim))-1D0
    ELSE
      matdorth_r=matdorth_r+SUM(aa(i,1:ndim)*aa(j,1:ndim))
    END IF
  END DO
END DO
END FUNCTION matdorth_r

!-------------------------------------------------------
! COMPLEX version
!-------------------------------------------------------
REAL(DP) FUNCTION matdorth_c(aa,n,ndim)

COMPLEX(DP), INTENT(IN)  :: aa(n,n)
INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: ndim
INTEGER                 :: i,j
matdorth_c=0D0
DO i=1,ndim
  DO j=i,ndim
    IF(i==j) THEN
      matdorth_c=matdorth_c+SUM(aa(i,1:ndim)*CONJG(aa(j,1:ndim)))-1D0
    ELSE
      matdorth_c=matdorth_c+SUM(aa(i,1:ndim)*CONJG(aa(j,1:ndim)))
    END IF
  END DO
END DO
END FUNCTION matdorth_c
!-------------------------------------------------------
! cmplxsic version
!-------------------------------------------------------

REAL(DP) FUNCTION matdorth_cmplxsic(aa,n,ndim)

COMPLEX(DP), INTENT(IN)  :: aa(n,n)
INTEGER, INTENT(IN OUT) :: n
INTEGER, INTENT(IN)     :: ndim
INTEGER                 :: i,j
matdorth_cmplxsic=0D0
DO i=1,ndim
  DO j=i,ndim
    IF(i==j) THEN
      matdorth_cmplxsic=matdorth_cmplxsic+SUM(aa(i,1:ndim)*aa(j,1:ndim))-1D0
    ELSE
      matdorth_cmplxsic=matdorth_cmplxsic+SUM(aa(i,1:ndim)*aa(j,1:ndim))
    END IF
  END DO
END DO
END FUNCTION matdorth_cmplxsic

!_________________________________________rMat Unite__________________________________________________
!                                         AA=II
!-------------------------------------------------------
! REAL version
!-------------------------------------------------------
SUBROUTINE matunite_r(aa,n,ndim)

REAL(DP), INTENT(OUT)  :: aa(n,n)
INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: ndim
INTEGER               :: i
aa(1:ndim,1:ndim)=0D0
DO i=1,ndim
  aa(i,i)=1D0
END DO

END SUBROUTINE matunite_r

!-------------------------------------------------------
! COMPLEX version
!-------------------------------------------------------
SUBROUTINE matunite_c(aa,n,ndim)

COMPLEX(DP), INTENT(OUT) :: aa(n,n)
INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: ndim
INTEGER                 :: i
aa(1:ndim, 1:ndim)=CMPLX(0D0,0D0,DP)
DO i=1,ndim
  aa(i,i)=CMPLX(1D0,0D0,DP)
END DO

END SUBROUTINE matunite_c

!_________________________________________rMat Exp__________________________________________________
!                                         BB=Exp(AA)
!-------------------------------------------------------
! REAL version
!-------------------------------------------------------
SUBROUTINE matexp_r(aa,bb,n,ndim)

REAL(DP), INTENT(IN OUT) :: aa(n,n)
REAL(DP), INTENT(IN OUT) :: bb(n,n)
INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN OUT) :: ndim
REAL(DP) cc(n,n), dd(n,n)

REAL(DP) eps, delta,rn

CALL matunite(bb,n,ndim)
CALL matunite(cc,n,ndim)
eps=1D-20
rn=1D0
delta=1D0
DO WHILE (delta > eps)
  dd(1:ndim,1:ndim) = MATMUL(aa(1:ndim,1:ndim),cc(1:ndim,1:ndim))
  cc(1:ndim,1:ndim) = dd(1:ndim,1:ndim)/rn
  delta = SQRT(SUM(cc(1:ndim,1:ndim)**2))
  bb(1:ndim,1:ndim) = cc(1:ndim,1:ndim) + bb(1:ndim,1:ndim) 
  rn=rn+1D0
END DO

END SUBROUTINE matexp_r

!-------------------------------------------------------
! REAL version
!-------------------------------------------------------
SUBROUTINE matexp_c(aa,bb,n,ndim)

COMPLEX(DP), INTENT(IN OUT)  :: aa(n,n)
COMPLEX(DP), INTENT(IN OUT)  :: bb(n,n)
INTEGER, INTENT(IN)         :: n
INTEGER, INTENT(IN OUT)     :: ndim

COMPLEX(DP) :: cc(n,n), dd(n,n)

!REAL(DP) matnorme
REAL(DP) eps, delta,rn

eps=1D-20
CALL matunite(bb,n,ndim)
CALL matunite(cc,n,ndim)
rn=CMPLX(1D0,0D0,DP)
delta= 1D0
!        call matprint('AAds exp',AA,N,Ndim)
!        call matprint('BB ds exp',BB,N,Ndim)

DO WHILE (delta > eps)
!  CALL matmult(aa,cc,dd,n,ndim)
!        call matprint('CC ds exp avant const',CC,N,Ndim)
!  CALL matconst(dd,cc, 1D0/rn,n,ndim)
  dd(1:ndim,1:ndim) = MATMUL(aa(1:ndim,1:ndim),cc(1:ndim,1:ndim))
  cc(1:ndim,1:ndim) =  dd(1:ndim,1:ndim)/rn
!        call matprint('CC ds exp',CC,N,Ndim)
!  delta= matnorme(cc, n, ndim)
  delta = SQRT(SUM(cc(1:ndim,1:ndim)*CONJG(cc(1:ndim,1:ndim))))
!  CALL matadd(bb,cc,bb,n,ndim)
  bb(1:ndim,1:ndim) =  cc(1:ndim,1:ndim) + bb(1:ndim,1:ndim)
  rn=rn+1
END DO
END SUBROUTINE matexp_c

!_________________________________________OrthNorm__________________________________________________
!     ortho-normalizes system of vectors 'vecs' with
!     dimension 'ndim'.
!-------------------------------------------------------
! REAL version
!-------------------------------------------------------
SUBROUTINE orthnorm_r(vecs,ndim,kdim)
USE params, ONLY: DP

INTEGER,INTENT(IN) :: ndim,kdim
INTEGER :: i,j
REAL(DP) :: acc2
REAL(DP) :: acc
REAL(DP),INTENT(IN OUT) :: vecs(kdim,kdim)

!-------------------------------------------------------
acc2 = 1D0/SQRT(vecnorm(vecs(:,1),ndim))
vecs(1:ndim,1) = vecs(1:ndim,1)*acc2
DO i=2,ndim
  DO j=1,i-1
    acc = vecovlp(vecs(:,j),vecs(:,i),ndim)
    vecs(1:ndim,i) = vecs(1:ndim,i)-acc*vecs(1:ndim,j)
  END DO
  acc2 = 1D0/SQRT(vecnorm(vecs(:,i),ndim))
  vecs(1:ndim,i) = vecs(1:ndim,i)*acc2
END DO

RETURN
END SUBROUTINE orthnorm_r
!-------------------------------------------------------
! COMPLEX version
!-------------------------------------------------------
SUBROUTINE orthnorm_c(vecs,ndim,kdim)
USE params, ONLY: DP

INTEGER,INTENT(IN) :: ndim,kdim
INTEGER :: i,j
REAL(DP) :: acc2
COMPLEX(DP) :: acc
COMPLEX(DP),INTENT(IN OUT) :: vecs(kdim,kdim)

!-------------------------------------------------------
acc2 = 1D0/SQRT(vecnorm(vecs(:,1),ndim))
vecs(1:ndim,1) = vecs(1:ndim,1)*acc2
DO i=2,ndim
  DO j=1,i-1
    acc = vecovlp(vecs(:,j),vecs(:,i),ndim)
    vecs(1:ndim,i) = vecs(1:ndim,i)-acc*vecs(1:ndim,j)
  END DO
  acc2 = 1D0/SQRT(vecnorm(vecs(:,i),ndim))
  vecs(1:ndim,i) = vecs(1:ndim,i)*acc2
END DO

RETURN
END SUBROUTINE orthnorm_c


!_________________________________________vecovlp__________________________________________________
!     Overlap 'vec1' with 'vec2' having dimension 'ndim'.
!-------------------------------------------------------
! REAL version
!-------------------------------------------------------
REAL(DP) FUNCTION vecovlp_r(vec1,vec2,ndim)

USE params, ONLY: DP
INTEGER,INTENT(IN) :: ndim
REAL(DP),INTENT(IN) :: vec1(ndim),vec2(ndim)

vecovlp_r = SUM(vec1(1:ndim)*vec2(1:ndim))

RETURN
END FUNCTION vecovlp_r
!-------------------------------------------------------
!COMPLEX version
!-------------------------------------------------------
COMPLEX(DP) FUNCTION vecovlp_c(vec1,vec2,ndim)

USE params, ONLY: DP
INTEGER,INTENT(IN) :: ndim
COMPLEX(DP),INTENT(IN) :: vec1(ndim),vec2(ndim)

vecovlp_c= SUM(CONJG( vec1(1:ndim))*vec2(1:ndim))

RETURN
END FUNCTION vecovlp_c

!_________________________________________vecnorm__________________________________________________
!     Calculates norm of vector 'vec' with
!     dimension 'ndim'.
!-------------------------------------------------------
! REAL version
!-------------------------------------------------------
REAL(DP) FUNCTION vecnorm_r(vec,ndim)
USE params, ONLY: DP
INTEGER,INTENT(IN) :: ndim
REAL(DP),INTENT(IN) :: vec(ndim)

vecnorm_r = SUM(vec(1:ndim)**2)

RETURN
END FUNCTION vecnorm_r
!-------------------------------------------------------
! COMPLEX version
!-------------------------------------------------------
REAL(DP) FUNCTION vecnorm_c(vec,ndim)

USE params, ONLY: DP
INTEGER,INTENT(IN) :: ndim
COMPLEX(DP),INTENT(IN) :: vec(ndim)

vecnorm_c=SUM(ABS(vec(1:ndim)**2))

RETURN
END FUNCTION vecnorm_c

END MODULE orthmat
