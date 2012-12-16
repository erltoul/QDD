!_________________________________________Mat Orth__________________________________________________
 

REAL(8) FUNCTION matdorth(aa,n,ndim)

COMPLEX(8), INTENT(IN)                      :: aa(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: ndim



matdorth=0D0
DO i=1,ndim
  DO j=i,ndim
    IF(i==j) THEN
      matdorth=matdorth+SUM(aa(i,1:ndim)*CONJG(aa(j,1:ndim)))-1D0
    ELSE
      matdorth=matdorth+SUM(aa(i,1:ndim)*CONJG(aa(j,1:ndim)))
    END IF
  END DO
END DO
END FUNCTION matdorth
!_________________________________________Vec Print__________________________________________________
!                                         print VV

SUBROUTINE vecprint(str,vv,n,ndim)

CHARACTER (LEN=*), INTENT(IN OUT)        :: str
COMPLEX(8), INTENT(IN OUT)                  :: vv(n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN OUT)                  :: ndim



WRITE (6,*) str

WRITE (6, '(200(2g13.6,2x))') (vv(i),i=1,ndim)
END SUBROUTINE vecprint
!_________________________________________Mat Print__________________________________________________
!                                         print AA

SUBROUTINE matprint(str,aa,n,ndim)

CHARACTER (LEN=*), INTENT(IN OUT)        :: str
COMPLEX(8), INTENT(IN OUT)                  :: aa(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN OUT)                  :: ndim



WRITE (6,*) str

WRITE (6, '(200(2g13.6,2x))') ((aa(i,j),j=1,ndim),i=1,ndim)
END SUBROUTINE matprint

!_________________________________________Mat cnjg__________________________________________________
!                                         BBij=T(AAji)

SUBROUTINE matcnjg(aa,bb,n,ndim)

COMPLEX(8), INTENT(IN OUT)                  :: aa(n,n)
COMPLEX(8), INTENT(OUT)                     :: bb(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: ndim
COMPLEX(8) :: cc(n,n)


DO i=1,ndim
  DO j=1,ndim
    bb(i,j)=CONJG(aa(j,i))
  END DO
END DO
END SUBROUTINE matcnjg
!_________________________________________Mat mult__________________________________________________
!                                         CC=AA*BB

SUBROUTINE matmult(aa,bb,cc,n,ndim)

COMPLEX(8), INTENT(IN)                      :: aa(n,n)
COMPLEX(8), INTENT(IN)                      :: bb(n,n)
COMPLEX(8), INTENT(OUT)                     :: cc(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: ndim

INTEGER :: i,j,k

DO i=1,ndim
  DO j=1,ndim
    cc(i,j)=0
    DO k=1,ndim
      cc(i,j)=cc(i,j)+aa(i,k)*bb(k,j)
    END DO
  END DO
END DO
END SUBROUTINE matmult
!_________________________________________MatABtoA__________________________________________________
!                                         AA=AA*BB

SUBROUTINE matabtoa(aa,bb,n,ndim)

COMPLEX(8), INTENT(IN OUT)                  :: aa(n,n)
COMPLEX(8), INTENT(IN OUT)                  :: bb(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN OUT)                  :: ndim
COMPLEX(8) :: cc(n,n)


CALL matmult(aa,bb,cc,n,ndim)
CALL matcopy(cc,aa,n,ndim)
RETURN
END SUBROUTINE matabtoa

!_________________________________________Mat add__________________________________________________
!                                         CC=AA+BB

SUBROUTINE matadd(aa,bb,cc,n,ndim)

COMPLEX(8), INTENT(IN)                      :: aa(n,n)
COMPLEX(8), INTENT(IN)                      :: bb(n,n)
COMPLEX(8), INTENT(OUT)                     :: cc(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: ndim

INTEGER :: i,j

DO i=1,ndim
  DO j=1,ndim
    cc(i,j)=aa(i,j)+bb(i,j)
  END DO
END DO
END SUBROUTINE matadd
!_________________________________________Mat Sub__________________________________________________
!                                         CC=AA-BB

SUBROUTINE matsub(aa,bb,cc,n,ndim)

COMPLEX(8), INTENT(IN)                      :: aa(n,n)
COMPLEX(8), INTENT(IN)                      :: bb(n,n)
COMPLEX(8), INTENT(OUT)                     :: cc(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: ndim

INTEGER :: i,j

DO i=1,ndim
  DO j=1,ndim
    cc(i,j)=aa(i,j)-bb(i,j)
  END DO
END DO
END SUBROUTINE matsub

!_________________________________________Mat Copy__________________________________________________
!                                          BB=CC

SUBROUTINE matcopy(aa,bb,n,ndim)
!   Copy AA into BB

COMPLEX(8), INTENT(IN)                      :: aa(n,n)
COMPLEX(8), INTENT(OUT)                     :: bb(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: ndim

INTEGER :: i,j

DO i=1,ndim
  DO j=1,ndim
    bb(i,j)=aa(i,j)
  END DO
END DO
END SUBROUTINE matcopy

!_________________________________________Mat Cons__________________________________________________
!                                         BB=x*CC

SUBROUTINE matconst(aa,bb,x,n,ndim)

COMPLEX(8), INTENT(IN)                      :: aa(n,n)
COMPLEX(8), INTENT(OUT)                     :: bb(n,n)
COMPLEX(8), INTENT(IN)                      :: x
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: ndim

INTEGER :: i,j

DO i=1,ndim
  DO j=1,ndim
    bb(i,j)=aa(i,j)*x
  END DO
END DO
END SUBROUTINE matconst

!_________________________________________Mat A*Exp(B)_______________________________________________
!                                         AA=AA*Exp(BB)

SUBROUTINE mataexpb(aa,bb,n,ndim)

COMPLEX(8), INTENT(IN OUT)                  :: aa(n,n)
COMPLEX(8), INTENT(IN OUT)                  :: bb(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN OUT)                  :: ndim

COMPLEX(8) :: cc(n,n),dd(n,n)


CALL matexp(bb,cc,n,ndim)
CALL matmult(aa,cc,dd,n,ndim)
CALL matcopy(dd,aa,n,ndim)
RETURN
END SUBROUTINE mataexpb

!_________________________________________Mat Exp__________________________________________________
!                                         BB=Exp(AA)

SUBROUTINE matexp(aa,bb,n,ndim)

COMPLEX(8), INTENT(IN OUT)                  :: aa(n,n)
COMPLEX(8), INTENT(IN OUT)                  :: bb(n,n)
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN OUT)                  :: ndim

COMPLEX(8) :: cc(n,n), dd(n,n)

REAL(8) matnorme
REAL(8) eps, delta,rn

eps=1D-20
CALL matunite(bb,n,ndim)
CALL matunite(cc,n,ndim)
rn=CMPLX(1D0,0D0)
delta= 1D0
!        call matprint('AAds exp',AA,N,Ndim)
!        call matprint('BB ds exp',BB,N,Ndim)

DO WHILE (delta > eps)
  CALL matmult(aa,cc,dd,n,ndim)
!        call matprint('CC ds exp avant const',CC,N,Ndim)
  CALL matconst(dd,cc, 1D0/rn,n,ndim)
!        call matprint('CC ds exp',CC,N,Ndim)
  delta= matnorme(cc, n, ndim)
  CALL matadd(bb,cc,bb,n,ndim)
  rn=rn+1
END DO
END SUBROUTINE matexp
!_________________________________________Mat Unite__________________________________________________
!                                         AA=II

SUBROUTINE matunite(aa,n,ndim)
!     BB=unite

COMPLEX(8), INTENT(OUT)                     :: aa(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: ndim



DO i=1,ndim
  DO j=1,ndim
    aa(i,j)=CMPLX(0D0,0D0)
  END DO
  aa(i,i)=CMPLX(1D0,0D0)
END DO
END SUBROUTINE matunite
!_________________________________________Mat Norme__________________________________________________
!                                         norme=|aa|**2

REAL(8) FUNCTION matnorme(aa,n,ndim)

COMPLEX(8), INTENT(IN)                      :: aa(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: ndim

REAL(8) norme


norme=0
DO i=1,ndim
  DO j=1,ndim
    norme=norme + aa(i,j)*CONJG(aa(i,j))
  END DO
END DO
norme=SQRT(norme)
matnorme=norme
RETURN
END FUNCTION matnorme


!****************    real arithmetic
!_________________________________________rMat Orth__________________________________________________
 
 
!                                         sigma(Vi*cjg(Vj))

REAL(8) FUNCTION rmatdorth(aa,n,ndim)

REAL(8), INTENT(IN)                       :: aa(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: ndim



rmatdorth=0D0
DO i=1,ndim
  DO j=i,ndim
    IF(i==j) THEN
      rmatdorth=rmatdorth+SUM(aa(i,1:ndim)*aa(j,1:ndim))-1D0
    ELSE
      rmatdorth=rmatdorth+SUM(aa(i,1:ndim)*aa(j,1:ndim))
    END IF
  END DO
END DO
END FUNCTION rmatdorth
!_________________________________________rMat Print__________________________________________________
!                                         print AA

SUBROUTINE rmatprint(str,aa,n,ndim)

CHARACTER (LEN=*), INTENT(IN OUT)        :: str
REAL(8), INTENT(IN OUT)                   :: aa(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN OUT)                  :: ndim



WRITE (6,*) str

WRITE (6, '(4(g13.6))') ((aa(i,j),j=1,ndim),i=1,ndim)
END SUBROUTINE rmatprint

!_________________________________________rMat cnjg__________________________________________________
!                                         BBij=T(AAji)

SUBROUTINE rmatcnjg(aa,bb,n,ndim)

REAL(8), INTENT(IN)                       :: aa(n,n)
REAL(8), INTENT(OUT)                      :: bb(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: ndim
REAL(8)  cc(n,n)


DO i=1,ndim
  DO j=1,ndim
    bb(i,j)=aa(j,i)
  END DO
END DO
END SUBROUTINE rmatcnjg
!_________________________________________rMat mult__________________________________________________
!                                         CC=AA*BB

SUBROUTINE rmatmult(aa,bb,cc,n,ndim)

REAL(8), INTENT(IN)                       :: aa(n,n)
REAL(8), INTENT(IN)                       :: bb(n,n)
REAL(8), INTENT(OUT)                      :: cc(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: ndim

INTEGER :: i,j,k

DO i=1,ndim
  DO j=1,ndim
    cc(i,j)=0
    DO k=1,ndim
      cc(i,j)=cc(i,j)+aa(i,k)*bb(k,j)
    END DO
  END DO
END DO
END SUBROUTINE rmatmult
!_________________________________________rMatABtoA__________________________________________________
!                                         AA=AA*BB

SUBROUTINE rmatabtoa(aa,bb,n,ndim)

REAL(8), INTENT(IN OUT)                   :: aa(n,n)
REAL(8), INTENT(IN OUT)                   :: bb(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN OUT)                  :: ndim
REAL(8)  cc(n,n)


CALL rmatmult(aa,bb,cc,n,ndim)
CALL rmatcopy(cc,aa,n,ndim)
RETURN
END SUBROUTINE rmatabtoa

!_________________________________________rMat add__________________________________________________
!                                         CC=AA+BB

SUBROUTINE rmatadd(aa,bb,cc,n,ndim)

REAL(8), INTENT(IN)                       :: aa(n,n)
REAL(8), INTENT(IN)                       :: bb(n,n)
REAL(8), INTENT(OUT)                      :: cc(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: ndim

INTEGER :: i,j

DO i=1,ndim
  DO j=1,ndim
    cc(i,j)=aa(i,j)+bb(i,j)
  END DO
END DO
END SUBROUTINE rmatadd
!_________________________________________rMat Sub__________________________________________________
!                                         CC=AA-BB

SUBROUTINE rmatsub(aa,bb,cc,n,ndim)

REAL(8), INTENT(IN)                       :: aa(n,n)
REAL(8), INTENT(IN)                       :: bb(n,n)
REAL(8), INTENT(OUT)                      :: cc(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: ndim

INTEGER :: i,j

DO i=1,ndim
  DO j=1,ndim
    cc(i,j)=aa(i,j)-bb(i,j)
  END DO
END DO
END SUBROUTINE rmatsub

!_________________________________________rMat Copy__________________________________________________
!                                          BB=CC

SUBROUTINE rmatcopy(aa,bb,n,ndim)
!   Copy AA into BB

REAL(8), INTENT(IN)                       :: aa(n,n)
REAL(8), INTENT(OUT)                      :: bb(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: ndim

INTEGER :: i,j

DO i=1,ndim
  DO j=1,ndim
    bb(i,j)=aa(i,j)
  END DO
END DO
END SUBROUTINE rmatcopy

!_________________________________________rMat Cons__________________________________________________
!                                         BB=x*CC

SUBROUTINE rmatconst(aa,bb,x,n,ndim)

REAL(8), INTENT(IN)                       :: aa(n,n)
REAL(8), INTENT(OUT)                      :: bb(n,n)
REAL(8), INTENT(IN)                       :: x
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: ndim

INTEGER :: i,j

DO i=1,ndim
  DO j=1,ndim
    bb(i,j)=aa(i,j)*x
  END DO
END DO
END SUBROUTINE rmatconst

!_________________________________________rMat A*Exp(B)_______________________________________________
!                                         AA=AA*Exp(BB)

SUBROUTINE rmataexpb(aa,bb,n,ndim)

REAL(8), INTENT(IN OUT)                   :: aa(n,n)
REAL(8), INTENT(IN OUT)                   :: bb(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN OUT)                  :: ndim

REAL(8) cc(n,n),dd(n,n)


CALL rmatexp(bb,cc,n,ndim)
CALL rmatmult(aa,cc,dd,n,ndim)
CALL rmatcopy(dd,aa,n,ndim)
RETURN
END SUBROUTINE rmataexpb

!_________________________________________rMat Exp__________________________________________________
!                                         BB=Exp(AA)

SUBROUTINE rmatexp(aa,bb,n,ndim)

REAL(8), INTENT(IN OUT)                   :: aa(n,n)
REAL(8), INTENT(IN OUT)                   :: bb(n,n)
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN OUT)                  :: ndim

REAL(8) cc(n,n), dd(n,n)

REAL(8) rmatnorme
REAL(8) eps, delta,rn

eps=1.e-20
CALL rmatunite(bb,n,ndim)
CALL rmatunite(cc,n,ndim)
rn=CMPLX(1D0,0D0)
delta= 1D0
!        call rMatprint('AAds exp',AA,N,Ndim)
!        call rMatprint('BB ds exp',BB,N,Ndim)

DO WHILE (delta > eps)
  CALL rmatmult(aa,cc,dd,n,ndim)
!        call rMatprint('CC ds exp avant const',CC,N,Ndim)
  CALL rmatconst(dd,cc, 1D0/rn,n,ndim)
!        call rMatprint('CC ds exp',CC,N,Ndim)
  delta= rmatnorme(cc, n, ndim)
  CALL rmatadd(bb,cc,bb,n,ndim)
  rn=rn+1
END DO
END SUBROUTINE rmatexp
!_________________________________________rMat Unite__________________________________________________
!                                         AA=II

SUBROUTINE rmatunite(aa,n,ndim)
!     BB=unite

REAL(8), INTENT(OUT)                      :: aa(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: ndim



DO i=1,ndim
  DO j=1,ndim
    aa(i,j)=0D0
  END DO
  aa(i,i)=1D0
END DO
END SUBROUTINE rmatunite
!_________________________________________rMat Norme__________________________________________________
!                                         norme=|aa|**2

REAL(8) FUNCTION rmatnorme(aa,n,ndim)

REAL(8), INTENT(IN)                       :: aa(n,n)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: ndim

REAL(8) norme


norme=0
DO i=1,ndim
  DO j=1,ndim
    norme=norme + aa(i,j)*aa(i,j)
  END DO
END DO
norme=SQRT(norme)
rmatnorme=norme
RETURN
END FUNCTION rmatnorme

