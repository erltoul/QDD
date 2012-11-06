#include "define.h"
#if(twostsic)

#ifdef REALSWITCH
MODULE twostr
USE params
USE kinetic
USE localize_rad
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP),PRIVATE :: qnewr(kdfull2,kstate)
REAL(DP),PRIVATE :: psirut(kdfull2,kstate)
!COMMON /twost/ qnewr,qnew,psirut,psiut

REAL(DP) :: step,precis,radmaxsym
INTEGER :: symutbegin,itut
!COMMON /twostut/ step,precis,symutbegin,radmaxsym,itut

INTEGER,PARAMETER,PRIVATE :: kdim=kstate
!     matrices of radial moments

REAL(DP),PRIVATE :: rrmatr(kdim,kdim,2)  ! matrix of r**2
REAL(DP),PRIVATE :: xxmatr(kdim,kdim,2)  ! matrix of x**2
REAL(DP),PRIVATE :: yymatr(kdim,kdim,2)  ! matrix of y**2
REAL(DP),PRIVATE :: zzmatr(kdim,kdim,2)  ! matrix of z**2
REAL(DP),PRIVATE :: xmatr(kdim,kdim,2)   ! matrix of x
REAL(DP),PRIVATE :: ymatr(kdim,kdim,2)   ! matrix of y
REAL(DP),PRIVATE :: zmatr(kdim,kdim,2)   ! matrix of z
REAL(DP) :: vecsr(kdim,kdim,2)    ! searched eigenvectors
INTEGER :: ndims(2)
!COMMON /radmatrix/ rrmatr,xxmatr,yymatr,zzmatr, xmatr,ymatr,zmatr,  &
!    vecsr
!COMMON /radmatrixn/ndims

#endif

#ifdef COMPLEXSWITCH
MODULE twost
USE params
USE kinetic
USE twostr, ONLY:step,precis,symutbegin,radmaxsym,itut,vecsr
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP) :: qnewut(kdfull2,kstate)
COMPLEX(DP) :: psiut(kdfull2,kstate)

INTEGER,PARAMETER,PRIVATE :: kdim=kstate

!     matrices of radial moments

COMPLEX(DP),PRIVATE :: rrmatr(kdim,kdim,2)  ! matrix of r**2
COMPLEX(DP),PRIVATE :: xxmatr(kdim,kdim,2)  ! matrix of x**2
COMPLEX(DP),PRIVATE :: yymatr(kdim,kdim,2)  ! matrix of y**2
COMPLEX(DP),PRIVATE :: zzmatr(kdim,kdim,2)  ! matrix of z**2
COMPLEX(DP),PRIVATE :: xmatr(kdim,kdim,2)   ! matrix of x
COMPLEX(DP),PRIVATE :: ymatr(kdim,kdim,2)   ! matrix of y
COMPLEX(DP),PRIVATE :: zmatr(kdim,kdim,2)   ! matrix of z
COMPLEX(DP) :: vecs(kdim,kdim,2)    ! searched eigenvectors
INTEGER,PRIVATE :: ndims(2)
!COMMON /radmatrix/ rrmatr,xxmatr,yymatr,zzmatr, xmatr,ymatr,zmatr,  &
!    vecs
!COMMON /radmatrixn/ndims

#endif

CONTAINS

#ifdef REALSWITCH
!-----init_fsic------------------------------------------------

SUBROUTINE init_fsic()

!     initializes fields for FSIC etc

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

!INCLUDE "twost.inc"
!INCLUDE 'radmatrixr.inc'
NAMELIST /fsic/step,precis,symutbegin,radmaxsym   !!! UT parameters

!----------------------------------------------------------------

IF(ifsicp < 6) RETURN

!     read SIC  specific parameters

OPEN(5,STATUS='old',FORM='formatted',FILE='for005.'//outnam)
READ(5,fsic)
WRITE(6,'(a,4(1pg13.5))')  &
    ' SIC running with step,precis,SymUtBegin,radmaxsym=',  &
    step,precis,symutbegin,radmaxsym
CLOSE(5)

!     dimensions of spin sup-spaces

ndims(1) = 0
ndims(2) = 0
DO na=1,nstate
  IF(ispin(nrel2abs(na)) == 1) THEN
    ndims(1) = na
  ELSE IF(ispin(nrel2abs(na)) == 2) THEN
    ndims(2) = na-ndims(1)
  END IF
END DO
WRITE(6,'(a,2i5)') ' dimension of sub-matrices:',ndims

!     initialize unitary matrix

DO is=1,2
  DO i=1,ndims(is)
    DO j=1,ndims(is)
      IF(i == j) THEN
        vecsr(i,j,is) = 1.D0
      ELSE
        vecsr(i,j,is) = 0.D0
      END IF
    END DO
  END DO
END DO

RETURN
END SUBROUTINE init_fsic
#endif




!     ******************************

#ifdef REALSWITCH
SUBROUTINE calc_fullsicr(q0,qsic)
#else
SUBROUTINE calc_fullsic(q0,qsic)
#endif

!     ******************************

!     full SIC:
!       input is set of wavefunctions on 'q0'
!       output are SIC s.p. wavefunctions on 'qsic'

USE params
USE kinetic
#if(symmcond)
USE symcond
!COMMON /sicsav/usicall(kdfull2,kstate)
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

#ifdef REALSWITCH
REAL(DP) :: q0(kdfull2,kstate)
REAL(DP) :: qsic(kdfull2,kstate)
#else
COMPLEX(DP) :: q0(kdfull2,kstate)
COMPLEX(DP) :: qsic(kdfull2,kstate)
#endif

!       workspaces

REAL(DP),ALLOCATABLE :: usicsp(:),rhosp(:)

LOGICAL,PARAMETER :: ttest=.false.

IF(numspin.NE.2) STOP "CALC_FULLSIC requires full spin"

!mb
enrearsave=enrear
enerpwsave=enerpw
enrear1=0D0
enrear2=0D0
enpw1  = 0D0
enpw2  = 0D0
encadd = 0D0
!mb/



!------------------------------------------------------------------

!     compute action of SIC potential
ALLOCATE(usicsp(2*kdfull2),rhosp(2*kdfull2))

DO nb=1,nstate
  IF(occup(nb) > small) THEN
    ishift = (ispin(nrel2abs(nb))-1)*nxyz ! store spin=2 in upper block
#ifdef REALSWITCH
    CALL calc_sicspr(rhosp,usicsp,q0(1,nb),nb)
#else
    CALL calc_sicsp(rhosp,usicsp,q0(1,nb),nb)
#endif
!#if(twostsic)
#if(symmcond)
    DO ind=1,nxyz
      usicall(ind,nb) = usicsp(ind+ishift)
    END DO
#endif
    
    IF (ispin(nrel2abs(nb)) == 1) THEN
      enrear1=enrear1+enrear*occup(nb)
      enpw1=enpw1+enerpw*occup(nb)
    ELSE
      enrear2=enrear2+enrear*occup(nb)
      IF(directenergy) THEN
        enpw2=enpw2+enerpw*occup(nb)
      END IF
    END IF
    IF(directenergy) THEN
      encadd=encadd+encoulsp*occup(nb)
    END IF
    
    IF (ispin(nrel2abs(nb)) == 1) THEN
      DO ind=1,nxyz
        qsic(ind,nb) = usicsp(ind)*q0(ind,nb)
      END DO
    ELSE
      DO ind=1,nxyz
        idx = ind+nxyz
        qsic(ind,nb) = usicsp(idx)*q0(ind,nb)
      END DO
    END IF
  END IF
END DO
encadd=encadd/2.0
enrear   = enrearsave-enrear1-enrear2
IF(directenergy) THEN
  enerpw   = enerpwsave-enpw1-enpw2-encadd
END IF
DEALLOCATE(usicsp,rhosp)

RETURN
#ifdef REALSWITCH
END SUBROUTINE calc_fullsicr
#else
END SUBROUTINE calc_fullsic
#endif


#ifdef COMPLEXSWITCH
!-----init_vecs-------------------------------------------------

SUBROUTINE init_vecs(vecs)

!     Initializes complex superposition coefficients for SIC.

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
!INCLUDE 'twost.inc'
!INCLUDE 'radmatrixr.inc'
COMPLEX(DP) :: vecs(kdim,kdim,2)    ! searched eigenvectors

!----------------------------------------------------------------

DO i=1,kstate
  DO j=1,kstate
    vecs(i,j,1) = CMPLX(vecsr(i,j,1),0D0,DP)
    vecs(i,j,2) = CMPLX(vecsr(i,j,2),0D0,DP)
  END DO
END DO

RETURN
END SUBROUTINE init_vecs
#endif

#ifdef REALSWITCH

!-----static_sicfield------------------------------------------------

SUBROUTINE static_sicfield(rho,aloc,psir,iter1)

!     The Coulomb part of the mean field.

!     Input:
!      rho    = electron density
!      psir   = real wavefunctions
!     Output:
!      aloc   = local mean-field potential

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
REAL(DP) :: rho(2*kdfull2)
REAL(DP) :: aloc(2*kdfull2)
REAL(DP) :: psir(kdfull2,kstate)
!INCLUDE "twost.inc"
!INCLUDE 'radmatrixr.inc'

!----------------------------------------------------------------

IF(ifsicp < 7) RETURN

!     computation of the mean fields

CALL calc_utwfr(psir,psirut,iter1)
!        write(*,*) ' UTWFR over. usew1=',usew1
IF(ifsicp == 7) THEN       ! Generalized Slater pot
  ifsicp=3
  CALL calc_sicr(rho,aloc,psirut)
  ifsicp=7
!          write(*,*) ' CALC_SICR over'
ELSE IF(ifsicp == 8) THEN   !    DSIC
  CALL calc_fullsicr(psirut,qnewr)
END IF

RETURN
END SUBROUTINE static_sicfield


!-----infor_sic------------------------------------------------

SUBROUTINE infor_sic(psir)

!     Compute and print variances


USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
REAL(DP) :: rho(2*kdfull2)
REAL(DP) :: aloc(2*kdfull2)
REAL(DP) :: psir(kdfull2,kstate)
!INCLUDE "twost.inc"
!INCLUDE 'radmatrixr.inc'

!----------------------------------------------------------------

IF(ifsicp < 7) RETURN

!     TESTS

IF(ifsicp > 6)THEN
  WRITE(6,*) 'DIAGONAL STATES :'
  CALL spmomsmatrixr(psir,1)       !!! to print the total variance
  WRITE(6,*) 'LOCALIZED STATES :'
  CALL spmomsmatrixr(psirut,1)
END IF

IF(ifsicp == 8) THEN   !!! to calculate the total
  DO is=1,2                       !!! violation of symmetry condition
    sum = 0.D0
    DO na=1,nstate
      IF(ispin(na) == is)THEN
        DO nb=1,nstate
          IF(ispin(nb) == is)THEN
            sum = ( rwfovlp(psirut(1,na),qnewr(1,nb)) -  &
                rwfovlp(psirut(1,nb),qnewr(1,na)) )**2 + sum
          END IF
        END DO
      END IF
    END DO
    WRITE(6,'(a,i3,a,(1pg12.4))')  &
        'For spin',is,'  Total violation of SymCond',SQRT(sum)
  END DO
END IF

RETURN
END SUBROUTINE infor_sic


#if(twostsic)
!-----diag_lagr------------------------------------------------

SUBROUTINE diag_lagr(psir)

!     Diagonalize the matrix of Lagrange parameters and print
!     information from that step.

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

!INCLUDE "twost.inc"
!INCLUDE 'radmatrixr.inc'
REAL(DP) :: psir(kdfull2,kstate)

!     vect1,2 must be declared in each spin subspace for correct
!     diagonalization of SIC with multipliers

REAL(DP) :: a(nstate*(nstate+1)/2),root(kdim,2),  &
            vect1(nstate-nspdw,nstate-nspdw),vect2(nspdw,nspdw)

!----------------------------------------------------------------

IF(ifsicp < 7) RETURN

DO is=1,2             !!! Diagonalization on the Lagrange matrix for FSIC
  
  i=1
  DO na=1,nstate
    IF(ispin(nrel2abs(na)) == is) THEN
      DO nb=1,na
        IF(ispin(nrel2abs(nb)) == is) THEN
          a(i) = hmatrix(nb,na)
          i=i+1
        END IF
      END DO
    END IF
  END DO
  
! Lagrange mult matrix diag (needs vect1,2 in each sp subs, instead bugs)
  IF(is == 1)THEN
    CALL givens(a,root(1,is),vect1(1,1),ndims(is),ndims(is),ndims(is))
  ELSE IF(is == 2)THEN
    CALL givens(a,root(1,is),vect2(1,1),ndims(is),ndims(is),ndims(is))
  END IF
  
  DO na=1,ndims(is)
    DO nb=1,ndims(is)
      IF(is == 1)THEN
        vecsr(na,nb,is) = vect1(nb,na)
      ELSE IF(is == 2)THEN
        vecsr(na,nb,is) = vect2(nb,na)
      END IF
    END DO
  END DO
  
  DO nb=1,nstate
    IF(is == ispin(nrel2abs(nb)))THEN
      nbeff = nb - (is-1)*ndims(1)
      IF(is == 1)THEN
        CALL superpose_stater(psir(1,nb),vect1(1,nbeff),psirut,is)
      ELSE IF(is == 2)THEN
        CALL superpose_stater(psir(1,nb),vect2(1,nbeff),psirut,is)
      END IF
    END IF
  END DO
  
  WRITE(6,*) 'Diag of the Lagrange matrix for spin subs',is
  DO j=1,ndims(is)
    WRITE(6,*) 'state',j,'single energies',root(j,is)
  END DO
  
END DO

ifsicp=8

WRITE(6,*) 'DIAGONAL STATES :'
CALL spmomsmatrixr(psir,1)  !!! to print the total variance
WRITE(6,*) 'LOCALIZED STATES :'
CALL spmomsmatrixr(psirut,1)

RETURN
END SUBROUTINE diag_lagr
#endif

!-----subtr_sicpot------------------------------------------------

SUBROUTINE subtr_sicpot(q1,nbe)

!     Diagonalize the matrix of Lagrange parameters and print
!     information from that step.

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
!INCLUDE "twost.inc"
!INCLUDE 'radmatrixr.inc'
REAL(DP) :: q1(kdfull2)

!----------------------------------------------------------------

IF(ifsicp /= 8) RETURN

!JM : subtract SIC potential for state NBE

is=ispin(nbe)
DO na=1,ndims(is)
  nb = nbe - (is-1)*ndims(1)
  nae = na + (is-1)*ndims(1)
  cf = vecsr(nb,na,is)
  DO i=1,nxyz
    q1(i)=q1(i)-qnewr(i,nae)*cf
  END DO
END DO

RETURN
END SUBROUTINE subtr_sicpot

#endif
!-----calc_utwf--------------------------------------------------

#ifdef REALSWITCH
SUBROUTINE calc_utwfr(q0,q0ut,iter1)
#else
SUBROUTINE calc_utwf(q0,q0ut,iter1)
#endif


!     computes localized wavefunctions
!       input is set of wavefunctions on 'q0'
!       output are localized wavefunctions 'q0UT'
!       the array 'qnewUT' is used as workspace

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

!     basic arrays and workspace

!INCLUDE 'twost.inc'
#ifdef REALSWITCH
!INCLUDE 'radmatrixr.inc'
REAL(DP) :: q0(kdfull2,kstate)
REAL(DP) :: q0ut(kdfull2,kstate)
#else
!INCLUDE 'radmatrix.inc'
COMPLEX(DP) :: q0(kdfull2,kstate)
COMPLEX(DP) :: q0ut(kdfull2,kstate)
COMPLEX(DP) :: wfovlp
#endif

!------------------------------------------------------------------

!     Compute matrix of radial moments and determine
!     optimally localizing transformation.
!     Results is unitary matrix 'vecs' communicated via
!     common /radmatrix/.

DO is=1,2
  IF(ndims(is) > 1)THEN
#ifdef REALSWITCH
    CALL utgradstepr(is,0,q0,iter1)   !!!!! gives the new vecs (gradient method)
#else
!      do nb=1,nstate
!        write(*,*) ' CALC_UTWF: nb,Q0 norm:',
!     &       nb,wfovlp(q0(1,nb),q0(1,nb))
!      enddo
    CALL utgradstep(is,0,q0,iter1)   !!!!! gives the new vecs (gradient method)
#endif
  END IF
END DO

DO nb=1,nstate
  is = ispin(nrel2abs(nb))
  nbeff = nb - (is-1)*ndims(1)
#ifdef REALSWITCH
  CALL superpose_stater(q0ut(1,nb),vecsr(1,nbeff,is),q0,is)
#else
  CALL superpose_state(q0ut(1,nb),vecs(1,nbeff,is),q0,is)
#endif
END DO

RETURN
#ifdef REALSWITCH
END SUBROUTINE calc_utwfr
#else
END SUBROUTINE calc_utwf
#endif



!-----utgradstep-----------------------------------------

#ifdef REALSWITCH
SUBROUTINE utgradstepr(is,iprint,q0,iter1)
#else
SUBROUTINE utgradstep(is,iprint,q0,iter1)
#endif

!c     Nonlinear gradient iteration to optimally localized states:
!c      'vecs'    system of eigen-vectors to be determined
!c      'is'      isospin
!c      'iprint'  print level: <0 --> no print at all
!c                             0  --> only final result
!c                             >0 --> print modulus
!c     The matrices and dimension are communicated via
!c     common /radmatrix/.

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)


!INCLUDE 'twost.inc'
#ifdef REALSWITCH
!INCLUDE 'radmatrixr.inc'   ! defines also 'KDIM'
REAL(DP) :: q0(kdfull2,kstate)
REAL(DP) :: xlambda(kdim,kdim)
REAL(DP) :: symm_mat(kdim,kdim)
REAL(DP) :: utcond(kdim,kdim)
REAL*8 sum
REAL(DP) :: w(kdim)             ! workspace
REAL*8 vecovlpr   ! function names
REAL*8 vecnormr   ! function names
#else
!INCLUDE 'radmatrix.inc'   ! defines also 'KDIM'
COMPLEX(DP) :: q0(kdfull2,kstate)
COMPLEX(DP) :: xlambda(kdim,kdim)
COMPLEX(DP) :: symm_mat(kdim,kdim)
COMPLEX(DP) :: utcond(kdim,kdim)
COMPLEX(DP) :: sum
COMPLEX(DP) :: w(kdim)             ! workspace
COMPLEX(DP) :: vecovlp   ! function names
REAL*8 vecnorm   ! function names
COMPLEX(DP) :: wfovlp
#endif

INTEGER :: itmax2
LOGICAL :: ttest
PARAMETER (ttest=.false.)

REAL*8 variance,variance2  ! variance of step
REAL*8 radmax              ! max. squared radius
REAL*8 actstep,radvary
REAL(DP) :: varstate(kdim),averstate(kdim)


!-------------------------------------------------------

IF(iprint > 0) WRITE(6,'(a)') '  iter  variance '

IF(iter1 <= symutbegin .AND. is == 1)THEN  !!! computed for is=1 only because also contains is=2
#ifdef REALSWITCH
  CALL spmomsmatrixr(q0,-1)
#else
  CALL spmomsmatrix(q0,-1)
#endif
END IF

IF(iter1 <= symutbegin)THEN
  itmax2=9000
ELSE
  itmax2=50
END IF

DO iter=1,itmax2
!      write(*,*) ' ITER=',iter
  
  IF(iter1 <= symutbegin)THEN
#ifdef REALSWITCH
    CALL locutr(is,utcond,radmax,radvary)
  ELSE
    CALL symut_enr(is,utcond,q0)
    radmax=radmaxsym       !!! radmax computed more precisely for locut
#else
    CALL locut(is,utcond,radmax,radvary)
  ELSE
!      do nb=1,nstate
!        write(*,*) ' UTgradstep: nb,Q0 norm:',
!     &       nb,wfovlp(q0(1,nb),q0(1,nb))
!      enddo
    CALL symut_en(is,utcond,q0)
    radmax=radmaxsym       !!! radmax computed more precisely for locut
!        write(*,'(a,i5/100(1pg12.4))') ' vecs at iter=',iter,
!     &    ((vecs(na,nb,is),na=1,ndims(is)),nb=1,ndims(is))
!        write(*,'(a,i5/100(1pg12.4))') ' utcond at iter=',iter,
!     &    ((utcond(na,nb),na=1,ndims(is)),nb=1,ndims(is))
#endif
  END IF
  
!       matrix of Lagrangian multipliers (xlambda)
  
  DO i=1,ndims(is)
    DO j=1,ndims(is)
#ifdef REALSWITCH
      sum=0.D0
      DO n=1,ndims(is)
        sum=sum+vecsr(n,j,is)*utcond(n,i)
#else
        sum=CMPLX(0.D0,0.D0,DP)
        DO n=1,ndims(is)
          sum=sum+CONJG(vecs(n,j,is))*utcond(n,i)
#endif
        END DO
        xlambda(i,j)=sum
      END DO
    END DO
    DO i=1,ndims(is)
      DO j=i+1,ndims(is)
#ifdef REALSWITCH
        symm_mat(i,j) = xlambda(i,j)-xlambda(j,i)
        symm_mat(j,i) = symm_mat(i,j)
        xlambda(i,j) = 0.5D0*(xlambda(i,j)+xlambda(j,i))
        xlambda(j,i) = xlambda(i,j)
#else
        symm_mat(i,j) = xlambda(i,j)-CONJG(xlambda(j,i))
        symm_mat(j,i) = CONJG(symm_mat(i,j))
        xlambda(i,j) = 0.5D0*(xlambda(i,j)+CONJG(xlambda(j,i)))
        xlambda(j,i) = CONJG(xlambda(i,j))
#endif
      END DO
    END DO
    IF(ttest) THEN
      WRITE(6,'(a)') 'xlambda:'
      WRITE(6,'(20(1pg12.4))') ((xlambda(i,j),i=1,ndims(is)),j=1,ndims(is))
    END IF
    
!       application of radii-matrix - subtract constraint
    
    variance = 0.D0
    variance2 = 0.D0
    DO i=1,ndims(is)
      DO na=1,ndims(is)
#ifdef REALSWITCH
        sum=0.D0
        DO j=1,ndims(is)
          sum = sum + xlambda(i,j)*vecsr(na,j,is)
#else
          sum=CMPLX(0.D0,0.D0,DP)
          DO j=1,ndims(is)
            sum = sum + xlambda(i,j)*vecs(na,j,is)
#endif
          END DO
          w(na) = utcond(na,i) - sum
        END DO
        
!         convergence criteria
        
#ifdef REALSWITCH
        averstate(i) =  vecovlpr(vecsr(1,i,is),w,ndims(1))
        varstate(i) =  vecnormr(w,ndims(1))
#else
        averstate(i) =  REAL(vecovlp(vecs(1,i,is),w,ndims(1)))
        varstate(i) =  vecnorm(w,ndims(1))
#endif
        variance = variance + varstate(i)-averstate(i)**2
        variance2 = variance2 + varstate(i)
        
!         step weighted by radmax
        
        actstep = step/radmax
        DO na=1,ndims(is)
#ifdef REALSWITCH
          vecsr(na,i,is) = vecsr(na,i,is)-actstep*w(na)
#else
          vecs(na,i,is) = vecs(na,i,is)-actstep*w(na)
#endif
        END DO
      END DO
      IF(ttest) WRITE(6,'(a,i5,2(1pg13.5))')  &
          ' iter,variances=',iter,variance,variance2
      
!       ortho-normalization
      
#ifdef REALSWITCH
      CALL orthnormr(vecsr(1,1,is),ndims(is),kdim)
#else
!        write(*,'(a/100(1pg12.4))') ' vecs before ortho:',
!     &    ((vecs(na,nb,is),na=1,ndims(is)),nb=1,ndims(is))
      CALL orthnorm(vecs(1,1,is),ndims(is),kdim)
#endif
      
      IF(variance < precis) GO TO 99
      
!       if(iter1.gt.SymUtBegin)then
!       write(6,*)variance    !!! for tests
!       endif
      
    END DO
    99   CONTINUE
    
    IF(iprint >= 0) THEN
      WRITE(6,'(2(a,i4),a,1pg12.4)') 'UT cv is reached at it nr.',iter,  &
          ' for spin =',is,' with variance=',variance
!#ifdef COMPLEXSWITCH
!        write(*,'(a/100(1pg12.4))') ' vecs:',
!     &    ((vecs(na,nb,is),na=1,ndims(is)),nb=1,ndims(is))
!#endif
    END IF
    
    
    
!ccccccccccccccccc
! TESTS FOR THE UT
!ccccccccccccccccc
!       if(iter1.gt.SymUtBegin)then
    
!       write(6,*)'UT delta relation'
!      do na=1,ndims(is)
!      do nb=1,ndims(is)
!#ifdef REALSWITCH
!       sum=0.D0
!#else
!       sum=dcmplx(0.D0,0.D0)
!#endif
!       do i=1,ndims(is)
!#ifdef REALSWITCH
!       sum=vecsr(na,i,is)*vecsr(nb,i,is)+sum
!#else
!       sum=vecs(na,i,is)*conjg(vecs(nb,i,is))+sum
!#endif
!       enddo
!       write(6,*)na,nb,sum
!      enddo
!      enddo
    
!       if(is.eq.1)then
!       write(6,*)'UT matrix'
!       do na=1,ndims(1)
!       if(ispin(nrel2abs(na)).eq.is)then
!       do nb=1,ndims(1)
!       if(ispin(nrel2abs(nb)).eq.is)then
!#ifdef REALSWITCH
!       write(6,*)na,nb,vecsr(na,nb,is)
!#else
!       write(6,*)na,nb,vecs(na,nb,is)
!#endif
!       endif
!       enddo
!       endif
!       enddo
!       endif
    
!       endif
!ccccccccccc
    
    
    
    RETURN
#ifdef REALSWITCH
  END SUBROUTINE utgradstepr
#else
  END SUBROUTINE utgradstep
#endif
  
!-----------locut--------------------------------------
  
#ifdef REALSWITCH
  SUBROUTINE locutr(is,utcond,radmax,radvary)
#else
  SUBROUTINE locut(is,utcond,radmax,radvary)
#endif
 
USE params 
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
!  INCLUDE 'com.inc'
  
  
#ifdef REALSWITCH
  !INCLUDE 'radmatrixr.inc'   ! defines also 'KDIM'
  REAL(DP) :: xaver(kdim),yaver(kdim),zaver(kdim)
  REAL(DP) :: utcond(kdim,kdim)
  REAL*8 sum
  REAL*8 avermatrixr          ! function names
#else
  !INCLUDE 'radmatrix.inc'   ! defines also 'KDIM'
  COMPLEX(DP) :: xaver(kdim),yaver(kdim),zaver(kdim)
  COMPLEX(DP) :: utcond(kdim,kdim)
  COMPLEX(DP) :: sum
  COMPLEX(DP) :: avermatrix          ! function names
#endif
  
  LOGICAL :: ttest
  PARAMETER (ttest=.false.)
  
  REAL*8 radvary             ! variance of radii
  REAL*8 radmax              ! max. squared radius
  REAL(DP) :: rvary(kdim),raver(kdim)
  
!-------------------------------------------------------
  
!       averages of x,y,z
  
  radvary = 0.D0
  radmax = 0.D0
  DO i=1,ndims(is)
#ifdef REALSWITCH
    xaver(i) = avermatrixr(xmatr(1,1,is),vecsr(1,i,is), ndims(is),kdim)
    yaver(i) = avermatrixr(ymatr(1,1,is),vecsr(1,i,is), ndims(is),kdim)
    zaver(i) = avermatrixr(zmatr(1,1,is),vecsr(1,i,is), ndims(is),kdim)
    raver(i) = avermatrixr(rrmatr(1,1,is),vecsr(1,i,is), ndims(is),kdim)
    rvary(i) = raver(i)-xaver(i)**2 -yaver(i)**2-zaver(i)**2
#else
    xaver(i) = avermatrix(xmatr(1,1,is),vecs(1,i,is), ndims(is),kdim)
    yaver(i) = avermatrix(ymatr(1,1,is),vecs(1,i,is), ndims(is),kdim)
    zaver(i) = avermatrix(zmatr(1,1,is),vecs(1,i,is), ndims(is),kdim)
    raver(i) = REAL(avermatrix(rrmatr(1,1,is),vecs(1,i,is), ndims(is),kdim))
    rvary(i) = raver(i)-REAL(xaver(i))**2 -REAL(yaver(i))**2-REAL(zaver(i))**2
#endif
    radvary = radvary + rvary(i)
    radmax = MAX(radmax,raver(i))
  END DO
  
  IF(ttest) THEN
    WRITE(6,'(a,i3)') ' initial for spin IS=',is,  &
        '   state  x y z variance '
    WRITE(6,'(4(1pg12.4))') (xaver(i),yaver(i),zaver(i),SQRT(rvary(i)),  &
        i=1,ndims(is))
    WRITE(6,'(36x,1pg12.4)') SQRT(radvary/ndims(is))
    WRITE(6,'(a)') 'aver x y z  rr  rvary:'
    WRITE(6,'(5(1pg12.4))') (xaver(i),yaver(i),zaver(i),raver(i),rvary(i),  &
        i=1,ndims(is))
  END IF
  
! variational result of the UT constraint
  
  DO i=1,ndims(is)
    DO na=1,ndims(is)
#ifdef REALSWITCH
      sum=0.D0
#else
      sum=CMPLX(0.D0,0.D0,DP)
#endif
      DO nb=1,ndims(is)
#ifdef REALSWITCH
        sum = (rrmatr(na,nb,is)  &
        -2.D0 * xaver(i)*xmatr(na,nb,is)    &
              -2.D0 * yaver(i)*ymatr(na,nb,is)  &
              -2.D0 * zaver(i)*zmatr(na,nb,is)) *vecsr(nb,i,is) &
#else
          sum = (rrmatr(na,nb,is)  &
          -2.D0 * xaver(i)*xmatr(na,nb,is)  &
                -2.D0 * yaver(i)*ymatr(na,nb,is)  &
                -2.D0 * zaver(i)*zmatr(na,nb,is)) *vecs(nb,i,is) &
#endif  
            + sum
          END DO
          utcond(na,i)=sum
        END DO
      END DO
      
      RETURN
#ifdef REALSWITCH
    END SUBROUTINE locutr
#else
    END SUBROUTINE locut
#endif
    
!-----------symut_en------------------- !!! sym cond calculated with minimization of energy
    
#ifdef REALSWITCH
    SUBROUTINE symut_enr(is,utcond,q0)
#else
    SUBROUTINE symut_en(is,utcond,q0)
#endif
    
USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
    
    LOGICAL :: ttest
    PARAMETER (ttest=.false.)
    
    REAL(DP) :: usicsp(2*kdfull2),rhosp(2*kdfull2)
    REAL*8 save1,save2,sum2
    
#ifdef REALSWITCH
    !INCLUDE 'radmatrixr.inc'   ! defines also 'KDIM'
    REAL(DP) :: q0(kdfull2,kstate)
    REAL(DP) :: qsym(kdfull2,kstate)
    REAL(DP) :: uqsym(kdfull2,kstate)
    REAL(DP) :: symcond(kstate,kstate)
    REAL(DP) :: utcond(kdim,kdim)
#else
    !INCLUDE 'radmatrix.inc'   ! defines also 'KDIM'
    COMPLEX(DP) :: q0(kdfull2,kstate)
    COMPLEX(DP) :: qsym(kdfull2,kstate)
    COMPLEX(DP) :: uqsym(kdfull2,kstate)
    COMPLEX(DP) :: symcond(kstate,kstate)
    COMPLEX(DP) :: utcond(kdim,kdim)
    COMPLEX(DP) :: wfovlp ! function declaration
#endif
    
!-------------------------------------------------------
    
    ishift = (is-1)*nxyz
    
    DO nb=1,nstate
      IF(ispin(nrel2abs(nb)) == is)THEN
        nbeff = nb - (is-1)*ndims(1)
        
#ifdef REALSWITCH
        CALL superpose_stater(qsym(1,nb),vecsr(1,nbeff,is),q0,is)
#else
        CALL superpose_state(qsym(1,nb),vecs(1,nbeff,is),q0,is)
#endif
        
!      write(*,*) ' Q0 norm:',wfovlp(q0(1,nb),q0(1,nb))
!      write(*,*) ' QSYM norm:',wfovlp(qsym(1,nb),qsym(1,nb))
        save1=enrear      !!! to deactivate accumulation of enrear, enerpw
        save2=enerpw      !!! (else wrong total energy)
#ifdef REALSWITCH
        CALL calc_sicspr(rhosp,usicsp,qsym(1,nb),nb)
#else
        CALL calc_sicsp(rhosp,usicsp,qsym(1,nb),nb)
#endif
        enrear=save1
        enerpw=save2
        
        DO ind=1,nxyz
          uqsym(ind,nb) = usicsp(ind+ishift)*qsym(ind,nb)
        END DO
        
! variational result of the UT constraint
        
        DO na=1,nstate
          IF(ispin(nrel2abs(na)) == is)THEN
            naa = na - (is-1)*ndims(1)
#ifdef REALSWITCH
            utcond(naa,nbeff) = -rwfovlp(q0(1,na),uqsym(1,nb))
#else
            utcond(naa,nbeff) = -wfovlp(q0(1,na),uqsym(1,nb))
#endif
          END IF
        END DO
        
!        do na=1,nstate
!        if(ispin(nrel2abs(na)).eq.is)then
!#ifdef REALSWITCH
!           symcond(nb,na) = rwfovlp(Uqsym(1,nb),qsym(1,na))
!     &                      - rwfovlp(qsym(1,nb),Uqsym(1,na))
        
!#else
!           symcond(nb,na) = wfovlp(Uqsym(1,nb),qsym(1,na))
!     &                      - wfovlp(qsym(1,nb),Uqsym(1,na))
        
!#endif
!        endif
!        enddo
        
      END IF
    END DO
    
!c
!c Sym Cond calculation
!c
!       sum2 = 0.D0
!       do na=1,nstate
!       if(ispin(na).eq.is)then
!         do nb=1,nstate
!         if(ispin(nb).eq.is .and. nb.ne.na)then
!#ifdef REALSWITCH
!         sum2 = symcond(na,nb)**2 + sum2
!#else
!         sum2 = abs(symcond(na,nb))**2 + sum2
!#endif
!         endif
!         enddo
!       endif
!       enddo
!       write(6,'(a,(1pg12.4))')
!     & 'Sym UT: evol of SymCond ',sqrt(sum2)
    
    RETURN
#ifdef REALSWITCH
  END SUBROUTINE symut_enr
#else
  END SUBROUTINE symut_en
#endif
  
!-----avermatrix-----------------------------------------
  
#ifdef REALSWITCH
  REAL(DP) FUNCTION avermatrixr(a,vec,ndim,kdim)
#else
  COMPLEX(DP) FUNCTION avermatrix(a,vec,ndim,kdim)
#endif
  USE params, ONLY: DP
  
!      implicit none
  
!     Average of matrix 'a' taken with state 'vec' of
!     actual length 'ndim', dimensioned with 'kdim'.
  
  INTEGER :: ndim,kdim
  INTEGER :: i,j
#ifdef REALSWITCH
  REAL(DP) :: a(kdim,kdim),vec(kdim)
  REAL(DP) :: aver
#else
  COMPLEX(DP) :: a(kdim,kdim),vec(kdim)
  COMPLEX(DP) :: aver
#endif
  
!-------------------------------------------------------
  
#ifdef REALSWITCH
  aver = 0.D0
#else
  aver = CMPLX(0.D0,0.D0,DP)
#endif
  DO i=1,ndim
    DO j=1,ndim
#ifdef REALSWITCH
      aver = aver + vec(j)*a(j,i)*vec(i)
#else
      aver = aver + CONJG(vec(j))*a(j,i)*vec(i)
#endif
    END DO
  END DO
#ifdef REALSWITCH
  avermatrixr = aver
#else
  avermatrix = aver
#endif
  
  RETURN
#ifdef REALSWITCH
END FUNCTION avermatrixr
#else
END FUNCTION avermatrix
#endif

!-----ovlpmatrix-----------------------------------------

#ifdef REALSWITCH
REAL(DP) FUNCTION ovlpmatrixr(a,vec1,vec2,ndim,kdim)
#else
COMPLEX(DP) FUNCTION ovlpmatrix(a,vec1,vec2,ndim,kdim)
#endif
USE params, ONLY: DP

!      implicit none

!     Matrix element of matrix 'a'
!     with respect to states 'vec1' and 'vec2'
!     having actual length 'ndim' and dimension 'kdim'.

INTEGER :: ndim,kdim
INTEGER :: i,j
#ifdef REALSWITCH
REAL(DP) :: a(kdim,kdim),vec1(kdim),vec2(kdim)
REAL(DP) :: ovlp
#else
COMPLEX(DP) :: a(kdim,kdim),vec1(kdim),vec2(kdim)
COMPLEX(DP) :: ovlp

#endif

!-------------------------------------------------------

#ifdef REALSWITCH
ovlp = 0.D0
#else
ovlp = CMPLX(0.D0,0.D0,DP)
#endif
DO i=1,ndim
  DO j=1,ndim
#ifdef REALSWITCH
    ovlp = ovlp + vec1(j)*a(j,i)*vec2(i)
#else
    ovlp = ovlp + CONJG(vec1(j))*a(j,i)*vec2(i)
#endif
  END DO
END DO
#ifdef REALSWITCH
ovlpmatrixr = ovlp
#else
ovlpmatrix = ovlp
#endif

RETURN
#ifdef REALSWITCH
END FUNCTION ovlpmatrixr
#else
END FUNCTION ovlpmatrix
#endif

!-----orthnorm-----------------------------------------

#ifdef REALSWITCH
SUBROUTINE orthnormr(vecs,ndim,kdim)
#else
SUBROUTINE orthnorm(vecs,ndim,kdim)
#endif
USE params, ONLY: DP

!      implicit none

!     ortho-normalizes system of vectors 'vecs' with
!     dimension 'ndim'.

INTEGER :: ndim,kdim
INTEGER :: i,j,n
REAL(DP) :: sum2
#ifdef REALSWITCH
REAL(DP) :: sum
REAL(DP) :: vecs(kdim,kdim)
REAL(DP) :: vecovlpr
REAL(DP) :: vecnormr
#else
COMPLEX(DP) :: sum
COMPLEX(DP) :: vecs(kdim,kdim)
COMPLEX(DP) :: vecovlp
REAL(DP) :: vecnorm
#endif

!-------------------------------------------------------

DO i=1,ndim
  IF(i > 1) THEN
    DO j=1,i-1
#ifdef REALSWITCH
      sum = vecovlpr(vecs(1,j),vecs(1,i),ndim)
#else
      sum = vecovlp(vecs(1,j),vecs(1,i),ndim)
#endif
      DO n=1,ndim
        vecs(n,i) = vecs(n,i)-sum*vecs(n,j)
      END DO
    END DO
  END IF
#ifdef REALSWITCH
  sum2 = 1.D0/SQRT(vecnormr(vecs(1,i),ndim))
#else
  sum2 = 1.D0/SQRT(vecnorm(vecs(1,i),ndim))
#endif
  DO n=1,ndim
    vecs(n,i) = vecs(n,i)*sum2
  END DO
END DO

RETURN
#ifdef REALSWITCH
END SUBROUTINE orthnormr
#else
END SUBROUTINE orthnorm
#endif

#ifdef REALSWITCH
REAL(DP) FUNCTION vecnormr(vec,ndim)
#else
REAL(DP) FUNCTION vecnorm(vec,ndim)
#endif
USE params, ONLY: DP

!      implicit none

!     ortho-normalizes system of vectors 'vecs' with
!     dimension 'ndim'.

INTEGER :: ndim
INTEGER :: n
#ifdef REALSWITCH
REAL(DP) :: vec(ndim)
#else
COMPLEX(DP) :: vec(ndim)
#endif

REAL(DP) :: sum

!-------------------------------------------------------

sum = 0.D0
DO n=1,ndim
#ifdef REALSWITCH
  sum = sum + vec(n)**2
#else
  sum = sum + ABS(vec(n))**2
#endif
END DO
#ifdef REALSWITCH
vecnormr = sum
#else
vecnorm = sum
#endif

RETURN
#ifdef REALSWITCH
END FUNCTION vecnormr
#else
END FUNCTION vecnorm
#endif

!-----vecovlp----------------------------------------------

#ifdef REALSWITCH
REAL(DP) FUNCTION vecovlpr(vec1,vec2,ndim)
#else
COMPLEX(DP) FUNCTION vecovlp(vec1,vec2,ndim)
#endif
USE params, ONLY: DP

!      implicit none

!     Overlap 'vec1' with 'vec2' having dimension 'ndim'.

INTEGER :: ndim
INTEGER :: n
#ifdef REALSWITCH
REAL(DP) :: vec1(ndim),vec2(ndim)
REAL(DP) :: sum
#else
COMPLEX(DP) :: vec1(ndim),vec2(ndim)
COMPLEX(DP) :: sum
#endif

!-------------------------------------------------------

#ifdef REALSWITCH
sum = 0.D0
#else
sum = CMPLX(0.D0,0.D0,DP)
#endif
DO n=1,ndim
#ifdef REALSWITCH
  sum = sum + vec1(n)*vec2(n)
#else
  sum = sum + CONJG(vec1(n))*vec2(n)
#endif
END DO
#ifdef REALSWITCH
vecovlpr = sum
#else
vecovlp = sum
#endif

RETURN
#ifdef REALSWITCH
END FUNCTION vecovlpr
#else
END FUNCTION vecovlp
#endif

!-----superpose_state--------------------------------------------------

#ifdef REALSWITCH
SUBROUTINE superpose_stater(wfsup,coeff,q0,is)
#else
SUBROUTINE superpose_state(wfsup,coeff,q0,is)
#endif

!     Superposition to new state:
!       wfsup     new single-particle state
!       coeff     superposition coefficients
!       q0        set of s.p.states to be combined
!       is        spin of states

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)


#ifdef REALSWITCH
!INCLUDE 'radmatrixr.inc'
REAL(DP) :: q0(kdfull2,kstate)
REAL(DP) :: wfsup(kdfull2)
REAL(DP) :: coeff(kstate)
#else
!INCLUDE 'radmatrix.inc'
COMPLEX(DP) :: q0(kdfull2,kstate)
COMPLEX(DP) :: wfsup(kdfull2)
COMPLEX(DP) :: coeff(kstate)
COMPLEX(DP) :: wfovlp,orbitaloverlap
#endif

!---------------------------------------------------------------------

!      write(*,*) ' SUPERPOSE: NDIMS=',ndims
!      write(*,*) ' SUPERPOSE: COEFF=',(coeff(i),i=1,ndims(is))
DO nas=1,ndims(is)
  na = nas + (is-1)*ndims(1)
!#ifdef COMPLEXSWITCH
!        write(*,*) ' Q0 norm at nas=',nas,
!     &     wfovlp(q0(1,na),q0(1,na))
!#endif
  IF(nas == 1) THEN
    DO i=1,nxyz
      wfsup(i) = coeff(nas)*q0(i,na)
    END DO
  ELSE
    DO i=1,nxyz
      wfsup(i) = coeff(nas)*q0(i,na) + wfsup(i)
    END DO
  END IF
END DO

RETURN
#ifdef REALSWITCH
END SUBROUTINE superpose_stater
#else
END SUBROUTINE superpose_state
#endif

!-----spmomsmatrix----------------------------------------------spmoms

#ifdef REALSWITCH
SUBROUTINE spmomsmatrixr(wfr,PRINT)
#else
SUBROUTINE spmomsmatrix(wfr,PRINT)
#endif

!     Matrix of spatial moments between single-particle states
!     from real  wf's:
!      wfr    = set of real single particle wavefunctions
!     The result is stored in common/radmatrix/ for further
!     use in localization transformation.

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)


INTEGER :: iunit,PRINT    ! print =1 : printing of the total variance
PARAMETER (iunit=0)    ! set zero to disable testprint

REAL :: var
LOGICAL :: tfirst
DATA tfirst/.true./

#ifdef REALSWITCH
!INCLUDE 'radmatrixr.inc'
REAL(DP) :: wfr(kdfull2,kstate)
#else
!INCLUDE 'radmatrix.inc'
COMPLEX(DP) :: s,wfmom,xmom,ymom,zmom,xxmom,yymom,zzmom
COMPLEX(DP) :: wfr(kdfull2,kstate)
#endif

!----------------------------------------------------------------------

OPEN(iunit,POSITION='append',FILE='pstatmomsmatrix.'//outnam)

!     check spin of states

!      ndims(1) = 0
!      ndims(2) = 0
!      do na=1,nstate
!        if(ispin(nrel2abs(na)).eq.1) then
!          if(ndims(2).ne.0) stop ' spins out of order'
!          ndims(1) = na
!        else
!          if(ndims(1).eq.0) stop ' spins out of order'
!          ndims(2) = na-ndims(1)
!        endif
!      enddo

noff = ndims(1)

!      if(iunit.gt.0) write(iunit,'(a,2i5)') 'NDIMS=',ndims

!     compute matrix elements and store

IF(iunit > 0) THEN
  WRITE(iunit,'(a)') 'matrix of s.p. moments:',  &
      '   na   nb      x     y      z     xx     yy     zz'
  tfirst = .false.
END IF
DO na=1,nstate
  
  DO nb=1,na
    IF(ispin(nrel2abs(nb)) == ispin(nrel2abs(na))) THEN
#ifdef REALSWITCH
      wfmom = 0D0
      xmom = 0D0
      ymom = 0D0
      zmom = 0D0
      xxmom = 0D0
      yymom = 0D0
      zzmom = 0D0
#else
      wfmom = CMPLX(0.D0,0.D0,DP)
      xmom = CMPLX(0.D0,0.D0,DP)
      ymom = CMPLX(0.D0,0.D0,DP)
      zmom = CMPLX(0.D0,0.D0,DP)
      xxmom = CMPLX(0.D0,0.D0,DP)
      yymom = CMPLX(0.D0,0.D0,DP)
      zzmom = CMPLX(0.D0,0.D0,DP)
#endif
      
      ind=0
      DO iz=minz,maxz
        z1=(iz-nzsh)*dz
        z2=z1*z1
        DO iy=miny,maxy
          y1=(iy-nysh)*dy
          y2=y1*y1
          DO ix=minx,maxx
            ind=ind+1
            IF((ix /= nx2).AND.(iy /= ny2).AND.(iz /= nz2)) THEN
              x1=(ix-nxsh)*dx
              x2=x1*x1
#ifdef REALSWITCH
              s=wfr(ind,na)*wfr(ind,nb)
#else
              s=CONJG(wfr(ind,na))*wfr(ind,nb)
#endif
              wfmom=wfmom+s       ! monopole
              xmom=xmom+s*x1      ! dipole
              ymom=ymom+s*y1
              zmom=zmom+s*z1
              xxmom=xxmom+s*x2    ! quadrupole
              yymom=yymom+s*y2
              zzmom=zzmom+s*z2
            END IF
          END DO
        END DO
      END DO
      
      wfmom = dvol*wfmom
      xmom = dvol*xmom
      ymom = dvol*ymom
      zmom = dvol*zmom
      xxmom = dvol*xxmom
      yymom = dvol*yymom
      zzmom = dvol*zzmom
      
      IF(iunit > 0) WRITE(iunit,'(3i4,6(1pg13.5))')  &
          na,nb,ispin(nrel2abs(nb)),xmom,ymom,zmom, xxmom,yymom,zzmom
      
      is = ispin(nrel2abs(nb))
      
      IF(is == 1) THEN
        ma = na
        mb = nb
      ELSE
        ma = na-noff
        mb = nb-noff
      END IF
      IF(iunit > 0) WRITE(iunit,*) 'ma,mb=',ma,mb,na,nb
      xmatr(ma,mb,is) = xmom
      ymatr(ma,mb,is) = ymom
      zmatr(ma,mb,is) = zmom
      xxmatr(ma,mb,is) = xxmom
      yymatr(ma,mb,is) = yymom
      zzmatr(ma,mb,is) = zzmom
#ifdef REALSWITCH
      xmatr(mb,ma,is) = xmom
      ymatr(mb,ma,is) = ymom
      zmatr(mb,ma,is) = zmom
      xxmatr(mb,ma,is) = xxmom
      yymatr(mb,ma,is) = yymom
      zzmatr(mb,ma,is) = zzmom
#else
      xmatr(mb,ma,is) = CONJG(xmom)
      ymatr(mb,ma,is) = CONJG(ymom)
      zmatr(mb,ma,is) = CONJG(zmom)
      xxmatr(mb,ma,is) = CONJG(xxmom)
      yymatr(mb,ma,is) = CONJG(yymom)
      zzmatr(mb,ma,is) = CONJG(zzmom)
#endif
    END IF
  END DO
END DO

!      write(6,'(10(/3i4,6(1pg13.5)))')
!     &   ((na,nb,is,
!     &    xmatr(na,nb,is),ymatr(na,nb,is),zmatr(na,nb,is),
!     &    xxmatr(na,nb,is),yymatr(na,nb,is),
!     &    zzmatr(na,nb,is),
!     &    na=1,ndims(is)),nb=1,ndims(is))

!     initialize eigen-vectors

DO is=1,2
  DO na=1,ndims(is)
    DO nb=1,ndims(is)
      rrmatr(na,nb,is) = xxmatr(na,nb,is) + yymatr(na,nb,is)  &
          + zzmatr(na,nb,is)
    END DO
  END DO
END DO

IF(PRINT == 1)THEN
  DO is=1,2
#ifdef REALSWITCH
    var=0.D0
#else
    var=CMPLX(0.D0,0.D0,DP)
#endif
    DO na=1,ndims(is)
      var = var + rrmatr(na,na,is) - xmatr(na,na,is)**2 -  &
          ymatr(na,na,is)**2 - zmatr(na,na,is)**2
    END DO
    WRITE(6,*)'For spin =',is,' Total spacial variance ',var
  END DO
END IF

IF(iunit > 0) CLOSE(iunit)

RETURN
#ifdef REALSWITCH
END SUBROUTINE spmomsmatrixr
#else
END SUBROUTINE spmomsmatrix
#endif


!#endif


#ifdef REALSWITCH
SUBROUTINE unitary_gradstepr(symm_mat,vect,ndima,norder)
#else
SUBROUTINE unitary_gradstep(symm_mat,vect,ndima,norder)
#endif
!
!  Unitary extrapolation of the gradient step.
!    symm_mat   = matrix of symmetry condition, rescaled by stepsize
!    vect       = matrix of state vectors to be propagated
!    ndima      = actual dimension of matrix
!    norder     = order of the Taylor expansion of the exponential
!
USE params
IMPLICIT REAL(DP) (A-H,O-Z)

!INTEGER,PARAMETER :: kdim=kstate

INTEGER,INTENT(IN) :: norder,ndima

#ifdef REALSWITCH
REAL(DP),INTENT(IN) :: symm_mat(kdim,kdim)
REAL(DP),INTENT(INOUT) :: vect(kdim,kdim)
REAL(DP) :: vectwork(kdim,kdim)
REAL(DP) :: symm_acc(kdim,kdim)
REAL(DP) :: sumv
REAL(DP) :: w(kdim)             ! workspace
#else
COMPLEX(DP),INTENT(IN) :: symm_mat(kdim,kdim)
COMPLEX(DP),INTENT(INOUT) :: vect(kdim,kdim)
COMPLEX(DP) :: vectwork(kdim,kdim)
COMPLEX(DP) :: symm_acc(kdim,kdim)
COMPLEX(DP) :: sumv
COMPLEX(DP) :: w(kdim)             ! workspace
#endif

INTEGER :: i,j,k,n      ! counters in DO loops
REAL(DP) :: facinv

DO i=1,ndima; DO j=1,ndima
  vectwork(i,j) = vect(i,j)
END DO; END DO

facinv = 1D0
DO n=1,norder

  facinv = facinv/n

  DO k=1,ndima
    w = 0D0
    DO i=1,ndima
      DO j=1,ndima
        w(i) = w(i) - symm_mat(i,j)*vectwork(k,j)
      END DO
    END DO
    DO i=1,ndima
      vectwork(k,i) = w(i)*facinv
    END DO
  END DO

  DO i=1,ndima; DO j=1,ndima
    vect(i,j) = vect(i,j)+vectwork(i,j)
  END DO; END DO

END DO

RETURN

#ifdef REALSWITCH
END SUBROUTINE unitary_gradstepr
#else
END SUBROUTINE unitary_gradstep
#endif



#ifdef REALSWITCH
END MODULE twostr
#endif

#ifdef COMPLEXSWITCH
END MODULE twost
#endif

#else
#ifdef REALSWITCH
!-----init_fsic------------------------------------------------

SUBROUTINE init_fsic()
STOP ' code not compiled for GSlat or double-set SIC'
RETURN
END SUBROUTINE init_fsic
#else
SUBROUTINE ccc  ! an unuseful subroutine so that the compilation does not abort IF twostsic=0
RETURN
END SUBROUTINE ccc  ! an unuseful subroutine so that
#endif

#endif

