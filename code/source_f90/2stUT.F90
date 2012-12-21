#include "define.h"
#if(twostsic)

#ifdef REALSWITCH
MODULE twostr
USE params
USE kinetic
USE localize_rad
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP),PRIVATE,ALLOCATABLE :: qnewr(:,:)
REAL(DP),PRIVATE,ALLOCATABLE :: psirut(:,:)
!COMMON /twost/ qnewr,qnew,psirut,psiut

REAL(DP) :: step=0.1D0,precis=1D-6,radmaxsym=1D0
REAL(DP) :: ener_2st(2)
INTEGER :: symutbegin=100,itut
!COMMON /twostut/ step,precis,symutbegin,radmaxsym,itut

INTEGER,PRIVATE :: kdim
!     matrices of radial moments

REAL(DP),ALLOCATABLE :: rExpDABold(:,:,:)!MV added
REAL(DP),PRIVATE,ALLOCATABLE :: rrmatr(:,:,:)  ! matrix of r**2
REAL(DP),PRIVATE,ALLOCATABLE :: xxmatr(:,:,:)  ! matrix of x**2
REAL(DP),PRIVATE,ALLOCATABLE :: yymatr(:,:,:)  ! matrix of y**2
REAL(DP),PRIVATE,ALLOCATABLE :: zzmatr(:,:,:)  ! matrix of z**2
REAL(DP),PRIVATE,ALLOCATABLE :: xmatr(:,:,:)   ! matrix of x
REAL(DP),PRIVATE,ALLOCATABLE :: ymatr(:,:,:)   ! matrix of y
REAL(DP),PRIVATE,ALLOCATABLE :: zmatr(:,:,:)   ! matrix of z
REAL(DP),ALLOCATABLE :: vecsr(:,:,:)    ! searched eigenvevtors
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

COMPLEX(DP),ALLOCATABLE :: qnewut(:,:)
COMPLEX(DP),ALLOCATABLE :: psiut(:,:)

!INTEGER,PRIVATE :: kdim=kstate
INTEGER,PRIVATE :: kdim

!     matrices of radial moments

COMPLEX(DP),ALLOCATABLE :: ExpDABold(:,:,:)!MV added
COMPLEX(DP),PRIVATE,ALLOCATABLE :: rrmatr(:,:,:)  ! matrix of r**2
COMPLEX(DP),PRIVATE,ALLOCATABLE :: xxmatr(:,:,:)  ! matrix of x**2
COMPLEX(DP),PRIVATE,ALLOCATABLE :: yymatr(:,:,:)  ! matrix of y**2
COMPLEX(DP),PRIVATE,ALLOCATABLE :: zzmatr(:,:,:)  ! matrix of z**2
COMPLEX(DP),PRIVATE,ALLOCATABLE :: xmatr(:,:,:)   ! matrix of x
COMPLEX(DP),PRIVATE,ALLOCATABLE :: ymatr(:,:,:)   ! matrix of y
COMPLEX(DP),PRIVATE,ALLOCATABLE :: zmatr(:,:,:)   ! matrix of z
COMPLEX(DP),ALLOCATABLE :: vecs(:,:,:)    ! searched eigenvevtors
!COMMON /radmatrix/ rrmatr,xxmatr,yymatr,zzmatr, xmatr,ymatr,zmatr,  &
!    vecs
!COMMON /radmatrixn/ndims

#endif

CONTAINS

#ifdef REALSWITCH
!-----init_fsicr------------------------------------------------

SUBROUTINE init_fsicr()

!     initializes fields for FSIC etc

USE params
USE kinetic
USE symcond
IMPLICIT REAL(DP) (A-H,O-Z)

LOGICAL,PARAMETER :: tconv=.true.       ! to print convergence

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

!     dimensions of spin sub-spaces

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


! initialize work space
ALLOCATE(usicall(kdfull2,kstate))
ALLOCATE(qnewr(kdfull2,kstate))
ALLOCATE(psirut(kdfull2,kstate))

kdim=kstate
ALLOCATE(rExpDABold(kstate, kstate, 2))!MV added
ALLOCATE(rrmatr(kdim,kdim,2))  ! matrix of r**2
ALLOCATE(xxmatr(kdim,kdim,2))  ! matrix of x**2
ALLOCATE(yymatr(kdim,kdim,2))  ! matrix of y**2
ALLOCATE(zzmatr(kdim,kdim,2))  ! matrix of z**2
ALLOCATE(xmatr(kdim,kdim,2))   ! matrix of x
ALLOCATE(ymatr(kdim,kdim,2))   ! matrix of y
ALLOCATE(zmatr(kdim,kdim,2))   ! matrix of z
ALLOCATE(vecsr(kdim,kdim,2))   ! searched eigenvevtors

!     initialize unitary matrix

DO is=1,2
  DO i=1,ndims(is)
    DO j=1,ndims(is)
      IF(i == j) THEN
        vecsr(i,j,is) = 1D0
      ELSE
        vecsr(i,j,is) = 0D0
      END IF
    END DO
  END DO
END DO

! initialize protocolo file
IF(tconv) THEN
OPEN(353,file='2st-stat-conv.res')
REWIND(353)
WRITE(353,'(a)') '# protocol of convergence if symmetry condition'
CLOSE(353)
END IF


RETURN
END SUBROUTINE init_fsicr
!-----end_fsicr------------------------------------------------

SUBROUTINE end_fsicr()

USE params
USE kinetic
USE symcond
IMPLICIT REAL(DP) (A-H,O-Z)

!   frees workspace for static SIC

DEALLOCATE(rExpDABold)!MV added
DEALLOCATE(rrmatr)  ! matrix of r**2
DEALLOCATE(xxmatr)  ! matrix of x**2
DEALLOCATE(yymatr)  ! matrix of y**2
DEALLOCATE(zzmatr)  ! matrix of z**2
DEALLOCATE(xmatr)   ! matrix of x
DEALLOCATE(ymatr)   ! matrix of y
DEALLOCATE(zmatr)   ! matrix of z
DEALLOCATE(vecsr)   ! searched eigenvevtors


RETURN
END SUBROUTINE end_fsicr

#endif



#ifdef COMPLEXSWITCH
!-----init_fsic------------------------------------------------

SUBROUTINE init_fsic()

!     initializes fields for FSIC etc

USE params
USE kinetic
USE symcond
IMPLICIT REAL(DP) (A-H,O-Z)

!INCLUDE "twost.inc"
!INCLUDE 'radmatrixr.inc'
NAMELIST /fsic/step,precis,symutbegin,radmaxsym   !!! UT parameters

!----------------------------------------------------------------

kdim=kstate
ALLOCATE(ExpDABold(kstate, kstate, 2))!MV added
ALLOCATE(rrmatr(kdim,kdim,2))  ! matrix of r**2
ALLOCATE(xxmatr(kdim,kdim,2))  ! matrix of x**2
ALLOCATE(yymatr(kdim,kdim,2))  ! matrix of y**2
ALLOCATE(zzmatr(kdim,kdim,2))  ! matrix of z**2
ALLOCATE(xmatr(kdim,kdim,2))   ! matrix of x
ALLOCATE(ymatr(kdim,kdim,2))   ! matrix of y
ALLOCATE(zmatr(kdim,kdim,2))   ! matrix of z
ALLOCATE(vecs(kdim,kdim,2))    ! searched eigenvevtors

ALLOCATE(qnewut(kdfull2,kstate))
ALLOCATE(psiut(kdfull2,kstate))


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

SUBROUTINE init_vecs()

!     Initializes complex superposition coefficients for SIC.

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
!INCLUDE 'twost.inc'
!INCLUDE 'radmatrixr.inc'
!COMPLEX(DP) :: vecs(kdim,kdim,2)    ! searched eigenvevtors

!----------------------------------------------------------------

WRITE(6,*) 'vecsr initvecs'!MV
WRITE (6,'(3f12.6)') ((vecsr(ii,jj,1), ii=1,3),jj=1,3)!MV

DO i=1,kstate
  DO j=1,kstate
    vecs(i,j,1) = CMPLX(vecsr(i,j,1),0D0,DP)
    vecs(i,j,2) = CMPLX(vecsr(i,j,2),0D0,DP)
  END DO
END DO

WRITE(6,*) 'vecsr fin initvecs'!MV
WRITE (6,'(6f12.6)') ((vecs(ii,jj,1), ii=1,3),jj=1,3)!MV

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
    acc = 0D0
    DO na=1,nstate
      IF(ispin(na) == is)THEN
        DO nb=1,nstate
          IF(ispin(nb) == is)THEN
            acc = ( rwfovlp(psirut(1,na),qnewr(1,nb)) -  &
                rwfovlp(psirut(1,nb),qnewr(1,na)) )**2 + acc
          END IF
        END DO
      END IF
    END DO
    WRITE(6,'(a,i3,a,(1pg12.4))')  &
        'For spin',is,'  Total violation of SymCond',SQRT(acc)
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
  
! lagrange mult matrix diag (necessits vect1,2 in each sp subs, instead bugs)
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
#ifdef COMPLEXSWITCH
!-----calc_utwfc--------------------------------------------------!MV

SUBROUTINE calc_utwfc(q0,q0ut,iter1)


!     computes localized wavefunctions
!       input is set of wavefunctions on 'q0'
!       output are localized wavefunctions 'q0UT'
!       the array 'qnewUT' is used as workspace

!#INCLUDE "all.inc"

!     basic arrays and workspace

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
!INCLUDE 'twost.inc'
!INCLUDE 'radmatrix.inc'
!INCLUDE 'vec.inc'!MV

COMPLEX(DP), INTENT(IN OUT) :: q0(kdfull2,kstate)
COMPLEX(DP), INTENT(IN OUT) :: q0ut(kdfull2,kstate)
COMPLEX(DP) :: wfovlp

!------------------------------------------------------------------

!     Compute matrix of radial moments and determine
!     optimally localizing transformation.
!     Results is unitary matrix 'vecs' communicated via
!     common /radmatrix/.

DO is=1,2
  IF(ndims(is) > 1)THEN
    CALL utgradstepc(is,0,q0,iter1)   !!!!! new vecs (gradient method)!MV
  END IF
END DO

DO nb=1,nstate
  is = ispin(nrel2abs(nb))
  nbeff = nb - (is-1)*ndims(1)
  CALL superpose_state(q0ut(1,nb),vecs(1,nbeff,is),q0,is)
END DO

RETURN
END SUBROUTINE calc_utwfc
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
!COMPLEX(DP) :: wfovlp
#endif

!------------------------------------------------------------------

!     Compute matrix of radial moments and determine
!     optimally localizing transformation.
!     Results is unitary matrix 'vecs' communicated via
!     common /radmatrix/.

DO is=1,2
  IF(ndims(is) > 1)THEN
#ifdef REALSWITCH
    CALL utgradstepr(is,0,q0,iter1)   !!!!! new vecs (gradient method)
#else
!      do nb=1,nstate
!        write(*,*) ' CALC_UTWF: nb,Q0 norm:',
!     &       nb,wfovlp(q0(1,nb),q0(1,nb))
!      enddo
    CALL utgradstepc(is,0,q0,iter1)   !!!!! new vecs (gradient method)
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

!-----utgradstepc------------------------------------------------------------------------

#ifdef COMPLEXSWITCH
SUBROUTINE utgradstepc(is,iprint,q0,iter1)

!c     Nonlinear gradient iteration to optmially localized states:
!c      'vecs'    system of eigen-vectors to be determined
!c      'is'      isospin
!c      'iprint'  print level: <0 --> no print at all
!c                             0  --> only final result
!c                             >0 --> print modulus

!#INCLUDE "all.inc"
USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
!INCLUDE 'twost.inc'
!INCLUDE 'radmatrix.inc'   ! defines also 'KDIM'
!INCLUDE 'vec.inc'!MV
COMPLEX(DP) :: q0(kdfull2,kstate)
COMPLEX(DP) :: xlambda(kdim,kdim)
COMPLEX(DP) :: acc
COMPLEX(DP) :: dab(kdim,kdim),expdab(kdim,kdim)           !MV! workspace
COMPLEX(DP) :: dabsto(kdim,kdim)          !MV! workspace
COMPLEX(DP) :: vecovlpc   ! function names
REAL(DP) :: vecnorm   ! function names
REAL(DP) :: matdorth !MV function
REAL(DP) :: matnorme !MV function
COMPLEX(DP) :: wfovlp

INTEGER :: itmax2,ni

REAL(DP) :: norm, ERR_c           ! variance of step, erreur du produit
!REAL(DP) :: radmax              ! max. squared radius
REAL(DP) :: actstep,radvary
!REAL(DP) :: varstate(kdim),averstate(kdim)

!-------------------------------------------------------
itmax2=symutbegin ! 50
ni=ndims(is)

! update unitary transformation 'vecs' with previous exponentional
CALL matmult (vecs(1,1, is),expdabold(1,1,is),dabsto,kdim,ni)
CALL matcopy(dabsto,vecs(1,1,is),kdim,ni)

actstep = step/radmaxsym     ! radmaxsym obsolete, set to 1D0

DO iter=1,itmax2
  CALL dalphabeta(is, dab, q0)     !actual symmetry condition on 'dab'
  norm=matnorme(dab,kdim,ni)       !norm of symm.cond. = error
  CALL matconst(dab,dab,CMPLX(-actstep,0),kdim,ni)  !mutiply 'dab' stepsize
  CALL matexp(dab,expdab,kdim,ni)  ! exponentiate 'dab' to 'ExpDab'
  CALL matabtoa(expdabold(1,1,is),expdab,kdim,ni)   !store 'expdab'
  CALL matabtoa(vecs(1,1,is), expdab,kdim,ni)       !update 'vecs'
  CALL matcnjg(expdabold(1,1,is),dab,kdim,ni)       ! conjugate on 'dab'
  CALL matmult(expdabold(1,1,is), dab,dabsto,kdim,ni)
  ERR_c=matnorme(dabsto,kdim,ni)
  WRITE(6,'(a,4f16.13)')' Ortho , variance, erreur, actstep',  &
      matdorth(vecs(1,1,is), kdim, ndims(is)), ABS(norm), ABS(ERR_c),actstep
  IF(ABS(norm) < precis) GO TO 99
END DO !iter
99   CONTINUE

RETURN
END SUBROUTINE utgradstepc
#endif

#ifdef REALSWITCH
SUBROUTINE utgradstepr(is,iprint,q0,iter1)

!c     Nonlinear gradient iteration to optmially localized states:
!c      'vecs'    system of eigen-vectors to be determined
!c      'is'      isospin
!c      'iprint'  print level: <0 --> no print at all
!c                             0  --> only final result
!c                             >0 --> print modulus
!c     The matrices and dimension are communicated via
!c     common /radmatrix/.

!#INCLUDE "all.inc"
USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
!INCLUDE 'twost.inc'
!INCLUDE 'radmatrixr.inc'   ! defines also 'KDIM'
!INCLUDE 'vec.inc'!MV
REAL(DP) :: q0(kdfull2,kstate)
REAL(DP) :: dab(kdim,kdim),expdab(kdim,kdim)
REAL(DP) :: dabsto(kdim,kdim)
REAL(DP) :: vecovlpr   ! function names
REAL(DP) :: vecnormr   ! function names
REAL(DP) :: rmathdorth   ! function names
REAL(DP) :: rmatnorme   ! function names

INTEGER :: itmax2,ni
LOGICAL,PARAMETER :: ttest=.false.
LOGICAL,PARAMETER :: tconv=.true.       ! to print convergence

REAL(DP) :: norm,ERR_r

!REAL(DP) :: variance,variance2  ! variance of step  ??
!REAL(DP) :: radmax              ! max. squared radius
REAL(DP) :: actstep,stepnew,enold_2st,actprecis   !,radvary
!REAL(DP) :: varstate(kdim),averstate(kdim)
REAL(DP) :: enstore(3),stepstore(3)        ! storage for optimized step
REAL(DP) :: dampopt=0.7D0,steplim=1.2D0    ! optimized stepsize
LOGICAL,PARAMETER :: topt=.false.       ! switch to optimized step


!-------------------------------------------------------

itmax2=symutbegin   ! 500    !Supprimé symutbegin

IF(ttest) THEN
  write(6,*) 'entree utgradstepr'!MV
  write (6,'(4f12.5)') &
    ((vecsr(ii,jj,is), ii=1,ndims(is)),jj=1,ndims(is))!MV
END IF

ni=ndims(is)

! update unitary transformation 'vecs' with previous exponentional
CALL rmatmult (vecsr(1,1, is),rexpdabold(1,1,is),dabsto,kdim,ni)
CALL rmatcopy(dabsto,vecsr(1,1,is),kdim,ni)

actstep = step/radmaxsym  ! radmaxsym obsolete, set to 1D0
actprecis = max(precis,1D-1*sumvar2)
!actprecis = precis
!WRITE(*,*) ' precis,actprecis,sumvar2=',precis,actprecis,sumvar2
IF(tconv) THEN
  OPEN(353,file='2st-stat-conv.res',POSITION='append')
  WRITE(353,*) '# convergence symmetry condition. is=',is
  WRITE(353,'(a)') '# Ortho , variance, erreur , actstep, Spin'
  WRITE(353,'(a,i4,1pg13.5)') &
    ' Iter,Ortho,variance,erreur,energy. Spin,precis=',is,actprecis
END IF
!IF(tconv .AND. ttest) THEN
IF(tconv) THEN
  write(353,*) 'entree utgradstepr. Spin=',is  !MV
  write (353,'(4f12.5)') &
    ((vecsr(ii,jj,is), ii=1,ndims(is)),jj=1,ndims(is))!MV
  CALL FLUSH(353)
END IF
WRITE(6,'(a,i4,1pg13.5)') &
 ' Iter,Ortho,variance,erreur,energy. Spin,precis=',is,actprecis

enold_2st=0D0
DO iter=1,itmax2

  CALL dalphabetar(is, dab, q0) !new DAB
  IF(topt) THEN
    IF(iter.LE.3) THEN
      enstore(iter)=ener_2st(is)
      IF(iter==1) THEN
         stepstore(1)=0D0
      ELSE
         stepstore(iter) = stepstore(iter-1)+actstep
      END IF
    ELSE
      enstore(1) = enstore(2)
      enstore(2) = enstore(3)
      enstore(3) = ener_2st(is)
      stepstore(1) = stepstore(2)
      stepstore(2) = stepstore(3)
      stepstore(3) = stepstore(3)+actstep
      e1der= (enstore(3)-enstore(2))/(stepstore(3)-stepstore(2))
      e1old= (enstore(2)-enstore(1))/(stepstore(2)-stepstore(1))
      e2der= (e1der-e1old)*2D0/(stepstore(3)-stepstore(1))
      stepnew = -dampopt*e1der/e2der
!WRITE(*,*) ' stepnew,actstep=',stepnew,actstep
      stepnew = stepnew/actstep
      IF(stepnew > steplim) stepnew=steplim
      IF(stepnew < 1D0/steplim) stepnew=1D0/steplim
      actstep = stepnew*actstep
!WRITE(*,*) ' stepnew,steplim=',stepnew,steplim,actstep
WRITE(*,*) '  enstore=',enstore
WRITE(*,*) 'stepstore=',stepstore
WRITE(6,*) 'e1der,e1old,e2der,actstep=',e1der,e1old,e2der,actstep
!WRITE(*,*) 'actstep=',actstep
    END IF    
  END IF
  norm=rmatnorme(dab,kdim,ni)
!  actstep = step/norm
  CALL rmatconst(dab,dab,-actstep, kdim,ni)  !MV mutiply DAB by eta
  CALL rmatexp(dab,expdab,kdim,ni) !MV exp in ExpDab
  CALL rmatabtoa (rexpdabold(1,1,is), expdab,kdim,ni)!MV update exp
  CALL rmatabtoa (vecsr(1,1,is), expdab,kdim,ni)!MV update Vecs
  CALL rmatcnjg(rexpdabold(1,1,is), dab,kdim,ni)
  CALL rmatmult(rexpdabold(1,1,is), dab,dabsto,kdim,ni)
  ERR_r=rmatnorme(dabsto,kdim,ni)
  IF(tconv) THEN
       WRITE(353,'(i4,6(1pg13.5))')   &
         iter,rmatdorth(vecsr(1,1,is),kdim,ndims(is)),&
         norm, ERR_r,ener_2st(is)-enold_2st,actstep,&
         actstep*norm**2/(ener_2st(is)-enold_2st)
       CALL FLUSH(353)
  ELSE
    WRITE(6,'(i4,5(1pg13.5))')   &
       iter,rmatdorth(vecsr(1,1,is),kdim,ndims(is)),&
       norm, ERR_r,ener_2st(is)-enold_2st,actstep
    CALL FLUSH(6)
  END IF
  IF(iter.GE.1) enold_2st=ener_2st(is)

  IF(iter>5 .AND. ABS(norm) < actprecis) GO TO 99
  
END DO
99   CONTINUE

IF(tconv) THEN
  WRITE(353,'(1x/1x)') 
  CLOSE(353)
END IF
! write(6,*) 'sortie de utgradstepr'!MV
! write (6,'(4f12.5)')
!     & ((vecsr(ii,jj,is), ii=1,ndims(is)),jj=1,ndims(is))!MV

IF(iprint >= 0) THEN
  WRITE(6,'(2(a,i4),a,1pg12.4)') 'UT cv is reached at it nr.',iter,  &
      ' for spin =',is,' with variance=',norm
END IF
RETURN
END SUBROUTINE utgradstepr
#endif

!-----------DAlphaBetar------------------- !MV!!! symmetry condition in antihermitian matrix

#ifdef REALSWITCH
SUBROUTINE dalphabetar(is,dab,q0)

!#INCLUDE "all.inc"
USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP),INTENT(IN) :: q0(kdfull2,kstate)
REAL(DP),INTENT(OUT) :: dab(kdim, kdim)

LOGICAL,PARAMETER :: ttest=.false.

REAL(DP) :: save1,save2,acc2

!INCLUDE 'radmatrix.inc'   ! defines also 'KDIM'
!INCLUDE 'vec.inc'!MV
REAL(DP),ALLOCATABLE :: usicsp(:),rhosp(:)
REAL(DP),ALLOCATABLE :: uqsym(:,:),symcond(:,:)
REAL(DP),ALLOCATABLE :: utcond(:,:),qsym(:,:)
REAL(DP) :: rwfovlp ! function declaration

!-------------------------------------------------------

ALLOCATE(usicsp(2*kdfull2),rhosp(2*kdfull2))
ALLOCATE(uqsym(kdfull2,kstate),symcond(kstate,kstate))
ALLOCATE(utcond(kdim,kdim),qsym(kdfull2,kstate) )

ishift = (is-1)*nxyz

ener_2st(is)=0D0
DO nb=1,nstate
  IF(ispin(nrel2abs(nb)) == is)THEN
    nbeff = nb - (is-1)*ndims(1)
    
    CALL superpose_stater(qsym(1,nb),vecsr(1,nbeff,is),q0,is)
    save1=enrear      !!! to deactivate cumulation of enrear, enerpw
    save2=enerpw      !!! (else wrong total energy)
    CALL calc_sicspr(rhosp,usicsp,qsym(1,nb),nb)
!WRITE(*,*) 'is,nb,encoulsp,enerpw=',is,nb,encoulsp,enerpw
    ener_2st(is)=ener_2st(is)+(encoulsp+enerpw)*occup(nb)
    enrear=save1
    enerpw=save2
    
    DO ind=1,nxyz
      uqsym(ind,nb) = usicsp(ind+ishift)*qsym(ind,nb)
    END DO
  END IF
END DO !nb

! variational result of the UT constraint

DO nb=1,nstate
  IF(ispin(nrel2abs(nb)) == is)THEN
    nbeff = nb - (is-1)*ndims(1)
    DO na=1,nstate
      IF(ispin(nrel2abs(na)) == is)THEN
        naa = na - (is-1)*ndims(1)
        utcond(naa,nbeff) = -rwfovlp(qsym(1,na),uqsym(1,nb))
      END IF !ispin
    END DO !na
  END IF
END DO !nb
!      call rMatPrint ('utCond', utcond, kstate, ndims(is))
DO nb=1,nstate
  IF(ispin(nrel2abs(nb)) == is)THEN
    nbeff = nb - (is-1)*ndims(1)
    DO na=1,nstate
      IF(ispin(nrel2abs(na)) == is)THEN
        naa = na - (is-1)*ndims(1)
        dab(naa,nbeff)=utcond(naa,nbeff)-utcond(nbeff,naa)
      END IF
    END DO
  END IF
END DO

IF(ttest) call rMatPrint ('DAB', DAB, kstate, ndims(is))

DEALLOCATE(usicsp,rhosp)
DEALLOCATE(uqsym,symcond)
DEALLOCATE(utcond,qsym)


RETURN
END SUBROUTINE dalphabetar
#endif


!-----------DAlphaBeta------------------- !MV!!! symmetry condition in antihermitian matrix

#ifdef COMPLEXSWITCH
SUBROUTINE dalphabeta(is,dab,q0)

!#INCLUDE "all.inc"
USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP),INTENT(IN) :: q0(kdfull2,kstate)
COMPLEX(DP),INTENT(IN OUT) :: dab(kdim, kdim) !  DAB

LOGICAL,PARAMETER :: ttest=.false.


REAL(DP) :: save1,save2,acc2
COMPLEX(DP) :: wfovlp ! function declaration 

!INCLUDE 'radmatrix.inc'   ! defines also 'KDIM'
!INCLUDE 'vec.inc'!MV
REAL(DP),ALLOCATABLE :: usicsp(:),rhosp(:)
COMPLEX(DP),ALLOCATABLE :: qsym(:,:),utcond(:,:)
COMPLEX(DP),ALLOCATABLE :: uqsym(:,:),symcond(:,:)

!-------------------------------------------------------

ALLOCATE(usicsp(2*kdfull2),rhosp(2*kdfull2))
ALLOCATE(qsym(kdfull2,kstate),utcond(kdim,kdim))
ALLOCATE(uqsym(kdfull2,kstate),symcond(kstate,kstate))

ishift = (is-1)*nxyz

DO nb=1,nstate
  IF(ispin(nrel2abs(nb)) == is)THEN
    nbeff = nb - (is-1)*ndims(1)
    
    CALL superpose_state(qsym(1,nb),vecs(1,nbeff,is),q0,is)
    save1=enrear      !!! to deactivate cumulation of enrear, enerpw
    save2=enerpw      !!! (else wrong total energy)
    CALL calc_sicspc(rhosp,usicsp,qsym(1,nb),nb)
    enrear=save1
    enerpw=save2
    
    DO ind=1,nxyz
      uqsym(ind,nb) = usicsp(ind+ishift)*qsym(ind,nb)
    END DO! nxyz
  END IF
END DO !nb

! variational result of the UT constraint

DO nb=1,nstate
  IF(ispin(nrel2abs(nb)) == is)THEN
    nbeff = nb - (is-1)*ndims(1)
    DO na=1,nstate
      IF(ispin(nrel2abs(na)) == is)THEN
        naa = na - (is-1)*ndims(1)
        utcond(naa,nbeff) = -wfovlp(qsym(1,na),uqsym(1,nb))
      END IF !ispin
    END DO !na
  END IF
END DO !nb
!      call MatPrint ('utCond', utcond, kstate, ndims(is))
DO nb=1,nstate
  IF(ispin(nrel2abs(nb)) == is)THEN
    nbeff = nb - (is-1)*ndims(1)
    DO na=1,nstate
      IF(ispin(nrel2abs(na)) == is)THEN
        naa = na - (is-1)*ndims(1)
        dab(naa,nbeff)=utcond(naa,nbeff)-CONJG (utcond(nbeff,naa))
      END IF
    END DO
  END IF
END DO
!      call MatPrint ('DAB', DAB, kstate, ndims(is))

DEALLOCATE(usicsp,rhosp)
DEALLOCATE(qsym,utcond)
DEALLOCATE(uqsym,symcond)


RETURN
END SUBROUTINE dalphabeta
#endif

  
!-----avermatrix-----------------------------------------
  
#ifdef REALSWITCH
  REAL(DP) FUNCTION avermatrixr(a,vec,ndim,kdim)
#else
  COMPLEX(DP) FUNCTION avermatrix(a,vec,ndim,kdim)
#endif
  USE params, ONLY: DP
  
!      implicit 
  
!     Average of matrix 'a' taken with state 'vec' of
!     actual length 'ndim', dimensioned with 'kdim'.
  
  INTEGER,INTENT(IN) :: ndim,kdim
  INTEGER :: i,j
#ifdef REALSWITCH
  REAL(DP),INTENT(IN) :: a(kdim,kdim),vec(kdim)
  REAL(DP) :: aver
#else
  COMPLEX(DP),INTENT(IN) :: a(kdim,kdim),vec(kdim)
  COMPLEX(DP) :: aver
#endif
  
!-------------------------------------------------------
  
#ifdef REALSWITCH
  aver = 0D0
#else
  aver = CMPLX(0D0,0D0,DP)
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
REAL(DP),INTENT(IN) :: a(kdim,kdim),vec1(kdim),vec2(kdim)
REAL(DP) :: ovlp
#else
COMPLEX(DP),INTENT(IN) :: a(kdim,kdim),vec1(kdim),vec2(kdim)
COMPLEX(DP) :: ovlp

#endif

!-------------------------------------------------------

#ifdef REALSWITCH
ovlp = 0D0
#else
ovlp = CMPLX(0D0,0D0,DP)
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

INTEGER,INTENT(IN) :: ndim,kdim
INTEGER :: i,j,n
REAL(DP) :: acc2
#ifdef REALSWITCH
REAL(DP) :: acc
REAL(DP),INTENT(IN OUT) :: vecs(kdim,kdim)
!REAL(DP) :: vecovlpr
!REAL(DP) :: vecnormr
#else
COMPLEX(DP) :: acc
COMPLEX(DP),INTENT(IN OUT) :: vecs(kdim,kdim)
!COMPLEX(DP) :: vecovlp
!REAL(DP) :: vecnorm
#endif

!-------------------------------------------------------

DO i=1,ndim
  IF(i > 1) THEN
    DO j=1,i-1
#ifdef REALSWITCH
      acc = vecovlpr(vecs(1,j),vecs(1,i),ndim)
#else
      acc = vecovlp(vecs(1,j),vecs(1,i),ndim)
#endif
      DO n=1,ndim
        vecs(n,i) = vecs(n,i)-acc*vecs(n,j)
      END DO
    END DO
  END IF
#ifdef REALSWITCH
  acc2 = 1D0/SQRT(vecnormr(vecs(1,i),ndim))
#else
  acc2 = 1D0/SQRT(vecnorm(vecs(1,i),ndim))
#endif
  DO n=1,ndim
    vecs(n,i) = vecs(n,i)*acc2
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
REAL(DP),INTENT(IN) :: vec(ndim)
#else
COMPLEX(DP),INTENT(IN) :: vec(ndim)
#endif

REAL(DP) :: acc

!-------------------------------------------------------

acc = 0D0
DO n=1,ndim
#ifdef REALSWITCH
  acc = acc + vec(n)**2
#else
  acc = acc + ABS(vec(n))**2
#endif
END DO
#ifdef REALSWITCH
vecnormr = acc
#else
vecnorm = acc
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
REAL(DP),INTENT(IN) :: vec1(ndim),vec2(ndim)
REAL(DP) :: acc
#else
COMPLEX(DP),INTENT(IN) :: vec1(ndim),vec2(ndim)
COMPLEX(DP) :: acc
#endif

!-------------------------------------------------------

#ifdef REALSWITCH
acc = 0D0
#else
acc = CMPLX(0D0,0D0,DP)
#endif
DO n=1,ndim
#ifdef REALSWITCH
  acc = acc + vec1(n)*vec2(n)
#else
  acc = acc + CONJG(vec1(n))*vec2(n)
#endif
END DO
#ifdef REALSWITCH
vecovlpr = acc
#else
vecovlp = acc
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
REAL(DP),INTENT(IN) :: q0(kdfull2,kstate)
REAL(DP),INTENT(IN OUT) :: wfsup(kdfull2)
REAL(DP),INTENT(IN) :: coeff(kstate)
#else
!INCLUDE 'radmatrix.inc'
COMPLEX(DP),INTENT(IN) :: q0(kdfull2,kstate)
COMPLEX(DP),INTENT(IN OUT) :: wfsup(kdfull2)
COMPLEX(DP),INTENT(IN) :: coeff(kstate)
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
!     The resuls is stored in common/radmatrix/ for further
!     use in localization transformation.

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)


INTEGER :: iunit,PRINT    ! print =1 : printing of the total variance
PARAMETER (iunit=0)    ! set zero to disable testprint

REAL(DP) :: var
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
      wfmom = CMPLX(0D0,0D0,DP)
      xmom = CMPLX(0D0,0D0,DP)
      ymom = CMPLX(0D0,0D0,DP)
      zmom = CMPLX(0D0,0D0,DP)
      xxmom = CMPLX(0D0,0D0,DP)
      yymom = CMPLX(0D0,0D0,DP)
      zzmom = CMPLX(0D0,0D0,DP)
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
    var=0D0
#else
    var=CMPLX(0D0,0D0,DP)
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
!-----init_fsic------------------------------------------------

SUBROUTINE init_fsicr()
STOP ' code not compiled for GSlat or double-set SIC'
RETURN
END SUBROUTINE init_fsicr
#else
SUBROUTINE ccc  ! a dummy subroutine so that the compilation does not abort IF twostsic=0
RETURN
END SUBROUTINE ccc 
#endif

#endif

