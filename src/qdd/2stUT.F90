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


#if(twostsic)

#ifdef REALSWITCH
MODULE twostr
USE params
USE kinetic
USE localize_rad
USE orthmat
IMPLICIT NONE


INTERFACE spmomsmatrix
! automated choice between real and complex versions of spmomsmatrix,
! with spmomsmatrix_r, spmomsmatrix_c both defined at first compilation of twostr.
  MODULE PROCEDURE spmomsmatrix_r, spmomsmatrix_c
END INTERFACE spmomsmatrix

INTERFACE dalphabeta
! same as before, for the real, complex and real+complex versions of dalphabeta.
  MODULE PROCEDURE dalphabeta_r ,dalphabeta_rc, dalphabeta_c
END INTERFACE dalphabeta

LOGICAL,PARAMETER :: tconv=.true.       ! to print convergence

#if(cmplxsic)
COMPLEX(DP),PRIVATE,ALLOCATABLE :: qnewr(:,:)
COMPLEX(DP),PRIVATE,ALLOCATABLE :: psirut(:,:)
#else
REAL(DP),PRIVATE,ALLOCATABLE :: qnewr(:,:)
REAL(DP),PRIVATE,ALLOCATABLE :: psirut(:,:)
#endif
REAL(DP),ALLOCATABLE :: usicall(:,:)
!COMMON /twost/ qnewr,qnew,psirut,psiut

REAL(DP) :: step=0.1D0,precis=1D-6,precisfact=1D-3,phiini=0D0
REAL(DP) :: ener_2st(2)
INTEGER :: symutbegin=100,itut
!COMMON /twostut/ step,precis,symutbegin,radmaxsym,itut
REAL(DP) :: dampopt=0.5D0,steplow=0.01D0,steplim=1.2D0    ! optimized stepsize
LOGICAL :: toptsicstep=.true.       ! switch to optimized step

INTEGER,PRIVATE :: kdim
!     matrices of radial moments

#if(cmplxsic)
COMPLEX(DP),ALLOCATABLE :: rExpDABold(:,:,:)
COMPLEX(DP),ALLOCATABLE,SAVE :: vecsr(:,:,:)    ! searched eigenvectors
#else
REAL(DP),ALLOCATABLE :: rExpDABold(:,:,:)!MV added
REAL(DP),ALLOCATABLE,SAVE :: vecsr(:,:,:)    ! searched eigenvectors
#endif

#endif

#ifdef COMPLEXSWITCH
MODULE twost
USE params
USE kinetic
USE twostr, ONLY:step,precis,symutbegin,precisfact,itut,vecsr,tconv,ener_2st,toptsicstep,steplim,dampopt
USE orthmat
IMPLICIT NONE

COMPLEX(DP),ALLOCATABLE :: qnewut(:,:)
COMPLEX(DP),ALLOCATABLE :: psiut(:,:)

!INTEGER,PRIVATE :: kdim=kstate
INTEGER,PRIVATE :: kdim

!     matrices of radial moments

COMPLEX(DP),ALLOCATABLE :: ExpDABold(:,:,:),wfrotate(:,:,:)
COMPLEX(DP),ALLOCATABLE,SAVE :: vecs(:,:,:)    ! searched eigenvectors

! This parameter 'tnearest' activates the computation of the
! rotation amongst occupuied states. The result ist stored on
! 'wfrotate' and so transferred to the SIC routines.
LOGICAL,PARAMETER :: tnearest=.true.

! The switch 'texpo' activates the "extrapolation" of the unitary
! transformation 'vecs' usign the array 'expdabold'. This array is
! initialzed as unit matrix and propagated in parallel to 'vecs'.
! This strategy seems to help. But it is not so obvious why.
LOGICAL,PARAMETER :: texpo=.true.       

#endif

CONTAINS

#ifdef REALSWITCH
!-----init_fsicr------------------------------------------------

SUBROUTINE init_fsicr()

!     initializes fields for FSIC etc

USE params
USE kinetic
IMPLICIT NONE

INTEGER :: i, j, is, na
!#if(cmplxsic)
INTEGER :: ifcmplxin, ni, nim, nip        ! oly used in complex SIC
!#endif
COMPLEX(DP) :: ccr,csi

!INCLUDE "twost.inc"
!INCLUDE 'radmatrixr.inc'
!
! Namelist FSIC contains numerical parameters for the solution
!               of the symmetry condition:
!  step        = step size (actually 'step/radmaxysm' is the step size)
!  precis      = final limit on precision,
!                early stages use actual variance*1D-4 as limit
!  symutbegin  = number of iterations in symmetry condition
!  toptsicstep = switch to tuning of stepsize by quadratic extrapolation 
!                of energy to a maximum using the last three steps;
!                this option is presently only relevant for the static case;
!                the following parameters are relevant only
!                if this option is set to .true.
!  precisfact  = factor applied to variance to compute actual precision limit
!  dampopt     = damping factor on the optimal step (default 0.7)
!  steplow     = minimum step size (tentative default 0.01)
!  steplim     = maximum step size (tentative default 1.2)
!
NAMELIST /fsic/step,precis,symutbegin,precisfact, &
               dampopt,steplow,steplim,toptsicstep,phiini     !!! UT parameters

!----------------------------------------------------------------

IF(ifsicp < 6) RETURN

!     read SIC  specific parameters

OPEN(5,STATUS='old',FORM='formatted',FILE='for005.'//outnam)
READ(5,fsic)
WRITE(6,'(a,4(1pg13.5))')  &
    ' SIC running with step,precis,SymUtBegin=',  &
    step,precis,symutbegin
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
#if(symmcond)
ALLOCATE(usicall(kdfull2,kstate))
#endif
ALLOCATE(qnewr(kdfull2,kstate))
ALLOCATE(psirut(kdfull2,kstate))

kdim=kstate
ALLOCATE(rExpDABold(kstate, kstate, 2))!MV added
ALLOCATE(vecsr(kstate,kstate,2))   ! searched eigenvectors

!     initialize unitary matrix

IF(istat==0) THEN
  DO is=1,2
    DO i=1,ndims(is)
      DO j=1,ndims(is)
        IF(i == j) THEN
!#if(cmplxsic)
!          vecsr(i,j,is) = CMPLX(0D0,1D0)
          vecsr(i,j,is) = 1D0
!#else
!          vecsr(i,j,is) = 1D0
!#endif
        ELSE
          vecsr(i,j,is) = 0D0
        END IF
      END DO
    END DO
#if(cmplxsic)
    ifcmplxin=2
    IF(ifcmplxin==1) THEN
!      phiini=PI*1D-3
      ni=ndims(is)
      nim=ni-1      
      vecsr(ni,ni,is)=CMPLX(COS(phiini),0D0,DP)
      vecsr(nim,nim,is)=vecsr(ni,ni,is)
      vecsr(ni,nim,is)=CMPLX(0D0,SIN(phiini),DP)
      vecsr(nim,ni,is)=CMPLX(0D0,-SIN(phiini),DP)
    ELSE IF(ifcmplxin==2) THEN
      csi = CMPLX(0D0,SIN(phiini)/SQRT(2D0),DP)
      ccr = CMPLX(SQRT(1D0-2D0*ABS(csi)**2),0D0,DP)
      DO ni=1,ndims(is)
        vecsr(ni,ni,is) = ccr
        nip=MOD(ni,ndims(is))+1
        vecsr(ni,nip,is) = csi
        nim=ni-1
        if(nim==0) nim=ndims(is)
        vecsr(ni,nim,is) = -csi
      END DO
      CALL orthnorm(vecsr(:,:,is),ndims(is),kdim)
    END IF
#endif    
  END DO
END IF

! initialize protocol file
IF(tconv) THEN
OPEN(353,file='2st-stat-conv-1.res')
REWIND(353)
WRITE(353,'(a)') '# protocol of convergence of symmetry condition for ispin=1'
CLOSE(353)
OPEN(353,file='2st-stat-conv-2.res')
REWIND(353)
WRITE(353,'(a)') '# protocol of convergence of symmetry condition for ispin=2'
CLOSE(353)
END IF


RETURN
END SUBROUTINE init_fsicr

!-----end_fsicr------------------------------------------------
SUBROUTINE end_fsicr()

USE params
USE kinetic
IMPLICIT NONE

!   frees workspace for static SIC

DEALLOCATE(rExpDABold)!MV added
DEALLOCATE(vecsr)   ! searched eigenvectors


RETURN
END SUBROUTINE end_fsicr

#endif



#ifdef COMPLEXSWITCH
!-----init_fsic------------------------------------------------

SUBROUTINE init_fsic()

!     initializes fields for FSIC etc

USE params
USE kinetic
IMPLICIT NONE

NAMELIST /fsic/step,precis,symutbegin   !!! UT parameters

!----------------------------------------------------------------

kdim=kstate
ALLOCATE(ExpDABold(kstate,kstate,2),wfrotate(kstate,kstate,2))
ALLOCATE(vecs(kstate,kstate,2))    ! searched eigenvectors

ALLOCATE(qnewut(kdfull2,kstate))
ALLOCATE(psiut(kdfull2,kstate))

IF(tconv) THEN
OPEN(353,file='2st-stat-conv-1.res')
REWIND(353)
WRITE(353,'(a)') '# protocol of convergence of symmetry condition for ispin=1'
CLOSE(353)
OPEN(353,file='2st-stat-conv-2.res')
REWIND(353)
WRITE(353,'(a)') '# protocol of convergence of symmetry condition for ispin=2'
CLOSE(353)
END IF

RETURN
END SUBROUTINE init_fsic

!-----ExpDabvol_rotatewf_init----------------------------------

SUBROUTINE expdabvol_rotate_init
IMPLICIT NONE
INTEGER::is
  DO is=1,2   !MV initialise ExpDabOld                              
     CALL MatUnite(ExpDabOld(:,:,is), kstate,ndims(is))
     CALL MatUnite(wfrotate(:,:,is), kstate,ndims(is))
  END DO
END SUBROUTINE expdabvol_rotate_init

!-----end_fsic------------------------------------------------

SUBROUTINE end_fsic()

!     terminates fields for FSIC etc

USE params
USE kinetic
IMPLICIT NONE

NAMELIST /fsic/step,precis,symutbegin   !!! UT parameters

!----------------------------------------------------------------

kdim=kstate
DEALLOCATE(ExpDABold,wfrotate,vecs)
DEALLOCATE(qnewut,psiut)

IF(tconv) CLOSE(353)

RETURN
END SUBROUTINE end_fsic
#endif

!     ******************************


#ifdef COMPLEXSWITCH
!-----init_vecs-------------------------------------------------

SUBROUTINE init_vecs()

!     Initializes complex superposition coefficients for SIC.

USE params
USE kinetic
IMPLICIT NONE 

INTEGER :: i, ii, j, jj
!----------------------------------------------------------------

WRITE(6,*) 'vecsr initvecs'!MV
WRITE (6,'(3f12.6)') ((vecsr(ii,jj,1), ii=1,3),jj=1,3)!MV

DO i=1,kstate
  DO j=1,kstate
!#if(cmplxsic)
    vecs(i,j,1) = vecsr(i,j,1)
    vecs(i,j,2) = vecsr(i,j,2)
!#else
!    vecs(i,j,1) = CMPLX(vecsr(i,j,1),0D0,DP)
!    vecs(i,j,2) = CMPLX(vecsr(i,j,2),0D0,DP)
!#endif
  END DO
END DO

WRITE(6,*) 'vecs fin initvecs'!MV
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
IMPLICIT NONE
REAL(DP) :: rho(2*kdfull2)
REAL(DP) :: aloc(2*kdfull2)
REAL(DP) :: psir(kdfull2,kstate)
INTEGER,INTENT(IN)::iter1

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
ELSE IF(ifsicp == 9) THEN   !    DSIC
  CALL calc_fullsic(psirut,qnewr)
END IF

RETURN
END SUBROUTINE static_sicfield


!-----infor_sic------------------------------------------------

SUBROUTINE infor_sic(psir)

!     Compute and print variances


USE params
USE kinetic
USE util, ONLY:wfovlp
IMPLICIT NONE
REAL(DP),INTENT(IN) :: psir(kdfull2,kstate)

INTEGER :: is, na, nb
#if(cmplxsic)
COMPLEX(DP) :: acc
#else
REAL(DP) :: acc
#endif

!----------------------------------------------------------------

IF(ifsicp < 7) RETURN

!     TESTS

WRITE(6,'(a)') 'DIAGONAL STATES :'
CALL spmomsmatrix(psir,1)       !!! to print the total variance
WRITE(6,'(a)') 'LOCALIZED STATES :'
CALL spmomsmatrix(psirut,1)

IF(ifsicp .GE. 8) THEN   !!! to calculate the total
  DO is=1,2                       !!! violation of symmetry condition
    acc = 0D0
    DO na=1,nstate
      IF(ispin(na) == is)THEN
        DO nb=1,nstate
          IF(ispin(nb) == is)THEN
            acc = ( wfovlp(psirut(:,na),qnewr(:,nb)) -  &
                wfovlp(psirut(:,nb),qnewr(:,na)) )**2 + acc
          END IF
        END DO
      END IF
    END DO
    WRITE(6,'(a,i3,a,2(1pg12.4))')  &
        'For spin',is,'  Total violation of SymCond',SQRT(acc)
  END DO
END IF

RETURN
END SUBROUTINE infor_sic


#if(twostsic)
#if(!cmplxsic)
!-----diag_lagr------------------------------------------------

SUBROUTINE diag_lagr(psir)

!     Diagonalize the matrix of Lagrange parameters and print
!     information from that step.

USE params
USE kinetic
USE localize_rad, ONLY:superpose_state
IMPLICIT NONE

REAL(DP),INTENT(IN OUT) :: psir(kdfull2,kstate)

!     vect1,2 must be declared in each spin subspace for correct
!     diagonalization of SIC with multipliers

INTEGER :: i, j, is, na, nb, nbeff
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
    CALL givens(a,root(:,is),vect1(:,:),ndims(is),ndims(is),ndims(is))
  ELSE IF(is == 2)THEN
    CALL givens(a,root(:,is),vect2(:,:),ndims(is),ndims(is),ndims(is))
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
        CALL superpose_state(psir(:,nb),vect1(1:kstate,nbeff),psirut,is)
      ELSE IF(is == 2)THEN
        CALL superpose_state(psir(:,nb),vect2(1:kstate,nbeff),psirut,is)
      END IF
    END IF
  END DO
  
  WRITE(6,*) 'Diag of the Lagrange matrix for spin subs',is
  DO j=1,ndims(is)
    WRITE(6,*) 'state',j,'single energies',root(j,is)
  END DO
  
END DO

!ifsicp=8

WRITE(6,*) 'DIAGONAL STATES :'
CALL spmomsmatrix(psir,1)  !!! to print the total variance
WRITE(6,*) 'LOCALIZED STATES :'
CALL spmomsmatrix(psirut,1)

RETURN
END SUBROUTINE diag_lagr
#endif
#endif

!-----subtr_sicpot------------------------------------------------

SUBROUTINE subtr_sicpot(q1,nbe)

! subtract the SIC potential acting on the acual 'q1'

USE params
USE kinetic
IMPLICIT NONE
REAL(DP), INTENT(IN OUT) :: q1(kdfull2)
INTEGER, INTENT(IN) :: nbe

INTEGER :: is, na, nb, nae
!----------------------------------------------------------------

IF(ifsicp < 8) RETURN

is=ispin(nbe)
nb = nbe - (is-1)*ndims(1)
DO na=1,ndims(is)
  nae = na + (is-1)*ndims(1)
#if(cmplxsic)
  q1(:) =q1(:)-qnewr(:,nae)*CONJG(vecsr(nb,na,is))
#else
  q1(:) =q1(:)-qnewr(:,nae)*vecsr(nb,na,is)
#endif
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


!     basic arrays and workspace

USE params
USE kinetic
USE util, ONLY:wfovlp
USE localize_rad, ONLY:superpose_state
IMPLICIT NONE

COMPLEX(DP), INTENT(IN OUT) :: q0(kdfull2,kstate)
COMPLEX(DP), INTENT(IN OUT) :: q0ut(kdfull2,kstate)
INTEGER,INTENT(IN) :: iter1

INTEGER :: is, nb, nbeff

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
  CALL superpose_state(q0ut(:,nb),vecs(:,nbeff,is),q0,is)
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
USE localize_rad, ONLY: superpose_state

IMPLICIT NONE

INTEGER :: is, nb, nbeff
!     basic arrays and workspace

#ifdef REALSWITCH
REAL(DP) :: q0(kdfull2,kstate)
#if(cmplxsic)
COMPLEX(DP) :: q0ut(kdfull2,kstate)
#else
REAL(DP) :: q0ut(kdfull2,kstate)
#endif
#else
COMPLEX(DP) :: q0(kdfull2,kstate)
COMPLEX(DP) :: q0ut(kdfull2,kstate)
#endif
INTEGER,INTENT(IN) :: iter1

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
    CALL utgradstepc(is,0,q0,iter1)   !!!!! new vecs (gradient method)
#endif
  END IF
END DO

DO nb=1,nstate
  is = ispin(nrel2abs(nb))
  nbeff = nb - (is-1)*ndims(1)
#ifdef REALSWITCH
  CALL superpose_state(q0ut(:,nb),vecsr(:,nbeff,is),q0,is)
#else
  CALL superpose_state(q0ut(:,nb),vecs(:,nbeff,is),q0,is)
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

!     Nonlinear gradient iteration to optimally localized states:
!      'vecs'    system of eigen-vectors to be determined
!      'is'      isospin
!      'iprint'  print level: <0 --> no print at all
!                             0  --> only final result
!                             >0 --> print modulus

USE params
USE kinetic
USE twostr, ONLY: dalphabeta
IMPLICIT NONE


INTEGER,INTENT(IN) :: is
INTEGER,INTENT(IN) :: iprint
COMPLEX(DP),INTENT(IN) :: q0(kdfull2,kstate)
INTEGER,INTENT(IN) :: iter1

INTEGER :: iter
REAL(DP) :: e1der, e2der, e1old, enold_2st, ERR_c,  stepnew, steplow
COMPLEX(DP) :: xlambda(kdim,kdim)
COMPLEX(DP) :: acc
COMPLEX(DP) :: dab(kdim,kdim),expdab(kdim,kdim)           !MV! workspace
COMPLEX(DP) :: dabsto(kdim,kdim)          !MV! workspace

INTEGER :: itmax2,ni
LOGICAL,PARAMETER :: ttest=.false.

REAL(DP) :: norm                         ! variance of step, erreur du produit
REAL(DP) :: actstep
REAL(DP) :: enstore(3),stepstore(3)      ! storage for optimized step


!-------------------------------------------------------

itmax2=symutbegin 
ni=ndims(is)

! update unitary transformation 'vecs' with previous exponentional
!IF(is==1) WRITE(*,*) ' vecs before:',vecs(1:ni,1:ni,is)
!IF(is==1) WRITE(*,*) ' wfrotate:',wfrotate(1:ni,1:ni,is)

IF(tnearest)  vecs(1:ni,1:ni,is) = MATMUL(wfrotate(1:ni,1:ni,is),vecs(1:ni,1:ni,is))

!IF(is==1) WRITE(*,*) ' vecs after:',vecs(1:ni,1:ni,is)

IF(texpo) THEN
  vecs(1:ni,1:ni,is) = MATMUL(vecs(1:ni,1:ni,is),expdabold(1:ni,1:ni,is))
END IF

actstep = step     ! radmaxsym obsolete, set to 1D0

IF(tconv) THEN
  IF(is==1) OPEN(353,file='2st-stat-conv-1.res',POSITION='append')
  IF(is==2) OPEN(353,file='2st-stat-conv-2.res',POSITION='append')
  WRITE(353,*) '# convergence symmetry condition. timestep,is=',iter1,is
!  WRITE(353,*) 'vecs before:'
!  DO na=1,ni
!    WRITE(353,'(3(2(1pg13.5),2x))') vecs(1:ni,na,is)
!  END DO
  WRITE(353,'(a)') '# iter Ortho , variance, erreur , actstep'
ELSE
  WRITE(6,'(a)')' iter Ortho , variance, non-linearity, actstep'
END IF

!WRITE(*,*) ' before: vecs=',vecs(1:ni,1:ni,is)

enold_2st=0D0
DO iter=1,itmax2
  CALL dalphabeta(is, dab, vecs, q0)     !actual symmetry condition on 'dab'
  IF(toptsicstep) THEN
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
      IF(stepnew < 0D0) stepnew=step
!      stepnew = stepnew/actstep
!      actstep = stepnew*actstep
!      IF(stepnew > steplim) stepnew=steplim
!      IF(stepnew < steplow) stepnew=steplow
      actstep = max(steplow,min(stepnew,steplim))
      IF(ttest) THEN
        WRITE(*,*) ' stepnew,steplim=',stepnew,steplim,actstep
        WRITE(*,'(a,30(1pg15.7))') '    enstore=',enstore
        WRITE(*,'(a,30(1pg15.7))') '  stepstore=',stepstore
        WRITE(6,'(a,3(1pg15.7))') 'e1d,e1o,e2d=',e1der,e1old,e2der
        WRITE(6,'(a,2(1pg15.7))') 'act/newstep=',actstep,stepnew
      END IF
    END IF    
  END IF
  norm = SQRT(SUM(dab(1:ni,1:ni)*CONJG(dab(1:ni,1:ni))))
  dab(1:ni,1:ni) = CMPLX(-actstep,0D0)*dab(1:ni,1:ni)
  CALL matexp(dab,expdab,kdim,ni)  ! exponentiate 'dab' to 'ExpDab'
  vecs(1:ni,1:ni,is) = MATMUL(vecs(1:ni,1:ni,is),expdab(1:ni,1:ni))
  IF(texpo) THEN
    expdabold(1:ni,1:ni,is) = MATMUL(expdabold(1:ni,1:ni,is),expdab(1:ni,1:ni))
    dab(1:ni,1:ni) = CONJG(expdabold(1:ni,1:ni,is))
    dabsto(1:ni,1:ni) = MATMUL(expdabold(1:ni,1:ni,is),dab(1:ni,1:ni))
    ERR_c = SQRT(ABS(SUM(dabsto(1:ni,1:ni)*CONJG(dabsto(1:ni,1:ni))-ndims(is))))
  END IF
  IF(tconv) THEN
    WRITE(353,'(i4,6(1pg13.5))') &
      iter,matdorth(vecs(:,:,is),kdim,ndims(is)),ABS(norm),&
         actstep*norm**2/(ener_2st(is)-enold_2st),actstep,&
         ener_2st(is)-enold_2st
       CALL FLUSH(353)

  ELSE
    WRITE(6,'(a,4f16.13)')' Ortho , variance, erreur, actstep',  &
      matdorth(vecs(:,:,is), kdim, ndims(is)), ABS(norm), ABS(ERR_c),actstep
  END IF
  IF(iter.GE.1) enold_2st=ener_2st(is)
  IF(ABS(norm) < precis) EXIT
END DO !iter


!WRITE(*,*) '  after: vecsr=',vecs(1:ni,1:ni,is)

IF(tconv) close(353)
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

USE params
USE kinetic
IMPLICIT NONE

REAL(DP),INTENT(IN) :: q0(kdfull2,kstate)
#if(cmplxsic)
COMPLEX(DP) :: dab(kdim,kdim),expdab(kdim,kdim)
COMPLEX(DP) :: dabsto(kdim,kdim)
#else
REAL(DP) :: dab(kdim,kdim),expdab(kdim,kdim)
REAL(DP) :: dabsto(kdim,kdim)
#endif
INTEGER,INTENT(IN) :: is,iprint,iter1

INTEGER :: i, iter, itmax2, jj,  ni
LOGICAL,PARAMETER :: ttest=.false.

REAL(DP) :: e1der, e1old, e2der, ERR_r, enorm, norm

!REAL(DP) :: variance,variance2  ! variance of step  ??
!REAL(DP) :: radmax              ! max. squared radius
REAL(DP) :: actstep,stepnew,enold_2st,actprecis   !,radvary
!REAL(DP) :: varstate(kdim),averstate(kdim)
REAL(DP) :: enstore(3),stepstore(3)        ! storage for optimized step

! The optimized step computes a quadratic from for the SIC energy
! from the previous three iterations and extrapolates this form
! to the maximum.
! The step size is limited from above by the parameter 'steplim'
! and from below by 'steplow'.
! If a negative step size emerges, the initial step size is used.
!

!-------------------------------------------------------

! tentative change of step size during iteration
IF(iter1 < 100000) THEN               
   itmax2=symutbegin 
ELSE
   itmax2=1
!   step=epswf
END IF

ni=ndims(is)

IF(ttest) THEN
  write(6,*) 'entree utgradstepr: is=',is
  DO jj=1,ni  
    write (6,'(6(1pg13.5))') vecsr(:,jj,is)
  END DO
END IF


! update unitary transformation 'vecs' with previous exponentional
!dabsto(1:ni,1:ni) = MATMUL(vecsr(1:ni,1:ni,is),rexpdabold(1:ni,1:ni,is))
!vecsr(1:ni,1:ni,is) = dabsto(1:ni,1:ni)
actstep = step  ! radmaxsym obsolete, set to 1D0
actprecis = max(precis,precisfact*sumvar2)
!actprecis = precis
!WRITE(*,*) ' precis,actprecis,sumvar2=',precis,actprecis,sumvar2
IF(tconv) THEN
  IF(is==1) OPEN(353,file='2st-stat-conv-1.res',POSITION='append')
  IF(is==2) OPEN(353,file='2st-stat-conv-2.res',POSITION='append')
  WRITE(353,*) '# convergence symmetry condition. iter1,is=',iter1,is
  WRITE(353,'(a)') '# Ortho , variance, erreur , actstep, Spin'
  IF(ttest) THEN
    WRITE(353,*) ' vecsr before: is=',is
    DO i=1,ni
      WRITE(353,'(20(1pg13.5))') vecsr(1:ni,i,is)
    END DO
    WRITE(353,*) ' rexpdabold before: is=',is
    DO i=1,ni
      WRITE(353,'(20(1pg13.5))') rexpdabold(1:ni,i,is)
    END DO
    WRITE(353,'(a,i4,1pg13.5)') &
      ' Iter,Ortho,variance,erreur,energy. Spin,precis=',is,actprecis
  END IF
  CALL FLUSH(353)
END IF
WRITE(6,'(a,i4,1pg13.5)') &
 ' Iter,Ortho,variance,erreur,energy. Spin,precis=',is,actprecis

!CALL test_symmcond(is,vecsr(1,1,is),q0)

enold_2st=0D0
DO iter=1,itmax2

  CALL dalphabeta(is, dab, vecsr, q0) !new DAB
  IF(iter==1) THEN
    WRITE(353,*) ' dab before: is=',is
    DO i=1,ni
      WRITE(353,'(20(1pg13.5))') dab(1:ni,i)
    END DO
  END IF
  IF(toptsicstep) THEN
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
      IF(stepnew < 0D0) stepnew=step
!      stepnew = stepnew/actstep
!      actstep = stepnew*actstep
!      IF(stepnew > steplim) stepnew=steplim
!      IF(stepnew < steplow) stepnew=steplow
      actstep = max(steplow,min(stepnew,steplim))
!      WRITE(*,'(i52(1pg13.5))') ' is,stepnew,actstep=',is,stepnew,actstep
      IF(ttest) THEN
        WRITE(*,*) ' stepnew,steplim=',stepnew,steplim,actstep
        WRITE(*,'(a,30(1pg15.7))') '    enstore=',enstore
        WRITE(*,'(a,30(1pg15.7))') '  stepstore=',stepstore
        WRITE(6,'(a,3(1pg15.7))') 'e1d,e1o,e2d=',e1der,e1old,e2der
        WRITE(6,'(a,2(1pg15.7))') 'act/newstep=',actstep,stepnew
      END IF
    END IF    
  END IF
#if(cmplxsic)
  norm=SQRT(SUM(dab(1:ni,1:ni)*CONJG(dab(1:ni,1:ni))))  ! rmatnorme(dab,kdim,ni)
#else
  norm=SQRT(SUM(dab(1:ni,1:ni)**2))  ! rmatnorme(dab,kdim,ni)
#endif
  dab(1:ni,1:ni) = -actstep*dab(1:ni,1:ni)
  CALL  matexp(dab,expdab,kdim,ni)            ! MV exp in ExpDab
  rexpdabold(1:ni,1:ni,is) = MATMUL(rexpdabold(1:ni,1:ni,is),expdab(1:ni,1:ni))
  vecsr(1:ni,1:ni,is) = MATMUL(vecsr(1:ni,1:ni,is),expdab(1:ni,1:ni))
#if(cmplxsic)
  dab(1:ni,1:ni) = CONJG(rexpdabold(1:ni,1:ni,is))
#else
  dab(1:ni,1:ni) = rexpdabold(1:ni,1:ni,is)
#endif
  dabsto(1:ni,1:ni) = MATMUL(rexpdabold(1:ni,1:ni,is),dab(1:ni,1:ni))
#if(cmplxsic)
  ERR_r=SQRT(ABS(SUM(dab(1:ni,1:ni)*CONJG(dab(1:ni,1:ni))-ndims(is))))
#else
  ERR_r=SQRT(ABS(SUM(dab(1:ni,1:ni)**2)-ndims(is))) 
#endif
  IF(tconv) THEN
    WRITE(353,'(i4,8(1pg13.5))')   &
     iter,matdorth(vecsr(:,:,is),kdim,ndims(is)),&
     norm, ERR_r,ener_2st(is)-enold_2st,actstep,&
     actstep*norm**2/(ener_2st(is)-enold_2st),ener_2st(is)
    CALL FLUSH(353)
  ELSE
!#if(cmplxsic)
!    WRITE(6,'(i4,5(1pg13.5))')   &
!       iter,matdorth_cmplxsic(vecsr(:,:,is),kdim,ndims(is)),&
!       norm, ERR_r,ener_2st(is)-enold_2st,actstep
!#else
    WRITE(6,'(i4,5(1pg13.5))')   &
       iter,matdorth(vecsr(:,:,is),kdim,ndims(is)),&
       norm, ERR_r,ener_2st(is)-enold_2st,actstep
!#endif
    CALL FLUSH(6)
  END IF
  IF(iter.GE.1) enold_2st=ener_2st(is)

  IF(iter>0 .AND. ABS(norm) < actprecis) EXIT
  
END DO

!CALL test_symmcond(is,vecsr(1,1,is),q0)

IF(tconv) THEN
  IF(ttest) THEN
    WRITE(353,*) '  after: vecsr='
    DO i=1,ni
      WRITE(353,'(20(1pg13.5))') vecsr(1:ni,i,is)
    END DO
    WRITE(353,*) '  after: rexpdabold='
    DO i=1,ni
      WRITE(353,'(20(1pg13.5))') rexpdabold(1:ni,i,is)
    END DO
  END IF
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

#if(!cmplxsic)
SUBROUTINE test_symmcond(is,vecact,q0)

USE params
USE kinetic
IMPLICIT NONE

REAL(DP),INTENT(IN) :: q0(kdfull2,kstate)
REAL(DP),INTENT(IN) :: vecact(kdim, kdim)
INTEGER,INTENT(IN) :: is

LOGICAL,PARAMETER :: ttest=.false.

REAL(DP) :: actstep,enold_2st,norm,enorm
INTEGER :: ni,i

REAL(DP),PARAMETER :: dactstep=0.1D0
INTEGER,PARAMETER :: numteststep=40

REAL(DP),ALLOCATABLE :: vecsav(:,:),dabstep(:,:),dab(:,:),expdab(:,:)

!-------------------------------------------------------

ALLOCATE(vecsav(kdim,kdim),dabstep(kdim,kdim),dab(kdim,kdim),expdab(kdim,kdim))

ni = ndims(is)
vecsav(1:ni,1:ni) = vecsr(1:ni,1:ni,is)
CALL dalphabeta(is, dab, vecsr, q0)
norm=SQRT(SUM(dab(1:ni,1:ni)**2))
enold_2st=ener_2st(is)

WRITE(353,'(a,i2)') ' test linearity of unitary transformation for is=',is
!WRITE(353,'(a,20(1pg13.5))') ' dab:',dab
WRITE(353,'(a)') ' step-size  energy old-energy  energy-ratio '


DO i=0,numteststep
  actstep = MAX(i*dactstep,2D-3)
  dabstep(1:ni,1:ni) = -actstep*dab(1:ni,1:ni)
  CALL matexp(dabstep,expdab,kdim,ni)
  vecsr(1:ni,1:ni,is) = MATMUL(vecsav(1:ni,1:ni),expdab(1:ni,1:ni))
  CALL dalphabeta(is, dabstep, vecsr, q0)
  WRITE(353,'(f8.3,6(1pg13.5))')   &
    actstep,ener_2st(is),ener_2st(is)-enold_2st, &
    actstep*norm**2/(ener_2st(is)-enold_2st),actstep*norm**2
END DO

vecsr(1:ni,1:ni,is) = vecsav(1:ni,1:ni)

DEALLOCATE(vecsav,dabstep,dab,expdab)


END SUBROUTINE test_symmcond
#endif
#endif

#ifdef REALSWITCH
SUBROUTINE dalphabeta_rc(is,dab,vec,q0)

USE params
USE kinetic
USE util, ONLY:wfovlp
USE localize_rad, ONLY:superpose_state
IMPLICIT NONE

INTEGER, INTENT(IN) :: is
COMPLEX(DP),INTENT(OUT) :: dab(kdim, kdim) !  DAB
COMPLEX(DP),INTENT(IN) :: vec(kstate,kstate,2)
REAL(DP),INTENT(IN) :: q0(kdfull2,kstate)

LOGICAL,PARAMETER :: ttest=.false.

INTEGER :: ind, ishift, na, naa, nb, nbeff
REAL(DP) :: save1,save2,acc2

REAL(DP),ALLOCATABLE :: usicsp(:),rhosp(:)
COMPLEX(DP),ALLOCATABLE :: qsym(:,:),utcond(:,:)
COMPLEX(DP),ALLOCATABLE :: uqsym(:,:),symcond(:,:)

!-------------------------------------------------------

ALLOCATE(usicsp(2*kdfull2),rhosp(2*kdfull2))
ALLOCATE(qsym(kdfull2,kstate),utcond(kdim,kdim))
ALLOCATE(uqsym(kdfull2,kstate),symcond(kstate,kstate))

ishift = (is-1)*nxyz

ener_2st(is)=0D0
DO nb=1,nstate
  IF(ispin(nrel2abs(nb)) == is)THEN
    nbeff = nb - (is-1)*ndims(1)
    CALL superpose_state(qsym(:,nb),vec(:,nbeff,is),q0,is)
    save1=enrear      !!! to deactivate cumulation of enrear, enerpw
    save2=enerpw      !!! (else wrong total energy)
!    WRITE(*,*) ' nb,wfnorm=',nb,wfnorm(qsym(1,nb))
    CALL calc_sicsp(rhosp,usicsp,qsym(1,nb),nb)
!    WRITE(*,*) ' nb,energs=',nb,encoulsp,enerpw,occup(nb)
    ener_2st(is)=ener_2st(is)+(encoulsp+enerpw)*occup(nb)
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
        utcond(naa,nbeff) = -wfovlp(qsym(:,na),uqsym(:,nb))
      END IF !ispin
    END DO !na
  END IF
END DO !nb
!WRITE(*,*) ' utcond:',utcond(1:nstate,a:nstate)
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
END SUBROUTINE dalphabeta_rc

!-----------DAlphaBetar------------------- !MV!!! symmetry condition in antihermitian matrix

SUBROUTINE dalphabeta_r(is,dab,vec,q0)

USE params
USE kinetic
USE util, ONLY:wfovlp
USE localize_rad, ONLY:superpose_state
IMPLICIT NONE

INTEGER,INTENT(IN) :: is
REAL(DP),INTENT(OUT) :: dab(kdim, kdim)
REAL(DP),INTENT(IN) :: vec(kstate,kstate,2)
REAL(DP),INTENT(IN) :: q0(kdfull2,kstate)

LOGICAL,PARAMETER :: ttest=.false.
INTEGER :: ind, ishift, na, naa, nb, nbeff
REAL(DP) :: save1,save2,acc2

REAL(DP),ALLOCATABLE :: usicsp(:),rhosp(:)
REAL(DP),ALLOCATABLE :: uqsym(:,:),symcond(:,:)
REAL(DP),ALLOCATABLE :: utcond(:,:),qsym(:,:)

!-------------------------------------------------------

ALLOCATE(usicsp(2*kdfull2),rhosp(2*kdfull2))
ALLOCATE(uqsym(kdfull2,kstate),symcond(kstate,kstate))
ALLOCATE(utcond(kdim,kdim),qsym(kdfull2,kstate) )

ishift = (is-1)*nxyz
dab=0D0
ener_2st(is)=0D0
DO nb=1,nstate
  IF(ispin(nrel2abs(nb)) == is)THEN
    nbeff = nb - (is-1)*ndims(1)
    
    CALL superpose_state(qsym(:,nb),vec(:,nbeff,is),q0,is)
    save1=enrear      !!! to deactivate cumulation of enrear, enerpw
    save2=enerpw      !!! (else wrong total energy)
    CALL calc_sicspr(rhosp,usicsp,qsym(1,nb),nb)
!    CALL printfield(6,rhosp,'rhosp')
!    CALL printfield(6,usicsp,'usicsp')
!WRITE(*,*) 'is,nb,encoulsp,enerpw=',is,nb,encoulsp,enerpw
!WRITE(*,*) 'is,nb,norm-wf=',is,nb,rwfnorm(qsym(1,nb)),&
! SUM(rhosp),SUM(usicsp)

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
        utcond(naa,nbeff) = -wfovlp(qsym(:,na),uqsym(:,nb))
      END IF !ispin
    END DO !na
  END IF
END DO !nb
!    WRITE(6,*) ' utcond: is=',is
!    DO i=1,ndims(is)
!      WRITE(6,'(20(1pg13.5))') utcond(1:ndims(is),i)
!    END DO
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

IF(ttest) THEN
  WRITE(*,*) 'dab:'
  WRITE(*,*) dab(1:ndims(is),1:ndims(is))
END IF

DEALLOCATE(usicsp,rhosp)
DEALLOCATE(uqsym,symcond)
DEALLOCATE(utcond,qsym)


RETURN
END SUBROUTINE dalphabeta_r

!-----------DAlphaBeta------------------- !MV!!! symmetry condition in antihermitian matrix

SUBROUTINE dalphabeta_c(is,dab,vec,q0)

USE params
USE kinetic
USE util, ONLY:wfovlp
USE localize_rad, ONLY:superpose_state
IMPLICIT NONE

INTEGER, INTENT(IN) :: is
COMPLEX(DP),INTENT(OUT) :: dab(kdim, kdim) !  DAB
COMPLEX(DP),INTENT(IN) :: vec(kstate,kstate,2)
COMPLEX(DP),INTENT(IN) :: q0(kdfull2,kstate)

LOGICAL,PARAMETER :: ttest=.false.
INTEGER :: ind, ishift, na, naa, nb, nbeff
REAL(DP) :: save1,save2,acc2

REAL(DP),ALLOCATABLE :: usicsp(:),rhosp(:)
COMPLEX(DP),ALLOCATABLE :: qsym(:,:),utcond(:,:)
COMPLEX(DP),ALLOCATABLE :: uqsym(:,:),symcond(:,:)

!-------------------------------------------------------

ALLOCATE(usicsp(2*kdfull2),rhosp(2*kdfull2))
ALLOCATE(qsym(kdfull2,kstate),utcond(kdim,kdim))
ALLOCATE(uqsym(kdfull2,kstate),symcond(kstate,kstate))

ishift = (is-1)*nxyz

ener_2st(is)=0D0
DO nb=1,nstate
  IF(ispin(nrel2abs(nb)) == is)THEN
    nbeff = nb - (is-1)*ndims(1)
    
    CALL superpose_state(qsym(:,nb),vec(:,nbeff,is),q0,is)
    save1=enrear      !!! to deactivate cumulation of enrear, enerpw
    save2=enerpw      !!! (else wrong total energy)
!    WRITE(*,*) 'DALPHAETA: nb,wfnorm=',nb,wfnorm(qsym(1,nb))
    CALL calc_sicsp(rhosp,usicsp,qsym(1,nb),nb)
!    WRITE(*,*) ' nb,energs=',nb,encoulsp,enerpw,occup(nb)
    ener_2st(is)=ener_2st(is)+(encoulsp+enerpw)*occup(nb)
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
        utcond(naa,nbeff) = -wfovlp(qsym(:,na),uqsym(:,nb))
      END IF !ispin
    END DO !na
  END IF
END DO !nb
!WRITE(*,*) ' utcond:',utcond(1:nstate,a:nstate)
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
END SUBROUTINE dalphabeta_c

!-----spmomsmatrix----------------------------------------------

!     Matrix of spatial moments between single-particle states
!     from real  wf's:
!      wfr    = set of single particle wavefunctions
!     The result is stored in common/radmatrix/ for further
!     use in localization transformation.
!----------------------------------------------------------------------
! REAL version
!----------------------------------------------------------------------
SUBROUTINE spmomsmatrix_r(wfr,PRINT)
USE params
USE kinetic
IMPLICIT NONE


INTEGER,INTENT(IN) :: PRINT    ! print =1 : printing of the total variance
REAL(DP),INTENT(IN) :: wfr(kdfull2,kstate)

INTEGER,PARAMETER :: iunit=0    ! set zero to disable testprint


LOGICAL :: tfirst
INTEGER :: ind, is, ix, iy, iz, ma, mb, na, nac, nb, noff 
REAL(DP) :: var, x1, y1, z1, x2, y2, z2

DATA tfirst/.true./

REAL(DP) :: s,wfmom,xmom,ymom,zmom,xxmom,yymom,zzmom
REAL(DP),ALLOCATABLE :: rrmatr(:,:,:)  ! matrix of r**2
REAL(DP),ALLOCATABLE :: xxmatr(:,:,:)  ! matrix of x**2
REAL(DP),ALLOCATABLE :: yymatr(:,:,:)  ! matrix of y**2
REAL(DP),ALLOCATABLE :: zzmatr(:,:,:)  ! matrix of z**2
REAL(DP),ALLOCATABLE :: xmatr(:,:,:)   ! matrix of x
REAL(DP),ALLOCATABLE :: ymatr(:,:,:)   ! matrix of y
REAL(DP),ALLOCATABLE :: zmatr(:,:,:)   ! matrix of z
!----------------------------------------------------------------------

kdim=kstate
ALLOCATE(rrmatr(kdim,kdim,2))  ! matrix of r**2
ALLOCATE(xxmatr(kdim,kdim,2))  ! matrix of x**2
ALLOCATE(yymatr(kdim,kdim,2))  ! matrix of y**2
ALLOCATE(zzmatr(kdim,kdim,2))  ! matrix of z**2
ALLOCATE(xmatr(kdim,kdim,2))   ! matrix of x
ALLOCATE(ymatr(kdim,kdim,2))   ! matrix of y
ALLOCATE(zmatr(kdim,kdim,2))   ! matrix of z
!
IF(iunit>0) OPEN(iunit,POSITION='append',FILE='pstatmomsmatrix.'//outnam)

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

write(*,'(a,2i5)') 'NDIMS=',ndims

!     compute matrix elements and store

IF(iunit > 0) THEN
  WRITE(iunit,'(a)') 'matrix of s.p. moments:',  &
      '   na   nb      x     y      z     xx     yy     zz'
  tfirst = .false.
END IF
DO na=1,nstate
  
  DO nb=1,na
    IF(ispin(nrel2abs(nb)) == ispin(nrel2abs(na))) THEN
      wfmom = 0D0
      xmom = 0D0
      ymom = 0D0
      zmom = 0D0
      xxmom = 0D0
      yymom = 0D0
      zzmom = 0D0
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
              s=wfr(ind,na)*wfr(ind,nb)
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
      xmatr(mb,ma,is) = xmom
      ymatr(mb,ma,is) = ymom
      zmatr(mb,ma,is) = zmom
      xxmatr(mb,ma,is) = xxmom

      yymatr(mb,ma,is) = yymom
      zzmatr(mb,ma,is) = zzmom
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
    var=0D0
    WRITE(6,'(a)') ' na  x   y   z     xx    yy     zz'
    DO na=1,ndims(is)
      var = var + rrmatr(na,na,is) - xmatr(na,na,is)**2 -  &
          ymatr(na,na,is)**2 - zmatr(na,na,is)**2
      nac = na+(is-1)*noff
      WRITE(6,'(i5,6(1pg13.5))')  &
       na,REAL(xmatr(na,na,is),DP),REAL(ymatr(na,na,is),DP), &
       REAL(zmatr(na,na,is),DP),REAL(xxmatr(na,na,is),DP), &
       REAL(yymatr(na,na,is),DP),REAL(zzmatr(na,na,is),DP)
    END DO
    WRITE(6,*)'For spin =',is,' Total spacial variance ',var
  END DO
END IF

IF(iunit > 0) CLOSE(iunit)
DEALLOCATE(rrmatr)  ! matrix of r**2
DEALLOCATE(xxmatr)  ! matrix of x**2
DEALLOCATE(yymatr)  ! matrix of y**2
DEALLOCATE(zzmatr)  ! matrix of z**2
DEALLOCATE(xmatr)   ! matrix of x
DEALLOCATE(ymatr)   ! matrix of y
DEALLOCATE(zmatr)   ! matrix of z

RETURN
END SUBROUTINE spmomsmatrix_r

!-----------------------------------------------------------------------
! COMPLEX version
!-----------------------------------------------------------------------

SUBROUTINE spmomsmatrix_c(wfr,PRINT)

USE params
USE kinetic
IMPLICIT NONE


INTEGER,INTENT(IN) :: PRINT    ! print =1 : printing of the total variance
INTEGER,PARAMETER :: iunit=0    ! set zero to disable testprint

LOGICAL :: tfirst
INTEGER :: ind, is, ix, iy, iz, ma, mb, na, nac, nb, noff 
REAL(DP) :: var, x1, y1, z1, x2, y2, z2

DATA tfirst/.true./

COMPLEX(DP) :: s,wfmom,xmom,ymom,zmom,xxmom,yymom,zzmom
COMPLEX(DP),INTENT(IN) :: wfr(kdfull2,kstate)
COMPLEX(DP),ALLOCATABLE :: rrmatr(:,:,:)  ! matrix of r**2
COMPLEX(DP),ALLOCATABLE :: xxmatr(:,:,:)  ! matrix of x**2
COMPLEX(DP),ALLOCATABLE :: yymatr(:,:,:)  ! matrix of y**2
COMPLEX(DP),ALLOCATABLE :: zzmatr(:,:,:)  ! matrix of z**2
COMPLEX(DP),ALLOCATABLE :: xmatr(:,:,:)   ! matrix of x
COMPLEX(DP),ALLOCATABLE :: ymatr(:,:,:)   ! matrix of y
COMPLEX(DP),ALLOCATABLE :: zmatr(:,:,:)   ! matrix of z

!----------------------------------------------------------------------

kdim=kstate
ALLOCATE(rrmatr(kdim,kdim,2))  ! matrix of r**2
ALLOCATE(xxmatr(kdim,kdim,2))  ! matrix of x**2
ALLOCATE(yymatr(kdim,kdim,2))  ! matrix of y**2
ALLOCATE(zzmatr(kdim,kdim,2))  ! matrix of z**2
ALLOCATE(xmatr(kdim,kdim,2))   ! matrix of x
ALLOCATE(ymatr(kdim,kdim,2))   ! matrix of y
ALLOCATE(zmatr(kdim,kdim,2))   ! matrix of z
!
IF(iunit>0) OPEN(iunit,POSITION='append',FILE='pstatmomsmatrix.'//outnam)

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

write(*,'(a,2i5)') 'NDIMS=',ndims

!     compute matrix elements and store

IF(iunit > 0) THEN
  WRITE(iunit,'(a)') 'matrix of s.p. moments:',  &
      '   na   nb      x     y      z     xx     yy     zz'
  tfirst = .false.
END IF
DO na=1,nstate
  
  DO nb=1,na
    IF(ispin(nrel2abs(nb)) == ispin(nrel2abs(na))) THEN
      wfmom = CMPLX(0D0,0D0,DP)
      xmom = CMPLX(0D0,0D0,DP)
      ymom = CMPLX(0D0,0D0,DP)
      zmom = CMPLX(0D0,0D0,DP)
      xxmom = CMPLX(0D0,0D0,DP)
      yymom = CMPLX(0D0,0D0,DP)
      zzmom = CMPLX(0D0,0D0,DP)
      
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
              s=CONJG(wfr(ind,na))*wfr(ind,nb)
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
      xmatr(mb,ma,is) = CONJG(xmom)
      ymatr(mb,ma,is) = CONJG(ymom)
      zmatr(mb,ma,is) = CONJG(zmom)
      xxmatr(mb,ma,is) = CONJG(xxmom)
      yymatr(mb,ma,is) = CONJG(yymom)
      zzmatr(mb,ma,is) = CONJG(zzmom)
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
    var=CMPLX(0D0,0D0,DP)
    WRITE(6,'(a)') ' na  x   y   z     xx    yy     zz'
    DO na=1,ndims(is)
      var = var + rrmatr(na,na,is) - xmatr(na,na,is)**2 -  &
          ymatr(na,na,is)**2 - zmatr(na,na,is)**2
      nac = na+(is-1)*noff
      WRITE(6,'(i5,6(1pg13.5))')  &
       na,REAL(xmatr(na,na,is),DP),REAL(ymatr(na,na,is),DP), &
       REAL(zmatr(na,na,is),DP),REAL(xxmatr(na,na,is),DP), &
       REAL(yymatr(na,na,is),DP),REAL(zzmatr(na,na,is),DP)
    END DO
    WRITE(6,*)'For spin =',is,' Total spacial variance ',var
  END DO
END IF

IF(iunit > 0) CLOSE(iunit)
DEALLOCATE(rrmatr)  ! matrix of r**2
DEALLOCATE(xxmatr)  ! matrix of x**2
DEALLOCATE(yymatr)  ! matrix of y**2
DEALLOCATE(zzmatr)  ! matrix of z**2
DEALLOCATE(xmatr)   ! matrix of x
DEALLOCATE(ymatr)   ! matrix of y
DEALLOCATE(zmatr)   ! matrix of z

RETURN
END SUBROUTINE spmomsmatrix_c

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


#  outside module to avoid cyclic referencing across modules

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
IMPLICIT NONE

#ifdef REALSWITCH
#if(cmplxsic)
COMPLEX(DP) :: q0(kdfull2,kstate)
COMPLEX(DP) :: qsic(kdfull2,kstate)
#else
REAL(DP) :: q0(kdfull2,kstate)
REAL(DP) :: qsic(kdfull2,kstate)
#endif
#else
COMPLEX(DP) :: q0(kdfull2,kstate)
COMPLEX(DP) :: qsic(kdfull2,kstate)
#endif

INTEGER :: ind, idx, ishift,  nb
REAL(DP) :: enrearsave, enerpwsave, enrear1, enrear2, enpw1, enpw2, encadd
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
  IF(ifsicp==9) THEN
    CALL calc_sicsp(rhosp,usicsp,q0(1,nb),nb)
  ELSE IF(ifsicp==8) THEN
    CALL calc_sicspr(rhosp,usicsp,q0(1,nb),nb)
  ELSE
    STOP "invalid IFSICP in CALC_FULLSIC"
  END IF
#else
    CALL calc_sicsp(rhosp,usicsp,q0(1,nb),nb)
#endif
    
    IF (ispin(nrel2abs(nb)) == 1) THEN
      enrear1=enrear1+enrear*occup(nb)
      IF(directenergy) THEN
        enpw1=enpw1+enerpw*occup(nb)
      END IF
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
!encadd=encadd/2.0
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

