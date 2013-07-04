
!-----utgradstep-----------------------------------------!PGR old version

#ifdef REALSWITCH
SUBROUTINE utgradstepoldr(is,iprint,q0,iter1)
#else
SUBROUTINE utgradstepold(is,iprint,q0,iter1)
#endif

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
IMPLICIT REAL(DP) (A-H,O-Z)


!INCLUDE 'twost.inc'
#ifdef REALSWITCH
!INCLUDE 'radmatrixr.inc'   ! defines also 'KDIM'
REAL(DP) :: q0(kdfull2,kstate)
REAL(DP) :: xlambda(kdim,kdim)
REAL(DP) :: symm_mat(kdim,kdim)
REAL(DP) :: utcond(kdim,kdim)
REAL(DP) acc
REAL(DP) :: w(kdim)             ! workspace
!REAL(DP) vecovlpr   ! function names
!REAL(DP) vecnormr   ! function names
#else
!INCLUDE 'radmatrix.inc'   ! defines also 'KDIM'
COMPLEX(DP) :: q0(kdfull2,kstate)
COMPLEX(DP) :: xlambda(kdim,kdim)
COMPLEX(DP) :: symm_mat(kdim,kdim)
COMPLEX(DP) :: utcond(kdim,kdim)
COMPLEX(DP) :: acc
COMPLEX(DP) :: w(kdim)             ! workspace
!COMPLEX(DP) :: vecovlp   ! function names
!REAL(DP) vecnorm   ! function names
COMPLEX(DP) :: wfovlp
#endif

INTEGER :: itmax2
LOGICAL,PARAMETER :: ttest=.false.

REAL(DP) :: variance,variance2  ! variance of step
REAL(DP) :: radmax              ! max. squared radius
REAL(DP) :: actstep,radvary
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
      acc=0.D0
      DO n=1,ndims(is)
        acc=acc+vecsr(n,j,is)*utcond(n,i)
#else
        acc=CMPLX(0.D0,0.D0,DP)
        DO n=1,ndims(is)
          acc=acc+CONJG(vecs(n,j,is))*utcond(n,i)
#endif
        END DO
        xlambda(i,j)=acc
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
    
!       application of radii-matrix - substract constraint
    
    variance = 0.D0
    variance2 = 0.D0
    DO i=1,ndims(is)
      DO na=1,ndims(is)
#ifdef REALSWITCH
        acc=0.D0
        DO j=1,ndims(is)
          acc = acc + xlambda(i,j)*vecsr(na,j,is)
#else
          acc=CMPLX(0.D0,0.D0,DP)
          DO j=1,ndims(is)
            acc = acc + xlambda(i,j)*vecs(na,j,is)
#endif
          END DO
          w(na) = utcond(na,i) - acc
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
!       acc=0.D0
!#else
!       acc=dcmplx(0.D0,0.D0)
!#endif
!       do i=1,ndims(is)
!#ifdef REALSWITCH
!       acc=vecsr(na,i,is)*vecsr(nb,i,is)+acc
!#else
!       acc=vecs(na,i,is)*conjg(vecs(nb,i,is))+acc
!#endif
!       enddo
!       write(6,*)na,nb,acc
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
  END SUBROUTINE utgradstepoldr
#else
  END SUBROUTINE utgradstepold
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
  REAL(DP) :: acc
!  REAL(DP) :: avermatrixr          ! function names
#else
  !INCLUDE 'radmatrix.inc'   ! defines also 'KDIM'
  COMPLEX(DP) :: xaver(kdim),yaver(kdim),zaver(kdim)
  COMPLEX(DP) :: utcond(kdim,kdim)
  COMPLEX(DP) :: acc
!  COMPLEX(DP) :: avermatrix          ! function names
#endif
  
  LOGICAL,PARAMETER :: ttest=.false.
  
  REAL(DP) :: radvary             ! variance of radii
  REAL(DP) :: radmax              ! max. squared radius
  REAL(DP) :: rvary(kdim),raver(kdim)
  
!-------------------------------------------------------
  
!       averages of x,y,z
  
  radvary = 0D0
  radmax = 0D0
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
      acc=0D0
#else
      acc=CMPLX(0D0,0D0,DP)
#endif
      DO nb=1,ndims(is)
#ifdef REALSWITCH
        acc = (rrmatr(na,nb,is)  &
        -2.D0 * xaver(i)*xmatr(na,nb,is)    &
              -2.D0 * yaver(i)*ymatr(na,nb,is)  &
              -2.D0 * zaver(i)*zmatr(na,nb,is)) *vecsr(nb,i,is) &
#else
          acc = (rrmatr(na,nb,is)  &
          -2.D0 * xaver(i)*xmatr(na,nb,is)  &
                -2.D0 * yaver(i)*ymatr(na,nb,is)  &
                -2.D0 * zaver(i)*zmatr(na,nb,is)) *vecs(nb,i,is) &
#endif  
            + acc
          END DO
          utcond(na,i)=acc
        END DO
      END DO
      
      RETURN
#ifdef REALSWITCH
    END SUBROUTINE locutr
#else
    END SUBROUTINE locut
#endif



!-----------symut_enc------------------- !MV!!! sym cond calculated with minimization of energy

#ifdef complexswitch
SUBROUTINE symut_enc(is,utcond,q0)

!#INCLUDE "all.inc"
USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

LOGICAL :: ttest
PARAMETER (ttest=.false.)

REAL(DP) :: usicsp(2*kdfull2),rhosp(2*kdfull2)
REAL(DP) :: save1,save2,acc2

!INCLUDE 'radmatrix.inc'   ! defines also 'KDIM'
!INCLUDE 'vec.inc'!MV
COMPLEX(DP) :: q0(kdfull2,kstate), qsym(kdfull2,kstate)
COMPLEX(DP) :: uqsym(kdfull2,kstate),symcond(kstate,kstate)
COMPLEX(DP) :: utcond(kdim,kdim), wfovlp ! function declaration

!-------------------------------------------------------

! write(6,*) 'passe par Sym_utc' !MV

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
    END DO
    
! variationnal result of the UT constraint
    
    DO na=1,nstate
      IF(ispin(nrel2abs(na)) == is)THEN
        naa = na - (is-1)*ndims(1)
        utcond(naa,nbeff) = -wfovlp(q0(1,na),uqsym(1,nb))
      END IF
    END DO
  END IF
END DO
RETURN
END SUBROUTINE symut_enc
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
    REAL(DP) save1,save2,acc2
    
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
        save1=enrear      !!! to deactivate cumulation of enrear, enerpw
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
!       acc2 = 0.D0
!       do na=1,nstate
!       if(ispin(na).eq.is)then
!         do nb=1,nstate
!         if(ispin(nb).eq.is .and. nb.ne.na)then
!#ifdef REALSWITCH
!         acc2 = symcond(na,nb)**2 + acc2
!#else
!         acc2 = abs(symcond(na,nb))**2 + acc2
!#endif
!         endif
!         enddo
!       endif
!       enddo
!       write(6,'(a,(1pg12.4))')
!     & 'Sym UT: evol of SymCond ',sqrt(acc2)
    
    RETURN
#ifdef REALSWITCH
  END SUBROUTINE symut_enr
#else
  END SUBROUTINE symut_en
#endif

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

!NTEGER,PARAMETER :: kdim=kstate

INTEGER,INTENT(IN) :: norder,ndima

#ifdef REALSWITCH
REAL(DP),INTENT(IN) :: symm_mat(kdim,kdim)
REAL(DP),INTENT(INOUT) :: vect(kdim,kdim)
REAL(DP) :: vectwork(kdim,kdim)
REAL(DP) :: symm_acc(kdim,kdim)
REAL(DP) :: accv
REAL(DP) :: w(kdim)             ! workspace
#else
COMPLEX(DP),INTENT(IN) :: symm_mat(kdim,kdim)
COMPLEX(DP),INTENT(INOUT) :: vect(kdim,kdim)
COMPLEX(DP) :: vectwork(kdim,kdim)
COMPLEX(DP) :: symm_acc(kdim,kdim)
COMPLEX(DP) :: accv
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

