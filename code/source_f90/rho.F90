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

#include "define.h"
 
!-----calcrho------------------------------------------------------

#ifdef REALSWITCH
SUBROUTINE calcrhor(rho,q0)
#else
SUBROUTINE calcrho(rho,q0)
#endif

!     density 'rho' for complex or real wavefunctions 'q0'

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
REAL(DP),DIMENSION(:),ALLOCATABLE :: rh
!DIMENSION rh(2*kdfull2)
#endif

#ifdef REALSWITCH
REAL(DP), INTENT(IN) :: q0(kdfull2,kstate)
#else
COMPLEX(DP), INTENT(IN) :: q0(kdfull2,kstate)         ! cPW
#endif
REAL(DP), INTENT(OUT) :: rho(2*kdfull2)

#if(parayes)
!EQUIVALENCE(rh(1),w1(1))
LOGICAL,PARAMETER :: ttestpara=.FALSE.
#endif

!-----------------------------------------------------------------


!     check workspace

#if(parayes)
  ALLOCATE(rh(2*kdfull2))
#endif

!k initialize densities:
#if(parayes)
  rh=0D0
#endif
#if(parano)
  rho=0D0
#endif



DO nb=1,nstate
  ishift = (ispin(nrel2abs(nb))-1)*nxyz ! store spin=2 in upper block
#if(parayes)
  IF(ttestpara) THEN
    WRITE(6,'(a,4i10)') ' RHO: myn,nb,is,ishift=',  &
        myn,nb,ispin(nrel2abs(nb)),ishift
  END IF
#endif
  DO ind=1,nxyz
#ifdef REALSWITCH
    
#if(parayes)
    rh(ind+ishift)=rh(ind+ishift)+ occup(nb)*q0(ind,nb)*q0(ind,nb)
#endif
#if(parano)
    rho(ind+ishift)=rho(ind+ishift)+ occup(nb)*q0(ind,nb)*q0(ind,nb)
#endif
    
#else
    
#if(parayes)
    rh(ind+ishift)=rh(ind+ishift)+ occup(nb)*(CONJG(q0(ind,nb)))*q0(ind,nb)
#endif
#if(parano)
    rho(ind+ishift)=rho(ind+ishift)+ occup(nb)*(CONJG(q0(ind,nb)))*q0(ind,nb)
#endif
    
#endif
  END DO
END DO

!     reorder to total density in lower block (1:nxyz)
!     and difference density in upper block (nxyz+1:2nxyz)

#if(parayes)
CALL mpi_barrier (mpi_comm_world, mpi_ierror)
!         mx=2*nxyz
CALL mpi_allreduce(rh,rho,2*nxyz,mpi_double_precision,  &
    mpi_sum,mpi_comm_world,ic)
CALL mpi_barrier (mpi_comm_world, mpi_ierror)
!         mx=2*nxyz
IF(ttestpara) WRITE(*,*) ' RHO: after allreduce'
DEALLOCATE(rh)
#endif


!GB      sum1=0D0
DO ind=1,nxyz
  IF(numspin==2) THEN
    rhotot      = rho(ind) + rho(ind+nxyz)
    rhodif      = rho(ind) - rho(ind+nxyz)
  ELSE
    rhotot      = rho(ind)
    rhodif      = 0D0
  END IF
  rho(ind)      = rhotot
  rho(ind+nxyz) = rhodif/MAX(rhotot,1D-8)
!GB        sum1 = sum1 + rho(ind)
END DO

CALL emoms(rho)
IF(eproj/=0) CALL projmoms(rho,q0)


#if(gridfft)
IF(istream == 1)  THEN
  CALL stream(rho,q0)
END IF
#endif


RETURN
#ifdef REALSWITCH
END SUBROUTINE calcrhor
#else
END SUBROUTINE calcrho
#endif

#ifdef COMPLEXSWITCH 
!-----calc_current------------------------------------------------------

SUBROUTINE calc_current(current,q0)

!  current 'current' for set of complex wavefunctions 'q0'

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN) :: q0(kdfull2,kstate)
REAL(DP), INTENT(OUT) :: current(kdfull2,3)

COMPLEX(DP), ALLOCATABLE :: dq0(:)

#if(parayes)
STOP "CALC_CURRENT presently not suited for parallel computing"              ! cPW
#endif

!-----------------------------------------------------------------

ALLOCATE(dq0(kdfull2))

! reset 
current=0D0

! accumulate
DO nb=1,nstate
  CALL xgradient_rspace(q0(1,nb),dq0)
  current(:,1) = current(:,1) + occup(nb)*AIMAG(CONJG(q0(:,nb))*dq0(:))
  CALL ygradient_rspace(q0(1,nb),dq0)
  current(:,2) = current(:,2) + occup(nb)*AIMAG(CONJG(q0(:,nb))*dq0(:))
  CALL zgradient_rspace(q0(1,nb),dq0)
  current(:,3) = current(:,3) + occup(nb)*AIMAG(CONJG(q0(:,nb))*dq0(:))
END DO

DEALLOCATE(dq0)

RETURN

END SUBROUTINE calc_current
#endif
!-----spmoms------------------------------------------------------------

#ifdef REALSWITCH
SUBROUTINE spmomsr(wfr,iunit)
#else
SUBROUTINE spmoms(wf,iunit)
#endif

!     spatial moments of single-particle densities from real  wf's:
!     input is
!      wfr    = set of real single particle wavefunctions
!      wf     = set of complex s.p. wavefunctions
!      iunit  = unit number for output

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

#ifdef REALSWITCH
REAL(DP), INTENT(IN)                     :: wfr(kdfull2,kstate)
#else
COMPLEX(DP), INTENT(IN)                  :: wf(kdfull2,kstate)
#endif

LOGICAL, PARAMETER :: ttest=.false.
REAL(DP), ALLOCATABLE :: qeorb(:,:)
INTEGER, INTENT(IN OUT)                  :: iunit
#if(parayes)
INTEGER  :: iprisav(kstate,2)     ! printing communication
#endif

!      logical tfirst
!      data tfirst/.true./

!----------------------------------------------------------------------

ALLOCATE(qeorb(kstate,11))

#if(parayes)
CALL mpi_barrier (mpi_comm_world, mpi_ierror)
CALL  mpi_comm_rank(mpi_comm_world,myn,icode)
IF(myn == 0) THEN
#endif
  WRITE(iunit,'(a)') 'protocol of s.p. moments:',  &
      '  state energy   x   y   z   variance  xx  yy  zz xy xz yz'
#if(parayes)
END IF
#endif

xcmel = 0D0
ycmel = 0D0
zcmel = 0D0
r2el = 0D0

DO nbe=1,nstate
  DO k=1,11
    qeorb(nbe,k)=0D0
  END DO
  
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
          s=wfr(ind,nbe)*wfr(ind,nbe)
#else
          s=wf(ind,nbe)*CONJG(wf(ind,nbe))
#endif
          qeorb(nbe,1)=amoy(nbe)
!                                                       monopole
          qeorb(nbe,2)=qeorb(nbe,2)+s
!                                                       dipole
          qeorb(nbe,3)=qeorb(nbe,3)+s*x1
          qeorb(nbe,4)=qeorb(nbe,4)+s*y1
          qeorb(nbe,5)=qeorb(nbe,5)+s*z1
!                                                       quadrupole
          qeorb(nbe,6)=qeorb(nbe,6)+s*x2
          qeorb(nbe,7)=qeorb(nbe,7)+s*y2
          qeorb(nbe,8)=qeorb(nbe,8)+s*z2

          qeorb(nbe,9)=qeorb(nbe,9)+s*x1*y1
          qeorb(nbe,10)=qeorb(nbe,10)+s*z1*x1
          qeorb(nbe,11)=qeorb(nbe,11)+s*z1*y1
        END IF
      END DO
    END DO
  END DO
 
  DO k=2,11
    qeorb(nbe,k)=qeorb(nbe,k)*dvol
  END DO
!    qeorb(nbe,6)=qeorb(nbe,6)-qeorb(nbe,3)*qeorb(nbe,3)    ?
!    qeorb(nbe,7)=qeorb(nbe,7)-qeorb(nbe,4)*qeorb(nbe,4)    ?
!    qeorb(nbe,8)=qeorb(nbe,8)-qeorb(nbe,5)*qeorb(nbe,5)    ?
!    qeorb(nbe,9)=qeorb(nbe,9)-qeorb(nbe,3)*qeorb(nbe,4)    ?
!    qeorb(nbe,10)=qeorb(nbe,10)-qeorb(nbe,3)*qeorb(nbe,5)  ?
!    qeorb(nbe,11)=qeorb(nbe,11)-qeorb(nbe,4)*qeorb(nbe,5)  ?

#if(parano)
  WRITE(iunit,'(i4,f7.3,4f6.2,2x,6f7.1)')  nbe,qeorb(nbe,1), &
   (qeorb(nbe,j),j=3,5), &
   SQRT(qeorb(nbe,6)+qeorb(nbe,7)+qeorb(nbe,8)),(qeorb(nbe,j),j=9,11)

  qeorb_all(nbe,:)=qeorb(nbe,:)
#endif
END DO

#if(parayes)
DO nbe=1,nstate
  iprisav(nbe,1) = nrel2abs(nbe)
  iprisav(nbe,2) = 3-2*ispin(nrel2abs(nbe))
END DO


IF(myn /= 0) THEN
  nod = myn
  CALL mpi_send(qeorb,11*kstate,mpi_double_precision,  &
      0,nod,mpi_comm_world,ic)
  CALL mpi_send(iprisav,2*kstate,mpi_integer, 0,nod,mpi_comm_world,ic)
  IF(ttest) WRITE(*,*) ' SPMOMS: sent at node:',myn
ELSE
  DO nod2=0,knode-1
    IF(nod2 > 0) THEN
      CALL mpi_recv(qeorb,11*kstate,mpi_double_precision,  &
          nod2,mpi_any_tag,mpi_comm_world,is,ic)
      CALL mpi_recv(iprisav,2*kstate,mpi_integer,  &
          nod2,mpi_any_tag,mpi_comm_world,is,ic)
      IF(ttest) WRITE(*,*)' SPMOMS: recv from  node=',nod2
    END IF
    DO nbe=1,nstate_node(nod2)
      nact = nrel2abs_other(nbe,nod2)
      qeorb_all(nact,:) = qeorb(nbe,:)
      WRITE(iunit,'(i4,f7.3,4f6.2,2x,6f7.1)')  iprisav(nbe,1),qeorb(nbe,1), &
       (qeorb(nbe,j),j=3,5), &
       SQRT(qeorb(nbe,6)+qeorb(nbe,7)+qeorb(nbe,8)),(qeorb(nbe,j),j=9,11)
    END DO
  END DO
END IF
#endif

DEALLOCATE(qeorb)

#if(parayes)
IF(myn == 0) THEN
#endif
  DO nbe=1,nstate_all
    xcmel = xcmel + qeorb_all(nbe,3)
    ycmel = ycmel + qeorb_all(nbe,4)
    zcmel = zcmel + qeorb_all(nbe,5)
    r2el  = r2el + qeorb_all(nbe,6)+qeorb_all(nbe,7)+qeorb_all(nbe,8)
  END DO

  xcmel = xcmel/nstate_all
  ycmel = ycmel/nstate_all
  zcmel = zcmel/nstate_all
  r2el  = SQRT(r2el/nstate_all)

  WRITE(iunit,'(a11,4f6.2)') 'average:   ',xcmel,ycmel,zcmel,r2el
#if(parayes)
END IF
#endif


RETURN
#ifdef REALSWITCH
END SUBROUTINE spmomsr
#else
END SUBROUTINE spmoms
#endif

!-----------------------------------------------------------------

