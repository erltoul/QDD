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
USE util, ONLY:emoms, projmoms
IMPLICIT NONE

#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
REAL(DP),DIMENSION(:),ALLOCATABLE :: rh
#endif

#ifdef REALSWITCH
REAL(DP), INTENT(IN) :: q0(kdfull2,kstate)
#else
COMPLEX(DP), INTENT(IN) :: q0(kdfull2,kstate)         ! cPW
#endif
INTEGER :: ind, ishift, nb
REAL(DP) :: rhodif, rhotot
REAL(DP), INTENT(OUT) :: rho(2*kdfull2)
REAL(DP)::rhoup(kdfull2),rhodown(kdfull2)
REAL(DP) :: rhouparrayfine(2*nx2-1,2*ny2-1,2*nz2-1),rhouparray(nx2,ny2,nz2)
REAL(DP) :: rhodownarrayfine(2*nx2-1,2*ny2-1,2*nz2-1),rhodownarray(nx2,ny2,nz2)
#if(parayes)
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

rhoup=0D0
rhodown=0D0
rhouparrayfine=0D0
rhouparray=0D0
rhodownarrayfine=0D0
rhodownarray=0D0

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
    mpi_sum,mpi_comm_world,icode)
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

CALL emoms(rho)                    ! moments for the whole system (in qe)
IF(eproj/=0) CALL projmoms(rho,q0) ! moments for the projectile and the target, (in qeproj and qetarget)


#if(gridfft)
IF(istream == 1)  THEN
  CALL stream(rho)
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
IMPLICIT NONE

COMPLEX(DP), INTENT(IN) :: q0(kdfull2,kstate)
REAL(DP), INTENT(OUT) :: current(kdfull2,3)

INTEGER :: nb
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

!     spatial moments of single-particle densities from wf's:
!     input is
!      wfr    = set of real single particle wavefunctions
!      wf     = set of complex s.p. wavefunctions
!      iunit  = unit number for output

USE params
IMPLICIT NONE
#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
INTEGER :: nact, nod, nod2
#endif

#ifdef REALSWITCH
REAL(DP), INTENT(IN)                     :: wfr(kdfull2,kstate)
#else
COMPLEX(DP), INTENT(IN)                  :: wf(kdfull2,kstate)
#endif
INTEGER, INTENT(IN)                  :: iunit

LOGICAL, PARAMETER :: ttest=.false.
INTEGER :: ind, ix, iy, iz, j, k, nbe
REAL(DP) :: r2el, s,  xcmel, ycmel, zcmel, x1, y1, z1, x2, y2, z2
REAL(DP), ALLOCATABLE :: qeorb(:,:)

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
      0,nod,mpi_comm_world,icode)
  CALL mpi_send(iprisav,2*kstate,mpi_integer, 0,nod,mpi_comm_world,icode)
  IF(ttest) WRITE(*,*) ' SPMOMS: sent at node:',myn
ELSE
  DO nod2=0,knode-1
    IF(nod2 > 0) THEN
      CALL mpi_recv(qeorb,11*kstate,mpi_double_precision,  &
          nod2,mpi_any_tag,mpi_comm_world,is,icode)
      CALL mpi_recv(iprisav,2*kstate,mpi_integer,  &
          nod2,mpi_any_tag,mpi_comm_world,is,icode)
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


!**************************from 1D to 3D************************
#ifdef REALSWITCH
SUBROUTINE from1Dto3D(va,a,xva,yva,zva)
#else
SUBROUTINE from1Dto3Dc(va,a,xva,yva,zva)
#endif

USE params
IMPLICIT NONE

INTEGER :: xva,yva,zva,ia,ja,ka,i0
#ifdef REALSWITCH
REAL(DP),INTENT (IN) :: va(xva*yva*zva)
REAL(DP),INTENT (OUT) :: a(xva,yva,zva)
#else
COMPLEX(DP),INTENT (IN) :: va(xva*yva*zva)
COMPLEX(DP),INTENT (OUT) :: a(xva,yva,zva)
#endif 

i0=0
DO ka=1,nz2
  DO ja=1,ny2
    DO ia=1,nx2
      i0=i0+1
      a(ia,ja,ka)=va(i0)
    END DO
  END DO
END DO

RETURN
#ifdef REALSWITCH
END SUBROUTINE from1Dto3D
#else
END SUBROUTINE from1Dto3Dc
#endif
!**************************from 3D to 1D************************
#ifdef REALSWITCH
SUBROUTINE from3Dto1D(va,a,xva,yva,zva)
#else
SUBROUTINE from3Dto1Dc(va,a,xva,yva,zva)
#endif

USE params
IMPLICIT NONE

INTEGER,INTENT (IN) :: xva
INTEGER,INTENT (IN) :: yva
INTEGER,INTENT (IN) :: zva
#ifdef REALSWITCH
REAL(DP),INTENT (OUT) :: va(xva*yva*zva)
REAL(DP),INTENT (IN) :: a(xva,yva,zva)
#else
COMPLEX(DP),INTENT (OUT) :: va(xva*yva*zva)
COMPLEX(DP),INTENT (IN) :: a(xva,yva,zva)
#endif 

INTEGER :: i0,ia,ja,ka

i0=0
DO ka=1,zva
  DO ja=1,yva
    DO ia=1,xva
      i0=i0+1
      va(i0)=a(ia,ja,ka)
    END DO
  END DO
END DO

RETURN
#ifdef REALSWITCH
END SUBROUTINE from3Dto1D
#else
END SUBROUTINE from3Dto1Dc
#endif

!----------------------------interpol3D(a)--------------------------------------
#ifdef REALSWITCH
SUBROUTINE interpol3Dn(a1,ao)
#else
SUBROUTINE interpol3Dcn(a1,ao)
#endif

USE params
IMPLICIT NONE
! INTEGER :: nx2,ny2,nz2!,nxyz=nx2*ny2*nz2

#ifdef REALSWITCH
!interpolation target array
! REAL(DP):: va(nxyz)
REAL(DP),INTENT (IN) :: a1(nx2,ny2,nz2)
!interpolation variables
! real(DP),ALLOCATABLE ::ai(:,:,:)
REAL(DP),INTENT (OUT) ::ao(2*nx2,2*ny2,2*nz2)
#else
! COMPLEX(DP):: va(nxyz)
COMPLEX(DP),INTENT (IN) :: a1(nx2,ny2,nz2)
!interpolation variables
! real(DP),ALLOCATABLE ::ai(:,:,:)
COMPLEX(DP),INTENT (OUT) ::ao(2*nx2,2*ny2,2*nz2)
#endif
COMPLEX(dp) :: ai(0:2*nx2,0:2*ny2,0:2*nz2)
COMPLEX(dp) :: a(0:nx2,0:ny2,0:nz2)
INTEGER:: nx2i,ny2i,nz2i
INTEGER :: ia,ja,ka,i0
INTEGER :: ia1,ja1,ka1

ai(:,:,:) = 0.0D0
a=0.0
!write (*,*)'hello1'
nx2i=2*nx2-1
ny2i=2*ny2-1
nz2i=2*nz2-1
! ALLOCATE(ai(nx2i,ny2i,nz2i))

!va is a vector no an array so the first step is to fill the corresponding 3D array a.
!#ifdef REALSWITCH
! call from1Dto3D(va,a,nx2,ny2,nz2) !***********
! #else
! call from1Dto3D(va,a,nx2,ny2,nz2) !***********
! #endif  

DO ia=1,nx2
  DO ja=1,ny2
    DO ka=1,ny2
      a(ia,ja,ka)=a1(ia,ja,ka)
    END DO
  END DO
END DO
!-------------------------------------------------------
DO ia=1,nx2i
  ia1=ia/2+1
  DO ja=1,ny2i
    ja1=ja/2+1
		DO ka=1,nz2i
      ka1=ka/2+1
      ai(ia,ja,ka)=(1.0D0/8.0D0)*(a(ia1-1,ja1-1,ka1-1)+a(ia1-1,ja1-1,ka1+1)+&
      a(ia1-1,ja1+1,ka1-1)+a(ia1-1,ja1+1,ka1+1)+&
      a(ia1+1,ja1-1,ka1-1)+a(ia1+1,ja1-1,ka1+1)+&
      a(ia1+1,ja1+1,ka1-1)+a(ia1+1,ja1+1,ka1+1))
		END DO
  END DO
END DO
DO ia=1,nx2i
  DO ja=1,ny2i
    DO ka=1,ny2i
      ao(ia,ja,ka)=ai(ia,ja,ka)
    END DO
  END DO
END DO
RETURN
#ifdef REALSWITCH
END SUBROUTINE interpol3Dn
#else
END SUBROUTINE interpol3Dcn
#endif
    
!----------------------------interpol3D(a)--------------------------------------
#ifdef REALSWITCH
SUBROUTINE interpol3D(a,ai)
#else
SUBROUTINE interpol3Dc(a,ai)
#endif
! 3D interpolation by a linear method
! with 1 interpolated point between 2 real point of each dimension
! Then there is 1 interpolated in a center of each elementary edge (2 real points),
! of each elementary face (4 real points) and of each elementary square (8 real points)

! This subroutine interpolates a nx2.ny2.nz2 dimension array to (2*nx2-1)(2*ny2-1)(2*nz2-1) array

!		illustration of the method in 1D
!
! 	    o--.--o--.--o--.--o--.--o--.--o--.--o--.--o		space discretisation
!
!	f:  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15		interpolated vector 
! 	F:  1     2     3     4     5     6     7     8		real vector 
! 	
!	f(n)=F(N) 			with an odd n=2N-1 (odd interpolated points correspond to real points)
!	f(n)=(1/2)[f(n-1)+f(n+1)] 	with an even n

USE params
IMPLICIT NONE
! INTEGER :: nx2,ny2,nz2!,nxyz=nx2*ny2*nz2

#ifdef REALSWITCH
!interpolation target array
! REAL(DP):: va(nxyz)
REAL(DP),INTENT (IN) :: a(nx2,ny2,nz2)
!interpolation variables
! REAL(DP),ALLOCATABLE ::ai(:,:,:)
REAL(DP),INTENT (OUT) ::ai(2*nx2,2*ny2,2*nz2)
#else
! COMPLEX(DP):: va(nxyz)
COMPLEX(DP),INTENT (IN) :: a(nx2,ny2,nz2)
!interpolation variables
! REAL(DP),ALLOCATABLE ::ai(:,:,:)
COMPLEX(DP),INTENT (OUT) ::ai(2*nx2,2*ny2,2*nz2)
#endif
INTEGER:: nx2i,ny2i,nz2i
INTEGER :: ia,ja,ka,i0

ai(:,:,:) = 0.0D0
!write (*,*)'hello1'
nx2i=2*nx2-1
ny2i=2*ny2-1
nz2i=2*nz2-1
! ALLOCATE(ai(nx2i,ny2i,nz2i))

!va is a vector no an array so the first step is to fill the corresponding 3D array a.
!#ifdef REALSWITCH
! call from1Dto3D(va,a,nx2,ny2,nz2) !***********
! #else
! call from1Dto3D(va,a,nx2,ny2,nz2) !***********
! #endif  

!-------------------------------------------------------
!filling the odd points of interpolated array with corresponding points of target array
DO ia=1,nx2
  DO ja=1,ny2
    DO ka=1,nz2
        ai(2*ia-1,2*ja-1,2*ka-1)=a(ia,ja,ka)
    END DO
  END DO 
END DO
    
!*******************************
!filling the even interpolated points in the center of elementary squares
DO ia=2,nx2i,2
  DO ja=2,ny2i,2
		DO ka=2,nz2i,2
      ai(ia,ja,ka)=(1.0D0/8.0D0)*(ai(ia-1,ja-1,ka-1)+ai(ia-1,ja-1,ka+1)+&
      ai(ia-1,ja+1,ka-1)+ai(ia-1,ja+1,ka+1)+&
      ai(ia+1,ja-1,ka-1)+ai(ia+1,ja-1,ka+1)+&
      ai(ia+1,ja+1,ka-1)+ai(ia+1,ja+1,ka+1))
    END DO
  END DO
END DO
    
!*******************************
!filling the even interpolated points in the center of elementary faces
DO ia=2,nx2i,2
	DO ja=2,ny2i,2
		DO ka=2,nz2i,2
      ai(ia,ja,ka)=(1.0D0/4.0D0)*(ai(ia,ja-1,ka-1)+ai(ia,ja-1,ka+1)+&
      ai(ia,ja+1,ka-1)+ai(ia,ja+1,ka+1))
    END DO
  END DO
END DO

DO ja=2,ny2i,2
	DO ka=2,nz2i,2
		DO ia=2,nx2i,2
      ai(ia,ja,ka)=(1.0D0/4.0D0)*(ai(ia-1,ja,ka-1)+ai(ia-1,ja,ka+1)+&
						    ai(ia+1,ja,ka-1)+ai(ia+1,ja,ka+1))
    END DO
  END DO
END DO

DO ka=2,nz2i,2
  DO ia=2,nx2i,2
		DO ja=2,ny2i,2
      ai(ia,ja,ka)=(1.0D0/4.0D0)*(ai(ia-1,ja-1,ka)+ai(ia-1,ja+1,ka)+&
			ai(ia+1,ja-1,ka)+ai(ia+1,ja+1,ka))
    END DO
  END DO
END DO

!*******************************
!filling the even interpolated points in the center of elementary edges
DO ja=2,ny2i,2
  DO ka=2,nz2i,2
		DO ia=2,nx2i,2
      ai(ia,ja,ka)=(1.0D0/2.0D0)*(ai(ia-1,ja,ka)+ai(ia+1,ja,ka))
    END DO
  END DO
END DO

DO ka=2,nz2i,2
	DO ia=2,nx2i,2
		DO ja=2,ny2i,2
      ai(ia,ja,ka)=(1.0D0/2.0D0)*(ai(ia,ja-1,ka)+ai(ia,ja+1,ka))
    END DO
  END DO
END DO

DO ia=2,nx2i,2
	DO ja=2,ny2i,2
		DO ka=2,nz2i,2
      ai(ia,ja,ka)=(1.0D0/2.0D0)*(ai(ia,ja,ka-1)+ai(ia,ja,ka+1))
    END DO
  END DO
END DO

!-------------------------------------------------------


!write (*,*)'end interpol'

RETURN
#ifdef REALSWITCH
END SUBROUTINE interpol3D
#else
END SUBROUTINE interpol3Dc
#endif

!----------------------------smoothing3D(a)--------------------------------------
#ifdef REALSWITCH
SUBROUTINE smoothing3D(a,as)
#else
SUBROUTINE smoothing3Dc(a,as)
#endif

!the smoothing method consists on the average of the 27 points. The 26 points around the real point count for the half of the contribution and the main point for the other half.
!
!how the 27 points are counted? e.g. in 2D, at z(k)
!                     y(j+1)   .     .     .
!
!                       y(j)   .     o     .
!
!                     y(j-1)   .     .     .
!
!                           x(i+1) x(i) x(i-1)

!   o is the real main point at  x(i),y(j),z(k)
!   . are the interpolated points so there are 8 at the z(k) plan, there are also 9 on each z(k-1) and z(k+1) plans

USE params
IMPLICIT NONE
#ifdef REALSWITCH
!smothing array
REAL(DP),INTENT (OUT) :: as(nx2,ny2,nz2)
!smoothing target array
INTEGER :: nx2i,ny2i,nz2i
REAL(DP),INTENT (IN) ::a(2*nx2,2*ny2,2*nz2)

#else
COMPLEX(DP),INTENT (OUT) :: as(nx2,ny2,nz2)
INTEGER :: nx2i,ny2i,nz2i
!smoothing target array
COMPLEX(DP),INTENT (IN) ::a(2*nx2,2*ny2,2*nz2)
#endif
INTEGER :: ia,ja,ka,i1,i2,i3
REAL(DP) :: fac
as=0.0

! nx2i=2*nx2-1
! ny2i=2*ny2-1
! nz2i=2*nz2-1
! 
!*******************************
DO ia=1,nx2
  IF ((ia==1).or.(ia==nx2))THEN ! for the orthogonal border plans to the x axis
    DO ja=1,ny2
      DO ka=1,nz2
        as(ia,ja,ka)=a(2*ia-1,2*ja-1,2*ka-1)
        END DO
      END DO
  END IF
END DO

DO ja=1,ny2
  IF ((ja==1).or.(ja==ny2))THEN ! for the orthogonal border plans to the y axis
    DO ka=1,nz2
      DO ia=1,nx2
    !    if((ia.ne.1).and.(ia.ne.nx2))then
        as(ia,ja,ka)=a(2*ia-1,2*ja-1,2*ka-1)
    !    end if
      END DO
    END DO
  END IF
END DO

DO ka=1,nz2
  IF ((ka==1).or.(ka==nz2))THEN ! for the orthogonal border plans to the z axis
    DO ia=1,nx2
      !if((ia.ne.1).and.(ia.ne.nx2))then
      DO ja=1,ny2
        IF((ja.ne.1).and.(ja.ne.ny2)) as(ia,ja,ka)=a(2*ia-1,2*ja-1,2*ka-1)
      END DO
      !end if
    END DO
  END IF
END DO
!*******************************
!for the inner points
DO ia=2,nx2-1
  DO ja=2,ny2-1
    DO ka=2,nz2-1
      DO i1=-2,0
        DO i2=-2,0
          DO i3=-2,0
            fac=1.0D0/27.D0
            IF(i1.eq.-1) THEN
              IF(i2.eq.-1) THEN
                IF(i3.eq.-1) THEN
                  fac=26.0D0/27.0D0
                END IF
              END IF
            END IF
            as(ia,ja,ka)=as(ia,ja,ka)+fac* a(2*ia+i1,2*ja+i2,2*ka+i3)
          END DO
        END DO
			END DO
    END DO
  END DO
END DO

! todo : Compare speed with this form :
!~ !REAL(DP):: fac(-2:1, -2:1, -2:1)

!~ fac = 1.0D0/27.D0 
!~ fac(-1,-1,-1)= 26.0D0/27.0D0

!~ DO ka=2,nz2-1
!~   DO ja=2,ny2-1
!~     DO ia=2,nx2-1
!~       DO i3=-2,0
!~         DO i2=-2,0
!~           DO i1=-2,0
!~             as(ia,ja,ka)=as(ia,ja,ka)+ fac(i1,i2,i3) * a(2*ia+i1,2*ja+i2,2*ka+i3)
!~           END DO
!~         END DO
!~ 			END DO
!~     END DO
!~   END DO
!~ END DO

RETURN
#ifdef REALSWITCH
END SUBROUTINE smoothing3D
#else
END SUBROUTINE smoothing3Dc
#endif
