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
SUBROUTINE CrankNicolson_exp(q0,aloc,rho,it,qwork)

!     Propagation with the Crank-Nicolson schem

!     The problem is written as a linear system
!       Ax=y
!     Where y = [1-iH(t)dt/2-...]*psi(t); A= 1+iH(t+dt)dt/2+...(taylor series); x=psi(t+dt)
!     at each time step:
!          1) compute y
!          2) predict H and psi at t+dt
!          3) compute the error Ax-y and verify that its norm tends to 0.threshold criterion.
!             other wise change x to x-err.r; where r is the chosen relaxation factor
!             predict H(t+dt)
!             do this step until norm(Ax-y)> threshold then do at the next time step.

!
!      q0          = s.p. wavefunctions to be propagated
!      aloc        = array for local mean field
!      rho         = array for local density
!      it          = nr. of time step
!      qwork       = work space for s.p. wavefunctions

USE params
#if(twostsic)
USE twost, ONLY:tnearest
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(INout)                  :: q0(kdfull2,kstate)
REAL(DP), INTENT(IN OUT)                 :: aloc(2*kdfull2)
!REAL(DP), INTENT(IN OUT)                :: akv(kdfull2)
REAL(DP), INTENT(IN OUT)                 :: rho(2*kdfull2)
INTEGER, INTENT(IN)                      :: it
COMPLEX(DP), INTENT(OUT)                 :: qwork(kdfull2,kstate)

COMPLEX(DP) :: q1(kdfull2)
COMPLEX(DP) :: cdtact
COMPLEX(DP)                     :: y(kdfull2,kstate),errcn(kdfull2,kstate),x(kdfull2,kstate)

! The parameter 'tnorotate' de-activates the subtraction of the
! Lagrangian matrix in the SIC step. The version of exponential
! evolution with subtraction of the Lagrangian matrix is found
! in 'exp_evolp'. The strategy needs yet testing. It was not
! really beneficial so far.
LOGICAL,PARAMETER :: tnorotate=.true.
#if(parayes||simpara||paraworld)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif
#if(paraworld)
COMPLEX(DP),ALLOCATABLE :: qarray (:,:,:),qarrayfine (:,:,:)
COMPLEX(DP),ALLOCATABLE :: qactfine(:)
ALLOCATE(qactfine(kdfull2fine))
ALLOCATE(qarray (nx2,ny2,nz2),qarrayfine (2*nx2,2*ny2,2*nz2))
#endif

  cdtact = CMPLX(dt1/2.0,0D0)

!----------------------------------------------------------------------



#if(paraworld)
CALL  mpi_comm_rank(mpi_comm_world,nrank,icode)
level=nrank
#endif


#if(parayes)
STOP 'exponential evolution not yet MPI parallelized'
#else
myn = 0
#endif


#if(raregas)
IF(nc+NE+nk > 0) STOP 'TSTEP_EXP not appropriate for rare gas'
#endif



!     one half time step to define new mean field
!     use exponential evolution to second order

DO nb=1,nstate
  CALL exp_evol(q0(1,nb),aloc,nb,1,cdtact,q1)
END DO
! save the (1-iH(t+dt/2)dt/2-...)*psi(t) for the right term of Crank-Nicolson scheme

y=q0

CALL dyn_mfield(rho,aloc,y,dt1)


DO nb=1,nstate
  !  qwork(:,nb) = q0(:,nb)
  CALL exp_evol(q0(1,nb),aloc,nb,1,cdtact,q1)
END DO

!      second half time step to estimate x and H(t+dt)
!estimate mean field of full time step
!CALL dyn_mfield(rho,aloc,qwork,dt1)
x=q0
!     One negative half time step to initialize  Ax
!DO nb=1,nstate
!  CALL exp_evol(qwork(1,nb),aloc,nb,1,-cdtact,q1)
!END DO

!error
rnorm_err=1.0
relax=0.25
threshold1=1e-1
threshold=1e-7
iter_CN=0





#if(paraworld)
CALL  mpi_comm_size(mpi_comm_world,nprocs,icode)
CALL  mpi_comm_rank(mpi_comm_world,level,icode)
kdfull28=(kxbox/2)*(kybox/2)*(kzbox/2)

        if(level.ne.0) then
   do nb=1,nstate
call mpi_recv(q0(1,nb),kdfull2,mpi_double_complex,level-1,1,mpi_comm_World,is,ic)
enddo
call mpi_recv(dt1,1,mpi_double_precision,level-1,3,mpi_comm_World,is,ic)
endif



#endif
! begin of self consistent Crank Nicolson relaxation loop 
do while (rnorm_err.gt.threshold1) 
      errcn= q0-y
      ro= rnorm_err
      rnorm_err= 0.0
   do nb=1,1
      rnorm_err=rnorm_err+wfnorm(errcn(1,nb))
   enddo
   write(6,*) 'err in crank step  1',rnorm_err,level,relax
   if(ro.lt.rnorm_err) then
               relax=relax*0.7
   else
               relax=relax*1.2
   endif
   x=x-relax*errcn
   q0=x
   iter_CN=iter_CN+1
   if (iter_CN.gt.25)then
     write(*,*) ' warning Crank-Nicolson 1',iter_CN
     rnorm_err=0
     else
   if (iter_CN.lt.1)then
     rnorm_err=1.0
     else
   if(mod(iter_CN,30).eq.0) CALL dyn_mfield(rho,aloc,x,dt1)

!     One negative half time step to compute Ax
   do nb=1,nstate
   CALL exp_evol(q0(1,nb),aloc,nb,1,-cdtact,q1)
   enddo
   endif
   endif
end Do

! end of first self consistent Crank Nicolson relaxation loop 

q0=x


#if(paraworld)
        if(level.ne.nprocs-1) then
do nb=1,nstate

call from1Dto3Dc(q0(1,nb),qarrayfine,nx2,ny2,nz2)
call smoothing3Dc(qarrayfine,qarray)
!from fine array to coarse array
call from3Dto1Dc(q1,qarray,nx2/2,ny2/2,nz2/2)
!from coarse array to coarse !vector
call mpi_send(q1,kdfull28,mpi_double_complex,level+1,1,mpi_comm_World,ic)
call mpi_send(dt1,1,mpi_double_precision,level+1,3,mpi_comm_World,ic)
enddo
endif




!begin of step 2 

if(level.ne.nprocs-1) then
do nb=1,nstate
call mpi_recv(q0(1,nb),kdfull2,mpi_double_complex,level+1,2,mpi_comm_World,is,ic)
enddo
endif


do while (rnorm_err.gt.threshold) 
      errcn= q0-y
      ro=rnorm_err
      rnorm_err= 0.0
   do nb=1,1
      rnorm_err=rnorm_err+wfnorm(errcn(1,nb))
   enddo
   write(6,*) 'err in crank step 2',rnorm_err,level,relax
   if(ro.lt.rnorm_err) then
               relax=relax*0.7
   else
               relax=relax*1.2
   endif
   x=x-relax*errcn
   q0=x
   iter_CN=iter_CN+1
   if(mod(iter_CN,30).eq.0) CALL dyn_mfield(rho,aloc,x,dt1)
   if (iter_CN.gt.150*(level+1))then
     write(*,*) ' warning Crank-Nicolson0',iter_CN
     rnorm_err=0
      
     if(level.eq.0) then
     dt1=dt1*1.00
     write(6,*) 'reduce dt to',dt1
     cdtact = CMPLX(dt1/2.0,0D0)
     endif
     else
   if (iter_CN.lt.1)then
     rnorm_err=1.0
     else

!     One negative half time step to compute Ax
   do nb=1,nstate
   CALL exp_evol(q0(1,nb),aloc,nb,1,-cdtact,q1)
   enddo
   endif
   endif
end Do
   if (iter_CN.lt.150*(level+1))then
     if(level.eq.0) then
     dt1=dt1*1.00
     write(6,*) 'increase dt to',dt1
     cdtact = CMPLX(dt1/2.0,0D0)
     endif
     endif

! end of second self consistent Crank Nicolson relaxation loop 

if(level.ne.0) then 
do nb=1,nstate
call from1Dto3Dc(q0(1,nb),qarray,nx2,ny2,nz2)    !from coarse vector to coarse array
call interpol3Dc(qarray,qarrayfine)               !from coarse array to fine array
call from3Dto1Dc(qactfine,qarrayfine,2*nx2,2*ny2,2*nz2)
!from coarse array to coarse !vector
call mpi_send(qactfine,kdfull2fine,mpi_double_complex,level-1,2,mpi_comm_World,ic)
enddo
endif



q0=x

#endif





!     compute mean field at new time

!CALL dyn_mfield(rho,aloc,qwork,dt1)

RETURN
END SUBROUTINE CrankNicolson_exp
