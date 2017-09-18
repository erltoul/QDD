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
SUBROUTINE cranknicolson_exp(q0,aloc,rho,it)

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
  USE util

#if(twostsic)
  USE twost, ONLY:tnearest
#endif

  IMPLICIT NONE

  COMPLEX(DP), INTENT(INOUT)               :: q0(kdfull2,kstate)
  REAL(DP)   , INTENT(INOUT)               :: aloc(2*kdfull2)
  REAL(DP)   , INTENT(INOUT)               :: rho(2*kdfull2)
  INTEGER    , INTENT(IN)                  :: it
!  COMPLEX(DP), INTENT(OUT)                 :: qwork(kdfull2,kstate)

  INTEGER :: iter_cn,nb
  REAL(DP) :: ro, rnorm_err, relax, threshold1, threshold

  COMPLEX(DP) :: q1(kdfull2)
  COMPLEX(DP) :: cdtact
  COMPLEX(DP),DIMENSION(KDFULL2,KSTATE)     :: y,errcn,x

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
  INTEGER :: kdfull28,level, nprocs, nrank, ilocbas
  COMPLEX(DP) :: hpsisave(kdfull2,kstate)
  COMPLEX(DP),ALLOCATABLE :: qarray (:,:,:),qarrayfine (:,:,:)
  COMPLEX(DP),ALLOCATABLE :: qactfine(:)
  ALLOCATE(qactfine(kdfull2fine))
  ALLOCATE(qarray (nx2,ny2,nz2),qarrayfine (2*nx2,2*ny2,2*nz2))
#endif

  cdtact = CMPLX(dt1/2D0,0D0)

!----------------------------------------------------------------------

#if(paraworld)
  CALL  mpi_comm_rank(mpi_comm_world,nrank,mpi_ierror)
  !CALL  MPI_COMM_SIZE(mpi_comm_world,nprocs,mpi_ierror)
  level=nrank
  !WRITE(2000*(level+1)+6,*) level,nprocs
  !CALL flush(2000*level)
  !CALL MPI_FINALIZE(mpi_ierror)
  !STOP
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

  DO nb = 1,nstate
     CALL exp_evol(q0(1,nb),aloc,nb,1,cdtact,q1)
  END DO
! save the (1-iH(t+dt/2)dt/2-...)*psi(t) for the right term of Crank-Nicolson scheme

  y = q0

  CALL dyn_mfield(rho,aloc,y,dt1,it)

  DO nb = 1,nstate
     !  qwork(:,nb) = q0(:,nb)
     CALL exp_evol(q0(1,nb),aloc,nb,1,cdtact,q1)
  END DO

!      second half time step to estimate x and H(t+dt)
!estimate mean field of full time step
!CALL dyn_mfield(rho,aloc,qwork,dt1,it)
  x = q0
!     One negative half time step to initialize  Ax
!DO nb = 1,nstate
!  CALL exp_evol(qwork(1,nb),aloc,nb,1,-cdtact,q1)
!END DO

!error
  rnorm_err = 1D0
  relax = 0.25D0
  threshold1 = 1D-2
  threshold = 1D-9
  iter_CN = 0

#if(paraworld)
  CALL  mpi_comm_size(mpi_comm_world,nprocs,mpi_ierror)
  CALL  mpi_comm_rank(mpi_comm_world,level,mpi_ierror)
  kdfull28 = (kxbox/2)*(kybox/2)*(kzbox/2)

  IF(level.NE.0) THEN
     DO nb = 1,nstate
        CALL mpi_recv(q0(1,nb),kdfull2,mpi_double_complex,level-1,1,mpi_comm_world,is,mpi_ierror)
     ENDDO
     CALL mpi_recv(dt1,1,mpi_double_precision,level-1,3,mpi_comm_world,is,mpi_ierror)
  ENDIF
#endif

! begin of self consistent Crank Nicolson relaxation loop 
  DO WHILE (rnorm_err.GT.threshold1) 
     errcn = q0-y
     ro = rnorm_err
     rnorm_err = 0D0
     DO nb = 1,1
        rnorm_err = rnorm_err+wfnorm(errcn(:,nb))
     ENDDO
!#if(paraworld)
     !   WRITE(6,*) 'err in crank step  1',rnorm_err,level,relax
!#endif
     IF(ro.LT.rnorm_err) THEN
        relax = relax*0.7D0
     ELSE
        relax = relax*1.2D0
     ENDIF
     x = x - relax*errcn
     q0 = x
     iter_CN = iter_CN + 1
     IF(iter_CN.GT.20)THEN
        WRITE(*,*) ' warning Crank-Nicolson 1',iter_CN
        rnorm_err = 0D0
     ELSE
        IF(iter_CN.LT.1)THEN
           rnorm_err = 1D0
        ELSE
           IF(MOD(iter_CN+1,30).EQ.0) CALL dyn_mfield(rho,aloc,x,dt1,it)

!     One negative half time step to compute Ax
           DO nb = 1,nstate
              CALL exp_evol(q0(1,nb),aloc,nb,1,-cdtact,q1)
           ENDDO
        ENDIF
     ENDIF
  ENDDO
! end of first self consistent Crank Nicolson relaxation loop 

  q0 = x

#if(paraworld)
  IF(level.NE.nprocs-1) THEN
     DO nb = 1,nstate
        CALL from1Dto3Dc(q0(1,nb),qarrayfine,nx2,ny2,nz2)
        CALL smoothing3Dc(qarrayfine,qarray)
!from fine array to coarse array
        CALL from3Dto1Dc(q1,qarray,nx2/2,ny2/2,nz2/2)
!from coarse array to coarse !vector
        CALL mpi_send(q1,kdfull28,mpi_double_complex,level+1,1,mpi_comm_world,mpi_ierror)
     ENDDO
     CALL mpi_send(dt1,1,mpi_double_precision,level+1,3,mpi_comm_world,mpi_ierror)
  ENDIF

  iter_CN = 0
  
!--------------------  
!begin of step 2 
!--------------------
  IF(level.NE.nprocs-1) THEN
     DO nb = 1,nstate
        CALL mpi_recv(q0(1,nb),kdfull2,mpi_double_complex,level+1,2,mpi_comm_world,is,mpi_ierror)
        CALL mpi_recv(hpsisave(1,nb),kdfull2,mpi_double_complex,level+1,12,mpi_comm_world,is,mpi_ierror)
     ENDDO
     CALL mpi_recv(aloc,kdfull2,mpi_double_precision,level+1,13,mpi_comm_world,is,mpi_ierror)
     CALL mpi_recv(aloc(1+nx2*ny2*nz2),kdfull2,mpi_double_precision,level+1,13,mpi_comm_world,is,mpi_ierror)
  ENDIF

  DO WHILE (rnorm_err.GT.threshold) 
     errcn =  q0 - y
     ro = rnorm_err
     rnorm_err =  0D0
     DO nb = 1,1
        rnorm_err = rnorm_err + wfnorm(errcn(:,nb))
     ENDDO
     WRITE(6,*) 'err in crank step 2',rnorm_err,level,relax,iter_CN
     IF(ro.LT.rnorm_err) THEN
        relax = relax*0.7D0
     ELSE
        relax = relax*1.2D0
     ENDIF
     x = x - relax*errcn
     q0 = x
     iter_CN = iter_CN + 1
     IF(MOD(iter_CN,30).EQ.0) CALL dyn_mfield(rho,aloc,x,dt1,it)
     IF(iter_CN.GT.150*(level+1))THEN
        WRITE(*,*) ' warning Crank-Nicolson0',iter_CN
        rnorm_err = 0
        IF(level.EQ.0) THEN
           dt1 = dt1*1D0
           WRITE(6,*) 'reduce dt to',dt1
           cdtact  =  CMPLX(dt1/2D0,0D0)
        ENDIF
     ELSE
        IF(iter_CN.LT.1)THEN
           rnorm_err = 1D0
        ELSE
!     One negative half time step to compute Ax
           DO nb = 1,nstate
              CALL exp_evol(q0(1,nb),aloc,nb,1,-cdtact,q1)
           ENDDO
        ENDIF
     ENDIF
  ENDDO
  IF(iter_CN.LT.150*(level+1))THEN
     IF(level.eq.0)THEN
        dt1 = dt1*1D0
        WRITE(6,*) 'increase dt to',dt1
        cdtact = CMPLX(dt1/2D0,0D0)
     ENDIF
  ENDIF

! end of second self consistent Crank Nicolson relaxation loop 

  IF(level.NE.0) THEN 
     DO nb = 1,nstate
        CALL from1Dto3Dc(q0(1,nb),qarray,nx2,ny2,nz2)    !from coarse vector to coarse array
        CALL interpol3Dc(qarray,qarrayfine)               !from coarse array to fine array
        CALL from3Dto1Dc(qactfine,qarrayfine,2*nx2,2*ny2,2*nz2)
        CALL mpi_send(qactfine,kdfull2fine,mpi_double_complex,level-1,2,mpi_comm_world,mpi_ierror)
!from coarse array to coarse !vector
        IF (ispin(nrel2abs(nb)) == 1) THEN
           ilocbas = 1
        ELSE IF (ispin(nrel2abs(nb)) == 2) THEN
           ilocbas = nx2*ny2*nz2 + 1
        ELSE
           STOP " EXPEVOL: spin index must be 1 or 2"
        END IF

        CALL hpsi(q0(1,nb),aloc(ilocbas),nb,1)

        CALL from1Dto3Dc(hpsisave(1,nb),qarray,nx2,ny2,nz2)  !from coarse vector to coarse array
        CALL interpol3Dc(qarray,qarrayfine)  !from coarse array to fine array
        CALL from3Dto1Dc(qactfine,qarrayfine,2*nx2,2*ny2,2*nz2)
        CALL mpi_send(qactfine,kdfull2fine,mpi_double_complex,level-1,12,mpi_comm_world,mpi_ierror)
        !from coarse array to coarse !vector

     ENDDO
     CALL from1Dto3D(aloc,qarray,nx2,ny2,nz2)  !from coarse vector to coarse array
     CALL interpol3D(qarray,qarrayfine)  !from coarse array to fine array
     CALL from3Dto1D(qactfine,qarrayfine,2*nx2,2*ny2,2*nz2)
     CALL mpi_send(qactfine,kdfull2fine,mpi_double_precision,level-1,13,mpi_comm_world,mpi_ierror)
     CALL from1Dto3D(aloc(nx2*ny2*nz2+1),qarray,nx2,ny2,nz2)  !from coarse vector to coarse array
     CALL interpol3D(qarray,qarrayfine)  !from coarse array to fine array
     CALL from3Dto1D(qactfine,qarrayfine,2*nx2,2*ny2,2*nz2)
     CALL mpi_send(qactfine,kdfull2fine,mpi_double_precision,level-1,13,mpi_comm_world,mpi_ierror)
ENDIF

q0 = x

#endif

!     compute mean field at new time
!CALL dyn_mfield(rho,aloc,qwork,dt1)
RETURN
END SUBROUTINE cranknicolson_exp
