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

!     reorder to total density in lower block (1-nxyz)
!     and difference density in upper block (nxyz+1-2nxyz)

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

COMPLEX(DP), INTENT(IN OUT) :: q0(kdfull2,kstate)
REAL(DP), INTENT(OUT) :: current(kdfull2,3)

COMPLEX(DP), ALLOCATABLE :: dq0(:)

#if(parayes)
STOP "CALC_CURRENT presently not suited for parallel computing"              ! cPW
#endif

!-----------------------------------------------------------------

ALLOCATE(dq0(kdfull2))

! reset 
current=0D0

! cumulate
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

