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
COMPLEX(DP), INTENT(IN) :: q0(kdfull2,kstate)
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
#if(fullspin)
  rhotot      = rho(ind) + rho(ind+nxyz)
  rhodif      = rho(ind) - rho(ind+nxyz)
#else
  rhotot      = rho(ind)
  rhodif      = 0D0
#endif
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

