#include "define.h"
 
!-----calcrhot------------------------------------------------------

SUBROUTINE calcrhot(rho,q0)

!     density 'rho' for complex or real wavefunctions 'q0'

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

CLASS(*),TARGET, INTENT(IN) :: q0(kdfull2,kstate)
REAL(DP),TARGET, INTENT(OUT) :: rho(2*kdfull2)

!REAL(DP),DIMENSION(:),POINTER :: rh


!-----------------------------------------------------------------

!rh => rho

rho=0D0

DO nb=1,nstate
  ishift = (ispin(nrel2abs(nb))-1)*nxyz ! store spin=2 in upper block
  SELECT TYPE(q0)
   TYPE IS(REAL(DP)) 
    DO ind=1,nxyz
      rho(ind+ishift)=rho(ind+ishift)+ occup(nb)*q0(ind,nb)**2
    END DO
   TYPE IS(COMPLEX(DP)) 
    DO ind=1,nxyz
      rho(ind+ishift)=rho(ind+ishift)+ occup(nb)*&
       (REAL(q0(ind,nb))**2+AIMAG(q0(ind,nb))**2)
    END DO
  END SELECT 
END DO

!     reorder to total density in lower block (1-nxyz)
!     and difference density in upper block (nxyz+1-2nxyz)
!DO ind=1,nxyz
!#if(fullspin)
!  rhootot      = rhoo(ind) + rhoo(ind+nxyz)
!  rhoodif      = rhoo(ind) - rhoo(ind+nxyz)
!#else
!  rhootot      = rhoo(ind)
!  rhoodif      = 0D0
!#endif
!  rhoo(ind)      = rhootot
!  rhoo(ind+nxyz) = rhoodif/MAX(rhootot,1D-8)
!END DO
!
!CALL emoms(rhoo)



RETURN
END SUBROUTINE calcrhot
