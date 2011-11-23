#include "define.h"
 
PROGRAM testcoul

!     test spped of Coulomb solver


!     *******************************************
!     declarations
!     *******************************************

USE params
USE kinetic
USE coulsolv
#if(twostsic)
USE twostr
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

!     non symmetrical dynamic fields

!     psir = real wavefunctions in coord space
!     psi  = complex wavefunctions in coord space (for dynamics)
!     rho  = electronic density in coord space
!     aloc = mean-field-potential (to be initialized before dynamics)
!     chpcoul = coulomb-potential of electrons
!     c$s/cp$s = auxiliary field to store coords and momenta in between
!              the trial protonic propagation
!     rhojel = jel.density

#if(parayes||simpara)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif


INTEGER :: getnearestgridpoint
INTEGER :: conv3to1


REAL(DP),ALLOCATABLE :: aloc(:),rho(:)
!REAL(DP),ALLOCATABLE :: rhotmp(:)

REAL(DP),ALLOCATABLE :: psir(:,:)
COMPLEX(DP),ALLOCATABLE :: psi(:,:)
COMPLEX(DP),ALLOCATABLE :: psiw(:,:)

!EQUIVALENCE (psir,psi)                    ! recycle storage space



LOGICAL :: tprints,tfinals,tmf

!REAL(4) tarray(2)
!REAL(4) etime



!     *******************************************

!     init

!     *******************************************


CALL cpu_time(time_absinit)

CALL init_parallele()

CALL checkoptions()

CALL initnamelists          ! read all input parameters

CALL init_baseparams()

CALL initisrtyp

!CALL iperio                     ! initializing the 'periodic table'

!CALL changeperio   ! overwrites default periodic system if necessary

CALL iparams()               ! check dynamic  parameters

CALL init_grid()

CALL init_fields()

!ALLOCATE(psir(kdfull2,kstate))
!psir=0D0
ALLOCATE(aloc(2*kdfull2),rho(2*kdfull2))
aloc=0D0
rho=0D0




IF(myn == 0) CALL ocoption(7)   ! output compiled options
IF(myn == 0) CALL ocoption(8)   ! output compiled options
CALL init_output()              ! headers for basic output files

CALL init_jellium()

rho=rhojel

DO i=1,100
CALL coul_mfield(rho)
WRITE(*,*) ' COULOMB nr. ',i
END DO


END PROGRAM testcoul


!#include "define.h"

#if(raregas)
!     ************************************

SUBROUTINE loc_mfield_dummy(rho,aloc)

!     ************************************

!     This routine is a dummy version -- usage still unclear.
!     Most probably obsolete!!!

!     The local part of the mean field
!     plus an update of the pseudo-potentials.

!     Input:
!      rho    = electron density
!      dt     = stepsize (in case of dynamics)
!      tdyn   = switch to dynamic case
!     Output:
!      aloc   = local mean field

USE params
USE kinetic
USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)


REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)



REAL(DP),ALLOCATABLE :: rhotmp(:)
COMPLEX(DP) :: psidummy(1)


ALLOCATE(rhotmp(2*kdfull2))
IF(idielec /= 1) THEN
  CALL calcpseudo()
  CALL calcrho(rho,psi)
  DO ii=1,2*kdfull2
    rhotmp(ii)=rho(ii)
  END DO
  CALL addimage(rho,1)
  
! Coulomb of the electronic density
#if(gridfft)
  CALL falr(rho,chpcoul,nx2,ny2,nz2,kdfull2)
#endif
#if(findiff|numerov)
  CALL solv_fft(rho,chpcoul,dx,dy,dz)
#endif
  CALL calclocal(rho,aloc)
END IF

IF(idielec == 1) THEN
  
  DO ii=1,2*kdfull2
    rho(ii)=rhotmp(ii)
  END DO
  DO ii=1,kdfull2
    rfieldtmp(ii)=chpcoul(ii)
  END DO
  CALL addimage(rho,0)
  
#if(gridfft)
!         if(ipseudo.eq.1)
  CALL falr(rho,chpcoul,nx2,ny2,nz2,kdfull2)
#endif
#if(findiff|numerov)
!         if(ipseudo.eq.1)
  CALL solv_fft(rho,chpcoul,dx,dy,dz)
#endif
  
  
  DO ii=1,kdfull2
    CALL conv1to3(ii)
    IF(iindtmp(1) > nint(xdielec/dx)+nx)THEN
      chpcoul(ii)=rfieldtmp(ii)
    END IF
  END DO
  DO ii=1,2*kdfull2
    rho(ii)=rhotmp(ii)
  END DO
  
END IF

CALL calclocal(rho,aloc)
IF(ifsicp > 0) CALL calc_sic(rho,aloc,psi)

DEALLOCATE(rhotmp)

RETURN
END SUBROUTINE loc_mfield_dummy
#endif
