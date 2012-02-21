#include "define.h"
 
!-----zeroforce-----------------------------------------------------------------

SUBROUTINE zeroforce(aloc,rho)

!     Corrects SIC-Slater and SIC-KLI potentials to meet the
!     zero-force condition.

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
#if(gridfft)

REAL(DP), INTENT(OUT)                        :: aloc(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
!REAL(DP), INTENT(IN)                         :: akv(kdfull2)



!       workspaces

COMPLEX(DP), ALLOCATABLE :: potk(:),dervk(:)     ! for Fourier transformed potentials
REAL(DP), ALLOCATABLE :: potwork(:)

!      equivalence (potwork(1),w1(1))
!      equivalence (potk(1),w2(1))              ! occupies also w3

!-----------------------------------------------------------------------------


!     check workspace

!      if(usew1) stop ' in SSTEP: workspace W1 already active '
!      if(usew2) stop ' in SSTEP: workspace W2 already active '
!      if(usew3) stop ' in SSTEP: workspace W3 already active '
!      usew1 = .true.
!      usew2 = .true.
!      usew3 = .true.
ALLOCATE(potk(kdfull2),dervk(kdfull2),potwork(kdfull2))

!     loop over spin-up and spin-down

DO is=1,numspin
    ishift = (is-1)*nxyz
    
!       Fourier transformation
    
    CALL rftf(aloc(1+ishift),potk)
    
!       Laplacian and integral
    
    DO  i=1,nxyz
      dervk(i) = potk(i)*akv(i)
    END DO
    CALL rfftback(dervk,potwork)
    denominator = rwfovlp(rho(1+ishift),potwork)
    
!       x-derivative
    
    ind=0
    DO i3=1,nz2
      DO i2=1,ny2
        DO i1=1,nx2
          IF(i1 >= (nx+1)) THEN
            zkx=(i1-nx2-1)*dkx
          ELSE
            zkx=(i1-1)*dkx
          END IF
!MB:
          IF(i1 >= (nx+1)) THEN
            zkx=(i1-nx2-1)*dkx
          ELSE
            zkx=(i1-1)*dkx
          END IF
!MB/
          ind=ind+1
          dervk(ind) = eye*zkx*potk(ind)
        END DO
      END DO
    END DO
    CALL rfftback(dervk,potwork)
    counter = rwfovlp(rho(1+ishift),potwork)
    xlambda = counter/denominator
    DO i=1,nxyz
      aloc(i+ishift) = aloc(i+ishift)-xlambda*potwork(i)
    END DO
    
!       y-derivative
    
    ind=0
    DO i3=1,nz2
      DO i2=1,ny2
        IF(i2 >= (ny+1)) THEN
          zky=(i2-ny2-1)*dky
        ELSE
          zky=(i2-1)*dky
        END IF
        DO i1=1,nx2
!MB:
          IF(i1 >= (nx+1)) THEN
            zkx=(i1-nx2-1)*dkx
          ELSE
            zkx=(i1-1)*dkx
          END IF
!MB/
          ind=ind+1
          dervk(ind) = eye*zky*potk(ind)
        END DO
      END DO
    END DO
    CALL rfftback(dervk,potwork)
    counter = rwfovlp(rho(1+ishift),potwork)
    ylambda = counter/denominator
    DO i=1,nxyz
      aloc(i+ishift) = aloc(i+ishift)-ylambda*potwork(i)
    END DO
    
!       z-derivative
    
    ind=0
    DO i3=1,nz2
      IF(i3 >= (nz+1)) THEN
        zkz=(i3-nz2-1)*dkz
      ELSE
        zkz=(i3-1)*dkz
      END IF
      DO i2=1,ny2
        DO i1=1,nx2
!MB:
          IF(i1 >= (nx+1)) THEN
            zkx=(i1-nx2-1)*dkx
          ELSE
            zkx=(i1-1)*dkx
          END IF
!MB/
          ind=ind+1
          dervk(ind) = eye*zkz*potk(ind)
        END DO
      END DO
    END DO
    CALL rfftback(dervk,potwork)
    counter = rwfovlp(rho(1+ishift),potwork)
    zlambda = counter/denominator
    DO i=1,nxyz
      aloc(i+ishift) = aloc(i+ishift)-zlambda*potwork(i)
    END DO
    
  END DO
  DEALLOCATE(potk,dervk,potwork)
#else
  STOP "ZERO-FORCE not yet implemented for finites differences"
#endif
  
  RETURN
END SUBROUTINE zeroforce

!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE checkzeroforce(rho,aloc)
!------------------------------------------------------------
USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)


#if(gridfft)

REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)

!REAL(DP) :: akv(kdfull2)

!       workspaces

COMPLEX(DP), ALLOCATABLE :: potk(:),dervk(:)     ! for Fourier transformed potentials
REAL(DP), ALLOCATABLE :: potwork(:)

ALLOCATE(potk(kdfull2),dervk(kdfull2),potwork(kdfull2))

DO is=1,numspin
    ishift = (is-1)*nxyz
    
!       Fourier transformation
    
    CALL rftf(aloc(1+ishift),potk)
    
    
    
!       x-derivative
    
    ind=0
    DO i3=1,nz2
      DO i2=1,ny2
        DO i1=1,nx2
          IF(i1 >= (nx+1)) THEN
            zkx=(i1-nx2-1)*dkx
          ELSE
            zkx=(i1-1)*dkx
          END IF
!MB:
          IF(i1 >= (nx+1)) THEN
            zkx=(i1-nx2-1)*dkx
          ELSE
            zkx=(i1-1)*dkx
          END IF
!MB/
          ind=ind+1
          dervk(ind) = eye*zkx*potk(ind)
        END DO
      END DO
    END DO
    CALL rfftback(dervk,potwork)
    zforcex = rwfovlp(rho(1+ishift),potwork)
    
    
!       y-derivative
    
    ind=0
    DO i3=1,nz2
      DO i2=1,ny2
        IF(i2 >= (ny+1)) THEN
          zky=(i2-ny2-1)*dky
        ELSE
          zky=(i2-1)*dky
        END IF
        DO i1=1,nx2
!MB:
          IF(i1 >= (nx+1)) THEN
            zkx=(i1-nx2-1)*dkx
          ELSE
            zkx=(i1-1)*dkx
          END IF
!MB/
          ind=ind+1
          dervk(ind) = eye*zky*potk(ind)
        END DO
      END DO
    END DO
    CALL rfftback(dervk,potwork)
    zforcey = rwfovlp(rho(1+ishift),potwork)
    
    
!       z-derivative
    
    ind=0
    DO i3=1,nz2
      IF(i3 >= (nz+1)) THEN
        zkz=(i3-nz2-1)*dkz
      ELSE
        zkz=(i3-1)*dkz
      END IF
      DO i2=1,ny2
        DO i1=1,nx2
!MB:
          IF(i1 >= (nx+1)) THEN
            zkx=(i1-nx2-1)*dkx
          ELSE
            zkx=(i1-1)*dkx
          END IF
!MB/
          ind=ind+1
          dervk(ind) = eye*zkz*potk(ind)
        END DO
      END DO
    END DO
    CALL rfftback(dervk,potwork)
    zforcez = rwfovlp(rho(1+ishift),potwork)
    
    
    
  END DO
  
  DEALLOCATE(potk,dervk,potwork)
#else
  STOP "ZERO-FORCE not yet implemented for finites differences"
#endif
  
  WRITE(6,'(a,3e17.7)') 'Checking Zero-Force Theorem: ',  &
      zforcex,zforcey,zforcez
  
  WRITE(500,'(1f17.5,3e17.7)') tfs,zforcex,zforcey,zforcez
  
  RETURN
END SUBROUTINE checkzeroforce
!------------------------------------------------------------


!------------------------------------------------------------

SUBROUTINE getgradfield(fieldfrom,fieldto,ishift,icoor)
!------------------------------------------------------------
USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)


#if(gridfft)


REAL(DP), INTENT(IN OUT)                     :: fieldfrom(2*kdfull2)
REAL(DP), INTENT(OUT)                        :: fieldto(kdfull2)
INTEGER, INTENT(IN OUT)                  :: ishift
INTEGER, INTENT(IN)                      :: icoor



COMPLEX(DP), ALLOCATABLE :: potk(:),dervk(:)     ! for Fourier transformed potentials
REAL(DP), ALLOCATABLE :: potwork(:)

ALLOCATE(potk(kdfull2),dervk(kdfull2),potwork(kdfull2))


!       Fourier transformation

CALL rftf(fieldfrom(1+ishift*kdfull2),potk)


IF (icoor == 1) THEN
  
!       x-derivative
  
  ind=0
  DO i3=1,nz2
    DO i2=1,ny2
      DO i1=1,nx2
        IF(i1 >= (nx+1)) THEN
          zkx=(i1-nx2-1)*dkx
        ELSE
          zkx=(i1-1)*dkx
        END IF
!MB:
        IF(i1 >= (nx+1)) THEN
          zkx=(i1-nx2-1)*dkx
        ELSE
          zkx=(i1-1)*dkx
        END IF
!MB/
        ind=ind+1
        dervk(ind) = eye*zkx*potk(ind)
      END DO
    END DO
  END DO
  CALL rfftback(dervk,potwork)
  
ELSE IF (icoor == 2) THEN
  
  
!       y-derivative
  
  ind=0
  DO i3=1,nz2
    DO i2=1,ny2
      IF(i2 >= (ny+1)) THEN
        zky=(i2-ny2-1)*dky
      ELSE
        zky=(i2-1)*dky
      END IF
      DO i1=1,nx2
!MB:
        IF(i1 >= (nx+1)) THEN
          zkx=(i1-nx2-1)*dkx
        ELSE
          zkx=(i1-1)*dkx
        END IF
!MB/
        ind=ind+1
        dervk(ind) = eye*zky*potk(ind)
      END DO
    END DO
  END DO
  CALL rfftback(dervk,potwork)
  
ELSE IF (icoor == 3) THEN
  
  
!       z-derivative
  
  ind=0
  DO i3=1,nz2
    IF(i3 >= (nz+1)) THEN
      zkz=(i3-nz2-1)*dkz
    ELSE
      zkz=(i3-1)*dkz
    END IF
    DO i2=1,ny2
      DO i1=1,nx2
!MB:
        IF(i1 >= (nx+1)) THEN
          zkx=(i1-nx2-1)*dkx
        ELSE
          zkx=(i1-1)*dkx
        END IF
!MB/
        ind=ind+1
        dervk(ind) = eye*zkz*potk(ind)
      END DO
    END DO
  END DO
  CALL rfftback(dervk,potwork)
  
ELSE
  STOP 'Error in getGradField'
  
END IF


#else
STOP "not yet implemented for finite differences"
#endif


DO ii=1,kdfull2
  fieldto(ii) = potwork(ii)
END DO
DEALLOCATE(potk,dervk,potwork)


RETURN
END SUBROUTINE getgradfield
!------------------------------------------------------------
