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

 
!-----zeroforce-----------------------------------------------------------------

SUBROUTINE zeroforce(aloc,rho)

!     Corrects SIC-Slater and SIC-KLI potentials to meet the
!     zero-force condition.

USE params
USE util, ONLY:wfovlp
USE kinetic
IMPLICIT NONE


REAL(DP), INTENT(OUT)                        :: aloc(2*kdfull2)
REAL(DP), INTENT(IN)                         :: rho(2*kdfull2)

INTEGER :: i, is, ishift
REAL(DP) :: denominator, counter, xlambda, ylambda, zlambda
!       workspaces

COMPLEX(DP), ALLOCATABLE :: potk(:),dervk(:)     ! for Fourier transformed potentials
REAL(DP), ALLOCATABLE :: potwork(:)

!      equivalence (potwork(1),w1(1))
!      equivalence (potk(1),w2(1))              ! occupies also w3

!-----------------------------------------------------------------------------

IF(.NOT.ALLOCATED(akv)) STOP "ZEROFORCE requires FFT"
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
!#if(netlib_fft|fftw_cpu)
    CALL rftf(aloc(1+ishift),potk)
    
!       Laplacian and integral
    
    DO  i=1,nxyz
      dervk(i) = potk(i)*akv(i)
    END DO
    CALL rfftback(dervk,potwork)
!#endif
    denominator = wfovlp(rho(1+ishift : kdfull2+ishift),potwork)
    
!       x-derivative      
!#if(netlib_fft|fftw_cpu)
!    ind=0
!    DO i3=1,nz2
!      DO i2=1,ny2
!        DO i1=1,nx2
!          IF(i1 >= (nx+1)) THEN
!            zkx=(i1-nx2-1)*dkx
!          ELSE
!            zkx=(i1-1)*dkx
!          END IF
!!MB:
!          IF(i1 >= (nx+1)) THEN
!            zkx=(i1-nx2-1)*dkx
!          ELSE
!            zkx=(i1-1)*dkx
!          END IF
!!MB/
!          ind=ind+1
          dervk = -akx*potk
!        END DO
!      END DO
!    END DO
    CALL rfftback(dervk,potwork)
!#endif
    counter = wfovlp(rho(1+ishift : kdfull2+ishift),potwork)
    xlambda = counter/denominator
    DO i=1,nxyz
      aloc(i+ishift) = aloc(i+ishift)-xlambda*potwork(i)
    END DO
    
!       y-derivative
!#if(netlib_fft|fftw_cpu)
!    ind=0
!    DO i3=1,nz2
!      DO i2=1,ny2
!        IF(i2 >= (ny+1)) THEN
!          zky=(i2-ny2-1)*dky
!        ELSE
!          zky=(i2-1)*dky
!        END IF
!        DO i1=1,nx2
!!MB:
!          IF(i1 >= (nx+1)) THEN
!            zkx=(i1-nx2-1)*dkx
!          ELSE
!            zkx=(i1-1)*dkx
!          END IF
!!MB/
!          ind=ind+1
          dervk = -aky*potk
!        END DO
!      END DO
!    END DO
    CALL rfftback(dervk,potwork)
!#endif
    counter = wfovlp(rho(1+ishift : kdfull2+ishift),potwork)
    ylambda = counter/denominator
    DO i=1,nxyz
      aloc(i+ishift) = aloc(i+ishift)-ylambda*potwork(i)
    END DO
    
!       z-derivative
!#if(netlib_fft|fftw_cpu)
!    ind=0
!    DO i3=1,nz2
!      IF(i3 >= (nz+1)) THEN
!        zkz=(i3-nz2-1)*dkz
!      ELSE
!        zkz=(i3-1)*dkz
!      END IF
!      DO i2=1,ny2
!        DO i1=1,nx2
!!MB:
!          IF(i1 >= (nx+1)) THEN
!            zkx=(i1-nx2-1)*dkx
!          ELSE
!            zkx=(i1-1)*dkx
!          END IF
!!MB/
!          ind=ind+1
!          dervk(ind) = eye*zkz*potk(ind)
          dervk = -akz*potk
!        END DO
!      END DO
!    END DO
    CALL rfftback(dervk,potwork)
!#endif
    counter = wfovlp(rho(1+ishift : kdfull2+ishift),potwork)
    zlambda = counter/denominator
    DO i=1,nxyz
      aloc(i+ishift) = aloc(i+ishift)-zlambda*potwork(i)
    END DO
    
  END DO
  DEALLOCATE(potk,dervk,potwork)
  
  RETURN
END SUBROUTINE zeroforce

!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE checkzeroforce(rho,aloc)
!------------------------------------------------------------
USE params
USE util, ONLY:wfovlp
USE kinetic
IMPLICIT NONE

REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
REAL(DP), INTENT(IN)                     :: aloc(2*kdfull2)

INTEGER :: is, ishift
REAL(DP) :: zforcex, zforcey, zforcez

!       workspaces

COMPLEX(DP), ALLOCATABLE :: potk(:),dervk(:)     ! for Fourier transformed potentials
REAL(DP), ALLOCATABLE :: potwork(:)

IF(.NOT.ALLOCATED(akv)) STOP "ZEROFORCE requires FFT"

ALLOCATE(potk(kdfull2),dervk(kdfull2),potwork(kdfull2))

DO is=1,numspin
    ishift = (is-1)*nxyz
    
!       Fourier transformation
    
!#if(netlib_fft|fftw_cpu)
    CALL rftf(aloc(1+ishift),potk)
!#endif
    
    
!       x-derivative
!#if(netlib_fft|fftw_cpu)
!    ind=0
!    DO i3=1,nz2
!      DO i2=1,ny2
!        DO i1=1,nx2
!          IF(i1 >= (nx+1)) THEN
!            zkx=(i1-nx2-1)*dkx
!          ELSE
!            zkx=(i1-1)*dkx
!          END IF
!!MB:
!          IF(i1 >= (nx+1)) THEN
!            zkx=(i1-nx2-1)*dkx
!          ELSE
!            zkx=(i1-1)*dkx
!          END IF
!!MB/
!          ind=ind+1
!          dervk(ind) = eye*zkx*potk(ind)
          dervk = -akx*potk
!        END DO
!      END DO
!    END DO
    CALL rfftback(dervk,potwork)
!#endif
    zforcex = wfovlp(rho(1+ishift : kdfull2+ishift),potwork)
    
    
!       y-derivative
!#if(netlib_fft|fftw_cpu)
!    ind=0
!    DO i3=1,nz2
!      DO i2=1,ny2
!        IF(i2 >= (ny+1)) THEN
!          zky=(i2-ny2-1)*dky
!        ELSE
!          zky=(i2-1)*dky
!        END IF
!        DO i1=1,nx2
!!MB:
!          IF(i1 >= (nx+1)) THEN
!            zkx=(i1-nx2-1)*dkx
!          ELSE
!            zkx=(i1-1)*dkx
!          END IF
!!MB/
!          ind=ind+1
!          dervk(ind) = eye*zky*potk(ind)
          dervk = -aky*potk
!        END DO
!      END DO
!    END DO
    CALL rfftback(dervk,potwork)
!#endif
    zforcey = wfovlp(rho(1+ishift : kdfull2+ishift),potwork)
    
    
!       z-derivative
!#if(netlib_fft|fftw_cpu)
!    ind=0
!    DO i3=1,nz2
!      IF(i3 >= (nz+1)) THEN
!        zkz=(i3-nz2-1)*dkz
!      ELSE
!        zkz=(i3-1)*dkz
!      END IF
!      DO i2=1,ny2
!        DO i1=1,nx2
!!MB:
!          IF(i1 >= (nx+1)) THEN
!            zkx=(i1-nx2-1)*dkx
!          ELSE
!            zkx=(i1-1)*dkx
!          END IF
!!MB/
!          ind=ind+1
!          dervk(ind) = eye*zkz*potk(ind)
          dervk = -akz*potk
!        END DO
!      END DO
!    END DO
    CALL rfftback(dervk,potwork)
!#endif
    zforcez = wfovlp(rho(1+ishift : kdfull2+ishift),potwork)
    
    
    
  END DO
  
  DEALLOCATE(potk,dervk,potwork)
  WRITE(6,'(a,3e17.7)') 'Checking Zero-Force Theorem: ',  &
      zforcex,zforcey,zforcez
  
  WRITE(500,'(1f17.5,3e17.7)') tfs,zforcex,zforcey,zforcez

  
  RETURN
END SUBROUTINE checkzeroforce
