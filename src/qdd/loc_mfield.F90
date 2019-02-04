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

 
!-----coul_mfield------------------------------------------------

SUBROUTINE coul_mfield(rho)

!     The Coulomb part of the mean field.

!     Input:
!      rho    = electron density
!     Output:
!      Coulomb field via 'chpcoul' in common

USE params
#if(netlib_fft|fftw_cpu)
USE coulsolv
#endif
#if(findiff|numerov)
USE coulsolv,ONLY:solv_poisson
#endif
IMPLICIT NONE

REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)

#if(raregas)
INTEGER :: ii
REAL(DP), ALLOCATABLE :: rhotmp(:)
#endif

!----------------------------------------------------------------

!     copy density on intermediate storage
!     add images for dielectric contribution
#if(raregas)
IF(idielec == 1) THEN
  
  ALLOCATE(rhotmp(2*kdfull2))
  DO ii=1,2*kdfull2
    rhotmp(ii)=rho(ii)
  END DO
  
  CALL addimage(rho,1)
  
END IF
#endif

! Coulomb of the electronic density

#if(gridfft)
CALL falr(rho,chpcoul,kdfull2)
#endif
#if(findiff|numerov)
CALL solv_poisson(rho,chpcoul,kdfull2)
#endif

!WRITE(6,'(/2a)') 'along x:  x  rho rho_image  coul'
!ind = 0
!DO jz=minz,maxz
!  DO jy=miny,maxy
!    DO jx=minx,maxx
!      ind    = 1 + ind
!      IF(jz == nzsh .AND. jy == nysh)  &
!          WRITE(6,'(1x,f6.2,3(1pg13.5))') (jx-nxsh)*dx,rhotmp(ind), &
!           rho(ind),chpcoul(ind)   
!    END DO
!  END DO
!END DO

! computation of the Coulomb potential from the electronic density
! before adjustdip or vstep for pseudodensity description

!     restore electron density, save Coulomb fields,
!     add images and solve Coulomb next round

#if(raregas)
IF(idielec == 1) THEN
  
  DO ii=1,2*kdfull2
    rho(ii)=rhotmp(ii)
  END DO
  DO ii=1,kdfull2
    rfieldtmp(ii)=chpcoul(ii)
  END DO
  
  CALL addimage(rho,0)
  
#if(gridfft)
  CALL falr(rho,chpcoul,kdfull2)
#endif
#if(findiff|numerov)
  CALL solv_fft(rho,chpcoul,dx,dy,dz)
#endif
  

!WRITE(6,'(/2a)') 'along x:  x  rho_image coul_imag'
!ind = 0
!DO jz=minz,maxz
!  DO jy=miny,maxy
!    DO jx=minx,maxx
!      ind    = 1 + ind
!      IF(jz == nzsh .AND. jy == nysh)  &
!          WRITE(6,'(1x,f6.2,3(1pg13.5))') (jx-nxsh)*dx, &
!           rho(ind),chpcoul(ind)   
!    END DO
!  END DO
!END DO
  
  DO ii=1,kdfull2
    CALL conv1to3(ii)
    IF(iindtmp(1) >= nint(xdielec/dx)+nxsh)THEN
      chpcoul(ii)=rfieldtmp(ii)
    END IF
  END DO
  DO ii=1,2*kdfull2
    rho(ii)=rhotmp(ii)
  END DO
  DEALLOCATE(rhotmp)
  
END IF
#endif

RETURN
END SUBROUTINE coul_mfield

!-----valence_step---------------------------------------------

SUBROUTINE valence_step(rho,dt,it,tdyn)

!     Propagates the substrates valence electron cloud.
!     Input:
!      rho    = electron density
!      dt     = stepsize (in case of dynamics)
!      it     = current iteration
!      tdyn   = switch to dynamic case
!     Output:
!      valence positions via common

!     This routine can still be used in statics as well as
!     in dynamics (see the switch 'tdyn'). It is presently
!     only called in dynamics.

USE params
IMPLICIT NONE

REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
REAL(DP), INTENT(IN)                     :: dt
INTEGER,INTENT(IN)                       :: it
LOGICAL, INTENT(IN)                      :: tdyn

#if(raregas)
COMPLEX(DP) :: psidummy(1)
#endif
LOGICAL,PARAMETER :: tspinprint=.true.


!---------------------------------------------------------------

#if(raregas)
IF (isurf /= 0 .AND. nc+NE+nk > 0) THEN    ! check condition ??
  
  IF(ifadiadip ==1) THEN 
!            instantaneous adjustment of substrate dipoles
    CALL adjustdip(rho,it)
  ELSE
!            dynamic propagation of substrate dipoles
    IF(tdyn) THEN
      IF(ipsptyp == 1) STOP ' VSTEP must not be used with Goedecker PsP'  ! ???
      IF(ionmdtyp==1) CALL vstep(rho,psidummy,it,dt)
      IF(ionmdtyp==2) CALL vstepv(rho,psidummy,it,dt)
    ELSE
      CALL adjustdip(rho,it)
    END IF
  END IF
  
!     update of the subgrids   ? here or outside substrate loop ?
  
  IF(.NOT.tdyn .AND. ipseudo == 1) CALL updatesubgrids  ! probably done in vstep
  CALL calcpseudo()                 !! ????????
  
  IF(nc > 0 .AND. ivdw == 1) CALL calc_frho(rho)
  
END IF
#endif

RETURN
END SUBROUTINE valence_step

!-----ldapsp_mfield------------------------------------------------

SUBROUTINE ldapsp_mfield(rho,aloc)

!     The LDA and pseudopotential part of the mean field.
!     Input:
!      rho    = electron density
!      ionic positions for PsP come in via common
!     Output:
!      aloc   = local mean field

USE params
IMPLICIT NONE

REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)



!------------------------------------------------------------------

CALL calcpseudo()                 ! update pseudo-potentials

CALL calclocal(rho,aloc)          ! LDA part of the potential

RETURN
END SUBROUTINE ldapsp_mfield

