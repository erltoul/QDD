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




!     ******************************

#ifdef REALSWITCH
SUBROUTINE calc_sicr(rho,aloc,q0)
#else
SUBROUTINE calc_sic(rho,aloc,q0)
#endif

!     ******************************

!     switchboard for the various brands of SIC

USE params
USE kinetic
IMPLICIT NONE
#ifdef REALSWITCH
REAL(DP), INTENT(IN) :: q0(kdfull2,kstate)
#else
COMPLEX(DP), INTENT(IN) :: q0(kdfull2,kstate)
#endif
!REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2),rho(2*kdfull2)
REAL(DP), INTENT(IN OUT) :: aloc(*),rho(*)

!------------------------------------------------------------




IF (ifreezekspot == 1 .AND. tfs > 0) RETURN


#ifdef REALSWITCH
SELECT CASE(ifsicp)
  CASE(2)
    CALL calc_adsicr(rho,aloc)
  CASE(3)
    CALL calc_slaterr(rho,aloc,q0)
  CASE(4)
    CALL calc_sicklir(rho,aloc,q0)
  CASE DEFAULT
    RETURN
END SELECT
#else
SELECT CASE(ifsicp)
  CASE(2)
    CALL calc_adsic(rho,aloc)
  CASE(3)
    CALL calc_slater(rho,aloc,q0)
  CASE(4)
    CALL calc_sickli(rho,aloc,q0)
  CASE DEFAULT
    RETURN
END SELECT
#endif
RETURN
#ifdef REALSWITCH
END SUBROUTINE calc_sicr
#else
END SUBROUTINE calc_sic
#endif







!     ******************************

#ifdef REALSWITCH
SUBROUTINE calc_adsicr(rho,aloc)
#else
SUBROUTINE calc_adsic(rho,aloc)
#endif

!     ******************************

!     computes local part of Hamiltonian
!     'time' <= 0  signal static iteration i.e. without laser field


USE params
USE kinetic
#if(netlib_fft|fftw_cpu)
USE coulsolv
#endif
#if(findiff|numerov)
USE coulsolv,ONLY:solv_poisson
#endif

IMPLICIT NONE

REAL(DP),PARAMETER :: sgnkli=-1D0

#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

REAL(DP) :: aloc(2*kdfull2),rho(2*kdfull2)
REAL(DP),DIMENSION(:),ALLOCATABLE :: rhosp,chpdftsp,coulsum,couldif
REAL(DP),DIMENSION(:),ALLOCATABLE :: rho1,rho2
REAL(DP),DIMENSION(:),ALLOCATABLE :: rhospu,rhospd,chpdftspu,chpdftspd
INTEGER :: ind, idx, npartdw, npartup, npartto, size
INTEGER :: indpri
REAL(DP) :: enrearsave, enerpwsave, enrear1, enrear2, enpw1, enpw2, encadd
REAL(DP) :: couldw, coulup, fac, facdw, facup, facdwh, facuph, factotal
!LOGICAL,PARAMETER :: ttest=.TRUE.
!INTEGER,EXTERNAL :: getnearestgridpoint
#if(parayes)
REAL(DP),DIMENSION(:),ALLOCATABLE :: rhonod1,rhonod2,rhotp
REAL(DP),DIMENSION(:),ALLOCATABLE :: tp1,tp2
REAL(DP),DIMENSION(:),ALLOCATABLE :: aloc1,aloc2
INTEGER :: nod
#endif

!-------------------------------------------------------------------------

IF(ifsicp /= 2) THEN
  STOP ' CALC_SIC called with wrong option IFSICP'
END IF

CALL act_part_num(npartup,npartdw,npartto)

#if(parayes)
CALL mpi_comm_rank(mpi_comm_world,nod,mpi_ierror)
size=lengnod(nod+1)
ALLOCATE(rhonod1(size))
ALLOCATE(rhonod2(size))
CALL pi_scatterv(rho,nxyz,rhonod1,size,mpi_ierror)
CALL pi_scatterv(rho(nxyz+1),nxyz,rhonod2,size,mpi_ierror)
ALLOCATE(tp1(size))
ALLOCATE(tp2(size))
#else
size=nxyz
#endif



enrearsave=enrear
enerpwsave=enerpw
enrear1=0D0
enrear2=0D0
enpw1  = 0D0
enpw2  = 0D0
encadd = 0D0
IF(numspin==2) THEN

  ALLOCATE(rhospu(2*kdfull2),rhospd(2*kdfull2),  &
                  chpdftspu(2*kdfull2),chpdftspd(2*kdfull2))

!     averaged spinup and spindown density

  IF(npartup > 0) THEN
    facuph = 0.5D0/npartup
  
!    DO ind=1,nxyz
    DO ind=1,size
#if(parayes)
      tp1(ind)=rhonod1(ind)*(1D0+rhonod2(ind))*facuph
      tp2(ind)=1D0
#else
      rhospu(ind)=rho(ind)*(1D0+rho(ind+nxyz))*facuph
      rhospu(ind+nxyz)=1D0
#endif
    END DO
  END IF
#if(parayes)
  CALL pi_allgatherv(tp1,size,rhospu,nxyz,mpi_ierror)
  CALL pi_allgatherv(tp2,size,rhospu(nxyz+1),nxyz,mpi_ierror)
#endif

  IF(npartdw > 0) THEN
    facdwh = 0.5D0/npartdw
  
!    DO ind=1,nxyz
    DO ind=1,size
#if(parayes)
      tp1(ind)=rhonod1(ind)*(1D0-rhonod2(ind))*facdwh
      tp2(ind)=-1D0
#else
      rhospd(ind)=rho(ind)*(1D0-rho(ind+nxyz))*facdwh
      rhospd(ind+nxyz)=-1D0
#endif
    END DO
  END IF
#if(parayes)
  CALL pi_allgatherv(tp1,size,rhospd,nxyz,mpi_ierror)
  CALL pi_allgatherv(tp2,size,rhospd(nxyz+1),nxyz,mpi_ierror)
#endif

!     DFT for averaged spinup and spindown density

  IF(npartup > 0) THEN
    CALL calc_lda(rhospu,chpdftspu)
    enrear1=enrear
  END IF
  IF(npartdw > 0) THEN
    CALL calc_lda(rhospd,chpdftspd)
    enrear2=enrear
  END IF

  enrear   = enrearsave-enrear1*npartup-enrear2*npartdw
  WRITE(6,*) ' enrear.s=',enrearsave, enrear1,enrear2,enrear

#if(parayes)
ALLOCATE(aloc1(size))
ALLOCATE(aloc2(size))
CALL pi_scatterv(chpdftspu,nxyz,tp1,size,mpi_ierror)
CALL pi_scatterv(chpdftspd(nxyz+1),nxyz,tp2,size,mpi_ierror)
CALL pi_scatterv(aloc,nxyz,aloc1,size,mpi_ierror)
CALL pi_scatterv(aloc(nxyz+1),nxyz,aloc2,size,mpi_ierror)
#endif

  IF(npartup > 0) THEN
!    DO ind=1,nxyz
    DO ind=1,size
#if(parayes)
      aloc1(ind) = aloc1(ind) - tp1(ind)
#else
      aloc(ind) = aloc(ind) - chpdftspu(ind)
#endif
    END DO
#if(parayes)
    CALL pi_allgatherv(aloc1,size,aloc,nxyz,mpi_ierror)
#endif
  END IF



  IF(npartdw > 0) THEN
#if(parayes)
    DO idx=1,size
      aloc2(idx) = aloc2(idx) - tp2(idx)
    END DO
#else
    DO idx=nxyz+1,2*nxyz
      aloc(idx) = aloc(idx) - chpdftspd(idx)
    END DO
#endif
#if(parayes)
    CALL pi_allgatherv(aloc2,size,aloc(nxyz+1),nxyz,mpi_ierror)
#endif
  END IF


  DEALLOCATE(rhospu,rhospd,chpdftspu,chpdftspd)

ELSE

  ALLOCATE(chpdftsp(2*kdfull2))
  ALLOCATE(rhosp(2*kdfull2))

  factotal=1D0/npartto

#if(parayes)
  ALLOCATE(aloc1(size))
#endif

!  DO ind=1,nxyz
  DO ind=1,size
#if(parayes)
    tp1(ind)=rhonod1(ind)*factotal
    tp2(ind)=1D0
#else
    rhosp(ind)=rho(ind)*factotal
    rhosp(ind+nxyz)=1D0
#endif
  END DO

#if(parayes)
CALL pi_allgatherv(tp1,size,rhosp,nxyz,mpi_ierror)
CALL pi_allgatherv(tp2,size,rhosp(nxyz+1),nxyz,mpi_ierror)
#endif

!     DFT for averaged s.p. state

  CALL calc_lda(rhosp,chpdftsp)

!     subtract from pure LDA potential

#if(parayes)
CALL pi_scatterv(chpdftsp,nxyz,tp1,size,mpi_ierror)
CALL pi_scatterv(aloc,nxyz,aloc1,size,mpi_ierror)
#endif
!  DO ind=1,nxyz
  DO ind=1,size
#if(parayes)
    aloc1(ind) = aloc1(ind) - tp1(ind)
#else
    aloc(ind) = aloc(ind) - chpdftsp(ind)
#endif
  END DO
#if(parayes)
  CALL pi_allgatherv(aloc1,size,aloc,nxyz,mpi_ierror)
#endif

  enrear   = enrearsave-enrear*npartto
  DEALLOCATE(chpdftsp)
  DEALLOCATE(rhosp)

END IF


!     correct Coulomb potential by 1/N

IF(numspin==2) THEN
  ALLOCATE(coulsum(kdfull2))
  ALLOCATE(couldif(kdfull2))
  ALLOCATE(rho1(kdfull2))
  ALLOCATE(rho2(kdfull2))

!       compute Coulomb


!  DO ind=1,nxyz
  DO ind=1,size
#if(parayes)
    tp2(ind)=  rhonod1(ind)
    tp1(ind)=  rhonod1(ind)*rhonod2(ind)
#else
    rho2(ind)=  rho(ind)
    rho1(ind)=  rho(ind)*rho(ind+nxyz)
#endif
  END DO
#if(parayes)
    CALL pi_allgatherv(tp1,size,rho1,nxyz,mpi_ierror)
    CALL pi_allgatherv(tp2,size,rho2,nxyz,mpi_ierror)
#endif
!  IF(ttest) THEN
!    indpri=getnearestgridpoint((maxx-nxsh)*dx,0D0,0D0)
!    WRITE(6,'(a,3(1pg13.5))') ' ADSIC before:',aloc(indpri),&
!        SUM(rho2(1:size)),SUM(rho1(1:size))
!  END IF

#if(gridfft)
  CALL falr(rho1,couldif,kdfull2)
  IF(idielec == 0 .AND. nion2 > 0) THEN
    coulsum = chpcoul
  ELSE
    CALL falr(rho2,coulsum,kdfull2)
  END IF
#endif
#if(findiff|numerov)
  CALL solv_poisson(rho1(1),couldif,kdfull2)
  CALL solv_poisson(rho2(1),coulsum,kdfull2)
#endif
  facup = 1D0/npartup

#if(parayes)
  CALL pi_scatterv(coulsum,nxyz,tp1,size,mpi_ierror)
  CALL pi_scatterv(couldif,nxyz,tp2,size,mpi_ierror)
  CALL pi_scatterv(aloc,nxyz,aloc1,size,mpi_ierror)
#endif
!  DO ind=1,nxyz
  DO ind=1,size
#if(parayes)
    coulup    = 0.5D0*(tp1(ind)+tp2(ind))
    aloc1(ind) = aloc1(ind) - coulup*facup
#else
    coulup    = 0.5D0*(coulsum(ind)+couldif(ind))
    aloc(ind) = aloc(ind) - coulup*facup
#endif
!  IF(ttest .AND. ind==indpri) &
!       WRITE(6,'(a,2(1pg13.5))') ' ADSIC up:',aloc(indpri),coulup*facup
  END DO
#if(parayes)
  CALL pi_allgatherv(aloc1,size,aloc,nxyz,mpi_ierror)
#endif


  IF (npartdw > 0) THEN
  
#if(parayes)
    CALL pi_scatterv(aloc(nxyz+1),nxyz,aloc2,size,mpi_ierror)
#endif
    facdw = 1D0/npartdw
!    DO ind=1,nxyz
    DO ind=1,size
#if(parayes)
      couldw    = 0.5D0*(tp1(ind)-tp2(ind))
      aloc2(ind) = aloc2(ind) - facdw*couldw
#else
      couldw    = 0.5D0*(coulsum(ind)-couldif(ind))
      idx       = ind + nxyz
      aloc(idx) = aloc(idx) - facdw*couldw
#endif
!      IF(ttest .AND. ind==indpri) &
!       WRITE(6,'(a,2(1pg13.5))') ' ADSIC dw:',aloc(indpri),couldw*facdw
    END DO
#if(parayes)
    CALL pi_allgatherv(aloc2,size,aloc(nxyz+1),nxyz,mpi_ierror)
    DEALLOCATE(aloc1,aloc2,tp1,tp2)
#endif
  
  END IF

  DEALLOCATE(coulsum)
  DEALLOCATE(couldif)
  DEALLOCATE(rho1)
  DEALLOCATE(rho2)

ELSE

!     recalculate Coulomb part for jellium

#if(gridfft)
  IF(nion2 == 0) THEN
    CALL falr(rho(1:kdfull2),chpcoul,kdfull2)
  END IF
#endif
#if(findiff|numerov)
  STOP ' ADSIC not yet ready for nospin and finite diff.'
#endif

  fac = 1D0/npartto
#if(parayes)
  CALL pi_scatterv(chpcoul,nxyz,tp1,size,mpi_ierror)
  CALL pi_scatterv(aloc,nxyz,aloc1,size,mpi_ierror)
!  WRITE(*,*) 'scattering of chpcoul and aloc done'
#endif
!  DO ind=1,nxyz
  DO ind=1,size
#if(parayes)
    aloc1(ind) = aloc1(ind) - fac*tp1(ind)
#else
    aloc(ind) = aloc(ind) - fac*chpcoul(ind)
#endif
  END DO
#if(parayes)
  CALL pi_allgatherv(aloc1,size,aloc,nxyz,mpi_ierror)
  DEALLOCATE(aloc1,tp1,tp2)
#endif

END IF


RETURN
#ifdef REALSWITCH
END SUBROUTINE calc_adsicr
#else
END SUBROUTINE calc_adsic
#endif





!     ******************************

#ifdef REALSWITCH
SUBROUTINE calc_slaterr(rho,aloc,q0)
#else
SUBROUTINE calc_slater(rho,aloc,q0)
#endif

!     ******************************

!     computes local part of Hamiltonian
!     'time' <= 0  signal static iteration i.e. without laser field


USE params
USE util, ONLY:prifld2
USE kinetic
IMPLICIT NONE

REAL(DP),PARAMETER :: sgnkli=-1D0
#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif
#ifdef REALSWITCH
REAL(DP) :: q0(kdfull2,kstate)
#else
COMPLEX(DP) :: q0(kdfull2,kstate)
#endif
REAL(DP) :: aloc(2*kdfull2),rho(2*kdfull2)

REAL(DP),DIMENSION(:),ALLOCATABLE :: usicsp,rhosp
REAL(DP),DIMENSION(:),ALLOCATABLE :: rhospu,rhospd
LOGICAL :: testprint
INTEGER :: ind, idx, ishift, nb, npartdw, npartup, npartto
REAL(DP) :: enrearsave, enerpwsave, enrear1, enrear2, enpw1, enpw2, encadd, reldw, relup

DATA testprint/.false./

!-------------------------------------------------------------------

IF(numspin.NE.2) STOP ' SIC-Slater requires full spin'
IF(ifsicp /= 3) STOP ' CALC_SLATER called with wrong option IFSICP'

CALL act_part_num(npartup,npartdw,npartto)


IF(testprint .AND. ifsicp == 3) THEN
  OPEN(11,FILE='testslater')
  CALL prifld2(11,aloc,'pot before u')
  CALL prifld2(11,aloc(1+nxyz),'pot before d')
END IF

! IF(usew1) STOP ' in CALC_SIC: workspace W1 already active '
! usew1 = .true.
ALLOCATE(rhosp(2*kdfull2))
ALLOCATE(usicsp(2*kdfull2))
ALLOCATE(rhospu(2*kdfull2),rhospd(2*kdfull2))

enrearsave=enrear
enerpwsave=enerpw
enrear1=0D0
enrear2=0D0
enpw1  = 0D0
enpw2  = 0D0
encadd = 0D0

!     recombine to spin up and spin down densities

IF(npartup > 0) THEN
  DO ind=1,nxyz
    rhospu(ind)=MAX(rho(ind)*(1D0+rho(ind+nxyz))*0.5D0,small)
  END DO
END IF
IF(npartdw > 0) THEN
  DO ind=1,nxyz
    rhospd(ind)=MAX(rho(ind)*(1D0-rho(ind+nxyz))*0.5D0,small)
  END DO
END IF

!     loop over s.p. states, compute and accumulate SIC-Slater potential

DO nb=1,nstate
  IF(occup(nb) > small) THEN
    ishift = (ispin(nrel2abs(nb))-1)*nxyz ! store spin=2 in upper block
#ifdef REALSWITCH
    CALL calc_sicspr(rhosp,usicsp,q0(1,nb),nb)
#else
    CALL calc_sicsp(rhosp,usicsp,q0(1,nb),nb)
#endif
    
    IF (ispin(nrel2abs(nb)) == 1) THEN
      enrear1=enrear1+enrear*occup(nb)
      enpw1=enpw1+enerpw*occup(nb)
    ELSE
      enrear2=enrear2+enrear*occup(nb)
      IF(directenergy) THEN
        enpw2=enpw2+enerpw*occup(nb)
      END IF
    END IF
    IF(directenergy) THEN
      encadd=encadd+encoulsp*occup(nb)
    END IF
    
    IF (ispin(nrel2abs(nb)) == 1) THEN
      DO ind=1,nxyz
        relup = occup(nb)*rhosp(ind)/MAX(rhospu(ind),small)
        aloc(ind)=aloc(ind)-relup*usicsp(ind)
      END DO
    ELSE IF (ispin(nrel2abs(nb)) == 2) THEN
      DO ind=1,nxyz
        idx = ind+nxyz
        reldw = occup(nb)*rhosp(ind)/MAX(rhospd(ind),small)
        aloc(idx)=aloc(idx)-reldw*usicsp(idx)
      END DO
    END IF
  END IF
END DO
encadd=encadd/2D0
enrear   = enrearsave-enrear1-enrear2
IF(directenergy) THEN
  enerpw   = enerpwsave-enpw1-enpw2-encadd
END IF
! usew1 = .false.
DEALLOCATE(rhosp)
DEALLOCATE(usicsp)
DEALLOCATE(rhospu,rhospd)


RETURN
#ifdef REALSWITCH
END SUBROUTINE calc_slaterr
#else
END SUBROUTINE calc_slater
#endif





!     ******************************

#ifdef REALSWITCH
SUBROUTINE calc_sicklir(rho,aloc,q0)
#else
SUBROUTINE calc_sickli(rho,aloc,q0)
#endif

!     ******************************

!     computes local part of Hamiltonian
!     'time' <= 0  signal static iteration i.e. without laser field


USE params
USE util, ONLY:prifld
USE kinetic
IMPLICIT NONE
REAL(DP),PARAMETER :: sgnkli=-1D0

#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

#ifdef REALSWITCH
REAL(DP), INTENT(IN) :: q0(kdfull2,kstate)
#else
COMPLEX(DP), INTENT(IN) :: q0(kdfull2,kstate)
#endif
REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2),rho(2*kdfull2)

REAL(DP),DIMENSION(:),ALLOCATABLE :: usicsp,rhosp,rho1,rho2
REAL(DP),DIMENSION(:),ALLOCATABLE :: rhospu,rhospd
REAL(DP),ALLOCATABLE :: rhokli(:,:),uslater(:),ukli(:)
LOGICAL,PARAMETER :: testprint=.false.
INTEGER :: itmaxkli=200
LOGICAL:: converged=.false.

INTEGER :: ind, idx, ishift, itkli, nb, nfdw, nfup, npartup, npartdw, npartto
REAL(DP) :: addn, addo, correct, efermup, efermdw
REAL(DP) :: enrearsave, enerpwsave, enrear1, enrear2, enpw1, enpw2, encadd, reldw, relup
REAL(DP) :: acc, avuslatdw, avuslatup, sumsldw, sumslup, sumdw, sumup
!--------------------------------------------------------------------

IF(numspin.NE.2) STOP ' SIC-KLI requires full spin code'
IF(ifsicp /= 4) STOP ' CALC_SICKLI called with wrong option IFSICP'

ALLOCATE(rhokli(kdfull2,kstate),uslater(2*kdfull2),ukli(2*kdfull2))

CALL act_part_num(npartup,npartdw,npartto)

!     determine Fermi levels

IF(temp /= 0D0) STOP 'KLI not yet compatible with temperature'
nfup = 0
nfdw = 0
IF(sgnkli == -1D0) THEN
  efermup = -1000D0
  efermdw = -1000D0
  DO nb=1,nstate
    IF (ispin(nrel2abs(nb)) == 1) THEN      ! M : 1=up 2=down
      IF(amoy(nb) > efermup .AND. occup(nb) > 0.5D0) THEN
        nfup = nb
        efermup = amoy(nb)
      END IF
    ELSE
      IF(amoy(nb) > efermdw .AND. occup(nb) > 0.5D0) THEN
        nfdw = nb
        efermdw = amoy(nb)
      END IF
    END IF
  END DO
  IF(nfup == 0) STOP ' Fermi state spin-up not found'
  IF(nfdw == 0) STOP ' Fermi state spin-down not found'
END IF


!     check workspace

! IF(usew1) STOP ' in CALC_SIC: workspace W1 already active '
! usew1 = .true.

enrearsave=enrear
enerpwsave=enerpw
enrear1=0D0
enrear2=0D0
enpw1  = 0D0
enpw2  = 0D0
encadd = 0D0
ALLOCATE(rhospu(2*kdfull2),rhospd(2*kdfull2))

!     spin up and spin down densities,

IF(npartup > 0) THEN
  DO ind=1,nxyz
    rhospu(ind)=MAX(rho(ind)*(1D0+rho(ind+nxyz))*0.5D0,small)
  END DO
END IF
IF(npartdw > 0) THEN
  DO ind=1,nxyz
    rhospd(ind)=MAX(rho(ind)*(1D0-rho(ind+nxyz))*0.5D0,small)
  END DO
END IF

!     reset accumulator

DO ind=1,2*nxyz
  uslater(ind)=0D0
END DO
sumslup = 0D0
sumsldw = 0D0

!     loop over s.p. states
ALLOCATE(rhosp(2*kdfull2))
ALLOCATE(usicsp(2*kdfull2))

DO nb=1,nstate
  IF(occup(nb) > small) THEN
    ishift = (ispin(nrel2abs(nb))-1)*nxyz ! store spin=2 in upper block
#ifdef REALSWITCH
    CALL calc_sicspr(rhosp,usicsp,q0(1,nb),nb)
#else
    CALL calc_sicsp(rhosp,usicsp,q0(1,nb),nb)
#endif
    
    IF (ispin(nrel2abs(nb)) == 1) THEN
      enrear1=enrear1+enrear*occup(nb)
      enpw1=enpw1+enerpw*occup(nb)
    ELSE
      enrear2=enrear2+enrear*occup(nb)
      IF(directenergy) THEN
        enpw2=enpw2+enerpw*occup(nb)
      END IF
    END IF
    IF(directenergy) THEN
      encadd=encadd+encoulsp*occup(nb)
    END IF
    IF (ispin(nrel2abs(nb)) == 1) THEN
      acc = 0D0
      DO ind=1,nxyz
        relup = rhosp(ind)/MAX(rhospu(ind),small)
        rhokli(ind,nb) = relup
        uslater(ind) = uslater(ind) + relup*usicsp(ind)
        acc = rhosp(ind)*usicsp(ind) + acc
      END DO
      IF(nb == nfup) avuslatup=acc*dvol
      acc = acc*dvol*sgnkli
      DO ind=1,nxyz
!old              uslater(ind) = uslater(ind)-acc*rhokli(ind,nb)
        uslater(ind) = uslater(ind)+acc*rhokli(ind,nb)
      END DO
      sumslup = sumslup+acc
    ELSE IF (ispin(nrel2abs(nb)) == 2) THEN
      acc = 0D0
      DO ind=1,nxyz
        idx = ind+nxyz
        reldw = rhosp(ind)/MAX(rhospd(ind),small)
        rhokli(ind,nb) = reldw
        uslater(idx) = uslater(idx) +reldw*usicsp(idx)
        acc = rhosp(ind)*usicsp(idx) + acc
      END DO
      IF(nb == nfdw) avuslatdw=acc*dvol
      acc = acc*dvol*sgnkli
      DO ind=1,nxyz
        idx = ind+nxyz
!old            uslater(idx) = uslater(idx)-acc*rhokli(ind,nb)
        uslater(idx) = uslater(idx)+acc*rhokli(ind,nb)
      END DO
      sumsldw = sumsldw+acc
    END IF
  END IF
END DO
encadd=encadd/2D0
enrear   = enrearsave-enrear1-enrear2
IF(directenergy) THEN
  enerpw   = enerpwsave-enpw1-enpw2-encadd
END IF

DEALLOCATE(rhosp)
DEALLOCATE(usicsp)
!ALLOCATE(rhospu(2*kdfull2),rhospd(2*kdfull2))

!     the first loop finished -  Slater potential now on 'uslater'


!     the KLI iteration

IF(testprint) WRITE(6,'(a,2(1pg12.4))') ' sumslp=',sumslup,sumsldw
sumslup=0D0
sumsldw=0D0
DO ind=1,2*nxyz
  ukli(ind)=uslater(ind)
END DO

ALLOCATE(rho1(kdfull2))
ALLOCATE(rho2(kdfull2))

DO itkli=1,itmaxkli
  DO ind=1,nxyz
    rho1(ind) = 0D0
!    w2(ind) = 0D0
  END DO
  DO nb=1,nstate
    IF(occup(nb) > small) THEN
      ishift = (ispin(nrel2abs(nb))-1)*nxyz
      acc = 0D0
      IF(ishift <= 0) THEN
        DO ind=1,nxyz
          acc = ukli(ind)*rhokli(ind,nb)*rhospu(ind) + acc
        END DO
        acc = acc*dvol*sgnkli
        DO ind=1,nxyz
          rho1(ind) = rhokli(ind,nb)*acc + rho1(ind)
        END DO
      ELSE
        DO ind=1,nxyz
          idx = ind+ishift
          acc = ukli(idx)*rhokli(ind,nb)*rhospd(ind) + acc
        END DO
        acc = acc*dvol*sgnkli
!        DO ind=1,nxyz
!          w2(ind) = rhokli(ind,nb)*acc + w2(ind)
!        END DO
      END IF
    END IF
  END DO
  addn = 0.3D0
  addo = 1D0-addn
  sumup = 0D0
  sumdw = 0D0
  DO ind=1,nxyz
    idx = ind+nxyz
!old          ukli(ind)=uslater(ind)+w1(ind)
!old          ukli(idx)=uslater(idx)+w2(ind)
    ukli(ind)=addo*ukli(ind)+addn*(uslater(ind)-rho1(ind))
    sumup = (ukli(ind)-(uslater(ind)-rho1(ind)))**2 *rhospu(ind) + sumup
    ukli(idx)=addo*ukli(idx)+addn*(uslater(idx)-rho2(ind))
    sumdw = (ukli(idx)-(uslater(idx)-rho2(ind)))**2 *rhospd(ind) + sumdw
!old          sumup = rhospu(ind)*w1(ind) +sumup
!old          sumdw = rhospd(ind)*w2(ind) +sumdw
  END DO
  sumup = dvol*sumup
  sumdw = dvol*sumdw
  IF(sgnkli == -1D0) THEN        ! asymptotic correction
    acc = 0D0
    DO ind=1,nxyz
      acc = rhospu(ind)*rhokli(ind,nfup)*ukli(ind)+acc
    END DO
    acc = acc*dvol
    correct = acc-avuslatup
    DO ind=1,nxyz
      ukli(ind) = ukli(ind)-correct
    END DO
    acc = 0D0
    DO ind=1,nxyz
      idx = ind+nxyz
      acc = rhospd(ind)*rhokli(ind,nfdw)*ukli(idx)+acc
    END DO
    acc = acc*dvol
    correct = acc-avuslatdw
    DO ind=nxyz+1,2*nxyz
      ukli(ind) = ukli(ind)-correct
    END DO
  END IF
  IF(testprint) WRITE(6,'(a,i4,3(1pg12.4))') ' itkli,sumup,sumdw=',itkli,  &
      sumup,sumdw,epsoro
!old   &      sumup-sumslup,sumdw-sumsldw,epsoro
!old        if(abs(sumup-sumslup)+abs(sumdw-sumsldw).lt.epsoro) goto 99
  IF(sumup+sumdw < epsoro**2*2D0)THEN
    converged=.true.
    EXIT
  END IF
  sumslup = sumup
  sumsldw = sumdw
END DO

IF(.NOT. converged) WRITE(6,'(a,2(1pg12.4))')  &
      ' KLI not converged: errors=',sumup-sumslup,sumdw-sumsldw
      
WRITE(6,'(a,i5,4(1pg12.4))')  &
    ' KLI itkli,sums=',itkli,sumup,sumslup,sumdw,sumsldw

DEALLOCATE(rho1)
DEALLOCATE(rho2)
!     KLI iteration terminated - modify local potential

DO ind=1,nxyz
  idx = ind+nxyz
  aloc(ind)=aloc(ind)-ukli(ind)
  aloc(idx)=aloc(idx)-ukli(idx)
END DO
IF(testprint) CALL prifld(ukli,' SIC-KLI  ')

DEALLOCATE(rhokli,uslater,ukli)

RETURN
#ifdef REALSWITCH
END SUBROUTINE calc_sicklir
#else
END SUBROUTINE calc_sickli
#endif


#ifdef REALSWITCH

!     ******************************

SUBROUTINE act_part_num(npartup,npartdw,npartto)

!     ******************************

!     computes actual particle number for spin-up and spin-down.

USE params
USE kinetic
IMPLICIT NONE
#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
REAL(DP) :: partdwp, partupp
#endif


INTEGER,INTENT(OUT) :: npartup
INTEGER,INTENT(OUT) :: npartdw
INTEGER,INTENT(OUT) :: npartto

INTEGER :: nb
REAL(DP) :: partup,partdw

!     compute number of states

partup=0D0
partdw=0D0
DO nb=1,nstate
  IF (ispin(nrel2abs(nb)) == 1) THEN      ! M : 1=up 2=down
    partup = occup(nb) + partup
  ELSE
    partdw = occup(nb) + partdw
  END IF
END DO

#if(parayes)
CALL mpi_barrier (mpi_comm_world, mpi_ierror)
CALL mpi_allreduce(partup,partupp,1,mpi_double_precision,  &
    mpi_sum,mpi_comm_world,mpi_ierror)
CALL mpi_allreduce(partdw,partdwp,1,mpi_double_precision,  &
    mpi_sum,mpi_comm_world,mpi_ierror)
CALL mpi_barrier (mpi_comm_world, mpi_ierror)
partup = partupp
partdw = partdwp
#endif



npartup = NINT(partup)
npartdw = NINT(partdw)
!npartup = partup+0.000001D0 ! M : integer ???
!npartdw = partdw+0.000001D0


npartto =  npartup+ npartdw


RETURN
END SUBROUTINE act_part_num
#endif

!     ******************************

#ifdef REALSWITCH
SUBROUTINE calc_sicspr(rhosp,usicsp,q0state,nb)
#else
SUBROUTINE calc_sicsp(rhosp,usicsp,q0state,nb)
#endif
!     ******************************

!     computes SIC potential for state 'nb'.
!     input is
!       q0state   = wavefunction for s.p. state
!       nb        = number of state
!     output is
!       usicsp   = the s.p. SIC potential
!     output via common
!       enrear,enerpw   = rearrangement energies for 'nb'
!       encoulsp        = Coulomb energy for 'nb'

USE params
USE util, ONLY:prifld2
USE kinetic
#if(netlib_fft|fftw_cpu)
USE coulsolv
#endif
#if(findiff|numerov)
USE coulsolv,ONLY:solv_poisson
#endif


IMPLICIT NONE

#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

#ifdef REALSWITCH
REAL(DP), INTENT(IN) :: q0state(kdfull2)
#else
COMPLEX(DP), INTENT(IN) :: q0state(kdfull2)
#endif
REAL(DP), INTENT(IN OUT) :: rhosp(2*kdfull2)
REAL(DP), INTENT(IN OUT) :: usicsp(2*kdfull2)


REAL(DP),DIMENSION(:),ALLOCATABLE :: chpdftsp,couldif,rho1

LOGICAL :: testprint=.false.
INTEGER :: ind, idx, ishift, nb
REAL(DP) :: enrearsave, enerpwsave, enrear1, enrear2, enpw1, enpw2, encadd, encsum
!--------------------------------------------------------------------

IF(numspin.NE.2) STOP ' CALC_SICSP requires fullspin code'

ALLOCATE(chpdftsp(2*kdfull2))
ALLOCATE(couldif(kdfull2))
ALLOCATE(rho1(kdfull2))

enrearsave=enrear
enerpwsave=enerpw
enrear1=0D0
enrear2=0D0
enpw1  = 0D0
enpw2  = 0D0
encadd = 0D0

IF(occup(nb) > small) THEN
  ishift = (ispin(nrel2abs(nb))-1)*nxyz ! store spin=2 in upper block
  
!         density for s.p. state
  
  DO ind=1,nxyz
#ifdef REALSWITCH
    rhosp(ind)  = q0state(ind)*q0state(ind)
#else
    rhosp(ind)  = REAL(CONJG(q0state(ind))*q0state(ind),DP)
#endif
    rhosp(ind+nxyz) =  3-2*ispin(nrel2abs(nb)) ! M : 1 if spinup -1 if spin down
    rho1(ind)        = rhosp(ind)
  END DO
!  WRITE(*,*) 'SICSP: nb,norm-dens=',nb,dvol*SUM(rhosp(1:nxyz))
!       DFT for s.p. state
  
  CALL calc_lda(rhosp,chpdftsp)    !  --> enrear,enerpw

  IF (ispin(nrel2abs(nb)) == 1) THEN
    DO ind=1,nxyz
      usicsp(ind) = chpdftsp(ind)
    END DO
  ELSE
    DO ind=1,nxyz
      idx = ind+nxyz
      usicsp(idx) = chpdftsp(idx)
    END DO
  END IF
  
!         s.p. Coulomb and subtract from LDA
  
  DO ind=1,nxyz
    rho1(ind) = rhosp(ind)
  END DO
#if(gridfft)
  CALL falr(rho1,couldif,kdfull2)
#endif
#if(findiff|numerov)
  CALL solv_poisson(rho1(1),couldif,kdfull2)
#endif
  IF(testprint) THEN
    WRITE(11,'(a)') '     ','     '
    WRITE(11,'(a,i3)') '# NB=',nb
    CALL prifld2(11,rhosp,' density')
    CALL prifld2(11,couldif,' Coulomb')
    CALL prifld2(11,chpdftsp,' P&W')
  END IF

 
  encsum = 0D0
  IF (ispin(nrel2abs(nb)) == 1) THEN
    DO ind=1,nxyz
      usicsp(ind) = usicsp(ind)+couldif(ind)
      IF(directenergy) THEN
        encsum=encsum+rhosp(ind)*couldif(ind)
      END IF
    END DO
    IF(testprint) CALL prifld2(11,usicsp,' SIC pot')
  ELSE
    DO ind=1,nxyz
      idx = ind+nxyz
      usicsp(idx) = usicsp(idx)+couldif(ind)
      IF(directenergy) THEN
        encsum=encsum+rhosp(ind)*couldif(ind)
      END IF
    END DO
    IF(testprint) CALL prifld2(11,usicsp(nxyz+1),' SIC pot')
  END IF
  encoulsp = encsum*dvol/2D0
ELSE
  usicsp(1:nxyz) = 0D0
  encoulsp = 0D0
END IF
DEALLOCATE(chpdftsp)
DEALLOCATE(couldif)
DEALLOCATE(rho1)


RETURN
#ifdef REALSWITCH
END SUBROUTINE calc_sicspr
#else
END SUBROUTINE calc_sicsp
#endif



!     ******************************

#ifdef REALSWITCH
SUBROUTINE exchgr(q0,qex)
#else
SUBROUTINE exchg(q0,qex,nbe)
#endif

!     ******************************

!     exact exchange:
!       input for static version is set of wavefunctions on 'q0'
!       input for the dynamic version are 
!             the wavefunctions 'q0' on which the exchange acts
!             and the wf's 'psisavex' for the density matrix
!             which are communicated via module PARAMS
!       output are accumulated exchange wavefunctions on 'qex'
!       only the exchange for the one state 'nbe' is handled
!            in case of complex wavefunctions

USE params
#if(netlib_fft|fftw_cpu)
USE coulsolv
#endif
#if(findiff|numerov)
USE coulsolv,ONLY:solv_poisson
#endif

IMPLICIT NONE

#ifdef REALSWITCH
REAL(DP), INTENT(IN)                         :: q0(kdfull2,kstate)
REAL(DP), INTENT(OUT)                        :: qex(kdfull2,kstate)
INTEGER :: i,nbe
REAL(DP) :: sump
#else
COMPLEX(DP), INTENT(IN)                         :: q0(kdfull2)
COMPLEX(DP), INTENT(OUT)                        :: qex(kdfull2)
INTEGER,INTENT(IN)                              :: nbe
REAL(DP),DIMENSION(:),ALLOCATABLE :: acli
REAL(DP),DIMENSION(:),ALLOCATABLE :: rhi
COMPLEX(DP) :: rhoc
#endif

!       workspaces

REAL(DP),DIMENSION(:),ALLOCATABLE :: rh
REAL(DP),DIMENSION(:),ALLOCATABLE :: acl

LOGICAL,PARAMETER :: ttest=.false.
INTEGER :: ind,nb2


IF(ifsicp /= 5) STOP ' in EXCHANGE: wrong option IFSICP'

ALLOCATE(acl(2*kdfull2))
ALLOCATE(rh(2*kdfull2))

#ifdef COMPLEXSWITCH
ALLOCATE(rhi(2*kdfull2))
ALLOCATE(acli(2*kdfull2))
#endif

qex = 0D0

!     double loop over states:
!     outer loop for states on which exchange acts
!     inner loop for states with which exchange acts

dvol=dx*dy*dz
#ifdef REALSWITCH
DO nbe=1,nstate
#endif
  DO nb2=1,nstate
    IF(ispin(nrel2abs(nbe)) == ispin(nrel2abs(nb2)) &
       .AND. occup(nrel2abs(nb2)) > 0.5D0) THEN
      IF(ttest) WRITE(*,*) ' NBE,NB2,SPINS,OCCUP:',  &
          nbe,nb2,ispin(nrel2abs(nbe)),ispin(nrel2abs(nb2)), &
       occup(nrel2abs(nb2))
      
!       compute transition density
      
      DO ind=1,nxyz
#ifdef REALSWITCH
        rh(ind)= q0(ind,nb2)*q0(ind,nbe)
#else
        rhoc= CONJG(psisavex(ind,nb2))*q0(ind)
        rh(ind)= REAL(rhoc,DP)
        rhi(ind)= AIMAG(rhoc)
#endif
      END DO
      
!       the Coulomb potential for the transition density
!         (warning : counet inserts the esquar factor)
      
#if(gridfft)
      CALL falr(rh,acl,kdfull2)
#ifdef COMPLEXSWITCH
      CALL falr(rhi,acli,kdfull2)
#endif
#endif
#if(findiff|numerov)
      CALL solv_poisson(rh,acl,kdfull2)
#ifdef COMPLEXSWITCH
      CALL solv_poisson(rhi,acli,kdfull2)
#endif
#endif
      
!       accumulate on wavefunction
      
#ifdef REALSWITCH
      IF(ttest) THEN
        sump = 0D0
        DO  i=1,nxyz
          sump  = q0(i,nbe)*q0(i,nb2)*acl(i) + sump
        END DO
        WRITE(6,'(a,2i5,1pg12.4)')  &
            ' EXCHANGE: nbe,nb2,overlap=',nbe,nb2,sump*dvol
      END IF
#endif

#ifdef REALSWITCH
      qex(:,nbe)=  qex(:,nbe) - q0(:,nb2)*acl(:)
#else
      qex(:)=  qex(:) - psisavex(:,nb2)*CMPLX(acl(:),acli(:),DP)
#endif

#ifdef REALSWITCH
      IF(ttest) THEN
        sump = 0D0
        DO  i=1,nxyz
          sump  = q0(i,nbe)*qex(i,nbe) + sump
        END DO
        WRITE(6,'(a,2i5,1pg12.4)')  &
            ' EXCHANGE: nbe,nb2,overlap=',nbe,nb2,sump*dvol
      END IF
#endif
      
    END IF
  END DO
#ifdef REALSWITCH
END DO
#endif

DEALLOCATE(acl)
DEALLOCATE(rh)

#ifdef COMPLEXSWITCH
DEALLOCATE(acli)
DEALLOCATE(rhi)
#endif

RETURN
#ifdef REALSWITCH
END SUBROUTINE exchgr
#else
END SUBROUTINE exchg
#endif



