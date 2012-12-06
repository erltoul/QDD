#include "define.h"

#if(fullsic)
#ifdef REALSWITCH
!-----sicstep_gen------------------------------------------

SUBROUTINE sicstep_gen(psir,qnew,aloc,iter1)

!     Switchboard for symmetry condition in static full SIC.
!     (? yet to be tested and developed ?)

USE params
USE kinetic
#if(symmcond)
USE symcond
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: psir(kdfull2,kstate)
REAL(DP), INTENT(IN OUT)                     :: qnew(kdfull2,kstate)
!REAL(DP), INTENT(IN OUT)                     :: akv(kdfull2)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)
INTEGER, INTENT(IN)                      :: iter1


!REAL(DP) ::  rho(2*kdfull2)

!------------------------------------------------------------


#if(symmcond)

!     tentative propagation of symmetry condition
!     the stepping parameters are set here at compile time.

itbegin = itersicp6+10000
itmodul = 100
symmstep = 3.0 ! 0.8

IF(iter1 > itbegin .AND. MOD(iter1/itmodul,2) == 0) THEN
!        call symmcond_step(psir,qnew,symmstep,symcond_add)
!        call symmcond_direct(psir,qnew,symmstep,symcond_add)
  CALL symmcond_emin(psir,qnew,symmstep,symcond_add)
  WRITE(7,'(a,i5,1pg13.5)') ' iter,symm.cond.=',iter1,SQRT(symcond_add/nstate)
  WRITE(11,'(a,i5,1pg13.5)')  &
      ' iter,symm.cond.=',iter1,SQRT(symcond_add/nstate)
  IF(MOD(iter1,itmodul) == 1) WRITE(11,'(5x)')
  symcon = SQRT(symcond_add/nstate)
ELSE
  CALL sstep_sic(psir,aloc,iter1,qnew)
END IF
#else

CALL sstep_sic(psir,aloc,iter1,qnew)

#endif

RETURN
END SUBROUTINE sicstep_gen

#endif
#endif





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
IMPLICIT REAL(DP) (A-H,O-Z)!
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



IF(ifsicp == 1) THEN
  CALL calc_sicgamr(rho,aloc,q0)
ELSE IF(ifsicp == 2) THEN
  CALL calc_adsicr(rho,aloc,q0)
ELSE IF(ifsicp == 3) THEN
  CALL calc_slaterr(rho,aloc,q0)
ELSE IF(ifsicp == 4) THEN
  CALL calc_sicklir(rho,aloc,q0)
ELSE
  RETURN
END IF
#else
IF(ifsicp == 1) THEN
  CALL calc_sicgam(rho,aloc,q0)
ELSE IF(ifsicp == 2) THEN
  CALL calc_adsic(rho,aloc,q0)
ELSE IF(ifsicp == 3) THEN
  CALL calc_slater(rho,aloc,q0)
ELSE IF(ifsicp == 4) THEN
  CALL calc_sickli(rho,aloc,q0)
ELSE
  RETURN
END IF
#endif
RETURN
#ifdef REALSWITCH
END SUBROUTINE calc_sicr
#else
END SUBROUTINE calc_sic
#endif





!     ******************************

#ifdef REALSWITCH
SUBROUTINE calc_sicgamr(rho,aloc,q0)
#else
SUBROUTINE calc_sicgam(rho,aloc,q0)
#endif

!     ******************************

!     computes local part of SIC-GAM Hamiltonian
!     'time' <= 0  signal static iteration i.e. without laser field


USE params
USE kinetic
USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)
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
REAL(DP), INTENT(IN OUT) :: aloc(*),rho(*)

REAL(DP),DIMENSION(:),ALLOCATABLE :: rhosp,chpdftsp,coulsum,couldif
REAL(DP),DIMENSION(:),ALLOCATABLE :: rho1,rho2
LOGICAL,PARAMETER :: testprint=.false.
INTEGER,PARAMETER :: itmaxkli=200

CALL act_part_num(npartup,npartdw,npartto)

ALLOCATE(chpdftsp(2*kdfull2))
ALLOCATE(rhosp(2*kdfull2))
ALLOCATE(rho1(kdfull2))



enrearsave=enrear
enerpwsave=enerpw
enrear1=0D0
enrear2=0D0
enpw1  = 0D0
enpw2  = 0D0
encadd = 0D0

!   first block:  Slater or non-Coulomb part of GAM

IF(ifsicp == 1 .OR. ifsicp >= 3) THEN
  DO nb=1,nstate
    IF(occup(nb) > small) THEN
      ishift = (ispin(nrel2abs(nb))-1)*nxyz ! store spin=2 in upper block
      
!         density for s.p. state
      
      DO ind=1,nxyz
#ifdef REALSWITCH
        rhosp(ind)       = q0(ind,nb)*q0(ind,nb)
#else
        rhosp(ind)       = REAL(CONJG(q0(ind,nb))*q0(ind,nb))
#endif
        rhosp(ind+nxyz) =  3-2*ispin(nrel2abs(nb)) ! M : 1 if spinup -1 if spin down
        rho1(ind)        = rhosp(ind)
      END DO
      
!         DFT for s.p. state
      
      CALL calc_lda(rhosp,chpdftsp)
      
      IF (ispin(nrel2abs(nb)) == 1) THEN
        enrear1=enrear1+enrear*occup(nb)
        enpw1=enpw1+enerpw*occup(nb)
      ELSE
        enrear2=enrear2+enrear*occup(nb)
        IF(directenergy) THEN
          enpw2=enpw2+enerpw*occup(nb)
        END IF
      END IF
      IF(ifsicp == 1) THEN
        
!           GAM: subtract with global weight from pure LDA potential
        
        IF(ishift == 0) THEN
          fac = occup(nb)/(npartup*1.0)
        ELSE
          fac = occup(nb)/(1.0*npartdw)
        END IF
        DO ind=1,nxyz
          idx = ind+ishift
          aloc(idx) = aloc(idx) - fac*chpdftsp(idx)
        END DO
      ELSE
      END IF
    END IF
  END DO
  encadd=encadd/2.0
  enrear   = enrearsave-enrear1-enrear2
  IF(directenergy) THEN
    enerpw   = enerpwsave-enpw1-enpw2-encadd
  END IF
!        write(6,'(a,1pg13.5)') ' ensic=',enpw1+enpw2
ELSE
  STOP ' CALC_SICGAMlled with wrong option IFSICP'
END IF

DEALLOCATE(chpdftsp)
DEALLOCATE(rhosp)

!     correct Coulomb potential by 1/N

IF(numspin==2) THEN
  ALLOCATE(coulsum(kdfull2))
  ALLOCATE(couldif(kdfull2))
  ALLOCATE(rho2(kdfull2))

!       compute Coulomb

  DO ind=1,nxyz
    rho2(ind)=  rho(ind)
    rho1(ind)=  rho(ind)*rho(ind+nxyz)
  END DO
#if(gridfft)
  CALL falr(rho1(1),couldif,nx2,ny2,nz2,kdfull2)
!                 new computation of total Coul not needed - check
!      call falr(rho2(1),coulsum,nx2,ny2,nz2,kdfull2)
#endif
#if(findiff|numerov)
  CALL solv_fft(rho1(1),couldif,dx,dy,dz)
!      call solv_FFT(rho2(1),coulsum,dx,dy,dz)
#endif
  facup = 1.0/npartup
  facdw = 1.0/npartdw


  DO ind=1,nxyz
!        coulup    = 0.5*(coulsum(ind)+couldif(ind))
!        couldw    = 0.5*(coulsum(ind)-couldif(ind))
    coulup    = 0.5*(chpcoul(ind)+couldif(ind))    ! reuse total Coul
    couldw    = 0.5*(chpcoul(ind)-couldif(ind))
    aloc(ind) = aloc(ind) - coulup*facup
    idx       = ind + nxyz
    aloc(idx) = aloc(idx) - facdw*couldw
  END DO
  DEALLOCATE(coulsum)
  DEALLOCATE(couldif)
  DEALLOCATE(rho2)

ELSE

!     recalculate Coulomb part for jellium

#if(gridfft)
  IF(nion2 == 0) CALL falr(rho(1),chpcoul,nx2,ny2,nz2,kdfull2)
#endif
#if(findiff|numerov)
  STOP ' SIC-GAM not yet ready for nospin and finite diff.'
#endif

  fac = 1.0/npartto
  DO ind=1,nxyz
    aloc(ind) = aloc(ind) - fac*chpcoul(ind)
  END DO

END IF

RETURN
#ifdef REALSWITCH
END SUBROUTINE calc_sicgamr
#else
END SUBROUTINE calc_sicgam
#endif





!     ******************************

#ifdef REALSWITCH
SUBROUTINE calc_adsicr(rho,aloc,q0)
#else
SUBROUTINE calc_adsic(rho,aloc,q0)
#endif

!     ******************************

!     computes local part of Hamiltonian
!     'time' <= 0  signal static iteration i.e. without laser field


USE params
USE kinetic
USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP),PARAMETER :: sgnkli=-1D0

#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

#ifdef REALSWITCH
REAL(DP) :: q0(kdfull2,kstate)
#else
COMPLEX :: q0(kdfull2,kstate)
#endif
REAL(DP) :: aloc(2*kdfull2),rho(2*kdfull2)
REAL(DP),DIMENSION(:),ALLOCATABLE :: rhosp,chpdftsp,coulsum,couldif
REAL(DP),DIMENSION(:),ALLOCATABLE :: rho1,rho2
REAL(DP),DIMENSION(:),ALLOCATABLE :: rhospu,rhospd,chpdftspu,chpdftspd
!EQUIVALENCE (coulsum,w1),(couldif,w2)
!EQUIVALENCE (rhosp,w1)
!EQUIVALENCE (chpdftsp,w3)
!EQUIVALENCE (rho2,rhosp)

!-------------------------------------------------------------------------

IF(ifsicp /= 2) THEN
  STOP ' CALC_SIC called with wrong option IFSICP'
END IF

CALL act_part_num(npartup,npartdw,npartto)




IF(numspin==2) ALLOCATE(rhospu(2*kdfull2),rhospd(2*kdfull2),  &
                        chpdftspu(2*kdfull2),chpdftspd(2*kdfull2))

enrearsave=enrear
enerpwsave=enerpw
enrear1=0D0
enrear2=0D0
enpw1  = 0D0
enpw2  = 0D0
encadd = 0D0
IF(numspin==2) THEN

!     averaged spinup and spindown density

  IF(npartup > 0) THEN
    facuph = 0.5/npartup
  
    DO ind=1,nxyz
      rhospu(ind)=rho(ind)*(1.0+rho(ind+nxyz))*facuph
      rhospu(ind+nxyz)=1.0
    END DO
  END IF
  IF(npartdw > 0) THEN
    facdwh = 0.5/npartdw
  
    DO ind=1,nxyz
      rhospd(ind)=rho(ind)*(1.0-rho(ind+nxyz))*facdwh
      rhospd(ind+nxyz)=-1.0
    END DO
  END IF

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

  IF(npartup > 0) THEN
    DO ind=1,nxyz
      aloc(ind) = aloc(ind) - chpdftspu(ind)
    END DO
  END IF



  IF(npartdw > 0) THEN
    DO idx=nxyz+1,2*nxyz
      aloc(idx) = aloc(idx) - chpdftspd(idx)
    END DO
  END IF


  DEALLOCATE(rhospu,rhospd,chpdftspu,chpdftspd)

ELSE

  ALLOCATE(chpdftsp(2*kdfull2))
  ALLOCATE(rhosp(2*kdfull2))

  factotal=1.0/npartto

  DO ind=1,nxyz
    rhosp(ind)=rho(ind)*factotal
  END DO

!     DFT for averaged s.p. state

  CALL calc_lda(rhosp,chpdftsp)

!     subtract from pure LDA potential

  DO ind=1,nxyz
    aloc(ind) = aloc(ind) - chpdftsp(ind)
  END DO
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

  DO ind=1,nxyz
    rho2(ind)=  rho(ind)
    rho1(ind)=  rho(ind)*rho(ind+nxyz)
  END DO
#if(gridfft)
  CALL falr(rho1(1),couldif,nx2,ny2,nz2,kdfull2)
  IF(idielec == 0 .AND. nion2 > 0) THEN
    coulsum = chpcoul
  ELSE
    CALL falr(rho2(1),coulsum,nx2,ny2,nz2,kdfull2)
  END IF
#endif
#if(findiff|numerov)
  CALL solv_fft(rho1(1),couldif,dx,dy,dz)
  CALL solv_fft(rho2(1),coulsum,dx,dy,dz)
#endif
  facup = 1.0/npartup

  DO ind=1,nxyz
    coulup    = 0.5*(coulsum(ind)+couldif(ind))
    aloc(ind) = aloc(ind) - coulup*facup
  END DO


  IF (npartdw > 0) THEN
  
    facdw = 1.0/npartdw
    DO ind=1,nxyz
      couldw    = 0.5*(coulsum(ind)-couldif(ind))
      idx       = ind + nxyz
      aloc(idx) = aloc(idx) - facdw*couldw
    END DO
  
  END IF

  DEALLOCATE(coulsum)
  DEALLOCATE(couldif)
  DEALLOCATE(rho1)
  DEALLOCATE(rho2)

ELSE

!     recalculate Coulomb part for jellium

#if(gridfft)
  IF(nion2 == 0) THEN
    CALL falr(rho(1),chpcoul,nx2,ny2,nz2,kdfull2)
  END IF
#endif
#if(findiff|numerov)
  STOP ' ADSIC not yet ready for nospin and finite diff.'
#endif

  fac = 1.0/npartto
  DO ind=1,nxyz
    aloc(ind) = aloc(ind) - fac*chpcoul(ind)
  END DO

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
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

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
DATA testprint/.false./

!-------------------------------------------------------------------

IF(numspin.NE.2) STOP ' SIC-Slater requires full spin code'
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
    rhospu(ind)=MAX(rho(ind)*(1.0+rho(ind+nxyz))*0.5,small)
  END DO
END IF
IF(npartdw > 0) THEN
  DO ind=1,nxyz
    rhospd(ind)=MAX(rho(ind)*(1.0-rho(ind+nxyz))*0.5,small)
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
encadd=encadd/2.0
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
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
PARAMETER (sgnkli=-1.0)
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
! EQUIVALENCE (rhosp,w1)
LOGICAL :: testprint
DATA testprint/.false./
DATA itmaxkli/200/

!--------------------------------------------------------------------

IF(numspin.NE.2) STOP ' SIC-KLI requires full spin code'
IF(ifsicp /= 4) STOP ' CALC_SICKLI called with wrong option IFSICP'

ALLOCATE(rhokli(kdfull2,kstate),uslater(2*kdfull2),ukli(2*kdfull2))

CALL act_part_num(npartup,npartdw,npartto)

!     determine Fermi levels

IF(temp /= 0D0) STOP 'KLI not yet compatible with temperature'
nfup = 0
nfdw = 0
IF(sgnkli == -1.0) THEN
  efermup = -1000D0
  efermdw = -1000D0
  DO nb=1,nstate
    IF (ispin(nrel2abs(nb)) == 1) THEN      ! M : 1=up 2=down
      IF(amoy(nb) > efermup .AND. occup(nb) > 0.5) THEN
        nfup = nb
        efermup = amoy(nb)
      END IF
    ELSE
      IF(amoy(nb) > efermdw .AND. occup(nb) > 0.5) THEN
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
    rhospu(ind)=MAX(rho(ind)*(1.0+rho(ind+nxyz))*0.5,small)
  END DO
END IF
IF(npartdw > 0) THEN
  DO ind=1,nxyz
    rhospd(ind)=MAX(rho(ind)*(1.0-rho(ind+nxyz))*0.5,small)
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
encadd=encadd/2.0
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
  addn = 0.3
  addo = 1.0-addn
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
  IF(sgnkli == -1.0) THEN        ! asymptotic correction
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
  IF(sumup+sumdw < epsoro**2*2.0) GO TO 99
  sumslup = sumup
  sumsldw = sumdw
END DO
WRITE(6,'(a,2(1pg12.4))')  &
    ' KLI not converged: errors=',sumup-sumslup,sumdw-sumsldw
99   CONTINUE
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


!  presently only static version of 'calc_fullsicr'
#ifdef REALSWITCH
#if(fullsic)

!     ******************************

SUBROUTINE calc_fullsicr(q0,qsic)

!     ******************************

!     full SIC:
!       input is set of wavefunctions on 'q0'
!       output are SIC s.p. wavefunctions on 'qsic'

USE params
USE kinetic
#if(symmcond)
USE symcond
!COMMON /sicsav/usicall(kdfull2,kstate)
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN) :: q0(kdfull2,kstate)
REAL(DP), INTENT(IN OUT) :: qsic(kdfull2,kstate)


!       workspaces

REAL(DP),DIMENSION(:),ALLOCATABLE :: usicsp,rhosp
! DIMENSION usicsp(2*kdfull2),rhosp(2*kdfull2)
! COMMON /sicwork/ rhospu(2*kdfull2),rhospd(2*kdfull2),  &
!     chpdftspu(2*kdfull2),chpdftspd(2*kdfull2)
! EQUIVALENCE (rhospu,rhosp),(rhospd,usicsp)

LOGICAL :: ttest
DATA ttest/.false./

IF(numspin.NE.2) STOP ' full SIC requires full spin code'

!mb
enrearsave=enrear
enerpwsave=enerpw
enrear1=0D0
enrear2=0D0
enpw1  = 0D0
enpw2  = 0D0
encadd = 0D0
!mb/

ALLOCATE(rhosp(2*kdfull2))
ALLOCATE(usicsp(2*kdfull2))


!------------------------------------------------------------------

!     compute action of SIC potential

DO nb=1,nstate
  IF(occup(nb) > small) THEN
    ishift = (ispin(nrel2abs(nb))-1)*nxyz ! store spin=2 in upper block
#ifdef REALSWITCH
    CALL calc_sicspr(rhosp,usicsp,q0(1,nb),nb)
#else
    CALL calc_sicsp(rhosp,usicsp,q0(1,nb),nb)
#endif
#if(symmcond)
    DO ind=1,nxyz
      usicall(ind,nb) = usicsp(ind+ishift)
    END DO
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
        qsic(ind,nb) = usicsp(ind)*q0(ind,nb)
      END DO
    ELSE
      DO ind=1,nxyz
        idx = ind+nxyz
        qsic(ind,nb) = usicsp(idx)*q0(ind,nb)
      END DO
    END IF
  END IF
END DO
encadd=encadd/2.0
enrear   = enrearsave-enrear1-enrear2
IF(directenergy) THEN
  enerpw   = enerpwsave-enpw1-enpw2-encadd
END IF


DEALLOCATE(rhosp)
DEALLOCATE(usicsp)

RETURN
END SUBROUTINE calc_fullsicr



!     ******************************

SUBROUTINE sstep_sic(q0,aloc,iter,qnew)

!     ******************************


!     one static step with full SIC - only serial version

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
REAL(DP), INTENT(IN OUT)   :: q0(kdfull2,kstate)
!REAL(DP), INTENT(IN)      :: akv(kdfull2)
REAL(DP), INTENT(IN OUT)   :: qnew(kdfull2,kstate)
REAL(DP), INTENT(IN)       :: aloc(2*kdfull2)

!       workspaces

REAL(DP),DIMENSION(:),ALLOCATABLE :: q1
!,q1(kdfull2)
#if(gridfft)
COMPLEX(DP),DIMENSION(:),ALLOCATABLE :: q2,psipr
! COMPLEX :: q2(kdfull2),psipr(kdfull2)
#else
REAL(DP),DIMENSION(:),ALLOCATABLE :: q2
! REAL :: q2(kdfull2)
#endif
! EQUIVALENCE (q1(1),w1(1))
! EQUIVALENCE (q2(1),w2(1))              ! occupies also w3

LOGICAL :: tocc,tcpu
DATA tocc,tcpu/.true.,.true./
LOGICAL :: tdiag,tprojec
PARAMETER (tdiag=.false.)   ! switch to diagonal subtraction
PARAMETER (tprojec=.false.) ! switch to full 1ph projection

!------------------------------------------------------------------------


!      write(*,*) ' SSTEP with IFSICP=',ifsicp
IF(tcpu) CALL cpu_time(time_init)
IF(ifsicp /= 6) STOP ' SSTEP_SIC only for full SIC'
!                       action of SIC potentials returned on 'qnew'

!MB:      call calc_fullsic(q0,qnew)


ALLOCATE(q1(kdfull2))
ALLOCATE(q2(kdfull2))
#if(gridfft)
ALLOCATE(psipr(kdfull2))
#endif

dvol=dx*dy*dz
DO nbe=1,nstate
  ishift = (ispin(nbe)-1)*nxyz        ! store spin=2 in upper block
  
!       action of mean-field Hamiltonian on 'hwfr'
  
  
  IF(ipsptyp == 1) THEN
    CALL nonlocalr(q0(1,nbe),q1)
    enonlo(nbe)= rwfovlp(q0(1,nbe),q1)
    DO  i=1,nxyz
      q1(i)=q1(i)+q0(i,nbe)*aloc(i+ishift)
    END DO
  ELSE
    DO  i=1,nxyz
      q1(i)=q0(i,nbe)*aloc(i+ishift)
    END DO
  END IF
!                             subtract SIC potential for state NBE
  DO i=1,nxyz
    q1(i)=q1(i)-qnew(i,nbe)
  END DO
  
!      optionally compute expectation value of potential energy
  
  IF(MOD(iter,istinf) == 0) epotsp(nbe) = rwfovlp(q0(1,nbe),q1)
  
  
#if(gridfft)
  
!      action of the kinetic energy in momentum space
  
  CALL rftf(q0(1,nbe),psipr)
  CALL rftf(q1,q2)
  DO  i=1,nxyz
    q2(i) = psipr(i)*akv(i)
  END DO
  IF(MOD(iter,istinf) == 0) THEN
    sum0 = 0D0
    sumk = 0D0
    DO  i=1,nxyz
      vol   = REAL(psipr(i))*REAL(psipr(i)) +AIMAG(psipr(i))*AIMAG(psipr(i))
      sum0  = vol + sum0
      sumk  = vol*akv(i) + sumk
    END DO
    ekinsp(nbe) = sumk/sum0
  END IF
  CALL rfftback(q2,qnew(1,nbe))
  DO i=1,nxyz
    qnew(i,nbe)=q1(i)+qnew(i,nbe)
  END DO
#endif
#if(findiff|numerov)
  STOP ' full SIC not yet working for finite diifferences'
#endif
  
!       compute 'hmatrix' for next step
  
  DO nbc=1,nstate
    IF(ispin(nrel2abs(nbe)) == ispin(nrel2abs(nbc))) THEN
      hmatrix(nbc,nbe) = rwfovlp(q0(1,nbc),qnew(1,nbe))
    ELSE
      hmatrix(nbc,nbe) = 0D0
    END IF
  END DO
  
!     first round loop completed
  
END DO

!       compute symmetry condition

symcon = 0D0
WRITE(6,'(a)') ' hmatrix:'
DO na=1,nstate
  DO nb=1,nstate
    symcon = (hmatrix(na,nb)-hmatrix(nb,na))**2 + symcon
  END DO
  WRITE(6,'(20(1pg13.5))') (hmatrix(nb,na),nb=1,nstate)
END DO
symcon = SQRT(symcon/nstate)

!       symmetrize Hamiltonian matrix

DO na=1,nstate
  DO nb=1,na-1
    hmatrix(na,nb) = 0.5*(hmatrix(na,nb)+hmatrix(nb,na))
    hmatrix(nb,na) = hmatrix(na,nb)
  END DO
END DO

!     second round loop

DO nbe=1,nstate
  DO ind=1,nxyz
    q1(ind) = qnew(ind,nbe)       ! copy to workspace
  END DO
  
!       augment with non-diagonal constraint correction
  
  IF(tprojec) THEN
    CALL orthogwfr(q1,ispin(nbe),q0)
  ELSE
    IF(tdiag) THEN
      DO ind=1,nxyz
        q1(ind) = q1(ind)-hmatrix(nbe,nbe)*q0(ind,nbe)
      END DO
    ELSE
      DO nbc=1,nstate
        IF(ispin(nrel2abs(nbc)) == ispin(nrel2abs(nbe))) THEN
          DO ind=1,nxyz
            q1(ind) = q1(ind)-hmatrix(nbe,nbc)*q0(ind,nbc)
          END DO
        END IF
      END DO
    END IF
  END IF
  
!       the action of the SIC Hamiltonian is completed
!       thus compute average and variance
  
  IF(MOD(iter,istinf) == 0) THEN
    evarsp(nbe) = SQRT(MAX(rwfnorm(q1),small))
    epotsp(nbe) = hmatrix(nbe,nbe)-ekinsp(nbe)
  END IF
  
  
  
!     perform the damped gradient step
  
  IF(idyniter /= 0 .AND. iter > 100) e0dmp = MAX(ABS(amoy(nbe)),0.5D0)

#if(gridfft)
  IF(e0dmp > small) THEN
    CALL rftf(q1,q2)
    DO i=1,nxyz
      q2(i)=epswf / (akv(i) + e0dmp )*q2(i)
    END DO
    CALL rfftback(q2,q1)
    DO i=1,nxyz
      qnew(i,nbe)=q0(i,nbe)-q1(i)
    END DO
  ELSE
    DO i=1,nxyz
      qnew(i,nbe)=q0(i,nbe)-epswf*q1(i)
    END DO
  END IF
  
  
#else
  STOP ' full SIC requires fast Fourier code '
#endif
  
END DO

! usew1 = .false.
! usew2 = .false.
! usew3 = .false.
DEALLOCATE(q1)
DEALLOCATE(q2)
#if(gridfft)
DEALLOCATE(psipr)
#endif

!     copy to 'q0'

DO nbe=1,nstate
  DO i=1,nxyz
    q0(i,nbe) = qnew(i,nbe)
  END DO
END DO


!     Schmidt ortho-normalisation

CALL schmidt(q0)

!     readjust occupations numbers to actual s.p. energies

IF(tocc .AND. nclust < nstate*nph) CALL reocc()

IF(tcpu) THEN
  CALL cpu_time(time_fin)
  time_cpu = time_fin-time_init
!        write(6,'(a,1pg13.5)') ' CPU time in SSTEP',time_cpu
!        write(7,'(a,1pg13.5)') ' CPU time in SSTEP',time_cpu
END IF
WRITE(6,'(a,i5,6(f10.4))') 'iter,up/down,CPU=',iter,se(4),se(5),time_cpu
WRITE(7,'(a,i5,6(f10.4))') 'iter,up/down,CPU=',iter,se(4),se(5),time_cpu


RETURN
END SUBROUTINE sstep_sic


#endif
#endif

#ifdef REALSWITCH

!     ******************************

SUBROUTINE act_part_num(npartup,npartdw,npartto)

!     ******************************

!     computes actual particle number for spin-up and spin-down.

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

REAL :: partup,partdw
INTEGER :: npartup,npartdw,npartto


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
    mpi_sum,mpi_comm_world,ic)
CALL mpi_allreduce(partdw,partdwp,1,mpi_double_precision,  &
    mpi_sum,mpi_comm_world,ic)
CALL mpi_barrier (mpi_comm_world, mpi_ierror)
partup = partupp
partdw = partdwp
#endif



!mb      npartup = nint(partup+0.000001) ! M : integer ???
!mb      npartdw = nint(partdw+0.000001)
npartup = partup+0.000001D0 ! M : integer ???
npartdw = partdw+0.000001D0


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
USE kinetic
USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

#ifdef REALSWITCH
REAL(DP), INTENT(IN) :: q0state(kdfull2)
#else
COMPLEX(DP), INTENT(IN) :: q0state(kdfull2)
#endif
REAL(DP), INTENT(IN OUT) :: rhosp(2*kdfull2),usicsp(2*kdfull2)


REAL(DP),DIMENSION(:),ALLOCATABLE :: chpdftsp,couldif,rho1
!DIMENSION chpdftsp(2*kdfull2)
!DIMENSION couldif(kdfull2)
!DIMENSION rho1(kdfull2)
!EQUIVALENCE (couldif,w2)
!EQUIVALENCE (chpdftsp,w3)
!EQUIVALENCE (rho1,w4)

LOGICAL :: testprint
DATA testprint/.false./

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
!      write(*,*) ' CALC_SICSPR: nb,occup(nb)',nb,occup(nb)

IF(occup(nb) > small) THEN
  ishift = (ispin(nrel2abs(nb))-1)*nxyz ! store spin=2 in upper block
  
!         density for s.p. state
  
  DO ind=1,nxyz
#ifdef REALSWITCH
    rhosp(ind)  = q0state(ind)*q0state(ind)
#else
    rhosp(ind)  = REAL(CONJG(q0state(ind))*q0state(ind))
#endif
!          if(ispin(nrel2abs(nb)).eq.1) then
    rhosp(ind+nxyz) =  3-2*ispin(nrel2abs(nb)) ! M : 1 if spinup -1 if spin down
!            rhosp(ind+nxyz) =  1.0
!          else
!            rhosp(ind+nxyz) =  -1.0
!          endif
    rho1(ind)        = rhosp(ind)
  END DO
  
!       DFT for s.p. state
  
  CALL calc_lda(rhosp,chpdftsp)    !  --> enrear,enerpw
!        do ind=1,nxyz+nxyz
!          usicsp(ind) = chpdftsp(ind)
!        enddo
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
  CALL falr(rho1(1),couldif,nx2,ny2,nz2,kdfull2)
#endif
#if(findiff|numerov)
  CALL solv_fft(rho1(1),couldif,dx,dy,dz)
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
    encoulsp = encsum*dvol
    IF(testprint) CALL prifld2(11,usicsp,' SIC pot')
  ELSE
    DO ind=1,nxyz
      idx = ind+nxyz
      usicsp(idx) = usicsp(idx)+couldif(ind)
      IF(directenergy) THEN
        encsum=encsum+rhosp(ind)*couldif(ind)
      END IF
    END DO
    encoulsp = encsum*dvol
    IF(testprint) CALL prifld2(11,usicsp(nxyz+1),' SIC pot')
  END IF
ELSE
  DO i=1,nxyz
    usicsp(ind) = 0D0
  END DO
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
!USE kinetic
USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

#ifdef REALSWITCH
REAL(DP), INTENT(IN)                         :: q0(kdfull2,kstate)
REAL(DP), INTENT(OUT)                        :: qex(kdfull2,kstate)
#else
COMPLEX(DP), INTENT(IN)                         :: q0(kdfull2)
COMPLEX(DP), INTENT(OUT)                        :: qex(kdfull2)
REAL(DP),DIMENSION(:),ALLOCATABLE :: acli
REAL(DP),DIMENSION(:),ALLOCATABLE :: rhi
COMPLEX(DP) :: rhoc
COMPLEX(DP) :: wfovlp
#endif

!       workspaces

REAL(DP),DIMENSION(:),ALLOCATABLE :: rh
REAL(DP),DIMENSION(:),ALLOCATABLE :: acl

LOGICAL,PARAMETER :: ttest=.false.


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
!WRITE(6,'(a,10i5)') ' ispin:',ispin
!WRITE(6,'(a,10f8.2)') ' occup:',occup
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
        rh(ind)= REAL(rhoc)
        rhi(ind)= AIMAG(rhoc)
#endif
      END DO
      
!       the Coulomb potential for the transition density
!         (warning : counet inserts the esquar factor)
      
#if(gridfft)
      CALL falr(rh,acl,nx2,ny2,nz2,kdfull2)
#ifdef COMPLEXSWITCH
      CALL falr(rhi,acli,nx2,ny2,nz2,kdfull2)
#endif
#endif
#if(findiff|numerov)
      CALL solv_fft(rh,acl,dx,dy,dz)
#ifdef COMPLEXSWITCH
      CALL solv_fft(rhi,acli,dx,dy,dz)
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



