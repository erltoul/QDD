#include "define.h"



!     ******************************

#ifdef REALSWITCH
SUBROUTINE calc_adsicr(rho,aloc,q0)
#else
SUBROUTINE calc_adsic(rho,aloc,q0)
#endif

!     ******************************

!     computes local part of hamiltonian
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
REAL(DP) :: aloc(*),rho(*)
REAL(DP),DIMENSION(:),ALLOCATABLE :: rhosp,chpdftsp,coulsum,couldif
REAL(DP),DIMENSION(:),ALLOCATABLE :: rho1,rho2
!DIMENSION rhosp(2*kdfull2),chpdftsp(2*kdfull2)
!DIMENSION coulsum(kdfull2),couldif(kdfull2)
!DIMENSION rho1(kdfull2),rho2(kdfull2)
#if(fullspin)
COMMON /sicwork/ rhospu(2*kdfull2),rhospd(2*kdfull2),  &
    chpdftspu(2*kdfull2),chpdftspd(2*kdfull2)
#endif
!EQUIVALENCE (coulsum,w1),(couldif,w2)
!EQUIVALENCE (rhosp,w1)
!EQUIVALENCE (chpdftsp,w3)
!EQUIVALENCE (rho2,rhosp)

!-------------------------------------------------------------------------

IF(ifsicp /= 2) THEN
  STOP ' CALC_SIC called with wrong option IFSICP'
END IF

CALL act_part_num(npartup,npartdw,npartto)




!     check workspace

!IF(usew1) STOP ' in CALC_SIC: workspace W1 already active '
!IF(usew2) STOP ' in CALC_SIC: workspace W2 already active '
!IF(usew3) STOP ' in CALC_SIC: workspace W3 already active '
!IF(usew4) STOP ' in CALC_SIC: workspace W4 already active '
!usew1 = .true.
!usew2 = .true.
!usew3 = .true.
!usew4 = .true.

enrearsave=enrear
enerpwsave=enerpw
enrear1=0D0
enrear2=0D0
enpw1  = 0D0
enpw2  = 0D0
encadd = 0D0
#if(fullspin)

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


#else


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
#endif

!usew1 = .false.
!usew2 = .false.
!usew3 = .false.
!usew4 = .false.

!     correct Coulomb potential by 1/N

#if(fullspin)
!IF(usew1) STOP ' in CALC_SIC: workspace W1 already active '
!usew1 = .true.
!IF(usew2) STOP ' in CALC_SIC: workspace W2 already active '
!usew2 = .true.
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

!usew1 = .false.
!usew2 = .false.
DEALLOCATE(coulsum)
DEALLOCATE(couldif)
DEALLOCATE(rho1)
DEALLOCATE(rho2)


#else
! IF(usew1) STOP ' in CALC_SIC: workspace W1 already in use'
! usew1 = .true.

!     recalculate coulomb part for jellium

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
! usew1 = .false.
#endif


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

!     computes local part of hamiltonian
!     'time' <= 0  signal static iteration i.e. without laser field


USE params
USE kinetic
!USE coulsolv
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
! DIMENSION usicsp(2*kdfull2),rhosp(2*kdfull2)
! EQUIVALENCE (rhosp,w1)
#if(fullspin)
COMMON /sicwork/ rhospu(2*kdfull2),rhospd(2*kdfull2),  &
    chpdftspu(2*kdfull2),chpdftspd(2*kdfull2)
#endif
LOGICAL :: testprint
DATA testprint/.false./

!-------------------------------------------------------------------

#if(fullspin)
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
#if(directenergy)
      enpw2=enpw2+enerpw*occup(nb)
#endif
    END IF
#if(directenergy)
    encadd=encadd+encoulsp*occup(nb)
#endif
    
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
#if(directenergy)
enerpw   = enerpwsave-enpw1-enpw2-encadd
#endif
! usew1 = .false.
DEALLOCATE(rhosp)
DEALLOCATE(usicsp)

#else
STOP ' SIC-Slater requires full spin code'
#endif

RETURN
#ifdef REALSWITCH
END SUBROUTINE calc_slaterr
#else
END SUBROUTINE calc_slater
#endif




#if(kli)

!     ******************************

#ifdef REALSWITCH
SUBROUTINE calc_sicklir(rho,aloc,q0)
#else
SUBROUTINE calc_sickli(rho,aloc,q0)
#endif

!     ******************************

!     computes local part of hamiltonian
!     'time' <= 0  signal static iteration i.e. without laser field


USE params
USE kinetic
!USE coulsolv
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
! DIMENSION rhosp(2*kdfull2),usicsp(2*kdfull2)
! DIMENSION rho1(kdfull2),rho2(kdfull2)
#if(fullspin)
COMMON /sicwork/ rhospu(2*kdfull2),rhospd(2*kdfull2),  &
    chpdftspu(2*kdfull2),chpdftspd(2*kdfull2)
COMMON /kliwork/ rhokli(kdfull2,kstate), uslater(2*kdfull2),ukli(2*kdfull2)
#endif
! EQUIVALENCE (rhosp,w1)
LOGICAL :: testprint
DATA testprint/.false./
DATA itmaxkli/200/

!--------------------------------------------------------------------

#if(fullspin)
IF(ifsicp /= 4) STOP ' CALC_SICKLI called with wrong option IFSICP'

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

!     loop oover s.p. states
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
#if(directenergy)
      enpw2=enpw2+enerpw*occup(nb)
#endif
    END IF
#if(directenergy)
    encadd=encadd+encoulsp*occup(nb)
#endif
    IF (ispin(nrel2abs(nb)) == 1) THEN
      sum = 0D0
      DO ind=1,nxyz
        relup = rhosp(ind)/MAX(rhospu(ind),small)
        rhokli(ind,nb) = relup
        uslater(ind) = uslater(ind) + relup*usicsp(ind)
        sum = rhosp(ind)*usicsp(ind) + sum
      END DO
      IF(nb == nfup) avuslatup=sum*dvol
      sum = sum*dvol*sgnkli
      DO ind=1,nxyz
!old              uslater(ind) = uslater(ind)-sum*rhokli(ind,nb)
        uslater(ind) = uslater(ind)+sum*rhokli(ind,nb)
      END DO
      sumslup = sumslup+sum
    ELSE IF (ispin(nrel2abs(nb)) == 2) THEN
      sum = 0D0
      DO ind=1,nxyz
        idx = ind+nxyz
        reldw = rhosp(ind)/MAX(rhospd(ind),small)
        rhokli(ind,nb) = reldw
        uslater(idx) = uslater(idx) +reldw*usicsp(idx)
        sum = rhosp(ind)*usicsp(idx) + sum
      END DO
      IF(nb == nfdw) avuslatdw=sum*dvol
      sum = sum*dvol*sgnkli
      DO ind=1,nxyz
        idx = ind+nxyz
!old            uslater(idx) = uslater(idx)-sum*rhokli(ind,nb)
        uslater(idx) = uslater(idx)+sum*rhokli(ind,nb)
      END DO
      sumsldw = sumsldw+sum
    END IF
  END IF
END DO
encadd=encadd/2.0
enrear   = enrearsave-enrear1-enrear2
#if(directenergy)
enerpw   = enerpwsave-enpw1-enpw2-encadd
#endif

DEALLOCATE(rhosp)
DEALLOCATE(usicsp)

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
    w2(ind) = 0D0
  END DO
  DO nb=1,nstate
    IF(occup(nb) > small) THEN
      ishift = (ispin(nrel2abs(nb))-1)*nxyz
      sum = 0D0
      IF(ishift <= 0) THEN
        DO ind=1,nxyz
          sum = ukli(ind)*rhokli(ind,nb)*rhospu(ind) + sum
        END DO
        sum = sum*dvol*sgnkli
        DO ind=1,nxyz
          rho1(ind) = rhokli(ind,nb)*sum + rho1(ind)
        END DO
      ELSE
        DO ind=1,nxyz
          idx = ind+ishift
          sum = ukli(idx)*rhokli(ind,nb)*rhospd(ind) + sum
        END DO
        sum = sum*dvol*sgnkli
        DO ind=1,nxyz
          w2(ind) = rhokli(ind,nb)*sum + w2(ind)
        END DO
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
    sum = 0D0
    DO ind=1,nxyz
      sum = rhospu(ind)*rhokli(ind,nfup)*ukli(ind)+sum
    END DO
    sum = sum*dvol
    correct = sum-avuslatup
    DO ind=1,nxyz
      ukli(ind) = ukli(ind)-correct
    END DO
    sum = 0D0
    DO ind=1,nxyz
      idx = ind+nxyz
      sum = rhospd(ind)*rhokli(ind,nfdw)*ukli(idx)+sum
    END DO
    sum = sum*dvol
    correct = sum-avuslatdw
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
usew1 = .false.
#else
STOP ' IFSICP=4: KLI requires fullspin code'
#endif

RETURN
#ifdef REALSWITCH
END SUBROUTINE calc_sicklir
#else
END SUBROUTINE calc_sickli
#endif
#endif

!  presently only static version of 'calc_fullsicr'
#ifdef REALSWITCH
#if(fullsic)
#if(fullspin)


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
!USE coulsolv
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
#if(directenergy)
      enpw2=enpw2+enerpw*occup(nb)
#endif
    END IF
#if(directenergy)
    encadd=encadd+encoulsp*occup(nb)
#endif
    
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
#if(directenergy)
enerpw   = enerpwsave-enpw1-enpw2-encadd
#endif


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
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)
REAL(DP), INTENT(IN) :: q0(kdfull2,kstate)
!REAL(DP), INTENT(IN) :: akv(kdfull2)
REAL(DP), INTENT(IN OUT) :: qnew(kdfull2,kstate)
REAL(DP), INTENT(IN) :: aloc(2*kdfull2)

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


!     check workspace

! IF(usew1) STOP ' in SSTEP: workspace W1 already active '
! IF(usew2) STOP ' in SSTEP: workspace W2 already active '
! IF(usew3) STOP ' in SSTEP: workspace W3 already active '
! usew1 = .true.
! usew2 = .true.
! usew3 = .true.
ALLOCATE(q1(kdfull2))
ALLOCATE(q2(kdfull2))
#if(gridfft)
ALLOCATE(psipr(kdfull2))
#endif

dvol=dx*dy*dz
DO nbe=1,nstate
  ishift = (ispin(nbe)-1)*nxyz        ! store spin=2 in upper block
  
!       action of mean-field hamiltonian on 'hwfr'
  
  
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
      vol   = REAL(psipr(i))*REAL(psipr(i)) +imag(psipr(i))*imag(psipr(i))
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

!       compute symmtry condition

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


#else

SUBROUTINE calc_fullsicr(q0,qsic)

STOP ' fullsic requires fullspin code'
RETURN
END SUBROUTINE calc_fullsicr
#endif
#endif
#endif

#ifdef REALSWITCH

!     ******************************

SUBROUTINE act_part_num(npartup,npartdw,npartto)

!     ******************************

!     computes actual particle number for spin-up and spin-down.

USE params
USE kinetic
!USE coulsolv
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

#if(fullspin)

!     check workspace

!IF(usew2) STOP ' in CALC_SIC: workspace W2 already active '
!IF(usew3) STOP ' in CALC_SIC: workspace W3 already active '
!IF(usew4) STOP ' in CALC_SIC: workspace W4 already active '
!usew2 = .true.
!usew3 = .true.
!usew4 = .true.
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
#if(directenergy)
      encsum=encsum+rhosp(ind)*couldif(ind)
#endif
    END DO
    encoulsp = encsum*dvol
    IF(testprint) CALL prifld2(11,usicsp,' SIC pot')
  ELSE
    DO ind=1,nxyz
      idx = ind+nxyz
      usicsp(idx) = usicsp(idx)+couldif(ind)
#if(directenergy)
      encsum=encsum+rhosp(ind)*couldif(ind)
#endif
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
!usew2 = .false.
!usew3 = .false.
!usew4 = .false.
DEALLOCATE(chpdftsp)
DEALLOCATE(couldif)
DEALLOCATE(rho1)
#else
STOP ' CALC_SICSP requires fullspin code'
#endif

RETURN
#ifdef REALSWITCH
END SUBROUTINE calc_sicspr
#else
END SUBROUTINE calc_sicsp
#endif


