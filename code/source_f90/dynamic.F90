#include "define.h"

!------------------------------------------------------------

SUBROUTINE init_dynwf(psi)
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

!     initializes dynamical wavefunctions from static solution



COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)
!REAL(DP), INTENT(IN OUT)                     :: psir(kdfull2,kstate)


IF(ifsicp.EQ.5 .AND. ifexpevol .NE. 1) &
   STOP ' exact exchange requires exponential evolution'
IF(ifsicp.EQ.5 .OR. jstateoverlap == 1) ALLOCATE(psisavex(kdfull2,kstate))

!     use of saved real wavefunctions, read and copy to complex

CALL restart2(psi,outnam,.true.)
IF(ifsicp.EQ.5 .OR. jstateoverlap == 1) psisavex = psi

!     optionally rotate a 1ph state

IF(ABS(phangle)+ABS(phphase) > small) CALL phstate(psi)

!     boost the wavefunctions

IF (ekin0pp > 0D0) CALL init_velocel(psi)

IF ( eproj > 0D0) CALL init_projwf(psi)

!     optionally add scattering electron


IF (iscatterelectron /= 0) CALL init_scattel(psi)

! optionally reset occupation number by hand 
! (run first static, save on 'rsave.<name>', check occupation numbers,
!  start dynamics with 'istat=1', reset occupation numbers here).
! (present example for Na-8 with 12  states):
ifreset = 00
IF(ifreset==1) THEN
  occup(1) = 1D0
  occup(2) = 1D0
  occup(3) = 1D0
  occup(4) = 1D0
  occup(5) = 1D0
  occup(6) = 0D0
  occup(7) = 1D0
  occup(8) = 1D0
  occup(9) = 1D0
  occup(10) = 0D0
  occup(11) = 1D0
  occup(12) = 0D0
END IF



rvectmp(1)=1D0             ! ??? what for ?

!JM
#if(twostsic)
CALL init_vecs(vecs)
#endif
!JM

!     optional initial excitations

IF(iexcit == 0) THEN
  IF(ABS(shiftinix)+ABS(shiftiniy)+ABS(shiftiniz) > small) THEN
    CALL dipole_qp(psi,shiftinix,shiftiniy,shiftiniz,centfx,centfy,centfz)
  ELSE   
    CALL boost(psi)     !dipole boost of electr. cloud
  END IF
ELSE IF(iexcit == 1) THEN
  IF(nion2 /= 0) THEN
!     rotate ionic coordinates
    CALL mrote
    CALL calcpseudo()
  ELSE IF(nion2 == 0) THEN
!     rotate jelliumdensity:
    palph=falph
    pbeta=fbeta
    phexe=fhexe
    sqtest=0D0
    xion=1.0*nion
    CALL jelbak(xion,palph,pbeta,phexe,sqtest,1)
!    acc=0D0
!    DO i=1,nxyz
!      acc=acc+rhojel(i)
!    END DO
!    acc=acc*dvol
    acc = dvol*SUM(rhojel)
    WRITE(6,*) 'jellium-density after rotation:',acc
    WRITE(7,*) 'jellium-density after rotation:',acc
    WRITE(7,*) ' '
  END IF
END IF

IF(jstateoverlap == 1) THEN
  CALL stateoverl(psi,psisavex)
!  DEALLOCATE(psisavex)
END IF

RETURN
END SUBROUTINE init_dynwf
 


!-----init_velocel-------------------------------------------------

SUBROUTINE init_velocel(psi)

!     Boosts all electronic wavefunctions by the same velocity.
!     The absolute value is associated with the boost kinetic
!     energy 'ekin0pp' and the direction is given by
!     'vxn0', 'vyn0','vzn0'.
!     This routine is similar to 'boost', but it is used in
!     connection with ionic initialization to produce a cluster
!     where electrons and ions move with the same velocity.


USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(OUT)                     :: psi(kdfull2,kstate)

COMPLEX(DP) :: cfac

!------------------------------------------------------------------

v0 = SQRT(2.*ekin0pp/(amu(np(nion))*1836.0*ame))
rnorm = vxn0**2 + vyn0**2+ vzn0**2
rnorm = SQRT(rnorm)

IF (rnorm == 0) STOP 'Velocity vector not normalizable'

vxn0 = vxn0/rnorm*v0*ame
vyn0 = vyn0/rnorm*v0*ame
vzn0 = vzn0/rnorm*v0*ame

ind = 0
DO iz=minz,maxz
  z1=(iz-nzsh)*dz
  DO iy=miny,maxy
    y1=(iy-nysh)*dy
    DO ix=minx,maxx
      x1=(ix-nxsh)*dx
      ind = ind + 1
      arg = x1*vxn0+y1*vyn0+z1*vzn0
      cfac = CMPLX(COS(arg),SIN(arg),DP)
      DO nbe=1,nstate
        psi(ind,nbe)=cfac*psi(ind,nbe)
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE init_velocel

!-----init_projwf-------------------------------------------------

SUBROUTINE init_projwf(psi)

!     Boosts the electronic wavefunctions of the projectile 
!     by the same velocity, which norm is given by 'eproj' and
!     the direction by 'vpx', 'vpy',and 'vpz'.
!     This routine is similar to 'init_velocel'.


USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(OUT)                     :: psi(kdfull2,kstate)

COMPLEX(DP) :: cfac


!------------------------------------------------------------------

v0 = SQRT(2.*eproj/(amu(np(nion))*1836.0*ame))
rnorm = vpx**2 + vpy**2+ vpz**2
rnorm = SQRT(rnorm)

IF (rnorm == 0) STOP 'Velocity vector not normalizable'

vpx = vpx/rnorm*v0*ame
vpy = vpy/rnorm*v0*ame
vpz = vpz/rnorm*v0*ame

IF(taccel>0D0) RETURN               ! boost is done adiabatically

IF(init_lcao.NE.1) THEN
  WRITE(*,*) ' instantaneous acceleration only for INIT_LCAO==1'
  STOP  ' instantaneous acceleration only for INIT_LCAO==1'
END IF
ind = 0
DO iz=minz,maxz
  z1=(iz-nzsh)*dz
  DO iy=miny,maxy
    y1=(iy-nysh)*dy
    DO ix=minx,maxx
      x1=(ix-nxsh)*dx
      ind = ind + 1
      arg = x1*vpx+y1*vpy+z1*vpz
      cfac = CMPLX(COS(arg),SIN(arg),DP)
      DO nbe=1,nmaxst(nion)
        psi(ind,nbe)=cfac*psi(ind,nbe)
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE init_projwf

!-----init_scattel-------------------------------------------------

SUBROUTINE init_scattel(psi)

!     Adds one electron in scattering state as Gaussian wavepacket
!     with a certain velocity.
!     The bookkeeping fiedls are extended accoprdingly.


USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(OUT)                     :: psi(kdfull2,kstate)

COMPLEX(DP) :: cfac

!------------------------------------------------------------------

pnorm = scatterelectronvxn**2+scatterelectronvyn**2+ scatterelectronvzn**2
pnorm = SQRT(pnorm)

IF (pnorm == 0D0) STOP 'Momentum of Scatt. Electron vanishes!'

scatterelectronvxn = scatterelectronvxn/pnorm
scatterelectronvyn = scatterelectronvyn/pnorm
scatterelectronvzn = scatterelectronvzn/pnorm

pnorm = SQRT(scatterelectronenergy*2.*ame)

scatterelectronvxn = scatterelectronvxn*pnorm
scatterelectronvyn = scatterelectronvyn*pnorm
scatterelectronvzn = scatterelectronvzn*pnorm

nstate = nstate + 1
IF(nstate > kstate) STOP ' insufficient KSTATE in INIT_SCATTEL'
occup(nstate)=1D0
nrel2abs(nstate)=nstate
nabs2rel(nstate)=nstate
ispin(nstate)=1


fac = 1D0/SQRT(pi**1.5*scatterelectronw**3)

ind = 0

DO iz=minz,maxz
  z1=(iz-nzsh)*dz
  DO iy=miny,maxy
    y1=(iy-nysh)*dy
    DO ix=minx,maxx
      x1=(ix-nxsh)*dx
      ind = ind + 1
      arg = (x1-scatterelectronx)*scatterelectronvxn+  &
          (y1-scatterelectrony)*scatterelectronvyn+  &
          (z1-scatterelectronz)*scatterelectronvzn
      rr = (x1 - scatterelectronx)**2 &
          + (y1 - scatterelectrony)**2 &
          + (z1 - scatterelectronz)**2
      fr = fac*EXP(-rr/2./scatterelectronw**2)
      psi(ind,nstate) = CMPLX(fr,0D0,DP)*CMPLX(COS(arg),SIN(arg),DP)
    END DO
  END DO
END DO

RETURN
END SUBROUTINE init_scattel

!-----tstep---------------------------------------------------------

SUBROUTINE tstep(q0,aloc,rho,it)

!     one electronic time step by TV splitting method.

!     'itsub' indicates the number of subiteration before
!     next analyzing step (routine 'info').
!     'itsub=1' is the first call in a series and
!     'itsub=ipasinf' is the last call.
!g     Now, 'ipasinf' gives the step of computation of the electronic force
!g     on the cluster ion.

!     For pure electronic propagation one has the option to
!     reduce the number of local unitary steps. The last half-step
!     is omitted (exept in case of the last call 'itsub=ipasinf')
!     and the first local step is doubled instead (except for the
!     first call 'itsub=1'). This option is switched on by the
!     run time switch 'iffastpropag'.

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN OUT)                  :: q0(kdfull2,kstate)
REAL(DP), INTENT(IN)                         :: aloc(2*kdfull2)
!COMPLEX(DP), INTENT(IN OUT)                  :: ak(kdfull2)
REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
INTEGER, INTENT(IN)                      :: it
COMPLEX(DP),DIMENSION(:),ALLOCATABLE :: q1,q2
COMPLEX(DP) :: rii,cfac

! CHARACTER (LEN=1) :: inttostring1
! CHARACTER (LEN=2) :: inttostring2
! 
! CHARACTER (LEN=1) :: str1
! CHARACTER (LEN=2) :: str2
! CHARACTER (LEN=3) :: str3
! CHARACTER (LEN=4) :: str4
! CHARACTER (LEN=5) :: str5
! CHARACTER (LEN=6) :: str6
! CHARACTER (LEN=7) :: str7
! CHARACTER (LEN=8) :: str8
! CHARACTER (LEN=9) :: str9

LOGICAL :: tenerg

#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)

CALL  mpi_comm_rank(mpi_comm_world,myn,icode)
#else
myn = 0
#endif


ALLOCATE(q1(2*kdfull2))
ALLOCATE(q2(kdfull2))

CALL cpu_time(time_init)
IF (ntref > 0 .AND. it > ntref) nabsorb = 0           ! is that the correct place?


itsub = MOD(it,ipasinf) + 1

ri = -dt1*0.5D0
dt = dt1*0.5D0





!     half time step in coordinate space
!     local phase field on workspace 'q1'

nlocact = numspin*nxyz
DO ind=1,nlocact
  pr=-dt*aloc(ind)
  q1(ind)=CMPLX(COS(pr),SIN(pr),DP)
END DO
DO nb=1,nstate
  ishift = (ispin(nrel2abs(nb))-1)*nxyz
!  CALL cmult3d(q0(1,nb),q1(1+ishift))
  q0(:,nb) = q1(ishift+1:ishift+kdfull2)*q0(:,nb)
END DO



!     half non-local step

IF(ipsptyp == 1 .AND. tnonlocany) THEN
  DO nb=1,nstate
    tenerg = itsub == ipasinf
    CALL nonlocstep(q0(1,nb),q1,q2,dt,tenerg,nb,6)   ! 4
  END DO
END IF


!       one full time step for the kinetic energy

DO nb=1,nstate
#if(gridfft)
  IF(iffastpropag == 1) THEN
    CALL kinprop(q0(1,nb),q1)
  ELSE
    CALL fftf(q0(1,nb),q1)
!    CALL cmult3d(q1,ak)
    q1 = ak*q1
    CALL fftback(q1,q0(1,nb))
  END IF
#endif
#if(findiff|numerov)
  CALL d3mixpropag (q0(1,nb),dt1)
#endif
  
END DO


!old       tfs = tfs + (dt - dt1)*0.0484/(2.*ame)






!     half non-local step

IF(ipsptyp == 1 .AND. tnonlocany) THEN
  DO nb = 1,nstate
    tenerg = .false. !   itsub.eq.ipasinf
    CALL nonlocstep(q0(1,nb),q1,q2,dt,tenerg,nb,6)   ! 4
  END DO
END IF



DEALLOCATE(q1)
DEALLOCATE(q2)



!     new density and local potential
!     (this is already the density at the end of the step,
!      because it is unchanged by unitary potential step)
!     propagation of substrate dipoles is done in 'dyn_mfield'.

CALL dyn_mfield(rho,aloc,q0,dt)




!     re-open workspace

ALLOCATE(q1(2*kdfull2))

!     half time step in coordinate space:

nup = numspin*nxyz
DO ind=1,nup
  pr=-dt*aloc(ind)
  q1(ind)=CMPLX(COS(pr),SIN(pr),DP)
END DO
DO nb=1,nstate
  ishift = (ispin(nrel2abs(nb))-1)*nxyz
  q0(:,nb) = q1(ishift+1:ishift+kdfull2)*q0(:,nb)
END DO

!     finally release workspace

DEALLOCATE(q1)


!      call manualPES(q0)



CALL cpu_time(time_fin)
time_cpu = time_fin-time_init
IF(myn == 0)THEN
  WRITE(6,'(a,1pg13.5)') ' CPU time in TSTEP',time_cpu
  WRITE(7,'(a,1pg13.5)') ' CPU time in TSTEP',time_cpu
  CALL FLUSH(6)
END IF

IF (izforcecorr /= -1) THEN
  CALL checkzeroforce(rho,aloc)
  CALL checkzeroforce(rho,chpdft)
END IF

IF ((jescmask > 0 .AND. MOD(it,jescmask) == 0) .OR. &
    (jescmaskorb > 0 .AND. MOD(it,jescmaskorb) == 0)  ) CALL  escmask(it)

RETURN
END SUBROUTINE tstep

!-----dyn_mfield---------------------------------------------------

SUBROUTINE dyn_mfield(rho,aloc,psi,dt)

!     The Coulomb part of the mean field.

!     Input:
!      rho    = electron density
!      psi    = complex wavefunctions
!      dt     = actual time step
!               dt=0D0 forces adiabatic computation of substrate dipoles
!     Output:
!      aloc   = local mean-field potential

USE params
#if(twostsic)
USE twost
#endif
IMPLICIT REAL(DP) (A-H,O-Z)


REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)
COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)
REAL(DP), INTENT(IN)                         :: dt



COMPLEX(DP) :: wfovlp

!----------------------------------------------------------------

CALL calcrho(rho,psi)
CALL coul_mfield(rho)
!      if (isurf.ne.0 .and. nc+ne+nk.gt.0)
IF (isurf /= 0 .AND. NE > 0) THEN
  IF(dt == 0D0) THEN
    CALL valence_step(rho,dt,.false.)
  ELSE
    CALL valence_step(rho,dt,.true.)
  END IF
END IF
CALL calclocal(rho,aloc)          ! LDA part of the potential
IF(ifsicp > 0 .AND.ifsicp <= 6) THEN
  CALL calc_sic(rho,aloc,psi)
!JM :  Generalized Slater and FSIC
#if(twostsic)
ELSE IF(ifsicp >= 7)THEN
  IF(symutbegin < itmax) itut = symutbegin+1  ! force symmetry condition
  CALL calc_utwf(psi,psiut,itut)
!ccccccJM     Generalized Slater pot
  IF(ifsicp == 7)THEN
    ifsicp=3
    CALL calc_sic(rho,aloc,psiut)
    ifsicp=7
!ccccccJM     DSIC
  ELSE IF(ifsicp == 8) THEN
    CALL calc_fullsic(psiut,qnewut)
  END IF
#endif
!JM
ELSE IF(ifsicp == 6) THEN
  STOP ' that kind of SIC not valid for dynamics'
END IF

RETURN
END SUBROUTINE dyn_mfield


! !-----cmult3D---------------------------------------------------------
! 
! SUBROUTINE cmult3d(q0,q1)
! 
! !     to multiply complex 3D array 'q0' by other array 'q1'
! !     resulting in new 'q0' on output
! 
! USE params
! !USE kinetic
! IMPLICIT REAL(DP) (A-H,O-Z)
! 
! COMPLEX(DP), INTENT(OUT)                     :: q0(kdfull2)
! COMPLEX(DP), INTENT(IN)                      :: q1(kdfull2)
! 
! 
! 
! !--------------------------------------------------------------------
! 
! DO i=1,kdfull2
!   q0(i) = q1(i)*q0(i)
! END DO
! 
! RETURN
! END SUBROUTINE cmult3d


!-----boost--------------------------------------------------------

SUBROUTINE boost(q0)        ! boost with given 'centfx,y,z'

!     boosts electronic wavefunctions by 'centfx', 'centfy',
!     and 'centfz'.

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)


COMPLEX(DP), INTENT(OUT)                     :: q0(kdfull2,kstate)
!, INTENT(IN OUT)                         :: ! boost wi
!REAL(DP), INTENT(IN OUT)                     :: y


!--------------------------------------------------------------------

DO nbe=1,nstate
!        write(6,*) ' before BOOST: nb,norm=',nbe,wfnorm(q0(1,nbe))
  IF(ispidi == 1) THEN
    aspin = 3-2*ispin(nbe)
    actsz = 0.5D0*centfz*aspin
    actsy = 0.5D0*centfy*aspin
    actsx = 0.5D0*centfx*aspin
    WRITE(6,*) ' state nr.',nbe,' : boost with spin=',aspin
    WRITE(7,*) ' state nr.',nbe,' : boost with spin=',aspin
  ELSE
    actsz = centfz
    actsy = centfy
    actsx = centfx
!          write(6,*) ' actsxyz=',actsx,actsy,actsz
  END IF
  ind = 0
  DO iz=minz,maxz
    z1=(iz-nzsh)*dz*actsz
    DO iy=miny,maxy
      y1=(iy-nysh)*dy*actsy
      DO ix=minx,maxx
        x1=(ix-nxsh)*dx*actsx
        ind=ind+1
        q0(ind,nbe) = EXP(CMPLX(0D0,x1+y1+z1,DP))*q0(ind,nbe)
      END DO
    END DO
  END DO
!        write(6,*) ' after BOOST: nb,norm=',nbe,wfnorm(q0(1,nbe))
END DO

RETURN
END SUBROUTINE boost


!-----info ------------------------------------------------------------

SUBROUTINE info(psi,rho,aloc,it)

!     information on energy observables: single particle energies,
!     kinetic energy, total energy, ionic energy, ...

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

COMPLEX(DP), INTENT(IN OUT)              :: psi(kdfull2,kstate)
REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
REAL(DP), INTENT(IN OUT)                 :: aloc(2*kdfull2)
!REAL(DP), INTENT(IN OUT)                 :: akv(kdfull2)
INTEGER, INTENT(IN)                      :: it

REAL(DP), PARAMETER :: alpha_ar=10.6             !  for VdW
REAL(DP) ::  en(kstate)
COMPLEX(DP),ALLOCATABLE :: qtmp(:)
REAL(DP),ALLOCATABLE :: current(:,:)
REAL(DP) ::  estar(2),estarETF(2)
COMPLEX(DP) :: wfovlp


LOGICAL :: topenf
LOGICAL,PARAMETER :: ttest=.FALSE.

!------------------------------------------------------------------


!DO nb=1,nstate
!  CALL testgradient(psi(1,nb))
!END DO


#if(parayes)
CALL  mpi_comm_rank(mpi_comm_world,myn,icode)
#else
myn = 0
#endif
tinfs=it*dt1*0.0484/2.0/ame
tfs=it*dt1*0.0484
IF(idielec == 1) THEN
  CALL energ_dielec(rho)
ELSE
  ecrhoimage = 0D0
END IF

!test      call prifld(aloc,'ALOC1       ')


!     compute single-particle energies and related quantities

IF(ifsicp==5)  psisavex = psi

eshell=0D0
esh1=0D0
enonlc = 0D0

DO nb=1,nstate
  spnorm(nb) = wfovlp(psi(1,nb),psi(1,nb))
  CALL  calc_ekin(psi(1,nb),ekin)
  ishift = (ispin(nrel2abs(nb))-1)*nxyz
  CALL calc_epot(psi(1,nb),aloc(ishift+1), epot,enonlo(nb),nb,it)
  ekinsp(nb) = ekin
  ehilf= epot
  epot =epot+ekin
  amoy(nb)=epot
  epot = epot + ekin
  en(nb)=epot*occup(nb)
  eshell=eshell+en(nb)
  esh1=esh1+ekin*occup(nb)
  enonlc   = enonlc + enonlo(nb)*occup(nb)
!  WRITE(6,'(a,i2,4(a,f9.5))') 'level:',nrel2abs(nb), '  ekin='  &
!      ,ekin,'  epot=',ehilf,'  esp=',amoy(nb) ,'  enonlo=', enonlo(nb)
#if(!parayes)
  WRITE(6,'(2(a,i2),5(a,f9.5))') 'level:',nrel2abs(nb), &
       '  ispin=',ispin(nrel2abs(nb)),' occup=',occup(nb), &
       '  ekin='  &
      ,ekin,'  epot=',ehilf,'  esp=',amoy(nb) ,'  enonlo=', enonlo(nb)
#endif
END DO

CALL spmoms(psi,6)

WRITE(441,'(f10.3,10(1pg13.5))') tfs,enonlo(1:nstate)
CALL FLUSH(441)


tstinf = jstinf > 0 .AND. MOD(it,jstinf)==0
IF(tstinf) then
  IF(ifsicp==5)  psisavex = psi
#if(!parayes)
  ALLOCATE(qtmp(kdfull2))
  DO nb=1,nstate
    qtmp = psi(:,nb)
    CALL hpsi(qtmp,aloc,nb,1)  
    CALL cproject(qtmp,qtmp,ispin(nb),psi)
    spvariancep(nb) = SQRT(REAL(wfovlp(qtmp,qtmp)))
  END DO
  DEALLOCATE(qtmp)
#endif

#if(!parayes)
  CALL safeopen(77,it,jstinf,'pspenergies')
  WRITE(77,'(1f15.6,500f12.6)') tfs,(amoy(nb),nb=1,nstate)
  CALL flush(77)
  CALL safeopen(76,it,jstinf,'pspvariances')
  WRITE(76,'(1f15.6,500f12.6)') tfs,(spvariance(nb),nb=1,nstate)
  CALL flush(76)
  CALL safeopen(75,it,jstinf,'pspvariancesp')
  WRITE(75,'(1f15.6,500f12.6)') tfs,(spvariancep(nb),nb=1,nstate)
  CALL flush(75)
!  WRITE(6,'(a,f15.6)') ' s.p. energies and variances written at tfs=',tfs
#endif
ENDIF

IF(jstboostinv>0 .AND. MOD(it,jstboostinv)==0) THEN
  ALLOCATE(current(kdfull2,3))
  ALLOCATE(qtmp(kdfull2))
  CALL calc_current(current,psi)
  DO ind=1,kdfull2
   current(ind,:) = current(ind,:)/MAX(rho(ind),1D-5)
  END DO
  DO nb=1,nstate
    qtmp = psi(:,nb)
    CALL hpsi_boostinv(qtmp,aloc,current,rho,nbe)
    spenergybi(nb) = REAL(wfovlp(psi(1,nb),qtmp))
    spvariancebi(nb) = SQRT(REAL(wfovlp(qtmp,qtmp))-spenergybi(nb)**2)
  END DO
  CALL safeopen(91,it,jstboostinv,'pspenergybi')
  WRITE(91,'(1f15.6,500f12.6)') tfs,(spenergybi(nb),nb=1,nstate)
  CLOSE(91)
  CALL safeopen(92,it,jstboostinv,'pspvariancesbi')
  WRITE(92,'(1f15.6,500f12.6)') tfs,(spvariancebi(nb),nb=1,nstate)
  CLOSE(92)
  DEALLOCATE(current,qtmp)
END IF

#if(parayes)

CALL  prispe_parallele(6,it)

IF(ttest) WRITE(*,*) ' INFO: before allreduce. myn=',myn
CALL mpi_barrier (mpi_comm_world, mpi_ierror)
CALL mpi_allreduce(eshell,eshellp,1,mpi_double_precision,  &
    mpi_sum,mpi_comm_world,ic)
CALL mpi_allreduce(esh1,esh1p,1,mpi_double_precision,  &
    mpi_sum,mpi_comm_world,ic)
CALL mpi_allreduce(enonlc,enonlcp,1,mpi_double_precision,  &
    mpi_sum,mpi_comm_world,ic)
!DO iss=2,1,-1
!   CALL mpi_allreduce(estar(iss),estarp(iss),1,mpi_double_precision,  &
!        mpi_sum,mpi_comm_world,ic)
!END DO
CALL mpi_barrier (mpi_comm_world, mpi_ierror)
IF(ttest) WRITE(*,*) ' INFO: after allreduce. myn=',myn
esh1=esh1p
eshell=eshellp
enonlc=enonlcp
!DO iss=2,1,-1
!   estar(iss)=estarp(iss)
!ENDDO
#endif

DO iss=2,1,-1
   CALL calc_estar(psi,iss,estar(iss),estarETF(iss))
END DO

WRITE(7,*) 'ekintot=',esh1
WRITE(7,*) 'estar(1)=',estar(1),'   estar(2)=',estar(2)



!     rearrangement and background Coulomb energy


ecback=0D0
ecrho=0D0
IF(nion2 /= 0) THEN
  DO ind=1,nxyz
    ecback=ecback-rho(ind)*potion(ind)
    ecrho=ecrho+rho(ind)*(chpcoul(ind)-potion(ind))
  END DO
!  WRITE(*,*) ' ECBACK loop:',ecback*dvol/2.0
ELSE ! jellium case
  DO ind=1,nxyz
    ecback=ecback-rhojel(ind)*chpcoul(ind)
    ecrho=ecrho+rho(ind)*chpcoul(ind)
  END DO
END IF
ecback=ecback*dvol/2.0
ecrho=ecrho*dvol/2.0

#if(raregas)
IF(nc > 0 .AND. ivdw == 1)THEN
  esub = 0D0
  IF(nion2 /= 0) THEN
    DO ind=1,nxyz
      esub = esub + rho(ind)*potvdw(ind)
    END DO
  END IF
  esub = esub*dvol
  evdw = esub
  DO iss=1,nc
    DO ico=1,3
      evdw = evdw - 0.5D0*e2*alpha_ar*frho(iss,ico)*frho(iss,ico)/(REAL(nclust))
    END DO
  END DO
  eshell = eshell - esub
END IF
#endif

esh1=esh1
eshell=eshell/2.0  !(=t+v/2)


!     ionic contributions to the energy


ecorr = energ_ions()

ekion=0D0        ! kinetic energy of Na cores
ekinion=0D0      ! kinetic energy of GSM cores
ekinel=0D0       ! kinetic energy of GSM shells
ekinkat=0D0      ! kinetic energy of kations
IF(ionmdtyp > 0) THEN
  DO ion=1,nion
    ek=cpx(ion)*cpx(ion)+cpy(ion)*cpy(ion)+cpz(ion)*cpz(ion)
    xm=1836.0D0*amu(np(ion))*ame
    ek=ek/2D0/xm
    ekion=ekion+ek
  END DO
#if(raregas)
  DO ion=1,nc
    ek=pxc(ion)*pxc(ion)+pyc(ion)*pyc(ion)+pzc(ion)*pzc(ion)
    xm=1836.0*mion*ame
    ek=ek/2.0/xm
    ekinion=ekinion+ek
  END DO
  DO ion=1,NE
    ek=pxe(ion)*pxe(ion)+pye(ion)*pye(ion)+pze(ion)*pze(ion)
    xm=1836.0*me*ame
    ek=ek/2.0/xm
    ekinel=ekinel+ek
  END DO
  DO ion=1,nk
    ek=pxk(ion)*pxk(ion)+pyk(ion)*pyk(ion)+pzk(ion)*pzk(ion)
    xm=1836.0*mkat*ame
    ek=ek/2.0/xm
    ekinkat=ekinkat+ek
  END DO
#endif
END IF


!     final composition and print

energy=eshell+enrear+ecback+ecorr+enonlc/2.+ecrhoimage
IF(ivdw == 1)energy = energy + evdw
ecoul=ecback+ecrho+ecorr+ecrhoimage
etot = energy + ekion + ekinion + ekinel + ekinkat

IF (myn == 0 .AND. jenergy > 0 .AND. MOD(it,jenergy) == 0 ) THEN
  CALL safeopen(163,it,jenergy,'penergies')
  WRITE(163,'(1f14.6,22e24.15)') tfs, &
     &                eshell*2.-esh1,     &
     &                enrear,             &
     &                ekion,              &
     &                ekinion,            &
     &                ekinel,             &
     &                ekinkat,            &
     &                ecorr,              &
     &        enii,                       &
     &        enig,                       &
     &        engg,                       &
     &                2*ecback,           &
     &                2.*ecback+ecorr,    &
     &                ecrho-ecback-ecrhoimage,&
     &                enonlc,             &
     &                2.*ecback+ecorr+enonlc,&
     &                energy,            &
     &                etot, &
     &                elaser, &
     &                estar,estarETF
  CALL flush(163)
  CLOSE(163)
END IF

#if(parayes)
IF(myn == 0) THEN
#endif
  WRITE(6,'(a)') ' '
!mb        write(6,*) 'tot sp energ     = ',eshell*2.-esh1/2.
  WRITE(6,*) 'tot sp energ     = ',eshell*2.-esh1
  WRITE(6,*) 'rearge. energ    = ',enrear
  WRITE(6,*) 'e_coul:ion-ion   = ',ecorr
  WRITE(6,*) 'e_coul:el-ion    = ',2.*ecback
  WRITE(6,*) 'extern. energy   = ',2.*ecback+ecorr
  WRITE(6,*) 'hartree energy   = ',ecrho-ecback-ecrhoimage
  WRITE(6,*) 'nonlocal energy  = ',enonlc
  WRITE(6,*) 'sim.ann.energy   = ',2.*ecback+ecorr+enonlc
  WRITE(6,*) 'laser energy     = ',elaser
  WRITE(6,*) 'internal exc. energy per spin    = ',estar(1),estar(2)
  IF(idielec == 1) WRITE(6,*) 'image energy     = ',ecrhoimage
  IF(ivdw == 1) WRITE(6,*) 'vdw energy       = ',evdw
  WRITE(6,*) 'binding energy   = ',energy
  WRITE(6,*) 'total energy     = ',etot
  IF (isurf /= 0) THEN
    WRITE(6,*) 'adsorb.energy = ',etot-enerinfty
  END IF
  
  
  IF(jinfo > 0 .AND. MOD(it,jinfo) == 0) THEN
    INQUIRE(17,OPENED=topenf)
    IF(.NOT.topenf) OPEN(17,POSITION='append',FILE='infosp.'//outnam)
    WRITE(17,'(a,i5,f8.3,1pg13.5)') 'time_step,time,energy=',it,tinfs,etot
    CLOSE(17)
!    IF(.NOT.topenf)  CLOSE(17)
  END IF
  
#if(parayes)
END IF
#endif

IF(myn==0) THEN
  WRITE(6,'(a,f8.4,a,f7.2,a,3f10.5)')  &
    'In info: t=',tfs,' moments: monop.=',qe(1), ' dip.: ',qe(2),qe(3),qe(4)
  WRITE(6,'(a,f8.4,a,f7.2,a,3f10.5)')  &
    'In info: t=',tfs,' moments: quad et al=',qe(5), ' dip.: ',qe(6),qe(7),qe(8)
END IF

tstinf = .false.

RETURN
END SUBROUTINE info

!-----calc_ekin------------------------------------------------------------

SUBROUTINE calc_ekin(psin,ekinout)

!     calculates kinetic energy for single particle state with
!     complex wavefunction 'psin'.

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN OUT)                  :: psin(kdfull2)
!REAL(DP), INTENT(IN)                         :: akv(kdfull2)
REAL(DP), INTENT(OUT)                        :: ekinout


COMPLEX(DP),DIMENSION(:),ALLOCATABLE :: psi2
!COMPLEX(DP) :: psi2(kdfull2)
COMPLEX(DP) :: wfovlp
!EQUIVALENCE(psi2(1),w1(1))

!------------------------------------------------------------------

ALLOCATE(psi2(kdfull2))

#if(gridfft)
CALL fftf(psin,psi2)
sum0 = 0D0
sumk = 0D0
DO ii=1,kdfull2
  vol   = REAL(psi2(ii))*REAL(psi2(ii)) +imag(psi2(ii))*imag(psi2(ii))
  sum0  = vol + sum0
  sumk  = vol*akv(ii) + sumk
END DO
ekinout = sumk/sum0
#endif
#if(findiff|numerov)

!     exp.value of kinetic energy

CALL ckin3d(psi(1,nb),psi2)
sum0 = 0D0
acc = 0D0
DO i=1,nxyz
  acc = REAL(psi(i,nb))*REAL(psi2(i)) + imag(psi(i,nb))*imag(psi2(i))  &
      + acc
  sum0 = REAL(psi(i,nb))*REAL(psi(i,nb)) + imag(psi(i,nb))*imag(psi(i,nb))  &
      + sum0
END DO
ekinout = REAL(wfovlp(psi(1,nb),psi2))
#endif

DEALLOCATE(psi2)

RETURN
END SUBROUTINE calc_ekin


!-----calc_epot------------------------------------------------------------

SUBROUTINE calc_epot(psin,alocact,epotout,enonlocout,nb,it)

!     Calculates total potential energy 'epotout' and non-local
!     part of potential energy 'enonlocout' for the
!     single particle state with complex wavefunction 'psin'.
!     The local potential comes in through 'alocact'.

USE params
!USE kinetic
#if(twostsic)
USE twost
USE twostr, ONLY: ndims
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN)                      :: psin(kdfull2)
REAL(DP), INTENT(IN)                         :: alocact(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: epotout
REAL(DP), INTENT(OUT)                        :: enonlocout
INTEGER, INTENT(IN)                      :: nb
INTEGER, INTENT(IN OUT)                  :: it


COMPLEX(DP),DIMENSION(:),ALLOCATABLE :: psi2,qex
COMPLEX(DP) :: wfovlp
COMPLEX(DP) :: cf

!------------------------------------------------------------------

ALLOCATE(psi2(kdfull2))

epot=0D0


!     non local part of ps

IF(ipsptyp == 1) THEN
  CALL nonlocalc(psin,psi2,0)
  sumnon = 0D0
  DO i=1,nxyz
    sumnon = sumnon + REAL(psin(i),DP)*REAL(psi2(i),DP) &
                    + AIMAG(psin(i))*AIMAG(psi2(i))
  END DO
  enonlocout = sumnon*dvol
  DO i=1,nxyz
    psi2(i)=psi2(i)+alocact(i)*psin(i)
  END DO
ELSE
  DO i=1,nxyz
    psi2(i)=alocact(i)*psin(i)
  END DO
END IF

!JM : subtract SIC potential for state NB
#if(twostsic)
IF(ifsicp == 8) THEN
  is=ispin(nrel2abs(nb))
  DO na=1,ndims(is)
    nbe = nb - (is-1)*ndims(1)
    nae = na + (is-1)*ndims(1)
    cf = CONJG(vecs(nbe,na,is))
!            write(*,*)   ' nb,is,na,nbe,cf=',nb,is,na,nbe,cf
    DO i=1,nxyz
      psi2(i)=psi2(i)-qnewut(i,nae)*cf
    END DO
  END DO
END IF
#endif
!JM

IF(ifsicp==5) THEN
  ALLOCATE(qex(kdfull2))
  CALL exchg(psin,qex,nb)
  psi2 = psi2 + qex
  DEALLOCATE(qex)
END IF

epot = REAL(wfovlp(psin,psi2))   ! develop real overlap

!      write(*,*) ' EPOT: nb,epot=',nb,epot

DEALLOCATE(psi2)

RETURN
END SUBROUTINE calc_epot

!-----calc_estar------------------------------------------------

SUBROUTINE calc_estar(psin,iss,excit,excitETF)
!SUBROUTINE calc_estar(psin,is,excit,ekint)

!     Compute excitation energy as E* = Ekin - Ekin(Thomas-Fermi)
!     of the whole cluster.

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

COMPLEX(DP), INTENT(IN)                      :: psin(kdfull2,kstate)
INTEGER, INTENT(IN)                          :: iss
REAL(DP), INTENT(OUT)                        :: excit,excitETF

REAL(DP),DIMENSION(:),ALLOCATABLE            :: arho,gradrho
COMPLEX(DP),DIMENSION(:),ALLOCATABLE         :: gradrhok
#if(parayes)
REAL(DP),DIMENSION(:),ALLOCATABLE            :: arhop
#endif
REAL(DP) :: fact,factgrad

LOGICAL,PARAMETER :: extendedTF=.TRUE.
REAL(DP),PARAMETER :: rholimit=1D-10

!------------------------------------------------------------

ALLOCATE(arho(kdfull2))
#if(parayes)
ALLOCATE(arhop(kdfull2))
#endif

#if(parayes)
CALL mpi_barrier (mpi_comm_world, mpi_ierror)
#endif

arho=0D0
ekintot=0D0
DO nb=1,nstate
   j=nrel2abs(nb)
   IF (ispin(j).NE.iss) THEN
      CYCLE
   ELSE
      ekintot = ekintot + occup(nb)*ekinsp(nb)
      DO i=1,kdfull2
         arho(i)=arho(i) + &
            occup(nb)*(REAL(psin(i,nb),DP)**2+AIMAG(psin(i,nb))**2)
      END DO
   END IF
END DO

#if(parayes)
CALL mpi_barrier (mpi_comm_world, mpi_ierror)
CALL mpi_allreduce(arho,arhop,kdfull2,mpi_double_precision,  &
                   mpi_sum,mpi_comm_world,ic)
CALL mpi_allreduce(ekintot,ekintotp,1,mpi_double_precision,  &
                   mpi_sum,mpi_comm_world,ic)
CALL mpi_barrier (mpi_comm_world, mpi_ierror)
#endif

#if(parayes)
arho = arhop
DEALLOCATE(arhop)
#endif


excit=0D0
DO i=1,kdfull2
   excit=excit+arho(i)**1.666666666667D0
END DO
!fact=dvol*0.6D0*(6.0D0*pi**2)**0.666666666667D0
excit=excit*dvol*h2m*0.6D0*(6.0D0*pi**2)**0.666666666667D0

IF(extendedTF) THEN

  ALLOCATE(gradrho(kdfull2),gradrhok(kdfull2))
  factgrad=dvol*h2m/18.0D0
!  factgrad=dvol*h2m/36.0D0
  excitETF = excit

! x derivative
  gradrho = log(arho)
!  CALL rftf(arho,gradrhok)
  CALL rftf(gradrho,gradrhok)
  CALL gradient(gradrhok,gradrhok,1)
  CALL rfftback(gradrhok,gradrho)
  DO i=1,kdfull2
!    IF (arho(i).ne.0D0) 
    IF (arho(i).gt.rholimit) &
      excitETF=excitETF +factgrad*(gradrho(i))**2*arho(i)
!      excitETF=excitETF +factgrad*(gradrho(i))**2/arho(i)
  END DO

! y derivative
  gradrho = log(arho)
!  CALL rftf(arho,gradrhok)
  CALL rftf(gradrho,gradrhok)
  CALL gradient(gradrhok,gradrhok,2)
  CALL rfftback(gradrhok,gradrho)
  DO i=1,kdfull2
!    IF (arho(i).ne.0D0) 
    IF (arho(i).gt.rholimit) &
      excitETF=excitETF +factgrad*(gradrho(i))**2*arho(i)
!      excitETF=excitETF +factgrad*(gradrho(i))**2/arho(i)
  END DO

! z derivative
  gradrho = log(arho)
!  CALL rftf(arho,gradrhok)
  CALL rftf(gradrho,gradrhok)
  CALL gradient(gradrhok,gradrhok,3)
  CALL rfftback(gradrhok,gradrho)
  DO i=1,kdfull2
!    IF (arho(i).ne.0D0) 
    IF (arho(i).gt.rholimit) &
      excitETF=excitETF +factgrad*(gradrho(i))**2*arho(i)
!      excitETF=excitETF +factgrad*(gradrho(i))**2/arho(i)
  END DO

END IF

anorm = SUM(arho)*dvol

DEALLOCATE(arho)
DEALLOCATE(gradrho)
DEALLOCATE(gradrhok)

WRITE(7,'(a,7(1pg13.5))') ' norm,TF,ETF=',anorm,excit,excitETF


#if(parayes)
excitETF=ekintotp-excitETF
excit=ekintotp-excit
#else
excitETF=ekintot-excitETF
excit=ekintot-excit
#endif

write(7,'(a,2i5,3(1pg13.5))') &
    ' CALC_ESTAR: myn,iss,ekintot,excit,excitETF',  &
                  myn,iss,ekintot,excit,excitETF

RETURN
END SUBROUTINE calc_estar


!-----mrote ------------------------------------------------------------

SUBROUTINE mrote

!      Rotates ionic configuration by 'phirot' around axis 'irotat:
!        irotat=1 -> rotate around x-axis
!        irotat=2 -> rotate around y-axis
!        irotat=3 -> rotate around z-axis
!        irotat=4 -> rotate around diagonal axis
!      The angle 'phirot' is to be given in degree.
!      Conifguration and parameters are communciated via 'common'

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
REAL(DP) :: vecin(3),vecout(3),vecalpha(3)

!------------------------------------------------------------------

! determine vector of rotation angles (in radian)
  IF(irotat == 1) THEN
    vecalpha(1)  = phirot*PI/180D0
    vecalpha(2)  = 0D0
    vecalpha(3)  = 0D0
  ELSE IF(irotat == 2) THEN
    vecalpha(1)  = 0D0
    vecalpha(2)  = phirot*PI/180D0
    vecalpha(3)  = 0D0
  ELSE IF(irotat == 3) THEN
    vecalpha(1)  = 0D0
    vecalpha(2)  = 0D0
    vecalpha(3)  = phirot*PI/180D0
  ELSE IF(irotat == 4) THEN
    vecalpha(1)  = phirot*PI/180D0
    vecalpha(2)  = vecalpha(1)
    vecalpha(3)  = vecalpha(1)
  END IF

  OPEN(79,FILE='rotated-config.res')
  WRITE(79,'(a)') '  ion     x    y    z'

  DO ion=1,nion

    vecin(1) = cx(ion)
    vecin(2) = cy(ion)
    vecin(3) = cz(ion)

    CALL rotatevec3D(vecin,vecout,vecalpha)

    cx(ion) = vecout(1) 
    cy(ion) = vecout(2) 
    cz(ion) = vecout(3) 

    WRITE(7,'(/a,i4)') 'ion no:',ion
    WRITE(7,'(a,f12.4)') 'cxnew=',cx(ion)
    WRITE(7,'(a,f12.4)') 'cynew=',cy(ion)
    WRITE(7,'(a,f12.4)') 'cznew=',cz(ion)
    WRITE(79,'(3f12.4)') cx(ion),cy(ion),cz(ion)

  END DO

  RETURN
  END

!  !-----mrote ------------------------------------------------------------
!  
!  SUBROUTINE mrote_old
!  
!  !      rotates ionic configuration by 'phirot' around axis 'irotat:
!  !        irotat=1 -> rotate around x-axis
!  !        irotat=2 -> rotate around y-axis
!  !        irotat=3 -> rotate around z-axis
!  !        irotat=4 -> rotate around diagonal axis
!  !      conifguration and parameters are communciated via 'common'
!  
!  USE params
!  USE kinetic
!
!  IMPLICIT REAL(DP) (A-H,O-Z)
!  REAL(DP) :: help(3),trm(3,3)
!  REAL(DP) :: oldion(20,3)
!  
!  !------------------------------------------------------------------
!  
!  phir1=phirot * pi / 180D0
!  
!  IF(irotat == 3) THEN
!  !      rotation around z-axis:
!    WRITE(6,'(a,f7.2,a)')  &
!        '-> ionic coordinates rotated around z-axis by ',phirot,' degrees'
!    WRITE(7,'(a,f7.2,a)')  &
!        '-> ionic coordinates rotated around z-axis by ',phirot,' degrees'
!    trm(1,1)=COS(phir1)
!    trm(1,2)=SIN(phir1)
!    trm(1,3)=0D0
!    trm(2,1)=-SIN(phir1)
!    trm(2,2)=COS(phir1)
!    trm(2,3)=0D0
!    trm(3,1)=0D0
!    trm(3,2)=0D0
!    trm(3,3)=1D0
!    
!  ELSE IF(irotat == 2) THEN
!  !      rotation around y-axis:
!    WRITE(6,'(a,f7.2,a)')  &
!        '-> ionic coordinates rotated around y-axis by ',phirot,' degrees'
!    WRITE(7,'(a,f7.2,a)')  &
!        '-> ionic coordinates rotated around y-axis by ',phirot,' degrees'
!    trm(1,1)=COS(phir1)
!    trm(1,2)=0D0
!    trm(1,3)=SIN(phir1)
!    trm(2,1)=0D0
!    trm(2,2)=1D0
!    trm(2,3)=0D0
!    trm(3,1)=-SIN(phir1)
!    trm(3,2)=0D0
!    trm(3,3)=COS(phir1)
!    
!  ELSE IF(irotat == 1) THEN
!  !      rotation around x-axis:
!    WRITE(6,'(a,f7.2,a)')  &
!        '-> ionic coordinates rotated around x-axis by ',phirot,' degrees'
!    WRITE(7,'(a,f7.2,a)')  &
!        '-> ionic coordinates rotated around x-axis by ',phirot,' degrees'
!    trm(1,1)=1D0
!    trm(1,2)=0D0
!    trm(1,3)=0D0
!    trm(2,1)=0D0
!    trm(2,2)=COS(phir1)
!    trm(2,3)=SIN(phir1)
!    trm(3,1)=0D0
!    trm(3,2)=-SIN(phir1)
!    trm(3,3)=COS(phir1)
!    
!  ELSE IF(irotat == 4) THEN
!  !      rotation around diagonal axis:
!    WRITE(6,'(a,f7.2,a)')  &
!        '-> ionic coordinates rotated around diag. by ',phirot,' degrees'
!    WRITE(7,'(a,f7.2,a)')  &
!        '-> ionic coordinates rotated around diag. by ',phirot,' degrees'
!    trm(1,1)=0.33333333+0.66666667*(COS(phir1)-SIN(phir1))
!    trm(1,2)=0.40824829*(1-COS(phir1)+SIN(phir1))
!    trm(1,3)=0.23570226*(1-COS(phir1)+SIN(phir1))
!    trm(2,1)=0.40824829*(1-COS(phir1)+SIN(phir1))
!    trm(2,2)=0.5D0*(1+COS(phir1))
!    trm(2,3)=0.28867513*(1-COS(phir1)-2*SIN(phir1))
!    trm(3,1)=0.23570226*(1+SIN(phir1)-COS(phir1))
!    trm(3,2)=0.28867513*(1-COS(phir1)-2*SIN(phir1))
!    trm(3,3)=0.16666667*(1+5*COS(phir1)+4*SIN(phir1))
!  END IF
!  
!  DO ion=1,nion
!    oldion(ion,1) = cx(ion)
!    oldion(ion,2) = cy(ion)
!    oldion(ion,3) = cz(ion)
!  END DO
!  
!  OPEN(79,FILE='rotated-config.res')
!  WRITE(79,'(a)') '  ion     x    y    z'
!  DO ion=1,nion
!    DO i=1,3
!      help(i)=0
!      DO j=1,3
!        help(i)=help(i) + trm(i,j)*oldion(ion,j)
!      END DO
!      cx(ion)=help(1)
!      cy(ion)=help(2)
!      cz(ion)=help(3)
!    END DO
!    WRITE(7,'(/a,i4)') 'ion no:',ion
!    WRITE(7,'(a,f12.4)') 'cxnew=',cx(ion)
!    WRITE(7,'(a,f12.4)') 'cynew=',cy(ion)
!    WRITE(7,'(a,f12.4)') 'cznew=',cz(ion)
!    WRITE(79,'(3f12.4)') cx(ion),cy(ion),cz(ion)
!  END DO
!  CLOSE(79)
!  
!  RETURN
!  END SUBROUTINE mrote




!-----instit----------------------------------------------------

SUBROUTINE instit (psi)

!   calculation of angular momentum

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)

COMPLEX(DP), ALLOCATABLE :: akx(:),aky(:),akz(:),q2(:)
COMPLEX(DP), ALLOCATABLE :: jtx(:),jty(:),jtz(:)
COMPLEX(DP) :: jalpha

!------------------------------------------------------------------

ALLOCATE(akx(kdfull2),q2(kdfull2),aky(kdfull2),akz(kdfull2), &
         jtx(kdfull2),jty(kdfull2),jtz(kdfull2))

dkx=pi/(dx*REAL(nx))
dky=pi/(dy*REAL(ny))
dkz=pi/(dz*REAL(nz))
!      eye=cmplx(0.0,1.0)
!      nxyf=nx2*ny2
!      nyf=nx2

ind=0
DO i3=1,nz2
  IF(i3 >= (nz+1)) THEN
    zkz=(i3-nz2-1)*dkz
  ELSE
    zkz=(i3-1)*dkz
  END IF
  
  DO i2=1,ny2
    IF(i2 >= (ny+1)) THEN
      zky=(i2-ny2-1)*dky
    ELSE
      zky=(i2-1)*dky
    END IF
    
    DO i1=1,nx2
      IF(i1 >= (nx+1)) THEN
        zkx=(i1-nx2-1)*dkx
      ELSE
        zkx=(i1-1)*dkx
      END IF
      
      ind=ind+1
      akx(ind)=-zkx*eye
      aky(ind)=-zky*eye
      akz(ind)=-zkz*eye
    END DO
  END DO
END DO

!  we going to calc. jx,jy,jz


DO i=1,kdfull2
  jtx(i)=CMPLX(0D0,0D0,DP)
  jty(i)=CMPLX(0D0,0D0,DP)
  jtz(i)=CMPLX(0D0,0D0,DP)
END DO


DO nb=1,nstate
  o=occup(nb)
  
  CALL fftf(psi(1,nb),q2)
  DO ind=1,kdfull2
    q2(ind)=q2(ind)*akx(ind)
  END DO
  CALL fftback(q2,q2)
  DO ind=1,kdfull2
    test=eye/2.0*(CONJG(psi(ind,nb))*q2(ind) -psi(ind,nb)*CONJG(q2(ind)))
    
    jalpha=test
    jtx(ind)=jtx(ind)-o*jalpha
  END DO
  
  
  CALL fftf(psi(1,nb),q2)
  DO ind=1,kdfull2
    q2(ind)=q2(ind)*aky(ind)
  END DO
  CALL fftback(q2,q2)
  DO ind=1,kdfull2
    test=eye/2.0*(CONJG(psi(ind,nb))*q2(ind) -psi(ind,nb)*CONJG(q2(ind)))
    jalpha=test
    jty(ind)=jty(ind)-o*jalpha
  END DO
  
  
  CALL fftf(psi(1,nb),q2)
  DO ind=1,kdfull2
    q2(ind)=q2(ind)*akz(ind)
  END DO
  CALL fftback(q2,q2)
  DO ind=1,kdfull2
    test=eye/2.0*(CONJG(psi(ind,nb))*q2(ind) -psi(ind,nb)*CONJG(q2(ind)))
    jalpha=test
    jtz(ind)=jtz(ind)-o*jalpha
  END DO
  
END DO
! end loop over state


ajx=0D0
ajy=0D0
ajz=0D0


ind=0
DO iz=minz,maxz
  z1=(iz-nzsh)*dz
  DO iy=miny,maxy
    y1=(iy-nysh)*dy
    DO ix=minx,maxx
      x1=(ix-nxsh)*dx
      ind=ind+1
      ajpx=y1*jtz(ind)-z1*jty(ind)
      ajpy=z1*jtx(ind)-x1*jtz(ind)
      ajpz=x1*jty(ind)-y1*jtx(ind)
      ajx=ajx+ajpx
      ajy=ajy+ajpy
      ajz=ajz+ajpz
    END DO
  END DO
END DO
ajx=ajx*dvol
ajy=ajy*dvol
ajz=ajz*dvol
WRITE(6,'(a,3f12.4)') 'moments',ajx,ajy,ajz
WRITE(6,*)

DEALLOCATE(akx,q2,aky,akz,jtx,jty,jtz)


RETURN
END SUBROUTINE instit



!-----nonlocstep---------------------------------------------------

SUBROUTINE nonlocstep(qact,q1,q2,ri,tenerg,nb,norder)

!     Computes one time step for the non-local potentials
!     using exponential evolution.
!     qact     = array for actual wavefunction to be propagated
!     q1,q2    = auxiliary wavefunction  arrays
!     ri       = size of time step
!     tenerg   = (logical) switch to cumulate non-local energy
!     nb       = number of state which is propagated
!     norder   = order of step (up to 6, 4 or 6 recommended)
!
!     !! Note: This routine should be replaced by truly exponential
!              propagation within the non-local plaquettes.
!
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(OUT)                     :: qact(kdfull2)
COMPLEX(DP), INTENT(IN)                      :: q1(kdfull2)
COMPLEX(DP), INTENT(IN)                      :: q2(kdfull2)
REAL(DP), INTENT(IN)                         :: ri
LOGICAL, INTENT(IN)                      :: tenerg
INTEGER, INTENT(IN OUT)                  :: nb
INTEGER, INTENT(IN)                      :: norder



COMPLEX(DP) :: rii,cfac


!---------------------------------------------------------------------


CALL nonlocalc(qact,q1,0)
IF(tenerg) THEN !  add nonloc.pot energy
  sumadd = 0D0
  DO  i=1,nxyz
    sumadd  = REAL(qact(i))*REAL(q1(i)) +imag(qact(i))*imag(q1(i))  + sumadd
  END DO
  enonlo(nb) = sumadd*dvol
  epotsp(nb) = sumadd*dvol + epotsp(nb)
END IF

!test            do ind=1,nxyz
!test               qact(ind) = qact(ind)
!test     &              - ri*eye*q1(ind) - ri2*q2(ind)
!test     &              + ri*ri*ri/6.0*eye*q3(ind) + ri2*ri2/6.0*q4(ind)
!test     &           -ri2*ri2*ri/30.*eye*q5(ind)-ri2*ri2*ri2/90.*q6(ind)
!test            enddo
rii=ri*eye
DO ind=1,nxyz
  qact(ind) = qact(ind) - rii*q1(ind)
END DO

CALL nonlocalc(q1,q2,0)
ri2=ri*ri/2.0
DO ind=1,nxyz
  qact(ind) = qact(ind) - ri2*q2(ind)
END DO
IF(norder <= 2) RETURN

CALL nonlocalc(q2,q1,0)
cfac = ri*ri*ri/6.0*eye
DO ind=1,nxyz
  qact(ind) = qact(ind) + cfac*q1(ind)
END DO
IF(norder <= 3) RETURN

CALL nonlocalc(q1,q2,0)
rfac = ri2*ri2/6.0
DO ind=1,nxyz
  qact(ind) = qact(ind) + rfac*q2(ind)
END DO
IF(norder <= 4) RETURN

CALL nonlocalc(q2,q1,0)
cfac = -ri2*ri2*ri/30D0*eye
DO ind=1,nxyz
  qact(ind) = qact(ind) + cfac*q1(ind)
END DO
IF(norder <= 5) RETURN

CALL nonlocalc(q1,q2,0)
rfac = -ri2*ri2*ri2/90D0
DO ind=1,nxyz
  qact(ind) = qact(ind) + rfac*q2(ind)
END DO

RETURN
END SUBROUTINE nonlocstep



!-----init_dynprotocol-------------------------------------------------

SUBROUTINE init_dynprotocol(rho,aloc,psi)

!     initializes file for protocol of dynamics

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)
!REAL(DP), INTENT(IN OUT)                     :: akv(kdfull2)
COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)


LOGICAL :: topenf

!---------------------------------------------------------------------


IF(irest <= 0) THEN                    !  write file headers
  
#if(simpara)
  IF(jdip /= 0) THEN
#else
  IF(myn == 0 .AND. jdip /= 0) THEN
#endif
    OPEN(8,STATUS='unknown',FORM='formatted',FILE='pdip.'//outnam)
    WRITE(8,'(a)') ' & '
    WRITE(8,'(a)') 'x:time'
    WRITE(8,'(a)') 'y:dipole-momenta'
    WRITE(8,'(a)') 's: dist(4.0) scal(0.8) thickn(0.6) xmin(0.0)'
    WRITE(8,'(a)') 'n: time in fs  dipole-momenta: x   y   z'
    WRITE(8,'(a)') 'H:   X        Yl              Yd             Yq'
    CLOSE(8)
  END IF
  
#if(simpara)
  IF(jdiporb /= 0) THEN
#else
  IF(myn == 0 .AND. jdiporb /= 0) THEN
#endif
    OPEN(810,STATUS='unknown',FORM='formatted',FILE='pdiporb.x.'//outnam)
    WRITE(810,'(a)') 'protocol of s.p. moments: time,x-dipole of orbitals'
    CLOSE(810)
    OPEN(811,STATUS='unknown',FORM='formatted',FILE='pdiporb.y.'//outnam)
    WRITE(811,'(a)') 'protocol of s.p. moments: time,y-dipole of orbitals'
    CLOSE(811)
    OPEN(812,STATUS='unknown',FORM='formatted',FILE='pdiporb.z.'//outnam)
    WRITE(812,'(a)') 'protocol of s.p. moments: time,z-dipole of orbitals'
    CLOSE(812)
  END IF

#if(simpara)
  IF(jmp /= 0) THEN
#else
  IF(jmp /= 0 .AND. myn == 0) THEN
#endif
    OPEN(803,STATUS='unknown',FILE='pMP.'//outnam)
    WRITE(803,'(a)') '# nr. of meas.points, nr. of state'
    WRITE(803,'(a,2i6)') '# ',nmps,nstate_all
    DO i=1,nmps
      WRITE(803,'(a,i6,a,3f12.1,3i4)') '# Point ',i , ' at : ',  &
          getxval(imps(i)),getyval(imps(i)),getzval(imps(i)), &
          nint(getxval(imps(i))/dx),nint(getyval(imps(i))/dy), &
          nint(getzval(imps(i))/dz)
    END DO
!    OPEN(803,STATUS='unknown',FORM='unformatted',FILE='pMP.'//outnam)    ! cPW
!    WRITE(803) nmps,nstate_all                                           
!    DO i=1,nmps                                                          
!      WRITE(803) i,  &                                                   
!          getxval(imps(i)),getyval(imps(i)),getzval(imps(i)), &
!          nint(getxval(imps(i))/dx),nint(getyval(imps(i))/dy), &
!          nint(getzval(imps(i))/dz)
!    END DO
    CLOSE(803)
  END IF
  
  IF(jnorms /= 0) THEN
    OPEN(806,STATUS='unknown',FILE='pescOrb.'//outnam)
    WRITE(806,'(a,i6,a,3f12.1)') '# tfs, 1.0-norms of orbitals'
    OPEN(808,STATUS='unknown',FILE='pproba.'//outnam)
    WRITE(808,'(a,i6,a,3f12.1)') '# tfs, charge state probabilities'
!    CLOSE(806)
  END IF
  
  IF(e0 /= 0) THEN
    OPEN(38,STATUS='unknown',FILE='plaser.'//outnam)
    WRITE(38,*) ' & '
    WRITE(38,*) ' tfs,Ex,Ey,Ez,power,elaser '
!    CLOSE(38)
  END IF
  
  IF(jgeomel /= 0) THEN
    OPEN(608,STATUS='unknown',FORM='formatted', FILE='pgeomel.'//outnam)
    WRITE(608,*) '&'
!    CLOSE(608)
  END IF
  
  IF(jquad /= 0) THEN
    OPEN(9,STATUS='unknown',FORM='formatted',FILE='pquad.'//outnam)
    WRITE(9,'(a)') ' & '
    WRITE(9,'(a)') 'x:time'
    WRITE(9,'(a)') 'y:quadrupole-momenta'
    WRITE(9,'(a)') 's: dist(4.0) scal(0.8) thickn(0.6) xmin(0.0)'
    WRITE(9,'(a)') 'n: time in fs  dipole-momenta: x   y   z'
    WRITE(9,'(a)') 'H:   X   Y1  Y2  Y3  Y4  Y5  Y6'
    CLOSE(9)
  END IF
  
  IF(jinfo /= 0) THEN
    INQUIRE(17,OPENED=topenf)
    IF(.NOT.topenf) OPEN(17,POSITION='append',FILE='infosp.'//outnam)
    WRITE(17,'(a)') ' & '
    IF(.NOT.topenf)  CLOSE(17)
  END IF
  
  IF(jspdp /= 0) THEN
    OPEN(78,STATUS='unknown',FORM='formatted', FILE='pspdip.'//outnam)
    WRITE(78,'(a)') ' & '
    WRITE(78,'(a)') 'x:time'
    WRITE(78,'(a)') 'y:spindipole-momenta'
    WRITE(78,'(a)') 's: dist(4.0) scal(0.8) thickn(0.6) xmin(0.0)'
    WRITE(78,'(a)') 'n: time in fs  spindipole-momenta: x   y   z'
    WRITE(78,'(a)') 'H:   X        Yl              Yd            Yq'
    CLOSE(78)
  END IF
  
  IF(jang /= 0) THEN
    OPEN(68,STATUS='unknown',FORM='formatted', FILE='pangmo.'//outnam)
    WRITE(68,'(a)') ' & '
    WRITE(68,'(a)') 'x:time'
    WRITE(68,'(a)') 'y:x,y,z-component of angular momentum'
    WRITE(68,'(a)') 's: dist(4.0) scal(0.8) thickn(0.6) xmin(0.0)'
    WRITE(68,'(a)') 'n: time in fs       <l_x>'
    WRITE(68,'(a)') 'n:   x      l_x       l_y        l_z'
    WRITE(68,'(a)') 'H:   X      Yl        Yd         Yq'
    CLOSE(68)
  END IF
  
  IF(jenergy /= 0) THEN
    OPEN(163,STATUS='unknown',FORM='formatted', FILE='penergies.'//outnam)
    WRITE(163,*) 'col 1: time (fs), col 2:  total sp en.'
    WRITE(163,*) 'col 3: rearr. en, col 4:  kin. en. Na ions'
    WRITE(163,*) 'col 5: kin. cores,col 6:  kin. en. shells'
    WRITE(163,*) 'col 7: kin. kations,col 8:  pot. energy of ions'
    WRITE(163,*) 'col 9: ion-ion pot., col 10: ion-surf pot'
    WRITE(163,*) 'col 11: intra-surf. pot'
    WRITE(163,*) 'col 12: el-ion energy,col 13:  external en.'
    WRITE(163,*) 'col 14: hartree en.,col 15:  nonloc. en'
    WRITE(163,*) 'col 16: sim.ann. en.,col 17:  binding en.'
    WRITE(163,*) 'col 18: total en. [Ry]'
    WRITE(163,*) 'col 19: energy absorbed from laser [Ry]'
    WRITE(163,*) 'col 20/21: internal exc. energy (spin up/down)'
    CLOSE(163)
  END IF
  
  IF(jesc /= 0) THEN
    OPEN(23,STATUS='unknown',FORM='formatted', FILE='pescel.'//outnam)
    WRITE(23,'(a)') ' & '
    CLOSE(23)
    
    IF (jovlp /= 0) THEN
      OPEN(805,STATUS='unknown',FILE='povlp.'//outnam)
    END IF
    
!??         endif
    
    IF(jelf /= 0) THEN
      OPEN(33,STATUS='unknown',FORM='formatted', FILE='pelf.'//outnam)
      WRITE(33,'(a)') '# TD electron localization function along axes'
      CLOSE(33)
      OPEN(33,STATUS='unknown',FORM='formatted', FILE='pelf2Dxy.'//outnam)
      WRITE(33,'(a)') '# TD electron localization function in z=0 plane'
      CLOSE(33)
      OPEN(33,STATUS='unknown',FORM='formatted', FILE='pelf2Dyz.'//outnam)
      WRITE(33,'(a)') '# TD electron localization function in x=0 plane'
      CLOSE(33)
      OPEN(33,STATUS='unknown',FORM='formatted', FILE='pelf2Dxz.'//outnam)
      WRITE(33,'(a)') '# TD electron localization function in y=0 plane'
      CLOSE(33)
    END IF
    
    IF(jstinf /= 0) then
      OPEN(77,status='unknown',form='formatted',file='pspenergies.'//outnam)
      WRITE(77,'(a)') '# Single particle energies'
      CLOSE(77)
      OPEN(76,status='unknown',form='formatted',file='pspvariances.'//outnam)
      WRITE(76,'(a)') '# s.p. energy variances'
      CLOSE(76)
      OPEN(75,status='unknown',form='formatted',file='pspvariancesp.'//outnam)
      WRITE(75,'(a)') '# s.p. energy variances (1ph-projected)'
      CLOSE(75)
    ENDIF
    IF(jstboostinv /= 0) then
      OPEN(91,status='unknown',form='formatted',file='pspenergybi.'//outnam)
      WRITE(91,'(a)') '# Single particle energies -- boost invariant'
      CLOSE(91)
      OPEN(92,status='unknown',form='formatted',file='pspvariancesbi.'//outnam)
      WRITE(92,'(a)') '# s.p. energy variances -- boost invariant'
      CLOSE(92)
    ENDIF
    
    IF(jcharges /= 0) THEN
      OPEN(323,STATUS='unknown',FILE='pcharges.'//outnam)
      WRITE(323,'(a,2i6)') '# Column 1 : time [fs]'
      WRITE(323,'(a,2i6)') '# Column 2 : total nr of electrons in box'
      DO ii=1,INT(nzsh*dz/drcharges)
        WRITE(323,'(a,i6,f8.3)') 'Column ',ii+2,ii*drcharges
      END DO
      CLOSE(323)
    END IF
    
    IF(iangabso /= 0 .AND. jangabso /= 0) THEN
      OPEN(47,STATUS='unknown',FORM='formatted', FILE='pangabso.'//outnam)
      WRITE(47,'(a)') ' & '
      CLOSE(47)
    END IF
    
    IF(ionmdtyp > 0) THEN
      
      IF(jpos /= 0) THEN
! Positions of cluster ions
        IF(nion > 0)THEN
          OPEN(21,STATUS='unknown',FORM='formatted', FILE='pposion.'//outnam)
          WRITE(21,'(a)') ' & '
          CLOSE(21)
          OPEN(621,STATUS='unknown',FILE='pgeomion.'//outnam)
          WRITE(621,'(a)') ' & '
          CLOSE(621)
        END IF
        
! Positions of GSM cores, clouds and kations
        IF(isurf /= 0) THEN
          OPEN(24,STATUS='unknown',FORM='formatted', FILE='pposcore.'//outnam)
          WRITE(24,'(a)') ' & '
          CLOSE(24)
          OPEN(25,STATUS='unknown',FORM='formatted', FILE='pposdip.'//outnam)
          WRITE(25,'(a)') ' & '
          CLOSE(25)
          OPEN(125,STATUS='unknown',FORM='formatted', FILE='pposkat.'//outnam)
          WRITE(125,'(a)') ' & '
          CLOSE(125)
        END IF
        
        OPEN(424,STATUS='unknown',FORM='formatted', FILE='penerSurf.'//outnam)
        WRITE(424,'(a)') ' & '
        WRITE(424,*) '#   tfs,ekinC,ekinE,ekinKf, ekinSurf'
        CLOSE(424)
        
      END IF
      
! Positions of center of mass
      IF(jposcm /= 0) THEN
        OPEN(281,STATUS='unknown',FORM='formatted', FILE='pposCM.'//outnam)
        WRITE(281,'(a)') ' & '
        CLOSE(281)
      END IF
      
      IF(jvel /= 0) THEN
        IF(nion > 0)THEN
! Velocities of cluster ions
          OPEN(22,STATUS='unknown',FORM='formatted', FILE='pvelion.'//outnam)
          WRITE(22,'(a)') ' & '
          CLOSE(22)
          OPEN(148,STATUS='unknown',FORM='formatted',  &
              FILE='pkinenion.'//outnam)
          WRITE(148,'(a)') ' & '
          CLOSE(148)
          OPEN(148,STATUS='unknown',FORM='formatted',  &
              FILE='ptempion.'//outnam)
          WRITE(148,'(a)') ' & '
          CLOSE(148)
        END IF
        IF(isurf /= 0)THEN
! Velocities of substrate cores
          OPEN(26,STATUS='unknown',FORM='formatted', FILE='pvelcore.'//outnam)
          WRITE(26,'(a)') ' & '
          CLOSE(26)
          OPEN(126,STATUS='unknown',FORM='formatted', FILE='pvelkat.'//outnam)
          WRITE(126,'(a)') ' & '
          CLOSE(126)
          OPEN(148,STATUS='unknown',FORM='formatted',  &
              FILE='pkinencore.'//outnam)
          WRITE(148,'(a)') ' & '
          CLOSE(148)
          OPEN(148,STATUS='unknown',FORM='formatted',  &
              FILE='pkinenkat.'//outnam)
          WRITE(148,'(a)') ' & '
          CLOSE(148)
          OPEN(148,STATUS='unknown',FORM='formatted',  &
              FILE='ptempsurf.'//outnam)
          WRITE(148,'(a)') ' & '
          CLOSE(148)
!                                   Velocities of Argon clouds (if available)
          IF(ifadiadip /= 1) THEN
            OPEN(27,STATUS='unknown',FORM='formatted', FILE='pveldip.'//outnam)
            WRITE(27,'(a)') ' & '
            CLOSE(27)
          END IF

        END IF
      END IF
      
      IF(jforce /= 0) THEN
        
        IF(nion > 0)THEN
! Forces on cluster ions
          OPEN(30,STATUS='unknown',FORM='formatted', FILE='pforce.3.'//outnam)
          WRITE(30,'(a)') ' & '
          CLOSE(30)
          
          IF(e0 /= 0D0) THEN
! Forces from laser
            OPEN(31,STATUS='unknown',FORM='formatted',  &
                FILE='plforce.3.'//outnam)
            WRITE(31,'(a)') ' & '
            CLOSE(31)
          END IF
          
        END IF
        IF(isurf /= 0) THEN
          IF(nc+NE+nk > 0)THEN
! Forces on Argon cores
            OPEN(28,STATUS='unknown',FORM='formatted',  &
                FILE='pforce.1.'//outnam)
            WRITE(28,'(a)') ' & '
            CLOSE(28)
            
            IF(e0 /= 0D0) THEN
! Forces from laser
              OPEN(32,STATUS='unknown',FORM='formatted',  &
                  FILE='plforce.1.'//outnam)
              WRITE(32,'(a)') ' & '
              CLOSE(32)
            END IF
            
            IF(nclust > 0)THEN
! Forces on Argon clouds
! (else, there have no forces)
              OPEN(29,STATUS='unknown',FORM='formatted',  &
                  FILE='pforce.2.'//outnam)
              WRITE(29,'(a)') ' & '
              CLOSE(29)
              
              IF(e0 /= 0D0) THEN
! Forces from laser
                OPEN(33,STATUS='unknown',FORM='formatted',  &
                    FILE='plforce.2.'//outnam)
                WRITE(33,'(a)') ' & '
                CLOSE(33)
              END IF
              
            END IF
          END IF
        END IF
      END IF
      IF(jener > 0) THEN
        
        IF(nion > 0)THEN
! Energies of cluster
          OPEN(34,STATUS='unknown',FORM='formatted', FILE='penerclu.'//outnam)
          WRITE(34,'(a)') ' & '
          CLOSE(34)
        END IF
        
        IF(nc+NE+nk > 0)THEN
! Energies of the matrix
          OPEN(35,STATUS='unknown',FORM='formatted', FILE='penermat.'//outnam)
          WRITE(35,'(a)') ' & '
          CLOSE(35)
        END IF
        
      END IF
    END IF
    
  END IF
  
  
END IF                 ! after jump over file headers


IF(nclust > 0)THEN
  IF(myn == 0)THEN
    WRITE(6,*)
    WRITE(7,*)
    WRITE(6,*) ' **** dynamical plasmon calculation - tdlda ****'
    WRITE(7,*) ' **** dynamical plasmon calculation - tdlda ****'
    WRITE(6,*)
    WRITE(7,*)
  END IF
ELSE
  IF(myn == 0)THEN
    WRITE(6,*)
    WRITE(7,*)
    WRITE(6,*)' **** molecular dynamic calculation ****'
    WRITE(7,*)' **** molecular dynamic calculation ****'
    WRITE(6,*)
    WRITE(7,*)
    
#if(raregas)
    IF (isurf /= 0) THEN
      CALL adjustdip(rho)
      DO i=1,NE
        xeinit(i)=xe(i)
        yeinit(i)=ye(i)
        zeinit(i)=ze(i)
      END DO
      CALL info(psi,rho,aloc,0)
      
      IF (myn == 0) THEN
        OPEN(308,POSITION='append',FILE='energies.res')
        dist = ABS(cz(1) - dmaximum(zc,nc))
        WRITE(308,'(1f20.10,2e25.14)') dist, energy-enerinfty,energy
        CLOSE(308)
      END IF
    END IF
#endif
    
    CALL getforces(rho,psi,0)
    
  END IF
END IF


RETURN
END SUBROUTINE init_dynprotocol



!-----print_densdiff--------------------------------------------------

SUBROUTINE print_densdiff(rho,it)

!     evaluation and print of density differences

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
INTEGER, INTENT(IN)                      :: it
REAL(DP),DIMENSION(:),ALLOCATABLE :: w1

CHARACTER (LEN=1) :: str1
CHARACTER (LEN=2) :: str2
CHARACTER (LEN=3) :: str3
CHARACTER (LEN=4) :: str4
CHARACTER (LEN=5) :: str5
CHARACTER (LEN=6) :: str6
CHARACTER (LEN=7) :: str7
CHARACTER (LEN=8) :: str8
CHARACTER (LEN=9) :: str9

!---------------------------------------------------------------------



IF (jplotdensitydiff /= 0 .AND. MOD(it,jplotdensitydiff) == 0) THEN
  
!  IF(usew1) STOP ' in PRINT_DENSDIFF:  W1 already active '
!  usew1 = .true.
  ALLOCATE(w1(kdfull2))
  OPEN(590,STATUS='unknown',FILE='densdiff')
  DO i=1,kdfull2
    READ(590,*) w1(i)
  END DO
  CLOSE(590)
  
  
  it0=it
  ndig=1
  
  DO i=1,10
    IF (it0/10D0 >= 1D0) THEN
      ndig=ndig+1
      it0=INT(it0/10)
    END IF
  END DO
  
  IF (ndig == 1) THEN
    WRITE (str1, '(I1)') it
    OPEN(689,STATUS='unknown',FILE='pdensdiff.'//str1//'.'//outnam)
  ELSE IF (ndig == 2) THEN
    WRITE (str2, '(I2)') it
    OPEN(689,STATUS='unknown',FILE='pdensdiff.'//str2//'.'//outnam)
  ELSE IF (ndig == 3) THEN
    WRITE (str3,'(I3)') it
    OPEN(689,STATUS='unknown',FILE='pdensdiff.'//str3//'.'//outnam)
  ELSE IF (ndig == 4) THEN
    WRITE (str4, '(I4)') it
    OPEN(689,STATUS='unknown',FILE='pdensdiff.'//str4//'.'//outnam)
  ELSE IF (ndig == 5) THEN
    WRITE (str5, '(I5)') it
    OPEN(689,STATUS='unknown',FILE='pdensdiff.'//str5//'.'//outnam)
  ELSE IF (ndig == 6) THEN
    WRITE (str6, '(I6)') it
    OPEN(689,STATUS='unknown',FILE='pdensdiff.'//str6//'.'//outnam)
  ELSE IF (ndig == 7) THEN
    WRITE (str7, '(I7)') it
    OPEN(689,STATUS='unknown',FILE='pdensdiff.'//str7//'.'//outnam)
  ELSE IF (ndig == 8) THEN
    WRITE (str8, '(I8)') it
    OPEN(689,STATUS='unknown',FILE='pdensdiff.'//str8//'.'//outnam)
  ELSE IF (ndig == 9) THEN
    WRITE (str9, '(I9)') it
    OPEN(689,STATUS='unknown',FILE='pdensdiff.'//str9//'.'//outnam)
  ELSE
    STOP '::too many time steps::'
  END IF
!  CALL addfields2(w1,rho,-1D0)
  w1=w1-rho
  CALL printfield(689,w1,'x')
  CLOSE(689)
!  usew1 = .false.
  DEALLOCATE(w1)  
  
END IF

!----------------------------------------
IF (jplotdensitydiff2d /= 0 .AND. MOD(it,jplotdensitydiff2d) == 0) THEN
  
!  usew1 = .true.
  ALLOCATE(w1(kdfull2))
  OPEN(590,STATUS='unknown',FILE='densdiff')
  DO i=1,kdfull2*2
    READ(590,*) w1(i)
  END DO
  CLOSE(590)
  it0=it
  ndig=1
  
  DO i=1,10
    IF (it0/10D0 >= 1D0) THEN
      ndig=ndig+1
      it0=INT(it0/10)
    END IF
  END DO
  
  IF (ndig == 1) THEN
    WRITE (str1, '(I1)') it
!    CALL addfields2(w1,rho,-1D0)
    w1=w1-rho
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dxy.'//str1//'.'//outnam)
    CALL pm3dcut(689,1,2,0D0,w1)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dxz.'//str1//'.'//outnam)
    CALL pm3dcut(689,1,3,0D0,w1)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dyz.'//str1//'.'//outnam)
    CALL pm3dcut(689,2,3,0D0,w1)
    CLOSE(689)
  ELSE IF (ndig == 2) THEN
    WRITE (str2, '(I2)') it
!    CALL addfields2(w1,rho,-1D0)
    w1=w1-rho
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dxy.'//str2//'.'//outnam)
    CALL pm3dcut(689,1,2,0D0,w1)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dxz.'//str2//'.'//outnam)
    CALL pm3dcut(689,1,3,0D0,w1)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dyz.'//str2//'.'//outnam)
    CALL pm3dcut(689,2,3,0D0,w1)
    CLOSE(689)
  ELSE IF (ndig == 3) THEN
    WRITE (str3,'(I3)') it
!    CALL addfields2(w1,rho,-1D0)
    w1=w1-rho
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dxy.'//str3//'.'//outnam)
    CALL pm3dcut(689,1,2,0D0,w1)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dxz.'//str3//'.'//outnam)
    CALL pm3dcut(689,1,3,0D0,w1)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dyz.'//str3//'.'//outnam)
    CALL pm3dcut(689,2,3,0D0,w1)
    CLOSE(689)
  ELSE IF (ndig == 4) THEN
    WRITE (str4, '(I4)') it
!    CALL addfields2(w1,rho,-1D0)
    w1=w1-rho
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dxy.'//str4//'.'//outnam)
    CALL pm3dcut(689,1,2,0D0,w1)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dxz.'//str4//'.'//outnam)
    CALL pm3dcut(689,1,3,0D0,w1)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dyz.'//str4//'.'//outnam)
    CALL pm3dcut(689,2,3,0D0,w1)
    CLOSE(689)
  ELSE IF (ndig == 5) THEN
    WRITE (str5, '(I5)') it
!    CALL addfields2(w1,rho,-1D0)
    w1=w1-rho
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dxy.'//str5//'.'//outnam)
    CALL pm3dcut(689,1,2,0D0,w1)
    CLOSE(689)
    WRITE (str5, '(I5)') it
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dxz.'//str5//'.'//outnam)
    CALL pm3dcut(689,1,3,0D0,w1)
    CLOSE(689)
    WRITE (str5, '(I5)') it
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dyz.'//str5//'.'//outnam)
    CALL pm3dcut(689,2,3,0D0,w1)
    CLOSE(689)
  ELSE IF (ndig == 6) THEN
    WRITE (str6, '(I6)') it
!    CALL addfields2(w1,rho,-1D0)
    w1=w1-rho
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dxy.'//str6//'.'//outnam)
    CALL pm3dcut(689,1,2,0D0,w1)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dxz.'//str6//'.'//outnam)
    CALL pm3dcut(689,1,3,0D0,w1)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dyz.'//str6//'.'//outnam)
    CALL pm3dcut(689,2,3,0D0,w1)
    CLOSE(689)
  ELSE IF (ndig == 7) THEN
    WRITE (str7, '(I7)') it
!    CALL addfields2(w1,rho,-1D0)
    w1=w1-rho
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dxy.'//str7//'.'//outnam)
    CALL pm3dcut(689,1,2,0D0,w1)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dxz.'//str7//'.'//outnam)
    CALL pm3dcut(689,1,3,0D0,w1)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dyz.'//str7//'.'//outnam)
    CALL pm3dcut(689,2,3,0D0,w1)
    CLOSE(689)
  ELSE IF (ndig == 8) THEN
    WRITE (str8, '(I8)') it
!    CALL addfields2(w1,rho,-1D0)
    w1=w1-rho
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dxy.'//str8//'.'//outnam)
    CALL pm3dcut(689,1,2,0D0,w1)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dxz.'//str8//'.'//outnam)
    CALL pm3dcut(689,1,3,0D0,w1)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dyz.'//str8//'.'//outnam)
    CALL pm3dcut(689,2,3,0D0,w1)
    CLOSE(689)
  ELSE IF (ndig == 9) THEN
    WRITE (str9, '(I9)') it
!    CALL addfields2(w1,rho,-1D0)
    w1=w1-rho
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dxy.'//str9//'.'//outnam)
    CALL pm3dcut(689,1,2,0D0,w1)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dxz.'//str9//'.'//outnam)
    CALL pm3dcut(689,1,3,0D0,w1)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdensdiff2Dyz.'//str9//'.'//outnam)
    CALL pm3dcut(689,2,3,0D0,w1)
    CLOSE(689)
  ELSE
    STOP '::too many time steps::'
  END IF
  
!  usew1 = .false.
  DEALLOCATE(w1)
  
  
END IF


!----------------------------------------



IF (jplotdensity2d /= 0 .AND. MOD(it,jplotdensity2d) == 0) THEN
  
  it0=it
  ndig=1
  
  DO i=1,10
    IF (it0/10D0 >= 1D0) THEN
      ndig=ndig+1
      it0=INT(it0/10)
    END IF
  END DO
  
  IF (ndig == 1) THEN
    WRITE (str1, '(I1)') it
    OPEN(689,STATUS='unknown', FILE='pdens2Dxy.'//str1//'.'//outnam)
    CALL pm3dcut(689,1,2,0D0,rho)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdens2Dxz.'//str1//'.'//outnam)
    CALL pm3dcut(689,1,3,0D0,rho)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdens2Dyz.'//str1//'.'//outnam)
    CALL pm3dcut(689,2,3,0D0,rho)
    CLOSE(689)
  ELSE IF (ndig == 2) THEN
    WRITE (str2, '(I2)') it
    OPEN(689,STATUS='unknown', FILE='pdens2Dxy.'//str2//'.'//outnam)
    CALL pm3dcut(689,1,2,0D0,rho)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdens2Dxz.'//str2//'.'//outnam)
    CALL pm3dcut(689,1,3,0D0,rho)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdens2Dyz.'//str2//'.'//outnam)
    CALL pm3dcut(689,2,3,0D0,rho)
    CLOSE(689)
  ELSE IF (ndig == 3) THEN
    WRITE (str3,'(I3)') it
    OPEN(689,STATUS='unknown', FILE='pdens2Dxy.'//str3//'.'//outnam)
    CALL pm3dcut(689,1,2,0D0,rho)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdens2Dxz.'//str3//'.'//outnam)
    CALL pm3dcut(689,1,3,0D0,rho)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdens2Dyz.'//str3//'.'//outnam)
    CALL pm3dcut(689,2,3,0D0,rho)
    CLOSE(689)
  ELSE IF (ndig == 4) THEN
    WRITE (str4, '(I4)') it
    OPEN(689,STATUS='unknown', FILE='pdens2Dxy.'//str4//'.'//outnam)
    CALL pm3dcut(689,1,2,0D0,rho)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdens2Dxz.'//str4//'.'//outnam)
    CALL pm3dcut(689,1,3,0D0,rho)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdens2Dyz.'//str4//'.'//outnam)
    CALL pm3dcut(689,2,3,0D0,rho)
    CLOSE(689)
  ELSE IF (ndig == 5) THEN
    WRITE (str5, '(I5)') it
    OPEN(689,STATUS='unknown', FILE='pdens2Dxy.'//str5//'.'//outnam)
    CALL pm3dcut(689,1,2,0D0,rho)
    CLOSE(689)
    WRITE (str5, '(I5)') it
    OPEN(689,STATUS='unknown', FILE='pdens2Dxz.'//str5//'.'//outnam)
    CALL pm3dcut(689,1,3,0D0,rho)
    CLOSE(689)
    WRITE (str5, '(I5)') it
    OPEN(689,STATUS='unknown', FILE='pdens2Dyz.'//str5//'.'//outnam)
    CALL pm3dcut(689,2,3,0D0,rho)
    CLOSE(689)
  ELSE IF (ndig == 6) THEN
    WRITE (str6, '(I6)') it
    OPEN(689,STATUS='unknown', FILE='pdens2Dxy.'//str6//'.'//outnam)
    CALL pm3dcut(689,1,2,0D0,rho)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdens2Dxz.'//str6//'.'//outnam)
    CALL pm3dcut(689,1,3,0D0,rho)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdens2Dyz.'//str6//'.'//outnam)
    CALL pm3dcut(689,2,3,0D0,rho)
    CLOSE(689)
  ELSE IF (ndig == 7) THEN
    WRITE (str7, '(I7)') it
    OPEN(689,STATUS='unknown', FILE='pdens2Dxy.'//str7//'.'//outnam)
    CALL pm3dcut(689,1,2,0D0,rho)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdens2Dxz.'//str7//'.'//outnam)
    CALL pm3dcut(689,1,3,0D0,rho)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdens2Dyz.'//str7//'.'//outnam)
    CALL pm3dcut(689,2,3,0D0,rho)
    CLOSE(689)
  ELSE IF (ndig == 8) THEN
    WRITE (str8, '(I8)') it
    OPEN(689,STATUS='unknown', FILE='pdens2Dxy.'//str8//'.'//outnam)
    CALL pm3dcut(689,1,2,0D0,rho)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdens2Dxz.'//str8//'.'//outnam)
    CALL pm3dcut(689,1,3,0D0,rho)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdens2Dyz.'//str8//'.'//outnam)
    CALL pm3dcut(689,2,3,0D0,rho)
    CLOSE(689)
  ELSE IF (ndig == 9) THEN
    WRITE (str9, '(I9)') it
    OPEN(689,STATUS='unknown', FILE='pdens2Dxy.'//str9//'.'//outnam)
    CALL pm3dcut(689,1,2,0D0,rho)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdens2Dxz.'//str9//'.'//outnam)
    CALL pm3dcut(689,1,3,0D0,rho)
    CLOSE(689)
    OPEN(689,STATUS='unknown', FILE='pdens2Dyz.'//str9//'.'//outnam)
    CALL pm3dcut(689,2,3,0D0,rho)
    CLOSE(689)
  ELSE
    STOP '::too many time steps::'
  END IF
  
  
  
END IF

RETURN
END SUBROUTINE print_densdiff



!-----open_protok_el----------------------------------------------

SUBROUTINE open_protok_el(it)

!     open protocol files for electronic properties

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

!----------------------------------------------------------------


IF (e0 /= 0) CALL safeopen(38,0,1,'plaser')
CALL safeopen(163,it,jenergy,'penergies')
CALL safeopen(323,it,jcharges,'pcharges')

#if(simpara)
IF(nclust > 0)THEN
#else
IF(nclust > 0 .AND. myn == 0)THEN
#endif
  
  CALL safeopen(8,it,jdip,'pdip')
  CALL safeopen(810,it,jdiporb,'pdiporb.x')
  CALL safeopen(811,it,jdiporb,'pdiporb.y')
  CALL safeopen(812,it,jdiporb,'pdiporb.z')
  CALL safeopen(9,it,jquad,'pquad')
  CALL safeopen(17,it,jinfo,'infosp')
  CALL safeopen(47,0,jangabso,'pangabso')
  CALL safeopen(68,it,jang,'pangmo')
  CALL safeopen(78,it,jspdp,'pspdip')
  CALL safeopen(608,it,jgeomel,'pgeomel')
  CALL safeopen(806,it,jnorms,'pescOrb')
  CALL safeopen(808,it,jnorms,'pproba')
!  CALL safeopen(803,it,jmp,'pMP')
  
END IF
  

RETURN
END SUBROUTINE open_protok_el

!-----analyze_ions-----------------------------------------------

SUBROUTINE analyze_ions(it)

!     open protocol files and print properties of ions

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

!----------------------------------------------------------------

!      tfs=it*dt1*0.0484

CALL getcm(1,0,0,0)                  !  c.m. now on 'rvectmp(1:3)'
IF (jposcm > 0 .AND. MOD(it,jposcm) == 0) THEN
  CALL  safeopen(281,it,jposcm,'pposCM')
!  OPEN(281,POSITION='append',FILE='pposCM.'//outnam)
  WRITE(281,'(1f15.6,3e17.8)') tfs, rvectmp(1), rvectmp(2),rvectmp(3)
  CLOSE(281)
!        call flush(281)
END IF

IF(((jpos > 0 .AND. MOD(it,jpos) == 0)  &
      .OR.(jvel > 0 .AND. MOD(it,jvel) == 0)  &
      .OR.(jforce /= 0 .AND. MOD(it,jforce) == 0)))THEN
  CALL  safeopen(21,it,jpos,'pposion')
!         open(21,position='append',file='pposion.'//outnam)
  CALL  safeopen(621,it,jgeomion,'pgeomion')
!         open(621,position='append',file='pgeomion.'//outnam)
  CALL  safeopen(22,it,jvel,'pvelion')
!         open(22,position='append',file='pvelion.'//outnam)
  CALL  safeopen(149,it,jvel,'pkinenion')
!         open(149,position='append',file='pkinenion.'//outnam)
  CALL  safeopen(151,it,jvel,'ptempion')
!         open(151,position='append',file='ptempion.'//outnam)

  r2ion = SQRT(SUM( (cx(1:nion)-rvectmp(1))**2 +(cy(1:nion)-rvectmp(2))**2 &
                   +(cz(1:nion)-rvectmp(3))**2)/nion)
  sumx=0D0
  sumy=0D0
  sumz=0D0
  DO ion=1,nion
    IF(MOD(it,jpos) == 0) WRITE(21,'(1f13.5,3e17.8,1pg13.5)')  &
        tfs,cx(ion),cy(ion),cz(ion),r2ion   !  ecorr
    IF(MOD(it,jvel) == 0) WRITE(22,'(1f13.5,3e17.8,1pg13.5)')  &
        tfs,cpx(ion),cpy(ion),cpz(ion),ekion
    sumx = sumx + (cpx(ion)**2)/amu(np(ion))/1836.
    sumy = sumy + cpy(ion)**2/amu(np(ion))/1836.
    sumz = sumz + cpz(ion)**2/amu(np(ion))/1836.
  END DO
  IF(MOD(it,jvel) == 0) THEN
    WRITE(149,'(1f13.5,4e17.8,i5)') tfs,sumx,sumy,sumz,sumx+sumy+sumz,nion
    CALL gettemperature(4)
  END IF
  
  IF(jgeomion > 0 .AND. MOD(it,jgeomion) == 0) THEN
    CALL getclustergeometry
    WRITE(621,'(12e15.5)') tfs,comx,comy,comz,rmsion,  &
        qtion(1,1),qtion(2,2),qtion(3,3),  &
        qtion(1,2),qtion(1,3),qtion(2,3),dmdistion
  END IF
  
  CALL flush(21)
  CALL flush(22)
  IF(MOD(it,jpos) == 0) CALL view3d()
END IF


!       Kinetic energy in the three directions for cluster and matrix
!               and total energy (both in penerclu and penermat)

IF(jener > 0 .AND. MOD(it,jener) == 0 .AND. (nion) > 0)THEN
! cluster
  ekx = 0D0
  eky = 0D0
  ekz = 0D0
  amfac = amu(np(nion))*1836.*ame*2.
  DO ion=1,nion
    ekx = ekx + cpx(ion)*cpx(ion)
    eky = eky + cpy(ion)*cpy(ion)
    ekz = ekz + cpz(ion)*cpz(ion)
  END DO
  ekx = ekx/amfac
  eky = eky/amfac
  ekz = ekz/amfac
  CALL  safeopen(34,it,jener,'penerclu')
!         open(34,position='append',file='penerclu.'//outnam)
  WRITE(34,'(5f13.5)')tfs,ekx,eky,ekz,etot
  CALL flush(34)
END IF

#if(raregas)
!    ?? not clear what that is doing

IF (jpos> 0 .AND. MOD(it,jpos) == 0 .AND. nclust == 0) THEN
!???            call info(psi,rho,aloc,akv,it)
  sumnc=0D0
  DO i=1,nc
    dist2 = getdistance2(i,i+nc)
    sumnc=sumnc+0.5D0*cspr*dist2
  END DO
  WRITE(773,'(1f12.5,1e17.7)') tfs,sumnc
END IF
#endif

RETURN
END SUBROUTINE analyze_ions

!-----analyze_surf------------------------------------------------

SUBROUTINE analyze_surf(it)

!     open protocol files and print properties of substrate

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

!----------------------------------------------------------------

!      tfs=it*dt1*0.0484
#if(raregas)
IF(nc+NE+nk > 0) THEN
  CALL safeopen(24,it,jpos,'pposcore')
  CALL safeopen(157,it,jpos,'ptempsurf')
  CALL safeopen(25,it,jpos,'pposdip')
  CALL safeopen(125,it,jpos,'pposkat')
  CALL safeopen(424,it,jpos,'penerSurf')
  CALL safeopen(26,it,jvel,'pvelcore')
  CALL safeopen(126,it,jvel,'pvelkat')
  IF(ifadiadip /= 1) CALL safeopen(27,it,jvel,'pveldip')
  CALL safeopen(152,it,jvel,'pkinencore')
  CALL safeopen(153,it,jvel,'pkinenkat')
  IF(((jpos > 0 .AND. MOD(it,jpos) == 0)  &
        .OR.(jvel > 0 .AND. MOD(it,jvel) == 0)))THEN
    
    IF (MOD(it,jpos) == 0) THEN
      CALL getsurfprops
      WRITE(424,'(f12.4,4e15.5)') tfs,ekincsurf,ekinesurf, ekinksurf, ekinsurf
    END IF
    
    IF (MOD(it,jsavesurf) == 0) THEN
      OPEN(464,STATUS='unknown',FILE='for006surf.'//outnam)
      DO i=1,nc
        WRITE(464,'(6e17.7,i6)') xc(i),yc(i),zc(i),  &
            xe(i),ye(i),ze(i),imobc(i),imobe(i)
      END DO
      DO i=1,nk
        WRITE(464,'(3e17.7,i6)') xk(i),yk(i),zk(i),imobk(i)
      END DO
      CLOSE(464)
    END IF
    
    
    sumcx=0D0
    sumcy=0D0
    sumcz=0D0
    nimobc=0
    DO ion=1,nc
      IF(MOD(it,jpos) == 0) THEN
        IF (iprintonlyifmob == 0 .OR. imobc(ion) /= 0)  &
            WRITE(24,'(1f13.5,3e17.8)') tfs,xc(ion),yc(ion),zc(ion)
        IF (iprintonlyifmob == 0 .OR. imobe(ion) /= 0)  &
            WRITE(25,'(1f13.5,3e17.8)') tfs,xe(ion),ye(ion),ze(ion)
      END IF
      
      IF(MOD(it,jvel) == 0) THEN
        IF (iprintonlyifmob == 0 .OR. imobc(ion) /= 0)  &
            WRITE(26,'(1f13.5,3e17.8)') tfs,pxc(ion),pyc(ion),pzc(ion)
        IF(ifadiadip /= 1 .AND. &
           (iprintonlyifmob == 0 .OR. imobe(ion) /= 0)) THEN
             WRITE(27,'(1f13.5,3e17.8)') tfs,pxe(ion),pye(ion),pze(ion)
             nimobc=nimobc+1
        END IF
      END IF
      
      IF (imobc(ion) /= 0) THEN
        sumcx = sumcx + (pxc(ion)**2)/mion/1836.
        sumcy = sumcy + (pyc(ion)**2)/mion/1836.
        sumcz = sumcx + (pzc(ion)**2)/mion/1836.
      END IF
      
    END DO
    
    WRITE(6,*) mion,mkat
    IF(MOD(it,jvel) == 0) THEN
      WRITE(152,'(1f13.5,4e17.8,i6)')  &
          tfs,sumcx,sumcy,sumcz,sumcx+sumcy+sumcz,nimobc
      CALL gettemperature(1)
    END IF
    
    sumkx=0D0
    sumky=0D0
    sumkz=0D0
    nimobk=0
    DO ion=1,nk
      IF(MOD(it,jpos) == 0 .AND.  &
          (iprintonlyifmob == 0 .OR. imobk(ion) /= 0))  &
          WRITE(125,'(1f13.5,3e17.8)') tfs,xk(ion),yk(ion),zk(ion)
      IF(MOD(it,jvel) == 0 .AND.  &
          (iprintonlyifmob == 0 .OR. imobk(ion) /= 0)) THEN
      WRITE(126,'(1f13.5,3e17.8)') tfs,pxk(ion),pyk(ion),pzk(ion)
      nimobk=nimobk+1
    END IF
    
    IF (imobk(ion) /= 0) THEN
      sumkx = sumkx + (pxk(ion)**2)/mkat/1836.
      sumky = sumky + (pyk(ion)**2)/mkat/1836.
      sumkz = sumkx + (pzk(ion)**2)/mkat/1836.
    END IF
  END DO
  
  IF(MOD(it,jvel) == 0) WRITE(153,'(1f13.5,4e17.8,i6)')  &
      tfs,sumkx,sumky,sumkz,nimobk
  
  CALL flush(24)
  CALL flush(25)
  CALL flush(26)
  IF(ifadiadip /= 1)   CALL flush(27)
END IF
END IF

IF(myn == 0 .AND. jener > 0 .AND. MOD(it,jener) == 0 .AND. nrare > 0 )THEN
! matrix
  ekcx = 0D0
  ekcy = 0D0
  ekcz = 0D0
  ekvx = 0D0
  ekvy = 0D0
  ekvz = 0D0
  amfac1 = amu(np(1))*1836.*ame*2.
  amfac2 = amu(np(nrare+1))*1836.*ame*2.
  DO ion=1,nrare
    ekcx = ekcx + pxc(ion)*pxc(ion)
    ekcy = ekcy + pyc(ion)*pyc(ion)
    ekcz = ekcz + pzc(ion)*pzc(ion)
    ekvx = ekvx + pxe(ion)*pxe(ion)
    ekvy = ekvy + pye(ion)*pye(ion)
    ekvz = ekvz + pze(ion)*pze(ion)
  END DO
  ekx = ekcx/amfac1 + ekvx/amfac2
  eky = ekcy/amfac1 + ekvy/amfac2
  ekz = ekcz/amfac1 + ekvz/amfac2
  CALL safeopen(35,it,jener,'penermat')
  WRITE(35,'(4f13.5)')tfs,ekx,eky,ekz,etot
  CALL flush(35)
END IF

IF(isurf == 1.AND.iforcecl2co /= 0 .AND.MOD(it,iforcecl2co) == 0) THEN
  
  CALL getforces_clust2cores(rho,psi,0)
  
  OPEN(256,STATUS='unknown',FILE='force_clust2cores.'//outnam)
  WRITE(256,'(a,f10.5)')'core pos val pos force on cores at t=' ,tfs
  
  DO ii=1,nc
    
    WRITE(256,'(9e20.10)') xc(ii),yc(ii),zc(ii),xe(ii),ye(ii),  &
        ze(ii),fxc(ii),fyc(ii),fzc(ii)
    
  END DO
  
  CLOSE(256)
END IF
#endif

RETURN
END SUBROUTINE analyze_surf


!-----analyze_elect-----------------------------------------------

SUBROUTINE analyze_elect(psi,rho,aloc,it)

!     Analysis and print of electronic properties during dynamics.

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)


COMPLEX(DP), INTENT(IN)                      :: psi(kdfull2,kstate)
REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)
!REAL(DP), INTENT(IN OUT)                     :: akv(kdfull2)
INTEGER, INTENT(OUT)                     :: it


LOGICAL,PARAMETER :: ttest=.FALSE.

COMPLEX(DP) :: orbitaloverlap
REAL(DP) :: rtmp(kstate,kstate)
COMPLEX(DP) :: cscal

!----------------------------------------------------------------

!      tfs=it*dt1*0.0484

IF(ABS(phangle) > small) CALL phoverl(psi)
!CALL phoverl(psi)
IF(jstateoverlap == 1) THEN
  CALL stateoverl(psi,psisavex)
!  DEALLOCATE(psisavex)
END IF



!     ***** dynamic density plot *****

IF(idenspl > 0 .AND. MOD(it,idenspl) == 0) CALL fmtv_fld(psi,rho,it)


!       Number of escaped electrons

IF(myn == 0 .AND. jesc > 0 .AND. MOD(it,jesc) == 0) THEN
  CALL nescape(it,rho)
END IF

IF (jcharges /= 0 .AND. MOD(it,jcharges) == 0) THEN
  CALL calcchargdist(it,rho)
END IF



!     Time-dependent Electron Localization Function (TD-ELF)

IF(jelf > 0 .AND. MOD(it,jelf) == 0) THEN
#if(!parayes) 
  CALL localize(rho,psi)
#else
  STOP ' LOCALIZE (switch JELF) should not be invoked in parallele code'
#endif
END IF


!     everything about single particles


IF((jstinf > 0 .AND. MOD(it,jstinf) == 0) &
   .OR. (jinfo > 0 .AND. MOD(it,jinfo)==0) &
   .OR. (jenergy > 0 .AND. MOD(it,jenergy)==0) & 
   .OR. (jesc > 0 .AND. jnorms>0 .AND. MOD(it,jnorms) == 0) &
   .OR. (jdiporb > 0 .AND. MOD(it,jdiporb)==0)) THEN
#if(parayes) 
  IF(ttest) WRITE(*,*) ' ANALYZE before INFO: myn=',myn
#endif
#if(parano)
  IF(ttest) WRITE(6,'(a,i4)') ' INFO from ANALYZE. it=',it
#endif
  CALL info(psi,rho,aloc,it)
#if(parayes) 
  IF(ttest) WRITE(*,*) ' ANALYZE after INFO: myn=',myn
#endif
END IF


!     escaped electrons analyzed orbital by orbital

#if(!parayes)
IF(jesc > 0 .AND. jnorms>0 .AND. MOD(it,jnorms) == 0) THEN
!  DO i=1,nstate
!    cscal=orbitaloverlap(psi(1,i),psi(1,i))
!    rtmp(i,1)=REAL(cscal)**2+imag(cscal)**2
!    rtmp(i,1)=1D0-SQRT(rtmp(i,1))
!  END DO
!call info(psi,rho,aloc,it)      !  move print 806 to 'pri_spe...'
  WRITE(806,'(500f12.8)') tfs,1D0-spnorm(1:nstate)
  CALL probab(psi)
END IF
#endif



IF(iangmo==1 .AND. iexcit == 1)  CALL instit (psi) !notes



IF(myn==0) THEN
  IF(MOD(it,100) == 0)THEN
    WRITE(6,'(a,i6,a,f12.6,a,f16.5)')  &
        'iter=',it,' tfs=',tfs,' total energy=',etot
  ELSE IF(nclust > 0)THEN
    WRITE(6,'(a,i6,a,f12.6,a,f16.5)')  &
        'iter=',it,' tfs=',tfs,' total energy=',etot
  END IF

  IF(jdip > 0 .AND. MOD(it,jdip) == 0) THEN
    WRITE(8,'(f10.5,3e17.8)') tfs,qe(2),qe(3),qe(4)
    CALL flush(8)
  END IF

  IF(jdiporb > 0 .AND. MOD(it,jdiporb) == 0) THEN
    WRITE(810,'(1000e17.8)') tfs,(qeorb_all(nbe,3),nbe=1,nstate_all)
    WRITE(811,'(1000e17.8)') tfs,(qeorb_all(nbe,4),nbe=1,nstate_all)
    WRITE(812,'(1000e17.8)') tfs,(qeorb_all(nbe,5),nbe=1,nstate_all)
    CALL flush(810)
    CALL flush(811)
    CALL flush(812)
  END IF
  
  IF(jquad > 0 .AND. MOD(it,jquad) == 0) THEN
    WRITE(9,'(f9.4,6f11.4)') tfs,qe(5),qe(6),qe(7), qe(8),qe(9),qe(10)
    CALL flush(9)
  END IF
END IF

IF(jmp > 0 .AND. MOD(it,jmp) == 0) CALL evalmp(803,psi)

IF(myn == 0 .AND. jangabso > 0 .AND. MOD(it,jangabso) == 0  &
    .AND. iangabso /= 0 ) CALL angabso(it)

IF(jgeomel > 0 .AND. MOD(it,jgeomel) == 0) THEN
  CALL getelectrongeometry(psi,0)
  WRITE(608,'(11e15.5)') tfs,codx,cody,codz,rmsel,  &
      qtel(1,1),qtel(2,2),qtel(3,3), qtel(1,2),qtel(1,3),qtel(2,3)
END IF

IF(jang > 0 .AND. MOD(it,jang) == 0)  THEN
  WRITE(68,'(f12.5,3f11.5)') tfs,ajx,ajy,ajz
  CALL flush(68)
END IF

IF(myn==0) THEN

  IF(jspdp > 0 .AND. MOD(it,jspdp) == 0 ) THEN
    WRITE(78,'(f10.5,3f13.5)') tfs,se(1),se(2),se(3)
    CALL flush(78)
  END IF

  WRITE(7,'(a,f8.4,a,f7.2,a,3f9.2)')  &
    't=',tfs,' moments: monop.=',qe(1),' dip.: ', qe(2),qe(3),qe(4)
  WRITE(7,'(a,i5,a/)') '** end of iter= ',it,'**'

  WRITE(6,'(a,f8.4,a,f7.2,a,3f9.2)')  &
    't=',tfs,' moments: monop.=',qe(1),' dip.: ', qe(2),qe(3),qe(4)

END IF

RETURN
END SUBROUTINE analyze_elect


!-----savings----------------------------------------------------

SUBROUTINE savings(psi,tarray,it)

!     Check status and optionally save wavefunctions

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

#if(simpara||parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)
INTEGER, INTENT(IN)                      :: it

REAL(4) :: tarray(2)

REAL(DP) :: trun

!----------------------------------------------------------------



IF(isave > 0 .AND. it /= 0) THEN
 IF(MOD(it,isave) == 0 .OR. it == itmax) THEN
  IF (irest /= 0 .AND. ABS(it-irest) <= 2) THEN
! do nothing, change later if needed
  ELSE
    CALL SAVE(psi,it,outnam)
  END IF
 END IF
END IF

!     if walltime has almost expired, save all relevant data to
!     prepare for restart and stop:

!      iiii=etime(tarray,trun)
IF(trequest > 0D0) THEN
 IF(myn == 0) THEN                                                        ! cPW
  !iiii=etime(tarray)
  CALL cpu_time(time_act)
  time_elapse=time_act-time_absinit
  !WRITE(6,'(a,5(1pg13.5))') &
  ! ' etime: trequest,timefrac,tarray=',trequest,timefrac,tarray,time_elapse
  CALL FLUSH(6)
  !IF ((tarray(1)+tarray(2)) > trequest*timefrac) THEN
 END IF
#if(parayes)
 CALL pi_scatter(time_elapse)
#endif
 IF (time_elapse > trequest*timefrac) THEN
  CALL SAVE(psi,it,outnam)
  IF(myn == 0) THEN
    OPEN(660,STATUS='unknown',FILE='progstatus')
    WRITE(660,*) '0'
    CLOSE(660)
  END IF
#if(parayes)
  CALL mpi_barrier (mpi_comm_world, mpi_ierror)
  CALL mpi_finalize (icode)
#endif
  STOP ' finish at TREQUEST'
 END IF
END IF                                                                    ! cPW
!     if program has received user message to shut down program, make it so!
IF (ishutdown /= 0) THEN
  CALL SAVE(psi,it,outnam)
#if(simpara)
  CALL mpi_barrier (mpi_comm_world, mpi_ierror)
#endif
  STOP 'Terminated by user message.'
END IF

RETURN
END SUBROUTINE savings

#if(raregas)
!-----vstep--------------------------------------------------------

SUBROUTINE vstep(rho,psi,it,dt)

!     propagation of rare gas valence clouds by leapfrog method

USE params
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)         :: rho(2*kdfull2)
REAL(DP), INTENT(IN OUT)         :: psi(kdfull2,kstate)
INTEGER, INTENT(IN OUT)          :: it
REAL(DP), INTENT(IN OUT)         :: dt





!------------------------------------------------------------------

!     optionally freeze valence momenta

DO i=1,NE
  IF (imobe(i) == 0) THEN
    pxe(i)=0D0
    pye(i)=0D0
    pze(i)=0D0
  END IF
END DO


!     propagation of positions first

!      xm=amu(np(nrare+1))*1836.0*ame

xm=amu(-18)*1836.0*ame

!      call leapfr(cx(nrare+1),cy(nrare+1),cz(nrare+1),
!     &     cpx(nrare+1),cpy(nrare+1),cpz(nrare+1),dt,xm,nrare)

CALL leapfr(xe(1),ye(1),ze(1), pxe(1),pye(1),pze(1),dt,xm,NE,2)


!     update subgrids in case of pseudo-densities

IF(ipseudo == 1) THEN
  CALL updatesubgrids
END IF


!     then propagation of momenta

CALL getforces(rho,psi,0) ! forces on valences with new positions


CALL leapfr(pxe(1),pye(1),pze(1), fxe(1),fye(1),fze(1),dt,1D0,NE,2)

RETURN
END SUBROUTINE vstep
!-----vstep--------------------------------------------------------

SUBROUTINE vstepv(rho,psi,it,dt)

!     propagation of rare gas valence clouds by velocity Verlet

USE params
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)         :: rho(2*kdfull2)
REAL(DP), INTENT(IN)             :: psi(kdfull2,kstate)
INTEGER, INTENT(IN OUT)          :: it
REAL(DP), INTENT(IN OUT)         :: dt





!------------------------------------------------------------------

!     optionally freeze valence momenta

DO i=1,NE
  IF (imobe(i) == 0) THEN
    pxe(i)=0D0
    pye(i)=0D0
    pze(i)=0D0
  END IF
END DO


!     propagation of positions first

!      xm=amu(np(nrare+1))*1836.0*ame

xm=amu(-18)*1836.0*ame
CALL velverlet1(xe(1),ye(1),ze(1),pxe(1),pye(1),pze(1), &
                fxe(1),fye(1),fze(1),dt,xm,ne,2)

!     update subgrids in case of pseudo-densities

IF(ipseudo == 1) THEN
  CALL updatesubgrids
END IF

! forces on valences with new positions

CALL getforces(rho,psi,0) 

!     then propagation of momenta

CALL velverlet2(pxe(1),pye(1),pze(1),fxe(1),fye(1),fze(1),dt,ne,2)

RETURN
END SUBROUTINE vstepv
#endif


