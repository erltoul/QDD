#include "define.h"
 
!-----tinit-----------------------------------------------------------

SUBROUTINE tinit(psir,psi)

!     to reshuffle real wavefunctions from 'psir' to the
!     complex wavefunction array 'psi'.
!     note that 'psir' and 'psi' may occupy the same storage.
!     thus the reshuffling expands the field while stepping
!     from highest state downwards.
!     intermediate storage is taken from workspace 'w1'.

USE params
USE kinetic
#if(twostsic)
USE twost
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN)                         :: psir(kdfull2,kstate)
COMPLEX(DP), INTENT(OUT)                     :: psi(kdfull2,kstate)

REAL(DP),DIMENSION(:),ALLOCATABLE :: w1


COMPLEX(DP) :: cfac

!--------------------------------------------------------------------



!     check workspace

!IF(usew1) THEN
!  STOP ' in TINIT: workspace W1 already in use'
!ELSE
!  usew1 = .true.
!END IF

!     copy wavefunctions from real to complex while
!     the wavefunction fields 'psir' (real) and 'psi' (complex)
!     may be overlayed

ALLOCATE(w1(kdfull2))
DO nbe=nstate,1,-1
  IF(nbe > 1) THEN
    DO i=1,kdfull2
      psi(i,nbe)=CMPLX(psir(i,nbe),0D0,DP)
    END DO
  ELSE
    DO i=1,kdfull2
      w1(i)=psir(i,nbe)
    END DO
    DO i=1,kdfull2
      psi(i,nbe)=CMPLX(w1(i),0D0,DP)
    END DO
  END IF
END DO

!     reset workspace

!usew1 = .false.
DEALLOCATE(w1)


!     boost the wavefunctions

IF (ekin0pp > 0D0) CALL init_velocel(psi)



!     optionally add scattering electron


IF (iscatterelectron /= 0) CALL init_scattel(psi)


rvectmp(1)=1D0             ! ??? what for ?

!JM
#if(twostsic)
CALL init_vecs(vecs)
#endif
!JM

RETURN
END SUBROUTINE tinit


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
USE kinetic
!USE coulsolv
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

!-----init_scattel-------------------------------------------------

SUBROUTINE init_scattel(psi)

!     Adds one electron in scattering state as Gaussian wavepacket
!     with a certain velocity.
!     The bookkeeping fiedls are extended accoprdingly.


USE params
USE kinetic
!USE coulsolv
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


fac = SQRT(pi**1.5*scatterelectronw**3)
fac = 1D0/fac

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
      cfac = CMPLX(COS(arg),SIN(arg),DP)
      rr = (x1 - scatterelectronx)**2
      rr = rr + (y1 - scatterelectrony)**2
      rr = rr + (z1 - scatterelectronz)**2
      fr = fac*EXP(-rr/2./scatterelectronw**2)
      psi(ind,nstate) = CMPLX(fr,0D0,DP)
      psi(ind,nstate) = psi(ind,nstate)*cfac
    END DO
  END DO
END DO

RETURN
END SUBROUTINE init_scattel

!-----tstep---------------------------------------------------------

SUBROUTINE tstep(q0,aloc,ak,rho,it)

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
!     compile time switch 'fastpropag'.

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN OUT)                  :: q0(kdfull2,kstate)
REAL(DP), INTENT(IN)                         :: aloc(2*kdfull2)
COMPLEX(DP), INTENT(IN OUT)                  :: ak(kdfull2)
REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
INTEGER, INTENT(IN)                      :: it
COMPLEX(DP),DIMENSION(:),ALLOCATABLE :: q1
!COMPLEX(DP) :: q1(2*kdfull2)


!DIMENSION  rhotmp(2*kdfull2)
!                                   workspace for non-local PsP
!                                   (to be revised!)
COMPLEX(DP),DIMENSION(:),ALLOCATABLE :: q2
!COMPLEX(DP) :: q2(kdfull2)
!      complex q3(kdfull2)
!      complex q4(kdfull2)
!      complex q5(kdfull2)
!      complex q6(kdfull2)
COMPLEX(DP) :: rii,cfac

CHARACTER (LEN=1) :: inttostring1
CHARACTER (LEN=2) :: inttostring2

CHARACTER (LEN=1) :: str1
CHARACTER (LEN=2) :: str2
CHARACTER (LEN=3) :: str3
CHARACTER (LEN=4) :: str4
CHARACTER (LEN=5) :: str5
CHARACTER (LEN=6) :: str6
CHARACTER (LEN=7) :: str7
CHARACTER (LEN=8) :: str8
CHARACTER (LEN=9) :: str9

LOGICAL :: tenerg

!EQUIVALENCE (q1(1),w1(1)) ! occupies also w3

!GB
#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)

CALL  mpi_comm_rank(mpi_comm_world,myn,icode)
#else
myn = 0
#endif
!GB


!     check workspace

!old      usew1=.false.

!IF(usew1) STOP ' in SSTEP: workspace W1 already active '
!IF(usew2) STOP ' in SSTEP: workspace W2 already active '
!IF(usew3) STOP ' in SSTEP: workspace W3 already active '
!IF(usew4) STOP ' in SSTEP: workspace W4 already active '
!usew1 = .true.
!usew2 = .true.
!usew3 = .true.
!usew4 = .true.
ALLOCATE(q1(2*kdfull2))
ALLOCATE(q2(kdfull2))

CALL cpu_time(time_init)
IF (it > ntref) nabsorb = 0           ! is that the correct place?


itsub = MOD(it,ipasinf) + 1

ri = -dt1*0.5D0
dt = dt1*0.5D0





!     half time step in coordinate space
!     local phase field on workspace 'q1'

#if(fullspin)
IF (nspdw == 0) THEN
  nup = nxyz
ELSE
  nup = 2*nxyz
END IF
DO ind=1,nup
  pr=-dt*aloc(ind)
  qr=COS(pr)
  rr=SIN(pr)
  q1(ind)=CMPLX(qr,rr,DP)
END DO
#else
DO ind=1,nxyz
  pr=-dt*aloc(ind)
  q1(ind)=CMPLX(COS(pr),SIN(pr),DP)
END DO
#endif
DO nb=1,nstate
  ishift = (ispin(nrel2abs(nb))-1)*nxyz
  CALL cmult3d(q0(1,nb),q1(1+ishift))
END DO



!     half non-local step


!old         if(ipsptyp.eq.1 .and. nrow(np(1)).ge.1) then     !??  .ge. ??
IF(ipsptyp == 1) THEN
  DO nb=1,nstate
!old            ind=0
    tenerg = itsub == ipasinf
    CALL nonlocstep(q0(1,nb),q1,q2,dt,tenerg,nb,4)   ! 6
  END DO
END IF


!       one full time step for the kinetic energy

DO nb=1,nstate
#if(gridfft)
#if(fastpropag)
  CALL kinprop(q0(1,nb),q1)
#else
  CALL fftf(q0(1,nb),q1)
  CALL cmult3d(q1,ak)
  CALL fftback(q1,q0(1,nb))
#endif
  
#endif
#if(findiff|numerov)
  CALL d3mixpropag (q0(1,nb),dt1)
#endif
  
END DO


!old       tfs = tfs + (dt - dt1)*0.0484/(2.*ame)




!old#if(fastpropag)
!old      if(itsub.ne.ipasinf) then
!old
!old         if(nabsorb.gt.0) then
!old            iComeFromAbso=0
!old            if (iangabso.ne.0) call calcrho(rho,q0)
!old                                ! necessary for angular distribution
!old
!old            if (ispherAbso.eq.0) then
!old               call abso(q0,it)
!old            else
!old               call spherAbso(q0,it)
!old            endif
!old         endif
!old      endif
!old#endif


!     half non-local step

IF(ipsptyp == 1 .AND. nrow(np(1)) >= 1) THEN
  DO nb = 1,nstate
    tenerg = .false. !   itsub.eq.ipasinf
    CALL nonlocstep(q0(1,nb),q1,q2,dt,tenerg,nb,4)   ! 6
  END DO
END IF



!     release workspace to make way for 'calcrho'
!usew1 = .false.
!usew2 = .false.
!usew3 = .false.
!usew4 = .false.
DEALLOCATE(q1)
DEALLOCATE(q2)



!     new density and local potential
!     (this is already the density at the end of the step,
!      because it is unchanged by unitary potential step)
!     propagation of substrate dipoles is done in 'dyn_mfield'.

CALL dyn_mfield(rho,aloc,q0,dt)




!     re-open workspace

!IF(usew1) STOP ' in SSTEP: workspace W1 already active '
!IF(usew2) STOP ' in SSTEP: workspace W2 already active '
!IF(usew3) STOP ' in SSTEP: workspace W3 already active '
!IF(usew4) STOP ' in SSTEP: workspace W4 already active '
!usew1 = .true.
!usew2 = .true.
!usew3 = .true.
!usew4 = .true.
ALLOCATE(q1(2*kdfull2))
ALLOCATE(q2(kdfull2))

!     half time step in coordinate space:

#if(fullspin)
IF (nspdw == 0) THEN
  nup = nxyz
ELSE
  nup = 2*nxyz
END IF
DO ind=1,nup
  pr=-dt*aloc(ind)
  qr=COS(pr)
  rr=SIN(pr)
  q1(ind)=CMPLX(qr,rr,DP)
END DO
#else
DO ind=1,nxyz
  pr=-dt*aloc(ind)
  q1(ind)=CMPLX(COS(pr),SIN(pr),DP)
END DO
#endif
DO nb=1,nstate
  ishift = (ispin(nrel2abs(nb))-1)*nxyz
  CALL cmult3d(q0(1,nb),q1(1+ishift))
END DO





!     finally release workspace

!usew1 = .false.
!usew2 = .false.
!usew3 = .false.
!usew4 = .false.
DEALLOCATE(q1)
DEALLOCATE(q2)



!old      tfs = tfs + (dt - dt1)*0.0484/(2.*ame)     ! ???



!      call manualPES(q0)



CALL cpu_time(time_fin)
time_cpu = time_fin-time_init
IF(myn == 0)THEN
  WRITE(6,'(a,1pg13.5)') ' CPU time in TSTEP',time_cpu
  WRITE(7,'(a,1pg13.5)') ' CPU time in TSTEP',time_cpu
END IF

IF (izforcecorr /= -1) THEN
  CALL checkzeroforce(rho,aloc)
  CALL checkzeroforce(rho,chpdft)
END IF

IF (MOD(it,jescmask) == 0 .OR. MOD(it,jescmaskorb) == 0) CALL  escmask(it)

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
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)
COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)
REAL(DP), INTENT(IN)                         :: dt



COMPLEX(DP) :: wfovlp
!JM
#if(twostsic)
INCLUDE 'radmatrix.inc'
INCLUDE 'twost.inc'
#endif
!JM

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
    CALL calc_fullsic(psiut,qnew)
  END IF
#endif
!JM
ELSE IF(ifsicp == 6) THEN
  STOP ' that kind of SIC not valid for dynamics'
END IF

RETURN
END SUBROUTINE dyn_mfield


!-----vstep--------------------------------------------------------

SUBROUTINE vstep(rho,psi,it,dt)

!     propagation of rare gas valence clouds by leapfrog method

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: psi(kdfull2,kstate)
INTEGER, INTENT(IN OUT)                  :: it
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


!-----cmult3D---------------------------------------------------------

SUBROUTINE cmult3d(q0,q1)

!     to multiply complex 3D array 'q0' by other array 'q1'
!     resulting in new 'q0' on output

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(OUT)                     :: q0(kdfull2)
COMPLEX(DP), INTENT(IN)                      :: q1(kdfull2)



!--------------------------------------------------------------------

DO i=1,kdfull2
  q0(i) = q1(i)*q0(i)
END DO

RETURN
END SUBROUTINE cmult3d


!-----booost--------------------------------------------------------

SUBROUTINE boost(q0)        ! boost with given 'centfx,y,z'

!     boosts electronic wavefunctions by 'centfx', 'centfy',
!     and 'centfz'.

USE params
USE kinetic
!USE coulsolv
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

SUBROUTINE info(psi,rho,aloc,akv,it)

!     information on energy observables: single particle energies,
!     kinetic energy, total energy, ionic energy, ...

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
!USE coulsolv
#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)
REAL(DP), INTENT(IN)                         :: rho(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: akv(kdfull2)
INTEGER, INTENT(IN)                      :: it

#if(vdw)
REAL(DP), PARAMETER :: alpha_ar=10.6
#endif
REAL(DP) ::  en(kstate)
COMPLEX(DP),ALLOCATABLE :: qtmp(:)
REAL(DP) ::  estar(2),estarp(2)


LOGICAL :: topenf
LOGICAL,PARAMETER :: ttest=.FALSE.

!------------------------------------------------------------------


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

eshell=0D0
esh1=0D0
enonlc = 0D0

DO nb=1,nstate
  CALL  calc_ekin(psi(1,nb),akv,ekin)
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
  WRITE(7,'(2(a,i2),4(a,f9.5))') 'level:',nrel2abs(nb), &
       '  ispin=',ispin(nrel2abs(nb)), &
       '  ekin='  &
      ,ekin,'  epot=',ehilf,'  esp=',amoy(nb) ,'  enonlo=', enonlo(nb)
  
END DO

IF(jstinf /= 0 .AND. mod(it,jstinf) == 0) then
  ALLOCATE(qtmp(kdfull2))
  DO nb=1,nstate
    qtmp = psi(:,nb)
    CALL hpsi(qtmp,aloc,akv,nb,1)  
  END DO
  CALL safeopen(77,it,jstinf,'pspenergies')
  WRITE(77,'(1f15.6,100f12.6)') tfs,(amoy(nb),nb=1,nstate)
  CALL flush(77)
  CALL safeopen(76,it,jstinf,'pspvariances')
  WRITE(76,'(1f15.6,100f12.6)') tfs,(spvariance(nb),nb=1,nstate)
  CALL flush(76)
!  WRITE(6,'(a,f15.6)') ' s.p. energies and variances written at tfs=',tfs
  DEALLOCATE(qtmp)
ENDIF

#if(parayes)
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
   CALL calc_estar(psi,iss,estar(iss))
END DO

WRITE(7,*) 'ekintot=',esh1
WRITE(7,*) 'estar(1)=',estar(1),'   estar(2)=',estar(2)



!     rearrangement and background Coulomb energy


ecback=0D0
ecrho=0D0
#if(vdw)
IF(nc > 0 .AND.ivdw == 1) esub = 0D0
#endif
DO ind=1,nxyz
  IF(nion2 /= 0) THEN
    ecback=ecback-rho(ind)*potion(ind)
    ecrho=ecrho+rho(ind)*(chpcoul(ind)-potion(ind))
#if(vdw)
    IF(nc > 0 .AND. ivdw == 1) esub = esub + rho(ind)*potvdw(ind)
#endif
  ELSE IF(nion2 == 0) THEN
    ecback=ecback-rhojel(ind)*chpcoul(ind)
    ecrho=ecrho+rho(ind)*chpcoul(ind)
  END IF
END DO
ecback=ecback*dvol/2.0
ecrho=ecrho*dvol/2.0

#if(vdw)
IF(nc > 0 .AND. ivdw == 1)THEN
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


!     ionic contributiosn to the energy


ecorr = energ_ions()

ekion=0D0        ! kinetic energy of Na cores
ekinion=0D0      ! kinetic energy of GSM cores
ekinel=0D0       ! kinetic energy of GSM shells
ekinkat=0D0      ! kinetic energy of kations
IF(ionmdtyp > 0) THEN
  DO ion=1,nion
    ek=0D0
    ek=ek+cpx(ion)*cpx(ion)
    ek=ek+cpy(ion)*cpy(ion)
    ek=ek+cpz(ion)*cpz(ion)
    xm=1836.0*amu(np(ion))*ame
    ek=ek/2.0/xm
    ekion=ekion+ek
  END DO
#if(raregas)
  DO ion=1,nc
    ek = 0D0
    ek=ek+pxc(ion)*pxc(ion)
    ek=ek+pyc(ion)*pyc(ion)
    ek=ek+pzc(ion)*pzc(ion)
    xm=1836.0*mion*ame
    ek=ek/2.0/xm
    ekinion=ekinion+ek
  END DO
  DO ion=1,NE
    ek = 0D0
    ek=ek+pxe(ion)*pxe(ion)
    ek=ek+pye(ion)*pye(ion)
    ek=ek+pze(ion)*pze(ion)
    xm=1836.0*me*ame
    ek=ek/2.0/xm
    ekinel=ekinel+ek
  END DO
  DO ion=1,nk
    ek = 0D0
    ek=ek+pxk(ion)*pxk(ion)
    ek=ek+pyk(ion)*pyk(ion)
    ek=ek+pzk(ion)*pzk(ion)
    xm=1836.0*mkat*ame
    ek=ek/2.0/xm
    ekinkat=ekinkat+ek
  END DO
#endif
END IF


!     final composition and print


energy=eshell+enrear+ecback+ecorr+enonlc/2.+ecrhoimage
#if(vdw)
IF(ivdw == 1)energy = energy + evdw
#endif
ecoul=ecback+ecrho+ecorr+ecrhoimage
etot = energy + ekion + ekinion + ekinel + ekinkat

!WRITE(*,*) ' before energy print: it,jenergy=',it,jenergy
IF (myn == 0 .AND. MOD(it,jenergy) == 0 ) THEN
!  WRITE(*,*) ' in energy print: it,jenergy=',it,jenergy
  CALL safeopen(163,it,jenergy,'penergies')
  WRITE(163,'(1f14.6,20e24.15)') tfs, &
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
     &                estar(1),estar(2)
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
#if(vdw)
  IF(ivdw == 1) WRITE(6,*) 'vdw energy       = ',evdw
#endif
  WRITE(6,*) 'binding energy   = ',energy
  WRITE(6,*) 'total energy     = ',etot
  IF (isurf /= 0) THEN
    WRITE(6,*) 'adsorb.energy = ',etot-enerinfty
  END IF
  
  
  IF(jinfo > 0 .AND. MOD(it,jinfo) == 0) THEN
    INQUIRE(17,OPENED=topenf)
    IF(.NOT.topenf) OPEN(17,POSITION='append',FILE='infosp.'//outnam)
    WRITE(17,'(a,i5,f8.3,1pg13.5)') 'time_step,time,energy=',it,tinfs,etot
    IF(.NOT.topenf)  CLOSE(17)
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


RETURN
END SUBROUTINE info

!-----calc_ekin------------------------------------------------------------

SUBROUTINE calc_ekin(psin,akv,ekinout)

!     calculates kinetic energy for single particle state with
!     complex wavefunction 'psin'.

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN OUT)                  :: psin(kdfull2)
REAL(DP), INTENT(IN)                         :: akv(kdfull2)
REAL(DP), INTENT(OUT)                        :: ekinout


COMPLEX(DP),DIMENSION(:),ALLOCATABLE :: psi2
!COMPLEX(DP) :: psi2(kdfull2)
COMPLEX(DP) :: wfovlp
!EQUIVALENCE(psi2(1),w1(1))

!------------------------------------------------------------------

!IF(usew1) STOP ' in CALC_EKIN: workspace W1 already active '
!IF(usew2) STOP ' in CALC_EKIN: workspace W2 already active '
!usew1 = .true.
!usew2 = .true.
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
sum = 0D0
DO i=1,nxyz
  sum = REAL(psi(i,nb))*REAL(psi2(i)) + imag(psi(i,nb))*imag(psi2(i))  &
      + sum
  sum0 = REAL(psi(i,nb))*REAL(psi(i,nb)) + imag(psi(i,nb))*imag(psi(i,nb))  &
      + sum0
END DO
ekinout = REAL(wfovlp(psi(1,nb),psi2))
#endif

!usew1 = .false.
!usew2 = .false.
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
USE kinetic
#if(twostsic)
USE twost
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN)                      :: psin(kdfull2)
REAL(DP), INTENT(IN)                         :: alocact(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: epotout
REAL(DP), INTENT(OUT)                        :: enonlocout
INTEGER, INTENT(IN)                      :: nb
INTEGER, INTENT(IN OUT)                  :: it


COMPLEX(DP),DIMENSION(:),ALLOCATABLE :: psi2
!COMPLEX(DP) :: psi2(kdfull2)
!EQUIVALENCE(psi2(1),w1(1))
COMPLEX(DP) :: wfovlp
!JM : dsic
!#if(twostsic)
COMPLEX(DP) :: cf
!#endif
!JM

!------------------------------------------------------------------

!IF(usew1) STOP ' in CALC_EPOT: workspace W1 already active '
!IF(usew2) STOP ' in CALC_EPOT: workspace W2 already active '
!usew1 = .true.
!usew2 = .true.
ALLOCATE(psi2(kdfull2))

epot=0D0


!     non local part of ps

!        if(ipsptyp.eq.1 .and. nrow(np(1)).gt.1) then
IF(ipsptyp == 1) THEN
  CALL nonlocalc(psin,psi2,0)
!??        if(it.LT.1) then
  sumnon = 0D0
  DO i=1,nxyz
    sumnon = sumnon + REAL(psin(i))*REAL(psi2(i)) +imag(psin(i))*imag(psi2(i))
  END DO
  enonlocout = sumnon*dvol
  DO i=1,nxyz
    psi2(i)=psi2(i)+alocact(i)*psin(i)
  END DO
!??        endif
ELSE
  DO i=1,nxyz
    psi2(i)=alocact(i)*psin(i)
  END DO
END IF
!JM : subtract SIC potential for state NB
#if(twostsic)
IF(ifsicp == 8) THEN
  is=ispin(nrel2abs(nb))
  DO na=1,ndim(is)
    nbe = nb - (is-1)*ndim(1)
    nae = na + (is-1)*ndim(1)
    cf = CONJG(vecs(nbe,na,is))
!            write(*,*)   ' nb,is,na,nbe,cf=',nb,is,na,nbe,cf
    DO i=1,nxyz
      psi2(i)=psi2(i)-qnew(i,nae)*cf
    END DO
  END DO
END IF
#endif
!JM
epot = REAL(wfovlp(psin,psi2))   ! develop real overlap

!      write(*,*) ' EPOT: nb,epot=',nb,epot

!usew1 = .false.
!usew2 = .false.
DEALLOCATE(psi2)

RETURN
END SUBROUTINE calc_epot

!-----calc_estar------------------------------------------------

SUBROUTINE calc_estar(psin,iss,excit)
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
REAL(DP), INTENT(OUT)                        :: excit

REAL(DP),DIMENSION(:),ALLOCATABLE            :: arho,arhop

!------------------------------------------------------------

ALLOCATE(arho(kdfull2))
ALLOCATE(arhop(kdfull2))

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
DO i=1,kdfull2
  CALL mpi_allreduce(arho(i),arhop(i),1,mpi_double_precision,  &
       mpi_sum,mpi_comm_world,ic)
ENDDO
CALL mpi_allreduce(ekintot,ekintotp,1,mpi_double_precision,  &
    mpi_sum,mpi_comm_world,ic)
CALL mpi_barrier (mpi_comm_world, mpi_ierror)
#endif

excit=0D0
DO i=1,kdfull2
   excit=excit+arhop(i)**1.666666666667D0
END DO
excit=excit*dvol*0.6D0*(6.0D0*pi**2)**0.666666666667D0


DEALLOCATE(arho)
DEALLOCATE(arhop)

excit=ekintotp-excit

! write(7,*) ' CALC_ESTAR: myn,iss,ekintot,excit',  &
!                          myn,iss,ekintotp,excit

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
USE kinetic
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
!  !USE coulsolv
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


#if(iangmo)

!-----instit----------------------------------------------------

SUBROUTINE instit (psi)

!   calculation of angular momentum

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)

COMPLEX(DP) :: akx(kdfull2),q2(kdfull2)
COMPLEX(DP) :: aky(kdfull2),akz(kdfull2)
!      complex eye
COMPLEX(DP) :: jtx(kdfull2),jalpha
COMPLEX(DP) :: jty(kdfull2),jtz(kdfull2)

!------------------------------------------------------------------

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

RETURN
END SUBROUTINE instit
#endif


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
USE kinetic
!USE coulsolv
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

SUBROUTINE init_dynprotocol(rho,aloc,akv,psi)

!     initializes file for protocol of dynamics

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: akv(kdfull2)
COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)




LOGICAL :: topenf

!---------------------------------------------------------------------


IF(irest <= 0) THEN                    !  write file headers
  
  IF(myn == 0 .AND. jdip /= 0) THEN
    OPEN(8,STATUS='unknown',FORM='formatted',FILE='pdip.'//outnam)
    WRITE(8,'(a)') ' & '
    WRITE(8,'(a)') 'x:time'
    WRITE(8,'(a)') 'y:dipole-momenta'
    WRITE(8,'(a)') 's: dist(4.0) scal(0.8) thickn(0.6) xmin(0.0)'
    WRITE(8,'(a)') 'n: time in fs  dipole-momenta: x   y   z'
    WRITE(8,'(a)') 'H:   X        Yl              Yd             Yq'
    CLOSE(8)
  END IF
  
  IF(jmp /= 0 .AND. myn == 0) THEN
    OPEN(803,STATUS='unknown',FILE='pMP.'//outnam)
    WRITE(803,'(a)') '# nr. of meas.points, nr. of state'
    WRITE(803,'(2i6)') nmps,nstate_all
    DO i=1,nmps
      WRITE(803,'(a,i6,a,3f12.1)') '# Point ',i , ' at : ',  &
          getxval(imps(i)),getyval(imps(i)),getzval(imps(i))
    END DO
!    CLOSE(803)
  END IF
  
  IF(jnorms /= 0) THEN
    OPEN(806,STATUS='unknown',FILE='pescOrb.'//outnam)
    WRITE(806,'(a,i6,a,3f12.1)') '# tfs, 1.0-norms of orbitals'
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
    
    IF(iangabso /= 0) THEN
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
#if(!adiadip)
! Velocities of Argon clouds (if available)
          OPEN(27,STATUS='unknown',FORM='formatted', FILE='pveldip.'//outnam)
          WRITE(27,'(a)') ' & '
          CLOSE(27)
#endif
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
    
    IF (isurf /= 0) THEN
      CALL adjustdip(rho)
      DO i=1,NE
        xeinit(i)=xe(i)
        yeinit(i)=ye(i)
        zeinit(i)=ze(i)
      END DO
      CALL info(psi,rho,aloc,akv,0)
      
      IF (myn == 0) THEN
        OPEN(308,POSITION='append',FILE='energies.res')
        dist = ABS(cz(1) - dmaximum(zc,nc))
        WRITE(308,'(1f20.10,2e25.14)') dist, energy-enerinfty,energy
        CLOSE(308)
      END IF
    END IF
    
    CALL getforces(rho,psi,0)
    
  END IF
END IF


RETURN
END SUBROUTINE init_dynprotocol



!-----print_densdiff--------------------------------------------------

SUBROUTINE print_densdiff(rho,it)

!     evaluation and print of density differences

USE params
USE kinetic
!USE coulsolv
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
  CALL addfields2(w1,rho,-1D0)
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
    CALL addfields2(w1,rho,-1D0)
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
    CALL addfields2(w1,rho,-1D0)
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
    CALL addfields2(w1,rho,-1D0)
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
    CALL addfields2(w1,rho,-1D0)
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
    CALL addfields2(w1,rho,-1D0)
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
    CALL addfields2(w1,rho,-1D0)
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
    CALL addfields2(w1,rho,-1D0)
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
    CALL addfields2(w1,rho,-1D0)
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
    CALL addfields2(w1,rho,-1D0)
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
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

!----------------------------------------------------------------


IF (e0 /= 0) CALL safeopen(38,0,1,'plaser')
CALL safeopen(163,it,jenergy,'penergies')
CALL safeopen(323,it,jcharges,'pcharges')

IF(nclust > 0 .AND. myn == 0)THEN
  
  CALL safeopen(8,it,jdip,'pdip')
  CALL safeopen(9,it,jquad,'pquad')
  CALL safeopen(17,it,jinfo,'infosp')
  CALL safeopen(47,0,iangabso,'pangabso')
  CALL safeopen(68,it,jang,'pangmo')
  CALL safeopen(78,it,jspdp,'pspdip')
  CALL safeopen(608,it,jgeomel,'pgeomel')
  CALL safeopen(806,it,jnorms,'pescOrb')
  
END IF
  
IF(nclust>0) CALL safeopen(803,it,jmp,'pMP')

RETURN
END SUBROUTINE open_protok_el

!-----analyze_ions-----------------------------------------------

SUBROUTINE analyze_ions(it)

!     open protocol files and print properties of ions

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

!----------------------------------------------------------------

!      tfs=it*dt1*0.0484

IF (jposcm > 0 .AND. MOD(it,jposcm) == 0) THEN
  CALL  safeopen(281,it,jposcm,'pposCM')
!  OPEN(281,POSITION='append',FILE='pposCM.'//outnam)
  CALL getcm(1,0,0,0)
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
  
  sumx=0D0
  sumy=0D0
  sumz=0D0
  DO ion=1,nion
    IF(MOD(it,jpos) == 0) WRITE(21,'(1f13.5,3e17.8,1pg13.5)')  &
        tfs,cx(ion),cy(ion),cz(ion),ecorr
    IF(MOD(it,jvel) == 0) WRITE(22,'(1f13.5,3e17.8,1pg13.5)')  &
        tfs,cpx(ion),cpy(ion),cpz(ion),ekion
    sumx = sumx + (cpx(ion)**2)/amu(np(ion))/1836.
    sumy = sumy + cpy(ion)**2/amu(np(ion))/1836.
    sumz = sumx + cpz(ion)**2/amu(np(ion))/1836.
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


!    ?? not clear what that is doing

IF (MOD(it,jpos) == 0 .AND. nclust == 0) THEN
!???            call info(psi,rho,aloc,akv,it)
  sum=0D0
  DO i=1,nc
    dist2 = getdistance2(i,i+nc)
    sum=sum+0.5D0*cspr*dist2
  END DO
  WRITE(773,'(1f12.5,1e17.7)') tfs,sum
END IF


RETURN
END SUBROUTINE analyze_ions

!-----analyze_surf------------------------------------------------

SUBROUTINE analyze_surf(it)

!     open protocol files and print properties of substrate

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

!----------------------------------------------------------------

!      tfs=it*dt1*0.0484

IF(nc+NE+nk > 0) THEN
  CALL safeopen(24,it,jpos,'pposcore')
  CALL safeopen(157,it,jpos,'ptempsurf')
  CALL safeopen(25,it,jpos,'pposdip')
  CALL safeopen(125,it,jpos,'pposkat')
  CALL safeopen(424,it,jpos,'penerSurf')
  CALL safeopen(26,it,jvel,'pvelcore')
  CALL safeopen(126,it,jvel,'pvelkat')
#if(!adiadip)
  CALL safeopen(27,it,jvel,'pveldip')
#endif
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
#if(!adiadip)
!                  if(nclust.gt.0) write(27,'(1f13.5,3e17.8)')
!     &              tfs,pxe(ion),pye(ion),pze(ion)
        IF (iprintonlyifmob == 0 .OR. imobe(ion) /= 0) THEN
          WRITE(27,'(1f13.5,3e17.8)') tfs,pxe(ion),pye(ion),pze(ion)
          nimobc=nimobc+1
        END IF
! else, there have no velocity
#endif
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
#if(!adiadip)
  CALL flush(27)
#endif
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


RETURN
END SUBROUTINE analyze_surf


!-----analyze_elect-----------------------------------------------

SUBROUTINE analyze_elect(psi,rho,aloc,akv,it)

!     Analysis and print of electronic properties during dynamics.

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)


COMPLEX(DP), INTENT(IN)                      :: psi(kdfull2,kstate)
REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: akv(kdfull2)
INTEGER, INTENT(OUT)                     :: it


LOGICAL,PARAMETER :: ttest=.FALSE.

COMPLEX(DP) :: orbitaloverlap
REAL(DP) :: rtmp(kstate,kstate)
COMPLEX(DP) :: cscal

!----------------------------------------------------------------

!      tfs=it*dt1*0.0484

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
  CALL localize(rho,psi)
END IF


!     escaped electrons analyzed orbital by orbital

IF(jesc > 0 .AND. MOD(it,jnorms) == 0) THEN
  DO i=1,nstate
    cscal=orbitaloverlap(psi(1,i),psi(1,i))
    rtmp(i,1)=REAL(cscal)**2+imag(cscal)**2
    rtmp(i,1)=1D0-SQRT(rtmp(i,1))
  END DO
  WRITE(806,'(12f12.8)') tfs,rtmp(1:nstate,1)
END IF


!     everything about single particles


IF(jstinf > 0 .AND. MOD(it,jstinf) == 0) THEN
#if(parayes) 
  IF(ttest) WRITE(*,*) ' ANALYZE before INFO: myn=',myn
#endif
#if(parano)
  IF(ttest) WRITE(6,'(a,i4)') ' INFO from ANALYZE. it=',it
#endif
  CALL info(psi,rho,aloc,akv,it)
#if(parayes) 
  IF(ttest) WRITE(*,*) ' ANALYZE after INFO: myn=',myn
#endif
END IF


#if(iangmo)
IF(iexcit == 1)  CALL instit (psi) !notes
#endif


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
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

#if(simpara)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)
INTEGER, INTENT(IN)                      :: it

REAL(4) :: tarray(2)

REAL(DP) :: trun

!----------------------------------------------------------------



IF(isave > 0 .AND. it /= 0 .AND. MOD(it,isave) == 0) THEN
  IF (irest /= 0 .AND. ABS(it-irest) <= 2) THEN
! do nothing, change later if needed
  ELSE
    CALL SAVE(psi,it,outnam)
  END IF
END IF

!     if walltime has almost expired, save all relevant data to
!     prepare for restart and stop:

!      iiii=etime(tarray,trun)
#if(elapse)
iiii=etime(tarray)
IF ((tarray(1)+tarray(2)) > trequest*timefrac) THEN
  CALL SAVE(psi,it,outnam)
  OPEN(660,STATUS='unknown',FILE='progstatus')
  WRITE(660,*) '0'
  CLOSE(660)
  STOP
END IF
#endif
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

