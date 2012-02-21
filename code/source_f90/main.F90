#include "define.h"
 
PROGRAM tdlda_m

!     mastercode "electronic dynamics in clusters"

!     first part : electronic static lda at fixed ions
!     second part: tdlda for electrons + classical MD for ions (optional)
!     including local (na, h, li) or nonlocal pseudo-potentials
!     for elements from h to al

!     excitations = electronic shift, fs laser, scissor excitation

!     *******************************************
!     declarations
!     *******************************************

USE params
USE kinetic
USE coulsolv
#if(fullsic)
USE localize_rad
#endif
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

!INTEGER :: getnearestgridpoint
!INTEGER :: conv3to1

REAL(DP),ALLOCATABLE :: aloc(:),rho(:)

REAL(DP),ALLOCATABLE :: psir(:,:)
COMPLEX(DP),ALLOCATABLE :: psi(:,:),psiw(:,:)

LOGICAL :: tmf

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

CALL iperio                     ! initializing the 'periodic table'

CALL changeperio   ! overwrites default periodic system if necessary

CALL iparams()               ! check dynamic  parameters

CALL init_grid()

CALL init_fields()

#if(fullsic)
IF(numspin==2) CALL init_radmatrix()
#endif

ALLOCATE(psir(kdfull2,kstate))
psir=0D0
ALLOCATE(aloc(2*kdfull2),rho(2*kdfull2))
aloc=0D0
rho=0D0



IF(myn == 0) CALL ocoption(7)   ! output compiled options
IF(myn == 0) CALL ocoption(8)   ! output compiled options
CALL init_output()              ! headers for basic output files

WRITE(*,*) ' nion2=',nion2
IF(nion2 == 1) THEN
  WRITE(*,*) ' ions switch'
  CALL initions()              ! reading ionic positions and inits.
ELSE IF(nion2 > 1) THEN
  WRITE(*,*) ' external background potential '
  CALL pseudo_external()
ELSE
  CALL init_jellium()
END IF

CALL initwf(psir)              ! init wf, jellium, static parameters

!                                     initialize surface
#if(raregas)
IF (isurf == 1) THEN
  CALL initfunctions
  CALL initsurface
END IF
#endif

!                                     initialize parameters for FSIC
IF(ifsicp >= 7) CALL init_fsic()

IF(ihome == 1) CALL init_homfield()  ! optional homogeneous E field

!      call priopt()                   !  print options


!   set timer

CALL timer(1)



!       *******************************************
!
!       static
!
!       *******************************************


!IF(nclust > 0 .AND. irest == 0 .AND. istat == 0 .AND. ismax > 0)  
IF(nclust > 0 .AND. irest == 0 .AND. ismax > 0)  THEN
    CALL statit(psir,rho,aloc)
END IF


!     using monte-carlo

IF(icooltyp == 3) THEN
  CALL simann(psir,rho,aloc)
  STOP ' Monte Carlo completed'
END IF

DEALLOCATE(psir)


!       *******************************************

!       dynamic

!       *******************************************

ALLOCATE(psi(kdfull2,kstate))
psi=CMPLX(0D0,0D0)

!     optionally initialize work arrays
IF(jescmaskorb /=0) ALLOCATE(rhoabsoorb(kdfull2,kstate))
IF(ifexpevol == 1) ALLOCATE(psiw(kdfull2,kstate))


!     initialise protocol files

IF(nclust > 0 .AND. nabsorb > 0) CALL init_absbc(rho)
IF (nclust > 0 .AND. jmp > 0) CALL initmeasurepoints
#if(simpara)
CALL init_dynprotocol(rho,aloc,psi)
#else
IF(myn == 0 .OR. knode == 1) CALL init_dynprotocol(rho,aloc,psi)
#endif

#if(raregas)
IF (surftemp > 0) CALL init_surftemp()
#endif

IF(nclust > 0)THEN
  
  IF(nabsorb > 0) CALL init_absbc(rho)
  IF(jmp > 0) CALL initmeasurepoints

!      if (nabsorb.gt.0 .and. ispherAbso.ne.0)
!     &    call spherAbso(psi,-1)
  
  
!   *** HOW TO USE THE RESTART ***
!   1. to start a dynamic run using saved real wave functions:
!       set irest=0 and istat=1 in the dynamic part of for005.outnam
!   2. to resume a dynamic run, set irest=1 in dynamic part
!      of for005.outnam
!   3. save data every n iterations by setting isave=n
!   ******************************
  
  IF(irest == 0) THEN
    CALL init_dynwf(psi)
    IF(nabsorb > 0) CALL  init_abs_accum()
  ELSE
    IF (ievaluate /= 0) CALL evaluate(rho,aloc,psi,iflag)
!                   ??: 'evaluate' should come after 'restart2' ??
    CALL restart2(psi,outnam,.false.)
    IF (iaddcluster /= 0) CALL addcluster(psi,outnam)
    WRITE(7,'(a,i3)') 'restart irest=',irest
!          if (irest.gt.1) irest=irest-1
  END IF
!                                           refresh potentials
  IF(irest > 0 .OR. istat > 0) THEN
    CALL calcpseudo()
    CALL calclocal(rho,aloc)                          !  ??
    IF(ifsicp > 0) CALL calc_sic(rho,aloc,psi)
    IF(ipsptyp == 1) THEN
      DO ion=1,nion
        CALL calc_proj(cx(ion),cy(ion),cz(ion),cx(ion),cy(ion),cz(ion),ion)
      END DO
    END IF
  END IF
  
#if(simpara)
  CALL mpi_barrier (mpi_comm_world, mpi_ierror)
#endif
  
  
!       initialization of mean field
  
  CALL dyn_mfield(rho,aloc,psi,0D0)
  
  IF(irest == 0) CALL info(psi,rho,aloc,0)
  IF(iangmo == 1 .AND. iexcit == 1)  CALL instit(psi)    !notes
  IF(ifrhoint_time == 1) THEN
    CALL rhointxy(rho,0)
    CALL rhointxz(rho,0)
    CALL rhointyz(rho,0)
  END IF
  
  
!     initialize accumulators for absorbed density
  
  
  
  time = irest*dt1*0.0484/(2.0*ame)
#if(!simpara)
  IF(myn == 0 .OR. knode == 1) THEN
#endif
    WRITE(7,'(a,i6,a,f12.6,a,f16.5)')  &
        'iter=',irest,' tfs=',time,' total energy=',etot
    WRITE(7,'(a,f8.4,a,f7.2,a,3f9.2)') 't=',time,' moments: monop.=',qe(1),  &
        ' dip.: ',qe(2),qe(3),qe(4)
    WRITE(6,'(a,i6,a,f12.6,a,f16.5)')  &
        'iter=',irest,' tfs=',time,' total energy=',etot
    WRITE(6,'(a,f8.4,a,f7.2,a,3f9.2)') 't=',time,' moments: monop.=',qe(1),  &
        ' dip.: ',qe(2),qe(3),qe(4)
#if(!simpara)
  END IF
#endif

  IF(myn == 0 .OR. knode == 1) CALL open_protok_el(0)

!  IF(myn == 0 .OR. knode == 1) THEN
!    IF(jdip /= 0)     CALL safeopen(8,0,1,'pdip')           ! cPW     ueberflue!ssig?
!    IF(jquad /= 0)    CALL safeopen(9,0,1,'pquad')          ! cPW            "
!    IF(jinfo /= 0)    CALL safeopen(17,0,1,'infosp')        ! cPW            "
!    IF(jspdp /= 0)    CALL safeopen(78,0,1,'pspdip')        ! cPW            "
!    IF(jang /= 0)     CALL safeopen(68,0,1,'pangmo')        ! cPW            "
!    IF(jenergy /= 0)  CALL safeopen(163,0,1,'penergies')    ! cPW            "
!    IF(jangabso /= 0) CALL safeopen(47,0,1,'pangabso')      ! cPW            "
!    IF(jnorms /= 0)   CALL safeopen(806,0,1,'pescOrb')      ! cPW            " !  
!    IF(jgeomel /= 0)  CALL safeopen(608,0,1,'pgeomel')      ! cPW            "
!    IF(jmp /= 0)      CALL safeopen(803,0,1,'pMP')          ! cPW            "
!  END IF
  
ELSE
  ecorr = energ_ions()
END IF                           ! end if nclust
  WRITE(*,*) '5.:cpx,y,z:',cpx(1:nion),cpy(1:nion),cpz(1:nion)

OPEN(660,STATUS='unknown',FILE='progstatus')
WRITE(660,*) 'electronic initialization done'
CLOSE(660)


IF(ionmdtyp == 1 .AND. irest == 0)THEN
! leap frog first step : propagation of momenta by half a time step
  CALL lffirststep(rho,psi)
  OPEN(660,STATUS='unknown',FILE='progstatus')
  WRITE(660,*) 'leap-frog initialized'
  CLOSE(660)
END IF

!GB
IF(nclust > 0) THEN
!      CALL info(psi,rho,aloc,irest)  ! ???
  IF(irest == 0)  CALL analyze_elect(psi,rho,aloc,0)
END IF
!GB

ekionold=0D0



!---           here starts true propagation  --------------

WRITE(*,*) 'before loop: cpx,y,z:',cpx(1:nion),cpy(1:nion),cpz(1:nion)
!cpx=0D0;cpy=0D0;cpz=0D0
DO it=irest,itmax   ! time-loop
  
  iterat = it      ! to communicate time step
  ijel=it
  tfs=it*dt1*0.0484  !/(2.0*ame)
  
  CALL print_densdiff(rho,it)       ! right place here ???
  
!         if(nclust.gt.0) call savings(psi,tarray,it)
  
!test           call calcrho(rho,psi)      !  ??????????
  
  
  
  IF(it > irest)THEN
    
    
!          pure electronic dynamics
    
    IF(nclust > 0) THEN
      
#if(raregas)
      tmf=.false.
      IF(tmf) CALL loc_mfield_dummy(rho,aloc)   ! ???
#endif
      
      
!     propagation of the wfs
      
      IF(ifexpevol == 1) THEN
        CALL tstep_exp(psi,aloc,rho,it,psiw)
      ELSE
        CALL tstep(psi,aloc,rho,it)
      END IF
!      WRITE(*,*) ' MAIN: nabsorb=',nabsorb
      IF(nabsorb > 0) CALL  absbc(psi,rho,it)
      
      
!            protocol of densities
      
      IF(ifrhoint_time == 1) THEN
        CALL rhointxy(rho,it)
        CALL rhointxz(rho,it)
        CALL rhointyz(rho,it)
      END IF
    ELSE
      IF(isurf /= 0 .AND. NE > 0) CALL valence_step(rho,dt,.true.)
    END IF
    IF(ionmdtyp >= 1) THEN
      
!            ionic propagation
      
      IF(ionmdtyp == 1) CALL itstep(rho,it,psi)
      IF(ionmdtyp == 2) CALL itstepv(rho,it,psi)
      ekion =  enerkin_ions()
      IF(icooltyp > 0) CALL  reset_ions()
      
!            new total mean field  ('calclocal' still necessary ?)
      
      IF(nclust > 0)THEN
        CALL calcpseudo()
        CALL calclocal(rho,aloc)                          !  ??
        IF(ifsicp > 0) CALL calc_sic(rho,aloc,psi)
        IF(ipsptyp == 1) THEN
          DO ion=1,nion
            CALL calc_proj(cx(ion),cy(ion),cz(ion),cx(ion),cy(ion),cz(ion),ion)
          END DO
        END IF
      END IF
    END IF
    
  END IF
  
!     ******** compute and write observables: ********
  
  CALL timer(2)
  
  IF(it > irest) THEN
    IF(myn == 0) THEN
      tfs=it*dt1*0.0484
      IF(nion2 > 0) CALL analyze_ions(it)
      IF(isurf > 0) CALL analyze_surf(it)
    END IF
    IF(nclust > 0) THEN
      CALL analyze_elect(psi,rho,aloc,it)
    ELSE IF(MOD(it,100) == 0)THEN
      ecorr = energ_ions()
      etot = ecorr + ekion
    END IF
    IF(nclust > 0) CALL savings(psi,tarray,it)
  END IF
  
  
#if(simpara)
  CALL mpi_barrier (mpi_comm_world, mpi_ierror)
  WRITE(7,*) ' After barrier. myn,it=',myn,it
#endif
  
  
END DO

!  ********************  end of dynamic loop ****************************



OPEN(660,STATUS='unknown',FILE='progstatus')
WRITE(660,*) 'dynamics finished'
CLOSE(660)

#if(simpara)
CALL mpi_barrier (mpi_comm_world, mpi_ierror)
#endif

CLOSE(805)

#if(raregas)
IF (iuselast == -3) THEN
  OPEN(89,STATUS='unknown',FILE='for005surf.cool')
  WRITE(89,'(a,i6,a,i6)') '# nc = ',nc,' nk= ',nk
  DO i=1,nc
    WRITE(89,'(6e20.9,2i6)') xc(i),yc(i),zc(i),xe(i),ye(i),ze(i),imobc(i),  &
        imobe(i)
  END DO
  DO i=1,nk
    WRITE(89,'(3e20.9,i6)') xk(i),yk(i),zk(i),imobk(i)
  END DO
  CLOSE(89)
END IF
#endif




IF (jcharges /= 0) CLOSE(323)


IF (isurf == 1) THEN
  
  CLOSE(121)
  CLOSE(122)
  CLOSE(123)
  CLOSE(124)
  CLOSE(126)
  CLOSE(24)
  CLOSE(128)
  CLOSE(130)
  CLOSE(157)
  
END IF

CLOSE(163)
CLOSE(68)
CLOSE(8)
CLOSE(9)
CLOSE(78)
CLOSE(806)

!  ********************  end of main program ****************************

DEALLOCATE(psi)

!                                       ! ends 'else' of 'if(ifscan)'
!#endif

#if(parayes)
CALL mpi_finalize(icode)
#endif
#if(simpara)
WRITE(7,*) ' before final barrier. myn=',myn
CALL mpi_barrier (mpi_comm_world, mpi_ierror)
WRITE(7,*) ' after final barrier. myn=',myn
CALL mpi_finalize(icode)
#endif

END PROGRAM tdlda_m


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
!USE kinetic
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
