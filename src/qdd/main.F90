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

 
PROGRAM tdlda_m

!     Main program of code package for "electronic dynamics in clusters"
!     according to TDLDA-MD with various options for a SIC.

!     *******************************************
!     declarations
!     *******************************************

USE params
USE kinetic
USE util
USE coulsolv
USE twostr
USE twost
USE twost_util
USE orthmat
IMPLICIT NONE

!     psir = real wavefunctions in coord space
!     psi  = complex wavefunctions in coord space (for dynamics)
!     rho  = electronic density in coord space
!     aloc = mean-field-potential (to be initialized before dynamics)
!     chpcoul = Coulomb-potential of electrons
!     c$s/cp$s = auxiliary field to store ionic coords and momenta in between
!              the trial ionic propagation
!     rhojel = jellium density

#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

REAL(DP),ALLOCATABLE :: aloc(:),rho(:)
REAL(DP),ALLOCATABLE :: psir(:,:)
COMPLEX(DP),ALLOCATABLE :: psi(:,:),psiw(:,:)

INTEGER :: ion, it,nspup
INTEGER :: i
#if(raregas)
REAL(DP):: dt
#endif
!REAL(DP):: totalprob,totalovlp
REAL(DP):: time_absfin

!REAL(DP), EXTERNAL:: energ_ions   ! declared in ion_md
!REAL(DP), EXTERNAL:: enerkin_ions ! declared in ion_md

LOGICAL,PARAMETER :: imaginary_time= .true.     ! activate static afterburn

!     *******************************************
!
!     initializations
!
!     *******************************************


CALL cpu_time(time_absinit)

CALL init_parallele() 

CALL initnamelists          ! read all input parameters

CALL checkoptions()         !check coherence of preprocessor option

IF(icooltyp == 3) CALL init_simann() ! initialize simulated annealing

CALL init_baseparams()      !init grid size, number of states ...

#if(raregas)
CALL check_isrtyp           ! check short range part of substrate
#endif

CALL iperio                 ! initializing the 'periodic table'

CALL changeperio            ! pseudo-potential parameters

CALL iparams()              ! check dynamic  parameters

CALL init_grid()            ! init coulomb solver, kinetic energy, grid properties

CALL init_fields()          ! allocate basic arrays in module params


ALLOCATE(psir(kdfull2,kstate))
psir=0D0
ALLOCATE(aloc(2*kdfull2),rho(2*kdfull2))
aloc=0D0
rho=0D0

IF(myn == 0) CALL ocoption(7)   ! output compiled options
IF(myn == 0) CALL ocoption(8)   ! output compiled options

CALL init_output()              ! headers for basic output files

! Choice of ionic background :
WRITE(*,*) ' nion2=',nion2
SELECT CASE(nion2)
  CASE(0)
    CALL init_jellium()         ! initialize jellium background
  CASE(1)
    WRITE(*,*) ' ions switch'
    CALL initions()             ! reading ionic positions and inits (for005ion.<name>)
  CASE(2)
    WRITE(*,*) ' external background potential '
    CALL pseudo_external()      ! read exernal ionic pseudopotential from a file
  CASE DEFAULT
    STOP  'nion2 has incorrect value. nion2 must be 0, 1 or 2' 
END SELECT

CALL initwf(psir)              ! init wf, jellium, static parameters

#if(raregas)
                               !      initialize substrate
IF (isurf == 1) THEN
  CALL initfunctions
  CALL initsurface
END IF
#endif

#if(parayes)
CALL init_boxpara()
WRITE(*,*) 'lengnod:',lengnod
#endif

IF(ihome == 1) CALL init_homfield()  ! optional homogeneous exernal E field

CALL timer(1)                        ! set timer

IF(nclust > 0 .AND. ifsicp > 7) THEN
  CALL init_fsicr()                  ! initialize parameters for full SIC
END IF


!       *******************************************
!
!       static
!
!       *******************************************

IF(nclust > 0 .AND. irest == 0 .AND. ismax > 0)  THEN
   CALL statit(psir,rho,aloc)
END IF

!     switch to Monte-Carlo optimization of ionic configurations

IF(icooltyp == 3) THEN
  CALL simann(psir,rho,aloc)
  STOP ' Monte Carlo completed'
END IF

DEALLOCATE(psir)
!outnam=outname

!          imaginary-time iteration (to improve statics)
IF(imaginary_time .AND. isitmax>0 .AND. ismax>0) THEN
  ifexpevol = 1
  CALL afterburn(psir,rho,aloc)
END IF


!       *******************************************
!
!                     dynamic
!
!       *******************************************

CALL  init_dynamic

CALL dyn_propag(psi,rho,aloc)

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


! epilog

IF (jcharges /= 0) CLOSE(323)

#if(raregas)
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
#endif

if (jattach>0) THEN
  DEALLOCATE(psi_target)
  DEALLOCATE(ispin_target)
  DEALLOCATE(occ_target)
  DEALLOCATE(spe_target)
  DEALLOCATE(match)
END IF
IF(nproj_states>0) DEALLOCATE(proj_states)


CLOSE(163)
CLOSE(68)
CLOSE(8)
CLOSE(9)
CLOSE(78)
IF (jnorms > 0) THEN
  CLOSE(806)
  CLOSE(808)
ENDIF

!  ********************  end of main program ****************************

DEALLOCATE(psi)
DEALLOCATE(psitophi)!MV


#if(parayes)
CALL mpi_finalize(mpi_ierror)
#endif

IF(myn==0) THEN
CALL cpu_time(time_absfin)
OPEN(123,ACCESS='append',STATUS='unknown',FILE='Time')
#if(netlib_fft)
WRITE(123,*)'NETLIB'
#endif
#if(fftw_cpu)
WRITE(123,*)'FFTW'
#endif
WRITE(123,*)'Box :',nx2,ny2,nz2
WRITE(123,*)'Walltime =',time_absfin-time_absinit
CLOSE(123)
WRITE(6,*)'Walltime =',time_absfin-time_absinit
ENDIF

CONTAINS

SUBROUTINE init_dynamic()
!------------------------------------------------------------
!USE params
!USE util, ONLY: safeopen,rhointxy,rhointyz,rhointxz
!USE twost, ONLY: init_fsic
!USE orthmat

IMPLICIT NONE
!COMPLEX(DP), ALLOCATABLE, INTENT(IN OUT) :: psi(:,:)
!COMPLEX(DP), ALLOCATABLE, INTENT(IN OUT) :: psiw(:,:)
!REAL(DP), INTENT(IN OUT)    :: aloc(2*kdfull2)
!REAL(DP), INTENT(IN OUT)    :: rho(2*kdfull2)

REAL(DP), EXTERNAL:: energ_ions   ! declared in ion_md
INTEGER,PARAMETER :: itest=0
INTEGER :: i,ion,nspup


!     optionally initialize wavefunction and work arrays
ALLOCATE(psi(kdfull2,kstate))
psi=CMPLX(0D0,0D0,DP)
IF(nabsorb > 0 .AND. jescmaskorb /=0) ALLOCATE(rhoabsoorb(kdfull2,kstate))
IF(ifexpevol == 1) ALLOCATE(psiw(kdfull2,kstate))

!     initialize protocol files
IF(nclust > 0 .AND. nabsorb > 0) CALL init_absbc(rho)
IF (nclust > 0 .AND. jmp > 0) CALL initmeasurepoints
IF(myn == 0 .OR. knode == 1) THEN
#if(raregas)
  CALL init_dynprotocol(rho,aloc,psi)
#else
  CALL init_dynprotocol()
#endif
  IF(nclust==0) CALL getforces(rho,psi,-1,0)  ! initial call, case of only MD
END IF

#if(raregas)
IF (surftemp > 0) CALL init_surftemp()
#endif

IF(nclust > 0) THEN

  IF(ifsicp >= 7) THEN
    CALL init_fsic()
!    CALL end_fsicr()                !??    check and correct
  END IF
  
  IF(nabsorb > 0) CALL init_absbc(rho)
  IF(jmp > 0) CALL initmeasurepoints

  
!   *** HOW TO USE THE RESTART ***
!   1. to start a dynamic run using saved real wave functions:
!       set irest=0 and istat=1 in the dynamic part of for005.outnam
!   2. to resume a dynamic run, set irest=1 in dynamic part
!      of for005.outnam
!   3. save data every n iterations by setting isave=n
!   ******************************
  

! initialize or recover wavefunctions
  IF(irest == 0) THEN
    CALL init_dynwf(psi)
    IF(nabsorb > 0) CALL  init_abs_accum()
  ELSE
!    IF (ievaluate /= 0) CALL evaluate(rho,aloc,psi)    ! ????
    IF (iscatterelectron /=0) CALL init_scattel(psi)
    outnam=outname                                     ! ???
    CALL restart2(psi,outnam,.false.)
    WRITE(7,'(a,i3)') 'restart irest=',irest
  END IF


  ! order states in contingent blocks of spin
  ALLOCATE (psitophi(nstate,nstate))
  psitophi=0D0
  DO i=1,nstate
    psitophi(i,i)=1D0
  END DO
  IF(irest==0)  then  !only for first pass. If irest <>0 then no temperature
    CALL ordo_per_spin(psi)!MV
    CALL calcpseudo()
    CALL dyn_mfield(rho,aloc,psi,0D0)  
    CALL info(psi,rho,aloc,0)  !to update the occupation numbers and amoy
    IF(rtatempinit > 1D-6) THEN
      CALL fermi_init(amoy,rtatempinit,occup,1)!occup changed in Fermi
      CALL fermi_init(amoy,rtatempinit,occup,2)!occup changed in Fermi
      CALL dyn_mfield(rho,aloc,psi,0D0)
    END IF
    CALL info(psi,rho,aloc,0)  ! to update the occupation numbers and amoy
  ELSE                         ! irest <>0 recover array 'ispin'
    nspup=0
    DO i=1,nstate
      IF(ispin(i)==1) nspup=nspup+1    !count nspup
    END DO
    DO i=1,nstate
      IF(i <= nspup) THEN
        ispin(i)=1
      ELSE
        ispin(i)=2
      END IF
    END DO
    eqstnspup=nspup
    eqstnspdw=nstate-eqstnspup
  END IF

!  OPEN(1003,status='replace',file='pPES.'//outnam)      !???
!  open(1004,status='replace',file='prhov.'//outnam)     !???

! optinally refresh (pseudo)potentials
  IF(irest > 0 .OR. istat > 0) THEN
    CALL calcpseudo()
    CALL calclocal(rho,aloc)                          !  ??
    IF(ifsicp > 0) CALL calc_sic(rho,aloc,psi)
    IF(ipsptyp == 1) THEN
      DO ion=1,nion
        IF (iswitch_interpol==1) THEN
          CALL calc_projFine(cx(ion),cy(ion),cz(ion),cx(ion),cy(ion),cz(ion),ion)
          CALL mergetabs
        ELSE
          CALL calc_proj(cx(ion),cy(ion),cz(ion),cx(ion),cy(ion),cz(ion),ion)
        END IF
      END DO
    END IF
  END IF
  CALL dyn_mfield(rho,aloc,psi,0D0,0)
  
! intialize protocols  
  IF(irest == 0) CALL info(psi,rho,aloc,0)
  IF(iangmo == 1 .AND. iexcit == 1)  CALL instit(psi)    ! notes
  IF(irhoint_time == 1) THEN
    CALL rhointxy(rho,0)
    CALL rhointxz(rho,0)
    CALL rhointyz(rho,0)
  END IF
  
  time = irest*dt1*0.0484D0/(2D0*ame)
  IF(myn == 0 .OR. knode == 1) THEN
    WRITE(7,'(a,i6,a,f12.6,a,f16.5)')  &
        'iter=',irest,' tfs=',time,' total energy=',etot
    WRITE(7,'(a,f8.4,a,f7.2,a,3f9.2)') 't=',time,' moments: monop.=',qe(1),  &
        ' dip.: ',qe(2),qe(3),qe(4)
    WRITE(6,'(a,i6,a,f12.6,a,f16.5)')  &
        'iter=',irest,' tfs=',time,' total energy=',etot
    WRITE(6,'(a,f8.4,a,f7.2,a,3f9.2)') 't=',time,' moments: monop.=',qe(1),  &
        ' dip.: ',qe(2),qe(3),qe(4)
  END IF

  IF(myn == 0 .OR. knode == 1) CALL open_protok_el(0)
 
ELSE                             ! case of ionic dynamics only
  ecorr = energ_ions()
END IF                           ! end if nclust
WRITE(*,*) '5.:cpx,y,z:',cpx(1:nion),cpy(1:nion),cpz(1:nion)
IF(nclust > 0 .AND. itest==0) CALL analyze_elect(psi,rho,aloc,0)

OPEN(660,STATUS='unknown',FILE='progstatus')
WRITE(660,*) 'electronic initialization done'
CLOSE(660)

! intialize ionic dynamics
IF(ionmdtyp == 1 .AND. irest == 0)THEN
! leap frog first step : propagation of momenta by half a time step
  CALL lffirststep(rho,psi)
  OPEN(660,STATUS='unknown',FILE='progstatus')
  WRITE(660,*) 'leap-frog initialized'
  CLOSE(660)
END IF
ekionold=0D0


IF(irest == 0) THEN
  totintegprob=0.D0
  reference_energy=etot
END IF


END SUBROUTINE init_dynamic



END PROGRAM tdlda_m



