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
USE util
#if(netlib_fft|fftw_cpu)
USE coulsolv
#endif
#if(twostsic)
USE twostr
USE twost
USE localize_rad
USE orthmat
#endif
IMPLICIT NONE

!     non symmetrical dynamic fields

!     psir = real wavefunctions in coord space
!     psi  = complex wavefunctions in coord space (for dynamics)
!     rho  = electronic density in coord space
!     aloc = mean-field-potential (to be initialized before dynamics)
!     chpcoul = Coulomb-potential of electrons
!     c$s/cp$s = auxiliary field to store coords and momenta in between
!              the trial protonic propagation
!     rhojel = jel.density

#if(parayes||simpara)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

REAL(DP),ALLOCATABLE :: aloc(:),rho(:)

REAL(DP),ALLOCATABLE :: psir(:,:)
COMPLEX(DP),ALLOCATABLE :: psi(:,:),psiw(:,:)

INTEGER :: ion, it
#if(raregas)
INTEGER :: i
#endif
REAL(DP):: dt
REAL(DP):: totalprob,totalovlp
REAL(DP):: time_absfin
LOGICAL :: imaginary_time= .true.


REAL(DP), EXTERNAL:: energ_ions   ! declared in ion_md
REAL(DP), EXTERNAL:: enerkin_ions ! declared in ion_md

!~ #if(raregas)
!~ LOGICAL :: tmf
!~ #endif

!     *******************************************

!     init

!     *******************************************


CALL cpu_time(time_absinit)

CALL init_parallele() 

#if(fftw_gpu)
CALL cuda_gpu_init()
#endif

CALL initnamelists          ! read all input parameters

CALL checkoptions()       !check coherence of preprocessor option
  
IF(icooltyp == 3) CALL init_simann() ! initialize and read parameters for simulated annealing
  
CALL init_baseparams()    !init grid size, number of states ...

CALL check_isrtyp         ! check short range interaction matrix

CALL iperio                     ! initializing the 'periodic table'

CALL changeperio   ! overwrites default periodic system if necessary

CALL iparams()               ! check dynamic  parameters

CALL init_grid()    ! init coulomb solver, kinetic energy, grid properties

CALL init_fields()  ! allocate params arrays

#if(lda_gpu)
CALL cuda_lda_init() 
#endif

#if(twostsic)
IF(numspin==2) CALL init_radmatrix()   ! initialize matrices of radial moments
#endif

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
    CALL init_jellium()       ! initialize jellium background
  CASE(1)
    WRITE(*,*) ' ions switch'
    CALL initions()              ! reading ionic positions and inits.
  CASE(2)
    WRITE(*,*) ' external background potential '
    CALL pseudo_external()      ! read pseudopotential from a file
  CASE DEFAULT
    STOP  'nion2 has incorrect value. nion2 must be 0, 1 or 2' 
END SELECT

CALL initwf(psir)              ! init wf, jellium, static parameters

!                                     initialize surface


#if(raregas)
IF (isurf == 1) THEN
  CALL initfunctions
  CALL initsurface
END IF
#endif

#if(parayes)
CALL init_boxpara()
WRITE(*,*) 'lengnod:',lengnod
#endif


IF(ihome == 1) CALL init_homfield()  ! optional homogeneous E field

!      call priopt()                   !  print options


!   set timer

CALL timer(1)



!       *******************************************
!
!       static
!
!       *******************************************
!SIC: Self Interaction Correction
!                                     initialize parameters for FSIC
IF(nclust > 0 .AND. ifsicp >= 7) THEN
  CALL init_fsicr()
END IF

!IF(nclust > 0 .AND. irest == 0 .AND. istat == 0 .AND. ismax > 0)  
IF(nclust > 0 .AND. irest == 0 .AND. ismax > 0)  THEN
  CALL statit(psir,rho,aloc)
END IF


!     using Monte-Carlo

IF(icooltyp == 3) THEN
  CALL simann(psir,rho,aloc)
  STOP ' Monte Carlo completed'
END IF

DEALLOCATE(psir)


!       ****************************************************
!          imaginary-time iteration (to improve statics)
!       ****************************************************

IF(isitmax>0 .AND. ismax>0) THEN
   WRITE(*,*) ' start afterburn. isitmax=',isitmax
   CALL flush(6)
   ALLOCATE(psi(kdfull2,kstate))
#if(twostsic)
   IF(ifsicp >= 7 .AND. nclust > 0 )  CALL init_fsic()
   IF(ifsicp == 8) CALL expdabvol_rotate_init
#endif
   CALL restart2(psi,outnam,.true.)     ! read static wf's
#if(twostsic)
   IF(ifsicp>=7) CALL init_vecs()
#endif
   CALL calcpseudo()
   CALL calclocal(rho,aloc)                          !  ??
   IF(ifsicp > 0) CALL calc_sic(rho,aloc,psi)    ! Not pure LDA : use some type of SIC (SIC-GAM, ADSIC, SIC-Slater, SIC-KLI ...)
   IF(ipsptyp == 1) THEN ! full Goedecker pseudopotentials require projectors
    DO ion=1,nion
      IF (iswitch_interpol==1) THEN
        CALL calc_projFine(cx(ion),cy(ion),cz(ion),cx(ion),cy(ion),cz(ion),ion)
        CALL calc_proj(cx(ion),cy(ion),cz(ion),cx(ion),cy(ion),cz(ion),ion)
      ELSE
        CALL calc_proj(cx(ion),cy(ion),cz(ion),cx(ion),cy(ion),cz(ion),ion)
      END IF
     END DO
     IF(iswitch_interpol==1) CALL mergetabs
   END IF
   CALL info(psi,rho,aloc,-1)
   
   IF(ifexpevol == 1) ALLOCATE(psiw(kdfull2,kstate))
   IF(ifcnevol == 1) THEN
    WRITE(6,*) 'allocate space for CN propagation'
    ALLOCATE(psiw(kdfull2,kstate))
   ENDIF
   
   ! Start imaginary time iterations
   DO it=1,isitmax
     WRITE(*,*) ' afterburn. iteration=',it
     IF(ifexpevol == 1) THEN
       CALL tstep_exp(psi,aloc,rho,it,psiw,imaginary_time)
     ELSE
       STOP 'imaginary-time step requires exponential evolution'
!       CALL tstep(psi,aloc,rho,it)
     END IF
     CALL cschmidt(psi)    ! schmidt normalization  of results
     IF(MOD(it,istinf)==0)   CALL info(psi,rho,aloc,it)
   END DO
   DEALLOCATE(psiw)
   CALL info(psi,rho,aloc,-1)
   CALL SAVE(psi,-1,outnam)
   DEALLOCATE(psi)
#if(twostsic)
   IF(nclust>0 .AND. ifsicp>=7) CALL end_fsic()
#endif
   IF(itmax <= 0) STOP ' terminate with afterburn '
END IF


!       *******************************************
!                     dynamic
!       *******************************************

ALLOCATE(psi(kdfull2,kstate))
psi=CMPLX(0D0,0D0,DP)

!     optionally initialize work arrays
IF(nabsorb > 0 .AND. jescmaskorb /=0) ALLOCATE(rhoabsoorb(kdfull2,kstate))
IF(ifexpevol == 1) ALLOCATE(psiw(kdfull2,kstate))
IF(ifcnevol == 1) THEN
  WRITE(6,*) 'allocate space for CN propagation'
  ALLOCATE(psiw(kdfull2,kstate))
ENDIF

!     initialize protocol files

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

  IF(ifsicp >= 7) THEN
    CALL init_fsic()
!    CALL end_fsicr()                !??    check and correct
  END IF
  
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
    IF (ievaluate /= 0) CALL evaluate(rho,aloc,psi)
!                   ??: 'evaluate' should come after 'restart2' ??
    IF (iscatterelectron /=0) CALL init_scattel(psi)
    CALL restart2(psi,outnam,.false.)
    IF (iaddcluster /= 0) CALL addcluster(psi,outnam)
    WRITE(7,'(a,i3)') 'restart irest=',irest
!          if (irest.gt.1) irest=irest-1
!      WRITE(*,*) 'ndims=',ndims,' , vecs after restart:'
!      DO n=1,ndims(1)
!        WRITE(*,'(5(2f10.5,2x))') vecs(1:ndims(1),n,1)
!      END DO
!      DO n=1,ndims(2)
!        WRITE(*,'(5(2f10.5,2x))') vecs(1:ndims(2),n,2)
!      END DO
  END IF
!                                           refresh potentials
  IF(irest > 0 .OR. istat > 0) THEN
    CALL calcpseudo()
    CALL calclocal(rho,aloc)                          !  ??
    IF(ifsicp > 0) CALL calc_sic(rho,aloc,psi)
    IF(ipsptyp == 1) THEN
      DO ion=1,nion
#if(paraworld)
#else
        IF (iswitch_interpol==1) CALL calc_projFine(cx(ion),cy(ion),cz(ion),cx(ion),cy(ion),cz(ion),ion)
#endif
        CALL calc_proj(cx(ion),cy(ion),cz(ion),cx(ion),cy(ion),cz(ion),ion)
      END DO
#if(paraworld)
#else
      IF (iswitch_interpol==1) CALL mergetabs
#endif
    END IF
!      WRITE(*,*) 'ndims=',ndims,' , vecs after mean field 1:'
!      DO n=1,ndims(1)
!        WRITE(*,'(5(2f10.5,2x))') vecs(1:ndims(1),n,1)
!      END DO
!      DO n=1,ndims(2)
!        WRITE(*,'(5(2f10.5,2x))') vecs(1:ndims(2),n,2)
!      END DO
  END IF
  
#if(simpara)
  CALL mpi_barrier (mpi_comm_world, mpi_ierror)
#endif
  
  
!       initialization of mean field
  
  CALL dyn_mfield(rho,aloc,psi,0D0,0)
  
  IF(irest == 0) CALL info(psi,rho,aloc,0)
  IF(iangmo == 1 .AND. iexcit == 1)  CALL instit(psi)    ! notes
  IF(ifrhoint_time == 1) THEN
    CALL rhointxy(rho,0)
    CALL rhointxz(rho,0)
    CALL rhointyz(rho,0)
  END IF
  
  
!     initialize accumulators for absorbed density
  
  
  
  time = irest*dt1*0.0484D0/(2D0*ame)
#if(paraworld)
  time = time+dt1*0.0484D0/(2D0*ame)
#endif
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


 IF(irest == 0) THEN
      totintegprob=0.D0
      reference_energy=etot
 END IF


!---           here starts true propagation  --------------

CALL flush(7)
CALL stimer(1)
WRITE(*,*) 'before loop: cpx,y,z:',cpx(1:nion),cpy(1:nion),cpz(1:nion)

!cpx=0D0;cpy=0D0;cpz=0D0

CALL stimer(1)

time=0
tfs=0
DO it=irest,itmax   ! time-loop

  ijel=it
#if(paraworld)
  tfs=tfs+dt1*0.0484D0
#else
  tfs=it*dt1*0.0484D0 !/(2.0*ame)
#endif
  CALL print_densdiff(rho,it)       ! right place here ???
  
!         if(nclust.gt.0) call savings(psi,it)
  
  IF(jattach/=0 .AND. it==irest) THEN
     CALL init_occ_target()
     WRITE(*,*) 'nstate_target, after init_occ_target:', nstate_target
     ALLOCATE(psi_target(kdfull2,nstate_target))
     CALL init_psitarget()
  END IF

  IF(it > irest)THEN
    
    
!          pure electronic dynamics
    
    IF(nclust > 0) THEN
      
!~ #if(raregas)
!~       tmf=.false.
!~       IF(tmf) CALL loc_mfield_dummy(rho,aloc)   ! ???
!~ #endif
      
      
!     propagation of the wfs
      WRITE(*,*) 'propagation of the wfs'
      IF(ifexpevol == 1) THEN
        CALL tstep_exp(psi,aloc,rho,it,psiw, .NOT. imaginary_time)
      ELSE IF(ifcnevol == 1) THEN
        WRITE(*,*) 'Crank call '
        CALL  CrankNicolson_exp(psi,aloc,rho,it,psiw)
      ELSE
        CALL tstep(psi,aloc,rho,it)
      END IF
      IF(nabsorb > 0) CALL  absbc(psi,rho) ! Existence of absorbing poins on boundary
      
      
!            protocol of densities
      
      IF(ifrhoint_time == 1) THEN
        CALL rhointxy(rho,it)
        CALL rhointxz(rho,it)
        CALL rhointyz(rho,it)
      END IF
#if(raregas)
    ELSE
      IF(isurf /= 0 .AND. NE > 0) CALL valence_step(rho,dt,.true.)
#endif
    END IF
    
    IF(ionmdtyp==1 .OR. (ionmdtyp==2 .AND. MOD(it,modionstep)==0)) THEN
      
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
            IF (iswitch_interpol==1) then
              CALL calc_projFine(cx(ion),cy(ion),cz(ion),cx(ion),cy(ion),cz(ion),ion)
            ELSE
              CALL calc_proj(cx(ion),cy(ion),cz(ion),cx(ion),cy(ion),cz(ion),ion)
            ENDIF
          END DO
        END IF
      END IF
    END IF
    
  END IF
  
!     ******** compute and write observables: ********
  
  CALL timer(2)
  CALL stimer(2)
  
  IF(it > irest) THEN
    IF(myn == 0) THEN
#if(paraworld)
!     tfs=tfs+dt1*0.0484D0
#else
      tfs=it*dt1*0.0484D0
#endif
      IF(nion2 > 0) CALL analyze_ions(it)
      IF(isurf > 0) CALL analyze_surf(it)
    END IF
    IF(nclust > 0) THEN
      CALL analyze_elect(psi,rho,aloc,it)
    ELSE IF(MOD(it,100) == 0)THEN
      ecorr = energ_ions()
      etot = ecorr + ekion
    END IF
    ! The calculation of an electron attachment on a water molecule should 
    ! proceed as follows:
    ! 1) Perform a static calculation for a water molecule alone, with a 
    !    'deocc' high enough, so that all bound states, occupied and empty, 
    !    are calculated.
    ! 2) Create the file 'occ_spe_target' containing the binding energy and 
    !    the total energy as a first line, and then the isopin, occupation 
    !    number and s.p. energie of each state.
    ! 3) Perform a dynamical calculation with itmax=1 to generate a save file 
    !    which should be moved to the parent directory and renamed 'wfs_target'.
    ! 4) Perform a static and dynamical calculation with 'iscatterelectron=1' 
    !    and a small deocc (e.g., 'deocc=0.1'), so that only occupied states 
    !    of the water molecule are considered. The incoming electron is added 
    !    just before the dynamics as the last (occupied) state, and thus 
    !    labeled by nstate.
    IF(jattach /= 0 .AND. it>irest .AND. MOD(it,jattach) == 0) then
       call attach_prob(totalprob,totalovlp,psi)
       totintegprob=totintegprob+dt1*0.0484D0*jattach*totalprob
       write(6,'(a,e12.5,1x,i8,3(1x,1pg13.5))') &
            'after ATTACHEMENT:',&
            tfs,nmatch,totalprob,totintegprob,totalovlp
       CALL safeopen(809,it,jattach,'pattach')
       write(809,'(e12.5,1x,i8,3(1x,1pg13.5))') & 
            tfs,nmatch,totalprob,totintegprob,totalovlp
       CALL FLUSH(809)
    END IF
  
    IF(nclust > 0) CALL savings(psi,it)
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

#if(fftw_cpu|fftw_gpu)
!CALL fft_end()
!CALL coulsolv_end()
#endif

#if(fftw_gpu)
!CALL cuda_end()
#endif
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

IF(myn==0) THEN
CALL cpu_time(time_absfin)
OPEN(123,ACCESS='append',STATUS='unknown',FILE='Time')
#if(netlib_fft)
WRITE(123,*)'NETLIB'
#endif
#if(fftw_cpu)
WRITE(123,*)'FFTW'
#endif
#if(fftw_gpu)
WRITE(123,*)'cuFFT'
#endif
WRITE(123,*)'Box :',nx2,ny2,nz2
WRITE(123,*)'Walltime =',time_absfin-time_absinit
CLOSE(123)
ENDIF

END PROGRAM tdlda_m


!~ !#include "define.h"

!~ #if(raregas)
!~ !     ************************************

!~ SUBROUTINE loc_mfield_dummy(rho,aloc)

!~ !     ************************************

!~ !     This routine is a dummy version -- usage still unclear.
!~ !     Most probably obsolete!!!

!~ !     The local part of the mean field
!~ !     plus an update of the pseudo-potentials.

!~ !     Input:
!~ !      rho    = electron density
!~ !      dt     = stepsize (in case of dynamics)
!~ !      tdyn   = switch to dynamic case
!~ !     Output:
!~ !      aloc   = local mean field

!~ USE params
!~ #if(netlib_fft|fftw_cpu)
!~ USE coulsolv
!~ #endif
!~ IMPLICIT REAL(DP) (A-H,O-Z)


!~ REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
!~ REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)



!~ REAL(DP),ALLOCATABLE :: rhotmp(:)
!~ COMPLEX(DP) :: psidummy(1)


!~ ALLOCATE(rhotmp(2*kdfull2))
!~ IF(idielec /= 1) THEN
!~   CALL calcpseudo()
!~   CALL calcrho(rho,psi)
!~   DO ii=1,2*kdfull2
!~     rhotmp(ii)=rho(ii)
!~   END DO
!~   CALL addimage(rho,1)
  
!~ ! Coulomb of the electronic density
!~ #if(gridfft)
!~   CALL falr(rho,chpcoul,kdfull2)
!~ #endif
!~ #if(findiff|numerov)
!~   CALL solv_fft(rho,chpcoul,dx,dy,dz)
!~ #endif
!~   CALL calclocal(rho,aloc)
!~ END IF

!~ IF(idielec == 1) THEN
  
!~   DO ii=1,2*kdfull2
!~     rho(ii)=rhotmp(ii)
!~   END DO
!~   DO ii=1,kdfull2
!~     rfieldtmp(ii)=chpcoul(ii)
!~   END DO
!~   CALL addimage(rho,0)
  
!~ #if(gridfft)
!~ !         if(ipseudo.eq.1)
!~   CALL falr(rho,chpcoul,kdfull2)
!~ #endif
!~ #if(findiff|numerov)
!~ !         if(ipseudo.eq.1)
!~   CALL solv_fft(rho,chpcoul,dx,dy,dz)
!~ #endif
  
  
!~   DO ii=1,kdfull2
!~     CALL conv1to3(ii)
!~     IF(iindtmp(1) > nint(xdielec/dx)+nx)THEN
!~       chpcoul(ii)=rfieldtmp(ii)
!~     END IF
!~   END DO
!~   DO ii=1,2*kdfull2
!~     rho(ii)=rhotmp(ii)
!~   END DO
  
!~ END IF

!~ CALL calclocal(rho,aloc)
!~ IF(ifsicp > 0) CALL calc_sic(rho,aloc,psi)

!~ DEALLOCATE(rhotmp)

!~ RETURN
!~ END SUBROUTINE loc_mfield_dummy
!~ #endif

SUBROUTINE mergetabs

USE params
IMPLICIT NONE

INTEGER :: i, ia, ion, ifinsp, ja, ka, nfin, nfinsp, nfinfinesp
INTEGER,ALLOCATABLE :: ialltabfine(:)

ALLOCATE(ialltabfine(kdfull2fine))

ialltabfine=0

DO ion=1,nion
  nfin=ifinfine(ion)
  DO i=1,nfin   
    ialltabfine(icountfine(i,ion))=1
  END DO
END DO


ifinsp=0
DO i=1,kdfull2fine
  IF(ialltabfine(i).ne.0) THEN
    ifinsp=ifinsp+1
    IF(ifinsp.gt.(ng*knl)) STOP 'increase KNL'
    icountfinesp(ifinsp)=i
  ENDIF
END DO
nfinfinesp=ifinsp
!write(6,*) 'nfinfinesp',nfinfinesp

ialltabfine=0
DO ion=1,nion
  nfin=ifin(ion)
  DO i=1,nfin   
    ialltabfine(icount(i,ion))=1
  END DO
END DO


ifinsp=0
DO i=1,kdfull2
  IF(ialltabfine(i).ne.0) THEN
    ifinsp=ifinsp+1
    IF(ifinsp.gt.(ng*knl)) STOP 'increase KNL'
    icountsp(ifinsp)=i
    ia=i/(nx2*ny2)
    ja=(i-ia*nx2*ny2)/ny2
    ka=i-ia*nx2*ny2-ja*ny2
    END IF
END DO
nfinsp=ifinsp

DEALLOCATE(ialltabfine)
RETURN
END SUBROUTINE mergetabs




!!!! ******************** CALCUL OF PROJECTORS ON FINE GRID ******************** !!!!

SUBROUTINE calc_projFine(cxa,cya,cza,cxg,cyg,czg,ion)

!     Computes the projectors for the non-local part of
!     the Goedecker PsP for ion 'ion' at positions
!     'cxa,'cya',cza'. The position 'cxg', 'cyg', 'czg'
!     define the reference point for the center of the
!     subgrid. These two vector are usually the same,
!     but differ when computing forces.

USE params
IMPLICIT NONE

REAL(DP),PARAMETER :: fac0_12=-0.387298334621D0 ! -0.5D0*SQRT(3D0/5D0)

REAL(DP),INTENT(IN) :: cxa
REAL(DP),INTENT(IN) :: cya
REAL(DP),INTENT(IN) :: cza
REAL(DP),INTENT(IN) :: cxg
REAL(DP),INTENT(IN) :: cyg
REAL(DP),INTENT(IN) :: czg
INTEGER, INTENT(IN) :: ion

REAL(DP) :: h0_12
!---------------------------------------------------------

write(6,*) 'nxfine,nx2fine,nxyfine'
write(6,*) nxfine,nx2fine,nxyfine

IF(h0_12g(np(ion)).gt.-1D10) THEN
  h0_12=h0_12g(np(ion))
ELSE
  h0_12 = fac0_12*h0_22g(np(ion))
ENDIF

IF(ABS(h2_11g(np(ion))) + ABS(h1_22g(np(ion)))+ ABS(h0_33g(np(ion))) > small) THEN
  CALL calpr4Fine(cxa,cya,cza,cxg,cyg,czg,ion)
ELSE IF(ABS(h1_11g(np(ion))) + ABS(h0_22g(np(ion))) + ABS(h0_12) > small) THEN
  CALL calpr3Fine(cxa,cya,cza,cxg,cyg,czg,ion)
ELSE IF(ABS(h0_11g(np(ion))) > small) THEN
  CALL calpr2Fine(cxa,cya,cza,cxg,cyg,czg,ion)
END IF

END SUBROUTINE calc_projFine

!!     ****************************************
SUBROUTINE calpr2Fine(cxact,cyact,czact,cxg,cyg,czg,ion)

!     ****************************************

USE params
IMPLICIT NONE

REAL(DP),INTENT(IN) :: cxact
REAL(DP),INTENT(IN) :: cyact
REAL(DP),INTENT(IN) :: czact
REAL(DP),INTENT(IN) :: cxg
REAL(DP),INTENT(IN) :: cyg
REAL(DP),INTENT(IN) :: czg
INTEGER, INTENT(IN) :: ion

INTEGER :: i, ii, il, in, ind, inn, i1, i2, i3, i1l, i2l, i3l, icrsx, icrsy, icrsz
REAL(DP) :: dvolfine, dxfine, dyfine, dzfine  ! should be global variables ? (in params.F90 ?)
REAL(DP) :: r0, r1, radion, rfac, rr, x, y, z, xion, yion, zion
REAL(DP) :: gamfac, proj, sum1, sum4

WRITE(*,*) ' calpr2Fine'
dxfine=dx/2
dyfine=dy/2
dzfine=dz/2
dvolfine=dxfine*dyfine*dzfine
write(6,*) 'dxfine',dxfine

r0=r0g(np(ion))
r1=r1g(np(ion))
radion=radiong(np(ion))
!fine grid

i1l = nint((cxg-radion)/dxfine)
i2l = nint((cyg-radion)/dyfine)
i3l = nint((czg-radion)/dzfine)
icrsx = 2*nint(radion/dxfine)
icrsy = 2*nint(radion/dyfine)
icrsz = 2*nint(radion/dzfine)

!  compute projectors on these auxiliary grid:
ind=0
DO i3=0,icrsz
  z = (i3l+i3)*dzfine
  zion = z-czact
  DO i2=0,icrsy
    y = (i2l+i2)*dyfine
    yion = y-cyact
    DO i1=0,icrsx
      x = (i1l+i1)*dxfine
      xion = x-cxact
      rr=SQRT(MAX(xion*xion+yion*yion+zion*zion,small))
      rr=rr
      IF(rr <= radion) THEN
        
        ind = ind + 1
        IF(ind > knl) STOP " CALCPR: subgrid exeeded. enhance KNL"
        p0_1fine(ind,ion) = 0D0
!        p1_1fine(ind,ion) = 0D0
!        p1_1xfine(ind,ion) = 0D0
!        p1_1yfine(ind,ion) = 0D0
!        p1_1zfine(ind,ion) = 0D0
        
!      counter of gridpoint in one-dimensional density array:
        
        ii = ((i3l+i3+nzfine)-1)*nxyfine+((i2l+i2+nyfine)-1)*nx2fine+(i1l+i1+nxfine)
        icountfine(ind,ion)=ii
!        write(6,*) i3l,i3,nzfine,nxyfine
!        write(6,*) i2l,i2,nyfine,nx2fine
!        write(6,*) i1l,i1,nxfine
!        write(6,*) ind,ion,ii

        
!        DO in=1,2
        DO in=1,2
          IF(in == 1) THEN
!     projectors p0_1:
            inn=in
            il=0
            rfac=r0
            gamfac=0.5D0*SQRT(pi)
!          ELSE IF(in == 2) THEN
!!     projectors p1_1:
!            il=1
!            rfac=r1
!            inn=1
!            gamfac=1.5*0.5*SQRT(pi)
          END IF
          proj=SQRT(2D0)*rr**(il+2*(inn-1))*  &
              EXP(-(rr*rr)/(2D0*rfac**2D0))/(rfac**(il+(4*inn-1)/2D0)*  &
              SQRT(gamfac))
          
!    compose the product of the radial and the angular components:
!    p_i -> p_i*Y_lm
          IF(in == 1) p0_1fine(ind,ion) = proj
!          IF(in == 2) THEN
!            p1_1fine(ind,ion)  = proj
!            p1_1xfine(ind,ion) = proj*xion/rr
!            p1_1yfine(ind,ion) = proj*yion/rr
!            p1_1zfine(ind,ion) = proj*zion/rr
!                  write(6,'(2i5,4(1pg12.4))')
!     &               ion, ind,rr,p1_1xfine(ind,ion),p1_1yfine(ind,ion),
!     &               p1_1zfine(ind,ion)
!          END IF
        END DO
      END IF
      
    END DO
  END DO
END DO

!     end of counter-array:

ifinfine(ion) = ind


!     re-normalize projectors on grid


sum1=0D0
sum4=0D0
DO i=1,ifinfine(ion)
  sum1=sum1 + p0_1fine(i,ion)*p0_1fine(i,ion)
!  sum4=sum4 + p1_1(i,ion)*p1_1(i,ion)
END DO
sum1=sum1*dvolfine/(4D0*pi)
write(6,*) 'in calcpr2 sum1',sum1,ifinfine(ion)
sum1 = 1D0/SQRT(sum1)
!sum4=sum4*dvol/(4D0*pi)
!sum4 = 1D0/SQRT(sum4)
DO i=1,ifinfine(ion)
  p0_1fine(i,ion)=p0_1fine(i,ion)*sum1
!  write(6,*)i,ion, p0_1(i,ion)
!  p1_1(i,ion)=p1_1(i,ion)*sum4
!  p1_1x(i,ion)=p1_1x(i,ion)*sum4
!  p1_1y(i,ion)=p1_1y(i,ion)*sum4
!  p1_1z(i,ion)=p1_1z(i,ion)*sum4
END DO

RETURN
END SUBROUTINE calpr2Fine

!     ********************
SUBROUTINE calpr3Fine(cxact,cyact,czact,cxg,cyg,czg,ion)

!     ****************************************

USE params
IMPLICIT NONE

REAL(DP),INTENT(IN) :: cxact
REAL(DP),INTENT(IN) :: cyact
REAL(DP),INTENT(IN) :: czact
REAL(DP),INTENT(IN) :: cxg
REAL(DP),INTENT(IN) :: cyg
REAL(DP),INTENT(IN) :: czg
INTEGER, INTENT(IN) :: ion

INTEGER :: i, ii, il, in, ind, inn, i1, i2, i3, i1l, i2l, i3l, icrsx, icrsy, icrsz
REAL(DP) :: dvolfine, dxfine, dyfine, dzfine  ! should be global variables ? (in params.F90 ?)
REAL(DP) :: r0, r1, radion, rfac, rr, x, y, z, xion, yion, zion
REAL(DP) :: gamfac, proj, xnorm

WRITE(*,*) ' in calpr3Fine'


r0=r0g(np(ion))
r1=r1g(np(ion))
radion=radiong(np(ion))

!  compute boundaries of auxiliary grid for each ion:
i1l = nint((cxg-radion)/dxfine)
i2l = nint((cyg-radion)/dyfine)
i3l = nint((czg-radion)/dzfine)
icrsx = 2*nint(radion/dxfine)
icrsy = 2*nint(radion/dyfine)
icrsz = 2*nint(radion/dzfine)

!  compute projectors on these auxiliary grid:
ind=0
DO i3=0,icrsz
  z = (i3l+i3)*dzfine
  zion = z-czact
  DO i2=0,icrsy
    y = (i2l+i2)*dyfine
    yion = y-cyact
    DO i1=0,icrsx
      x = (i1l+i1)*dxfine
      xion = x-cxact
      rr=SQRT(MAX(xion*xion+yion*yion+zion*zion,small))
      IF(rr <= radion) THEN
        
        ind = ind + 1
        p0_1fine(ind,ion) = 0D0
        p0_2fine(ind,ion) = 0D0
        p1_1fine(ind,ion) = 0D0
        p1_1xfine(ind,ion) = 0D0
        p1_1yfine(ind,ion) = 0D0
        p1_1zfine(ind,ion) = 0D0
        
!      counter of gridpoint in one-dimensional density array:
        
        ii = ((i3l+i3+nzfine)-1)*nxyfine+((i2l+i2+nyfine)-1)*nx2fine+(i1l+i1+nxfine)
        icountfine(ind,ion)=ii
        
        DO in=1,3
          IF(in == 1 .OR. in == 2) THEN
!     projectors p0_1,p0_2:
            inn=in
            il=0
            rfac=r0
            IF(in == 1) THEN
              gamfac=0.5D0*SQRT(pi)
            ELSE IF(in == 2) THEN
              gamfac=2.5D0*1.5D0*0.5D0*SQRT(pi)
            END IF
          ELSE IF(in == 3) THEN
!     projectors p1_1,p1_2:
            il=1
            rfac=r1
            inn=1
            gamfac=1.5D0*0.5D0*SQRT(pi)
          END IF
          proj=SQRT(2D0)*rr**(il+2*(inn-1))*  &
              EXP(-(rr*rr)/(2D0*rfac**2D0))/(rfac**(il+(4*inn-1)/2D0)*  &
              SQRT(gamfac))
          
!    compose the product of the radial and the angular components:
!    p_i -> p_i*Y_lm
          IF(in == 1) p0_1fine(ind,ion) = proj
          IF(in == 2) p0_2fine(ind,ion) = proj
          IF(in == 3) THEN
            p1_1fine(ind,ion)  = proj
            p1_1xfine(ind,ion) = proj*xion/rr
            p1_1yfine(ind,ion) = proj*yion/rr
            p1_1zfine(ind,ion) = proj*zion/rr
!                  write(6,'(2i5,4(1pg12.4))')
!     &               ion, ind,rr,p1_1xfine(ind,ion),p1_1yfine(ind,ion),
!     &               p1_1zfine(ind,ion)
          END IF
        END DO
      END IF
      
    END DO
  END DO
END DO

!     end of counter-array:

! normalize  

xnorm = 1D0/SQRT(dvolfine*SUM(p0_1(:,ion)**2)/(4D0*pi))
p0_1(:,ion) = xnorm*p0_1(:,ion)
!WRITE(6,'(a,1pg13.5)') ' norm of 1s projector:',xnorm**(-2)

xnorm = 1D0/SQRT(dvolfine*SUM(p0_2(:,ion)**2)/(4D0*pi))
p0_2(:,ion) = xnorm*p0_2(:,ion)
!WRITE(6,'(a,1pg13.5)') ' norm of 2s projector:',xnorm**(-2)

xnorm = 1D0/SQRT(dvolfine*SUM(p1_1(:,ion)**2)/(4D0*pi))
p1_1(:,ion) = xnorm*p1_1(:,ion)
p1_1x(:,ion) = xnorm*p1_1x(:,ion)
p1_1y(:,ion) = xnorm*p1_1y(:,ion)
p1_1z(:,ion) = xnorm*p1_1z(:,ion)
!WRITE(6,'(a,1pg13.5)') ' norm of 1p projector:',xnorm**(-2)


ifinfine(ion) = ind
!         write(*,*) ' ion nr:',ion,' has ifinfine=',ind


END SUBROUTINE calpr3Fine

!     ********************
SUBROUTINE calpr4Fine(cxact,cyact,czact,cxg,cyg,czg,ion)

!     ****************************************

USE params
IMPLICIT NONE

REAL(DP),INTENT(IN) :: cxact
REAL(DP),INTENT(IN) :: cyact
REAL(DP),INTENT(IN) :: czact
REAL(DP),INTENT(IN) :: cxg
REAL(DP),INTENT(IN) :: cyg
REAL(DP),INTENT(IN) :: czg
INTEGER, INTENT(IN) :: ion

INTEGER :: i, ii, il, in, ind, inn, i1, i2, i3, i1l, i2l, i3l, icrsx, icrsy, icrsz
REAL(DP) :: dvolfine, dxfine, dyfine, dzfine  ! should be global variables ? (in params.F90 ?)
REAL(DP) :: r0, r1, r2, radion, rfac, rr, x, y, z, xion, yion, zion
REAL(DP) :: gamfac, proj, xnorm

dxfine=dx/2
dyfine=dy/2
dzfine=dz/2
r0=r0g(np(ion))
r1=r1g(np(ion))
r2=r2g(np(ion))
radion=3D0

!  compute boundaries of auxiliary grid for each ion:
i1l = nint((cxg-radion)/dxfine)
i2l = nint((cyg-radion)/dyfine)
i3l = nint((czg-radion)/dzfine)
icrsx = 2*nint(radion/dxfine)
icrsy = 2*nint(radion/dyfine)
icrsz = 2*nint(radion/dzfine)

!  compute projectors on these auxiliary grid:
ind=0
DO i3=0,icrsz
  z = (i3l+i3)*dzfine
  zion = z-czact
  DO i2=0,icrsy
    y = (i2l+i2)*dyfine
    yion = y-cyact
    DO i1=0,icrsx
      x = (i1l+i1)*dxfine
      xion = x-cxact
      rr=SQRT(MAX(xion*xion+yion*yion+zion*zion,small))
      IF(rr <= radion) THEN
        
        ind = ind + 1
        p0_1fine(ind,ion)   = 0D0
        p0_2fine(ind,ion)   = 0D0
        p0_3fine(ind,ion)   = 0D0
        p1_1fine(ind,ion)   = 0D0
        p1_1xfine(ind,ion)  = 0D0
        p1_1yfine(ind,ion)  = 0D0
        p1_1zfine(ind,ion)  = 0D0
        p1_2fine(ind,ion)   = 0D0
        p1_2xfine(ind,ion)  = 0D0
        p1_2yfine(ind,ion)  = 0D0
        p1_2zfine(ind,ion)  = 0D0
        p2_1fine(ind,ion)   = 0D0
        p2_xyfine(ind,ion)  = 0D0
        p2_xzfine(ind,ion)  = 0D0
        p2_yzfine(ind,ion)  = 0D0
        p2_xy2fine(ind,ion) = 0D0
        p2_z2fine(ind,ion)  = 0D0
        
!      counter of gridpoint in one-dimensional density array:
        
        ii = ((i3l+i3+nzfine)-1)*nxyfine+((i2l+i2+nyfine)-1)*nx2fine+(i1l+i1+nxfine)
        icountfine(ind,ion)=ii
        
        DO in=1,6
          IF(in == 1 .OR. in == 2 .OR. in == 3) THEN
!     projectors p0_1,p0_2,p0_3:
            inn=in
            il=0
            rfac=r0
            IF(in == 1) THEN
              gamfac=0.5D0*SQRT(pi)
            ELSE IF(in == 2) THEN
              gamfac=2.5D0*1.5D0*0.5D0*SQRT(pi)
            ELSE IF(in == 3) THEN
              gamfac=4.5D0*3.5D0*2.5D0*1.5D0*0.5D0*SQRT(pi)
            END IF
          ELSE IF(in == 4 .OR. in == 5) THEN
!     projectors p1_1,p1_2:
            il=1
            rfac=r1
            IF(in == 4) THEN
              inn=1
              gamfac=1.5D0*0.5D0*SQRT(pi)
            ELSE IF(in == 5) THEN
              inn=2
              gamfac=3.5D0*2.5D0*1.5D0*0.5D0*SQRT(pi)
            END IF
          ELSE IF(in == 6) THEN
!     projector p2_1:
            il=2
            inn=1
            rfac=r2
            gamfac=2.5D0*1.5D0*0.5D0*SQRT(pi)
          END IF
          proj=SQRT(2D0)*rr**(il+2*(inn-1))*  &
              EXP(-(rr*rr)/(2D0*rfac**2D0))/(rfac**(il+(4*inn-1)/2D0)*  &
              SQRT(gamfac))
          
!    compose the product of the radial and the angular components:
!    p_i -> p_i*Y_lm
          IF(in == 1) p0_1fine(ind,ion) = proj
          IF(in == 2) p0_2fine(ind,ion) = proj
          IF(in == 3) p0_3fine(ind,ion) = proj
          IF(in == 4) THEN
            p1_1fine(ind,ion)  = proj
            p1_1xfine(ind,ion) = proj*xion/rr
            p1_1yfine(ind,ion) = proj*yion/rr
            p1_1zfine(ind,ion) = proj*zion/rr
          END IF
          IF(in == 5) THEN
            p1_2fine(ind,ion)  = proj
            p1_2xfine(ind,ion) = proj*xion/rr
            p1_2yfine(ind,ion) = proj*yion/rr
            p1_2zfine(ind,ion) = proj*zion/rr
          END IF
          IF(in == 6) THEN
            p2_1fine(ind,ion)  = proj
            p2_xyfine(ind,ion) = proj*xion*yion/(rr*rr)
            p2_xzfine(ind,ion) = proj*xion*zion/(rr*rr)
            p2_yzfine(ind,ion) = proj*yion*zion/(rr*rr)
            p2_xy2fine(ind,ion)= proj*SQRT(3D0)*(xion*xion - yion*yion)/(rr*rr)
            p2_z2fine(ind,ion) = proj*(2D0*zion*zion  &
                - xion*xion-yion*yion)/(rr*rr)
          END IF
        END DO
      END IF
      
    END DO
  END DO
END DO

xnorm = 1D0/SQRT(dvol*SUM(p0_1fine(:,ion)**2)/(4D0*pi))
p0_1fine(:,ion) = xnorm*p0_1fine(:,ion)

!     end of counter-array:

ifinfine(ion) = ind
RETURN
END SUBROUTINE calpr4Fine

!     ********************
