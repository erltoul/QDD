#include "define.h"
 
!-----statit---------------------------------------------------------

SUBROUTINE statit(psir,rho,aloc)

!     master routine for static iteration

USE params
! USE kinetic
USE coulsolv, ONLY:falr
#if(twostsic)
USE twostr
#endif
#if(fullsic)
USE localize_rad
#endif
IMPLICIT REAL(DP) (A-H,O-Z)
#if(parayes)
INCLUDE 'mpif.h'
REAL(DP) :: is(mpi_status_size)
#endif

REAL(DP), INTENT(IN OUT)                     :: psir(kdfull2,kstate)
REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)
!REAL(DP), INTENT(IN OUT)                     :: akv(kdfull2)

!REAL(DP) ::  rhotmp(2*kdfull2)
LOGICAL,PARAMETER :: tcpu=.true.
LOGICAL,PARAMETER :: tspinprint=.true.
LOGICAL,PARAMETER :: tp_prints=.false.

!MB:
#if(fullsic)
REAL(DP) :: qaux(kdfull2,kstate)
#else
REAL(DP) :: qaux(1,1)       ! dummy array
#endif
!MB/

! test Coulomb
CALL calcrhor(rho,psir)
WRITE(*,*) 'for charge=',SUM(rho)*dvol
CALL falr(rho,chpcoul,nx2,ny2,nz2,kdfull2)
CALL prifld(chpcoul,'coulomb pot')

!     Number of pre-iterations for static solution with IFSICP=6.
!     This parameter is to be set here "by hand" such that it can
!     be communicated by 'all.inc'

WRITE(*,*) ' in STATIT'

itersicp6=20


    IF(istat == 1) THEN
!        CALL sstep(psir,akv,aloc,0)    ! also for GSlat, DSIC
        CALL resume(psir,outnam)
        CALL calcpseudo()                 ! initialize pseudo-potentials
        CALL static_mfield(rho,aloc,psir,qaux,0)
        CALL pricm(rho)
        CALL infor(psir,rho,0)
    END IF


!     initial computation of mean field, protocol prints, headers


!      call priCM_state(psir)
CALL calcpseudo()                 ! initialize pseudo-potentials
CALL static_mfield(rho,aloc,psir,qaux,0)
CALL pricm(rho)
IF(myn == 0) WRITE(6,'(a,3e17.7)') 'rhoCM: ',  &
    rvectmp(1),rvectmp(2),rvectmp(3)
CALL infor(psir,rho,0)
#if(twostsic)
IF(ifsicp >= 7) CALL infor_sic(psir)
#endif


IF(myn == 0)THEN
  CALL prifld(rho,'density    ')
  CALL prifld(aloc,'potential   ')
  IF(nion2 /= 0) CALL prifld(potion,'potential_io')
  WRITE(7,'(f8.4,a,4f12.4)') 0.0,' initial moments',(qe(j),j=1,4)
  WRITE(6,'(f8.4,a,4f12.4)') 0.0,' initial moments',(qe(j),j=1,4)
  WRITE(7,*)
  
  WRITE(6,'(a)')'+++ start of static iteration +++'
  WRITE(7,'(a)')'+++ start of static iteration +++'
  
  IF(dpolx*dpolx+dpoly*dpoly*dpolz*dpolz > 0D0) WRITE(7,'(a,3f8.4)')  &
      ' static dipole potential: dpolx,dpoly,dpolz=', dpolx,dpoly,dpolz
  WRITE(7,*)
  WRITE(6,*) 'ismax=',ismax
END IF



!     the static iteration starts here



ifsicpsav = ifsicp       ! save for later recycling
sumvar = 100D0           ! initial value for terminator

DO iter1=1,ismax

!  WRITE(7,*) ' myn,iter1=',myn,iter1
!  CALL prifld2(7,rho,'density    ')
!  CALL prifld2(7,aloc,'KS potent. ')
!  IF(nion2 /= 0) CALL prifld2(7,potion,'potential_io')


  
!?         iterat = -iter1
  IF(tcpu) CALL cpu_time(time_init)
  
  IF(ifsicpsav == 4) THEN      ! switch safe pre-iterations for KLI
    IF(iter1 < 40) THEN
      ifsicp = 3
    ELSE
      ifsicp = ifsicpsav
    END IF
  END IF
  
  IF(ifsicpsav == 6) THEN
#if(fullsic)
    IF(iter1 <= itersicp6 .AND. sumvar > 10D0*epsoro) THEN
      ifsicp = 3   !   or 0 ?
    ELSE
      ifsicp = 6
!             if(iter1.le.itersicp6+10) then
      IF(iter1 <= 0) THEN
        WRITE(*,*) ' before LOCWF'
        CALL calc_locwf(psir,qaux)
        CALL calc_fullsicr(psir,qaux)
        WRITE(*,*) ' after LOCWF'
      END IF
    END IF
#else
    STOP "IFSICP=6 requires LOCSIC or FULLSIC"
#endif
!           write(*,*) ' ITERSICP6,IFSICP=',itersicp6,ifsicp   ! <--  strange line
  END IF
  
  IF(ifsicp /= 6) THEN
    CALL sstep(psir,aloc,iter1)    ! also for GSlat, DSIC
#if(fullsic)
  ELSE
    CALL  sicstep_gen(psir,qaux,aloc,iter1)
#endif
  END IF
  
  CALL static_mfield(rho,aloc,psir,qaux,iter1)
  CALL pricm(rho)
  IF(myn == 0) WRITE(6,'(a,3e17.7)') 'rhoCM: ',  &
      rvectmp(1),rvectmp(2),rvectmp(3)
  
  IF(MOD(iter1,istinf) == 0) THEN
!test          call prifld(aloc,'KS-potential')
    CALL infor(psir,rho,iter1)
#if(twostsic)
    IF(ifsicp >= 7) CALL infor_sic(psir)
#endif
    
    IF(myn == 0)THEN
      WRITE(7,'(a,i5)') 'iter= ',iter1
      WRITE(7,'(a,f12.4,a,2(/5x,5f12.4))') 'binding energy=',binerg,  &
          ', moments: monop.,dip,quad=', qe(1),qe(2),qe(3),qe(4),  &
          qe(5),qe(6),qe(7),qe(8),qe(9),qe(10)

IF(numspin==2)  WRITE(7,'(a,3f10.4)') 'spindipole',se(1),se(2),se(3)
!            write(7,*) ' sumvar,epsoro=',sumvar,epsoro
!            write(6,*) ' sumvar,epsoro=',sumvar,epsoro
    END IF
    IF(sumvar2 < epsoro) GO TO 99
  END IF
  

CALL flush(6)

END DO                                      ! end iteration loop
99   CONTINUE


IF(myn == 0)THEN
  WRITE(7,*) ' static iteration terminated with ', iter1,' iterations'
  
  CALL prifld(rho,'density    ')
  CALL prifld(aloc,'potential   ')
  CALL prifld(chpcoul,'Coul-potent.')
END IF




!     compute and print localization

IF(iflocaliz == 1) THEN                                 
#if(!parayes)
  CALL localizer(rho,psir)          
#else
  STOP ' LOCALIZER (switch IFLOCALIZ) should not be invoked in parallele code'   ! cPW
#endif
END IF

#if(twostsic)

!     diagonalize Lagrange parameters

IF(ifsicp >= 7) THEN
  CALL diag_lagr(psir)
!        itut=iter1         !!! for TD switch MaxLoc/SymCond   ???
END IF
#endif


!     final protocol on file 'pstat.<name>''

CALL pri_pstat(psir,iter1,rho)

!     save real wavefunctions for further applications

IF(tspinprint) CLOSE(12)          ! ???



!       call printSurfPot(592)

IF(tp_prints .AND. (myn == 0 .OR. knode == 1)) THEN
  CALL printfield(491,aloc,'tp.aloc')
  CALL printfield(492,rho,'tp.density')
  CALL printfield(496,chpcoul,'tp.coulomb')
  CALL printfield(497,potion,'tp.potion')
END IF

IF (iplotorbitals /= 0) THEN
  OPEN(522,STATUS='unknown',FILE='pOrbitals.'//outnam)
  DO i=1,nstate
    WRITE(522,'(a,i3)') '# state nr: ',i
    WRITE(522,'(a,f12.5)') '# occupation: ',occup(i)
    WRITE(522,'(a,f12.5)') '# s.p. energy: ',amoy(i)
    CALL printfield(522,psir(1,i),'tp.psir')
    WRITE(522,*)  ! separate blocks for gnuplot
    WRITE(522,*)  !
  END DO
  CLOSE(522)
END IF

!       call calcrhor(rhor,psir)
CALL calcrhor(rho,psir)


!old       if (jPlotDensityDiff.le.1e7) then
IF (jplotdensitydiff /= 0) THEN
  OPEN(590,STATUS='unknown',FILE='densdiff')
  DO i=1,kdfull2*2
    WRITE(590,*) rho(i)
  END DO
  CLOSE(590)
END IF


CALL rsave(psir,iter1,outnam)
istat=1                        ! prepare for reading in time step
IF(itmax == 0 .AND. isave > 0) THEN
         write(*,*) ' CALL RSAVE  1. case'
  CALL infor(psir,rho,-1)
  
#if(parayes)
  CALL mpi_finalize(icode)
#endif
  
  STOP ' terminate with static iteration '
END IF

!     file for groundstate at the end of static iteration

!     optionally print KS potential along axes
!     use 'jforce' as switch

IF(jforce > 0 .AND. myn == 0) CALL prifld(aloc,'potential   ')



IF(icooltyp == 1 .AND. itmax > 0) THEN
  WRITE(*,*) ' CALL RSAVE  2. case'
  CALL rsave(psir,iter1,outnam)
END IF

!k check 1p-h transition matrix elements to unoccpied state:
!k attention: nstate must be bigger than number of electrons

#if(parano)
IF(iftransme==1 .AND. (nstate > nclust)) CALL transel(psir)
#endif



#if(raregas)
IF (isurf /= 0) THEN
  IF (iuselast == -1) THEN
    IF (myn == 0) THEN
      OPEN(308,STATUS='unknown',FILE='for005surf.init')
      WRITE(308,*) nc,nk
      DO i=1,nc
        WRITE(308,'(6e17.7,2i6)') xc(i),yc(i),zc(i),xe(i),ye(i),  &
            ze(i),imobc(i),imobe(i)
      END DO
      DO i=1,nk
        WRITE(308,'(3e17.7,i6)') xk(i),yk(i),zk(i),imobk(i)
      END DO
      CLOSE(308)
    END IF
    STOP
  ELSE IF (iuselast == -2) THEN
    WRITE(*,*) ' ADJUSTDIP from STATIC uselast'
    CALL adjustdip(rho)
    IF (myn == 0) THEN
      OPEN(308,STATUS='unknown',FILE='for005surf.init')
      WRITE(308,*) nc,nk
      DO i=1,nc
        WRITE(308,'(6e17.7,2i6)') xc(i),yc(i),zc(i),xe(i),ye(i),  &
            ze(i),imobc(i),imobe(i)
      END DO
      DO i=1,nk
        WRITE(308,'(3e17.7,i6)') xk(i),yk(i),zk(i),imobk(i)
      END DO
      CLOSE(308)
    END IF
    STOP
  END IF
END IF
#endif


CALL pricm(rho)
IF(myn == 0) WRITE(6,'(a,3e17.7)') 'rhoCM: ',  &
    rvectmp(1),rvectmp(2),rvectmp(3)

xcm=rvectmp(1)
ycm=rvectmp(2)
zcm=rvectmp(3)

IF(idenspl /= 0) CALL mtv_fld(rho,1)

IF(itmax <= 0) STOP ' terminate with static iteration '

RETURN
END SUBROUTINE statit
!-----static_mfield------------------------------------------------

SUBROUTINE static_mfield(rho,aloc,psir,psiaux,iter1)

!     The Coulomb part of the mean field.

!     Input:
!      rho    = electron density
!      psir   = real wavefunctions
!     Output:
!      aloc   = local mean-field potential

USE params
#if(fullsic)
USE localize_rad
#endif
! USE kinetic
#if(twostsic)
USE twostr
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: psir(kdfull2,kstate)
REAL(DP), INTENT(IN OUT)                     :: psiaux(kdfull2,kstate)
INTEGER, INTENT(IN OUT)                  :: iter1




!----------------------------------------------------------------

CALL calcrhor(rho,psir)
CALL coul_mfield(rho)
#if(raregas)
IF (NE > 0) THEN
!   WRITE(*,*) ' ADJUSTDIP from STATIC mfield'
  CALL adjustdip(rho)
  CALL calcpseudo()                 ! update pseudo-potentials   ??
END IF
#endif
CALL calclocal(rho,aloc)          ! LDA part of the potential

!      call infor(psir,rho,0)

IF(ifsicp > 0 .AND.ifsicp < 6) THEN
  CALL calc_sicr(rho,aloc,psir)
#if(fullsic)
ELSE IF(ifsicp == 6) THEN
  CALL calc_fullsicr(psir,psiaux)
#endif
#if(twostsic)
ELSE IF(ifsicp >= 7 .AND. iter1 > 0) THEN
  CALL static_sicfield(rho,aloc,psir,iter1)
#endif
END IF

RETURN
END SUBROUTINE static_mfield


!#if(parano)
!-----sstep----------------------------------------------------------

SUBROUTINE sstep(q0,aloc,iter)

!     Performs one static step for all wavefunctions and for given
!     mean fields.
!     The step involves: action of H->psi, some analysis, and damping.

!     Optionally the mean field Hamiltonian is diagonalized in the
!     space of occupied states.

!     Input is the old real wavefunction 'q0'
!              the kinetic energy on 'akv'
!              the local mean field on 'aloc'
!              the iteration number 'iter' (for switching analysis)
!                 (only analysis, no stepping for 'iter=0') 
!     Ouput is new new real wavefunction 'q0'

!     This part for the serial version.

USE params
USE kinetic
!USE localize_rad
#if(fullsic)
USE localize_rad
#endif
#if(twostsic)
USE twostr
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: q0(kdfull2,kstate)
!REAL(DP), INTENT(IN)                         :: akv(kdfull2)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)
INTEGER, INTENT(IN)                      :: iter


REAL(DP) :: occold(kstate),ocwork(kstate)
#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif
!#if(hamdiag && parano)
INTEGER :: ktridig  !=(kstate+kstate*kstate)/2
REAL(DP),ALLOCATABLE :: hmatr(:,:)
REAL(DP),ALLOCATABLE :: heigen(:)
REAL(DP),ALLOCATABLE :: vect(:,:)
REAL(DP),ALLOCATABLE :: psistate(:)
INTEGER,ALLOCATABLE :: npoi(:,:)
INTEGER :: ntridig(2),nstsp(2)
LOGICAL, PARAMETER :: tprham=.false.
!#endif

COMPLEX(DP) :: csum



!REAL(DP) :: ph(ksttot)             ! degeneracy of wavefunction, for
INTEGER :: nph
LOGICAL :: tocc,tcpu
#if(parano)
DATA tocc,tcpu/.false.,.true./
#endif
#if(parayes)
DATA tocc,tcpu/.false.,.true./     ! no reoccupation in parallel
#endif
!       workspaces

REAL(DP),DIMENSION(:),ALLOCATABLE :: q1,w4
REAL(DP),ALLOCATABLE :: qex(:,:)
#if(gridfft)
COMPLEX(DP),DIMENSION(:),ALLOCATABLE :: psipr
COMPLEX(DP),DIMENSION(:),ALLOCATABLE :: q2
#if(fftw_gpu)
COMPLEX(C_DOUBLE_COMPLEX),DIMENSION(:),ALLOCATABLE :: q3
LOGICAL,PARAMETER :: copyback=.false.,recopy=.false.,recopy2=.false.
#endif
#else
REAL(DP),DIMENSION(:),ALLOCATABLE :: q2
#endif

!-------------------------------------------------------------------------

!      write(*,*) ' SSTEP with IFSICP=',ifsicp

nph=3-numspin

!     set timer

IF(tcpu) CALL cpu_time(time_init)


!     exact exchange, to be completed before wavefunctions are modified

IF(ifsicp == 5) THEN
  ALLOCATE(qex(kdfull2,kstate))
  CALL exchgr(q0,qex)
END IF


ALLOCATE(q1(kdfull2))
ALLOCATE(q2(kdfull2))
#if(fftw_gpu)
ALLOCATE(q3(kdfull2))
#endif


#if(parano)
IF(ifhamdiag == 1) THEN
  ktridig=(kstate+kstate*kstate)/2
  ALLOCATE(hmatr(ktridig,2),heigen(kstate),vect(kstate,kstate))
  ALLOCATE(psistate(kstate),npoi(kstate,2))
  ntridig(1) = 0
  ntridig(2) = 0
  nstsp(1) = 0
  nstsp(2) = 0
END IF
#endif


!?      dvol=dx*dy*dz
DO nbe=1,nstate
  ishift = (ispin(nbe)-1)*nxyz        ! store spin=2 in upper block
  
  
!       action of the potential in coordinate space
!       plus non local part of ps
!       plus optionally exchange part
  
  
  IF(ipsptyp == 1) THEN
    CALL nonlocalr(q0(1,nbe),q1)
    enonlo(nbe)= rwfovlp(q0(1,nbe),q1)
    DO  i=1,nxyz
      q1(i)=q1(i)+q0(i,nbe)*(aloc(i+ishift)-amoy(nbe))
    END DO
  ELSE
    DO  i=1,nxyz
      q1(i)=q0(i,nbe)*(aloc(i+ishift)-amoy(nbe))
    END DO
  END IF
  
  IF(ifsicp == 5) THEN
    q1=q1+qex(:,nbe)
  END IF
  

!JM : subtract SIC potential for state NBE
#if(twostsic)
  IF(ifsicp == 8) CALL subtr_sicpot(q1,nbe)
#endif
!JM
  
  
!       optionally compute expectation value of potential energy
  
  IF(MOD(iter,istinf) == 0) epotsp(nbe) = rwfovlp(q0(1,nbe),q1) + amoy(nbe)
  
  
#if(gridfft)
ALLOCATE(psipr(kdfull2))
  
!        action of the kinetic energy in momentum space
#if(netlib_fft|fftw_cpu)
  CALL rftf(q0(1,nbe),psipr)
!    ALLOCATE(w4(kdfull2))
!    CALL rfftback(psipr,w4)
!    WRITE(*,*) 'test RFTF: norms=',SUM(w4),SUM(q0(:,nbe))
!    CALL prifld(q0(1,nbe),'wavefct. in ')
!    CALL prifld(w4,'wavefct. out')
!    DEALLOCATE(w4)
  
!        FFT of V*psi
  CALL rftf(q1,q2)
#endif
#if(fftw_gpu)
  CALL rftf(q0(1,nbe),psipr,ffta,gpu_ffta,copyback)

  CALL rftf(q1,q2,ffta2,gpu_ffta2,copyback)
#endif
  
!       compose to h|psi>
#if(netlib_fft|fftw_cpu)
  q2 = psipr*akv+q2
#endif
#if(fftw_gpu)
!  q2 = psipr*akv+q2
  CALL hpsi_cuda(gpu_ffta,gpu_ffta2,gpu_akvfft,kdfull2)
#endif
  
!       Optionally compute expectation value of kinetic energy
!       and variance of mean field Hamiltonian.
!       This is done in Fourier space to save FFT's.
!       Variance 'evarsp2' excludes non-diagonal elements within
!       occupied space.
  
  IF(MOD(iter,istinf) == 0 .AND. ifsicp /= 6) THEN
    
#if(parano)
    
    sum0 = 0D0
    sumk = 0D0
    sume = 0D0
    sum2 = 0D0
#if(netlib_fft|fftw_cpu)
    DO  i=1,nxyz
      vol   = REAL(psipr(i))*REAL(psipr(i)) +imag(psipr(i))*imag(psipr(i))
      sum0  = vol + sum0
      sumk  = vol*akv(i) + sumk
      sume =  REAL(q2(i))*REAL(psipr(i)) +imag(q2(i))*imag(psipr(i))  + sume
      sum2 =  REAL(q2(i))*REAL(q2(i)) +imag(q2(i))*imag(q2(i))  + sum2
    END DO
#endif
#if(fftw_gpu)
    CALL sum_calc(sum0,sumk,sume,sum2,gpu_ffta,gpu_ffta2,gpu_akvfft,nxyz)
#endif

    ekinsp(nbe) = sumk/sum0
    sume = sume/sum0
    sum2 = sum2/sum0
!          write(6,*) ' norm,spe=',sum0,sume
!          amoy(nbe)   = sume
    evarsp(nbe) = SQRT(MAX(sum2-sume**2,small))

    ALLOCATE(w4(kdfull2))

#if(netlib_fft|fftw_cpu)
    CALL rfftback(psipr,w4)
    CALL rfftback(q2,w4)
#endif
#if(fftw_gpu)
    CALL gpu_to_gpu(gpu_ffta,gpu_ffta_int,kdfull2) !save gpu_ffta (psipr on GPU) for later
    CALL rfftback(psipr,w4,ffta,gpu_ffta_int,recopy)
    CALL gpu_to_gpu(gpu_ffta2,gpu_ffta_int,kdfull2) !save gpu_ffta2 (q2 on GPU) for later
    CALL rfftback(q2,w4,ffta2,gpu_ffta_int,recopy)
#endif
    CALL project(w4,w4,ispin(nbe),q0)
    evarsp2(nbe) =  SQRT(rwfovlp(w4,w4))
    DEALLOCATE(w4)
#endif

#if(parayes)
    evarsp2(nbe) = evarsp(nbe)
#endif
  END IF
  
#if(parano)
  IF(ifhamdiag == 1) THEN
!       accumulate mean-field Hamiltonian within occupied states,
!       for later diagonalization

    IF(iter > 0) THEN  
#if(netlib_fft|fftw_cpu)
      CALL rfftback(q2,q1)
#endif
#if(fftw_gpu)
      CALL gpu_to_gpu(gpu_ffta2,gpu_ffta_int,kdfull2) !save gpu_ffta2 (q2 on GPU) for later
      CALL rfftback(q2,q1,ffta2,gpu_ffta_int,recopy)
#endif
      iactsp = ispin(nbe)
      nstsp(iactsp) = 1+nstsp(iactsp)
      npoi(nstsp(iactsp),iactsp) = nbe
      DO nbc=1,nbe
        IF(iactsp == ispin(nbc)) THEN
          ntridig(iactsp) = 1+ntridig(iactsp)
          hmatr(ntridig(iactsp),iactsp) = rwfovlp(q0(1,nbc),q1)
          IF(nbc == nbe) hmatr(ntridig(iactsp),iactsp) =  &
              hmatr(ntridig(iactsp),iactsp) + amoy(nbe)
          IF(tprham) WRITE(6,'(a,2i5,1pg13.5)') ' nbe,nbc,hmatr=',nbe,nbc,  &
              hmatr(ntridig(iactsp),iactsp)
        END IF
      END DO
    END IF

  END IF
#endif
  
!     perform the damped gradient step and orthogonalise the new basis
  
  IF(idyniter /= 0 .AND. iter > 100) e0dmp = MAX(ABS(amoy(nbe)),0.5D0)
#if(netlib_fft|fftw_cpu)    
!  IF(iter > 0) THEN  
    IF(e0dmp > small) THEN
!       DO i=1,nxyz
!         psipr(i)=psipr(i)-q2(i)*epswf / (akv(i) + e0dmp )
!       END DO
      psipr = psipr - q2*epswf/(akv+e0dmp)
    ELSE
      psipr = psipr - epswf*q2
    END IF
    CALL rfftback(psipr,q0(1,nbe))
#endif
#if(fftw_gpu)
    IF(e0dmp > small) THEN
!      psipr = psipr - q2*epswf/(akv+e0dmp)
      CALL d_grad1(gpu_ffta,gpu_ffta2,gpu_akvfft,epswf,e0dmp,kdfull2)
    ELSE
!      psipr = psipr - epswf*q2
      CALL d_grad2(gpu_ffta,gpu_ffta2,epswf,kdfull2)
    END IF
    CALL rfftback(psipr,q0(1,nbe),ffta,gpu_ffta,recopy)
#endif
!  END IF
  
DEALLOCATE(psipr)
  
#endif
! end of FFT switch
  
  
  
#if(findiff|numerov)
  
!      action of kinetic energy, exp.values, and gradient step
!      (this version of finite differences needs yet development)
  
  IF(e0dmp > small) STOP 'damped gradient not yet for finite differences'
  CALL rkin3d(q0(1,nbe),q2)
  sumk = 0D0
  sume = 0D0
  sum2 = 0D0
  DO i=1,nxyz
    wfstep = (q2(i)+q1(i))
    wf0    = q0(i,nbe)
    sume   = wf0*wfstep + sume
    sum2   = wfstep*wfstep + sum2
    sumk   = wf0*q2(i) + sumk
    q0(i,nbe)=q0(i,nbe)-epswf*wfstep
  END DO
  sume = sume*dvol
  sum2 = sum2*dvol
  sumk = sumk*dvol
  ekinsp(nbe) = sumk
  evarsp(nbe) = SQRT(MAX(sum2-sume**2,small))
  amoy(nbe) = ekinsp(nbe)+epotsp(nbe)
#endif
  
  
END DO                                            ! end loop over states

DEALLOCATE(q1)
DEALLOCATE(q2)
#if(fftw_gpu)
DEALLOCATE(q3)
#endif
IF(ifsicp == 5)  DEALLOCATE(qex)


#if(parano)
IF(ifhamdiag == 1) THEN
!     diagonalize mean-field Hamiltonian
!     and transform occupied states correspondingly



  IF(tprham) THEN
    WRITE(6,'(a,2i5)')  ' nstsp=',nstsp
    WRITE(6,'(a)') 'npoi:'
    DO iactsp=1,2
      WRITE(6,'(10i5)') (npoi(nbcs,iactsp),nbcs=1,nstsp(iactsp))
    END DO
  END IF

  DO iactsp=1,2
    CALL givens(hmatr(1,iactsp),heigen,vect, nstsp(iactsp),nstsp(iactsp),kstate)
    IF(tprham) WRITE(6,'(a/20(1pg13.5))') ' eigenvalues:',  &
        (heigen(nbe),nbe=1,nstsp(iactsp))
    DO ii=1,nxyz
      psistate = 0D0
      DO nbes=1,nstsp(iactsp)
        nbe = npoi(nbes,iactsp)
!      IF(tprham .AND. ii==1) WRITE(6,'(a,3i5,a)')  &
!          ' iactsp,nbes,nbe=',iactsp,nbes,nbe,'  vect:'
!      IF(tprham) WRITE(6,'(10(1pg12.4))') (vect(nbc,nbe),nbc=1,nstsp(iactsp))
        DO nbcs=1,nstsp(iactsp)
          nbc = npoi(nbcs,iactsp)
          psistate(nbe) = psistate(nbe) + q0(ii,nbc)*vect(nbcs,nbes)
        END DO
      END DO
      DO nbes=1,nstsp(iactsp)
        nbe = npoi(nbes,iactsp)
        q0(ii,nbe) = psistate(nbe)
      END DO
    END DO
  END DO
  DEALLOCATE(hmatr,heigen,vect)
  DEALLOCATE(psistate,npoi)

END IF
#endif


!     Schmidt ortho-normalisation


CALL schmidt(q0)


!     readjust occupations numbers to actual s.p. energies


IF(tocc .AND. nclust < nstate*nph .AND. iter > istinf) CALL reocc()



!     protocol of SSTEP on standard output files

IF(tcpu) THEN
  CALL cpu_time(time_fin)
  time_cpu = time_fin-time_init
!        write(6,'(a,1pg13.5)') ' CPU time in SSTEP',time_cpu
!        write(7,'(a,1pg13.5)') ' CPU time in SSTEP',time_cpu
END IF
WRITE(6,'(a,i5,6(f10.4))') 'iter,up/down,CPU=',iter,se(4),se(5),time_cpu
WRITE(7,'(a,i5,6(f10.4))') 'iter,up/down,CPU=',iter,se(4),se(5),time_cpu



RETURN
END SUBROUTINE sstep
!#endif

!-----infor--------------------------------------------------------

SUBROUTINE infor(psi,rho,i)

!     Computes observables (energies, radii, ...)
!     and prints to standard output.

USE params
#if(fullsic)
USE localize_rad
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: psi(kdfull2,kstate)
REAL(DP), INTENT(IN)                         :: rho(2*kdfull2)
INTEGER, INTENT(IN)                      :: i
REAL(DP), PARAMETER :: alpha_ar=10.6
REAL(DP),SAVE :: energyold=0D0

COMPLEX(DP),DIMENSION(:),ALLOCATABLE :: psipr
COMPLEX(DP),DIMENSION(:),ALLOCATABLE :: psi2
!COMPLEX(DP) :: psipr(kdfull2)
!COMPLEX(DP) :: psi2(kdfull2)

#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif
!EQUIVALENCE (psi2(1),w1(1))
!EQUIVALENCE (psipr(1),w3(1))
REAL(DP) :: en(kstate)

ALLOCATE(psipr(kdfull2))
ALLOCATE(psi2(kdfull2))


!   compute h*psi

eshell=0D0
esh1=0D0

sumvar = 0D0
sumvar2 = 0D0
espnb  = 0D0
enonlc = 0D0
DO nb=1,nstate
  
  ekin     = ekinsp(nb)
  epot     = epotsp(nb)
  ehilf    = epot
  epot     = epot+ekin
  espnb    = espnb + (epotsp(nb) + ekinsp(nb))*occup(nb)
  amoy(nb) = epot
  epot     = epot+ekin
  en(nb)   = epot*occup(nb)
  eshell   = eshell+en(nb)
  esh1     = esh1+ekin*occup(nb)
  sumvar   = sumvar + occup(nb)*evarsp(nb)**2
  sumvar2   = sumvar2 + occup(nb)*evarsp2(nb)**2
  enonlc   = enonlc + enonlo(nb)*occup(nb)
!        if(nc+nk+ne.gt.0) ecorr=energ_ions()
END DO

ecorr=energ_ions()

#if(parano)
DO nb=1,nstate
  ekin = ekinsp(nb)
  IF(numspin==2) THEN
    WRITE(6,'(a,i2,a,i3,3f9.5,2(1pg12.4))')  &
      'level:',nrel2abs(nb),'  spin,occup,ekin,esp,var=',  &
      3-2*ispin(nrel2abs(nb)),occup(nb),ekin,amoy(nb), evarsp(nb),evarsp2(nb)
  ELSE
    WRITE(6,'(a,i2,a,3f9.5,1pg12.4)')  &
      'level:',nrel2abs(nb),'  occup,ekin,esp,var=',  &
      occup(nb),ekin,amoy(nb),evarsp(nb)
  END IF
END DO
#endif

#if(parayes)

CALL mpi_comm_rank(mpi_comm_world,myn,icode)
CALL mpi_barrier (mpi_comm_world, mpi_ierror)
CALL mpi_allreduce(espnb,espnbp,1,mpi_double_precision,  &
    mpi_sum,mpi_comm_world,ic)
CALL mpi_allreduce(eshell,eshellp,1,mpi_double_precision,  &
    mpi_sum,mpi_comm_world,ic)
CALL mpi_allreduce(esh1,esh1p,1,mpi_double_precision,  &
    mpi_sum,mpi_comm_world,ic)
CALL mpi_allreduce(enonlc,enonlcp,1,mpi_double_precision,  &
    mpi_sum,mpi_comm_world,ic)
CALL mpi_allreduce(sumvar,sumvarp,1,mpi_double_precision,  &
    mpi_sum,mpi_comm_world,ic)
CALL mpi_allreduce(sumvar2,sumvar2p,1,mpi_double_precision,  &
    mpi_sum,mpi_comm_world,ic)
CALL mpi_barrier (mpi_comm_world, mpi_ierror)
espnb=espnbp
esh1=esh1p
eshell=eshellp
enonlc=enonlcp
sumvar = sumvarp
sumvar2 = sumvar2p
CALL prispe_parallele(6,-1)
#endif
sumvar = SQRT(sumvar/REAL(nclust))
sumvar2 = SQRT(sumvar2/REAL(nclust))
IF(myn == 0) WRITE(6,'(a,2(1pg12.4))') ' total variance  =',sumvar,sumvar2


ecback=0D0
ecrho=0D0
DO ind=1,nxyz
  IF(nion2 /= 0) THEN
    ecback=ecback-rho(ind)*potion(ind)
    ecrho=ecrho+rho(ind)*(chpcoul(ind)-potion(ind))
  ELSE IF(nion2 == 0) THEN
    ecback=ecback-rhojel(ind)*chpcoul(ind)
    ecrho=ecrho+rho(ind)*chpcoul(ind)
  END IF
END DO
ecback=ecback*dvol/2D0
ecrho=ecrho*dvol/2D0

#if(raregas)
IF(nc > 0 .AND. ivdw == 1)THEN
  elnum = nclust*1D0
  esub = 0D0
  IF(nion2 /= 0) THEN
    DO ind=1,nxyz
      esub = esub + rho(ind)*potvdw(ind)
    END DO
  END IF
  esub = esub*dvol
  evdw = esub
  DO ion=1,nc
    DO ico=1,3
      evdw = evdw - 0.5*e2*alpha_ar*frho(ion,ico)*frho(ion,ico)/elnum
    END DO
  END DO
  espnb = espnb - esub
  eshell = eshell - esub
END IF
#endif
esh1=esh1
eshell=eshell/2D0  !(=t+v/2)

CALL energ_dielec(rho)
energy = espnb/2.+esh1/2.+enrear+ecback+ecorr+enonlc/2. -ecrhoimage
IF(directenergy) &
     energ2 = esh1+enerpw+ecrho+ecback+ecorr+enonlc -ecrhoimage
IF(ivdw == 1) energy = energy + evdw
binerg = energy
#if(parayes)
IF(myn == 0) THEN
#endif
  WRITE(6,*)  'sp pot. energy  =',espnb-esh1
  WRITE(6,*)  'sp kin. energy  =',esh1
  WRITE(6,*)  'tot sp energy   =',espnb
  WRITE(6,*)  'e_coul:ion-ion  =',ecorr
  WRITE(6,*)  'e_coul:el-ion   =',2.*ecback
  WRITE(6,*)  'e_coul:el-el    =',ecrho-ecback
  WRITE(6,*)  'e_coul:total    =',ecback+ecrho+ecorr
  WRITE(6,*)  'rearge. energy  =',enrear
  WRITE(6,*)  'nonlocal energy =',enonlc
  IF(idielec == 1) WRITE(6,*) 'dielect.Coul.e. =',ecrhoimage
  IF(ivdw == 1) WRITE(6,*)  'vdw energy      =',evdw
  WRITE(6,*)  'binding energy  =',energy
  IF (isurf /= 0) THEN
    WRITE(6,*)  'adsorption energy = ', energy-enerinfty
  END IF
  IF(directenergy) THEN
    WRITE(6,*)  'binding energy2 =',energ2
    WRITE(6,*)  'potential energ =',enerpw
    ensav  = energ2
    energ2 = energy
    energy = ensav
  END IF
  WRITE(6,'(a,i5,a,f12.6)') 'iter= ',i,'  binding energy',binerg
  WRITE(6,'(a)') ' '
  
  IF(jinfo > 0 .AND. MOD(i,jinfo) == 0) THEN
    CALL cleanfile(17)
    OPEN(17,POSITION='append',FILE='infosp.'//outnam)
    WRITE(17,'(a,i5,4(1pg13.5))') 'iteration,energy,variances=',  &
        i,energy,(energy-energyold)/jinfo,sumvar,sumvar2
    CALL flush(17)
    energyold = energy
  END IF
  
#if(parayes)
END IF
#endif

DEALLOCATE(psipr)
DEALLOCATE(psi2)

IF(myn == 0) WRITE(6,'(a,f8.4,a,f7.2,a,3f9.3)')  &
    'In info: t=',tfs,' moments: monop.=',qe(1), ' dip.: ',qe(2),qe(3),qe(4)


IF (idebug == 1) THEN
  
  WRITE(6,*) 'Printing KS potential...'
  
  WRITE(6,*) 'Printing electron density...'
  CALL  printfieldx(2000+i,rho,0,0)
  CALL  printfieldz(2100+i,rho,0,0)
  
  WRITE(6,*) 'Printing ion potential...'
  CALL  printfieldx(3000+i,potion,0,0)
  CALL  printfieldz(3100+i,potion,0,0)
  
END IF


!     protocol of dipoles and polarizability

IF(myn == 0 .AND. dpolx*dpolx+dpoly*dpoly+dpolz*dpolz > 0D0) THEN
  WRITE(7,'(a,3f8.4)') ' external dipole fields: dpolx,dpoly,dpolz=',  &
      dpolx,dpoly,dpolz
  IF(dpolx > 0D0) WRITE(7,'(a,1pg13.5)')  &
      ' dipole polarizability in x=',-qe(1)*qe(2)/dpolx
  IF(dpoly > 0D0) WRITE(7,'(a,1pg13.5)')  &
      ' dipole polarizability in y=',-qe(1)*qe(3)/dpoly
  IF(dpolz > 0D0) WRITE(7,'(a,1pg13.5)')  &
      ' dipole polarizability in z=',-qe(1)*qe(4)/dpolz
  CALL flush(7)
END IF


!mb for creating the potential energy curve of Na on MgO
#if(raregas)
IF (i == -1 .AND. myn == 0) THEN ! static iteration ended, print out binding energies
  OPEN(221,POSITION='append',FILE='penerstat')
  WRITE(221,'(1f15.5,8e20.10)') zc(1),energy,ecorr,2*ecback, ecrho-ecback,  &
      amoy(1)
! ,amoy(2),amoy(3),amoy(4)
  CLOSE(221)
END IF
#endif
!mb/


RETURN
END SUBROUTINE infor

!-----pri_pstat----------------------------------------------------

SUBROUTINE pri_pstat(psi,i,rho)

!     print short protocol on file 'pstat.*'


USE params
! USE kinetic
USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)             :: psi(kdfull2,kstate)
INTEGER, INTENT(IN OUT)              :: i
REAL(DP), INTENT(IN OUT)             :: rho(kdfull2)

COMPLEX(DP),DIMENSION(:),ALLOCATABLE :: psipr
COMPLEX(DP),DIMENSION(:),ALLOCATABLE :: psi2
!COMPLEX(DP) :: psipr(kdfull2)
!COMPLEX(DP) :: psi2(kdfull2)
!EQUIVALENCE (psi2(1),w1(1))
!EQUIVALENCE (psipr(1),w3(1))
REAL(DP) :: en(kstate)
LOGICAL,PARAMETER :: mumpri=.false.

#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

omegam = omega_mieplasmon(rho)

ALLOCATE(psipr(kdfull2))
ALLOCATE(psi2(kdfull2))


!      eshell=0.0
!      esh1=0.0
!      sumvar = 0.0
!      espnb  = 0.0
!      enonlc = 0.0

#if(parayes)
CALL mpi_barrier (mpi_comm_world, mpi_ierror)
CALL  mpi_comm_rank(mpi_comm_world,myn,icode)
IF(myn == 0) THEN
#endif
  OPEN(42,POSITION='append',FILE='pstat.'//outnam)
  IF(ifsicp == 6) THEN
    WRITE(42,'(a,i3,a,i5)') 'final protocol of static for IFSICP=',ifsicp,  &
        ', pre-iterations with Slater=',itersicp6
  ELSE
    WRITE(42,'(a,i3)') 'final protocol of static for IFSICP=',ifsicp
  END IF
  IF(idielec /= 0) WRITE(42,'(a,i3,a,f10.5)')  &
      'iDielec=',idielec,'   ,xDielec=',xdielec
#if(parayes)
END IF
#endif

#if(coufou)
CALL mulmws(rho,p00,p10,p11r,p11i,p20,p21r,p21i,  &
    p22r,p22i,p30,p31r,p31i,p32r,p32i,p33r,p33i,  &
    p40,p41r,p41i,p42r,p42i,p43r,p43i,p44r,p44i, pr2,mumpri)
#endif


#if(parano)
DO nb=1,nstate
  IF(numspin==2) THEN
    WRITE(42,'(a,i2,a,i3,3f9.5,1pg12.4)')  &
      'level:',nrel2abs(nb),'  spin,occup,ekin,esp,variance =',  &
      3-2*ispin(nrel2abs(nb)),occup(nb),ekinsp(nb), amoy(nb),evarsp(nb)
  ELSE
    WRITE(42,'(a,i2,a,3f9.5,1pg12.4)')  &
      'level:',nrel2abs(nb),'  occup,ekin,esp,variance=',  &
      occup(nb),ekinsp(nb),amoy(nb),evarsp(nb)
  END IF
END DO
#endif
#if(parayes)
CALL prispe_parallele(42,-1)
#endif

IF(myn == 0) THEN
  IF(temp == 0D0) THEN
    IF(directenergy) THEN
      WRITE(42,'(a,2f13.7)') 'binding energy  =',energ2,energy
    ELSE
      WRITE(42,'(a,f13.7)') 'binding energy  =',energy
    END IF
  ELSE
    IF(directenergy) THEN
      energy = energ2
    END IF
    WRITE(42,'(a,4f13.7)')  'energies:E,TS,F,T =',  &
        energy,temp*entrop,energy-temp*entrop,temp
  END IF
  WRITE(42,'(a,1pg12.4)') 'total variance  =',sumvar
  WRITE(42,'(a,4f11.5)') 'sp pot, sp kin, rearr, nonlocal=',  &
      espnb-esh1,esh1,enrear,enonlc
  IF(idielec == 1) THEN
    WRITE(42,'(a,5f11.5)') 'e_coul: i-i , e-i , e-e , e-eimage, total=',  &
        ecorr,2D0*ecback,ecrho-ecback-ecrhoimage,2D0*ecrhoimage,  &
        ecback+ecrho+ecorr+ecrhoimage
  ELSE
    WRITE(42,'(a,4f11.5)') 'e_coul: i-i , e-i , e-e , total=',  &
        ecorr,2D0*ecback,ecrho-ecback,ecback+ecrho+ecorr
  END IF
  WRITE(42,'(a,f7.2)')    'mon.:',qe(1)
  WRITE(42,'(a,3f11.5)')  'dip.in  :',dpolx,dpoly,dpolz
  WRITE(42,'(a,3f11.5)')  'dip.out :',qe(2),qe(3),qe(4)
  WRITE(42,'(a)')         'quadrupol moments:'
  WRITE(42,'(a,3f11.4)')  'xx,yy,zz:',qe(5),qe(6),qe(7)
  WRITE(42,'(a,3f11.4)')  'xy,zx,zy:',qe(8),qe(9),qe(10)
  rms = SQRT(qe(5)+qe(6)+qe(7))
  WRITE(42,'(a,3f11.4)')  'q20,30,4:',p20,p30,p40
  WRITE(42,'(a,3f11.4)')  'renorm. :',  &
      p20/(qe(1)*rms**2),p30/(qe(1)*rms**3),p40/(qe(1)*rms**4)
  WRITE(42,'(a,3f11.4)')  'spindip.:',se(1),se(2),se(3)
  WRITE(42,'(a,4f11.4)') 'omegam,rhops,N_el,rhomix:',  &
      omegam,rhopss,apnum,rhomix
END IF

IF(ifspemoms == 1) CALL spmomsr(psi,42)

IF(myn == 0) THEN
  WRITE(42,'(1x)')
  CLOSE(42)
END IF

DEALLOCATE(psipr)
DEALLOCATE(psi2)

RETURN
END SUBROUTINE pri_pstat



!-----transel-------------------------------------------------transel--

SUBROUTINE transel(psir)

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)


REAL(DP), INTENT(IN)                         :: psir(kdfull2,kstate)
REAL(DP) ::  te(kstate,kstate,3)

!      nx=nx2/2
!      ny=ny2/2
!      nz=nz2/2

OPEN(34,STATUS='unknown',FILE='mte-xyz')

DO nb1=1,nclust
  DO nb2=nclust+1,nstate
    te(nb1,nb2,1) = 0D0
    te(nb1,nb2,2) = 0D0
    te(nb1,nb2,3) = 0D0
    
    ind=0
    DO iz=minz,maxz
      z1=(iz-nzsh)*dz
      DO iy=miny,maxy
        y1=(iy-nysh)*dy
        DO ix=minx,maxx
          ind=ind+1
          IF((ix /= nx2).AND.(iy /= ny2).AND.(iz /= nz2))THEN
            x1=(ix-nxsh)*dx
            
            te(nb1,nb2,1)=te(nb1,nb2,1)+psir(ind,nb1)*x1*psir(ind,nb2)
            te(nb1,nb2,2)=te(nb1,nb2,2)+psir(ind,nb1)*y1*psir(ind,nb2)
            te(nb1,nb2,3)=te(nb1,nb2,3)+psir(ind,nb1)*z1*psir(ind,nb2)
            
          END IF
        END DO
      END DO
    END DO
    te(nb1,nb2,1)=ABS(te(nb1,nb2,1)*dvol)**2D0
    te(nb1,nb2,2)=ABS(te(nb1,nb2,2)*dvol)**2D0
    te(nb1,nb2,3)=ABS(te(nb1,nb2,3)*dvol)**2D0
    
    IF(te(nb1,nb2,1) >= 1D0 .OR. te(nb1,nb2,2) >= 1D0  &
          .OR. te(nb1,nb2,3) >= 1D0) THEN
      WRITE(34,'(a,i3,a,i3,a,3f12.5,1e17.7)') 'wf1:',nb1,' wf2:',nb2,  &
          ' mte-x,y,z=',(te(nb1,nb2,i),i=1,3),amoy(nb1)-amoy(nb2)
    END IF
  END DO
END DO


RETURN
END SUBROUTINE transel





!     ******************************

SUBROUTINE reocc()

!     ******************************


!     readjust occupations numbers to actual s.p. energies

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
REAL(DP) :: occold(kstate),ocwork(kstate)
REAL(DP) :: ph(ksttot)             ! degeneracy of wavefunction, for
INTEGER :: nelecs(2)

!--------------------------------------------------------------------

DO nbe=1,nstate
  amoy(nbe)   = ekinsp(nbe)+epotsp(nbe)
  occold(nbe) = occup(nbe)
END DO
epstmp = 1D-8
occo = 1D0-occmix
IF(numspin==2) THEN
  nelecs(1) = nclust-nspdw
  nelecs(2) = nspdw
  DO is=1,numspin
    DO nbe=1,nstate
      IF(ispin(nbe) == is) THEN
        ph(nbe) = 1D0
      ELSE
        ph(nbe) = 0D0
      END IF
    END DO
    eferm  = 0D0         ! restart search of eferm from scratch
    CALL pair(amoy,ocwork,ph,nelecs(is),nstate,gp,eferm,temp,  &
        partnm,200,4,epstmp,-1,ksttot)
    DO nbe=1,nstate
      IF(ispin(nbe) == is) THEN
        occup(nbe) = occmix*ocwork(nbe)*ph(nbe) + occo*occold(nbe)
      END IF
    END DO
  END DO
ELSE
  eferm  = 0D0         ! restart search of eferm from scratch
  CALL pair(amoy,ocwork,ph,nclust,nstate,gp,eferm,temp,partnm,  &
      60,4,epstmp,-1,ksttot)
  DO nbe=1,nstate
    occup(nbe) = occmix*ocwork(nbe)*ph(nbe) + occo*occold(nbe)
  END DO
END IF
RETURN
END SUBROUTINE reocc


