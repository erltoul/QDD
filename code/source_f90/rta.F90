!MV
!_____________________________rta_________________________________________
SUBROUTINE rta(psi,aloc,rho,iterat)

USE params, ONLY: DP, kdfull2,kstate,nstate,dvol,ispin,ifsicp,outnam,&
                  nyf,nxyf,centfx,nx2,ny,nz,enonlo,&
                  amoy,energy,rtasumvar2max,occup,hbar,&
                  jrtaint,dt1,eqstnspup,eqstnspdw,ekinsp,psitophi,tfs, &
                  elasermask,rtasigee,rtars
USE kinetic
IMPLICIT NONE

COMPLEX(DP), INTENT(IN OUT)  :: psi(kdfull2,kstate)
REAL(DP), INTENT(IN OUT)     :: aloc(2*kdfull2),rho(2*kdfull2)
INTEGER, INTENT(IN)          :: iterat

COMPLEX(DP)       :: psiorth(kdfull2,kstate),psieq(kdfull2,kstate)
COMPLEX(DP)       :: rhopsieq(nstate,nstate),scal(nstate,nstate)
REAL(DP)          :: occuppsi(kstate),occuporth(kstate),Estar,trace, sigee, rs,oneontaurelax,Nel0,Nel
REAL(DP)          :: Occspin(2),Occref(2),eta,Eperspinref(2),Eperspintarget(2),entropy,Eref
REAL(DP)          :: mu,T
INTEGER           :: i,j,ii,nspup,nspdw

sigee=rtasigee!6.5d0 (rta2)
rs=rtars!3.7(rta2)
        nspup=eqstnspup
        nspdw=eqstnspdw
        CALL dyn_mfield(rho,aloc,psi,0D0) 
        WRITE(*,*)'_____________before orth'
        CALL info(psi,rho,aloc,iterat)!to update parameters as energy per mode, etc
        occuppsi=occup!to store it
	CALL OccupPerSpin('start rta',Occref)!memorize for final correction
!        CALL calc_Eref(occup,ispin,amoy,Eperspinref)!memorize for final correction
	CALL srhomat(psi,aloc,psiorth,occuporth)! new basis is psi orth with occuporth
	occup=occuporth!because dyn_mfield needs occup
	Nel0=sum(occup(1:nstate))!number of electrons
        CALL dyn_mfield(rho,aloc,psiorth,0D0) !  to check after transfo
        WRITE(*,*)'_____________before eq: iterat=',iterat
        CALL info(psiorth,rho,aloc,iterat)
        Eref=energy
        WRITE(*,'(a,40f10.6)')'occup    ',(occup(i),i=1,nstate)
	WRITE(*,'(a,40f10.6)')'occuporth',(occuporth(i),i=1,nstate)
	WRITE(*,'(a,40i5)')'ispin',(ispin(i),i=1,nstate)
	CALL OccupPerSpin('start rta',Occref)!memorize for final correction
	
	!compute the T distribution of occupation numbers
 CALL eqstate(psiorth,aloc,rho,psieq,occuporth,iterat)!find density current constained state of same energy(T increased)
 CALL OccupPerSpin('end eqstate',Occspin)

 CALL occupT0(occup(1:nstate),amoy(1:nstate),Estar)! computes Estar, occup is unchanged

WRITE(*,*)'_____________________ after eqstate_______________________________'
        CALL dyn_mfield(rho,aloc,psieq,0D0) 
        CALL analyze_elect(psieq,rho,aloc,iterat)!update field and print
Nel=0.d0;DO ii=1,nstate;Nel=Nel+occup(ii);ENDDO!sum(occup)! compute new number of electron, used in relaxation time calculation
!compose the dcmf state with the psiorth
oneontaurelax=0.4d0*sigee/rs**2*Estar/Nel!inverse of relaxation time
eta=dt1*oneontaurelax*jrtaint !eta for rta scheme
WRITE(*,*)'dt1,Nel,oneontaurelax,jrtaint',dt1,Nel,oneontaurelax,jrtaint

!compose rhoorth and rhoeq
!spin up
 CALL calcrhoeq(psiorth(:,1:nspup),psieq(:,1:nspup),psi(:,1:nspup), occuporth(1:nspup)&
               ,occup(1:nspup),nspup)!compute the rho eq operator for next loop of iterations, modifies psiorth->psi, occup
               WRITE(*,*)'just after rhoeq'
!spin down
 CALL calcrhoeq(psiorth(:,nspup+1:nstate),psieq(:,nspup+1:nstate),psi(:,nspup+1:nstate), occuporth(nspup+1:nstate)&
               ,occup(nspup+1:nstate),nspdw)!compute the rho eq operator for next loop of iterations, modifies psiorth->psi, occup
 CALL OccupPerSpin('end calcrhoeq',Occspin)
!   WRITE(*,*)'electron loss',(sum(occup(1:nstate))-Nel0)/Nel0 !unnecessary, number of electrons adjusted in calcrhoeq
trace=0.d0
DO i=1,nstate
 trace=trace+real(scal(i,i))
ENDDO
WRITE(*,'(a,10f12.8)') 'Estar,trace, oneontaurelax Nel eta',Estar,trace,oneontaurelax,Nel,eta
!compute energies, etc
WRITE(*,*)'_____________________energies after rhoeq_______________________________'
       CALL dyn_mfield(rho,aloc,psi,0D0) !update field to new density
       CALL analyze_elect(psi,rho,aloc,iterat)!to have energy
WRITE(*,*) 'occup after rta',sum(occup(1:nstate))
         !CALL CorrectEnergy(energy,Eref)! one first CALL to correct energy
        CALL calc_Eref(occup,ispin,amoy,Eperspinref)
        Eperspintarget=Eperspinref+(Eref-energy)*Eperspinref/sum(Eperspinref)
        CALL CorrectEnergy2(occref(1),Eperspintarget(1),occup(1:nspup),amoy(1:nspup),occup(1:nspup),nspup)
        CALL CorrectEnergy2(occref(2),Eperspintarget(2),occup(nspup+1:nstate),amoy(nspup+1:nstate),occup(nspup+1:nstate),nspdw)

WRITE(*,*)'_____________________energies after Correct_Energy________________________'
       CALL dyn_mfield(rho,aloc,psi,0D0) !update field to new density
       CALL analyze_elect(psi,rho,aloc,iterat)!to have energy
         !CALL CorrectEnergy(energy,Eref)! one second CALL to correct energy
        CALL calc_Eref(occup,ispin,amoy,Eperspinref)
        Eperspintarget=Eperspinref+(Eref-energy)*Eperspinref/sum(Eperspinref)
        CALL CorrectEnergy2(occref(1),Eperspintarget(1),occup(1:nspup),amoy(1:nspup),occup(1:nspup),nspup)
        CALL CorrectEnergy2(occref(2),Eperspintarget(2),occup(nspup+1:nstate),amoy(nspup+1:nstate),occup(nspup+1:nstate),nspdw)
        CALL OccupPerSpin('first correct nergy',Occspin)
WRITE(*,*)'_____________________energies after second Correct_Energy________________________'
       CALL dyn_mfield(rho,aloc,psi,0D0) !update field to new density
       CALL analyze_elect(psi,rho,aloc,iterat)
       CALL OccupPerSpin('second correct nergy',Occspin)
!compute and store entropy
entropy=0.d0
DO i=1,nstate
  IF((occup(i)>1.d-15).and.(occup(i).lt.(1.d0-1.d-15)))THEN
   entropy=entropy+occup(i)*log(occup(i))+(1-occup(i))*log(1-occup(i))
 ENDIF
 !compute temperature
ENDDO
 CALL temperature(mu,T)
      OPEN(1002,POSITION='append',FILE='prta')
      WRITE(1002,'(11f14.7)') tfs,entropy, elasermask,mu,T
      close(1002)

!stop!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
!_________________________SUBROUTINE calcrhoeq()_____________________________
SUBROUTINE calcrhoeq(psiorthloc,psieqloc,psiloc,occuporthloc,occuploc,nstateloc)
    !to make a sum of psiorth(occuporth)*(1-eta) with psieq(occup)*eta projected 
    !and diagonalize
    !this is done in a basis of 2*nstate
    !then conserve only the nstate most populated states
    USE util, ONLY:wfovlp
    IMPLICIT NONE

    INTEGER::nstateloc,k,l
    COMPLEX(DP)::psiorthloc(kdfull2,nstateloc),psieqloc(kdfull2,nstateloc),psiloc(kdfull2,nstateloc)
    REAL(DP)::occuporthloc(nstateloc),occuploc(nstateloc),coeffoccup

    

    COMPLEX(DP)::scal0(2*nstateloc,2*nstateloc),scal1(2*nstateloc,2*nstateloc),tab(kdfull2,2*nstateloc)
    COMPLEX(DP)::taborth(kdfull2,2*nstateloc)
    COMPLEX(DP)::occloc(2*nstateloc,2*nstateloc)
    COMPLEX(DP)::uloc(2*nstateloc,2*nstateloc),vloc(2*nstateloc,2*nstateloc),vploc(2*nstateloc,2*nstateloc)
    COMPLEX(DP)::wloc(2*nstateloc,2*nstateloc),D1onrholoc(2*nstateloc,2*nstateloc)!for the diagonalisation process
    REAL(DP)::eigenloc(2*nstateloc)
    
    INTEGER indx(2*nstateloc),i,j,ncall
 WRITE(*,*) 'in calcrhoeq, nstateloc',nstateloc
    !ALLOCATE matrices
    !tab is the table of wave functions
    tab(:,1:nstateloc)=psiorthloc(:,1:nstateloc)
    tab(:,nstateloc+1:2*nstateloc)=psieqloc(:,1:nstateloc)
 	!occloc diagonal matrix of occupation numbers
	occloc=0.d0
	DO i=1,nstateloc
	     occloc(i,i)=(1-eta)*occuporthloc(i)
	     occloc(nstateloc+i,nstateloc+i)=eta*occuploc(i)
	ENDDO
	!compute scal product <psi|psi>
	DO j=1,2*nstateloc
	  DO i=1,2*nstateloc
	      scal1(i,j)=wfovlp(tab(:,i),tab(:,j))!no spin involved, each call to rhoeq in a subspace of spin
	  ENDDO
	ENDDO
    !compute occup in the basis psiorth, psi
    occloc=matmul(transpose(conjg(scal1) ),occloc)
    occloc=matmul(occloc,scal1) 
 WRITE(*,*) 'in calcrhoeq bef diag'
   !diagonalize scalar product 
    CALL cdiagmat(scal1, eigenloc, vloc, 2*nstateloc)
 WRITE(*,*) 'in calcrhoeq aft diag'
   !compute diagonal matrices D1on_rho
    D1onrholoc=0.D0
    DO i=1,2*nstateloc
		D1onrholoc(i,i)=1/sqrt(eigenloc(i))
    ENDDO
   !scale v into vp to have a unity matrix
    vploc=matmul(vloc,D1onrholoc)
   !compute Vp T *occloc*Vp
    scal0=matmul(transpose(conjg(Vploc) ),occloc)
    scal0=matmul(scal0,vploc )
   !diagonalize Vp T *rhomat*Vp: eigenvectors in u
 WRITE(*,*) 'in calcrhoeq bef diag2'
    CALL cdiagmat(scal0,eigenloc,uloc,2*nstateloc)
   !compute the new change of basis matix
 WRITE(*,*) 'in calcrhoeq aft diag2'
    wloc=matmul(vploc,uloc)
	!compute the phi_j in the r basis: wloc are the taborth in the tab basis
	!here phi is in table tab
	taborth=0.d0
	DO j=1,2*nstateloc
	  DO i=1,2*nstateloc
	  taborth(:,j)=taborth(:,j)+wloc(i,j)*tab(:,i)
	  ENDDO
	ENDDO
   
   coeffoccup=sum(eigenloc(1:2*nstateloc))
   !ordonnate the occupation numbers
    ncall=2*nstateloc
    CALL indexx(ncall,eigenloc,indx) 
 WRITE(*,*) 'in calcrhoeq aft indxx'
   WRITE(*,*)'eigensorted'
   WRITE(*,'(10f14.10)') (eigenloc(indx(i)),i=2*nstateloc,1,-1)
!    WRITE(*,*) 'ispin '
!    WRITE(*,'(40i3)')(ispinloc(i),i=1,2*nstateloc)
!   WRITE(*,*) 'ispin sorted'
!   WRITE(*,'(40i3)')(ispinloc(indx(i)),i=2*nstateloc,nstateloc+1,-1)
   !select and store in psi the nstateloc modes with higher occupation numbers
   j=0
   DO i=2*nstateloc,nstateloc+1,-1
          j=j+1
          psiloc(:,j)=taborth(:,indx(i))
          occuploc(j)=eigenloc(indx(i))
   ENDDO
   l=0
!   coeffoccup=coeffoccup/sum(occuploc(1:nstateloc))!compute a scaling factor for occupation numbers
!   occup(1:nstateloc)=occup(1:nstateloc)*coeffoccup!scale occupation numbers by this factor
!   WRITE(*,*)'1-coeffoccup, occuptot'
!   WRITE(*,'(12f14.10)')1.d0-coeffoccup,sum(occup(1:nstateloc))

END SUBROUTINE calcrhoeq

END SUBROUTINE rta

!________________________________________cdiagmat________________________________________________
SUBROUTINE cdiagmat(mat, eigen, vect, N)
!mat	:matrice to be diagonalised.
!eigen 	:vector of eigenvalues
!vect  	:eigenvectors written in columns
!N 		:dimension, assuming that all states are occupied
IMPLICIT NONE
INTEGER,PARAMETER :: DP=KIND(1D0)  ! precision  setting
COMPLEX(DP), INTENT(IN)                  	:: mat(N,N)
REAL(DP), INTENT(OUT)                   :: eigen(N)
COMPLEX(DP),INTENT(OUT)					:: Vect(N,N)
COMPLEX(DP)::mat0(N,N)
INTEGER ::N,i,j
vect=0.d0
eigen=0.d0
mat0=mat
!CALL CDIAG(mat(1:nspup,1:nspup),eigen(1:nspup),vect(1:nspup,1:nspup),nspup,nspup)   
!CALL CDIAG(mat(nspup+1:N,nspup+1:N),eigen(nspup+1:N),vect(nspup+1:N,nspup+1:N),nspdw,nspdw)   
CALL HEigensystem(N,mat0,N,eigen,vect,N,0) 
vect=transpose(conjg(vect))
   
END SUBROUTINE cdiagmat


!_________________________________Calc-Eref__________________________________________
!compute the single particule energy per spin, for use in the fermi convergence
SUBROUTINE calc_Eref(occup,ispin,Ei,Eref)
use params,only:DP,nstate
IMPLICIT NONE
REAL(DP), intent (in):: occup(nstate),Ei(nstate)
REAL(DP), intent(out):: Eref(2)
INTEGER, intent (in) ::ispin(nstate)
INTEGER:: i
Eref=0.d0
 DO i=1,nstate
    Eref(ispin(i))=Eref(ispin(i))+occup(i)*Ei(i)
 ENDDO
END SUBROUTINE calc_Eref
!_________________________________eqstate____________________________________________

SUBROUTINE eqstate(psi,aloc,rho,psi1,occuporth,iterat)
!update psi into psi1, and occuporth(entrance) into occup(in params)
USE params, only: DP, kdfull2,kstate,nstate,dvol,ispin,ifsicp,outnam,&
                  nyf,nxyf,centfx,nx2,ny,nz,enonlo,&
                  amoy,energy,rtamu,rtamuj,rtasumvar2max,&
                  occup,eqstnspup,eqstnspdw,rtatempinit
USE kinetic
USE util, ONLY: prifld
IMPLICIT NONE
COMPLEX(DP), INTENT(IN OUT)  :: psi(kdfull2,kstate),psi1(kdfull2,kstate)
REAL(DP), INTENT(IN OUT)     :: aloc(2*kdfull2),rho(2*kdfull2),occuporth(kstate)
INTEGER, INTENT(IN) :: iterat

COMPLEX(DP)                  :: dpsi(kdfull2)

REAL(DP)::rhotot0(kdfull2,2),rhotot1(kdfull2,2),err,ma0,ma1,err1
REAL(DP)::lambda(kdfull2,2),mu,time0,time1
REAL(DP)::curr0(kdfull2,3),curr1(kdfull2,3),lambdaj(kdfull2,3),errj,muj!parameters for current
REAL(DP)::sumvar2,eal,fac,EspPerMod(nstate),EspTotRef
REAL(DP)::EspPerSpinRef(2),EspPerSpinAchieved(2),EspPerSpinTarget(2),temp(2),mut,EspTotAchieved
INTEGER::i,j,ishift,ii,nspup,nspdw

WRITE(*,*) "EQSTATE: iterat=",iterat
nspup=eqstnspup
nspdw=eqstnspdw
err=1.d3
!initialize psi, rhotot,lambda,j lamdaj
 sumvar2=1.d0!a big value, to start minimisation
 psi1=psi
 lambda=0.d0!a set of lagrange coefficients, for spin up and down, associated with density
 lambdaj=0.d0!a set of lagrange coefficients, associated to j
 mu=rtamu
 muj=rtamuj
 CALL calcrhotot(rhotot0,psi)
 rhotot1=rhotot0
 CALL calc_current(curr0,psi)
 curr1=curr0
 ma0=dvol*sum(rhotot0)
 CALL prifld(curr0(:,1),'j_1')
 CALL prifld(curr0(:,2),'j_2')
 CALL prifld(curr0(:,3),'j_3')
WRITE(*,*) '  j      centfx       err           errj          mu        muj           &
            ma0           ma1'
j=0
!once for all compute aloc, without the local part
CALL coul_mfield(rho)


    CALL calc_Eref(occuporth,ispin,amoy,EspPerSpinRef)!compute a reference sp energy per spin,with occuporth
    EspTotRef=EspPerSpinRef(1)+EspPerSpinRef(2)
 WRITE(*,*) 'appel initial fermi'
    CALL fermi1(amoy,EspPerSpinRef(1),occup,1,rtatempinit/10,1.0,temp(1),mut)!occup changed in Fermi
    CALL fermi1(amoy,EspPerSpinRef(2),occup,2,rtatempinit/10,1.0,temp(2),mut)!occup changed in Fermi
    EspPerSpinTarget=EspPerSpinRef!initial value of EspPerSpinTarget
    EspTotAchieved=0.d0!to avoid exit at first pass
!The big issue:
! when converging, some modes of low weght are crossing, nd getting low energy.
!at next fermi the density is changed a lot and err and sumvar are exploding. Also energy is ver out of balance
!in order to control this this version introduces an additionnalloop, with little authority (0.001 on Temp)
! and evry tenth step at most, in order to reorder the Occ vs Esp     
   WRITE(177,'(a,i5/a)') 'RTA protocol at iterat=',iterat,&
    '  j  err   errj    diffrho  sumj0 sumj1  diffE  sumvar2 '
   DO while(sumvar2>rtasumvar2max.or.abs(EspTotRef-EspTotAchieved)>1.d-4.or.err>0.2d0)!main loop
        IF((sumvar2<1.d-1.and.err<1.0).and.mod(j,10).eq.0)then !update the occupation numbers without changing the target energy
                !idea: keep the occ in line with  the energies
                !at high temp and temp rate, the energy level are crossing
    		CALL fermi1(EspPerMod,EspPerSpinTarget(1),occup,1,max(temp(1)-0.001,rtatempinit/10.),temp(1)+0.001,temp(1),mut)!
    		CALL fermi1(EspPerMod,EspPerSpinTarget(2),occup,2,max(temp(2)-0.001,rtatempinit/10.),temp(2)+0.001,temp(2),mut)!
        ENDIF

	IF((sumvar2<1.d-2.and.err<.5).or.err<1.d-1) THEN!if there is enough convergence update fermi occupation numbers
 		EspPerSpinAchieved=0.d0
  		DO ii=1,nstate
   			EspPerSpinAchieved(ispin(ii))=EspPerSpinAchieved(ispin(ii))+occup(ii)*EspPerMod(ii)!err on ener per sp tot per spin with ekmod computed after cschmidt (in calc psi1)
  		ENDDO
    		EspPerSpinTarget=EspPerSpinTarget-(EspPerSpinAchieved-EspPerSpinRef)/2.0d0
    		WRITE(*,'(a,20f14.7)') 'EspPerSpinAchieved-EspPerSpinRef',EspPerSpinAchieved(1)-EspPerSpinRef(1)&
    		,EspPerSpinAchieved(2)-EspPerSpinRef(2)
    		WRITE(*,'(a,20f14.7)') 'EspPerSpinTarget-EspPerSpinRef',EspPerSpinTarget(1)-EspPerSpinRef(1)&
    		,EspPerSpinTarget(2)-EspPerSpinRef(2)
   		WRITE(*,*) 'fermi dans la boucle'
    		CALL fermi1(EspPerMod,EspPerSpinTarget(1),occup,1,max(temp(1)-0.005,rtatempinit/10.),temp(1)+0.005,temp(1),mut)!call with some phase advance, using error on ek occup changed in Fermi
    		CALL fermi1(EspPerMod,EspPerSpinTarget(2),occup,2,max(temp(2)-0.005,rtatempinit/10.),temp(2)+0.005,temp(2),mut)!call with some phase    	ELSE
   		EspPerSpinTarget=EspPerSpinRef
   	ENDIF
 j=j+1! counts the number of loops
!next value of psi1
  CALL calc_psi1(psi1,aloc,rhotot0,rhotot1,curr0,curr1,j,lambda,mu,lambdaj,muj,sumvar2,eal,EspPerMod)! one step of damped gradient
  EspTotAchieved=0.d0; DO ii=1,nstate; EspTotAchieved= EspTotAchieved+&
               EspPerMod(ii)*occup(ii);ENDDO
  err=sum(abs(rhotot1-rhotot0))
  errj=sum(abs(curr1-curr0))
  ma1=dvol*sum(rhotot1)
  lambda=lambda+2.d0*(rhotot1-rhotot0)*mu!update lambda(Augmented Lagrangian rule)
  lambdaj=lambdaj+2.d0*(curr1-curr0)*muj
  IF(mod(j,1).eq.0) THEN
      WRITE(*,'(a,i6,20f14.7)') ' RTA-iteration:',j,centfx,err,errj,mu,muj,ma0,ma1,EspTotRef,EspTotAchieved&
      ,EspTotRef-EspTotAchieved,sumvar2,EspPerSpinAchieved(1)-EspPerSpinAchieved(2)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     CALL flush(6)
     WRITE(177,'(a,i6,20(1pg14.6))') ' RTA-iteration:',j,err,errj,ma0-ma1,&
      sum(abs(curr0)),sum(abs(curr1)),EspTotRef-EspTotAchieved,sumvar2
     CALL flush(177) 
  ENDIF 
  IF(j>100) STOP
ENDDO!DO j=1...
OPEN(1000,STATUS='unknown',POSITION='append',FORM='formatted',FILE='pspeed.'//outnam)
        WRITE (1000,'(6g14.5)') (rhotot0(ny*(nyf+nxyf)+ii,1),&
          rhotot0(ny*(nyf+nxyf)+ii,2),rhotot1(ny*(nyf+nxyf)+ii,1),rhotot1(ny*(nyf+nxyf)+ii,2),&
          curr0(ny*(nyf+nxyf)+ii,1),curr1(ny*(nyf+nxyf)+ii,1),ii=1,nx2)
          WRITE(1000,*)
 close(1000)
        CALL cpu_time(time1)
        WRITE(*,*)'j,time', j, (time1-time0)
        CALL dyn_mfield(rho,aloc,psi1,0D0) 
       WRITE(*,*)'_____________in eqstate 2: iterat=',iterat
        CALL info(psi1,rho,aloc,iterat)
      WRITE(*   ,'(2i4,11f14.7)') iterat,j,sumvar2,err,errj,mu,muj,energy
      OPEN(1001,POSITION='append',FILE='peqstate')
      WRITE(1001,'(2i6,11f14.7)') iterat,j,sumvar2,err,errj,mu,muj,energy
      close(1001)

END SUBROUTINE eqstate
!_________________________________OccupT0______________________________________________
SUBROUTINE OccupT0(occloc,esploc,Estar)
!this routine computes occupation numbers at temperature=0
!for each spin, ordonate energies, putting occupation numbers to 1, until only a value between 0,1 is left for the last occupation number
! then computes Estar. The calling occupation numbers are not changed
USE params, only:DP,ispin,nstate
IMPLICIT NONE
INTEGER::isp,i,n,indx(nstate),order(nstate)
REAL(DP)::occloc(nstate),esploc(nstate),Estar
REAL(DP)::Etot,OccTot,EspSpin(nstate),OccSpin(nstate),EspSpinT0,occup
! first compute tables of energies per spi
!here there is no assumption thar spin are in order
EspSpinT0=0.d0!the ergy at T=0
Etot=dot_product(Occloc(1:nstate),Esploc(1:nstate))!cumulate enrgies of both spins
DO isp=1,2! for each spin
   n=0
   DO i=1,nstate
      IF (ispin(i).eq.isp) THEN    
      n=n+1!count spin isp elements
        OccSpin(n)=occloc(i)!sp occupation per spin
        EspSpin(n)=EspLoc(i)!sp energy per spin
        indx(N)=i
      ENDIF
   ENDDO
   OccTot=sum(Occspin(1:N))
   !now we have energies per spin in ener
   CALL indexx (n,EspSpin,order) !order contains the indices of energies in increasing order 
   !and a second loop to populate the occupation numbers
   DO i=1,n!a loop for this acting spin
     IF(OccTot>=1.d0) THEN
       occup=1.d0
       OccTot=OccTot-1.d0
       ELSEIF(Occtot<1.d0.and.OccTot>=0) THEN
       occup=OccTot
       OccTot=OccTot-1.d0
       ELSE
       occup=0.d0
       OccTot=OccTot-1.d0
     ENDIF
     EspSpinT0=EspSpinT0+occup*EspLoc(indx(order(i)))
     WRITE(*,'(a,4f14.8)') 'occup, Esp', occup,occloc(indx(order(i))),EspLoc(indx(order(i)))
   ENDDO
ENDDO
Estar=Etot-EspSpinT0!stop
WRITE(*,'(a,3f14.7)') 'estar,Etot,EspSpinT0',estar,Etot,EspSpinT0
END SUBROUTINE OccupT0
!_______________________________Calc_psi1____________________________________________
SUBROUTINE Calc_psi1(psi1,aloc,rhotot0,rhototloc,curr0,curr1,j,lambda,mu,lambdaj,muj,sumvar2,eal,ekmod)
USE params, only: DP, kdfull2,kstate,nstate,dvol,ispin,eye,ipsptyp,h2m,hbar,ame,pi,nx2,ny2,nz2,&
                  rtaeps,rtae0dmp,occup,enonlo,nxyz
USE kinetic
USE util, ONLY: wfovlp
IMPLICIT NONE
COMPLEX(DP), INTENT(IN OUT)                  :: psi1(kdfull2,kstate)
REAL(DP), INTENT(IN OUT)                     :: rhotot0(kdfull2,2),ekmod(nstate)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: rhototloc(kdfull2,2)
REAL(DP), INTENT(IN OUT)                     :: lambda(kdfull2,2),mu
REAL(DP), INTENT(IN OUT)                     ::curr0(kdfull2,3),curr1(kdfull2,3)
REAL(DP), INTENT(IN OUT)                     ::lambdaj(kdfull2,3),muj,eal!parameters for current

INTEGER, INTENT (IN)::j
COMPLEX(DP) :: psik(kdfull2),ALpsi(kdfull2),hpsi(kdfull2,kstate),psi1new(kdfull2,kstate)
COMPLEX(DP) :: vpsi(kdfull2),HpsiN(kdfull2),Apsi(kdfull2)
COMPLEX(DP) ::cnl,scal(nstate,nstate),hpsip(kdfull2,kstate)
REAL(DP) ::epswf,e0dmp,deltalambda,emod(nstate),epot(nstate)
REAL(DP) ::sumvar2,dkvol,eee,ekin
INTEGER N,ii,ishift,ii1

dkvol=((2*pi)**3)*dvol!I do not understand this coefficient, it is used in calc nn, and works
 epswf=rtaeps
 !WRITE (*,*)epswf
 e0dmp=rtae0dmp
 eal=0.d0
    !compute real space part of augmented hamiltonien,transfer in k space
    DO N=1,nstate
!IF(occup(N)>1.d-8) THEN! this is added since I had bad convergenece when non ocupied states where used. Should work as no states do not contribute to rho or j
   CALL calc_hamiltonien
     !in k space
      CALL fftf(psi1(1,N),psik)!psik=psi in k space
      CALL fftf(ALpsi,ALpsi)!Alpsi here in fact all the potential + penalty in k space
      CALL fftf(vpsi,vpsi)!vpsi the potential in the k space
           !penalty +real energy
                Apsi=h2m*akv*psik!kin energy in ksce
                vpsi=vpsi+Apsi!T+V hpsi in k space
                ekmod(N)=real(dot_product(vpsi,psik))*dkvol!energy per mode
                Apsi=Apsi+ALpsi!Apsi unused here to compute ek + tot penalty +v
                emod(N)=real(dot_product(Apsi,psik))*dkvol!energy +AL per mode
                hpsip(:,N)=Apsi!store hpsi total in k space for calc of sumvar
                Apsi=Apsi-emod(N)*psik

     !gradient step in k space
      psik=psik-epswf/(e0dmp+akv)*(Apsi)
      CALL fftback(psik,psi1(1,N))!here psinew just to DO printings with homogeneous psi, hpsi
      eal=eal+emod(N)*occup(N)
!ENDIF
    ENDDO!N
	CALL calc_var(hpsip,psi1,sumvar2 )!hppsip in k space
IF ((mod(j,1).eq.0)) THEN
	WRITE(*,'(a,i6,(20f10.6))')'emod1 ener+Lagrangi. j=',j,(emod(N),N=1,nstate)
	WRITE(*,'(a,(21f10.6))')'ekmod1 ener per mode',(ekmod(N),N=1,nstate)
	WRITE(*,'(a,(21f10.6))')'occup',(occup(N),N=1,nstate)
!	WRITE(*,'(a,f14.9)')'sumvar2 penalty',sumvar2
ENDIF                     
 CALL cschmidt(psi1)
 CALL calcrhotot(rhototloc,psi1)
 CALL calc_current(curr1,psi1)
 DO N=1,nstate
  CALL calc_ekin(psi1(1,N),ekmod(N))
  ishift = (ispin(N)-1)*nxyz
!  WRITE(*,*) ' before CALC_EPOT',N,1,ishift
!  CALL calc_epot(psi1(1,N),aloc(ishift+1),epot(N),enonlo(N),N,1)
vpsi=0.d0           
IF (ipsptyp==1) THEN!use local potential and non local only if there is a non local potential
             CALL  nonlocalc(psi1(1,N),vpsi,0)  !vpsi is the term Vnonloc *psi  
           ENDIF
             !add local potential
             DO ii=1,kdfull2
                vpsi(ii)=vpsi(ii)+aloc(ishift+ii) *psi1(ii,N)
             ENDDO
epot(N)=0.d0
DO ii=1,kdfull2
 epot(N)=epot(N)+conjg(psi1(ii,N))*Vpsi(ii)*dvol
ENDDO
  ekmod(N)=ekmod(N)+epot(N)!ekmod is total sp energy
 ENDDO
contains 
!________________________________________Calc_hamiltonien______________________________________________________________________
 SUBROUTINE calc_hamiltonien
   COMPLEX(DP),ALLOCATABLE::A(:),grad(:),gradApsi(:),jgradpsi(:)
   REAL(DP)::hbarm
   ALLOCATE(A(kdfull2),grad(kdfull2),gradApsi(kdfull2),jgradpsi(kdfull2))
      ishift=(ispin(N)-1)*kdfull2
      hbarm=hbar/ame
      vpsi=0.d0

     ! compute lambda*psi
           ALpsi=lambda(:,ispin(N))*psi1(:,N)!real space
     ! compute rho*psi
           ALpsi=ALpsi+mu*(rhototloc( :,ispin(N) )-rhotot0( :,ispin(N) ) )*psi1(:,N)!real space
     ! compute v*psi
           IF (ipsptyp==1) THEN!use local potential and non local only if there is a non local potential
             CALL  nonlocalc(psi1(1,N),vpsi,0)  !vpsi is the term Vnonloc *psi  
           ENDIF
            !aloc is a pure local potential 
            !add local potential
             DO ii=1,kdfull2
                vpsi(ii)=vpsi(ii)+aloc(ishift+ii) *psi1(ii,N)
             ENDDO
 !           ENDIF
     !compute (muj (j-j0)+lambda) nabla phi+nabla((muj (j-j0)+lambda) phi)
        A=eye*(2*muj* (curr1(:,1)-curr0(:,1))+lambdaj(:,1) )
        Apsi=A(:)*psi1(:,N)
        CALL xgradient_rspace(psi1(1,N),grad)
        CALL xgradient_rspace(Apsi,gradApsi)
        jgradpsi=         A(:)*grad(:)+gradApsi
        A=eye*(2*muj* (curr1(:,2)-curr0(:,2))+lambdaj(:,2) )
        Apsi=A(:)*psi1(:,N)
        CALL ygradient_rspace(psi1(1,N),grad)
        CALL ygradient_rspace(Apsi,gradApsi)
        jgradpsi=jgradpsi+A(:)*grad(:)+gradApsi
        A=eye*(2*muj* (curr1(:,3)-curr0(:,3))+lambdaj(:,3) )
        Apsi=A(:)*psi1(:,N)
        CALL zgradient_rspace(psi1(1,N),grad)
        CALL zgradient_rspace(Apsi,gradApsi)
        jgradpsi=jgradpsi+A(:)*grad(:)+gradApsi
        !total hamiltonian term
        ALpsi=vpsi+ALpsi-hbarm*jgradpsi!total potential+penalty term in real space
 DEALLOCATE(A,grad,gradApsi,jgradpsi)
 END SUBROUTINE calc_hamiltonien
END SUBROUTINE Calc_psi1
!____________________________________calcrhospin______________________________
SUBROUTINE calcrhotot(rho,q0)

!     density 'rho' for complex or real wavefunctions 'q0'
!experimental, only for non parallel scheme. To be rewritten
USE params
IMPLICIT NONE


COMPLEX(DP), INTENT(IN) :: q0(kdfull2,kstate)         ! cPW
REAL(DP), INTENT(OUT) :: rho(kdfull2,2)
INTEGER nb,ishift,ind
!k initialize densities:
  rho=0D0

DO nb=1,nstate
  ishift = ispin(nb) ! store spin=2 in upper block
  DO ind=1,nxyz
    rho(ind,ishift)=rho(ind,ishift)+ occup(nb)*(CONJG(q0(ind,nb)))*q0(ind,nb)
  END DO
END DO
RETURN
END SUBROUTINE calcrhotot
!___________________________________calc_var__________________________________
SUBROUTINE calc_var(hpsi,psi1,sumvar2)
  use params, only :DP,kdfull2,nstate,kstate,ispin,occup
  USE kinetic
  USE util, ONLY: project,wfovlp,c_project
  IMPLICIT NONE
  COMPLEX(DP), INTENT(IN)                  :: psi1(kdfull2,kstate)
  COMPLEX(DP), INTENT(IN OUT)              :: hpsi(kdfull2,kstate)
  REAL(DP), INTENT(OUT)                    :: sumvar2

  REAL(DP)                                 :: evarsp2,evarsp2b
  INTEGER:: N
  COMPLEX(DP)                              :: qout(kdfull2),ovl
  ! hpsi is the term h psi in k space
  ! psi1 the wavefunctions (in real space)
!                WRITE(*,*)'in var'
sumvar2=0.d0
!pass hpsi in real space
DO N=1,nstate
    CALL fftback(hpsi(1,N),hpsi(1,N))
ENDDO

!do the calc     !?? PGR: is 'c_project' the right thing to call ??
DO N=1,nstate
    CALL cproject(hpsi(:,N),qout,ispin(N),psi1)
!    ovl=wfovlp(psi1(:,N),hpsi(:,N))
!    qout=hpsi(:,N)-psi1(:,N)*ovl
    evarsp2 =  SQRT( real( wfovlp(qout,qout) ) )
    sumvar2   = sumvar2 + occup(N)*evarsp2**2
!WRITE(*,'(a,2i4,4(1pg13.5))') "RTA CALC_VAR:",N,ispin(N),occup(N),evarsp2,&
!  SQRT( real( wfovlp(hpsi(:,N),hpsi(:,N)) ) )
ENDDO
!WRITE(*,'(a, 40i3)') 'ispin in sumvar2', (ispin(N),N=1,nstate)
!WRITE(*,'(a, 40f10.6)') 'occup in sumvar2', (occup(N),N=1,nstate)


END SUBROUTINE calc_var
!_______________________________Fermi1______________________________________________________
SUBROUTINE fermi1(ekmod,eref,occup,ispinact,T0i,T1i,t2,mu)
!changes occup to achieve the total sp energy of eref, and the right total occupation
use params, only: DP, kstate,nstate,eqstnspup,eqstnspdw,rtatempinit
IMPLICIT NONE
   REAL(DP)::ekmod(kstate),occup(kstate),Eref,occref(2),T(2)
   REAL(DP)::Et,mu,T0i,T1i,T0,T1,t2,occtot(2)
   INTEGER ispinact,nspup,nspdw
   nspup=eqstnspup
   nspdw=eqstnspdw
   occref(1)=sum(occup(1:nspup))!total occupation spin1
   occref(2)=sum(occup(nspup+1:nstate))
   !T0, T1 are low and high boudaries of temp
   T0=T0i
   T1=T1i
   !dichotomy on T
   DO while(abs(T1-T0)>1.d-10)
       T2=(T1+T0)/2.d0
       IF (ispinact.eq.1) THEN
           CALL occT1(occref(1),ekmod(1:nspup),et,mu,occtot(1),nspup,T2,occup(1:nspup))! change mu at constatnt T2 to have the required total occupation
         ELSE
           CALL occT1(occref(2),ekmod(nspup+1:nstate),et,mu,occtot(2),nspdw,T2,occup(nspup+1:nstate))
       ENDIF
       IF(et>eref) THEN !et is the energy at T2, with the mu giving the right occupation number
          T1=T2
       ELSE 
          T0=T2
       ENDIF
   ENDDO
   WRITE(*,*) 'fermi1, T mu',T2,mu
END SUBROUTINE fermi1
!_______________________________________forceTemp____________________________________
SUBROUTINE forceTemp(amoy,occup,n,temp,mu)
use params,only:DP
IMPLICIT NONE
INTEGER n
REAL(DP)::amoy(n),occup(n),temp,mu,occtot,Etot,occtotloc
occtot=sum(occup(1:n))
Etot=dot_product(occup(1:n),amoy(1:n))
WRITE(*,*) 'occtot,Etot',occtot,etot
CALL OccT1(occtot,amoy,Etot,mu,occtotloc,n,temp,occup)
WRITE(*,*)'force:T,mu,occtot,occtotloc',temp,mu,occtot,occtotloc
END SUBROUTINE forceTemp
!_________________________________OccT1______________________________________________
SUBROUTINE OccT1(occrefloc,enerloc,Etotloc,muloc,occtotloc,n,T,occuploc)
!optimize mu at given temperature to have 
!	the right occupation
!occrefloc: total occupation
!enerloc:    sp energy per mode
!Etotloc: resulting sp total energy
!muloc: energy chemical potential
!occtotloc: taget total occupation
! t: temperature
!occloc: return the occupation numbers
use params, only: DP, kstate,nstate
IMPLICIT NONE
   INTEGER::n,i
   REAL(DP)::occrefloc,enerloc(n),Etotloc,muloc,T,occuploc(n)
   REAL(DP)::muloc0,muloc1,occtotloc,fac
   INTEGER::orderloc(n)
   muloc0=minval(enerloc)
   muloc1=maxval(enerloc)
   DO while(abs(muloc0-muloc1)>1.d-12)!dichotomy on mu to have the right ocupation number at T2
      muloc=(muloc0+muloc1)/2.d0
      occtotloc=0.d0
      Etotloc=0.d0
      DO i=1,n
         occuploc(i)=1.d0/(1.d0+exp((enerloc(i)-muloc)/T))
         occtotloc=occtotloc+occuploc(i)
         Etotloc=Etotloc+occuploc(i)*enerloc(i)!Etotloc is the corresponding energy
      ENDDO
      IF(occtotloc<occrefloc) THEN
         muloc0=muloc
      ELSE
         muloc1=muloc
      ENDIF
   ENDDO    
      
END SUBROUTINE occT1

!________________________________Fermi_init_________________________________________________   
SUBROUTINE fermi_init(ekmod,T,occup,ispinact)
use params, only: DP, kdfull2,kstate,nstate,ispin
IMPLICIT NONE
       REAL(DP), INTENT(in) :: ekmod(nstate),T
       REAL(DP), INTENT(in out) :: occup(kstate)
       INTEGER, intent (in)::ispinact
       REAL(DP)::mu1,mu2,mu3,o1,o2,o3,O,Occ
       INTEGER::i
       O=0.d0
       DO i=1,nstate
          IF (ispin(i).eq.ispinact) O=O+occup(i)
       ENDDO
       
       mu1=-2
       mu2=0
       o1=occ(mu1,ekmod,T,occup,O,ispinact)
       o2=occ(mu2,ekmod,T,occup,O,ispinact)
       DO while(abs(mu1-mu2)>1.d-5)
          mu3=(mu1+mu2)/2.d0
          o3=occ(mu3,ekmod,T,occup,O,ispinact)
          IF(o3<0)THEN
             mu1=mu3
             o1=o3
          ELSE  
             mu2=mu3
             o2=o3
          ENDIF
       ENDDO
       DO i=1,nstate
         IF(ispin(i)==ispinact) occup(i)=1.d0/(1.d0+exp((ekmod(i)-mu3)/T))
       ENDDO
  END SUBROUTINE fermi_init
       real function occ(mu,ekmod,T,occup,O,ispinact)
       use params, only: DP,nstate,ispin,kstate
      IMPLICIT NONE
       REAL(DP), INTENT(in) :: ekmod(nstate),O
       REAL(DP), INTENT(in ) :: occup(kstate)
       REAL(DP)::mu,T
       INTEGER, intent (in)::ispinact
       REAL(DP)::mu1,mu2,mu3,o1,o2,o3
       INTEGER::i
       OCC=0.d0
       DO i=1,nstate
          IF(ispin(i)==ispinact) Occ=Occ+1.d0/(1.d0+exp((ekmod(i)-mu)/T))
       ENDDO
       Occ=Occ-O
       end function Occ


!__________________________________________srhomat____________________________________________________________________
SUBROUTINE srhomat(psi,aloc,psiorth,occuporth)
USE params, only: DP, kdfull2,kstate,nstate,ispin,nrel2abs,nxyz,occup,apnum,psitophi
IMPLICIT NONE
REAL(DP), INTENT(IN)            :: aloc(2*kdfull2)
COMPLEX(DP), INTENT(IN OUT)     :: psi(kdfull2,kstate),psiorth(kdfull2,kstate)
REAL(DP),intent(out)	        :: occuporth(kstate)


INTEGER:: i,j
COMPLEX(DP):: scal(nstate,nstate),scal0(nstate,nstate)
      COMPLEX(DP)::D1on_rho(nstate,nstate),rhomat(nstate,nstate)
      REAL(DP) :: eigen(nstate) ,eigen0(nstate)     ! eigenvalues
      COMPLEX(8) :: v(nstate,nstate),vp(nstate,nstate),u(nstate,nstate),w(nstate,nstate)  ! eigenvectors
      COMPLEX(8) :: winv(nstate,nstate)!inverse of w
      REAL(DP)    :: occup1(nstate)





!Calculate scalar product of the <psi_i |psi_j>
!not hermitian
 CALL scalar(psi,psi,scal,ispin, 'scal psi(i) psi(j)')
!compute the density matrix in the basis psi
   rhomat=0.d0
   DO i=1,nstate
       rhomat(i,i)=occup(i)
   ENDDO!i
scal0=matmul(transpose(conjg(scal) ),rhomat)
rhomat=matmul(scal0,scal )

!diagonalize scal
 CALL cdiagmat(scal, eigen0, v, nstate)
 WRITE(*,'(a,30 f12.7)')'eigen0',(eigen0(i),i=1,nstate)
!compute diagonal matrices D1on_rho
D1on_rho=0.D0
DO i=1,nstate
    D1on_rho(i,i)=1/sqrt(eigen0(i))
ENDDO
!scale v into vp to have a unity matrix
vp=matmul(v,D1on_rho)


!compute Vp T *rhomat*Vp
scal0=matmul(transpose(conjg(Vp) ),rhomat)
scal0=matmul(scal0,vp )



!diagonalize Vp T *rhomat*Vp: eigenvectors in u
 CALL cdiagmat(scal0,eigen,u,nstate)
!compute the new change of basis matix
w=matmul(vp,u)
winv=matmul (transpose ( conjg(u) ) ,transpose(conjg(v)))

!compute the phi_j in the r basis: Vect are the phi_j in the psi basis
psiorth=0.d0
DO j=1,nstate
  DO i=1,nstate
  psiorth(:,j)=psiorth(:,j)+w(i,j)*psi(:,i)
  ENDDO
ENDDO
 WRITE(*,'(a,30 f12.7)')'eigen',(eigen(i),i=1,nstate)
!uptage the change of basis matrix
psitophi=matmul(winv,psitophi)!matmul(psitophi,winv)

occuporth(1:nstate)=eigen

END SUBROUTINE srhomat

!______________________________scalar product___________________________________
SUBROUTINE scalar(tab1,tab2,scal,ispin, mess)
USE params, only: DP, kdfull2,nstate,kstate,nrel2abs

 USE util, ONLY: wfovlp

 IMPLICIT NONE

 COMPLEX(DP), INTENT(IN)       :: tab1(kdfull2,kstate),tab2(kdfull2,kstate)
 COMPLEX(DP), INTENT(OUT)      :: scal(nstate,nstate)


 INTEGER::ispin(nstate)
 CHARACTER(*)::mess
INTEGER:: i,j
DO j=1,nstate
  DO i=1, nstate
    IF( ispin(nrel2abs(i)) == ispin(nrel2abs(j) )) THEN
      scal(i,j)=wfovlp(tab1(:,i),tab2(:,j))
    ELSE
      scal(i,j)=0.d0
    ENDIF    
  ENDDO
ENDDO
WRITE(*,*) mess
DO i=1,nstate
    WRITE(*, '(20f12.9)')(scal(i,j),j=1,nstate)
ENDDO
END SUBROUTINE scalar

!________________________________________cdiagspin________________________________________________
SUBROUTINE cdiagspin(mat, eigen, vect, N)
!mat	:matrice to be diagonalised.
!eigen 	:vector of eigenvalues
!vect  	:eigenvectors written in columns
!N 		:dimension, assuming that all states are occupied
! It assumes that all are with spins in order: nspup of spin up, THEN nspdw
USE params, only: DP,eqstnspup, eqstnspdw
IMPLICIT NONE
COMPLEX(DP), INTENT(IN)                  	:: mat(N,N)
REAL(DP), INTENT(OUT)                   :: eigen(N)
COMPLEX(DP),INTENT(OUT)					:: Vect(N,N)
INTEGER ::N,nspup,nspdw,i,j
nspup=eqstnspup
nspdw=eqstnspdw
vect=0.d0
eigen=0.d0
WRITE (*,*) 'nspup, nspdw', nspup, nspdw
CALL CDIAG(mat(1:nspup,1:nspup),eigen(1:nspup),vect(1:nspup,1:nspup),nspup,nspup)   
CALL CDIAG(mat(nspup+1:N,nspup+1:N),eigen(nspup+1:N),vect(nspup+1:N,nspup+1:N),nspdw,nspdw)      
END SUBROUTINE cdiagspin

!_______________________________________________indexx___________________________________
SUBROUTINE indexx (n,arrin,indx)

!c     Use the Heapsort algorithm to index an array arrin of length n.
!c     Output the array indx such that arrin(indx(j)) is in ascending
!c     order for j = 1,2,...,n.  The input quantities n and arrin are
!c     not changed.

!c     This is a Numerical Recipes routine, but modified by one
!c     line to work if n equals 1.

      INTEGER i, n, indx(n), indxt, ir, j, l
      real*8 arrin(n), q
!WRITE(*,'(20f12.6)')(arrin(i),i=1,n)

      DO 11 j=1,n
        indx(j)=j
11    continue
      IF (n .eq. 1) return
      l=n/2+1
      ir=n
10    continue
        IF(l.gt.1)THEN
          l=l-1
          indxt=indx(l)
          q=arrin(indxt)
        ELSE
          indxt=indx(ir)
          q=arrin(indxt)
          indx(ir)=indx(1)
          ir=ir-1
          IF(ir.eq.1)THEN
            indx(1)=indxt
            return
          ENDIF
        ENDIF
        i=l
        j=l+l
20      IF(j.le.ir)THEN
          IF(j.lt.ir)THEN
            IF(arrin(indx(j)).lt.arrin(indx(j+1)))j=j+1
          ENDIF
          IF(q.lt.arrin(indx(j)))THEN
            indx(i)=indx(j)
            i=j
            j=j+j
          ELSE
            j=ir+1
          ENDIF
        go to 20
        ENDIF
        indx(i)=indxt
      go to 10
      END SUBROUTINE indexx
!__________________________________________OccupPerSpin_________________________________________________
SUBROUTINE occupPerSpin(mess,Occ)
use params,only: DP, nstate, ispin,occup
IMPLICIT NONE
 CHARACTER(*)::mess
REAL(DP)::Occ(2)
INTEGER::i
Occ=0.d0
DO i=1,nstate
 Occ(ispin(i))=Occ(ispin(i))+occup(i)
ENDDO
WRITE (*,'(a,a,2f12.8)')'OccupPerSpin',mess,Occ(1),Occ(2)
END SUBROUTINE occupPerSpin


!_____________________________CorrectEnergy2________________________________________________________
SUBROUTINE CorrectEnergy2(Wref,Eref,w,E,Wout,nloc)
use params,only: DP
REAL(DP):: Wref,Eref, W(nloc),E(nloc),Wout(nloc)
REAL(DP):: Wloc,Eloc,dE, dW,fi(Nloc),A0,A1,A2,D,lam,mu
INTEGER i
!compute Wloc, Eloc
Wloc=sum(W)
Eloc=dot_product(W,E)
!compute de,dW
dw=Wref-Wloc
dE=Eref-Eloc
A0=0.d0
A1=0.d0
A2=0.D0
fi=0.d0
DO i=1,nloc
 IF(w(i).eq.0.5d0) THEN
     fi(i)=1.d20
 ELSEIF(w(i)>1.d-20 .and. w(i)<(1-1.d-20) ) THEN
     fi(i)=1/( log(w(i))- Log(1-w(i)) )**2
 ENDIF
 A0=A0+fi(i)
 A1=A1+fi(i)*E(i)
 A2=A2+fi(i)*E(i)**2
ENDDO
D=A0*A2-A1*A1
lam=(+dE*A1-dW*A2)/D
mu=(dW*A1-dE*A0)/D
DO i=1,nloc
  Wout(i)=W(i)-lam*fi(i)-mu*fi(i)*E(i)
ENDDO
    WRITE (*,'(a,10f15.8)')'wr, er w, e', Wref, Eref, sum(Wout), dot_product(Wout, E)
    WRITE (*,'(a,12f15.8)')'E', E

END SUBROUTINE CorrectEnergy2


!_____________________________ordo_per_spin_________________________________
SUBROUTINE ordo_per_spin(psi)

USE params, only: DP,nspdw,nstate,ispin,occup,kdfull2,kstate,eqstnspup,eqstnspdw
IMPLICIT NONE

COMPLEX(DP), INTENT(IN OUT)    :: psi(kdfull2,kstate)

REAL(DP),ALLOCATABLE     :: occuploc(:)
COMPLEX(DP), ALLOCATABLE :: psiloc(:,:)
INTEGER,ALLOCATABLE      :: ispinloc(:)
INTEGER                  :: i,j,ind(kstate),compt,nspup

ALLOCATE(occuploc(kstate),ispinloc(kstate),psiloc(kdfull2,kstate))
compt=0
DO i=1,nstate
 
 IF (ispin(i).eq.1)THEN
    compt=compt+1
    ind(compt)=i
 ENDIF
 WRITE(*,*) i,compt, ind(compt)
ENDDO
eqstnspup=compt
!IF (.not.(compt.eq.nspup)) stop'inconsistent number of spins up'
DO i=1,nstate
 IF(ispin(i).eq.2)THEN
    compt=compt+1
    ind(compt)=i
 ENDIF
 WRITE(*,*) i,compt, ind(compt)
ENDDO
eqstnspdw=nstate-eqstnspup
!IF (.not.(compt.eq.nstate)) stop'inconsistent number of spins tot'
psiloc=psi
ispinloc=ispin
occuploc=occup
DO i=1,nstate
 psi(:,i)=psiloc(:,ind(i))
 ispin(i)=ispinloc(ind(i))
 occup(i)=occuploc(ind(i))
ENDDO
DEALLOCATE(occuploc,ispinloc,psiloc)
END SUBROUTINE ordo_per_spin
!___________________________________________Temperature_______________________________________________________________
SUBROUTINE temperature(mu,T)
!compute a sort atemperature by fitting occup
USE params,only: DP,occup,amoy,nstate
 IMPLICIT NONE
 double precision mu,T
 DOUBLE PRECISION TOL
 INTEGER INFO,lwa,i
 REAL(DP), ALLOCATABLE:: x(:),fvec(:),wa(:)
 INTEGER,ALLOCATABLE::iwa(:)
 INTEGER,parameter::n=2!number of independant variables(mu,t)
 EXTERNAL ff
 lwa= nstate*n+5*n+nstate
 ALLOCATE(iwa(nstate),wa(lwa),x(n),fvec(nstate))
       mu=(amoy(1)+amoy(nstate))/2.d0
       T=0.1d0
       x(1)=mu
       x(2)=T
       TOL=1.d-5
       WRITE(*,*)'bef lmdi', nstate,x(1),x(2),mu,T
       CALL lmdif1(ff,nstate,2,x,fvec,tol,info,iwa,wa,lwa)
       mu=x(1)
       T=x(2)
       DO i=1,nstate
         WRITE(*,*) 'amoy(i), occup(i)',amoy(i), occup(i),fvec(i)+occup(i)
       ENDDO
 DEALLOCATE(iwa,wa,x,fvec)
END SUBROUTINE temperature
!________________________________ff__________________________________________
!fcn(m,n,x,fvec,iflag)
SUBROUTINE ff(m,n,X,FVEC,IFLAG)
USE params,only: DP,occup,amoy
IMPLICIT NONE
INTEGER:: m,n,i,iflag
REAL(DP):: mu,T
REAL(DP) ::X(n),Fvec(m)
mu=X(1)
T=x(2)
WRITE(*,*) 'in ff m,n,mu, T',m,n,mu,T
DO i=1,m
  fvec(i)=1.d0/(1.d0+exp((amoy(i)-mu)/T))-occup(i)! fvec is a difference!!!
ENDDO

!WRITE(*,'(a,i3,25f9.6)') 'in ff: iflag, mu T fvec',iflag,mu,T,fvec

END SUBROUTINE ff
SUBROUTINE cproject(qin,qout,ispact,q0)

!     ******************************
!     projects all (!) states 'q0' out of 'qin'.
!      q0     = set of s.p. wavefunctions (complex)
!      qin    = wavefunction from which 'q0' are to be removed (complex)
!      qout   = resulting wavefunction
!      ispact = spin of 'qin'
!
!     old version! 

USE params
IMPLICIT NONE

COMPLEX(DP), INTENT(IN)   :: q0(kdfull2,kstate)
COMPLEX(DP), INTENT(IN)   :: qin(kdfull2)
COMPLEX(DP), INTENT(OUT)  :: qout(kdfull2)
INTEGER, INTENT(IN)    :: ispact

INTEGER :: nbe
COMPLEX(DP) :: ovl

!*********************************************************

qout=qin

DO nbe=1,nstate
  IF(ispin(nbe) == ispact) THEN
    ovl=dvol*SUM(CONJG(q0(:,nbe))*qout)
    qout(:)=qout(:)-q0(:,nbe)*ovl
  END IF
END DO

RETURN
END SUBROUTINE cproject
