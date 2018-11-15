 !MV
!___________________________________rta______________________________________________
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
INTEGER, INTENT(IN) :: iterat

COMPLEX(DP)                  :: psiorth(kdfull2,kstate),psieq(kdfull2,kstate),scal(nstate,nstate)
COMPLEX(DP)                  :: rhopsieq(nstate,nstate)
REAL(DP)                     :: occuppsi(kstate),occuporth(kstate),Estar,trace, sigee, rs,oneontaurelax,Nel0,Nel
REAL(DP)                     :: Occspin(2),Occref(2),eta,Eperspinref(2),Eperspintarget(2),entropy,Eref
REAL(DP)                     :: mu,T
integer::i,j,ii,nspup,nspdw
sigee=rtasigee!6.5d0 (rta2)
rs=rtars!3.7(rta2)
        nspup=eqstnspup
        nspdw=eqstnspdw
        CALL dyn_mfield(rho,aloc,psi,0D0) 
        write(*,*)'_____________before orth'
        call info(psi,rho,aloc,iterat)!to update parameters as energy per mode, etc
        occuppsi=occup!to store it
	call OccupPerSpin('start rta',Occref)!memorize for final correction
!        call calc_Eref(occup,ispin,amoy,Eperspinref)!memorize for final correction
	call srhomat(psi,aloc,psiorth,occuporth)! new basis is psi orth with occuporth
	occup=occuporth!because dyn_mfield needs occup
	Nel0=sum(occup(1:nstate))!number of electrons
        CALL dyn_mfield(rho,aloc,psiorth,0D0) 
        write(*,*)'_____________before eq'
        call info(psiorth,rho,aloc,iterat)
        Eref=energy
        write(*,'(a,40f10.6)')'occup    ',(occup(i),i=1,nstate)
	write(*,'(a,40f10.6)')'occuporth',(occuporth(i),i=1,nstate)
	write(*,'(a,40i5)')'ispin',(ispin(i),i=1,nstate)
	call OccupPerSpin('start rta',Occref)!memorize for final correction
	
	!compute the T distribution of occupation numbers
 call eqstate(psiorth,aloc,rho,psieq,occuporth)!find density current constained state of same energy(T increased)
 call OccupPerSpin('end eqstate',Occspin)

 call occupT0(occup(1:nstate),amoy(1:nstate),Estar)! computes Estar, occup is unchanged

write(*,*)'_____________________ after eqstate_______________________________'
        CALL dyn_mfield(rho,aloc,psieq,0D0) 
        CALL analyze_elect(psieq,rho,aloc,iterat)!update field and print
Nel=0.d0;do ii=1,nstate;Nel=Nel+occup(ii);enddo!sum(occup)! compute new number of electron, used in relaxation time calculation
!compose the dcmf state with the psiorth
oneontaurelax=0.4d0*sigee/rs**2*Estar/Nel!inverse of relaxation time
eta=dt1*oneontaurelax*jrtaint !eta for rta scheme
write(*,*)'dt1,Nel,oneontaurelax,jrtaint',dt1,Nel,oneontaurelax,jrtaint

!compose rhoorth and rhoeq
!spin up
 call calcrhoeq(psiorth(:,1:nspup),psieq(:,1:nspup),psi(:,1:nspup), occuporth(1:nspup)&
               ,occup(1:nspup),nspup)!compute the rho eq operator for next loop of iterations, modifies psiorth->psi, occup
               write(*,*)'just after rhoeq'
!spin down
 call calcrhoeq(psiorth(:,nspup+1:nstate),psieq(:,nspup+1:nstate),psi(:,nspup+1:nstate), occuporth(nspup+1:nstate)&
               ,occup(nspup+1:nstate),nspdw)!compute the rho eq operator for next loop of iterations, modifies psiorth->psi, occup
 call OccupPerSpin('end calcrhoeq',Occspin)
!   write(*,*)'electron loss',(sum(occup(1:nstate))-Nel0)/Nel0 !unnecessary, number of electrons adjusted in calcrhoeq
trace=0.d0
do i=1,nstate
 trace=trace+real(scal(i,i))
enddo
write(*,'(a,10f12.8)') 'Estar,trace, oneontaurelax Nel eta',Estar,trace,oneontaurelax,Nel,eta
!compute energies, etc
write(*,*)'_____________________energies after rhoeq_______________________________'
       CALL dyn_mfield(rho,aloc,psi,0D0) !update field to new density
       CALL analyze_elect(psi,rho,aloc,iterat)!to have energy
write(*,*) 'occup after rta',sum(occup(1:nstate))
         !call CorrectEnergy(energy,Eref)! one first call to correct energy
        call calc_Eref(occup,ispin,amoy,Eperspinref)
        Eperspintarget=Eperspinref+(Eref-energy)*Eperspinref/sum(Eperspinref)
        call CorrectEnergy2(occref(1),Eperspintarget(1),occup(1:nspup),amoy(1:nspup),occup(1:nspup),nspup)
        call CorrectEnergy2(occref(2),Eperspintarget(2),occup(nspup+1:nstate),amoy(nspup+1:nstate),occup(nspup+1:nstate),nspdw)

write(*,*)'_____________________energies after Correct_Energy________________________'
       CALL dyn_mfield(rho,aloc,psi,0D0) !update field to new density
       CALL analyze_elect(psi,rho,aloc,iterat)!to have energy
         !call CorrectEnergy(energy,Eref)! one second call to correct energy
        call calc_Eref(occup,ispin,amoy,Eperspinref)
        Eperspintarget=Eperspinref+(Eref-energy)*Eperspinref/sum(Eperspinref)
        call CorrectEnergy2(occref(1),Eperspintarget(1),occup(1:nspup),amoy(1:nspup),occup(1:nspup),nspup)
        call CorrectEnergy2(occref(2),Eperspintarget(2),occup(nspup+1:nstate),amoy(nspup+1:nstate),occup(nspup+1:nstate),nspdw)
        call OccupPerSpin('first correct nergy',Occspin)
write(*,*)'_____________________energies after second Correct_Energy________________________'
       CALL dyn_mfield(rho,aloc,psi,0D0) !update field to new density
       CALL analyze_elect(psi,rho,aloc,iterat)
       call OccupPerSpin('second correct nergy',Occspin)
!compute and store entropy
entropy=0.d0
do i=1,nstate
  if((occup(i)>1.d-15).and.(occup(i).lt.(1.d0-1.d-15)))then
   entropy=entropy+occup(i)*log(occup(i))+(1-occup(i))*log(1-occup(i))
 endif
 !compute temperature
enddo
 call temperature(mu,T)
      OPEN(1002,POSITION='append',FILE='prta')
      write(1002,'(11f14.7)') tfs,entropy, elasermask,mu,T
      close(1002)

!stop!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    contains
!_________________________subroutine calcrhoeq()_______________________________________________
subroutine calcrhoeq(psiorthloc,psieqloc,psiloc,occuporthloc,occuploc,nstateloc)
    !to make a sum of psiorth(occuporth)*(1-eta) with psieq(occup)*eta projected 
    !and diagonalize
    !this is done in a basis of 2*nstate
    !then conserve only the nstate most populated states
    USE util, ONLY:wfovlp
    implicit none

    integer::nstateloc,k,l
    complex(DP)::psiorthloc(kdfull2,nstateloc),psieqloc(kdfull2,nstateloc),psiloc(kdfull2,nstateloc)
    real(DP)::occuporthloc(nstateloc),occuploc(nstateloc),coeffoccup

    

    complex(DP)::scal0(2*nstateloc,2*nstateloc),scal1(2*nstateloc,2*nstateloc),tab(kdfull2,2*nstateloc)
    complex(DP)::taborth(kdfull2,2*nstateloc)
    complex(DP)::occloc(2*nstateloc,2*nstateloc)
    complex(DP)::uloc(2*nstateloc,2*nstateloc),vloc(2*nstateloc,2*nstateloc),vploc(2*nstateloc,2*nstateloc)
    complex(DP)::wloc(2*nstateloc,2*nstateloc),D1onrholoc(2*nstateloc,2*nstateloc)!for the diagonalisation process
    real(DP)::eigenloc(2*nstateloc)
    
    integer indx(2*nstateloc),i,j,ncall
 write(*,*) 'in calcrhoeq, nstateloc',nstateloc
    !allocate matrices
    !tab is the table of wave functions
    tab(:,1:nstateloc)=psiorthloc(:,1:nstateloc)
    tab(:,nstateloc+1:2*nstateloc)=psieqloc(:,1:nstateloc)
 	!occloc diagonal matrix of occupation numbers
	occloc=0.d0
	do i=1,nstateloc
	     occloc(i,i)=(1-eta)*occuporthloc(i)
	     occloc(nstateloc+i,nstateloc+i)=eta*occuploc(i)
	enddo
	!compute scal product <psi|psi>
	do j=1,2*nstateloc
	  do i=1,2*nstateloc
	      scal1(i,j)=wfovlp(tab(:,i),tab(:,j))!no spin involved, each call to rhoeq in a subspace of spin
	  enddo
	enddo
    !compute occup in the basis psiorth, psi
    occloc=matmul(transpose(conjg(scal1) ),occloc)
    occloc=matmul(occloc,scal1) 
 write(*,*) 'in calcrhoeq bef diag'
   !diagonalize scalar product 
    call cdiagmat(scal1, eigenloc, vloc, 2*nstateloc)
 write(*,*) 'in calcrhoeq aft diag'
   !compute diagonal matrices D1on_rho
    D1onrholoc=0.D0
    do i=1,2*nstateloc
		D1onrholoc(i,i)=1/sqrt(eigenloc(i))
    enddo
   !scale v into vp to have a unity matrix
    vploc=matmul(vloc,D1onrholoc)
   !compute Vp T *occloc*Vp
    scal0=matmul(transpose(conjg(Vploc) ),occloc)
    scal0=matmul(scal0,vploc )
   !diagonalize Vp T *rhomat*Vp: eigenvectors in u
 write(*,*) 'in calcrhoeq bef diag2'
    call cdiagmat(scal0,eigenloc,uloc,2*nstateloc)
   !compute the new change of basis matix
 write(*,*) 'in calcrhoeq aft diag2'
    wloc=matmul(vploc,uloc)
	!compute the phi_j in the r basis: wloc are the taborth in the tab basis
	!here phi is in table tab
	taborth=0.d0
	do j=1,2*nstateloc
	  do i=1,2*nstateloc
	  taborth(:,j)=taborth(:,j)+wloc(i,j)*tab(:,i)
	  enddo
	enddo
   
   coeffoccup=sum(eigenloc(1:2*nstateloc))
   !ordonnate the occupation numbers
    ncall=2*nstateloc
    call indexx(ncall,eigenloc,indx) 
 write(*,*) 'in calcrhoeq aft indxx'
   write(*,*)'eigensorted'
   write(*,'(10f14.10)') (eigenloc(indx(i)),i=2*nstateloc,1,-1)
!    write(*,*) 'ispin '
!    write(*,'(40i3)')(ispinloc(i),i=1,2*nstateloc)
!   write(*,*) 'ispin sorted'
!   write(*,'(40i3)')(ispinloc(indx(i)),i=2*nstateloc,nstateloc+1,-1)
   !select and store in psi the nstateloc modes with higher occupation numbers
   j=0
   do i=2*nstateloc,nstateloc+1,-1
          j=j+1
          psiloc(:,j)=taborth(:,indx(i))
          occuploc(j)=eigenloc(indx(i))
   enddo
   l=0
!   coeffoccup=coeffoccup/sum(occuploc(1:nstateloc))!compute a scaling factor for occupation numbers
!   occup(1:nstateloc)=occup(1:nstateloc)*coeffoccup!scale occupation numbers by this factor
!   write(*,*)'1-coeffoccup, occuptot'
!   write(*,'(12f14.10)')1.d0-coeffoccup,sum(occup(1:nstateloc))

end subroutine calcrhoeq

end subroutine rta

!________________________________________cdiagmat________________________________________________
subroutine cdiagmat(mat, eigen, vect, N)
!mat	:matrice to be diagonalised.
!eigen 	:vector of eigenvalues
!vect  	:eigenvectors written in columns
!N 		:dimension, assuming that all states are occupied
implicit none
INTEGER,PARAMETER :: DP=KIND(1D0)  ! precision  setting
COMPLEX(DP), INTENT(IN)                  	:: mat(N,N)
REAL(DP), INTENT(OUT)                   :: eigen(N)
COMPLEX(DP),INTENT(OUT)					:: Vect(N,N)
COMPLEX(DP)::mat0(N,N)
integer ::N,i,j
vect=0.d0
eigen=0.d0
mat0=mat
!call CDIAG(mat(1:nspup,1:nspup),eigen(1:nspup),vect(1:nspup,1:nspup),nspup,nspup)   
!call CDIAG(mat(nspup+1:N,nspup+1:N),eigen(nspup+1:N),vect(nspup+1:N,nspup+1:N),nspdw,nspdw)   
call HEigensystem(N,mat0,N,eigen,vect,N,0) 
vect=transpose(conjg(vect))
   
end subroutine cdiagmat


!_________________________________Calc-Eref__________________________________________
!compute the single particule energy per spin, for use in the fermi convergence
subroutine calc_Eref(occup,ispin,Ei,Eref)
use params,only:DP,nstate
implicit none
Real(DP), intent (in):: occup(nstate),Ei(nstate)
Real(DP), intent(out):: Eref(2)
integer, intent (in) ::ispin(nstate)
integer:: i
Eref=0.d0
 do i=1,nstate
    Eref(ispin(i))=Eref(ispin(i))+occup(i)*Ei(i)
 enddo
end subroutine calc_Eref
!_________________________________eqstate____________________________________________

subroutine eqstate(psi,aloc,rho,psi1,occuporth,iterat)
!update psi into psi1, and occuporth(entrance) into occup(in params)
USE params, only: DP, kdfull2,kstate,nstate,dvol,ispin,ifsicp,outnam,&
                  nyf,nxyf,centfx,nx2,ny,nz,enonlo,&
                  amoy,energy,rtamu,rtamuj,rtasumvar2max,&
                  occup,eqstnspup,eqstnspdw,rtatempinit
USE kinetic
implicit none
COMPLEX(DP), INTENT(IN OUT)  :: psi(kdfull2,kstate),psi1(kdfull2,kstate)
REAL(DP), INTENT(IN OUT)     :: aloc(2*kdfull2),rho(2*kdfull2),occuporth(kstate)
INTEGER, INTENT(IN) :: iterat

COMPLEX(DP)                  :: dpsi(kdfull2)

real(DP)::rhotot0(kdfull2,2),rhotot1(kdfull2,2),err,ma0,ma1,err1
real(DP)::lambda(kdfull2,2),mu,time0,time1
real(DP)::curr0(kdfull2,3),curr1(kdfull2,3),lambdaj(kdfull2,3),errj,muj!parameters for current
real(DP)::sumvar2,eal,fac,EspPerMod(nstate),EspTotRef
REAL(DP)::EspPerSpinRef(2),EspPerSpinAchieved(2),EspPerSpinTarget(2),temp(2),mut,EspTotAchieved
integer::i,j,ishift,ii,nspup,nspdw
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
 call calcrhotot(rhotot0,psi)
 rhotot1=rhotot0
 call calc_current(curr0,psi)
 curr1=curr0
 ma0=dvol*sum(rhotot0)
write(*,*) '  j      centfx       err           errj          mu        muj           &
            ma0           ma1'
j=0
!once for all compute aloc, without the local part
CALL coul_mfield(rho)


    call calc_Eref(occuporth,ispin,amoy,EspPerSpinRef)!compute a reference sp energy per spin,with occuporth
    EspTotRef=EspPerSpinRef(1)+EspPerSpinRef(2)
 write(*,*) 'appel initial fermi'
    call fermi1(amoy,EspPerSpinRef(1),occup,1,rtatempinit/10,1.0,temp(1),mut)!occup changed in Fermi
    call fermi1(amoy,EspPerSpinRef(2),occup,2,rtatempinit/10,1.0,temp(2),mut)!occup changed in Fermi
    EspPerSpinTarget=EspPerSpinRef!initial value of EspPerSpinTarget
    EspTotAchieved=0.d0!to avoid exit at first pass
!The big issue:
! when converging, some modes of low weght are crossing, nd getting low energy.
!at next fermi the density is changed a lot and err and sumvar are exploding. Also energy is ver out of balance
!in order to control this this version introduces an additionnalloop, with little authority (0.001 on Temp)
! and evry tenth step at most, in order to reorder the Occ vs Esp     
   do while(sumvar2>rtasumvar2max.or.abs(EspTotRef-EspTotAchieved)>1.d-4.or.err>0.2d0)!main loop
        if((sumvar2<1.d-1.and.err<1.0).and.mod(j,10).eq.0)then !update the occupation numbers without changing the target energy
                !idea: keep the occ in line with  the energies
                !at high temp and temp rate, the energy level are crossing
    		call fermi1(EspPerMod,EspPerSpinTarget(1),occup,1,max(temp(1)-0.001,rtatempinit/10.),temp(1)+0.001,temp(1),mut)!
    		call fermi1(EspPerMod,EspPerSpinTarget(2),occup,2,max(temp(2)-0.001,rtatempinit/10.),temp(2)+0.001,temp(2),mut)!
        endif

	if((sumvar2<1.d-2.and.err<.5).or.err<1.d-1) then!if there is enough convergence update fermi occupation numbers
 		EspPerSpinAchieved=0.d0
  		do ii=1,nstate
   			EspPerSpinAchieved(ispin(ii))=EspPerSpinAchieved(ispin(ii))+occup(ii)*EspPerMod(ii)!err on ener per sp tot per spin with ekmod computed after cschmidt (in calc psi1)
  		enddo
    		EspPerSpinTarget=EspPerSpinTarget-(EspPerSpinAchieved-EspPerSpinRef)/2.0d0
    		write(*,'(a,20f14.7)') 'EspPerSpinAchieved-EspPerSpinRef',EspPerSpinAchieved(1)-EspPerSpinRef(1)&
    		,EspPerSpinAchieved(2)-EspPerSpinRef(2)
    		write(*,'(a,20f14.7)') 'EspPerSpinTarget-EspPerSpinRef',EspPerSpinTarget(1)-EspPerSpinRef(1)&
    		,EspPerSpinTarget(2)-EspPerSpinRef(2)
   		write(*,*) 'fermi dans la boucle'
    		call fermi1(EspPerMod,EspPerSpinTarget(1),occup,1,max(temp(1)-0.005,rtatempinit/10.),temp(1)+0.005,temp(1),mut)!call with some phase advance, using error on ek occup changed in Fermi
    		call fermi1(EspPerMod,EspPerSpinTarget(2),occup,2,max(temp(2)-0.005,rtatempinit/10.),temp(2)+0.005,temp(2),mut)!call with some phase    	else
   		EspPerSpinTarget=EspPerSpinRef
   	endif
 j=j+1! counts the number of loops
!next value of psi1
  call calc_psi1(psi1,aloc,rhotot0,rhotot1,curr0,curr1,j,lambda,mu,lambdaj,muj,sumvar2,eal,EspPerMod)! one step of damped gradient
  EspTotAchieved=0.d0; do ii=1,nstate; EspTotAchieved= EspTotAchieved+&
               EspPerMod(ii)*occup(ii);enddo
  err=sum(abs(rhotot1-rhotot0))
  errj=sum(abs(curr1-curr0))
  ma1=dvol*sum(rhotot1)
  lambda=lambda+2.d0*(rhotot1-rhotot0)*mu!update lambda(Augmented Lagrangian rule)
  lambdaj=lambdaj+2.d0*(curr1-curr0)*muj
  if(mod(j,1).eq.0) then
      write(*,'(i6,20f14.7)') j,centfx,err,errj,mu,muj,ma0,ma1,EspTotRef,EspTotAchieved&
      ,EspTotRef-EspTotAchieved,sumvar2,EspPerSpinAchieved(1)-EspPerSpinAchieved(2)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  endif 
enddo!do j=1...
OPEN(1000,STATUS='old',POSITION='append',FORM='formatted',FILE='pspeed.'//outnam)
        write (1000,'(6g14.5)') (rhotot0(ny*(nyf+nxyf)+ii,1),&
          rhotot0(ny*(nyf+nxyf)+ii,2),rhotot1(ny*(nyf+nxyf)+ii,1),rhotot1(ny*(nyf+nxyf)+ii,2),&
          curr0(ny*(nyf+nxyf)+ii,1),curr1(ny*(nyf+nxyf)+ii,1),ii=1,nx2)
          write(1000,*)
 close(1000)
        CALL cpu_time(time1)
        write(*,*)'j,time', j, (time1-time0)
        CALL dyn_mfield(rho,aloc,psi1,0D0) 
       write(*,*)'_____________in eqstate 2'
        call info(psi1,rho,aloc,iterat)
      write(*   ,'(2i4,11f14.7)') iterat,j,sumvar2,err,errj,mu,muj,energy
      OPEN(1001,POSITION='append',FILE='peqstate')
      write(1001,'(2i6,11f14.7)') iterat,j,sumvar2,err,errj,mu,muj,energy
      close(1001)

end subroutine eqstate
!_________________________________OccupT0______________________________________________
Subroutine OccupT0(occloc,esploc,Estar)
!this routine computes occupation numbers at temperature=0
!for each spin, ordonate energies, putting occupation numbers to 1, until only a value between 0,1 is left for the last occupation number
! then computes Estar. The calling occupation numbers are not changed
USE params, only:DP,ispin,nstate
implicit none
integer::isp,i,n,indx(nstate),order(nstate)
real(DP)::occloc(nstate),esploc(nstate),Estar
Real(DP)::Etot,OccTot,EspSpin(nstate),OccSpin(nstate),EspSpinT0,occup
! first compute tables of energies per spi
!here there is no assumption thar spin are in order
EspSpinT0=0.d0!the ergy at T=0
Etot=dot_product(Occloc(1:nstate),Esploc(1:nstate))!cumulate enrgies of both spins
do isp=1,2! for each spin
   n=0
   do i=1,nstate
      if (ispin(i).eq.isp) then    
      n=n+1!count spin isp elements
        OccSpin(n)=occloc(i)!sp occupation per spin
        EspSpin(n)=EspLoc(i)!sp energy per spin
        indx(N)=i
      endif
   enddo
   OccTot=sum(Occspin(1:N))
   !now we have energies per spin in ener
   call indexx (n,EspSpin,order) !order contains the indices of energies in increasing order 
   !and a second loop to populate the occupation numbers
   do i=1,n!a loop for this acting spin
     if(OccTot>=1.d0) then
       occup=1.d0
       OccTot=OccTot-1.d0
       elseif(Occtot<1.d0.and.OccTot>=0) then
       occup=OccTot
       OccTot=OccTot-1.d0
       else
       occup=0.d0
       OccTot=OccTot-1.d0
     endif
     EspSpinT0=EspSpinT0+occup*EspLoc(indx(order(i)))
     write(*,'(a,4f14.8)') 'occup, Esp', occup,occloc(indx(order(i))),EspLoc(indx(order(i)))
   enddo
enddo
Estar=Etot-EspSpinT0!stop
write(*,'(a,3f14.7)') 'estar,Etot,EspSpinT0',estar,Etot,EspSpinT0
end Subroutine OccupT0
!_______________________________Calc_psi1____________________________________________
subroutine Calc_psi1(psi1,aloc,rhotot0,rhototloc,curr0,curr1,j,lambda,mu,lambdaj,muj,sumvar2,eal,ekmod)
USE params, only: DP, kdfull2,kstate,nstate,dvol,ispin,eye,ipsptyp,h2m,hbar,ame,pi,nx2,ny2,nz2,&
                  rtaeps,rtae0dmp,occup,enonlo,nxyz
USE kinetic
USE util, ONLY: wfovlp
implicit none
COMPLEX(DP), INTENT(IN OUT)                  :: psi1(kdfull2,kstate)
REAL(DP), INTENT(IN OUT)                     :: rhotot0(kdfull2,2),ekmod(nstate)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: rhototloc(kdfull2,2)
REAL(DP), INTENT(IN OUT)                     :: lambda(kdfull2,2),mu
real(DP), INTENT(IN OUT)                     ::curr0(kdfull2,3),curr1(kdfull2,3)
real(DP), INTENT(IN OUT)                     ::lambdaj(kdfull2,3),muj,eal!parameters for current

INTEGER, INTENT (IN)::j
COMPLEX(DP) :: psik(kdfull2),ALpsi(kdfull2),hpsi(kdfull2,kstate),psi1new(kdfull2,kstate)
COMPLEX(DP) :: vpsi(kdfull2),HpsiN(kdfull2),Apsi(kdfull2)
COMPLEX(DP) ::cnl,scal(nstate,nstate),hpsip(kdfull2,kstate)
Real(DP) ::epswf,e0dmp,deltalambda,emod(nstate),epot(nstate)
Real(DP) ::sumvar2,dkvol,eee,ekin
integer N,ii,ishift,ii1

dkvol=((2*pi)**3)*dvol!I do not understand this coefficient, it is used in calc nn, and works
 epswf=rtaeps
 !write (*,*)epswf
 e0dmp=rtae0dmp
 eal=0.d0
    !compute real space part of augmented hamiltonien,transfer in k space
    do N=1,nstate
!if(occup(N)>1.d-8) then! this is added since I had bad convergenece when non ocupied states where used. Should work as no states do not contribute to rho or j
   call calc_hamiltonien
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
      call fftback(psik,psi1(1,N))!here psinew just to do printings with homogeneous psi, hpsi
      eal=eal+emod(N)*occup(N)
!endif
    enddo!N
	call calc_var(hpsip,psi1,sumvar2 )!hppsip in k space
if ((mod(j,1).eq.0)) then
	write(*,'(a,(20f10.6))')'emod1 ener + Lagrangi',(emod(N),N=1,nstate)
	write(*,'(a,(21f10.6))')'ekmod1 ener per mode',(ekmod(N),N=1,nstate)
	write(*,'(a,(21f10.6))')'occup',(occup(N),N=1,nstate)
!	write(*,'(a,f14.9)')'sumvar2 penalty',sumvar2
endif                     
 call cschmidt(psi1)
 call calcrhotot(rhototloc,psi1)
 call calc_current(curr1,psi1)
 do N=1,nstate
  call calc_ekin(psi1(1,N),ekmod(N))
  ishift = (ispin(N)-1)*nxyz
!  WRITE(*,*) ' before CALC_EPOT',N,1,ishift
!  call calc_epot(psi1(1,N),aloc(ishift+1),epot(N),enonlo(N),N,1)
vpsi=0.d0           
if (ipsptyp==1) then!use local potential and non local only if there is a non local potential
             call  nonlocalc(psi1(1,N),vpsi,0)  !vpsi is the term Vnonloc *psi  
           endif
             !add local potential
             do ii=1,kdfull2
                vpsi(ii)=vpsi(ii)+aloc(ishift+ii) *psi1(ii,N)
             enddo
epot(N)=0.d0
do ii=1,kdfull2
 epot(N)=epot(N)+conjg(psi1(ii,N))*Vpsi(ii)*dvol
enddo
  ekmod(N)=ekmod(N)+epot(N)!ekmod is total sp energy
 enddo
contains 
!________________________________________Calc_hamiltonien______________________________________________________________________
 subroutine calc_hamiltonien
   COMPLEX(DP),allocatable::A(:),grad(:),gradApsi(:),jgradpsi(:)
   real(DP)::hbarm
   allocate(A(kdfull2),grad(kdfull2),gradApsi(kdfull2),jgradpsi(kdfull2))
      ishift=(ispin(N)-1)*kdfull2
      hbarm=hbar/ame
      vpsi=0.d0

     ! compute lambda*psi
           ALpsi=lambda(:,ispin(N))*psi1(:,N)!real space
     ! compute rho*psi
           ALpsi=ALpsi+mu*(rhototloc( :,ispin(N) )-rhotot0( :,ispin(N) ) )*psi1(:,N)!real space
     ! compute v*psi
           if (ipsptyp==1) then!use local potential and non local only if there is a non local potential
             call  nonlocalc(psi1(1,N),vpsi,0)  !vpsi is the term Vnonloc *psi  
           endif
            !aloc is a pure local potential 
            !add local potential
             do ii=1,kdfull2
                vpsi(ii)=vpsi(ii)+aloc(ishift+ii) *psi1(ii,N)
             enddo
 !           endif
     !compute (muj (j-j0)+lambda) nabla phi+nabla((muj (j-j0)+lambda) phi)
        A=eye*(2*muj* (curr1(:,1)-curr0(:,1))+lambdaj(:,1) )
        Apsi=A(:)*psi1(:,N)
        call xgradient_rspace(psi1(1,N),grad)
        call xgradient_rspace(Apsi,gradApsi)
        jgradpsi=         A(:)*grad(:)+gradApsi
        A=eye*(2*muj* (curr1(:,2)-curr0(:,2))+lambdaj(:,2) )
        Apsi=A(:)*psi1(:,N)
        call ygradient_rspace(psi1(1,N),grad)
        call ygradient_rspace(Apsi,gradApsi)
        jgradpsi=jgradpsi+A(:)*grad(:)+gradApsi
        A=eye*(2*muj* (curr1(:,3)-curr0(:,3))+lambdaj(:,3) )
        Apsi=A(:)*psi1(:,N)
        call zgradient_rspace(psi1(1,N),grad)
        call zgradient_rspace(Apsi,gradApsi)
        jgradpsi=jgradpsi+A(:)*grad(:)+gradApsi
        !total hamiltonian term
        ALpsi=vpsi+ALpsi-hbarm*jgradpsi!total potential+penalty term in real space
 deallocate(A,grad,gradApsi,jgradpsi)
 end subroutine calc_hamiltonien
end subroutine Calc_psi1
!____________________________________calcrhospin______________________________
SUBROUTINE calcrhotot(rho,q0)

!     density 'rho' for complex or real wavefunctions 'q0'
!experimental, only for non parallel scheme. To be rewritten
USE params
IMPLICIT none


COMPLEX(DP), INTENT(IN) :: q0(kdfull2,kstate)         ! cPW
REAL(DP), INTENT(OUT) :: rho(kdfull2,2)
integer nb,ishift,ind
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
!_________________________________________calc_var_____________________________________
subroutine calc_var(hpsi,psi1,sumvar2)
  use params, only :DP,kdfull2,nstate,kstate,ispin,occup
  USE kinetic
  USE util, ONLY: project
  implicit none
  COMPLEX(DP), INTENT(IN)                  :: psi1(kdfull2,kstate)
  COMPLEX(DP), INTENT(IN OUT)              :: hpsi(kdfull2,kstate)
  REAL(DP)                                 :: evarsp2,sumvar2
  Integer:: N
  COMPLEX(DP)                              ::wfovlp,qout(kdfull2)
  ! hpsi is the term h psi in k space
  ! psi1 the wavefunctions (in real space)
!                write(*,*)'in var'
sumvar2=0.d0
!pass hpsi in real space
do N=1,nstate
    call fftback(hpsi(1,N),hpsi(1,N))
enddo

!do the calc
do N=1,nstate
    CALL c_project(hpsi(:,N),qout,ispin(N),psi1)
    evarsp2 =  SQRT( real( wfovlp(qout,qout) ) )
    sumvar2   = sumvar2 + occup(N)*evarsp2**2
enddo
!write(*,'(a, 40i3)') 'ispin in sumvar2', (ispin(N),N=1,nstate)
!write(*,'(a, 40f10.6)') 'occup in sumvar2', (occup(N),N=1,nstate)

end subroutine calc_var
!_______________________________Fermi1______________________________________________________
subroutine fermi1(ekmod,eref,occup,ispinact,T0i,T1i,t2,mu)
!changes occup to achieve the total sp energy of eref, and the right total occupation
use params, only: DP, kstate,nstate,eqstnspup,eqstnspdw,rtatempinit
implicit none
   Real(DP)::ekmod(kstate),occup(kstate),Eref,occref(2),T(2)
   real(DP)::Et,mu,T0i,T1i,T0,T1,t2,occtot(2)
   integer ispinact,nspup,nspdw
   nspup=eqstnspup
   nspdw=eqstnspdw
   occref(1)=sum(occup(1:nspup))!total occupation spin1
   occref(2)=sum(occup(nspup+1:nstate))
   !T0, T1 are low and high boudaries of temp
   T0=T0i
   T1=T1i
   !dichotomy on T
   do while(abs(T1-T0)>1.d-10)
       T2=(T1+T0)/2.d0
       if (ispinact.eq.1) then
           call occT1(occref(1),ekmod(1:nspup),et,mu,occtot(1),nspup,T2,occup(1:nspup))! change mu at constatnt T2 to have the required total occupation
         else
           call occT1(occref(2),ekmod(nspup+1:nstate),et,mu,occtot(2),nspdw,T2,occup(nspup+1:nstate))
       endif
       if(et>eref) then !et is the energy at T2, with the mu giving the right occupation number
          T1=T2
       else 
          T0=T2
       endif
   enddo
   write(*,*) 'fermi1, T mu',T2,mu
end subroutine fermi1
!_______________________________________forceTemp____________________________________
subroutine forceTemp(amoy,occup,n,temp,mu)
use params,only:DP
implicit none
integer n
real(DP)::amoy(n),occup(n),temp,mu,occtot,Etot,occtotloc
occtot=sum(occup(1:n))
Etot=dot_product(occup(1:n),amoy(1:n))
write(*,*) 'occtot,Etot',occtot,etot
call OccT1(occtot,amoy,Etot,mu,occtotloc,n,temp,occup)
write(*,*)'force:T,mu,occtot,occtotloc',temp,mu,occtot,occtotloc
end subroutine forceTemp
!_________________________________OccT1______________________________________________
subroutine OccT1(occrefloc,enerloc,Etotloc,muloc,occtotloc,n,T,occuploc)
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
implicit none
   integer::n,i
   real(DP)::occrefloc,enerloc(n),Etotloc,muloc,T,occuploc(n)
   real(DP)::muloc0,muloc1,occtotloc,fac
   integer::orderloc(n)
   muloc0=minval(enerloc)
   muloc1=maxval(enerloc)
   do while(abs(muloc0-muloc1)>1.d-12)!dichotomy on mu to have the right ocupation number at T2
      muloc=(muloc0+muloc1)/2.d0
      occtotloc=0.d0
      Etotloc=0.d0
      do i=1,n
         occuploc(i)=1.d0/(1.d0+exp((enerloc(i)-muloc)/T))
         occtotloc=occtotloc+occuploc(i)
         Etotloc=Etotloc+occuploc(i)*enerloc(i)!Etotloc is the corresponding energy
      enddo
      if(occtotloc<occrefloc) then
         muloc0=muloc
      else
         muloc1=muloc
      endif
   enddo    
      
end subroutine occT1

!________________________________Fermi_init_________________________________________________   
subroutine fermi_init(ekmod,T,occup,ispinact)
use params, only: DP, kdfull2,kstate,nstate,ispin
implicit none
       REAL(DP), INTENT(in) :: ekmod(nstate),T
       REAL(DP), INTENT(in out) :: occup(kstate)
       INTEGER, intent (in)::ispinact
       REAL(DP)::mu1,mu2,mu3,o1,o2,o3,O,Occ
       integer::i
       O=0.d0
       do i=1,nstate
          if (ispin(i).eq.ispinact) O=O+occup(i)
       enddo
       
       mu1=-2
       mu2=0
       o1=occ(mu1,ekmod,T,occup,O,ispinact)
       o2=occ(mu2,ekmod,T,occup,O,ispinact)
       do while(abs(mu1-mu2)>1.d-5)
          mu3=(mu1+mu2)/2.d0
          o3=occ(mu3,ekmod,T,occup,O,ispinact)
          if(o3<0)then
             mu1=mu3
             o1=o3
          else  
             mu2=mu3
             o2=o3
          endif
       enddo
       do i=1,nstate
         if(ispin(i)==ispinact) occup(i)=1.d0/(1.d0+exp((ekmod(i)-mu3)/T))
       enddo
  end subroutine fermi_init
       real function occ(mu,ekmod,T,occup,O,ispinact)
       use params, only: DP,nstate,ispin,kstate
      implicit none
       REAL(DP), INTENT(in) :: ekmod(nstate),O
       REAL(DP), INTENT(in ) :: occup(kstate)
       real(DP)::mu,T
       INTEGER, intent (in)::ispinact
       REAL(DP)::mu1,mu2,mu3,o1,o2,o3
       integer::i
       OCC=0.d0
       do i=1,nstate
          if(ispin(i)==ispinact) Occ=Occ+1.d0/(1.d0+exp((ekmod(i)-mu)/T))
       enddo
       Occ=Occ-O
       end function Occ


!__________________________________________srhomat____________________________________________________________________
subroutine srhomat(psi,aloc,psiorth,occuporth)
USE params, only: DP, kdfull2,kstate,nstate,ispin,nrel2abs,nxyz,occup,apnum,psitophi
implicit none
REAL(DP), INTENT(IN)                     :: aloc(2*kdfull2)
COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate),psiorth(kdfull2,kstate)
real(DP),intent(out)				:: occuporth(kstate)
integer:: i,j
complex(DP):: scal(nstate,nstate),scal0(nstate,nstate)
 


      complex(DP)::D1on_rho(nstate,nstate),rhomat(nstate,nstate)
      REAL(8) :: eigen(nstate) ,eigen0(nstate)     ! eigenvalues
      COMPLEX(8) :: v(nstate,nstate),vp(nstate,nstate),u(nstate,nstate),w(nstate,nstate)  ! eigenvectors
      COMPLEX(8) :: winv(nstate,nstate)!inverse of w
      REAL(8)    :: occup1(nstate)





!Calculate scalar product of the <psi_i |psi_j>
!not hermitian
 call scalar(psi,psi,scal,ispin, 'scal psi(i) psi(j)')
!compute the density matrix in the basis psi
   rhomat=0.d0
   do i=1,nstate
       rhomat(i,i)=occup(i)
   enddo!i
scal0=matmul(transpose(conjg(scal) ),rhomat)
rhomat=matmul(scal0,scal )

!diagonalize scal
 call cdiagmat(scal, eigen0, v, nstate)
 write(*,'(a,30 f12.7)')'eigen0',(eigen0(i),i=1,nstate)
!compute diagonal matrices D1on_rho
D1on_rho=0.D0
do i=1,nstate
    D1on_rho(i,i)=1/sqrt(eigen0(i))
enddo
!scale v into vp to have a unity matrix
vp=matmul(v,D1on_rho)


!compute Vp T *rhomat*Vp
scal0=matmul(transpose(conjg(Vp) ),rhomat)
scal0=matmul(scal0,vp )



!diagonalize Vp T *rhomat*Vp: eigenvectors in u
 call cdiagmat(scal0,eigen,u,nstate)
!compute the new change of basis matix
w=matmul(vp,u)
winv=matmul (transpose ( conjg(u) ) ,transpose(conjg(v)))

!compute the phi_j in the r basis: Vect are the phi_j in the psi basis
psiorth=0.d0
do j=1,nstate
  do i=1,nstate
  psiorth(:,j)=psiorth(:,j)+w(i,j)*psi(:,i)
  enddo
enddo
 write(*,'(a,30 f12.7)')'eigen',(eigen(i),i=1,nstate)
!uptage the change of basis matrix
psitophi=matmul(winv,psitophi)!matmul(psitophi,winv)

occuporth=eigen

end subroutine srhomat

!______________________________scalar product________________________________________
subroutine scalar(tab1,tab2,scal,ispin, mess)
USE params, only: DP, kdfull2,nstate,kstate,nrel2abs
implicit none
 COMPLEX(DP), INTENT(IN)       :: tab1(kdfull2,kstate),tab2(kdfull2,kstate)
 COMPLEX(DP), INTENT(OUT)      :: scal(nstate,nstate)

 USE util, ONLY: wfovlp

 integer::ispin(nstate)
 CHARACTER(*)::mess
integer:: i,j
do j=1,nstate
  do i=1, nstate
    if( ispin(nrel2abs(i)) == ispin(nrel2abs(j) )) then
      scal(i,j)=wfovlp(tab1(:,i),tab2(:,j))
    else
      scal(i,j)=0.d0
    endif    
  enddo
enddo
write(*,*) mess
do i=1,nstate
    write(*, '(20f12.9)')(scal(i,j),j=1,nstate)
enddo
end subroutine scalar

!________________________________________cdiagspin________________________________________________
subroutine cdiagspin(mat, eigen, vect, N)
!mat	:matrice to be diagonalised.
!eigen 	:vector of eigenvalues
!vect  	:eigenvectors written in columns
!N 		:dimension, assuming that all states are occupied
! It assumes that all are with spins in order: nspup of spin up, then nspdw
USE params, only: DP,eqstnspup, eqstnspdw
implicit none
COMPLEX(DP), INTENT(IN)                  	:: mat(N,N)
REAL(DP), INTENT(OUT)                   :: eigen(N)
COMPLEX(DP),INTENT(OUT)					:: Vect(N,N)
integer ::N,nspup,nspdw,i,j
nspup=eqstnspup
nspdw=eqstnspdw
vect=0.d0
eigen=0.d0
write (*,*) 'nspup, nspdw', nspup, nspdw
call CDIAG(mat(1:nspup,1:nspup),eigen(1:nspup),vect(1:nspup,1:nspup),nspup,nspup)   
call CDIAG(mat(nspup+1:N,nspup+1:N),eigen(nspup+1:N),vect(nspup+1:N,nspup+1:N),nspdw,nspdw)      
end subroutine cdiagspin

!_______________________________________________indexx___________________________________
subroutine indexx (n,arrin,indx)

!c     Use the Heapsort algorithm to index an array arrin of length n.
!c     Output the array indx such that arrin(indx(j)) is in ascending
!c     order for j = 1,2,...,n.  The input quantities n and arrin are
!c     not changed.

!c     This is a Numerical Recipes routine, but modified by one
!c     line to work if n equals 1.

      integer i, n, indx(n), indxt, ir, j, l
      real*8 arrin(n), q
!write(*,'(20f12.6)')(arrin(i),i=1,n)

      do 11 j=1,n
        indx(j)=j
11    continue
      if (n .eq. 1) return
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          indxt=indx(l)
          q=arrin(indxt)
        else
          indxt=indx(ir)
          q=arrin(indxt)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir.eq.1)then
            indx(1)=indxt
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(arrin(indx(j)).lt.arrin(indx(j+1)))j=j+1
          endif
          if(q.lt.arrin(indx(j)))then
            indx(i)=indx(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        go to 20
        endif
        indx(i)=indxt
      go to 10
      end subroutine indexx
!__________________________________________OccupPerSpin_________________________________________________
subroutine occupPerSpin(mess,Occ)
use params,only: DP, nstate, ispin,occup
implicit none
 CHARACTER(*)::mess
real(DP)::Occ(2)
integer::i
Occ=0.d0
do i=1,nstate
 Occ(ispin(i))=Occ(ispin(i))+occup(i)
enddo
Write (*,'(a,a,2f12.8)')'OccupPerSpin',mess,Occ(1),Occ(2)
end subroutine occupPerSpin


!_____________________________CorrectEnergy2________________________________________________________
Subroutine CorrectEnergy2(Wref,Eref,w,E,Wout,nloc)
use params,only: DP
Real(DP):: Wref,Eref, W(nloc),E(nloc),Wout(nloc)
Real(DP):: Wloc,Eloc,dE, dW,fi(Nloc),A0,A1,A2,D,lam,mu
integer i
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
do i=1,nloc
 if (w(i).eq.0.5d0) then
     fi(i)=1.d20
 elseif (w(i)>1.d-20 .and. w(i)<(1-1.d-20) ) then
     fi(i)=1/( log(w(i))- Log(1-w(i)) )**2
 endif
 A0=A0+fi(i)
 A1=A1+fi(i)*E(i)
 A2=A2+fi(i)*E(i)**2
enddo
D=A0*A2-A1*A1
lam=(+dE*A1-dW*A2)/D
mu=(dW*A1-dE*A0)/D
do i=1,nloc
  Wout(i)=W(i)-lam*fi(i)-mu*fi(i)*E(i)
enddo
    write (*,'(a,10f15.8)')'wr, er w, e', Wref, Eref, sum(Wout), dot_product(Wout, E)
    write (*,'(a,12f15.8)')'E', E

end subroutine CorrectEnergy2


!_____________________________ordo_per_spin_________________________________________________________
subroutine ordo_per_spin(psi)
USE params, only: DP,nspdw,nstate,ispin,occup,kdfull2,kstate,eqstnspup,eqstnspdw
implicit none
COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)
REAL(DP),ALLOCATABLE :: occuploc(:)
COMPLEX(DP), allocatable                  :: psiloc(:,:)
integer,allocatable:: ispinloc(:)
integer i,j,ind(kstate),compt,nspup
allocate(occuploc(kstate),ispinloc(kstate),psiloc(kdfull2,kstate))
compt=0
do i=1,nstate
 
 if (ispin(i).eq.1)then
    compt=compt+1
    ind(compt)=i
 endif
 write(*,*) i,compt, ind(compt)
enddo
eqstnspup=compt
!if (.not.(compt.eq.nspup)) stop'inconsistent number of spins up'
do i=1,nstate
 if (ispin(i).eq.2)then
    compt=compt+1
    ind(compt)=i
 endif
 write(*,*) i,compt, ind(compt)
enddo
eqstnspdw=nstate-eqstnspup
!if (.not.(compt.eq.nstate)) stop'inconsistent number of spins tot'
psiloc=psi
ispinloc=ispin
occuploc=occup
do i=1,nstate
 psi(:,i)=psiloc(:,ind(i))
 ispin(i)=ispinloc(ind(i))
 occup(i)=occuploc(ind(i))
enddo
deallocate(occuploc,ispinloc,psiloc)
end subroutine ordo_per_spin
!___________________________________________Temperature_______________________________________________________________
subroutine temperature(mu,T)
!compute a sort atemperature by fitting occup
USE params,only: occup,amoy,nstate
 implicit none
 double precision mu,T
 DOUBLE PRECISION TOL
 integer INFO,lwa,i
 Real(8), allocatable:: x(:),fvec(:),wa(:)
 Integer,allocatable::iwa(:)
 integer,parameter::n=2!number of independant variables(mu,t)
 EXTERNAL ff
 lwa= nstate*n+5*n+nstate
 allocate(iwa(nstate),wa(lwa),x(n),fvec(nstate))
       mu=(amoy(1)+amoy(nstate))/2.d0
       T=0.1d0
       x(1)=mu
       x(2)=T
       TOL=1.d-5
       write(*,*)'bef lmdi', nstate,x(1),x(2),mu,T
       call lmdif1(ff,nstate,2,x,fvec,tol,info,iwa,wa,lwa)
       mu=x(1)
       T=x(2)
       do i=1,nstate
         write(*,*) 'amoy(i), occup(i)',amoy(i), occup(i),fvec(i)+occup(i)
       enddo
 deallocate(iwa,wa,x,fvec)
end subroutine temperature
!________________________________ff___________________________________________
!fcn(m,n,x,fvec,iflag)
subroutine ff(m,n,X,FVEC,IFLAG)
USE params,only: occup,amoy
implicit none
integer:: m,n,i,iflag
real(8):: mu,T
real(8) ::X(n),Fvec(m)
mu=X(1)
T=x(2)
write(*,*) 'in ff m,n,mu, T',m,n,mu,T
do i=1,m
  fvec(i)=1.d0/(1.d0+exp((amoy(i)-mu)/T))-occup(i)! fvec is a difference!!!
enddo

!write(*,'(a,i3,25f9.6)') 'in ff: iflag, mu T fvec',iflag,mu,T,fvec

end subroutine ff
