!-----afterburn---------------------------------------------------------

SUBROUTINE afterburn(psir,rho,aloc)

! dynamic propagation for electrons and ions

USE params
!USE kinetic
!USE util, ONLY: stimer,timer,safeopen,testcurrent,rhointxy,rhointyz,rhointxz
#if(twostsic)
USE twost, ONLY:tnearest,init_fsic,init_vecs,end_fsic,expdabvol_rotate_init
#endif
IMPLICIT NONE

REAL(DP), INTENT(IN OUT)                  :: psir(kdfull2,kstate)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)

COMPLEX(DP),ALLOCATABLE :: psiw(:,:),psi(:,:)

INTEGER  :: it,ion
!REAL(DP):: totalprob,totalovlp
!REAL(DP), EXTERNAL:: energ_ions   ! declared in ion_md
!REAL(DP), EXTERNAL:: enerkin_ions ! declared in ion_md



!       ****************************************************
!          imaginary-time iteration (to improve statics)
!       ****************************************************

 WRITE(*,*) ' start afterburn. isitmax=',isitmax
 CALL flush(6)
 ALLOCATE(psi(kdfull2,kstate))
#if(twostsic)
 IF(ifsicp >= 7 .AND. nclust > 0 )  CALL init_fsic()
 IF(ifsicp >= 8) CALL expdabvol_rotate_init
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
     CALL tstep_exp(psi,aloc,rho,it,psiw,.TRUE.)
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

RETURN

END SUBROUTINE afterburn
