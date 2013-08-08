#include"define.h"                                                                      

!-----init_occ_target--------------------------------------------
SUBROUTINE init_occ_target()
USE params

REAL(DP) :: xxdum,delta_etrgt
INTEGER  :: iidum
!------------------------------------------------------------------------
! Reads the file 'occ_eps_target' that should be in the working directory.
!------------------------------------------------------------------------

WRITE(6,*) 'opening occ_eps_target file'
OPEN(UNIT=91,STATUS='unknown',FORM='formatted', FILE='occ_spe_target')
!IF(istat == 0.AND.irest == 0) THEN
!   READ(91,*) xxdum,xxdum
!ELSE IF(istat == 1.AND.irest == 0) THEN
!IF(istat == 1.AND.irest == 0) THEN
!   READ(91,*) binerg,xxdum
!ELSE IF(istat == 0.AND.irest > 0) THEN
!ELSE IF(istat == 1.AND.irest > 0) THEN
   READ(91,*) binerg,reference_energy
!ELSE
!   STOP ' invalid option in ATTACH_PROB'
!END IF
WRITE(*,*) ' istat,irest,binerg,reference_energy:',istat,irest,binerg,reference_energy

!      aver_estar  = etot - binerg                                                            
aver_estar  = reference_energy - binerg
delta_etrgt = 3.d0*h2m/4.d0/scatterelectronw**2
emin_target = aver_estar - delta_etrgt
emax_target = aver_estar + delta_etrgt
WRITE(6,*)'aver_estar,delta_etrgt,emin_target,emax_target',  &
     aver_estar,delta_etrgt,emin_target,emax_target

DO i=1,10000
   READ(91,*,END=899) iidum,iidum,xxdum,xxdum,xxdum,xxdum
   nstate_target=nstate_target+1
END DO
899 CONTINUE

ALLOCATE(ispin_temp(nstate_target))
ALLOCATE(occ_target(nstate_target))
ALLOCATE(spe_target(nstate_target))

REWIND(91)
READ(91,*) xxdum,xxdum
DO i=1,nstate_target
   READ(91,*)iidum,ispin_temp(i),occ_target(i),xxdum,spe_target(i),xxdum
END DO
CLOSE(91)

WRITE(*,*) 'ispin_target:',ispin_temp(1:nstate_target)
WRITE(*,*) 'occ_target:',occ_target(1:nstate_target)
WRITE(*,*) 'nclust,nstate,nstate_target:',nclust,nstate,nstate_target

RETURN
END SUBROUTINE init_occ_target

!-----init_psitarget------------------------------------------------                        

SUBROUTINE init_psitarget()
USE params

INTEGER    :: idum,nstatedum,nclust_target,nion_target,nspdw_target
REAL(DP)   :: xxdum

!-------------------------------------------------------------------
! reading of w.f. from wfs_target
!-------------------------------------------------------------------
WRITE(*,*) 'enter init_psitarget'
OPEN(UNIT=90,STATUS='unknown',FORM='unformatted',FILE='../wfs_target')

READ(90) idum,nstatedum,nclust_target,nion_target,nspdw_target
write(6,*) 'target: nstate_t,nclust_t,nion_t,nspdw_t'                 
write(6,*) nstatedum,nclust_target,nion_target,nspdw_target
write(6,*) 'actual: nstate,nclust,nion,nspdw,kstate'
write(6,*) nstate,nclust,nion,nspdw,kstate

IF (nstatedum.ne.nstate_target) STOP 'number of w.f. in wfs_target not consistent with occ_spe_target file'
IF((nion_target /= nion).OR.(nspdw_target /= nspdw)) THEN
   WRITE(6,*)'attachment calculation impossible:'
   WRITE(6,*)'mismatch betw # of ions or spins'
   STOP ' ATTACHEMENT: mismatch betw # of ions or spins'
END IF

IF(nclust > 0)THEN
   DO nb=1,nstate_target
      READ(90) xxdum,psitemp(1:nxyz,nb)
!      WRITE(*,*) 'nb,xxdum',nb,xxdum
   END DO
END IF

678 CONTINUE
CLOSE(90)
write(6,*)'target wfs read'

RETURN
END SUBROUTINE init_psitarget

!-----attach_prob------------------------------------------------                        

SUBROUTINE attach_prob(totalprob,psi)
USE params

REAL(DP), INTENT(OUT)           :: totalprob
COMPLEX(DP), INTENT(IN)         :: psi(kdfull2,kstate)

COMPLEX(DP),ALLOCATABLE :: psitarget(:,:)

COMPLEX(DP) :: overlaps(kstate,kstate),submatr(kstate,kstate)
COMPLEX(DP) :: determinant,tbelement,det

COMPLEX(DP) :: psip(kdfull2),psipp(kdfull2)

COMPLEX(DP) :: wfovlp
INTEGER :: indx(nstate)             ! index field for LU decomp.                          
INTEGER :: index(nstate-2)          ! index field for LU decomp.                          
COMPLEX(DP) :: temp1,temp2
REAL(DP) :: xxdum,d
INTEGER :: iidum,nexpand,ierror
LOGICAL,PARAMETER :: ttest=.true.

!--------------------------------------------------------------                           

vcoll=1D0      ! strength of collisional pot.                                             

!WRITE(*,*) 'ispin(nstate)',ispin(nstate)
WRITE(*,*) 'enter attach_prob'

totalprob=0D0
nmatchenergy=0

ALLOCATE(ispin_target(nstate_target))
ispin_target(1:nstate_target)=ispin_temp(1:nstate_target)

!WRITE(*,'(200i3)') ispin_target(1:nstate_target)
!WRITE(*,'(200i3)') ispin(1:nstate)

ALLOCATE(psitarget(kdfull2,nstate))
DO ih=1,nstate-1
   IF(ispin_target(ih) == ispin(nstate)) THEN
!      DO ip1=nstate,nstate_target-1
      DO ip1=1,nstate_target-1
         DO ip2=ip1+1,nstate_target
!            WRITE(*,'(3i5,3(1pg13.5))') ih,ip1,ip2, &
!                  occ_target(ih),occ_target(ip1),occ_target(ip1)
            IF(ispin_target(ip1) == ispin(nstate).AND.ispin_target(ip2) == ispin(nstate)) THEN
               IF((occ_target(ih) > 0.5).AND.(occ_target(ip1) < 0.5).AND.(occ_target(ip2) < 0.5)) THEN
                  delta_e = spe_target(ip2)+spe_target(ip1)-spe_target(ih)
                  IF(ttest) WRITE(*,'(a,3i5,4(1pg13.5))') 'ih etc:',ih,ip1,ip2,&
                    spe_target(ih),spe_target(ip1),spe_target(ip2),delta_e
                  IF(delta_e > emin_target.AND.delta_e < emax_target) THEN
!                     write(6,*)'ih,ip1,ip2,delta_e',ih,ip1,ip2,delta_e
                     nmatchenergy=nmatchenergy+1
                     IF(ttest) WRITE(*,*) 'nmatch:',nmatchenergy

                     DO nb=1,nstate-1
                        psitarget(1:nxyz,nb)=psitemp(1:nxyz,nb)
                     END DO
                     psitarget(1:nxyz,ih)=psitemp(1:nxyz,ip1)
                     ispin_target(ih)=ispin_target(ip1)
                     psitarget(1:nxyz,nstate)=psitemp(1:nxyz,ip2)
                     ispin_target(nstate)=ispin_target(ip2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             
!     build matrix overlaps of s.p. states                                            

                     DO i=1,kstate
                        DO j=1,kstate
                           overlaps(i,j)=CMPLX(0.d0,0.d0)
                           submatr(i,j)=CMPLX(0.d0,0.d0)
                        END DO
                     END DO
                     
                     nexpand = nstate

                     IF(ttest) write(6,*) &
                        'entering overlap comp. loop with nexpand=', nexpand
                     DO i=1,nexpand
                        DO j=1,nexpand
                           IF(ispin_target(i) == ispin(j)) THEN
                              overlaps(i,j) = wfovlp(psitarget(1,i),psi(1,j))
                           ELSE
                              overlaps(i,j) = CMPLX(0.d0,0.d0)
                           END IF
                        END DO
                     END DO

                     IF(ttest) WRITE(*,*) 'overlap(5,:):',overlaps(5,:)

!                     WRITE(*,*) 'LU decomposition of overlaps matrix'
!                     CALL cludcmp_d(overlaps,nexpand,kstate,indx,d,det,ierror)
                     
                     CALL cludcmp_d(overlaps,nexpand,indx,d,det,ierror)

                     IF(ierror == 99) det = CMPLX(0D0,0D0)
!                     WRITE(6,'(f12.5,4i5,3(1pg13.5))') tfs,nmatchenergy,ih,ip1,ip2,&
!                          delta_e,REAL(det),imag(det)
!                     WRITE(811,'(f12.5,4i8,3(1pg13.5))') tfs,nmatchenergy,ih,ip1,ip2,&
!                          delta_e,REAL(det),imag(det)

!     accumulate total transition matrix element                                     

                     tbelement = CMPLX(0D0,0D0)
                     DO i1=1,nexpand
                        DO i2=i1+1,nexpand
                           IF(ttest) WRITE(*,'(a,4i3)') &
                              ' submatrix outer loops: i1,i2=',i1,i2
                           DO j1=1,nexpand
                              DO j2=j1+1,nexpand
                                 i1nn = i1
                                 i2nn = i2
                                 j1nn = j1
                                 j2nn = j2
                                 IF((ispin_target(i1)+ispin_target(i2))  &
                                      == ispin(j1)+ispin(j2)) THEN
!                                     extract submatrix                                  
                        
                                    ishift = 0
                                    DO i=1,nexpand
                                       IF(i == i1) ishift = 1+ishift
                                       IF(i == i2) ishift = 1+ishift
                                       jshift = 0
                                       DO j=1,nexpand
                                          IF(j == j1) jshift = 1+jshift
                                          IF(j == j2) jshift = 1+jshift
                                          IF(i /= i1 .AND. i /= i2 .AND. j /= j1 .AND. j /= j2) THEN
                                             submatr(i-ishift,j-jshift) = overlaps(i,j)
                                          END IF
                                       END DO
                                    END DO
                                    CALL cludcmp_d(submatr,nexpand-2,indx,d,det,ierror)

                                    IF(ierror == 99) det = CMPLX(0D0,0D0)
                                    IF(ierror == 0) THEN
                                       DO ind=1,kdfull2
                                          temp1=psitarget(ind,i1nn)*psitarget(ind,i2nn)
                                          temp2=psi(ind,j1nn)*psi(ind,j2nn)
                                          tbelement = CONJG(temp1)*temp2  + tbelement
                                       END DO
                                       tbelement=tbelement*det
!                                       tbelement = tbelement + wfovlp(psip,psipp)*det
                                    END IF
                                 END IF
                              END DO
                           END DO
                        END DO
                     END DO
!     final composition                                                                   
                    tbelement=tbelement*dvol
                    totalprob = totalprob + ABS(tbelement*vcoll)**2

                  END IF    ! test on 2p1h energies                                       
               END IF    ! test on 2p1h occupations                                       
            END IF    ! test on ip1 and ip2 spins                                         
         END DO      ! loop on ip2                                                        
      END DO      ! loop on ip1                                                           
   END IF        ! test on ih spin                                                        
END DO! loop on ih   
                                                               
WRITE(*,'(a,200i3)') 'ispin_target:',ispin_target(1:nstate_target)

DEALLOCATE(psitarget)
DEALLOCATE(ispin_target)

RETURN
END SUBROUTINE attach_prob

!-----cludcmp_d----------------------------------------------                        


SUBROUTINE cludcmp_d(a,n,indx,d,det,ierror)
USE params

!     Lower upper decomposition for a complex matrix:                                     
!     a        matrix to be decomposed, at the end decomposed matrix                       
!     n        actual dimension of matrix                                                  
!!     np       physical dimension of matrix                                               
!     indx     array keeping the permutations from pivoting                               
!     d        sign  of number of row interchanges                                        
!     det      determinant of 'a'                                                          


COMPLEX(DP), INTENT(IN OUT)  :: a(kstate,kstate)
INTEGER, INTENT(IN)          :: n
INTEGER, INTENT(OUT)         :: indx(n)
REAL(DP), INTENT(OUT)        :: d
COMPLEX(DP), INTENT(OUT)     :: det
INTEGER, INTENT(OUT)         :: ierror

!      COMPLEX*8 a(np,np),sum,dum2,newd,det                                                
COMPLEX(DP) :: sum,dum2,newd

INTEGER, PARAMETER :: nmax=500
DOUBLE PRECISION, PARAMETER :: tiny=1.0D-20
INTEGER :: i,imax,j,k
REAL(DP) :: aamax,dum,vv(nmax)
!      write(6,*)'enter LU decomp'            

!WRITE(*,*) 'kstate,n,a:',kstate,n,a(:,:)
                                                
IF(kstate > nmax) STOP 'LUDCMP: too large matrix'
ierror = 0
d=1D0
DO i=1,n
  aamax=0.D0
  DO j=1,n
    IF (ABS(a(i,j)) > aamax) aamax=ABS(a(i,j))
  END DO
!  WRITE(*,*) 'aamax:',aamax
  if (aamax.eq.0D0) stop 'singular matrix in ludcmp'
  IF (aamax == 0.D0) THEN
    ierror = 99
!    WRITE(6,*)'ierror=99'
    RETURN
  END IF
  vv(i)=1D0/aamax
END DO

DO j=1,n
  DO i=1,j-1
    sum=a(i,j)
    DO k=1,i-1
      sum=sum-a(i,k)*a(k,j)
    END DO
    a(i,j)=sum
  END DO
  aamax=0D0
  DO i=j,n
    sum=a(i,j)
    DO k=1,j-1
      sum=sum-a(i,k)*a(k,j)
    END DO
    a(i,j)=sum
    dum=vv(i)*ABS(sum)
    IF (dum >= aamax) THEN
      imax=i
      aamax=dum
    END IF
  END DO
  IF (j /= imax)THEN
    DO k=1,n
      dum2=a(imax,k)
      a(imax,k)=a(j,k)
      a(j,k)=dum2
    END DO
    d=-d
    vv(imax)=vv(j)
  END IF
  indx(j)=imax
  IF(a(j,j) == CMPLX(0D0,0D0)) a(j,j)=tiny
  IF(j /= n)THEN
    dum2=1D0/a(j,j)
    DO i=j+1,n
      a(i,j)=a(i,j)*dum2
    END DO
  END IF
END DO

newd = CMPLX(d,0D0)
!      write(6,*) newd                                                                        
DO i=1,n
!        write(6,*)'i,a(i,i)',i,a(i,i)                                                        
  newd = newd*a(i,i)
END DO
det = newd

RETURN
END SUBROUTINE cludcmp_d


