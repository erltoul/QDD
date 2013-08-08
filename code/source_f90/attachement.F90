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
ALLOCATE(occ_target(nstate))
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

!COMPLEX(DP),ALLOCATABLE :: psitarget(:,:)

COMPLEX(DP) :: overlaps(kstate,kstate),submatr(kstate,kstate)
COMPLEX(DP) :: determinant,tbelement,det,tbacc

COMPLEX(DP) :: psip(kdfull2),psipp(kdfull2)

COMPLEX(DP) :: wfovlp
REAL(DP),ALLOCATABLE :: occ_act(:)
INTEGER :: indx(nstate)             ! index field for LU decomp.                INTEGER :: index(nstate-2)          ! index field for LU decomp.                
INTEGER,ALLOCATABLE :: ipoint(:),ipoi_target(:),ipoi_act(:)
          
COMPLEX(DP) :: temp1,temp2
REAL(DP) :: xxdum,d
INTEGER :: iidum,nexpand,ierror,npoi
LOGICAL,PARAMETER :: ttestb=.false.   ! large test output
LOGICAL,PARAMETER :: ttest=.true.   ! compact test output

!--------------------------------------------------------------                           
         
!WRITE(*,*) 'ispin(nstate)',ispin(nstate)
IF(ttest) WRITE(*,'(a)') 'enter attach_prob'

vcoll=1D0      ! strength of collisional pot.                                   totalprob=0D0
nmatchenergy=0

! prepare actual number of active states
nexpand=0
DO i=1,nstate
  IF(occup(i) > 0.5D0) nexpand=1+nexpand
END DO
IF(ttest) WRITE(*,*) 'number of active states=',nexpand

!ALLOCATE(psitarget(kdfull2,nstate))
ALLOCATE(ispin_target(nstate))
ALLOCATE(occ_act(nstate))
ALLOCATE(ipoint(nexpand),ipoi_target(nexpand),ipoi_act(nexpand))

ispin_target(1:nstate_target)=ispin_temp(1:nstate_target)

! prepare pointers
npoi=0
DO i=1,nstate
  IF(occup(i) > 0.5D0) THEN
    npoi=1+npoi
    ipoint(npoi) = i
  END IF
END DO
IF(npoi .NE. nexpand) STOP 'mismatch in nr. of active TDHF states'
IF(ttestb) WRITE(*,'()') 'IPOINT:',ipoint

npoi=0
DO i=1,nstate-1
  IF(occ_target(i) > 0.5D0) THEN
    npoi=1+npoi
    ipoi_target(npoi) = i
  END IF
END DO
IF(npoi+1 .NE. nexpand) STOP 'mismatch in nr. of target states'
IF(ttestb) WRITE(*,'()') 'IPOI_TARGET:',ipoi_target

IF(ttestb) WRITE(*,'(a,200i3)') 'ispin_target:',ispin_target(1:nstate_target)
!WRITE(*,'(200i3)') ispin(1:nstate)

IF(ttestb) THEN
  WRITE(*,'(a)') 'occups:'
  WRITE(*,'(i3,2f5.1)') (i,occup(i),occ_target(i),i=1,nstate)
  CALL FLUSH(6)
END IF


DO ih0=1,nexpand-1
   ih=ipoi_target(ih0)
   IF(ispin_target(ih) == ispin(nstate)) THEN
!      DO ip1=nstate,nstate_target-1
      DO ip1=1,nstate_target-1
         DO ip2=ip1+1,nstate_target
            IF(ispin_target(ip1) == ispin(nstate).AND.ispin_target(ip2) == ispin(nstate)) THEN
               IF((occ_target(ih) > 0.5).AND.(occ_target(ip1) < 0.5).AND.(occ_target(ip2) < 0.5)) THEN
                  delta_e = spe_target(ip2)+spe_target(ip1)-spe_target(ih)
                  IF(ttestb) WRITE(*,'(a,3i5,4(1pg13.5))') 'ih etc:',ih,ip1,ip2,&
                    spe_target(ih),spe_target(ip1),spe_target(ip2),delta_e
!                  IF(ttestb) WRITE(*,'(3(1pg13.5))') &
!                       occ_target(ih),occ_target(ip1),occ_target(ip1)
                  IF(delta_e > emin_target.AND.delta_e < emax_target) THEN
!                     write(6,*)'ih,ip1,ip2,delta_e',ih,ip1,ip2,delta_e
                     nmatchenergy=nmatchenergy+1
                     IF(ttest) WRITE(*,*) 'nmatch:',nmatchenergy

                     ipoi_act=ipoi_target
                     ipoi_act(nexpand) = ip2
                     ipoi_act(ih0) = ip1
                     IF(ttestb) WRITE(*,'(a,200i3)') 'IPOI_ACT:',ipoi_act

!                     DO nb=1,nstate-1
!                        psitarget(1:nxyz,nb)=psitemp(1:nxyz,nb)
!                     END DO
!                     psitarget(1:nxyz,ih)=psitemp(1:nxyz,ip1)
!                     ispin_target(ih)=ispin_target(ip1)
!                     psitarget(1:nxyz,nstate)=psitemp(1:nxyz,ip2)
!                     ispin_target(nstate)=ispin_target(ip2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             
!     build matrix overlaps of s.p. states                                            

                     overlaps=CMPLX(0D0,0D0)
                     submatr=CMPLX(0D0,0D0)

                     IF(ttestb) write(6,*) &
                        'entering overlap comp. loop'

                     DO i=1,nexpand
                        DO j=1,nexpand
                           IF(ispin_target(i) == ispin(j)) THEN
                              overlaps(i,j) = wfovlp(psitemp(1,ipoi_act(i)),psi(1,ipoint(j)))
                           ELSE
                              overlaps(i,j) = CMPLX(0D0,0D0)
                           END IF
                        END DO
                     END DO

                     IF(ttestb) THEN
                        WRITE(*,*) 'overlaps:'
                        DO i=1,nexpand
                           WRITE(*,'(20(1pg13.5))') overlaps(1:nexpand,i)
                        END DO
                        CALL FLUSH(6)
                     END IF

!                     WRITE(*,*) 'LU decomposition of overlaps matrix'
!                     CALL cludcmp_d(overlaps,nexpand,kstate,indx,d,det,ierror)
                     
                     CALL cludcmp_d(overlaps,nexpand,indx,d,det,ierror)

                     IF(ierror == 99) det = CMPLX(0D0,0D0)
!                     WRITE(6,'(f12.5,4i5,3(1pg13.5))') tfs,nmatchenergy,ih,ip1,ip2,&
!                          delta_e,REAL(det),imag(det)
!                     WRITE(811,'(f12.5,4i8,3(1pg13.5))') tfs,nmatchenergy,ih,ip1,ip2,&
!                          delta_e,REAL(det),imag(det)

!     accumulate total transition matrix element                                     
                     IF(ttestb) WRITE(*,*) ' accumulate transition matrix'
                     tbelement = CMPLX(0D0,0D0)
                     DO i1=1,nexpand
                        i1nn = ipoi_act(i1)
                        DO i2=i1+1,nexpand
                           i2nn = ipoi_act(i2)
!                           IF(ttestb) WRITE(*,'(a,4i3)') &
!                              ' submatrix outer loops: i1,i2=',i1,i2
                           DO j1=1,nexpand
                              j1nn = ipoint(j1)
                              DO j2=j1+1,nexpand
                                 j2nn = ipoint(j2)
                                 IF((ispin_target(i1nn)+ispin_target(i2nn))  &
                                      == ispin(j1nn)+ispin(j2nn)) THEN
!      extract submatrix                     
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
                                       tbacc = CMPLX(0D0,0D0)
                                       DO ind=1,kdfull2
                                          temp1=psitemp(ind,i1nn)*psitemp(ind,i2nn)
                                          temp2=psi(ind,j1nn)*psi(ind,j2nn)
                                          tbacc = CONJG(temp1)*temp2  + tbacc
                                       END DO
                                       tbelement=tbacc*dvol*det+tbelement
                                    END IF
                                 END IF
                              END DO
                           END DO
                        END DO
                     END DO
!     final composition                                                                   
                    IF(ttest) THEN
                       WRITE(*,'(a,3i3,2(1pg13.5))') &
                         ' ih,ip1,ip2,tbelement=',ih,ip1,ip2,tbelement
                       CALL FLUSH(6)
                    END IF
                    totalprob = totalprob + ABS(tbelement*vcoll)**2

                  END IF    ! test on 2p1h energies                                       
               END IF    ! test on 2p1h occupations                                       
            END IF    ! test on ip1 and ip2 spins                                         
         END DO      ! loop on ip2                                                        
      END DO      ! loop on ip1                                                           
   END IF        ! test on ih spin                                                        
END DO! loop on ih   

IF(ttest) THEN
  WRITE(*,'(a,i4,2(1pg13.5))') &
    'nmatchenergy,totalprob=',nmatchenergy,totalprob
  CALL FLUSH(6)
END IF

!DEALLOCATE(psitarget)
DEALLOCATE(ispin_target)
DEALLOCATE(occ_act)
DEALLOCATE(ipoint,ipoi_target,ipoi_act)

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


