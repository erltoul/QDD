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

#include"define.h"                                                                      

!-----init_occ_target--------------------------------------------
SUBROUTINE init_occ_target()
USE params

REAL(DP) :: xxdum,delta_etrgt
INTEGER  :: iidum
LOGICAL,PARAMETER :: ttest=.false.
INTEGER,PARAMETER :: Nmatchmax=500

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
   READ(91,*) binerg
!ELSE
!   STOP ' invalid option in ATTACH_PROB'
!END IF
WRITE(*,'(a,2i5,2(1pg15.7))')  &
  ' istat,irest,binerg,reference_energy:',istat,irest,binerg,reference_energy

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

ALLOCATE(ispin_target(nstate_target))
ALLOCATE(occ_target(nstate_target))
ALLOCATE(spe_target(nstate_target))

REWIND(91)
READ(91,*) xxdum
DO i=1,nstate_target
   READ(91,*)iidum,ispin_target(i),occ_target(i),xxdum,spe_target(i),xxdum
END DO
CLOSE(91)

WRITE(*,*) 'ispin_target:',ispin_target(1:nstate_target)
WRITE(*,*) 'occ_target:',occ_target(1:nstate_target)
WRITE(*,*) 'nclust,nstate,nstate_target:',nclust,nstate,nstate_target
WRITE(*,*) 'ispin:',ispin(1:nstate)
!commPG
!
! The following loops make sense only if the target states are sorted
! such that the 'nstate-1' occupied states come first. For a robust
! scheme, it would be better to loop over all target states and to
! pick particle- and hole-states through 'occup_target'.
!
!commPG

! Calculation of the triples (ih,ip1,ip2) corresponding to 2p1h transitions
! which are possible according to the matching energy condition.
nmatch=0
!Nmatchmax=(nstate-1)*(nstate_target-nstate)*(nstate_target-nstate-1)
ALLOCATE(match(Nmatchmax,3))
DO ih=1,nstate_target ! loop over the holes
   IF(ispin_target(ih) == ispin(nstate) .AND. (occ_target(ih) > 0.5D0)) THEN

      DO ip1=1,nstate_target-1  ! loop over the first particles
         IF(ispin_target(ip1) == ispin(nstate) .AND. (occ_target(ip1) < 0.5D0)) THEN

            DO ip2=ip1+1,nstate_target ! loops over the second particles
               IF(ispin_target(ip2) == ispin(nstate).AND.(occ_target(ip2) < 0.5D0)) THEN

                  delta_e = spe_target(ip2)+spe_target(ip1)-spe_target(ih)
                  IF(ttest) WRITE(*,'(a,3i5,4(1pg13.5))') 'ih etc:',ih,ip1,ip2,&
                    spe_target(ih),spe_target(ip1),spe_target(ip2),delta_e
                  IF(delta_e > emin_target.AND.delta_e < emax_target) THEN
                     nmatch=nmatch+1
                     IF(nmatch>Nmatchmax) STOP "too many matching 2p1h states"
                     match(nmatch,1)=ih  ! index of the created hole
                     match(nmatch,2)=ip1 ! index of the 1st particle state
                     match(nmatch,3)=ip2 ! index of the 2nd particle state, corresponding to the incoming
                                         ! electron, labeled by nstate
                     IF (ttest) WRITE(*,*) 'in match array: ',nmatch,ih,ip1,ip2
                  END IF
               END IF
            END DO
         END IF
      END DO
   END IF
END DO

IF (ttest) THEN
   DO iener=1,nmatch
      write(*,*) 'iener,match: ',iener,match(iener,:)
   END DO
END IF
WRITE(*,*) 'nmatch:',nmatch
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
      READ(90) xxdum,psi_target(1:nxyz,nb)
!      WRITE(*,*) 'nb,xxdum',nb,xxdum
   END DO
END IF

678 CONTINUE
CLOSE(90)
write(6,*)'target wfs read'

RETURN
END SUBROUTINE init_psitarget

!-----attach_prob------------------------------------------------                        

SUBROUTINE attach_prob(totalprob,totalovlp,psi)
USE params
USE util, ONLY:wfovlp
REAL(DP), INTENT(OUT)           :: totalprob,totalovlp
COMPLEX(DP), INTENT(IN)         :: psi(kdfull2,kstate)

!COMPLEX(DP),ALLOCATABLE :: psitarget(:,:)

COMPLEX(DP) :: overlaps(kstate,kstate),submatr(kstate,kstate)
COMPLEX(DP) :: tbelement,det,tbacc,testovlp

!COMPLEX(DP) :: psip(kdfull2),psipp(kdfull2)

INTEGER :: indx(nstate)             ! index field for LU decomp.                
INTEGER :: index(nstate-2)          ! index field for LU decomp.                
INTEGER,ALLOCATABLE :: ipoint(:),ipoi_act(:)
          
COMPLEX(DP) :: temp1,temp2
REAL(DP) :: d
INTEGER :: nexpand,ierror,npoi
LOGICAL,PARAMETER :: ttestb=.false.   ! large test output
LOGICAL,PARAMETER :: ttest=.true.   ! compact test output

!--------------------------------------------------------------                           
         
!WRITE(*,*) 'ispin(nstate)',ispin(nstate)
IF(ttest) WRITE(*,'(a)') 'enter attach_prob'

vcoll=1D0      ! strength of collisional pot.                                   
totalprob=0D0
totalovlp=0D0
totaltestovlp=0D0

IF(ttestb) THEN
  WRITE(*,'(a,200f8.4)') 'OCCUP:',occup(1:nstate)
  WRITE(*,'(a,200f8.4)') 'OCC_TARGET:',occ_target(1:nstate_target)
END IF

! test actual number of active states, only occcupied states allowed
nexpand=0
DO i=1,nstate
  IF(occup(i) > 0.5D0) nexpand=1+nexpand
END DO
IF(ttest) WRITE(*,*) 'number of active states=',nexpand
IF(nexpand .NE. nstate) STOP 'mismatch in nr. of active TDHF states'

!ALLOCATE(psitarget(kdfull2,nstate))
ALLOCATE(ipoint(nstate),ipoi_act(nstate))

! prepare pointers
npoi=0
DO i=1,nstate_target
  IF(occ_target(i) > 0.5D0) THEN        !modPG
    npoi=1+npoi
    ipoint(npoi) = i
 END IF
END DO
IF(npoi+1 .NE. nstate) STOP 'mismatch in nr. of active TDHF states'
ipoint(nstate)=0
IF(ttestb) WRITE(*,'(a,200i3)') 'IPOINT:',ipoint

IF(ttestb) WRITE(*,'(a,200i3)') 'ispin_target:',ispin_target(1:nstate_target)
!WRITE(*,'(200i3)') ispin(1:nstate)

! Loop over the possible 2p1h transitions
DO iener=1,nmatch
   ipoi_act = ipoint
   ih = match(iener,1)
   ip1 = match(iener,2)
   ip2 = match(iener,3)
   ipoi_act(ih) = ip1
   ipoi_act(nstate) = ip2

   IF (ttest) write(*,*) 'ipoi_act: ',ipoi_act

!     build matrix overlaps of s.p. states                                            
   overlaps=CMPLX(0D0,0D0,DP)
   submatr=CMPLX(0D0,0D0,DP)
   
   IF(ttestb) write(6,*) 'entering overlap comp. loop'
   
   DO i=1,nstate
      inn=ipoi_act(i)
      DO j=1,nstate
         IF(ispin_target(inn) == ispin(j)) THEN
            overlaps(i,j) = wfovlp(psi_target(:,inn),psi(:,j))
         ELSE
            overlaps(i,j) = CMPLX(0D0,0D0,DP)
         END IF
      END DO
   END DO

   IF(ttestb) THEN
      WRITE(*,*) 'overlaps:'
      DO i=1,nstate
         WRITE(*,'(20(1pg13.5))') overlaps(1:nstate,i)
      END DO
      CALL FLUSH(6)
   END IF

!
!  as a test, compute determinant of full overlap
!
   submatr(1:nstate,1:nstate) = overlaps(1:nstate,1:nstate)
   CALL cludcmp_d(submatr,nstate,indx,d,det,ierror)
   IF(ierror == 99) det = CMPLX(0D0,0D0,DP)
   totalovlp=ABS(det)**2+totalovlp

!     accumulate total transition matrix element                                     
   IF(ttestb) WRITE(*,*) ' accumulate transition matrix'
   tbelement = CMPLX(0D0,0D0,DP)
   testovlp =  CMPLX(0D0,0D0,DP)
!   WRITE(*,'(a,i5)') ' for iener=',iener
!   WRITE(*,'(a)') ' i1,i2,j1,j2,tbacc,det,tbacc*dvol*det:'
   DO i1=1,nstate
      i1nn = ipoi_act(i1)
      DO i2=i1+1,nstate
         i2nn = ipoi_act(i2)
         !                           IF(ttestb) WRITE(*,'(a,4i3)') &
         !                              ' submatrix outer loops: i1,i2=',i1,i2
         DO j1=1,nstate
            DO j2=j1+1,nstate
               IF((ispin_target(i1nn)+ispin_target(i2nn))  &
                    == ispin(j1)+ispin(j2) ) THEN
                  !      extract submatrix                     
                  ishift = 0
                  DO i=1,nstate
                     IF(i == i1) ishift = 1+ishift
                     IF(i == i2) ishift = 1+ishift
                     jshift = 0
                     DO j=1,nstate
                        IF(j == j1) jshift = 1+jshift
                        IF(j == j2) jshift = 1+jshift
                        IF(i /= i1 .AND. i /= i2 .AND. j /= j1 .AND. j /= j2) THEN
                           submatr(i-ishift,j-jshift) = overlaps(i,j)
                        END IF
                     END DO
                  END DO
                  CALL cludcmp_d(submatr,nstate-2,index,d,det,ierror)
                  IF(ierror == 99) det = CMPLX(0D0,0D0,DP)
                  
                  testovlp = det*(2*MOD(i1,2)-1)*(2*MOD(i2,2)-1) &
                                *(2*MOD(j1,2)-1)*(2*MOD(j2,2)-1) &
                             *(overlaps(i1,j1)*overlaps(i2,j2) &
                               -overlaps(i1,j2)*overlaps(i2,j1)) &
                            + testovlp
                  IF(ispin_target(i1nn).NE.ispin_target(i2nn) ) THEN
!                       IF(mod(j1-j2,2)==0) det = -det
                     det = det*(2*MOD(i1,2)-1)*(2*MOD(i2,2)-1) &
                              *(2*MOD(j1,2)-1)*(2*MOD(j2,2)-1)
                     IF(ispin_target(i1nn).NE.ispin(j1)) det = -det
                     tbacc = CMPLX(0D0,0D0,DP)
                     DO ind=1,kdfull2
                        temp1=psi_target(ind,i1nn)*psi_target(ind,i2nn)
                        temp2=psi(ind,j1)*psi(ind,j2)
                        tbacc = CONJG(temp1)*temp2  + tbacc
                     END DO
                     tbelement=tbacc*dvol*det+tbelement
!                     add = REAL(ABS(tbacc*dvol*det))
!IF(add>1D-5) WRITE(*,'(4i5,6(1pg13.5))') i1,i2,j1,j2,REAL(ABS(tbacc)),REAL(ABS(det)),add
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
   totaltestovlp=ABS(testovlp)**2+totaltestovlp
   
END DO! loop on 2p1h transitions

totaltestovlp=SQRT(totaltestovlp)
totalovlp=SQRT(totalovlp)

IF(ttest) THEN
  WRITE(*,'(a,i4,4(1pg13.5))') &
    'nmatch,totalprob=',nmatch,totalprob
  WRITE(*,'(a,4(1pg13.5))') &
    'test determinants=',totalovlp,totaltestovlp,totaltestovlp/totalovlp
  CALL FLUSH(6)
END IF

DEALLOCATE(ipoint,ipoi_act)

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
  aamax=0D0
  DO j=1,n
    IF (ABS(a(i,j)) > aamax) aamax=ABS(a(i,j))
  END DO
!  WRITE(*,*) 'aamax:',aamax
  if (aamax.eq.0D0) stop 'singular matrix in ludcmp'
  IF (aamax == 0D0) THEN
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
  IF(a(j,j) == CMPLX(0D0,0D0,DP)) a(j,j)=tiny
  IF(j /= n)THEN
    dum2=1D0/a(j,j)
    DO i=j+1,n
      a(i,j)=a(i,j)*dum2
    END DO
  END IF
END DO

newd = CMPLX(d,0D0,DP)
!      write(6,*) newd                                                                        
DO i=1,n
!        write(6,*)'i,a(i,i)',i,a(i,i)                                                        
  newd = newd*a(i,i)
END DO
det = newd

RETURN
END SUBROUTINE cludcmp_d


