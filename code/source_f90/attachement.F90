!#include"define.h"
 
!-----attach_prob------------------------------------------------

SUBROUTINE attach_prob(nmatchenergy,totalprob,psi)
USE params

INTEGER, INTENT(IN OUT)         :: nmatchenergy
REAL(DP), INTENT(OUT)           :: totalprob
COMPLEX(DP), INTENT(IN)         :: psi(kdfull2,kstate)
!IMPLICIT REAL(DP) (a-h,o-z)

!     Package to compute the transition matrix element for
!     electron attachement. The wavefunctions of the
!     target system are read in from file 'save'
!     List parameters are:
!       totalprob    = resulting transition probability

!     The two-body interaction for is taken as zero-range interaction
!     and the overall strength is given as parameter 'Vcoll' below.

!     First reads  data of target wavefunctions for overlap calculations.
!     target wfs stored in 'wfs_target' which is the 'save'
!     file of the target system (caution with
!     dispacement of cluster in box --> computation
!     requires 1 dynamical step of target
!     system to ensure proper location of target in computational box
!     NB: routine usuable in sequential mode only
!c
!     CAUTION: a proper scan of possibly attaching states (2p1h excitations)
!     requires to compute a large span of sp excited levels. In order to save space
!     and avoid expansive copying the set of excited states and occupation numbers
!     are stored in the parent directory.

!INCLUDE "all.inc"


COMPLEX(DP) :: psitarget(kdfull2,kstate)
INTEGER :: ispin_target(10*kstate)
REAL(DP) :: xdum(kstate),idum(kstate)

REAL(DP) :: occ_target(10*kstate)
REAL(DP) :: eps_target(10*kstate)

!      complex*8 overlaps(kstate,kstate),submatr(kstate,kstate)
!      complex*8 determinant,tbelement,det
COMPLEX(DP) :: overlaps(kstate,kstate),submatr(kstate,kstate)
COMPLEX(DP) :: determinant,tbelement,det
!      complex determinant,tbelement,det
!      complex overlaps(kstate,kstate),submatr(kstate,kstate)
COMPLEX(DP) :: psip(kdfull2),psipp(kdfull2)

COMPLEX(DP) :: psidum(kdfull2)

EXTERNAL wfovlp
COMPLEX(DP) :: wfovlp
INTEGER :: indx(nstate)             ! index field for LU decomp.
COMPLEX(DP) :: cdum

!--------------------------------------------------------------

WRITE(6,*) 'Entering attach_prob'

vcoll=1D0      ! strength of collisional pot.

!      complex psitarget(kdfull2,kstate)
!      dimension ispin_target(kstate)
!      dimension xdum(kstate),idum(kstate)

WRITE(6,*) 'opening unit 91'
OPEN(UNIT=91,STATUS='unknown',FORM='formatted', FILE='occ_eps_target')
IF(istat == 0.AND.irest == 0) READ(91,*) xdum,xdum
IF(istat == 1.AND.irest == 0) READ(91,*) binerg,xdum
IF(istat == 0.AND.irest > 0) READ(91,*) binerg,reference_energy

!      aver_estar  = etot - binerg
aver_estar  = reference_energy - binerg
delta_etrgt = 3.d0*h2m/4.d0/scatterelectronw**2
emin_target = aver_estar - delta_etrgt
emax_target = aver_estar + delta_etrgt
WRITE(6,*)'aver_estar,delta_etrgt,emin_target,emax_target',  &
    aver_estar,delta_etrgt,emin_target,emax_target

!      write(6,*) 'opening unit 91'
!        open(unit=91,status='unknown',form='formatted',
!    &        file='occ_eps_target')

nstate_target=0
DO i=1,10000
  READ(91,*,END=899)xxdum,iidum,ispin_target(i), occ_target(i),eps_target(i)
!            write(6,'(a,2(1x,i2),2(1x,e12.5) )')
!     &      'target state #, spin, occup,energy',
!     &      i,ispin_target(i),occ_target(i),eps_target(i)
  nstate_target=nstate_target+1
END DO
899     CONTINUE
CLOSE(91)
!         write(6,*)'nstate,nstate_target',nstate,nstate_target
!         do i=1,nstate
!         write(6,*)i,ispin(i),occup(i)
!         enddo

totalprob=0.
nmatchenergy=0

IF(nstate <= 1) RETURN

!         if(nstate.gt.1) then
DO ih=1,nstate-1
  
  IF(ispin_target(ih) == ispin(nstate)) THEN
    
    DO ip1=nstate,nstate_target-1
      DO ip2=ip1+1,nstate_target
        IF(ispin_target(ip1) == ispin(nstate).AND.  &
              ispin_target(ip2) == ispin(nstate)) THEN
          
          IF((occ_target(ih) > 0.5).AND. (occ_target(ip1) < 0.5).AND.  &
                (occ_target(ip2) < 0.5)) THEN
            
            delta_e = eps_target(ip2)+eps_target(ip1)-eps_target(ih)
!               write(6,*)'ih,ip1,ip2,delta_e',ih,ip1,ip2,delta_e
            
            IF(delta_e > emin_target.AND.delta_e < emax_target) THEN
              nmatchenergy=nmatchenergy+1
              
!          write(6,*)'e*,de,emin,emax',
!     &       aver_estar,delta_etrgt,emin_target,emax_target
!          write(6,*)'ih,ip1,ip2,delta_e',ih,ip1,ip2,delta_e
!          write(6,*) 'opening unit 90'
              OPEN(UNIT=90,STATUS='unknown',FORM='unformatted',  &
                  FILE='../wfs_target')
              
!  read the iteration where the data has been saved last:
              READ(90) iact_t,nstate_target,nclust_target,  &
                  nion_target,nspdw_target
!         write(6,*) 'target: nstate_t,nclust_t,nion_t,nspdw_t'
!         write(6,*) nstate_target,nclust_target,
!     &              nion_target,nspdw_target
!         write(6,*) 'actual: nstate,nclust,nion,nspdw,kstate'
!         write(6,*) nstate,nclust,nion,nspdw,kstate
              IF((nion_target /= nion).OR.(nspdw_target /= nspdw)) THEN
                WRITE(6,*)'attachment calculation impossible:'
                WRITE(6,*)'mismatch betw # of ions or spins'
                STOP
              END IF
              
              DO ipp=1,nstate_target
!            read(60) occ_target(i)
                READ(90) xdum(i)
              END DO
              
!  read wavefunctions:
              
              IF(nclust > 1)THEN
                DO nb=1,nstate-1
                  READ(90) (psitarget(i,nb),i=1,nxyz)
                END DO
                DO ipp=nstate,nstate_target
                  READ(90) (psidum(i),i=1,nxyz)
                  IF(ipp == ip1) THEN
                    DO ixyz=1,nxyz
                      psitarget(ixyz,ih)=psidum(ixyz)
                    END DO
                    ispin_target(ih)=ispin_target(ipp)
                  END IF
                  IF(ipp == ip2) THEN
                    DO ixyz=1,nxyz
                      psitarget(ixyz,nstate)=psidum(ixyz)
                    END DO
                    ispin_target(nstate)=ispin_target(ipp)
                  END IF
                END DO
                
!   read spins
                
                GO TO 678
                DO i=1,nstate_target
                  ispin_target(i) = 0
                END DO
                DO nb=1,nstate-1
                  READ(90) ispin_target(nb),iidum,iddum,iidum
                END DO
                DO ipp=nstate,nstate_target
                  READ(90) idum(ipp),iidum,iddum,iidum
                  IF(ipp == ip1) THEN
                    DO ixyz=1,nxyz
                      ispin_target(ih)=idum(ipp)
                    END DO
                  END IF
                  IF(ipp == ip2) THEN
                    DO ixyz=1,nxyz
                      ispin_target(nstate)=idum(ipp)
                    END DO
                  END IF
                END DO
              END IF
              
              678  CONTINUE
              CLOSE(90)
!      write(6,*)'target wfs and occ/energ read'
              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     build matrix overlaps of s.p. states
              
              DO i=1,kstate
                DO j=1,kstate
                  overlaps(i,j)=CMPLX(0.d0,0.d0)
                  submatr(i,j)=CMPLX(0.d0,0.d0)
                END DO
              END DO
              
              nexpand = nstate
              
!      write(6,*)'entering overlap comp. loop with nexpand=', nexpand
              DO i=1,nexpand
                DO j=1,nexpand
!          write(6,*)'i,j,ispin_target_i,ispin_j',
!     &               i,j,ispin_target(i),ispin(j)
                  IF(ispin_target(i) == ispin(j)) THEN
                    overlaps(i,j) = wfovlp(psitarget(1,i),psi(1,j))
                  ELSE
                    overlaps(i,j) = CMPLX(0.d0,0.d0)
                  END IF
                  
!          write(6,'(a,2i5,2(1pg13.5))')' overlap=',i,j,overlaps(i,j)
!          write(66,'(a,i1,a,i1,a,1e12.5,a,1e12.5),a')
!     &    '       overlaps(',i,',',j,')=cmplx(',
!     &            real(overlaps(i,j)),',',
!     &            imag(overlaps(i,j)),')'
                END DO
              END DO

!      write(6,*)'Vf nexpand bef. entering LU', nexpand
!old      call cludcmp_d(overlaps,nexpand,nexpand,indx,d,det,ierror)
              CALL cludcmp_d(overlaps,nexpand,kstate,indx,d,det,ierror)
              
              IF(ierror == 99) det = CMPLX(0D0,0D0)
              WRITE(6,'(f12.5,4i5,3(1pg13.5))') tfs,nmatchenergy,ih,ip1,ip2,&
!     & ispin_target(ip1),ispin_target(ip2),  &
              delta_e,REAL(det),imag(det)
              WRITE(811,'(f12.5,4i5,3(1pg13.5))') tfs,nmatchenergy,ih,ip1,ip2,&
!     & ispin_target(ip1),ispin_target(ip2),  &
              delta_e,REAL(det),imag(det)
              
!     accumulate total transition matrix element
              
              tbelement = CMPLX(0D0,00D0)
              DO i1=1,nexpand
                DO i2=i1+1,nexpand
                  DO j1=1,nexpand
                    DO j2=j1+1,nexpand
                      i1nn = i1
                      i2nn = i2
                      j1nn = j1
                      j2nn = j2
!        if(mtarget(i1)+mtarget(i2).eq.mpsiexp(j1)+mpsiexp(j2)) then
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
!                write(6,'(8i3)')
!     &            i1,i2,j1,j2,i,j,i-ishift,j-jshift
                            END IF
                          END DO
                        END DO
!        write(6,'(a,6i5)') ' submatrix composed: coeffs,shifts=',
!     &    i1,i2,j1,j2,ishift,jshift
!old          call cludcmp_d(submatr,nexpand-2,nexpand,indx,d,det,ierror)
                        CALL cludcmp_d(submatr,nexpand-2,kstate,indx,d,det,ierror)
                        IF(ierror == 99) det = CMPLX(0D0,0D0)
!          write(6,'(a,2(1pg13.5))') ' det(submatrix)=',det
!                                    compute matrix element of "delta"
                        IF(ierror == 0) THEN
                          DO ind=1,kdfull2
                            psip(ind)  = psitarget(ind,i1nn)*psitarget(ind,i2nn)
                            psipp(ind) = psi(ind,j1nn)*psi(ind,j2nn)
                          END DO
                          
                          tbelement = tbelement + wfovlp(psip,psipp)*det
!            write(6,'(a,2(1pg13.5))') ' tbelement=',tbelement
                        END IF
                      END IF
                    END DO
                  END DO
                END DO
              END DO
!      write(6,'(a)') ' fourfold loop done'
              
!     final composition
              
              totalprob = totalprob + ABS(tbelement*vcoll)**2
              
            END IF    ! test on 2p1h energies
          END IF    ! test on 2p1h occupations
        END IF    ! test on ip1 and ip2 spins
      END DO      ! loop on ip2
    END DO      ! loop on ip1
  END IF        ! test on ih spin
END DO          ! loop on ih
!      endif      ! test on nstate>1

RETURN
END SUBROUTINE attach_prob

!-----cludcmp_d----------------------------------------------


SUBROUTINE cludcmp_d(a,n,np,indx,d,det,ierror)

!     Lowwer upper decomposition for a complex matrix:
!     a        matrix to be decomposed, at the end decomposed matrix
!     n        actual dimension of matrix
!     np       physical dimension of matrix
!     indx     array keeping the permutations from pivoting
!     d        sign  of number of row interchanges
!     det      determinant of 'a'


COMPLEX, INTENT(IN OUT)                  :: a(np,np)
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: np
INTEGER, INTENT(OUT)                     :: indx(n)
DOUBLE PRECISION, INTENT(OUT)            :: d
COMPLEX, INTENT(OUT)                     :: det
INTEGER, INTENT(OUT)                     :: ierror

!      COMPLEX*8 a(np,np),sum,dum2,newd,det
COMPLEX :: sum,dum2,newd

INTEGER, PARAMETER :: nmax=500
DOUBLE PRECISION, PARAMETER :: tiny=1.0D-20
INTEGER :: i,imax,j,k
DOUBLE PRECISION :: aamax,dum,vv(nmax)
!      write(6,*)'enter LU decomp'
IF(np > nmax) STOP 'LUDCMP: too large matrix'
ierror = 0
d=1D0
DO i=1,n
  aamax=0D0
  DO j=1,n
    IF (ABS(a(i,j)) > aamax) aamax=ABS(a(i,j))
  END DO
!        if (aamax.eq.0D0) stop 'singular matrix in ludcmp'
  IF (aamax == 0D0) THEN
    ierror = 99
    WRITE(6,*)'ierror=99'
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
