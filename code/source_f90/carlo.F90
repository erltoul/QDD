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
 

SUBROUTINE init_simann()

!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
INTEGER :: ifall

    ionsin    = 3
    iknow     = 0
    facann    = 1.5D0
    nrun      = 100
    nloop1    = 100
    nloop2    = 200
    cptemp    = 0.008D0
    delpos    = 0.2D0
    ERR       = 6.0D-5
    ifall     = 121466
    trfac2    = 0.03D0
    prfac2    = 0.05D0
    errsim    = 0.001D0
    ncsim     = 10
    errtot    = 0.002D0
    ncon      = 8
    erfac1    = 0.1D0
    trfac1    = 0.08D0
    prfac1    = 0.05D0
    errks0    = 5.0D-3

! Initialize random seed.
    CALL RANDOM_SEED()
    CALL RANDOM_SEED(SIZE=isize_seed)
    ALLOCATE(rand_seed(isize_seed))

    IF(myn == 0)THEN

      OPEN(UNIT=44,STATUS='old',FORM='formatted', FILE='ionen-in.'//outnam)
      WRITE(6,'(/a)') 'PARAMETERS FOR METROPOLIS:'
      READ(44,*) ionsin,iknow,facann
      WRITE(6,'(a,2i2,f7.3)') ' ionsin,iknow,facann: ',ionsin, iknow,facann
      READ(44,*) nrun,nloop1,nloop2
      WRITE(6,'(a,3i4)') ' nrun,nloop1,nloop2: ',nrun,nloop1 ,nloop2
      READ(44,*) cptemp,delpos,ERR,errks0
      WRITE(6,'(a,4f12.6)') ' cptemp,delpos,err,errks0: ',cptemp  &
          ,delpos,ERR,errks0
      READ(44,*) ifall
      WRITE(6,'(a,i10)') ' ifall (argument no longer used):',ifall
      READ(44,*) trfac2,prfac2,errsim,ncsim
      WRITE(6,'(a,3f8.3,i3)') ' trfac2,prfac2,errsim,ncsim: '  &
          ,trfac2,prfac2,errsim,ncsim
      READ(44,*) errtot,ncon
      WRITE(6,'(a,f8.3,i3)') ' errtot,ncon: ' ,errtot,ncon
      READ(44,*) erfac1,trfac1,prfac1
      WRITE(6,'(a,3f8.3)') ' erfac1,trfac1,prfac1: ' ,erfac1,trfac1,prfac1
      WRITE(6,*) ' '
      WRITE(6,*) ' '
      CLOSE(44)

! Take random seed from first node.
      CALL RANDOM_SEED(GET=rand_seed)

    END IF
    
#if(parayes)
    CALL comm_simann()
#endif

! Start all nodes with the same random seed.
    CALL RANDOM_SEED(PUT=rand_seed)

END SUBROUTINE init_simann

!     *************************************

SUBROUTINE simann(psir,rho,aloc)

!     *************************************

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: psir(kdfull2,kstate)
REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)
!REAL(DP), INTENT(IN OUT)                     :: akv(kdfull2)





!     save starting values of sim. annealing for further use:

delps0 = delpos
cptem0 = cptemp

!     loop over various runs of the metropolis part:

IF(istat == 1) THEN
  CALL resume(psir,outnam)
END IF
DO jrun = 1, nrun
  
!     reset annealing schedule for further sim. annealing:
  
  cptemp = cptem0
  delpos = delps0
!         nloop2 = nloop2/(2.0*facann)
  errks  = errks0
  nyes = 0               !counter for convergence in main loop
  ebold = 1D20
  hph2 = 1D0
  
  DO loop1 = 1,nloop1
    
!     reach convergence in Kohn-Sham loop
    
!g            call calcpseudo(rho)
    CALL calcpseudo()
    
!old            sumion = 0.0
!old            do ion1=2,nion
!old               do ion2=1,ion1-1
!old                  dist   = (cx(ion1)-cx(ion2))**2+(cy(ion1)-cy(ion2))**2
!old     &                 +(cz(ion1)-cz(ion2))**2
!old                  sumion = e2*ch(np(ion1))**2/sqrt(dist) + sumion
!old               enddo
!old            enddo
!old            ecorr=sumion
    ecorr = energ_ions()
    
!     calculate the projectors around the new ionic positions
    
    IF(ipsptyp == 1) THEN
      DO ion=1,nion
        CALL calc_proj(cx(ion),cy(ion),cz(ion),cx(ion),cy(ion),cz(ion),ion)
!        IF(nrow(np(1)) == 2) CALL calpr2(cx(ion),cy(ion),cz(ion),ion)
!        IF(nrow(np(1)) == 3) CALL calpr3(cx(ion),cy(ion),cz(ion),ion)
!        IF(nrow(np(1)) == 4) CALL calpr4(cx(ion),cy(ion),cz(ion),ion)
        WRITE(6,'(a,i2,a,3f9.4)') 'ion',ion,'=', cx(ion),cy(ion),cz(ion)
      END DO
      CALL checkproj(ion)
    END IF
    ismax = 200
    CALL statit(psir,rho,aloc)
    
    CALL fmtv_fld(psir,rho,0)
    
    delbin = ABS(ebold-binerg)
    IF(loop1 >= 2) THEN
      WRITE(6,'(/a,f8.5,a,f8.5)') 'new binding energy =',binerg,  &
          ' energy difference with last loop =',binerg-ebold
    END IF
    
    IF(delbin < errtot) THEN
      nyes = nyes + 1
      IF(nyes == ncon) GO TO 999 !leave metropolis
    ELSE
      nyes = 0
    END IF
    WRITE(6,'(a,f9.4,i3)') 'deltaE, nyes =',binerg-ebold,nyes
    ebold = binerg
    
!     if convergence is not achieved, reduce annealing schedule
!     and Kohn-Sham accuracy 'errks':
    
    errks = errks*(1D0-erfac1)
    IF(errks < ERR) errks = ERR !maximum accuracy
    cptemp = cptemp*(1D0-trfac1)
    delpos = delpos*(1D0-prfac1)
    
    WRITE(6,'(a)') '-->NOW STARTING TO OPTIMIZE THE IONIC POSITIONS keep waiting ...'
    
    CALL minpos(rho,psir,aloc)
    
    CALL cenmass()
    CALL view3d()
    
    OPEN(27,POSITION='append',FILE='pminpos.'//outnam)
    DO ion=1,nion
      WRITE(27,'(i3,4f13.5)') loop1,cx(ion),cy(ion),cz(ion),binerg
    END DO
    CALL flush(27)
    
    WRITE(6,'(A/)') '-->ROLLING THE DICE COMPLETED'
  END DO
  
!     !if you come here, energetic accuracy is finally achieved
  
  999     CONTINUE
  WRITE(6,'(A)') ' NORMAL TERMINATION IN METROPOLIS'
END DO
CALL rsave(psir,outnam)

RETURN
END SUBROUTINE simann
!     *************************************

SUBROUTINE minpos(rho,psimc,aloc)

!     *************************************

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN)                         :: rho(kdfull2)
REAL(DP), INTENT(IN)                         :: psimc(kdfull2,kstate)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)
!REAL(DP), INTENT(IN OUT)                     :: akv(kdfull2)

REAL(DP),ALLOCATABLE ::  q1(:)
REAL(DP),ALLOCATABLE :: rhoion(:)

REAL(DP) :: rand0, rand3(3)

LOGICAL :: ans

ALLOCATE(q1(kdfull2),rhoion(kdfull2))

eloctot = 0D0
enoloctot = 0D0
eiontot = 0D0
DO ion=1,nion
  
  xion=cx(ion)
  yion=cy(ion)
  zion=cz(ion)
  
  CALL rhopsg(xion,yion,zion,rhoion,ion)
  eloc(ion) = 0D0
  DO i = 1, nxyz
    eloc(ion) = eloc(ion) + rho(i)*rhoion(i)
  END DO
  eloc(ion)=eloc(ion)*dvol
  eloctot = eloctot + eloc(ion)
  
  IF(ipsptyp == 1) THEN
    CALL calc_proj(cx(ion),cy(ion),cz(ion),cx(ion),cy(ion),cz(ion),ion)
!    IF(nrow(np(1)) == 2) CALL calpr2(xion,yion,zion,ion)
!    IF(nrow(np(1)) == 3) CALL calpr3(xion,yion,zion,ion)
!    IF(nrow(np(1)) == 4) CALL calpr4(xion,yion,zion,ion)
    enoloc(ion) = 0D0
    DO nb=1,nstate
      CALL nonlocalr(psimc(1,nb),q1)
      sumnl = 0D0
      DO i=1,ifin(ion)
        ii = icount(i,ion)
        sumnl  = sumnl + psimc(ii,nb)*q1(ii)
      END DO
      enoloc(ion) = enoloc(ion) + sumnl*occup(nb)*dvol
    END DO
    enoloctot = enoloctot + enoloc(ion)
  END IF
  
  DO j= ion+1, nion
    dist   = SQRT((cx(j)-xion)**2+(cy(j)-yion)**2 +(cz(j)-zion)**2)
    eion(ion,j) = e2*ch(np(ion))**2/dist
    eiontot     = eiontot + eion(ion,j)
  END DO
END DO

etot = -eloctot + enoloctot + eiontot
eold = etot

delps3 = delpos
teml3  = cptemp
nmin   = 0        ! counts loops with energy change < errsim

!     NOW THE FUN BEGINS
!     the outer dice loop:

diftot = 0D0              ! total energy gain through moving ions
difii  = 0D0              ! ionic energy gain through moving ions

DO loop2 = 1,nloop2
  
  ndown  = 0             ! steps downwards
  nbut   = 0             ! steps upwards although energy gets worse
  nup    = 0             ! steps without any success
  oh     = 0.5D0
  
!     the inner loop:
  
  DO loop3 = 1,nion
    
!     roll the dice to decide which ion is to be moved:

    CALL RANDOM_NUMBER(rand0)
    ionvar = INT(nion*rand0)+1
    
!     vary the coordinates of this ion:
    
    5          CALL RANDOM_NUMBER(rand3)
    deltax = 2.*delps3*(rand3(1)-oh)
    deltay = 2.*delps3*(rand3(2)-oh)
    deltaz = 2.*delps3*(rand3(3)-oh)
    delta2 = deltax*deltax + deltay*deltay + deltaz*deltaz
    IF (delta2 > delps3*delps3) GO TO 5
    
    IF(mxforce == 1) deltax = 0D0
    IF(myforce == 1) deltay = 0D0
    IF(mzforce == 1) deltaz = 0D0
    
    xvar   = cx(ionvar) + deltax
    yvar   = cy(ionvar) + deltay
    zvar   = cz(ionvar) + deltaz
    
!     Ion within sensible grid area?
    
    radone = 3D0
    facmo =  2D0
    
    xred2 = xvar+facmo*radone
    yred2 = yvar+facmo*radone
    zred2 = zvar+facmo*radone
    xled2 = xvar-facmo*radone
    yled2 = yvar-facmo*radone
    zled2 = zvar-facmo*radone
    
    gxr = REAL(nx,DP)*dx
    gyr = REAL(ny,DP)*dy
    gzr = REAL(nz,DP)*dz
    gxl = -gxr
    gyl = -gyr
    gzl = -gzr
    
    IF((gxr <= xred2).OR.(gxl >= xled2)) GO TO 5
    IF((gyr <= yred2).OR.(gyl >= yled2)) GO TO 5
    IF((gzr <= zred2).OR.(gzl >= zled2)) GO TO 5
    
!     calculate new one-ion-electron  energy of ion no. 'ionvar'
    
    CALL rhopsg(xvar,yvar,zvar,rhoion,ionvar)
    eloc(0) = 0D0
    DO i = 1, nxyz
      eloc(0) = eloc(0) + rho(i)*rhoion(i)
    END DO
    eloc(0)=eloc(0)*dvol
    
    IF(ipsptyp == 1) THEN
      CALL calc_proj(cx(ion),cy(ion),cz(ion),cx(ion),cy(ion),cz(ion),ion)
!      IF(nrow(np(1)) == 2) CALL calpr2(xvar,yvar,zvar,0)
!      IF(nrow(np(1)) == 3) CALL calpr3(xvar,yvar,zvar,0)
!      IF(nrow(np(1)) == 4) CALL calpr4(xvar,yvar,zvar,0)
      enoloc(0) = 0D0
      DO nb=1,nstate
        CALL nonlocalr(psimc(1,nb),q1)
        sumnl = 0D0
        DO i=1,ifin(0)
          ii = icount(i,0)
          sumnl  = sumnl + psimc(ii,nb)*q1(ii)
        END DO
        enoloc(0) = enoloc(0) + sumnl*occup(nb)*dvol
      END DO
    END IF
    
!     compute energy-difference diffen=etvar-etold caused by the move:
    
    CALL move(ionvar,xvar,yvar,zvar,diffen,dii)
    
!     The mysterious decision: go down or call the ORACLE
    
    ans = (diffen < 0D0)
    IF (ans) THEN
      ndown = ndown + 1 !good call
    ELSE
      iknowl=iknow
      CALL metrop(diffen,teml3,iknowl,ans)
      IF (ans)   nbut = nbut + 1
    END IF
    IF(ans) THEN
      diftot = diftot + diffen
      difii  = difii  + dii
    END IF
    
!     If 'yes', swap coordinates and energies
    
    IF (ans) THEN
      cx(ionvar) = xvar
      cy(ionvar) = yvar
      cz(ionvar) = zvar
      eloc(ionvar)   = eloc(0)
      enoloc(ionvar) = enoloc(0)
      DO ion = 1,ionvar-1
        eion(ion,ionvar) = eiinew(ion)
      END DO
      DO ion = ionvar+1,nion
        eion(ionvar,ion) = eiinew(ion)
      END DO
    ELSE
      nup = nup + 1    !bad call
    END IF
  END DO                  ! end of dice loop at given T
  
!     NEW ONE-ION ENERGY AFTER ORACLE:
  
  eloctot = 0D0           !reset energies for new calculation
  eiontot = 0D0
  enoloctot  = 0D0
  
  DO i = 1, nion
    eloctot = eloctot + eloc(i)
    enoloctot = enoloctot + enoloc(i)
    DO j = i+1, nion
      eiontot = eiontot + eion(i,j)
    END DO
  END DO
  eionic = -eloctot + enoloctot + eiontot
  
!     lower annealing parameters for new dice loop:
  
  delps3 = delps3*(1D0-prfac2)
  teml3  = teml3*(1D0-trfac2)
  
!     Criteria for convergence: end of loop if ionic energy has been
!     changed for 'ncsim' consecutive times by less than 'errsim':
  
  IF (ABS(eold-eionic) < errsim) THEN
    nmin = nmin + 1
  ELSE
    nmin = 0
  END IF
  IF (nmin == ncsim) GO TO 99
  
  eold = eionic
END DO                     !end of annealing

99   CONTINUE

!     history of dice loops:

WRITE(6,'(/a,i3)')'no. of dice loops: ',loop2
WRITE(6,'(a,f12.4)') 'total energy gain= ',diftot
WRITE(6,'(a,f12.4)') 'ionic energy gain= ',difii
WRITE(6,'(3x,a,/,5x,a,2f8.5,/,5x,a,3g15.6)') 'Leaving minpos with:',  &
    'cptemp, delpos = ',teml3,delps3, 'eloctot,enoloctot,eiontot= ',  &
    eloctot,enoloctot,eiontot
WRITE(6,'(5x,a,f12.4)') 'eionic= ',eionic

DEALLOCATE(q1,rhoion)

RETURN
END SUBROUTINE minpos


!     *************************************************

SUBROUTINE move(ionvar,xvar,yvar,zvar,diffen,dii)

!     *************************************************

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

!     Computes energy-difference 'diffen' resulting from
!     moving 'ion' to (xvar,yvar,zvar).
!     Prepare energies on fieldindex '0' for a possible swapping
!     of coordinates after oracle.
!     eiinew = vector of inverse distances

!     calculate new ion-ion-energy of ion no. 'ionvar' with
!     all the other ions

eiivar = 0D0
DO ion = 1,nion
  IF (ion /= ionvar) THEN
    dist   = SQRT((cx(ion)-xvar)**2+(cy(ion)-yvar)**2 +(cz(ion)-zvar)**2)
    eiinew(ion) = e2*ch(np(ion))**2/dist
    eiivar = eiivar + eiinew(ion)
  END IF
END DO

!     calculate old ion-ion-energy of ion no. 'ionvar' with
!     all the other ions

eionold = 0D0
DO ion = 1,ionvar-1
  eionold = eionold + eion(ion,ionvar)
END DO
DO ion = ionvar+1,nion
  eionold = eionold + eion(ionvar,ion)
END DO

!     compare the energies for the oracle

etold  = eionold - eloc(ionvar) + enoloc(ionvar)
etvar  = eiivar - eloc(0) + enoloc(0)
diffen = etvar - etold
dii    = eiivar - eionold
dieloc = eloc(0) - eloc(ionvar)
dienl  = enoloc(0) - enoloc(ionvar)

RETURN
END SUBROUTINE move

!     *******************

REAL(DP) FUNCTION ran0(idum)
USE params, ONLY: DP
IMPLICIT REAL(DP) (A-H,O-Z)

!     *******************

!     numerical recipes function ran0 (p. 270):
!     creates a random number from input 'idum'


INTEGER, INTENT(OUT)                     :: idum

!REAL(8):: ran0              !danger, do not declare ran0 as real*8
INTEGER, PARAMETER :: ia=16807
INTEGER, PARAMETER :: im=2147483647
REAL(DP), PARAMETER :: am=1D0/im
INTEGER, PARAMETER :: iq=127773
INTEGER, PARAMETER :: ir=2836
INTEGER, PARAMETER :: mask=123459876
INTEGER :: k

idum=IEOR(idum,mask)
k=idum/iq
idum=ia*(idum-k*iq)-ir*k
IF (idum < 0) idum=idum+im
ran0=am*idum
idum=IEOR(idum,mask)

RETURN
END FUNCTION ran0

!     **************************************

SUBROUTINE metrop(diffen,t,iknowi,ans)

!     **************************************

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)



REAL(DP), INTENT(IN)                         :: diffen
REAL(DP), INTENT(IN)                         :: t
INTEGER, INTENT(IN)                      :: iknowi
LOGICAL, INTENT(OUT)                     :: ans
REAL(DP) :: rand0


IF (t > 0D0) THEN
  expon = diffen/t
  IF (expon < 80D0) THEN
    bolfac = EXP(-expon)
  ELSE
    bolfac = 0D0
  END IF
ELSE
  bolfac = 0D0
END IF

IF(iknowi == 0) THEN
  CALL RANDOM_NUMBER(rand0)
  ans = (rand0 < bolfac)
  
!     proper choice whenever the user has not the
!     faintest idea about the shape of the cluster!
!     -> provides the possibility to escape a local
!     minimum in order to reach the global one.
  
ELSE IF(iknowi == 1) THEN
  IF(loop1 == 0) THEN
    ans = (0.95D0 < bolfac) ! provides good starting point
  ELSE IF(loop1 > 0) THEN
    ans = ((1.8D0 - 0.955D0**(LOG(loop1+0.00001D0))) < bolfac)
    
!     this choice cares for a fast iteration ->
!     DANGER: use only if starting point is chosen well,
!     meaning the starting configuration is not too
!     far away from the desired final result.
!     If this is not the case, being trapped in a local
!     minimum is unavoidable!!
    
  END IF
END IF

RETURN
END SUBROUTINE metrop

!     ********************

SUBROUTINE cenmass()

!     ********************

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP) :: cenmas(3)
REAL(DP) :: pos(ng,3)

DO ion=1,nion
  pos(ion,1) = cx(ion)
  pos(ion,2) = cy(ion)
  pos(ion,3) = cz(ion)
END DO

!   compute center of mass

cenmas(1) = zero
cenmas(2) = zero
cenmas(3) = zero
gamu      = zero

DO ion = 1,nion
  DO icoo = 1,3
    cenmas(icoo) = cenmas(icoo) + amu(np(ion))*pos(ion,icoo)
  END DO
  gamu = gamu + amu(np(ion))
END DO

DO icoo = 1,3
  cenmas(icoo) = cenmas(icoo)/gamu
END DO

!     transform to center of mass coordinates

DO ion = 1,nion
  DO icoo = 1,3
    pos(ion,icoo) = pos(ion,icoo) - cenmas(icoo)
  END DO
END DO

DO ion=1,nion
  cx(ion) = pos(ion,1)
  cy(ion) = pos(ion,2)
  cz(ion) = pos(ion,3)
END DO

RETURN
END SUBROUTINE cenmass
