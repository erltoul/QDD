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
 
!------------------------------------------------------------

SUBROUTINE heatcluster(totemp)
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
!     This routine heats up the cluster.
!     The momenta of the cluster ions
!     are reset to values such that the kinetic energy
!     corresponds to an instantaneous "thermal" energy toTemp.
!     toTemp is in Kelvin!!





RETURN
END SUBROUTINE heatcluster
!------------------------------------------------------------


!-----reset_ions--------------------------------------------

SUBROUTINE reset_ions()

!     Resets all ionic momento to zero, is used in cooling.

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

!-----------------------------------------------------------

IF ((ekionold-ekion) > 0D0) THEN
  jekion=jekion+1
  DO ion=1,nion
    cpx(ion)=0D0
    cpy(ion)=0D0
    cpz(ion)=0D0
  END DO
#if(raregas)
  DO ion=1,nc
    pxc(ion)=0D0
    pyc(ion)=0D0
    pzc(ion)=0D0
    pxe(ion)=0D0
    pye(ion)=0D0
    pze(ion)=0D0
  END DO
  DO ion=1,nk
    pxk(ion)=0D0
    pyk(ion)=0D0
    pzk(ion)=0D0
  END DO
#endif
  icooltimes=icooltimes+1
  WRITE(6,*)'jekion=',jekion
  ekionold=0D0
END IF
IF(ekionold < ekion) ekionold = ekion

RETURN
END SUBROUTINE reset_ions

!-----energkin_ions--------------------------------------------

FUNCTION enerkin_ions()

!     Kinetic energy of ions

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

!     calculate kinetic energy of ions

sumekinion=0D0
!      if (nion.gt.0) then
DO ion=1,nion
  ek=cpx(ion)*cpx(ion)+cpy(ion)*cpy(ion)+cpz(ion)*cpz(ion)
  ek=ek/(2D0*1836D0*amu(np(ion))*ame)
  sumekinion=sumekinion+ek
END DO
!      endif
#if(raregas)
DO  ion=1,nc
  ek=pxc(ion)*pxc(ion)
  ek=ek+pyc(ion)*pyc(ion)
  ek=ek+pzc(ion)*pzc(ion)
  ek=ek/(2D0*1836D0*mion*ame)
  sumekinion=sumekinion+ek
END DO

DO  ion=1,NE
  ek=pxe(ion)*pxe(ion)
  ek=ek+pye(ion)*pye(ion)
  ek=ek+pze(ion)*pze(ion)
  ek=ek/(2D0*1836D0*me*ame)
  sumekinion=sumekinion+ek
END DO

DO  ion=1,nk
  ek=pxk(ion)*pxk(ion)
  ek=ek+pyk(ion)*pyk(ion)
  ek=ek+pzk(ion)*pzk(ion)
  ek=ek/(2D0*1836D0*mkat*ame)
  sumekinion=sumekinion+ek
END DO
#endif
enerkin_ions = sumekinion

RETURN
END FUNCTION enerkin_ions

!-----energy_ions--------------------------------------------

REAL(DP) FUNCTION energ_ions()

!     Potential energy of ions

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)


INTEGER,EXTERNAL :: iptyp

IF(nion2 == 2) THEN
  energ_ions = 0D0
  RETURN
END IF

srenergy=0.0D0
sumion = 0.0D0
enii=0.0D0
enig=0.0D0
engg=0.0D0

IF(ipsptyp == 0) THEN
  IF(nion2 /= 0) THEN
    sumion = 0D0
    DO ion=2,nion
      DO ion1=1,ion-1
        dist2  = (cx(ion)-cx(ion1))**2+(cy(ion)-cy(ion1))**2  &
            +(cz(ion)-cz(ion1))**2
        dist   = SQRT(dist2)
        
        sumion = sumion + v_ion_ion(dist,ion,ion1)
      END DO
    END DO
  END IF
ELSE IF(ipsptyp >= 1) THEN
  IF(nion2 /= 0) THEN
    sumion = 0D0
    DO ion1=2,nion
      DO ion2=1,ion1-1
        dist2 = (cx(ion1)-cx(ion2))**2+(cy(ion1)-cy(ion2))**2  &
            +(cz(ion1)-cz(ion2))**2
        dist = SQRT(dist2)
        sumion = sumion + e2*ch(np(ion1))*ch(np(ion2))/dist
        
        
      END DO
    END DO
    
  END IF
ELSE
  STOP 'this sort of PsP not yet implemented'
END IF


IF(idielec /= 0) THEN
  
!     each ion interacts with each image ion!!
  
  IF(ipsptyp == 0) THEN
    IF(nion2 /= 0) THEN
      DO ion=1,nion
        DO ion1=1,nion
          dist2  = (cx(ion)+cx(ion1)-2D0*xdielec)**2 +(cy(ion)-cy(ion1))**2  &
              +(cz(ion)-cz(ion1))**2
          dist   = SQRT(dist2)
          sumion = sumion - (epsdi-1)/(epsdi+1)* v_ion_ion(dist,ion,ion1)
        END DO
      END DO
    END IF
  ELSE IF(ipsptyp >= 1) THEN
    IF(nion2 /= 0) THEN
      DO ion1=1,nion
        DO ion2=1,nion1
          dist2 = (cx(ion1)+cx(ion2)-2D0*xdielec)**2 +(cy(ion1)-cy(ion2))**2  &
              +(cz(ion1)-cz(ion2))**2
          dist = SQRT(dist2)
          sumion = sumion - (epsdi-1)/(epsdi+1)*  &
              e2*ch(np(ion1))*ch(np(ion2))/dist
          
          
        END DO
      END DO
      
    END IF
  ELSE
    STOP 'this sort of PsP not yet implemented'
  END IF
  
  
  
END IF


enii=sumion
            !WRITE(*,*) ' IONS: sumion=',sumion,enii


#if(raregas)
!     now the ion-GSM energy
IF (isurf /= 0) THEN
  IF (ipsptyp == 0) THEN
    
    DO ion=1,nion
      
      DO ion1=1,nc+NE+nk
        
        CALL getparas(ion1)
        
        dist2 = (rvectmp(1)-cx(ion))**2+(rvectmp(2)-cy(ion))**2  &
            +(rvectmp(3)-cz(ion))**2
        dist = SQRT(dist2)
        dist = MAX(dist,small)
        
        sss1 = sigsig(sigtmp,sgm1(np(ion)))
        sss2 = sigsig(sigtmp,sgm2(np(ion)))
        
                !WRITE(*,*) ' ion1 etc:',ion,ion1,dist, &
                !    e2*chgtmp*chg1(np(ion))* v_soft(dist,sss1*sq2), &
                !    e2*chgtmp*chg2(np(ion))* v_soft(dist,sss2*sq2),sumion

        sumion = sumion + e2*chgtmp*chg1(np(ion))* v_soft(dist,sss1*sq2)
        sumion = sumion + e2*chgtmp*chg2(np(ion))* v_soft(dist,sss2*sq2)

                  !WRITE(*,*) ' ion1 etc:',ion,ion1,dist, &
                !    e2*chgtmp*chg1(np(ion))* v_soft(dist,sss1*sq2), &
                !    e2*chgtmp*chg2(np(ion))* v_soft(dist,sss2*sq2),sumion

        IF(idielec /= 0) THEN
          
          dist2 = (-rvectmp(1)+2D0*xdielec-cx(ion))**2+  &
              (rvectmp(2)-cy(ion))**2+ (rvectmp(3)-cz(ion))**2
          dist = SQRT(dist2)
          dist = MAX(dist,small)
          
          sumion = sumion - (epsdi-1)/(epsdi+1)* e2*chgtmp*chg1(np(ion))*  &
              v_soft(dist,sss1*sq2)
          sumion = sumion - (epsdi-1)/(epsdi+1)* e2*chgtmp*chg2(np(ion))*  &
              v_soft(dist,sss2*sq2)
          
          
          
        END IF
        
        
        ion11=iconvlongtoshort(ion1)
        CALL getshortenergy(iptyp(ion1),4,ion11,ion)
        
             !WRITE(*,*) ' ion11:',ion11,rscaltmp
        sumion = sumion + rscaltmp
        
      END DO
      
      
    END DO
    
  ELSE
    
    STOP 'Goedecker still to be re-implemented with isurf.ne.1'
    
  END IF
  
  enig = sumion - enii
  
!     now the GSM-GSM
  
  
  DO ii = 1,nc+NE+nk-1
    
    ccc1 = getcharge(ii)
    
    
    DO jj = ii+1,nc+NE+nk
      
      dist2 = getdistance2(ii,jj)
      
      
      
      IF (ii <= nc .AND. jj == ii+nc) THEN ! spring force
        sumion = sumion + 0.5D0*cspr*dist2
         ! WRITE(*,*) ' cspr:',ii,jj,0.5D0*cspr*dist2
      ELSE ! Coulombian erf force
        
        dist = SQRT(dist2)
        
        
        ccc2 = getcharge(jj)
        sss = getmixedwidth(ii,jj)
        
        
!     electrostatic part
        sumion = sumion + ccc1*ccc2*e2*v_soft(dist,sss*sq2)
        
        
        ii2= iconvlongtoshort(ii)
        jj2= iconvlongtoshort(jj)
        
        CALL getshortenergy(iptyp(ii),iptyp(jj),ii2,jj2)
        
        sumion = sumion + rscaltmp
          !WRITE(*,*) ' iijj:',ii,jj,ccc1,ccc2,dist, &
          !         ccc1*ccc2*e2*v_soft(dist,sss*sq2),rscaltmp
        
      END IF
      
!               write(6,'(2a,3f15.8)') 'Achtung############sumion=',
!     &             'rScalTmp=    dist^2=',
!     &            sumion,rScalTmp,dist2
    END DO
  END DO

WRITE(*,*) ' GSM: sumion=',sumion,idielec
  
  IF(idielec /= 0) THEN
    
    DO ii = 1,nc+NE+nk
      
      CALL getparas(ii)
      
      ccc1 = chgtmp
      xtmp = rvectmp(1)
      ytmp = rvectmp(2)
      ztmp = rvectmp(3)
      
      DO jj = 1,nc+NE+nk
        
        CALL getparas(jj)
        
        dist2 = (xtmp+rvectmp(1)-2D0*xdielec)**2+ (ytmp-rvectmp(2))**2+  &
            (ztmp-rvectmp(3))**2
        
        dist = dist2**0.5D0
        
        
        ccc2 = chgtmp
        sss = getmixedwidth(ii,jj)
        
        
!     electrostatic part
        sumion = sumion - (epsdi-1D0)/(epsdi+1D0)*  &
            ccc1*ccc2*e2*v_soft(dist,sss*sq2)
        
      END DO
    END DO
  END IF
  
  
  
END IF ! isurf eq 0
#endif

engg = sumion - enii - enig

          !WRITE(*,*) ' sumion,enii,enig=',sumion,enii,enig

energ_ions = sumion

RETURN
END FUNCTION energ_ions
!------------------------------------------------------------


!-----itstep-------------------------------------------------

SUBROUTINE itstep(rho,it,psi)

!     ionic time step, propagates all heavy particles (not valence cloud)


USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                 :: rho(2*kdfull2)
INTEGER, INTENT(IN)                      :: it
COMPLEX(DP), INTENT(IN OUT)              :: psi(kdfull2,kstate)
REAL(DP),ALLOCATABLE :: xm(:)       ! array for actual particle mass
REAL(DP),ALLOCATABLE :: tfac(:)     ! array acceleration factors

LOGICAL,PARAMETER :: tsmooth = .true.  ! switch to smooth termination of acceleration
!INTEGER, PARAMETER :: mm=ng
!REAL(DP) :: px(ng),py(ng),pz(ng)


! what time is it ?
tfs = it*dt1*0.0484D0/(2D0*ame)
!WRITE(*,*) ' ITSTEP entered. tfs=',tfs

#if(raregas)
IF (isurf /= 0) THEN
  DO i=1,nc
    IF (imobc(i) == 0) THEN
      pxc(i)=0D0
      pyc(i)=0D0
      pzc(i)=0D0
    END IF
  END DO
!val         do i=1,ne
!val            if (imobe(i).eq.0) then
!val               pxe(i)=0.0
!val               pye(i)=0.0
!val               pze(i)=0.0
!val            endif
!val         enddo
  DO i=1,nk
    IF (imobk(i) == 0) THEN
      pxk(i)=0D0
      pyk(i)=0D0
      pzk(i)=0D0
    END IF
  END DO
END IF
#endif



!     PROPAGATION OF POSITIONS

#if(raregas)
IF (isurf /= 0 .AND. nc > 0) THEN
  IF(ALLOCATED(xm)) DEALLOCATE(xm)
  ALLOCATE(xm(1:nc))
  xm = mion*1836D0*ame
  CALL leapfr(xc(1),yc(1),zc(1),pxc(1),pyc(1),pzc(1),dt1,xm,nc,1)
  DEALLOCATE(xm)
END IF

IF (nk > 0) THEN
! correction of velocities to fox c.m.
  IF(ifixcmion == 1) THEN
    pxcmion = SUM(pxk(1:nion))/nion
    pycmion = SUM(pyk(1:nion))/nion
    pzcmion = SUM(pzk(1:nion))/nion
    pxk = pxk - pxcmion
    pyk = pyk - pycmion
    pzk = pzk - pzcmion
  END IF
  IF(ALLOCATED(xm)) DEALLOCATE(xm)
  ALLOCATE(xm(1:nk))
  xm = mkat*1836D0*ame
  CALL leapfr(xk(1),yk(1),zk(1),pxk(1),pyk(1),pzk(1),dt1,xm,nk,3)
  DEALLOCATE(xm)
END IF
#endif


IF (nion > 0 .AND.imob /= 0) THEN

! compute c.m. of ions
   IF(ifixcmion == 1) THEN
     xcmion = SUM(cx(1:nion))/nion
     ycmion = SUM(cy(1:nion))/nion
     zcmion = SUM(cz(1:nion))/nion
   END IF

!     propagation of cluster ions
  IF(ALLOCATED(xm)) DEALLOCATE(xm)
  ALLOCATE(xm(1:nion))
  DO i=1,nion
     xm(i)=amu(np(i))*1836D0*ame
  END DO
!  xm(:)=amu(np(:))*1836D0*ame
  CALL leapfr(cx(1),cy(1),cz(1), cpx(1),cpy(1),cpz(1),  &
      dt1,xm,nion,4)
  DEALLOCATE(xm)

! correct ionic c.m. to restore value before step
   IF(ifixcmion == 1) THEN
     xcmion = xcmion-SUM(cx(1:nion))/nion
     ycmion = ycmion-SUM(cy(1:nion))/nion
     zcmion = zcmion-SUM(cz(1:nion))/nion
     cx = cx+xcmion
     cy = cy+ycmion
     cz = cz+zcmion
   END IF

  
END IF


!MB      if(nclust.eq.0 .and. isurf.ne.0) call adjustdip(rho)

! for pure rare gas cluster in uniform translation,
! the cloud must be adjust before computation of forces



!     PROPAGATION OF MOMENTA


IF(ipseudo == 1)THEN
  CALL updatesubgrids
END IF


!mb      if(nc+ne+nk.gt.0 .and. isurf.ne.0)then

!     compute forces on rare gas cores with new positions
!WRITE(*,*) ' before GETFORCES, time=',time
CALL getforces(rho,psi,0)
  IF(tfs < taccel-1D-6) THEN
!    WRITE(*,*) 'FZ:',fz(1:nion)
    ALLOCATE(tfac(1:nion))
    tfac(:) = ame*amu(np(:))*1836D0*0.048D0/taccel
    IF(tsmooth) THEN
      dtaccel=taccel/10D0
      tfac = tfac*taccel/(taccel-dtaccel/2D0)
      IF(tfs>taccel-dtaccel) tfac=tfac*sin(PI*(taccel-tfs)/(dtaccel*2D0))**2
    END IF
    WRITE(*,*) ' tfs,tfac=',tfs,tfac
    fx(nion) = fx(nion)+vpx*tfac(nion)
    fy(nion) = fy(nion)+vpy*tfac(nion)
    fz(nion) = fz(nion)+vpz*tfac(nion)
    IF(nion>1) THEN
      fx(1:nion-1)=fx(1:nion-1)-vpx*tfac(1:nion-1)/(nion-1)
      fy(1:nion-1)=fy(1:nion-1)-vpy*tfac(1:nion-1)/(nion-1)
      fz(1:nion-1)=fz(1:nion-1)-vpz*tfac(1:nion-1)/(nion-1)
    END IF
!    WRITE(*,*) 'FZ:',fz(1:nion)
    DEALLOCATE(tfac)
  END IF


!mb      endif

!mb      call calcforce(rho,it,psi,3)
!     compute forces on cluster ions with new positions

#if(raregas)
IF (isurf /= 0) THEN
!     propagation of rare gas cores
  IF(nc > 0)THEN
    IF(ALLOCATED(xm)) DEALLOCATE(xm)
    ALLOCATE(xm(1:nc))
    xm = 1D0               ! setting for propagation of momenta
    CALL leapfr(pxc(1),pyc(1),pzc(1), fxc(1),fyc(1),fzc(1),dt1,xm,nc,1)
    DEALLOCATE(xm)
  END IF
  
  IF (nk > 0) THEN
    IF(ALLOCATED(xm)) DEALLOCATE(xm)
    ALLOCATE(xm(1:nk))
    xm = 1D0               ! setting for propagation of momenta
    CALL leapfr(pxk(1),pyk(1),pzk(1), fxk(1),fyk(1),fzk(1),dt1,xm,nk,3)
    DEALLOCATE(xm)
  END IF
END IF
#endif


!     propagation of cluster ions
IF (nion > 0 .AND. imob /= 0) THEN
  IF(ALLOCATED(xm)) DEALLOCATE(xm)
  ALLOCATE(xm(1:nion))
  xm = 1D0               ! setting for propagation of momenta
  CALL leapfr(cpx(1),cpy(1),cpz(1),fx(1),fy(1),fz(1),dt1,xm,nion,4)
  DEALLOCATE(xm)
END IF

!WRITE(6,'(a,3f15.6)') 'fionx,fiony,fionz',fx(1),fy(1),fz(1)
!WRITE(6,'(a,3f15.6)') 'pix,piy,piz',cpx(1),cpy(1),cpz(1)

RETURN
END SUBROUTINE itstep


!     ****************************************************

SUBROUTINE leapfr(x,y,z,xprop,yprop,zprop,ddt,xm,n,ityp)

!     ****************************************************

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

!INTEGER, PARAMETER :: mm=maxpar
INTEGER, INTENT(IN)                      :: n
REAL(DP), INTENT(OUT)                    :: x(n)
REAL(DP), INTENT(OUT)                    :: y(n)
REAL(DP), INTENT(OUT)                    :: z(n)
REAL(DP), INTENT(IN)                     :: xprop(n)
REAL(DP), INTENT(IN)                     :: yprop(n)
REAL(DP), INTENT(IN)                     :: zprop(n)
REAL(DP), INTENT(IN)                     :: ddt
REAL(DP), INTENT(IN)                     ::  xm(1:n)   
INTEGER, INTENT(IN)                      :: ityp



!   computes a leap frog step for n particles of masses xm
!   propagation of positions and momenta are separate
!   if x,y,z are positions, then xprop,... are momenta : ddt=dt, xm=mass
!   if x,y,z are momenta, then xprop,... are forces : ddt=dt, xm=1.


#if(raregas)
DO i = 1,n
  
  IF ((ityp == 4) .OR. (ityp == 1 .AND. imobc(i) == 1) .OR.  &
        (ityp == 2 .AND. imobe(i) == 1) .OR.  &
        (ityp == 3 .AND. imobk(i) == 1)) THEN
    
    x(i) = x(i) + ddt * xprop(i) / xm(i)
    y(i) = y(i) + ddt * yprop(i) / xm(i)
    z(i) = z(i) + ddt * zprop(i) / xm(i)
    
  END IF
END DO
#else
IF (ityp == 4) THEN
  DO i = 1,n
    x(i) = x(i) + xprop(i) * (ddt/xm(i))
    y(i) = y(i) + yprop(i) * (ddt/xm(i))
    z(i) = z(i) + zprop(i) * (ddt/xm(i))
  END DO
END IF
#endif


RETURN
END SUBROUTINE leapfr


!-----itstepv-------------------------------------------------

SUBROUTINE itstepv(rho,it,psi)

!  ionic time step, propagates all heavy particles (not valence cloud)
!  by velocity Verlet

USE params
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                 :: rho(2*kdfull2)
INTEGER, INTENT(IN)                      :: it
COMPLEX(DP), INTENT(IN OUT)              :: psi(kdfull2,kstate)

REAL(DP),ALLOCATABLE :: xm(:)          ! array for particle masses
REAL(DP),ALLOCATABLE :: tfac(:)        ! acceleration factors
LOGICAL,PARAMETER :: tsmooth = .true.  ! switch to smooth termination of acceleration


! what time is it ?
tfs = it*dt1*0.0484D0/(2D0*ame)
WRITE(*,*) 'enter ITSTEPV: tfs,cpx,cpy,cpz=',tfs,cpx(nion),cpy(nion),cpz(nion)

#if(raregas)
!                    zero in mobile atoms
IF (isurf /= 0) THEN
  DO i=1,nc
    IF (imobc(i) == 0) THEN
      pxc(i)=0D0
      pyc(i)=0D0
      pzc(i)=0D0
    END IF
  END DO
  DO i=1,nk
    IF (imobk(i) == 0) THEN
      pxk(i)=0D0
      pyk(i)=0D0
      pzk(i)=0D0
    END IF
  END DO
END IF
#endif



!     PROPAGATION OF POSITIONS

#if(raregas)
IF (isurf /= 0 .AND. nc > 0) THEN
  IF(ALLOCATED(xm)) DEALLOCATE(xm)
  ALLOCATE(xm(1:nc))
  xm = mion*1836D0*ame
  CALL velverlet1(xc(1),yc(1),zc(1),pxc(1),pyc(1),pzc(1), &
                  fxc(1),fyc(1),fzc(1),dt1,xm,nc,1)
  DEALLOCATE(xm)
END IF
IF (nk > 0) THEN
! correction of velocities to fix c.m.
   IF(ifixcmion == 1) THEN
     pxcmion = SUM(pxk(1:nk))/nk
     pycmion = SUM(pyk(1:nk))/nk
     pzcmion = SUM(pzk(1:nk))/nk
     pxk = pxk - pxcmion
     pyk = pyk - pycmion
     pzk = pzk - pzcmion
   END IF
  IF(ALLOCATED(xm)) DEALLOCATE(xm)
  ALLOCATE(xm(1:nk))
  xm = mkat*1836D0*ame
  CALL velverlet1(xk(1),yk(1),zk(1),pxk(1),pyk(1),pzk(1), &
                  fxk(1),fyk(1),fzk(1),dt1*modionstep,xm,nk,3)
  DEALLOCATE(xm)
END IF
#endif


IF (nion > 0 .AND.imob /= 0) THEN

! compute c.m. of ions
   IF(ifixcmion == 1) THEN
     xcmion = SUM(cx(1:nion))/nion
     ycmion = SUM(cy(1:nion))/nion
     zcmion = SUM(cz(1:nion))/nion
   END IF

!   propagation of cluster ions, first half of momenta
  IF(ALLOCATED(xm)) DEALLOCATE(xm)
  ALLOCATE(xm(1:nion))
  xm(:)=amu(np(:))*1836D0*ame
  CALL velverlet1(cx(1),cy(1),cz(1), cpx(1),cpy(1),cpz(1), &
                  fx(1),fy(1),fz(1),dt1*modionstep,xm,nion,4)
  DEALLOCATE(xm)

! correct ionic c.m. to restore value before step
   IF(ifixcmion == 1) THEN
     xcmion = xcmion-SUM(cx(1:nion))/nion
     ycmion = ycmion-SUM(cy(1:nion))/nion
     zcmion = zcmion-SUM(cz(1:nion))/nion
     cx = cx+xcmion
     cy = cy+ycmion
     cz = cz+zcmion
   END IF

  
END IF


!MB      if(nclust.eq.0 .and. isurf.ne.0) call adjustdip(rho)

! for pure rare gas cluster in uniform translation,
! the cloud must be adjust before computation of forces


!     update forces


IF(ipseudo == 1)THEN
  CALL updatesubgrids
END IF


!     compute forces on rare gas cores with new positions
CALL getforces(rho,psi,0)
!WRITE(*,*) ' ITSTEPV: after GETFORCES'
  IF(tfs < taccel-1D-6) THEN
!    WRITE(*,*) 'FZ:',fz(1:nion)
    ALLOCATE(tfac(1:nion))
    tfac(:) = ame*amu(np(:))*1836D0*0.048D0/taccel
    IF(tsmooth) THEN
      dtaccel=taccel/10D0
      tfac = tfac*taccel/(taccel-dtaccel/2D0)
      IF(tfs>taccel-dtaccel) tfac=tfac*sin(PI*(taccel-tfs)/(dtaccel*2D0))**2
    END IF
    WRITE(*,*) ' tfs,tfac=',tfs,tfac
    fx(nion) = vpx*tfac(nion)
    fy(nion) = vpy*tfac(nion)
    fz(nion) = vpz*tfac(nion)
    IF(nion>1) THEN
      fx(1:nion-1)=-vpx*tfac(1:nion-1)/(nion-1)
      fy(1:nion-1)=-vpy*tfac(1:nion-1)/(nion-1)
      fz(1:nion-1)=-vpz*tfac(1:nion-1)/(nion-1)
    END IF
    DEALLOCATE(tfac)
!    WRITE(*,*) 'FZ:',fz(1:nion)
  END IF


!     second half of propagation of momenta


#if(raregas)
IF (isurf /= 0) THEN
!     propagation of rare gas cores
  IF(nc > 0)THEN
    CALL velverlet2(pxc(1),pyc(1),pzc(1), &
                    fxc(1),fyc(1),fzc(1),dt1*modionstep,nc,1)
  END IF
  
  IF (nk > 0) THEN
    CALL velverlet2(pxk(1),pyk(1),pzk(1), &
                    fxk(1),fyk(1),fzk(1),dt1*modionstep,nc,3)
  END IF
END IF
#endif


!     propagation of cluster ions
IF (nion > 0 .AND. imob /= 0) CALL velverlet2(cpx(1),cpy(1),cpz(1),  &
    fx(1),fy(1),fz(1), dt1*modionstep,nion,4)

WRITE(6,'(a,3f15.6)') 'fionx,fiony,fionz',fx(1),fy(1),fz(1)
WRITE(6,'(a,3f15.6)') 'pix,piy,piz',cpx(1),cpy(1),cpz(1)

RETURN
END SUBROUTINE itstepv



SUBROUTINE velverlet1(x,y,z,px,py,pz,fox,foy,foz,ddt,xm,n,ityp)

!     ****************************************************

! first part of ionic time step with velocity-Verlet
!  x,y,z         = arrays of positions
!  px,py,pz      = arrays of momenta  
!  fox,foy,foz   = arrays of forces, before the step
!  ddt           = size of time step
!  xm            = ionic mass
!  n             = number of ions
!  ityp          = type of ions/atoms
!
!  note that momenta are only stepped by a half step,
!  waiting for a second half step after computing the new
!  forces

USE params
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)        :: n
REAL(DP), INTENT(IN OUT)   :: x(n),y(n),z(n)
REAL(DP), INTENT(IN OUT)   :: px(n),py(n),pz(n)
REAL(DP), INTENT(IN)       :: fox(n),foy(n),foz(n)
REAL(DP), INTENT(IN)       :: ddt
REAL(DP), INTENT(IN)       :: xm(1:n)
INTEGER, INTENT(IN)        :: ityp


ddth = 0.5D0*ddt
#if(raregas)
DO i = 1,n
  
  IF ((ityp == 1 .AND. imobc(i) == 1) .OR.  &
        (ityp == 2 .AND. imobe(i) == 1) .OR.  &
        (ityp == 3 .AND. imobk(i) == 1) .OR. (ityp == 4)           ) THEN
    
    x(i) = x(i) + (px(i)+ddth*fox(i))*ddt/xm(i)
    y(i) = y(i) + (py(i)+ddth*foy(i))*ddt/xm(i)
    z(i) = z(i) + (pz(i)+ddth*foz(i))*ddt/xm(i)
    px(i) = px(i) + ddth*fox(i)
    py(i) = py(i) + ddth*foy(i)
    pz(i) = pz(i) + ddth*foz(i)
    
  END IF
END DO
#else
IF (ityp == 4) THEN
  DO i = 1,n
    x(i) = x(i) + (px(i)+ddth*fox(i))*(ddt/xm(i))
    y(i) = y(i) + (py(i)+ddth*foy(i))*(ddt/xm(i))
    z(i) = z(i) + (pz(i)+ddth*foz(i))*(ddt/xm(i))
    px(i) = px(i) + ddth*fox(i)
    py(i) = py(i) + ddth*foy(i)
    pz(i) = pz(i) + ddth*foz(i)
  END DO
END IF
#endif

RETURN
END SUBROUTINE velverlet1

SUBROUTINE velverlet2(px,py,pz,fox,foy,foz,ddt,n,ityp)

!     ****************************************************

! second part of ionic time step with velocity-Verlet
!  px,py,pz      = arrays of momenta  
!  fox,foy,foz   = arrays of forces, after the step
!  ddt           = size of time step
!  xm            = ionic mass
!  n             = number of ions
!  ityp          = type of ions/atoms

USE params
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)        :: n
REAL(DP), INTENT(IN OUT)   :: px(n),py(n),pz(n)
REAL(DP), INTENT(IN)       :: fox(n),foy(n),foz(n)
REAL(DP), INTENT(IN)       :: ddt
INTEGER, INTENT(IN)        :: ityp


ddth = 0.5D0*ddt

#if(raregas)
DO i = 1,n
  
  IF ((ityp == 1 .AND. imobc(i) == 1) .OR.  &
        (ityp == 2 .AND. imobe(i) == 1) .OR.  &
        (ityp == 3 .AND. imobk(i) == 1) .OR. (ityp == 4)           ) THEN
    
    px(i) = px(i) + ddth*fox(i)
    py(i) = py(i) + ddth*foy(i)
    pz(i) = pz(i) + ddth*foz(i)
    
  END IF
END DO
#else
IF (ityp == 4) THEN
  DO i = 1,n
    px(i) = px(i) + ddth*fox(i)
    py(i) = py(i) + ddth*foy(i)
    pz(i) = pz(i) + ddth*foz(i)
  END DO
END IF
#endif

RETURN
END SUBROUTINE velverlet2


!     *******************************

SUBROUTINE lffirststep(rho,psi)

!     *******************************
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)

REAL(DP),ALLOCATABLE :: xm(:)

tfs = 0D0

IF(ipseudo == 1)THEN
  CALL updatesubgrids
END IF




!GB      if(nc+ne+nk.gt.0 .and. isurf.ne.0)then
CALL getforces(rho,psi,0)
!GB      endif

#if(raregas)
IF (isurf /= 0) THEN
  IF(nc+NE > 0)THEN
!     propagation of cores
    IF(ALLOCATED(xm)) DEALLOCATE(xm)
    ALLOCATE(xm(1:nc))
    xm = 1D0
    CALL leapfr(pxc(1),pyc(1),pzc(1), fxc(1),fyc(1),fzc(1),dt1/2D0,xm,nc,1)
    DEALLOCATE(xm)
!     propagation of clouds
    IF(ifadiadip /= 1) THEN
      ALLOCATE(xm(1:ne))    ! possible buf with gfortran = debian64
      xm = 1D0
      CALL leapfr(pxe(1),pye(1),pze(1), fxe(1),fye(1),fze(1),dt1/4D0,xm,ne,2)
      DEALLOCATE(xm)
    END IF
  END IF
!     propagation of cation
  IF(nk > 0) THEN
    IF(ALLOCATED(xm)) DEALLOCATE(xm)
    ALLOCATE(xm(1:nk))
    xm = 1D0
    CALL leapfr(pxk(1),pyk(1),pzk(1),fxk(1),fyk(1),fzk(1),dt1/2D0,xm,nk,3)
    DEALLOCATE(xm)
  END IF
  
END IF
#endif


!     propagation of cluster ions
IF(nion > 0) THEN
  IF(ALLOCATED(xm)) DEALLOCATE(xm)
  ALLOCATE(xm(1:nion))
  xm = 1D0
  CALL leapfr(cpx(1),cpy(1),cpz(1),  &
              fx(1),fy(1),fz(1),dt1/2.,xm,nion,4)
  DEALLOCATE(xm)
END IF

WRITE(6,*)'initial momenta and kinetic energy'
ekion =  enerkin_ions()

RETURN
END SUBROUTINE lffirststep

!GB



!     ******************************

SUBROUTINE spheric(r,x,y,z,n)
USE params, ONLY: DP
IMPLICIT REAL(DP) (A-H,O-Z)

!     ******************************

INTEGER, PARAMETER :: mm=50000
REAL(DP), INTENT(IN)                         :: r
REAL(DP), INTENT(OUT)                        :: x(mm)
REAL(DP), INTENT(OUT)                        :: y(mm)
REAL(DP), INTENT(OUT)                        :: z(mm)
INTEGER, INTENT(IN)                      :: n


IF(n > mm) THEN
  WRITE(6,*)'warning !! number of particles exceeds array s size'
  WRITE(6,*)'in routine spheric'
  RETURN
END IF

!   samples a spherical distribution of n particles  on
!   a sphere of radius r. results are x(i), y(i),z(i)

DATA im,ia,ic/259200,7141,54773/
DATA pi/3.1415927D0/
jran = 12345

WRITE(6,*)r,n
DO i = 1,n
  rr = r
  jran = MOD(jran*ia + ic,im)
  xx = REAL(jran,DP) / REAL(im,DP)
  cth = 1D0 - 2D0 * xx
  sth = SQRT ( 1D0 - cth * cth )
  jran = MOD(jran*ia + ic,im)
  xx = REAL(jran,DP) / REAL(im,DP)
  phi = 2D0 * pi * xx
  cosp = COS (phi)
  sinp = SIN (phi)
  x(i) = rr * sth * cosp
  y(i) = rr * sth * sinp
  z(i) = rr *cth
END DO
RETURN
END SUBROUTINE spheric

!     ************************************

SUBROUTINE conslw(x,y,z,px,py,pz,n)
USE params, ONLY: DP
IMPLICIT REAL(DP) (A-H,O-Z)

!     ************************************

!INTEGER, PARAMETER :: mm=50000
INTEGER, INTENT(IN)                      :: n
REAL(DP), INTENT(IN OUT)                     :: x(n)
REAL(DP), INTENT(IN OUT)                     :: y(n)
REAL(DP), INTENT(IN OUT)                     :: z(n)
REAL(DP), INTENT(IN OUT)                     :: px(n)
REAL(DP), INTENT(IN OUT)                     :: py(n)
REAL(DP), INTENT(IN OUT)                     :: pz(n)


!    renormalization of c.o.m, of momentum and of angular momentum

xcm = 0D0
ycm = 0D0
zcm = 0D0
xkcm = 0D0
ykcm = 0D0
zkcm = 0D0
rlx = 0D0
rly = 0D0
rlz = 0D0
rixx = 0D0
riyy = 0D0
rizz = 0D0
rixy = 0D0
rixz = 0D0
riyz = 0D0
npart = n
reno = 1D0/REAL(n,DP)

DO i =1,n
  rx = x(i)
  ry = y(i)
  rz = z(i)
  rkx = px(i)
  rky = py(i)
  rkz = pz(i)
  rx2 = rx * rx
  ry2 = ry * ry
  rz2 = rz * rz
  
!         c . o . m. in r space
  
  xcm = xcm + rx
  ycm = ycm + ry
  zcm = zcm + rz
  
!         c . o . m. in p space
  
  xkcm = xkcm + rkx
  ykcm = ykcm + rky
  zkcm = zkcm + rkz
  
!         angular momentum
  
  rlx = rlx + ry * rkz - rz * rky
  rly = rly + rz * rkx - rx * rkz
  rlz = rlz + rx * rky - ry * rkx
  
!         tensor of inertia
  
  rixx = rixx + ry2 + rz2
  riyy = riyy + rx2 + rz2
  rizz = rizz + rx2 + ry2
  rixy = rixy - rx * ry
  rixz = rixz - rx * rz
  riyz = riyz - ry * rz
END DO

!        renormalization of c . o . m . in r and p space
!               and of the total angular momentum

xcm = xcm * reno
ycm = ycm * reno
zcm = zcm * reno
xkcm = xkcm * reno
ykcm = ykcm * reno
zkcm = zkcm * reno

rixx = rixx * reno - ycm * ycm - zcm * zcm
riyy = riyy * reno - xcm * xcm - zcm * zcm
rizz = rizz * reno - xcm * xcm - ycm * ycm
rixy = rixy * reno + xcm * ycm
rixz = rixz * reno + xcm * zcm
riyz = riyz * reno + ycm * zcm

rlxcm = ycm * zkcm - zcm * ykcm
rlycm = zcm * xkcm - xcm * zkcm
rlzcm = xcm * ykcm - ycm * xkcm

rlx = rlx * reno - rlxcm
rly = rly * reno - rlycm
rlz = rlz * reno - rlzcm

det = rixx * riyy * rizz + 2. * rixy * rixz * riyz  &
    - rixx * riyz ** 2 - riyy * rixz ** 2 - rizz * rixy **2
wx  = (rlx * (riyy * rizz - riyz ** 2) + rly * (rixz * riyz - rizz * rixy)  &
    + rlz * (rixy * riyz - riyy * rixz)) / det
wy  = (rlx * (rixz * riyz - rixy * rizz) + rly * (rixx * rizz - rixz ** 2)  &
    + rlz * (rixy * rixz - rixx * riyz)) / det
wz  = (rlx * (rixy * riyz - riyy * rixz)  &
    + rly * (rixz * rixy - rixx * riyz) + rlz * (rixx * riyy - rixy ** 2)) / det

WRITE(6,*)' before...'
WRITE(6,*)' r c.m....',xcm,ycm,zcm
WRITE(6,*)' k c.m....',xkcm,ykcm,zkcm
WRITE(6,*)' l........',rlx,rly,rlz

xxcm = 0D0
yycm = 0D0
zzcm = 0D0
xxkcm = 0D0
yykcm = 0D0
zzkcm = 0D0
rlx = 0D0
rly = 0D0
rlz = 0D0

DO i = 1,n
  x(i) = x(i) - xcm
  y(i) = y(i) - ycm
  z(i) = z(i) - zcm
  xxcm = xxcm + x(i)
  yycm = yycm + y(i)
  zzcm = zzcm + z(i)
  px(i) = px(i) - xkcm - (wy * z(i) - wz * y(i))
  py(i) = py(i) - ykcm - (wz * x(i) - wx * z(i))
  pz(i) = pz(i) - zkcm - (wx * y(i) - wy * x(i))
  xxkcm = xxkcm + px(i)
  yykcm = yykcm + py(i)
  zzkcm = zzkcm + pz(i)
  rx = x(i)
  ry = y(i)
  rz = z(i)
  rkx = px(i)
  rky = py(i)
  rkz = pz(i)
  rlx = rlx + ry * rkz - rz * rky
  rly = rly + rz * rkx - rx * rkz
  rlz = rlz + rx * rky - ry * rkx
  rixx = rixx + ry * ry + rz * rz
  rixy = rixy + rx * ry
  rixz = rixz + rx * rz
  riyy = riyy + rx * rx + rz * rz
  riyz = riyz + ry * rz
  rizz = rizz + rx * rx + ry * ry
END DO

xxcm = xxcm * reno
yycm = yycm * reno
zzcm = zzcm * reno
xxkcm = xxkcm * reno
yykcm = yykcm * reno
zzkcm = zzkcm * reno
rlx = rlx * reno
rly = rly * reno
rlz = rlz * reno

WRITE(6,*)' after....'
WRITE(6,*)' r c.m....',xxcm,yycm,zzcm
WRITE(6,*)' k c.m....',xxkcm,yykcm,zzkcm
WRITE(6,*)' l .......',rlx,rly,rlz

!        end of the renormalizations ( r , p , l )


RETURN
END SUBROUTINE conslw
#if(gridfft)

!     *******************

SUBROUTINE stream(rho,psi)

!     *******************

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN)                         :: rho(2*kdfull2)
COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)


REAL(DP), ALLOCATABLE :: rhos(:),drho(:)

ALLOCATE(rhos(kdfull2),drho(kdfull2))


query2=qe(2)*qold2
qold2=qe(2)
query3=qe(3)*qold3
qold3=qe(3)
query4=qe(4)*qold4
qold4=qe(4)
!test
query2=1.0D0
query3=1.0D0
!test
IF(tfs == 0D0) THEN
  iquery4=0
  WRITE(6,*) 'iquery4=',iquery4
  DO i=1,nxyz
    rhos(i)=0D0
  END DO
END IF
IF(tfs > 0D0) THEN
  IF(query2 < 0D0) THEN
    iquery4=iquery4+1
    WRITE(6,*) 'iquery4=',iquery4
    IF(MOD(iquery4,10) == 0) THEN
      DO i=1,nxyz
        drho(i)=rho(i)-rhos(i)
      END DO
      CALL write_density(rhos,'sdensi')
      CALL write_density(drho,'sderiv')
!           stop 'calculated a derivative of density after x-stimulation'
    END IF
  END IF
  IF(query3 < 0D0) THEN
    iquery4=iquery4+1
    WRITE(6,*) 'iquery4=',iquery4
    IF(MOD(iquery4,10) == 0) THEN
      DO i=1,nxyz
        drho(i)=rho(i)-rhos(i)
      END DO
      CALL write_density(rhos,'sdensi')
      CALL write_density(drho,'sderiv')
!           stop 'calculated a derivative of density after y-stimulation'
    END IF
  END IF
  IF(query4 < 0D0) THEN
    iquery4=iquery4+1
    WRITE(6,*) 'iquery4=',iquery4
    IF(MOD(iquery4,10) == 0) THEN
      DO i=1,nxyz
        drho(i)=rho(i)-rhos(i)
      END DO
      CALL write_density(rhos,'sdensi')
      CALL write_density(drho,'sderiv')
!           stop 'calculated a derivative of density after z-stimulation'
    END IF
  END IF
END IF
rhos = rho

DEALLOCATE(rhos,drho)

RETURN
END SUBROUTINE stream
! ---eltherm----------------------------------------------eltherm-------

SUBROUTINE eltherm(rho,psi,expjx,expjy,expjz,eeth,iask)

!  computes electronic thermal energy via the expectation value of the
!  the current operator to check electronic thermalization

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN)                         :: rho(2*kdfull2)
COMPLEX(DP), INTENT(IN)                      :: psi(kdfull2,kstate)
REAL(DP), INTENT(OUT)                        :: expjx
REAL(DP), INTENT(OUT)                        :: expjy
REAL(DP), INTENT(OUT)                        :: expjz
REAL(DP), INTENT(OUT)                        :: eeth
INTEGER, INTENT(IN OUT)                  :: iask

COMPLEX(DP),DIMENSION(:),ALLOCATABLE :: p,q2
REAL(DP),DIMENSION(:),ALLOCATABLE :: ajtx,ajty,ajtz
!REAL(DP),DIMENSION(:),ALLOCATABLE :: akkx,akky,akkz
!DIMENSION akkx(kdfull2),akky(kdfull2),akkz(kdfull2)
!EQUIVALENCE (akkx(1),w1(1))
!EQUIVALENCE (akky(1),w2(1))
!EQUIVALENCE (akkz(1),w3(1))

REAL(DP) :: exjx(ksttot),exjy(ksttot),exjz(ksttot)
COMPLEX(DP) :: test

CHARACTER (LEN=6) :: ext
#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
REAL(DP), ALLOCATABLE :: aj(:)
ALLOCATE(aj(kdfull2))
!EQUIVALENCE (aj(1),q2(1))
#endif

!ALLOCATE(akkx(kdfull2))
!ALLOCATE(akky(kdfull2))
!ALLOCATE(akkz(kdfull2))
ALLOCATE(p(kdfull2),q2(kdfull2))
ALLOCATE(ajtx(kdfull2),ajty(kdfull2),ajtz(kdfull2))

!   init derivative

dkx=pi/(dx*REAL(nx,DP))
dky=pi/(dy*REAL(ny,DP))
dkz=pi/(dz*REAL(nz,DP))

!ind=0
!DO i3=1,nz2
!  IF(i3 >= (nz+1)) THEN
!    zkz=(i3-nz2-1)*dkz
!  ELSE
!    zkz=(i3-1)*dkz
!  END IF
!  DO i2=1,ny2
!    IF(i2 >= (ny+1)) THEN
!      zky=(i2-ny2-1)*dky
!    ELSE
!      zky=(i2-1)*dky
!    END IF
!    DO i1=1,nx2
!      IF(i1 >= (nx+1)) THEN
!        zkx=(i1-nx2-1)*dkx
!      ELSE
!        zkx=(i1-1)*dkx
!      END IF
!      ind=ind+1
!      akkx(ind)=-zkx
!      akky(ind)=-zky
!      akkz(ind)=-zkz
!    END DO
!  END DO
!END DO

!   compute the currents for each direction: ajtx/y/z
!   and their expectation values: exjx/y/z

q1=0D0
tel=0D0
DO i=1,kdfull2
  ajtx(i)=0D0
  ajty(i)=0D0
  ajtz(i)=0D0
END DO

DO nb=1,nstate
  
  exjx(nb)=0D0
  exjy(nb)=0D0
  exjz(nb)=0D0
  
#if(netlib_fft|fftw_cpu)
  CALL fftf(psi(1,nb),q2)

  DO ind=1,kdfull2
    q2(ind)=q2(ind)*akx(ind)
  END DO

  CALL fftback(q2,p)
#endif
#if(fftw_gpu)
  CALL fftf(psi(1,nb),q2,ffta,gpu_ffta)

  CALL multiply_ak2(gpu_ffta,gpu_akxfft,kdfull2)

  CALL fftback(q2,p,ffta,gpu_ffta)
#endif

  DO ind=1,kdfull2
    test=eye/2D0*(CONJG(psi(ind,nb))*p(ind) -psi(ind,nb)*CONJG(p(ind)))
    ajalpha=test
    exjx(nb)=exjx(nb)+ajalpha*CONJG(psi(ind,nb))*psi(ind,nb)
    ajtx(ind)=ajtx(ind)+ajalpha
  END DO
  
#if(netlib_fft|fftw_cpu)
  CALL fftf(psi(1,nb),q2)

  DO ind=1,kdfull2
    q2(ind)=q2(ind)*aky(ind)
  END DO

  CALL fftback(q2,p)
#endif

#if(fftw_gpu)
  CALL fftf(psi(1,nb),q2,ffta,gpu_ffta)

  CALL multiply_ak2(gpu_ffta,gpu_akyfft,kdfull2)

  CALL fftback(q2,p,ffta,gpu_ffta)
#endif

  DO ind=1,kdfull2
    test=eye/2D0*(CONJG(psi(ind,nb))*p(ind) -psi(ind,nb)*CONJG(p(ind)))
    ajalpha=test
    exjy(nb)=exjy(nb)+ajalpha*CONJG(psi(ind,nb))*psi(ind,nb)
    ajty(ind)=ajty(ind)+ajalpha
  END DO
  
#if(netlib_fft|fftw_cpu)
  CALL fftf(psi(1,nb),q2)

  DO ind=1,kdfull2
    q2(ind)=q2(ind)*akz(ind)
  END DO

  CALL fftback(q2,p)
#endif

#if(fftw_gpu)
  CALL fftf(psi(1,nb),q2,ffta,gpu_ffta)

  CALL multiply_ak2(gpu_ffta,gpu_akzfft,kdfull2)

  CALL fftback(q2,p,ffta,gpu_ffta)
#endif

  DO ind=1,kdfull2
    test=eye/2D0*(CONJG(psi(ind,nb))*p(ind) -psi(ind,nb)*CONJG(p(ind)))
    ajalpha=test
    exjz(nb)=exjz(nb)+ajalpha*CONJG(psi(ind,nb))*psi(ind,nb)
    ajtz(ind)=ajtz(ind)+ajalpha
  END DO
  
END DO                     !loop over states

!  compute sum of current expectation values in either direction:

expjx=0D0
expjy=0D0
expjz=0D0
DO n=1,nstate
  expjx=expjx+exjx(n)
  expjy=expjy+exjy(n)
  expjz=expjz+exjz(n)
END DO
expjx=expjx*dvol
expjy=expjy*dvol
expjz=expjz*dvol
#if(parayes)
DO ind=1,kdfull2
  aj(ind)=ajtx(ind)
END DO
CALL pi_allreduce(aj,ajtx,kdfull2,mpi_double_precision,mpi_sum,  &
    mpi_comm_world,ic)
DO ind=1,kdfull2
  aj(ind)=ajty(ind)
END DO
CALL pi_allreduce(aj,ajty,kdfull2,mpi_double_precision,mpi_sum,  &
    mpi_comm_world,ic)
DO ind=1,kdfull2
  aj(ind)=ajtz(ind)
END DO
CALL pi_allreduce(aj,ajtz,kdfull2,mpi_double_precision,mpi_sum,  &
    mpi_comm_world,ic)

CALL pi_allreduce(expjx,ex,1,mpi_double_precision,mpi_sum, mpi_comm_world,ic)
CALL pi_allreduce(expjy,ey,1,mpi_double_precision,mpi_sum, mpi_comm_world,ic)
CALL pi_allreduce(expjz,ez,1,mpi_double_precision,mpi_sum, mpi_comm_world,ic)
expjx=ex
expjy=ey
expjz=ez
#endif

!  compute average velocities from currents:

DO ind=1,kdfull2
  ajtx(ind)=ajtx(ind)/rho(ind)
  ajty(ind)=ajty(ind)/rho(ind)
  ajtz(ind)=ajtz(ind)/rho(ind)
END DO


! finally compute thermal energy:
! e_th = sum_i int[d^3r 0.5D0*rho_i*(v_av - v_i)**2]

tel=0D0
DO nb=1,nstate
  
#if(netlib_fft|fftw_cpu)
  CALL fftf(psi(1,nb),q2)

  DO ind=1,kdfull2
    q2(ind)=q2(ind)*akx(ind)
  END DO

  CALL fftback(q2,p)
#endif

#if(fftw_gpu)
  CALL fftf(psi(1,nb),q2,ffta,gpu_ffta)

  CALL multiply_ak2(gpu_ffta,gpu_akxfft,kdfull2)

  CALL fftback(q2,p,ffta,gpu_ffta)
#endif

  DO ind=1,kdfull2
    test=eye/2D0*(CONJG(psi(ind,nb))*p(ind) -psi(ind,nb)*CONJG(p(ind)))
    ajalpha=test
    rhoalpha=CONJG(psi(ind,nb))*psi(ind,nb)
    ajalpha=ajalpha/rhoalpha
    tel=tel+0.5D0*rhoalpha*(ajalpha-ajtx(ind))**2
  END DO
  
#if(netlib_fft|fftw_cpu)
  CALL fftf(psi(1,nb),q2)

  DO ind=1,kdfull2
    q2(ind)=q2(ind)*aky(ind)
  END DO

  CALL fftback(q2,p)
#endif

#if(fftw_gpu)
  CALL fftf(psi(1,nb),q2,ffta,gpu_ffta)

  CALL multiply_ak2(gpu_ffta,gpu_akyfft,kdfull2)

  CALL fftback(q2,p,ffta,gpu_ffta)
#endif

  DO ind=1,kdfull2
    test=eye/2D0*(CONJG(psi(ind,nb))*p(ind) -psi(ind,nb)*CONJG(p(ind)))
    ajalpha=test
    rhoalpha=CONJG(psi(ind,nb))*psi(ind,nb)
    ajalpha=ajalpha/rhoalpha
    tel=tel+0.5D0*rhoalpha*(ajalpha-ajty(ind))**2
  END DO
  
#if(netlib_fft|fftw_cpu)
  CALL fftf(psi(1,nb),q2)

  DO ind=1,kdfull2
    q2(ind)=q2(ind)*akz(ind)
  END DO

  CALL fftback(q2,p)
#endif

#if(fftw_gpu)
  CALL fftf(psi(1,nb),q2,ffta,gpu_ffta)

  CALL multiply_ak2(gpu_ffta,gpu_akzfft,kdfull2)

  CALL fftback(q2,p,ffta,gpu_ffta)
#endif

  DO ind=1,kdfull2
    test=eye/2D0*(CONJG(psi(ind,nb))*p(ind) -psi(ind,nb)*CONJG(p(ind)))
    ajalpha=test
    rhoalpha=CONJG(psi(ind,nb))*psi(ind,nb)
    ajalpha=ajalpha/rhoalpha
    tel=tel+0.5D0*rhoalpha*(ajalpha-ajtz(ind))**2
  END DO
  
END DO                     !loop over states

eeth=tel*dvol*0.5D0*hbar*hbar/2D0/ame !atomic units (h/2m)**2
#if(parayes)
tel=eeth
CALL pi_allreduce(tel,eeth,1,mpi_double_precision, mpi_sum,mpi_comm_world,ic)
#endif

WRITE(6,'(a,f12.5)') ' electronic thermal energy: ',eeth
!      write(17,'(a,f12.5)') 'electronic thermal energy: ',eeth

!DEALLOCATE(akkx)
!DEALLOCATE(akky)
!DEALLOCATE(akkz)
DEALLOCATE(p,q2)
DEALLOCATE(ajtx,ajty,ajtz)
#if(parayes)
DEALLOCATE(aj)
#endif

RETURN
END SUBROUTINE eltherm

!     **************************

SUBROUTINE write_density(drho,filename)
!  write density or derivative of density
!
!     **************************
USE params
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: drho(kdfull2)
CHARACTER (LEN=6), INTENT(IN) :: filename
CHARACTER (LEN=4) :: ext

WRITE(ext,'(a,i3.3)') ".", iquery4/10    ! e.g:  iquery4 = 1234  --->  ext = ".123"
                                         !       iquery4 = 568   --->  ext = ".056"
OPEN(UNIT=60,STATUS='unknown',FORM='unformatted',FILE= filename//'1'//ext)


DO i=1,nxyz
  WRITE(60) drho(i)
END DO

CLOSE(UNIT=60,STATUS='keep')

RETURN
END SUBROUTINE write_density


#endif
