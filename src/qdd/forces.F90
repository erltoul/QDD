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

 


!------------------------------------------------------------

SUBROUTINE getforces(rho,psi,it,iflag)
!------------------------------------------------------------
USE params
IMPLICIT NONE
!     calculates forces on all particles except DFT-electrons

!     GSM means a particle described by the Gaussian Shell Model, 
!     i.e. cores, valence shells and kations.

!     - calculation of the GSM-GSM forces
!     - calculation of the forces between GSM particles and
!       cluster cores, labelled Na
!     - calculation of the forces *between* cluster cores
!     - calculation of the force on GSM and cluster cores due
!       to the electronic density
!     - force from external potentials (laser, Madelung, etc.)

!     iflag   --->   0       calculate forces on all particles
!                    1                 forces on valence shells only
!                    2                 forces on valence shells only, but
!                                          skipping the spring force
!                                          (used in adjustdip)




REAL(DP), INTENT(IN OUT)   :: rho(2*kdfull2)
COMPLEX(DP), INTENT(IN)    :: psi(kdfull2,kstate)
INTEGER, INTENT(IN)        :: it
INTEGER, INTENT(IN)        :: iflag

INTEGER :: i
#if(raregas)
INTEGER :: ii
#endif

WRITE(6,*) 'Entering getforces: ipsptyp=',ipsptyp


!     In case of Goedecker PsP, switch to the corresponding routine
!     in the old force.F file

IF(ipsptyp == 1 .AND. tnonlocany) THEN
  CALL calcf_goenonl(rho,it,psi)
  CALL laserf(rho)
  CALL forceproject()
  OPEN(772,FILE='forces.'//outnam)
  DO i=1,nion
    WRITE(772,*) fx(i)
    WRITE(772,*) fy(i)
    WRITE(772,*) fz(i)
  ENDDO
  CALL flush(772)
  CLOSE(772)
  RETURN
ELSE IF((ipsptyp==1  .AND. .NOT.tnonlocany) .OR. ipsptyp == 2) THEN
  CALL calcf_goeloc(rho,it)
  CALL forceproject()
  OPEN(772,FILE='forces.'//outnam)
  DO i=1,nion
    WRITE(772,*) fx(i)
    WRITE(772,*) fy(i)
    WRITE(772,*) fz(i)
  ENDDO
  CALL FLUSH(772)
  CLOSE(772)
  RETURN
END IF

IF (iflag == 0) THEN
  fx(1:nion)=0D0;fy(1:nion)=0D0;fz(1:nion)=0D0
#if(raregas)
  DO ii=1,nc
    fxc(ii)=0.0D0
    fyc(ii)=0.0D0
    fzc(ii)=0.0D0
  END DO
  DO ii=1,nk
    fxk(ii)=0.0D0
    fyk(ii)=0.0D0
    fzk(ii)=0.0D0
  END DO
#endif
END IF

#if(raregas)
DO ii=1,NE
  fxe(ii)=0.0D0
  fye(ii)=0.0D0
  fze(ii)=0.0D0
END DO
#endif

!     pure cluster part

IF (iflag == 0) THEN
  
  IF (nclust > 0) CALL getforceelna(rho)
  
  CALL getforcenana()
  
END IF



!     GSM relevant parts
#if(raregas)
IF (isurf /= 0) THEN
  
!  IF (iswforce == 0) THEN
    CALL getforcegsmgsm(iflag)
!  ELSE
!    CALL getforcegsmgsm1(iflag)                !  ???
!  END IF
  
  CALL getforcenagsm(iflag)
  
  IF (nclust > 0) CALL getforceelgsm(iflag,rho)
  
  
  
!     van der Waals part
  IF (ivdw == 1) THEN
    CALL getvdwforce(rho)
  END IF  ! if ivdw.eq.2 vdw is done implicitely by vArElCore
  
END IF
#endif


OPEN(772,FILE='forces.'//outnam)
DO i=1,nion
  WRITE(772,*) fx(i)
  WRITE(772,*) fy(i)
  WRITE(772,*) fz(i)
ENDDO
CLOSE(772)


CALL laserf(rho)



!      stop
RETURN
END SUBROUTINE getforces
!------------------------------------------------------------



!------------------------------------------------------------

SUBROUTINE getforcenana()

! Force between ionic cores
! Input and output communicated via module 'params'

USE params
IMPLICIT NONE


!REAL(DP), INTENT(IN)                     :: rho(2*kdfull2)

INTEGER :: ii, jj
REAL(DP) :: dist, dist2, radfor, forcex, forcey, forcez
REAL(DP) :: xi, yi, zi, xr, yr, zr


!     Na(core)-Na(core) forces

DO ii=1,nion-1
  xi = cx(ii)
  yi = cy(ii)
  zi = cz(ii)
  
  DO jj=ii+1,nion
    
    xr = xi - cx(jj)
    yr = yi - cy(jj)
    zr = zi - cz(jj)
    
    dist2 = xr**2+yr**2+zr**2
    dist = SQRT(dist2)
    
    radfor = -e2*ch(np(ii))*ch(np(jj))/dist2
    
    forcex = radfor*xr/dist
    forcey = radfor*yr/dist
    forcez = radfor*zr/dist
    
    fx(ii) = fx(ii) - forcex
    fy(ii) = fy(ii) - forcey
    fz(ii) = fz(ii) - forcez
    
    fx(jj) = fx(jj) + forcex
    fy(jj) = fy(jj) + forcey
    fz(jj) = fz(jj) + forcez
    
!    CALL getshortforce(4,4,ii,jj,rho,0)
    
    
  END DO
  
!          write(6,'(a,i,3e15.5)') 'Na-Na:',ii,fx(ii),fy(ii),fz(ii)
  
END DO


#if(raregas)
!     Na(core)-Na(core)image forces
IF(idielec /= 0) THEN
  DO ii=1,nion
    xi = cx(ii)
    yi = cy(ii)
    zi = cz(ii)
    
    DO jj=1,nion
      xr = xi + cx(jj) - 2D0*xdielec
      yr = yi - cy(jj)
      zr = zi - cz(jj)
      
      dist2 = xr**2+yr**2+zr**2
      dist = SQRT(dist2)
      
      radfor = (epsdi-1D0)/(epsdi+1D0)*e2*ch(np(ii))* ch(np(jj))/dist2
      
      forcex = radfor*xr/dist
      forcey = radfor*yr/dist
      forcez = radfor*zr/dist
      
      fx(ii) = fx(ii) - forcex
      fy(ii) = fy(ii) - forcey
      fz(ii) = fz(ii) - forcez
      
    END DO
    
  END DO
  
END IF
#endif

RETURN
END SUBROUTINE getforcenana
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE getforceelna(rho)

! Force from electron cloud on ions.
!
! Input:
!   rho   = electronic local density
!   Other I/O is communicated via module 'params'

USE params
IMPLICIT NONE


REAL(DP), INTENT(IN OUT)   :: rho(2*kdfull2)


INTEGER :: ii, ind
REAL(DP) :: prefac
#if(raregas)
REAL(DP),ALLOCATABLE ::  rhotmp(:)
#endif

EXTERNAL v_soft,gauss
EXTERNAL dgaussdr


DO ii=1,nion
  
  IF (ipseudo == 0) THEN
#if(raregas)    
    IF(idielec == 1) THEN
      
      ALLOCATE(rhotmp(2*kdfull2))
      DO ind=1,2*kdfull2
        rhotmp(ind)=rho(ind)
      END DO
      
      CALL addimage(rho,1)
      
    END IF
#endif    
    CALL foldgradfunc(rho,v_soft,cx(ii),cy(ii),cz(ii), sgm1(np(ii))*sq2)
! contribution from first gaussian
! the plus sign in the forces is really a plus because
! rho is positive but the electronic charge is -rho
    
    
    fx(ii) = fx(ii) + e2*chg1(np(ii))*rvectmp(1)
    fy(ii) = fy(ii) + e2*chg1(np(ii))*rvectmp(2)
    fz(ii) = fz(ii) + e2*chg1(np(ii))*rvectmp(3)
    
    
    CALL foldgradfunc(rho,v_soft,cx(ii),cy(ii),cz(ii), sgm2(np(ii))*sq2)
! contribution from second gaussian
    
    fx(ii) = fx(ii) + e2*chg2(np(ii))*rvectmp(1)
    fy(ii) = fy(ii) + e2*chg2(np(ii))*rvectmp(2)
    fz(ii) = fz(ii) + e2*chg2(np(ii))*rvectmp(3)
    
#if(rearegas)
    IF(idielec == 1) THEN
      
      DO ind=1,2*kdfull2
        rho(ind)=rhotmp(ind)
      END DO
      DEALLOCATE(rhotmp)

    END IF
#endif
    
  ELSE  ! ipseudo = 1
    
    CALL foldgradfunconsubgrid(chpcoul,gauss,cx(ii),  &
        cy(ii),cz(ii),sgm1(np(ii))*sq2)
    
    prefac = chg1(np(ii))/(2D0*pi*sgm1(np(ii))**2)**1.5D0
! no factor e2, because it is already in the field chpcoul()
    
    fx(ii) = fx(ii) + prefac*rvectmp(1)
    fy(ii) = fy(ii) + prefac*rvectmp(2)
    fz(ii) = fz(ii) + prefac*rvectmp(3)
    
    CALL foldgradfunconsubgrid(chpcoul,gauss,cx(ii),  &
        cy(ii),cz(ii),sgm2(np(ii))*sq2)
    
    prefac = chg2(np(ii))/(2*pi*sgm2(np(ii))**2)**1.5D0
! no factor e2, because it is already in the field chpcoul()
    
    fx(ii) = fx(ii) + prefac*rvectmp(1)
    fy(ii) = fy(ii) + prefac*rvectmp(2)
    fz(ii) = fz(ii) + prefac*rvectmp(3)
    
  END IF
  
#if(raregas)
  CALL getshortforce(4,5,ii,0,rho,0)
#endif
  
END DO

RETURN
END SUBROUTINE getforceelna
!------------------------------------------------------------



!     ******************************

SUBROUTINE calcf_goeloc(rho,it)

! Force from local part of Goedecketr PsP on ions
! Input:
!   rho   = electronic local density
!   it    = time step nr. in calling routine
!   Other I/O is communicated via module 'params'

USE params
IMPLICIT NONE

REAL(DP), INTENT(IN)    :: rho(2*kdfull2)
INTEGER, INTENT(IN)     :: it

CHARACTER (LEN=1) :: ext
INTEGER :: ind, ion, ion1, is, ix, iy, iz
REAL(DP) :: c1, c2, chpddr, dist2, dist3, pch
REAL(DP) :: r2, rdn, rloc, zion
REAL(DP) :: rr, rr3,  rx, ry, rz, x1, y1, z1

REAL(DP),EXTERNAL :: v_ion_el_lgoed

REAL(DP),PARAMETER :: rder=1D-1    ! radius step for finite difference


DO ion=1,nion
  fx(ion)=0D0
  fy(ion)=0D0
  fz(ion)=0D0
END DO

!     derivatives of the n goedecker pseudos

DO is=1,nion
  c1=cc1(np(is))
  c2=cc2(np(is))
  rloc=crloc(np(is))*SQRT(2D0)
  ind=0
  DO iz=minz,maxz
    z1=(iz-nzsh)*dz
    DO iy=miny,maxy
      y1=(iy-nysh)*dy
      DO ix=minx,maxx
        x1=(ix-nxsh)*dx
        ind=ind+1
        IF((ix /= nx2).AND.(iy /= ny2).AND.(iz /= nz2))THEN
          rdn=2D0/SQRT(pi)*e2
          zion=ch(np(is))*e2
          rx=x1-cx(is)
          ry=y1-cy(is)
          rz=z1-cz(is)
          r2=rx*rx+ry*ry+rz*rz+1D-12
          rr=SQRT(r2)
          rr3=rr*r2
          chpddr =-(v_ion_el_lgoed(rr+rder,rloc,c1,c2,zion)  &
                  -v_ion_el_lgoed(rr-rder,rloc,c1,c2,zion)) /(rder+rder)/rr
          fx(is)=fx(is)+(chpddr*rho(ind))*rx
          fy(is)=fy(is)+(chpddr*rho(ind))*ry
          fz(is)=fz(is)+(chpddr*rho(ind))*rz
        END IF
      END DO
    END DO
  END DO
END DO
DO ion=1,nion
  fx(ion)=fx(ion)*dvol
  fy(ion)=fy(ion)*dvol
  fz(ion)=fz(ion)*dvol
END DO

!     force from ion-ion interactions

DO ion=1,nion
  DO ion1=1,nion
    IF(ion1 /= ion) THEN
      dist2   = (cx(ion1)-cx(ion))**2+(cy(ion1)-cy(ion))**2  &
          +(cz(ion1)-cz(ion))**2
      dist3=-SQRT(dist2)*dist2
      pch=e2*ch(np(ion))*ch(np(ion1))
      fx(ion)=fx(ion)-pch*(cx(ion)-cx(ion1))/dist3
      fy(ion)=fy(ion)-pch*(cy(ion)-cy(ion1))/dist3
      fz(ion)=fz(ion)-pch*(cz(ion)-cz(ion1))/dist3
    END IF
  END DO
  
  IF(jforce /= 0 .AND. MOD(it,jforce) == 0 .AND. it >= 0) THEN
    WRITE(ext,'(i1)') ion 
    OPEN(24,POSITION='append',FILE='pforce.'//ext//'.'//outnam)
    WRITE(24,'(4f13.5)') tfs,fx(ion),fy(ion),fz(ion)
  END IF
END DO


!     force from laser

CALL laserf(rho)

!     protocol

IF(jforce /= 0 .AND. MOD(it,jforce) == 0 .AND. it >= 0) THEN
  DO ion=1,nion
    WRITE(ext,'(i1)') ion 
    OPEN(25,POSITION='append',FILE='plforce.'//ext//'.'//outnam)
    WRITE(25,'(4f13.5)') tfs,flx(ion),fly(ion),flz(ion)
  END DO
END IF

RETURN
END SUBROUTINE calcf_goeloc


SUBROUTINE forceproject()

! Force from colliding projectile (if any).
! Input:
!   rho   = electronic local density
!   Other I/O is communicated via module 'params'


USE params
IMPLICIT NONE

INTEGER :: ion
REAL(DP) :: dist2, dist3, pch

IF (ABS(projcharge) <= 1.0D-5) RETURN


OPEN(8886,POSITION='append',FILE='projforce'//outnam)


DO ion=1,nion
  
  dist2   = (cx(ion)-(projvelx*tfs+projinix))**2 + &
       (cy(ion)-(projvely*tfs+projiniy))**2 + &
       (cz(ion)-(projvelz*tfs+projiniz))**2
  dist3=-SQRT(dist2)*dist2
  pch=e2*ch(np(ion))*projcharge
  
  fprojx(ion) = -pch*(cx(ion)-(projvelx*tfs+projinix))/dist3
  fprojy(ion) = -pch*(cy(ion)-(projvely*tfs+projiniy))/dist3
  fprojz(ion) = -pch*(cz(ion)-(projvelz*tfs+projiniz))/dist3
  
  WRITE(8886,'(4f13.5)') tfs,fprojx(ion),fprojy(ion),fprojz(ion)
  fx(ion) = fx(ion) + fprojx(ion)
  fy(ion) = fy(ion) + fprojy(ion)
  fz(ion) = fz(ion) + fprojz(ion)
  
END DO

RETURN

END SUBROUTINE forceproject


SUBROUTINE laserf(rho)

! Force from external laser pulse on ionic cores.

USE params
USE util, ONLY:laserp
IMPLICIT NONE

REAL(DP), INTENT(IN)                         :: rho(2*kdfull2)

INTEGER :: ind
REAL(DP) :: ascal, foft
REAL(DP) :: tstart, tend, tvend, tmax, try
REAL(DP) :: snorm, sx, sy, sz


IF (ABS(e0) <= 1D-12) THEN
  RETURN
END IF


ascal=e1x*e2x+e1y*e2y+e1z*e2z
sx=e1x+e2x*COS(phi)
sy=e1y+e2y*COS(phi)
sz=e1z+e2z*COS(phi)
snorm=sx*sx+sy*sy+sz*sz
snorm=SQRT(snorm)
e1x=e1x/snorm
e1y=e1y/snorm
e1z=e1z/snorm
e2x=e2x/snorm
e2y=e2y/snorm
e2z=e2z/snorm

!     prepare time profile

IF(itft == 1) THEN
  tstart = tnode  + tpeak
  tend   = tstart + deltat
  tvend  = tend   + tpeak
  
  IF(tfs <= tnode.OR.tfs >= tvend) foft = 0D0
  IF(tfs > tnode.AND.tfs <= tstart) foft = SIN(0.5D0*pi*(tfs-tnode)/tpeak)
  IF(tfs > tstart.AND.tfs <= tend) foft = 1D0
  IF(tfs > tend.AND.tfs < tvend) foft = SIN (0.5D0 * pi * (tfs-tend) / tpeak)
END IF

IF(itft == 2) THEN
  tmax = tnode + tpeak
  
  IF(tfs <= tnode) foft = 0D0
  IF(tfs > tnode) foft = EXP (- ((tfs - tmax) / deltat)**2)
END IF

IF(itft == 3) THEN
  tvend = tnode + 2D0*tpeak + deltat
  IF(tfs <= tnode.OR.tfs >= tvend) THEN
    foft = 0D0
  ELSE
    foft = COS((-0.5D0+(tfs-tnode)/(2*tpeak+deltat))*pi)**2
  END IF
END IF

IF(itft == 4) THEN
  tvend = tnode + 2*tpeak + deltat
  IF(tfs <= tnode.OR.tfs >= tvend) THEN
    foft = 0.
  ELSE
    foft = COS((-0.5D0+(tfs-tnode)/(2*tpeak+deltat))*pi)**4
  END IF
END IF

IF(itft <= 0 .OR. itft >= 5) THEN
  STOP ' this pulse profile not yet implemented'
END IF
power=e0*e0*foft*foft


!     calculate force of the laser field


try = tfs/0.0484D0

DO ind=1,nion
  flx(ind) = - e0*e1x*COS(omega*try+phi)*e2*ch(np(ind))*foft
  fly(ind) = - e0*e1y*COS(omega*try+phi)*e2*ch(np(ind))*foft
  flz(ind) = - e0*e1z*COS(omega*try+phi)*e2*ch(np(ind))*foft
  
!wzp         write(8887,'(4f13.5)') tfs,flx(ind),fly(ind),flz(ind)
  fx(ind) = fx(ind)+flx(ind)
  fy(ind) = fy(ind)+fly(ind)
  fz(ind) = fz(ind)+flz(ind)
END DO

#if(raregas)
DO ind=1,nc
  fxc(ind) = fxc(ind) - e0*e1x*COS(omega*try+phi)*e2*chgc(ind)*foft
  fyc(ind) = fyc(ind) - e0*e1y*COS(omega*try+phi)*e2*chgc(ind)*foft
  fzc(ind) = fzc(ind) - e0*e1z*COS(omega*try+phi)*e2*chgc(ind)*foft
  
  fxe(ind) = fxe(ind) - e0*e1x*COS(omega*try+phi)*e2*chge(ind)*foft
  fye(ind) = fye(ind) - e0*e1y*COS(omega*try+phi)*e2*chge(ind)*foft
  fze(ind) = fze(ind) - e0*e1z*COS(omega*try+phi)*e2*chge(ind)*foft
END DO

DO ind=1,nk
  fxk(ind) = fxk(ind) - e0*e1x*COS(omega*try+phi)*e2*chgk(ind)*foft
  fyk(ind) = fyk(ind) - e0*e1y*COS(omega*try+phi)*e2*chgk(ind)*foft
  fzk(ind) = fzk(ind) - e0*e1z*COS(omega*try+phi)*e2*chgk(ind)*foft
END DO
#endif

IF (nclust == 0) THEN
  CALL laserp(potion,rho) ! dummy field
END IF

RETURN
END SUBROUTINE laserf
!     ********************************

SUBROUTINE calcf_goenonl(rho,it,psi)

! Force from non-local Goedecker PsP on ionic cores.
! Input:
!   rho   = electronic local density
!   Other I/O is communicated via module 'params'

USE params
USE util, ONLY:realoverlap, realovsubgrid
IMPLICIT NONE

#if(parayes)
INCLUDE 'mpif.h'
#endif

REAL(DP), INTENT(IN)                         :: rho(2*kdfull2)
INTEGER, INTENT(IN)                  :: it
COMPLEX(DP), INTENT(IN)                  :: psi(kdfull2,kstate)


INTEGER :: i, ion, ion1, nb
REAL(DP) :: dist2, dist3, forceloc, forcenonl, pch
REAL(DP) :: sumfor,sumfulp, sumfulm,  sumslp, sumslm
REAL(DP) :: xion, yion, zion 
REAL(DP) :: zshift=0.001D0,yshift=0.001D0,xshift=0.001D0
REAL(DP) :: zshinv=1D3,yshinv=1D3,xshinv=1D3              ! must be inverse of zshift
COMPLEX(DP),ALLOCATABLE :: q1(:)
REAL(DP),ALLOCATABLE :: rhoslp(:),rhoslm(:)
REAL(DP),ALLOCATABLE :: fxnl(:),fynl(:),fznl(:)
#if(parayes)
REAL(DP),ALLOCATABLE :: ftemp(:)
#endif

CHARACTER (LEN=1) :: ext

ALLOCATE(q1(kdfull2),rhoslp(kdfull2),rhoslm(kdfull2))
ALLOCATE(fxnl(nion),fynl(nion),fznl(nion))

!     force on the ions

DO ion=1,nion
  fx(ion)=0D0
  fy(ion)=0D0
  fz(ion)=0D0
END DO
fznl=0D0; fynl=0D0; fxnl=0D0

DO ion = 1,nion
  
  xion = cx(ion)
  yion = cy(ion)
  zion = cz(ion)
  
  CALL rhopsg(xion,yion,zion+zshift*0.5D0,rhoslp,ion)
  CALL rhopsg(xion,yion,zion-zshift*0.5D0,rhoslm,ion)
  sumfor = 0D0
  DO i = 1, nxyz
    sumfor = sumfor + rho(i)*(rhoslp(i)-rhoslm(i))
  END DO
  fz(ion) = sumfor*dvol * zshinv
    
  IF(tnonloc(ion)) THEN
    CALL calc_proj(xion,yion,zion+zshift*0.5D0,xion,yion,zion,ion)
    sumslp = 0D0
    DO nb=1,nstate
      CALL nonlocalc(psi(1,nb),q1,ion)
      sumslp = realovsubgrid(psi(1,nb),q1,ion) + sumslp
    END DO
    
    CALL calc_proj(xion,yion,zion-zshift*0.5D0,xion,yion,zion,ion)
    sumslm = 0D0
    DO nb=1,nstate
      CALL nonlocalc(psi(1,nb),q1,ion)
      sumslm = realovsubgrid(psi(1,nb),q1,ion) + sumslm
    END DO
    fznl(ion) = (sumslm - sumslp) * zshinv
  END IF

  CALL rhopsg(xion,yion+yshift*0.5D0,zion,rhoslp,ion)
  CALL rhopsg(xion,yion-yshift*0.5D0,zion,rhoslm,ion)
  sumfor = 0D0
  DO i = 1, nxyz
    sumfor = sumfor + rho(i)*(rhoslp(i)-rhoslm(i))
  END DO
  fy(ion) = sumfor*dvol * yshinv
    
  IF(tnonloc(ion)) THEN
    CALL calc_proj(xion,yion+yshift*0.5D0,zion,xion,yion,zion,ion)
    sumslp = 0D0
    DO nb=1,nstate
      CALL nonlocalc(psi(1,nb),q1,ion)
      sumslp = realovsubgrid(psi(1,nb),q1,ion) + sumslp
    END DO
    CALL calc_proj(xion,yion-yshift*0.5D0,zion,xion,yion,zion,ion)
    sumslm = 0D0
    DO nb=1,nstate
      CALL nonlocalc(psi(1,nb),q1,ion)
      sumslm = realovsubgrid(psi(1,nb),q1,ion) + sumslm
    END DO
    fynl(ion) = (sumslm - sumslp) * yshinv
  END IF
    
  CALL rhopsg(xion+xshift*0.5D0,yion,zion,rhoslp,ion)
  CALL rhopsg(xion-xshift*0.5D0,yion,zion,rhoslm,ion)
  sumfor = 0D0
  DO i = 1, nxyz
    sumfor = sumfor + rho(i)*(rhoslp(i)-rhoslm(i))
  END DO
  fx(ion) = sumfor*dvol * xshinv
  forceloc = sumfor * xshinv              ! for testing
  IF(tnonloc(ion)) THEN
    CALL calc_proj(xion+xshift*0.5D0,yion,zion,xion,yion,zion,ion)
    sumslp = 0D0
    sumfulp = 0D0
    DO nb=1,nstate
      CALL nonlocalc(psi(1,nb),q1,ion)
      sumslp = realovsubgrid(psi(1,nb),q1,ion) + sumslp
      sumfulp = realoverlap(psi(1,nb),q1) + sumfulp
    END DO
    CALL calc_proj(xion-xshift*0.5D0,yion,zion,xion,yion,zion,ion)
    sumslm = 0D0
    sumfulm = 0D0
    DO nb=1,nstate
      CALL nonlocalc(psi(1,nb),q1,ion)
      sumslm = realovsubgrid(psi(1,nb),q1,ion) + sumslm
      sumfulm = realoverlap(psi(1,nb),q1) + sumfulm
    END DO
    fxnl(ion) = (sumslm - sumslp) * xshinv
!           sudiff = (sumfulm - sumfulp)
!      sumfor = sumfor + sudiff
  END IF
!    fx(ion) = sumfor * xshinv
  forcenonl = fxnl(ion) 

  sumslm =  sumslm * xshinv
  sumslp =  sumslp * xshinv
  sumfulp =  sumfulp * xshinv
  sumfulm =  sumfulm * xshinv
  
!       optional prints
  
  IF(jforce /= 0 .AND. MOD(it,jforce) == 0 .AND. it >= 0) THEN
    WRITE(ext,'(i1)') ion 
!            write(*,*) ' open 25: ion,ext=',ion,ext
    OPEN(25,POSITION='append',FILE='pfrcel.'//ext//'.'//outnam)
    WRITE(25,'(10f13.5)')  tfs,fx(ion),fy(ion),fz(ion),  &
        forceloc,forcenonl,sumslm,sumslp,sumfulm,sumfulp
    CLOSE(25)
  END IF
  
!      restore  projectors for original ionic positions
  CALL calc_proj(xion,yion,zion,xion,yion,zion,ion)

END DO

! compose total force

#if(parayes)
  ALLOCATE(ftemp(nion))
  CALL mpi_barrier (mpi_comm_world, mpi_ierror)
  CALL mpi_allreduce(fxnl,ftemp,nion,mpi_double_precision,  &
                     mpi_sum,mpi_comm_world,mpi_ierror)
  fx(1:nion) = fx(1:nion)+ftemp(1:nion)
  CALL mpi_barrier (mpi_comm_world, mpi_ierror)
  CALL mpi_allreduce(fynl,ftemp,nion,mpi_double_precision,  &
                     mpi_sum,mpi_comm_world,mpi_ierror)
  fy(1:nion) = fy(1:nion)+ftemp(1:nion)
  CALL mpi_barrier (mpi_comm_world, mpi_ierror)
  CALL mpi_allreduce(fznl,ftemp,nion,mpi_double_precision,  &
                     mpi_sum,mpi_comm_world,mpi_ierror)
  fz(1:nion) = fz(1:nion)+ftemp(1:nion)
  CALL mpi_barrier (mpi_comm_world, mpi_ierror)
  DEALLOCATE(ftemp)
#else
  fx(1:nion) = fx(1:nion)+fxnl(1:nion)
  fy(1:nion) = fy(1:nion)+fynl(1:nion)
  fz(1:nion) = fz(1:nion)+fznl(1:nion)
#endif

DEALLOCATE(fxnl,fynl,fznl)


DO ion=1,nion
  DO ion1=1,nion
    IF(ion1 /= ion) THEN
      dist2   = (cx(ion1)-cx(ion))**2+(cy(ion1)-cy(ion))**2  &
          +(cz(ion1)-cz(ion))**2
      dist3=-SQRT(dist2)*dist2
      pch=e2*ch(np(ion))*ch(np(ion1))
      fx(ion)=fx(ion)-pch*(cx(ion)-cx(ion1))/dist3
      fy(ion)=fy(ion)-pch*(cy(ion)-cy(ion1))/dist3
      fz(ion)=fz(ion)-pch*(cz(ion)-cz(ion1))/dist3
    END IF
  END DO
  
  IF(jforce /= 0 .AND. MOD(it,jforce) == 0 .AND. it >=0) THEN
    WRITE(ext,'(i1)') ion 
    OPEN(24,POSITION='append',FILE='pforce.'//ext//'.'//outnam)
    WRITE(24,'(4f13.5)') tfs,fx(ion),fy(ion),fz(ion)
    CLOSE(24)
  END IF
  
END DO
DEALLOCATE(q1,rhoslp,rhoslm)

RETURN
END SUBROUTINE calcf_goenonl

!-----rhopsg-----------------------------------------------------------

!     **********************************************

SUBROUTINE rhopsg(cxact,cyact,czact,rhopsp,is)

! Pseudo-density of ion 'is' related to a local PsP.
!
! Input:
!   cxact,cyact,czact = coordinates of ion
!   is                = nr. of ion
! Output:
!   rhopsp            = emerging pseudo-density

USE params
IMPLICIT NONE


INTEGER, INTENT(IN)    :: is
REAL(DP), INTENT(IN)   :: cxact
REAL(DP), INTENT(IN)   :: cyact
REAL(DP), INTENT(IN)   :: czact
REAL(DP), INTENT(OUT)  :: rhopsp(kdfull2)

INTEGER :: ind, ix, iy, iz
REAL(DP) :: c1, c2, f1, f2, pt1, pt2, rloc, zion
REAL(DP) :: rr, rx, ry, rz, x1, y1, z1

REAL(DP),EXTERNAL :: v_ion_el_lgoed
REAL(DP),EXTERNAL :: v_soft

DO ind=1,nxyz
  rhopsp(ind)=0D0
END DO


!     soft local PsP

IF(ipsptyp == 0) THEN
  ind=0
  DO iz=minz,maxz
    z1=(iz-nzsh)*dz
    DO iy=miny,maxy
      y1=(iy-nysh)*dy
      DO ix=minx,maxx
        x1=(ix-nxsh)*dx
        ind=ind+1
        rx=x1-cxact
        ry=y1-cyact
        rz=z1-czact
        rr=SQRT(rx*rx+ry*ry+rz*rz)
        rr=rr+1D-10             ! avoid zero
        f1=(v_soft(rr,(sgm1(np(is))*SQRT(2D0))))
        f2=(v_soft(rr,(sgm2(np(is))*SQRT(2D0))))
        pt1=e2*(f1)
        pt2=e2*(f2)
        rhopsp(ind)=rhopsp(ind)+chg1(np(is))*pt1+chg2(np(is))*pt2
      END DO
    END DO
  END DO
  
  
!     local part of Goedecker
  
ELSE IF(ipsptyp >= 1) THEN
  c1=cc1(np(is))
  c2=cc2(np(is))
  rloc=crloc(np(is))
  zion=ch(np(is))
  
  ind=0
  DO iz=minz,maxz
    z1=(iz-nzsh)*dz
    DO iy=miny,maxy
      y1=(iy-nysh)*dy
      DO ix=minx,maxx
        x1=(ix-nxsh)*dx
        ind=ind+1
        rx=x1-cxact
        ry=y1-cyact
        rz=z1-czact
        rr=SQRT(rx*rx+ry*ry+rz*rz)
        rr=rr+1D-6             ! avoid zero
        rhopsp(ind)=rhopsp(ind)+v_ion_el_lgoed(rr,rloc,c1,c2,zion)
      END DO
    END DO
  END DO
  
ELSE
  STOP 'this type of PsP not yet provided'
END IF

RETURN
END SUBROUTINE rhopsg


REAL(DP) FUNCTION vgian(r)
USE params, ONLY: DP
IMPLICIT NONE

! returns the Gianozzi hydrogen PP in Rydberg units.
! Reference: F. Gygi, PRB 48, 11692 (1993).
! This PP was used for the calculations in
! K"ummel, Kronik, Perdew, PRL 93, 213002 (2004)

REAL(DP), INTENT(IN)     :: r

REAL(DP), PARAMETER :: rc1=0.25D0
REAL(DP), PARAMETER :: rc2=0.284D0
REAL(DP), PARAMETER :: a=-1.9287D0*2D0
REAL(DP), PARAMETER :: b=0.3374D0*2D0

vgian=-2.d0*erf(r/rc1)/r + (a+b*r**2)*EXP(-(r/rc2)**2)

END FUNCTION vgian



 
!------------------------------------------------------------
!
REAL(DP) FUNCTION gauss(r,s)
!------------------------------------------------------------
USE params, ONLY:DP
IMPLICIT NONE
REAL(DP),INTENT(IN) :: r
REAL(DP),INTENT(IN) :: s
REAL(DP)::rs
rs = r/s
gauss = EXP(-rs*rs)
!      endif
RETURN
END FUNCTION gauss
!------------------------------------------------------------

REAL(DP) FUNCTION dgaussdr(r,s)
!------------------------------------------------------------
! Derivative of a Gaussian of width 's' at coordinate 'r'.
USE params, ONLY: DP
IMPLICIT NONE

REAL(DP),INTENT(IN)  :: r
REAL(DP),INTENT(IN)  :: s
REAL(DP)::ra

ra = r/s
dgaussdr = - 2D0*ra*EXP(-(ra*ra))/s
RETURN
END FUNCTION dgaussdr

