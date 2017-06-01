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

#include "define.h"
 
!------------------------------------------------------------

SUBROUTINE calcchargdist(field)
!------------------------------------------------------------
USE params
USE util, ONLY:getcm
IMPLICIT NONE

REAL(DP), INTENT(IN)                     :: field(kdfull2)

INTEGER :: ii, ind, ix, iy, iz, n
REAL(DP) :: chtotal, r, rr, x1, y1, z1
REAL(DP) :: chfld(100)


chtotal = 0.0D0
DO ii = 1,INT(nzsh*dz/drcharges)
  chfld(ii)=0.0D0
END DO

ind = 0

CALL getcm(1,0,0)

DO iz=minz,maxz
  z1=(iz-nzsh)*dz
  DO iy=miny,maxy
    y1=(iy-nysh)*dy
    DO ix=minx,maxx
      x1=(ix-nxsh)*dx
      
      rr = (x1-rvectmp(1))**2+(y1-rvectmp(2))**2+ (z1-rvectmp(3))**2
      
      ind = ind + 1
      
      chtotal = chtotal + field(ind)
      
      DO ii = 1,INT(nzsh*dz/drcharges)
        r = ii * drcharges
        IF (rr <= r*r) THEN
          chfld(ii)=chfld(ii)+field(ind)
        END IF
      END DO
      
    END DO
  END DO
END DO

DO ii = 1,INT(nzsh*dz/drcharges)
  chfld(ii)=chfld(ii)*dvol
END DO

chtotal = chtotal * dvol

WRITE(323,'(f15.5,101e17.7)') tfs, chtotal,  &
    (chfld(n),n=1,INT(nzsh*dz/drcharges))


RETURN
END SUBROUTINE calcchargdist

!------------------------------------------------------------

SUBROUTINE evaluate(rho,aloc,psi)
!------------------------------------------------------------
USE params
IMPLICIT NONE


REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
REAL(DP), INTENT(IN OUT)                     :: aloc(2*kdfull2)
COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)
#if(raregas)
INTEGER :: i
#endif
REAL(DP) :: dummy, enpol, enpol0
REAL(DP), EXTERNAL :: energ_ions

WRITE(6,*) 'Doing some postrun evaluation only...'


CALL calcrho(rho,psi)
CALL calclocal(rho,aloc)
CALL calcpseudo()
CALL calcrho(rho,psi)
IF(ifsicp > 0) CALL calc_sic(rho,aloc,psi)
jenergy=1

CALL info(psi,rho,aloc,10)

WRITE(6,*) 'Reading restart File'
nion=0
!     call energ_ions() ??? BF
dummy=energ_ions()

OPEN(834,STATUS='unknown',FILE='penerinfty')
#if(raregas)
enerinfty=engg
WRITE(834,*) enerinfty
#else
WRITE(834,*) engg
#endif

CLOSE(834)

enpol0=0D0
#if(raregas)
DO i=1,nc
  enpol0=enpol0+(xc(i)-xe(i))**2+(yc(i)-ye(i))**2+ (zc(i)-ze(i))**2
END DO
enpol0=0.5D0*cspr*enpol0
#endif

CALL restart2(psi,outnam,.false.)
WRITE(6,*) 'Done.'

CLOSE(163)

CALL calcrho(rho,psi)
CALL calclocal(rho,aloc)
CALL calcpseudo()
CALL calcrho(rho,psi)
IF(ifsicp > 0) CALL calc_sic(rho,aloc,psi)
jenergy=1
CALL info(psi,rho,aloc,10)

enpol=0D0
#if(raregas)
DO i=1,nc
  enpol=enpol+(xc(i)-xe(i))**2+(yc(i)-ye(i))**2+ (zc(i)-ze(i))**2
END DO
enpol=0.5D0*cspr*enpol
WRITE(6,'(a,3f17.8)') 'enii,enig,engg',enii,enig ,engg-enerinfty
#else
WRITE(6,'(a,3f17.8)') 'enii,enig,engg',enii,enig ,engg
#endif
WRITE(6,'(a,2f17.8)') 'enpol0,enpol',enpol0,enpol

WRITE(6,*) 'Postrun evaluation done.'

STOP 'Postrun evaluation done.'


!old      if (iflag.eq.1) then
!old
!old         call restart2(psi,outnam)
!old
!oldc         nsav=nion
!oldc         nion=0
!old         call calcrho(rho,psi)
!old         call calcpseudo()
!old
!old         eMgO=energ_ions()
!old         call info(psi,rho,aloc,akv,7)
!oldc        nion=nsav
!old
!old         write(6,*) eMgO
!old
!old      elseif (iflag.eq.2) then
!old      endif


RETURN
END SUBROUTINE evaluate
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE getsurfprops
!------------------------------------------------------------
USE params
IMPLICIT NONE

! calculate kinetic energy of surface ions
#if(raregas)
INTEGER :: i

ekincsurf=0D0
ekinesurf=0D0
ekinksurf=0D0
ekinsurf=0D0

DO i=1,nc
  ekincsurf=ekincsurf+(pxc(i)*pxc(i)+pyc(i)*pyc(i)+pzc(i)*pzc(i))/  &
      2D0/mion/1836D0/ame
END DO
DO i=1,NE
  ekinesurf=ekinesurf+(pxe(i)*pxe(i)+pye(i)*pye(i)+pze(i)*pze(i))/  &
      2D0/me/1836D0/ame
END DO
DO i=1,NE
  ekinksurf=ekinksurf+(pxk(i)*pxk(i)+pyk(i)*pyk(i)+pzk(i)*pzk(i))/  &
      2D0/mkat/1836D0/ame
END DO

ekinsurf = ekincsurf + ekinesurf + ekinksurf
#endif

RETURN
END SUBROUTINE getsurfprops
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE getclustergeometry
!------------------------------------------------------------
USE params
IMPLICIT NONE
!     calculates characteristic observables for describing
!     the geometry of the cluster ions

INTEGER :: i, ico, ion, j
REAL(DP) :: dist, rr, sumx, sumy, sumz
REAL(DP) :: cw(ng,3)

! get center of mass first

sumx=0D0
sumy=0D0
sumz=0D0

DO i=1,nion
  sumx = sumx + cx(i)
  sumy = sumy + cy(i)
  sumz = sumz + cz(i)
END DO

comx = sumx/nion
comy = sumy/nion
comz = sumz/nion

! now use center of mass as origin

DO i=1,nion
  cw(i,1)=cx(i)-comx
  cw(i,2)=cy(i)-comy
  cw(i,3)=cz(i)-comz
END DO


!     ! calculate quadrupole tensor of mass distribution

DO i=1,3
  DO j=i,3
    qtion(i,j)=0D0
    DO ion=1,nion
      qtion(i,j)=qtion(i,j) + 3*cw(ion,i)*cw(ion,j)
      IF (i == j) THEN
        rr=cw(ion,1)**2+cw(ion,2)**2+cw(ion,3)**2
        qtion(i,j)=qtion(i,j)-rr
      END IF
    END DO
    qtion(i,j)=qtion(i,j)/nion
    qtion(j,i)=qtion(i,j)
  END DO
END DO

!      ! calculate root mean square radius

rmsion=0D0

DO ion=1,nion
  rmsion=rmsion+cw(ion,1)**2+cw(ion,2)**2+cw(ion,3)**2
END DO

rmsion=SQRT(rmsion/nion)

! determine maximum distance; allows to see if
! cluster is dissociating

dmdistion=0D0
DO i=1,nion-1
  DO j=i+1,nion
    dist=0D0
    DO ico=1,3
      dist = dist + (cw(i,ico)-cw(j,ico))**2
    END DO
    IF (dist > dmdistion) dmdistion=dist
  END DO
END DO

dmdistion=SQRT(dmdistion)


RETURN
END SUBROUTINE getclustergeometry
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE getelectrongeometry(q0,nbe)
!------------------------------------------------------------
USE params
IMPLICIT NONE


COMPLEX(DP), INTENT(IN)                      :: q0(kdfull2,kstate)
INTEGER, INTENT(IN)                      :: nbe

INTEGER :: i, ind, ix, iy, iz, j, k, ncount
REAL(DP) :: rr, sum,  x1, y1, z1
REAL(DP) :: xic(3)
REAL(DP),ALLOCATABLE :: q1(:)

ALLOCATE(q1(kdfull2))

!     calculates geometry observables of electron orbital density
!     number nbe.
!     for total electronic density, set nbe=0
!     for total spin-up density, set nbe=-1
!     for total spin-down density, set nbe=-2


q1=0D0

IF (nbe == 0) THEN ! total density
  ncount=nclust
  DO i=1,nstate
    DO ind=1,kdfull2
      q1(ind)=q1(ind)+q0(ind,i)*CONJG(q0(ind,i))
    END DO
  END DO
ELSE IF (nbe == -1) THEN ! total spin-up density
  ncount=nclust-nspdw
  STOP 'analyse.F: to be implemented'
ELSE IF (nbe == -2) THEN   ! total spin-down density
  ncount=nspdw
  STOP 'analyse.F: to be implemented'
ELSE ! single orbital density
  ncount=1
  DO ind=1,kdfull2
    q1(ind)=q1(ind)+q0(ind,nbe)*CONJG(q0(ind,nbe))
  END DO
END IF

! get center of density

codx=0D0
cody=0D0
codz=0D0

ind=0
DO iz=minz,maxz
  z1=(iz-nzsh)*dz
  DO iy=miny,maxy
    y1=(iy-nysh)*dy
    DO ix=minx,maxx
      x1=(ix-nxsh)*dx
      
      ind=ind+1
      codx = codx + q1(ind)*x1
      cody = cody + q1(ind)*y1
      codz = codz + q1(ind)*z1
      
    END DO
  END DO
END DO

codx = codx * dvol / ncount
cody = cody * dvol / ncount
codz = codz * dvol / ncount

! get quadrupole tensor

DO i=1,3
  DO j=i,3
    qtel(i,j)=0D0
    
    ind=0
    DO iz=minz,maxz
      xic(3)=(iz-nzsh)*dz-codz
      DO iy=miny,maxy
        xic(2)=(iy-nysh)*dy-cody
        DO ix=minx,maxx
          xic(1)=(ix-nxsh)*dx-codx
          
          ind=ind+1
          sum = 3*xic(i)*xic(j)
          IF (i == j) THEN
            rr=0.
            DO k=1,3
              rr=xic(k)**2+rr
            END DO
            sum=sum-rr
          END IF
          
          qtel(i,j)=qtel(i,j)+sum*q1(ind)
          
        END DO
      END DO
    END DO
    
    qtel(i,j)=qtel(i,j)*dvol/ncount
    qtel(j,i)=qtel(i,j)
    
  END DO
END DO


! get root mean square radius

ind=0
rmsel=0D0
DO iz=minz,maxz
  xic(3)=(iz-nzsh)*dz-codz
  DO iy=miny,maxy
    xic(2)=(iy-nysh)*dy-cody
    DO ix=minx,maxx
      xic(1)=(ix-nxsh)*dx-codx
      
      ind=ind+1
      rr=0
      DO k=1,3
        rr=rr+xic(k)*xic(k)
      END DO
      rmsel = rmsel + q1(ind)*rr
      
    END DO
  END DO
END DO

rmsel = SQRT(rmsel*dvol/ncount)

DEALLOCATE(q1)


RETURN
END SUBROUTINE getelectrongeometry
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE evalprops
!------------------------------------------------------------
USE params
IMPLICIT NONE
!     evaluates properties at a given iteration, such like
!     kinetic energy of ions, cores, clouds and cations
!     as well as coupling energies etc.

INTEGER :: i
REAL(DP) :: kinenergy
#if(raregas)
REAL(DP) :: kinenergyc, kinenergye, kinenergyk
#endif
REAL(DP) :: vsum, vvsum, vsumc, vvsumc, vsume, vvsume, vsumk, vvsumk
REAL(DP) :: vx, vy, vz, xm

vsum = 0D0
vvsum = 0D0
vsumc = 0D0
vvsumc = 0D0
vsume = 0D0
vvsume = 0D0
vsumk = 0D0
vvsumk = 0D0



!     kinetic energies first

xm = amu(np(i)) ! is that correct??? check with ionmd.F

DO i=1,nion
  
! reset velocities by half a time step because of
! the use of leapfrog
  vx = (cpx(i) - 0.5D0*fx(i)*dt1)
  vy = (cpy(i) - 0.5D0*fy(i)*dt1)
  vz = (cpz(i) - 0.5D0*fz(i)*dt1)
  vsum = vsum + vx + vy + vz
  vvsum = vvsum + vx*vx + vy*vy + vz*vz
END DO

kinenergy = 0.5D0* vvsum/xm

#if(raregas)
xm = mion ! is that correct??? check with ionmd.F

DO i=1,nc
  
! reset velocities by half a time step because of
! the use of leapfrog
  vx = (pxc(i) - 0.5D0*fxc(i)*dt1)
  vy = (pyc(i) - 0.5D0*fyc(i)*dt1)
  vz = (pzc(i) - 0.5D0*fzc(i)*dt1)
  vsumc = vsumc + vx + vy + vz
  vvsumc = vvsumc + vx*vx + vy*vy + vz*vz
END DO

kinenergyc = 0.5D0* vvsumc/xm / nc

xm = me ! is that correct??? check with ionmd.F

DO i=1,NE
  
! reset velocities by half a time step because of
! the use of leapfrog
  vx = (pxe(i) - 0.5D0*fxe(i)*dt1)
  vy = (pye(i) - 0.5D0*fye(i)*dt1)
  vz = (pze(i) - 0.5D0*fze(i)*dt1)
  vsume = vsume + vx + vy + vz
  vvsume = vvsume + vx*vx + vy*vy + vz*vz
END DO

kinenergye = 0.5D0* vvsume/xm / NE

xm = mkat ! is that correct??? check with ionmd.F

DO i=1,nk
  
! reset velocities by half a time step because of
! the use of leapfrog
  vx = (pxk(i) - 0.5D0*fxk(i)*dt1)
  vy = (pyk(i) - 0.5D0*fyk(i)*dt1)
  vz = (pzk(i) - 0.5D0*fzk(i)*dt1)
  vsumk = vsumk + vx + vy + vz
  vvsumk = vvsumk + vx*vx + vy*vy + vz*vz
END DO

kinenergyk = 0.5D0* vvsumk/xm / nk
#endif

RETURN
END SUBROUTINE evalprops
!------------------------------------------------------------


!~ !------------------------------------------------------------

!~ SUBROUTINE accumprops(code)
!~ !------------------------------------------------------------
!~ USE params
!~ IMPLICIT NONE
!~ !     accumulates properties for averaging over successive
!~ !     iterations


!~ INTEGER, INTENT(IN) :: code

!~ SELECT CASE(code)
!~   CASE(0) ! (re-)set accum. variables to zero
!~   skinenergy = 0D0
!~   sskinenergy = 0D0
!~   skinenergyc = 0D0
!~   sskinenergyc = 0D0
!~   skinenergye = 0D0
!~   sskinenergye = 0D0
!~   skinenergyk = 0D0
!~   sskinenergyk = 0D0
!~   CASE(1) ! accumulate variables
!~   skinenergy =  skinenergy + kinenergy
!~   sskinenergy = sskinenergy + kinenergy*kinenergy
!~   skinenergyc =  skinenergyc + kinenergyc
!~   sskinenergyc = sskinenergyc + kinenergyc*kinenergyc
!~   skinenergye =  skinenergye + kinenergye
!~   sskinenergye = sskinenergye + kinenergye*kinenergye
!~   skinenergyk =  skinenergyk + kinenergyk
!~   sskinenergyk = sskinenergyk + kinenergyk*kinenergyk
!~   CASE(2) ! evaluate averages and variances
!~   skinenergy = skinenergy / istepavg
!~   sskinenergy = SQRT(sskinenergy / istepavg-skinenergy)
!~   skinenergyc = skinenergyc / istepavg
!~   sskinenergyc = SQRT(sskinenergyc / istepavg-skinenergyc)
!~   skinenergye = skinenergye / istepavg
!~   sskinenergye = SQRT(sskinenergye / istepavg-skinenergye)
!~   skinenergyk = skinenergyk / istepavg
!~   sskinenergyk = SQRT(sskinenergyk / istepavg-skinenergyk)
!~ ! now the ss... variables are the standard deviations
!~   CASE DEFAULT
!~   STOP 'Error in AccumProps: wrong argument'
!~ END SELECT

!~ RETURN
!~ END SUBROUTINE accumprops
!~ !------------------------------------------------------------


