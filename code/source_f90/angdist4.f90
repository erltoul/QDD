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

program angdist4
IMPLICIT REAL(KIND(1D0)) (A-H,O-Z)
INTEGER,PARAMETER                         ::  DP=KIND(1D0)
COMMON dx,kxbox,kybox,kzbox,nabsorb,ispherabso
! plots the angular distribution of ionization as density plot;
! requires pescmask file
! collects electron loss by using grid point represented by 
! box/tent distribution 

INTEGER,PARAMETER                         ::  nstepstheta = 80
INTEGER,PARAMETER                         ::  nstepsphi   = 80
REAL(DP),DIMENSION(nstepstheta,nstepsphi) ::  angdist
REAL(DP),DIMENSION(nstepstheta)           ::  angdistp,wt,theta

REAL(DP),ALLOCATABLE                      ::  rhoabso(:)
REAL(DP),PARAMETER                        ::  pi          = 3.141592653589793D0

! ---------------- Settings ----------------------------------------------------
dx          = 0.8D0    !dx          = 0.824187016040328D0
kxbox       = 64
kybox       = 64
kzbox       = 64
ispherabso  = 2        ! (0) cartesian (1) spherical (2) ellipsoidal mask
nabsorb     = 8
drho        = dx/8D0   ! experience value
!xo          = 0D0
!yo          = 0D0
!zo          = 0D0

deltatheta = pi/(nstepstheta-1)
kdfull2    = kxbox*kybox*kzbox
bcrad      = nabsorb*dx
ALLOCATE(rhoabso(kdfull2))
dmax2     = MAX(kxbox,kybox,kzbox)/2D0*dx
nmaxrho   = dmax2/drho

do nph=1,nstepsphi
do nth=1,nstepstheta
 angdist(nth,nph) = 0D0
 angdistp(nth)    = 0D0
enddo
enddo

! read in mask
do ind=1,kdfull2
 read(5,*) dum,dum,dum,rhoabso(ind)
enddo

write(6,'(1a,1i4)')   ' nmaxrho     = ',nmaxrho
write(6,'(1a,1i4)')   ' nstepstheta = ',nstepstheta
write(6,'(1a,1i4)')   ' nstepsphi   = ',nstepsphi
write(6,'(1a,1i4)')   ' nabsorb     = ',nabsorb
write(6,'(1a,1f8.4)') ' dx          = ',dx
write(6,'(1a,1f8.4)') ' dmax2       = ',dmax2
write(6,'(1a,1f8.4)') ' drho        = ',drho

open(28,status='unknown',file='pangdist4-1.res')
open(29,status='unknown',file='pangdist4-2.res')

! cccccccccccccccccccccccccccccccccccccccccccccccccc


! initialize theta and phi
do nth=1,nstepstheta
 theta(nth) = deltatheta*(nth-1)
 if(nth == 1 .OR. nth == nstepstheta) then
  wt(nth)   = 2D0*pi*(1D0 - cos(deltatheta/2D0))
 else
  wt(nth)   = 2D0*pi*(dcos(theta(nth)-deltatheta/2.)-dcos(theta(nth)+deltatheta/2.))
 endif
 wt(nth)    = dabs(wt(nth))
enddo


! calculate angdist
do nth=1,nstepstheta
do nph=1,nstepsphi
do nrh=1,nmaxrho
   phi   = 2D0*pi*(nph-1)/nstepsphi
   rho   = drho*nrh

   rvecx = rho*dcos(phi)*dsin(theta(nth))
   rvecy = rho*dsin(phi)*dsin(theta(nth))
   rvecz = rho*dcos(theta(nth))

   IF(rvecx > kxbox/2D0*dx .OR. rvecx < -(kxbox-2)/2D0*dx) goto 23  ! out of box
   IF(rvecy > kybox/2D0*dx .OR. rvecy < -(kybox-2)/2D0*dx) goto 23  ! out of box
   IF(rvecz > kzbox/2D0*dx .OR. rvecz < -(kzbox-2)/2D0*dx) goto 23  ! out of box
!   IF(icheckabsozone(rvecx,rvecy,rvecz) == 0) write(777,'(3f,i)') rvecx,rvecy,rvecz,0    ! out of abso zone
   IF(icheckabsozone(rvecx,rvecy,rvecz) == 0) goto 23               ! out of abso zone

!write(777,'(3f,i)') rvecx,rvecy,rvecz,1

ind = 0
do iz=1,kzbox
do iy=1,kybox
do ix=1,kxbox
 z1=(iz-kzbox/2)*dx
 y1=(iy-kybox/2)*dx
 x1=(ix-kxbox/2)*dx

 ind = ind + 1

 rxrel = x1!-xo
 ryrel = y1!-xo
 rzrel = z1!-xo

 rxdiff=rvecx-rxrel
 rydiff=rvecy-ryrel
 rzdiff=rvecz-rzrel

 func=gtent(rxdiff,rydiff,rzdiff)  !  tent function
 angdist(nth,nph)=angdist(nth,nph)+rhoabso(ind)*rho**2D0*func*drho
enddo
enddo
enddo

23 continue
enddo   !  nmaxrho
 if(nth == 1 .OR. nth == nstepstheta) exit 
enddo   !  nstepsphi
enddo   !  nstepstheta


! write out
do nth=1,nstepstheta
do nph=1,nstepsphi
 if(nth == 1)           angdist(nth,nph) = angdist(1,1)
 if(nth == nstepstheta) angdist(nth,nph) = angdist(nstepstheta,1)
 write(28,'(2f17.4,1e17.7)') (nph-1)*360D0/nstepsphi,theta(nth)/pi*180D0,angdist(nth,nph)
enddo
 write(28,'(2f17.4,1e17.7)') 360D0,theta(nth)/pi*180D0,angdist(nth,1)
 write(28,*)
enddo


! calculate phi average and print out pangdist3-2
total=0.0
wg   =0.0
do nth=1,nstepstheta
do nph=1,nstepsphi
 angdistp(nth)=angdistp(nth)+angdist(nth,nph)/nstepsphi
enddo
 total=total+angdistp(nth)*wt(nth)
 wg=wg+wt(nth)
enddo

do nth=1,nstepstheta
 write(29,'(1f17.4,2e17.7)') theta(nth)/pi*180D0,angdistp(nth)/(total/4D0/pi), &
                                                angdistp(nth)/4D0/pi
enddo


! finish
write(*,'(1a,1e17.7)')  ' total       = ',total
write(*,'(1a,1f17.14)') ' wg          = ',wg
stop " finished"

end

! cccccccccccccccccccccccccccccccccccccccccccccccccc

REAL FUNCTION gbox(xx,yy,zz)
IMPLICIT REAL(KIND(1D0)) (A-H,O-Z)
COMMON dx
gbox = 1D0
if(abs(xx).gt.dx/2D0) gbox=0D0
if(abs(yy).gt.dx/2D0) gbox=0D0
if(abs(zz).gt.dx/2D0) gbox=0D0
end function

REAL FUNCTION gtent(xx,yy,zz)
IMPLICIT REAL(KIND(1D0)) (A-H,O-Z)
COMMON dx
gtent=max(dx-dabs(xx),0D0)*max(dx-dabs(yy),0D0)*max(dx-dabs(zz),0D0)
gtent=gtent/dx/dx/dx
end function

! cccccccccccccccccccccccccccccccccccccccccccccccccc

INTEGER FUNCTION icheckabsozone(xx,yy,zz)
IMPLICIT REAL(KIND(1D0)) (A-H,O-Z)
COMMON dx,kxbox,kybox,kzbox,nabsorb,ispherabso

icheckabsozone = 0

IF(ispherabso == 0) THEN
if(xx < (-(kxbox-2)/2D0 + nabsorb)*dx .OR. xx > (kxbox/2D0 - nabsorb)*dx) icheckabsozone = 1
if(yy < (-(kybox-2)/2D0 + nabsorb)*dx .OR. yy > (kybox/2D0 - nabsorb)*dx) icheckabsozone = 1
if(zz < (-(kzbox-2)/2D0 + nabsorb)*dx .OR. zz > (kzbox/2D0 - nabsorb)*dx) icheckabsozone = 1
ELSE IF(ispherabso == 1) THEN
 dmin2 = MIN(kxbox,kybox,kzbox)/2D0*dx
 dmin1 = dmin2 - nabsorb*dx
 dist2 = xx*xx + yy*yy + zz*zz
 dmin12 = dmin1*dmin1
 dmin22 = dmin2*dmin2
 if(dist2 > dmin12 .AND. dist2 < dmin22) icheckabsozone = 1
ELSE IF(ispherabso == 2) THEN
 dmin2  = MIN(kxbox,kybox,kzbox)/2D0*dx
 dmin2x = kxbox/2D0*dx
 dmin2y = kybox/2D0*dx
 dmin2z = kzbox/2D0*dx
 dmin1x = dmin2x - dmin2x/dmin2*nabsorb*dx
 dmin1y = dmin2y - dmin2y/dmin2*nabsorb*dx
 dmin1z = dmin2z - dmin2z/dmin2*nabsorb*dx
 dmin12x = dmin1x*dmin1x
 dmin12y = dmin1y*dmin1y
 dmin12z = dmin1z*dmin1z
 dmin22x = dmin2x*dmin2x
 dmin22y = dmin2y*dmin2y
 dmin22z = dmin2z*dmin2z
 ellips1 = xx*xx/dmin12x+yy*yy/dmin12y+zz*zz/dmin12z
 ellips2 = xx*xx/dmin22x+yy*yy/dmin22y+zz*zz/dmin22z
 if(ellips1 > 1D0 .AND. ellips2 <= 1D0) icheckabsozone = 1  
ENDIF

RETURN
END FUNCTION








!ELSE IF(ispherabso == 1) THEN
! dmin2 = MIN(kxbox,kybox,kzbox)/2D0*dx
! dmin1 = dmin2 - nabsorb*dx
! dist2 = xx*xx + yy*yy + zz*zz
!! circ1 = xx*xx/dmin1/dmin1 + yy*yy/dmin1/dmin1 + zz*zz/dmin1/dmin1
!! circ2 = xx*xx/dmin2/dmin2 + yy*yy/dmin2/dmin2 + zz*zz/dmin2/dmin2
!! if(circ1 > 1D0 .AND. circ2 <= 1D0) icheckabsozone = 1
! dmin12 = dmin1*dmin1
! dmin22 = dmin2*dmin2
! if(dist2*dist2*dist2 > dmin12*dmin12*dmin12 .AND. dist2*dist2*dist2 <= dmin22*dmin22*dmin22) icheckabsozone = 1
! circ1  = xx*xx*dmin12*dmin12+yy*yy*dmin12*dmin12+zz*zz*dmin12*dmin12
! circ2  = xx*xx*dmin22*dmin22+yy*yy*dmin22*dmin22+zz*zz*dmin22*dmin22
!! if(circ1 > dmin12*dmin12*dmin12 .AND. circ1 <= dmin22*dmin22*dmin22) icheckabsozone = 1
!
! if(dist2 > dmin12 .AND. dist2 < dmin22) icheckabsozone = 1
!ELSE IF(ispherabso == 2) THEN
! dmin2  = MIN(kxbox,kybox,kzbox)/2D0*dx
! dmin2x = kxbox/2D0*dx
! dmin2y = kybox/2D0*dx
! dmin2z = kzbox/2D0*dx
! dmin1x = dmin2x - dmin2x/dmin2*nabsorb*dx
! dmin1y = dmin2y - dmin2y/dmin2*nabsorb*dx
! dmin1z = dmin2z - dmin2z/dmin2*nabsorb*dx
!! ellips1 = xx*xx/dmin1x/dmin1x+yy*yy/dmin1y/dmin1y+zz*zz/dmin1z/dmin1z
!! ellips2 = xx*xx/dmin2x/dmin2x+yy*yy/dmin2y/dmin2y+zz*zz/dmin2z/dmin2z
!! if(ellips1 > 1D0 .AND. ellips2 <= 1D0) icheckabsozone = 1
! dmin12x = dmin1x*dmin1x
! dmin12y = dmin1y*dmin1y
! dmin12z = dmin1z*dmin1z
! dmin22x = dmin2x*dmin2x
! dmin22y = dmin2y*dmin2y
! dmin22z = dmin2z*dmin2z
! ellips1 = xx*xx*dmin12y*dmin12z+yy*yy*dmin12x*dmin12z+zz*zz*dmin12x*dmin12y
! ellips2 = xx*xx*dmin22y*dmin22z+yy*yy*dmin22x*dmin22z+zz*zz*dmin22x*dmin22y
! if(ellips1 > dmin12x*dmin12y*dmin12z .AND. ellips2 <= dmin22x*dmin22y*dmin22z) icheckabsozone = 1  
!ENDIF
!
!RETURN
!END FUNCTION

