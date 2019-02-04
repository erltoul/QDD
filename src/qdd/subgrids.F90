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
!------------------------------------------------------------

!SUBROUTINE addfields2(field1,field2,factor2)
!!------------------------------------------------------------
!USE params
!IMPLICIT REAL(DP) (A-H,O-Z)
!
!
!
!REAL(DP), INTENT(IN OUT)                        :: field1(kdfull2)
!REAL(DP), INTENT(IN)                         :: field2(kdfull2)
!REAL(DP), INTENT(IN)                         :: factor2
!
!
!DO ind=1,kdfull2
!  field1(ind)=field1(ind)+factor2*field2(ind)
!END DO
!
!RETURN
!END SUBROUTINE addfields2
!------------------------------------------------------------

#if(raregas)
!------------------------------------------------------------

!SUBROUTINE addfields(field1,field2,isignum)
!!------------------------------------------------------------
!USE params
!IMPLICIT REAL(DP) (A-H,O-Z)
!
!
!
!REAL(DP), INTENT(IN)                         :: field1(kdfull2)
!REAL(DP), INTENT(IN)                         :: field2(kdfull2)
!INTEGER, INTENT(IN)                      :: isignum
!
!
!DO ind=1,kdfull2
!  rfieldtmp(ind)=field1(ind)+(-1D0)**isignum*field2(ind)
!END DO
!
!RETURN
!END SUBROUTINE addfields




!------------------------------------------------------------

SUBROUTINE initsubgrids
!------------------------------------------------------------
USE params
IMPLICIT NONE

INTEGER :: i
!     give every particle a subgrid which surrounds it
DO i=1,nc
  CALL makesubgrid(i,xc(i),yc(i),zc(i))
!         call debuggrid
!         stop 'aslkdj'
END DO
DO i=1,ne
  CALL makesubgrid(i+nc,xe(i),ye(i),ze(i))
END DO
DO i=1,nk
  CALL makesubgrid(i+nc+NE,xk(i),yk(i),zk(i))
END DO

END SUBROUTINE initsubgrids
!------------------------------------------------------------
#endif

!------------------------------------------------------------

SUBROUTINE addfunctofield(field,func,x0,y0,z0,fact)
!------------------------------------------------------------
USE params
IMPLICIT NONE
!     adds function func(x-x0,y-y0,z-z0) to field(x,y,z)
!     version with zero parameters in func

REAL(DP), INTENT(IN OUT)                        :: field(kdfull2)
REAL(DP), INTENT(IN)                         :: x0
REAL(DP), INTENT(IN)                         :: y0
REAL(DP), INTENT(IN)                         :: z0
REAL(DP), INTENT(IN)                         :: fact
! Function name passed as an argument :
INTERFACE 
  REAL(DP) FUNCTION func(r)
    USE params
    REAL(DP),INTENT(IN) :: r
  END FUNCTION func
END INTERFACE

INTEGER :: ind, ix, iy, iz, max 
REAL(DP) :: rr, rx, ry, rz, x1, y1, z1

ind = 0

DO iz=minz,maxz
  z1=(iz-nzsh)*dz
  rz=z1-z0
  DO iy=miny,maxy
    y1=(iy-nysh)*dy
    ry=y1-y0
    DO ix=minx,maxx
      x1=(ix-nxsh)*dx
      rx=x1-x0
      rr=SQRT(rx*rx+ry*ry+rz*rz)
      
      rr=MAX(rr,small)
      
      ind = ind + 1
      
      field(ind) = field(ind) + func(rr)*fact
      
    END DO
  END DO
END DO


RETURN
END SUBROUTINE addfunctofield
!------------------------------------------------------------


!------------------------------------------------------------

SUBROUTINE addfunctofield1(field,func,x0,y0,z0,fact,para)
!------------------------------------------------------------
USE params
IMPLICIT NONE
!     adds function func(x-x0,y-y0,z-z0) to field(x,y,z)
!     version with 1 parameter in func

REAL(DP), INTENT(IN OUT)                        :: field(kdfull2)
REAL(DP), INTENT(IN)                         :: x0
REAL(DP), INTENT(IN)                         :: y0
REAL(DP), INTENT(IN)                         :: z0
REAL(DP), INTENT(IN)                         :: fact
REAL(DP), INTENT(IN)                     :: para
! Function name passed as an argument :
INTERFACE 
  REAL(DP) FUNCTION func(r,s)
    USE params
    REAL(DP),INTENT(IN) :: r
    REAL(DP),INTENT(IN) :: s
  END FUNCTION func
END INTERFACE

INTEGER :: ind, ix, iy, iz
REAL(DP) ::rr, rx, ry, rz, x1, y1, z1 

ind = 0

DO iz=minz,maxz
  z1=(iz-nzsh)*dz
  rz=z1-z0
  DO iy=miny,maxy
    y1=(iy-nysh)*dy
    ry=y1-y0
    DO ix=minx,maxx
      x1=(ix-nxsh)*dx
      rx=x1-x0
      rr=SQRT(rx*rx+ry*ry+rz*rz)
      
      rr=MAX(rr,small)
      
      ind = ind + 1
      
      field(ind) = field(ind) + func(rr,para)*fact
      
    END DO
  END DO
END DO


RETURN
END SUBROUTINE addfunctofield1
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE addfunctofield3(field,func,x0,y0,z0,fact,para1,para2, para3)
!------------------------------------------------------------
USE params
IMPLICIT NONE
!     adds function func(x-x0,y-y0,z-z0) to field(x,y,z)
!     version with 3 parameter in func


REAL(DP), INTENT(IN OUT)           :: field(kdfull2)
REAL(DP), INTENT(IN)               :: x0
REAL(DP), INTENT(IN)               :: y0
REAL(DP), INTENT(IN)               :: z0
REAL(DP), INTENT(IN)               :: fact
REAL(DP), INTENT(IN)               :: para1
REAL(DP), INTENT(IN)               :: para2
REAL(DP), INTENT(IN)               :: para3
! Function name passed as an argument :
INTERFACE 
  REAL(DP) FUNCTION func(r,a,b,c)
    USE params
    REAL(DP), INTENT(IN) :: r
    REAL(DP), INTENT(IN) :: a
    REAL(DP), INTENT(IN) :: b
    REAL(DP), INTENT(IN) :: c
  END FUNCTION func
END INTERFACE

INTEGER :: ind, ix, iy, iz
REAL(DP) :: rr, rx, ry, rz, x1, y1, z1
ind = 0

DO iz=minz,maxz
  z1=(iz-nzsh)*dz
  rz=z1-z0
  DO iy=miny,maxy
    y1=(iy-nysh)*dy
    ry=y1-y0
    DO ix=minx,maxx
      x1=(ix-nxsh)*dx
      rx=x1-x0
      rr=SQRT(rx*rx+ry*ry+rz*rz)
      
      rr=MAX(rr,small)
      
      ind = ind + 1
      
      field(ind) = field(ind) + func(rr,para1,para2, para3)*fact
      
    END DO
  END DO
END DO


RETURN
END SUBROUTINE addfunctofield3
!------------------------------------------------------------


!------------------------------------------------------------

SUBROUTINE addfunctofieldonsubgrid(field,func,x0,y0,z0,fact, nsgsize)
!------------------------------------------------------------
USE params
IMPLICIT NONE
!     same as above routine but with subgrids

REAL(DP), INTENT(IN OUT)                        :: field(kdfull2)
REAL(DP), INTENT(IN)                         :: x0
REAL(DP), INTENT(IN)                         :: y0
REAL(DP), INTENT(IN)                         :: z0
REAL(DP), INTENT(IN)                     :: fact
INTEGER, INTENT(IN)                      :: nsgsize
! Function name passed as an argument :
INTERFACE 
  REAL(DP) FUNCTION func(r)
      USE params
      REAL(DP),INTENT(IN) :: r
  END FUNCTION func
END INTERFACE

INTEGER :: ind, ix, iy, iz
REAL(DP) :: rr, rx, ry, rz, x1, y1, z1

INTEGER,EXTERNAL :: getnearestgridpoint
INTEGER,EXTERNAL :: conv3to1


ind = getnearestgridpoint(x0,y0,z0)


CALL conv1to3(ind)

DO iz=iindtmp(3)-nsgsize,iindtmp(3)+nsgsize
  IF (iz >= minz .AND. iz <= maxz) THEN
    z1=(iz-nzsh)*dz
    rz=z1-z0
    
    
    DO iy=iindtmp(2)-nsgsize,iindtmp(2)+nsgsize
      IF (iy >= miny .AND. iy <= maxy) THEN
        y1=(iy-nysh)*dy
        ry=y1-y0
        
        DO ix=iindtmp(1)-nsgsize,iindtmp(1)+nsgsize
          
          IF (ix >= minx .AND. ix <= maxx) THEN
            
            x1=(ix-nxsh)*dx
            rx=x1-x0
            
            rr=SQRT(rx*rx+ry*ry+rz*rz)
            
            rr=MAX(rr,small) ! avoid zero
            
            ind=conv3to1(ix,iy,iz)
            
            
            
            field(ind) = field(ind) + func(rr)*fact
            
            
            
          END IF
        END DO
      END IF
    END DO
  END IF
  
END DO


RETURN
END SUBROUTINE addfunctofieldonsubgrid
!------------------------------------------------------------


!------------------------------------------------------------

SUBROUTINE addfunctofieldonsubgrid1(field,func,x0,y0,z0,fact, para,nsgsize)
!------------------------------------------------------------
USE params
IMPLICIT NONE
!     same as above routine but with subgrids


REAL(DP), INTENT(IN OUT)                        :: field(kdfull2)
REAL(DP), INTENT(IN)                         :: x0
REAL(DP), INTENT(IN)                         :: y0
REAL(DP), INTENT(IN)                         :: z0
REAL(DP), INTENT(IN)                     :: fact
REAL(DP), INTENT(IN)                     :: para
INTEGER, INTENT(IN)                      :: nsgsize
! Function name passed as an argument :
INTERFACE 
  REAL(DP) FUNCTION func(r,s)
    USE params
    REAL(DP),INTENT(IN) :: r
    REAL(DP),INTENT(IN) :: s
  END FUNCTION func
END INTERFACE

INTEGER :: ind, ix, iy, iz
REAL(DP) :: rr, rx, ry, rz, x1, y1, z1

INTEGER,EXTERNAL :: getnearestgridpoint
INTEGER,EXTERNAL :: conv3to1


ind = getnearestgridpoint(x0,y0,z0)


CALL conv1to3(ind)

DO iz=iindtmp(3)-nsgsize,iindtmp(3)+nsgsize
  IF (iz >= minz .AND. iz <= maxz) THEN
    z1=(iz-nzsh)*dz
    rz=z1-z0
    
    
    DO iy=iindtmp(2)-nsgsize,iindtmp(2)+nsgsize
      IF (iy >= miny .AND. iy <= maxy) THEN
        y1=(iy-nysh)*dy
        ry=y1-y0
        
        DO ix=iindtmp(1)-nsgsize,iindtmp(1)+nsgsize
          
          IF (ix >= minx .AND. ix <= maxx) THEN
            
            x1=(ix-nxsh)*dx
            rx=x1-x0
            
            rr=SQRT(rx*rx+ry*ry+rz*rz)
            
            rr=MAX(rr,small) ! avoid zero
            
            ind=conv3to1(ix,iy,iz)
            
            
            
            field(ind) = field(ind) + func(rr,para)*fact
            
            
            
          END IF
        END DO
      END IF
    END DO
  END IF
  
END DO


RETURN
END SUBROUTINE addfunctofieldonsubgrid1
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE addfunctofieldonsubgrid3(field,func,x0,y0,z0,fact,  &
    para1,para2,para3,nsgsize)
!------------------------------------------------------------
USE params
IMPLICIT NONE
!     same as above routine but with subgrids


REAL(DP), INTENT(IN OUT)                        :: field(kdfull2)
REAL(DP), INTENT(IN)                         :: x0
REAL(DP), INTENT(IN)                         :: y0
REAL(DP), INTENT(IN)                         :: z0
REAL(DP), INTENT(IN)                     :: fact
REAL(DP), INTENT(IN)                     :: para1
REAL(DP), INTENT(IN)                     :: para2
REAL(DP), INTENT(IN)                     :: para3
INTEGER, INTENT(IN)                      :: nsgsize
! Function name passed as an argument :
INTERFACE 
  REAL(DP) FUNCTION func(r,a,b,c)
    USE params
    REAL(DP), INTENT(IN) :: r
    REAL(DP), INTENT(IN) :: a
    REAL(DP), INTENT(IN) :: b
    REAL(DP), INTENT(IN) :: c
  END FUNCTION func
END INTERFACE

INTEGER :: ind, ix, iy, iz
REAL(DP) :: rr, rx, ry, rz, x1, y1, z1


INTEGER,EXTERNAL :: getnearestgridpoint
INTEGER,EXTERNAL :: conv3to1


ind = getnearestgridpoint(x0,y0,z0)


CALL conv1to3(ind)

DO iz=iindtmp(3)-nsgsize,iindtmp(3)+nsgsize
  IF (iz >= minz .AND. iz <= maxz) THEN
    z1=(iz-nzsh)*dz
    rz=z1-z0
    
    
    DO iy=iindtmp(2)-nsgsize,iindtmp(2)+nsgsize
      IF (iy >= miny .AND. iy <= maxy) THEN
        y1=(iy-nysh)*dy
        ry=y1-y0
        
        DO ix=iindtmp(1)-nsgsize,iindtmp(1)+nsgsize
          
          IF (ix >= minx .AND. ix <= maxx) THEN
            
            x1=(ix-nxsh)*dx
            rx=x1-x0
            
            rr=SQRT(rx*rx+ry*ry+rz*rz)
            
            rr=MAX(rr,small) ! avoid zero
            
            ind=conv3to1(ix,iy,iz)
            
            
            
            field(ind) = field(ind) + func(rr,para1,para2,para3)*fact
            
            
            
          END IF
        END DO
      END IF
    END DO
  END IF
  
END DO


RETURN
END SUBROUTINE addfunctofieldonsubgrid3
!------------------------------------------------------------

#if(raregas)
!------------------------------------------------------------

SUBROUTINE addtabtofield(field,table,x0,y0,z0,fact)
!------------------------------------------------------------
USE params
IMPLICIT NONE
!     Adds function f(sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2))
!     to field(x,y,z).
!     The function is interpolated from tabulated values
!     in array 'table'.
REAL(DP), INTENT(IN OUT)                        :: field(kdfull2)
REAL(DP), INTENT(IN)                     :: table(kfermi)
REAL(DP), INTENT(IN)                         :: x0
REAL(DP), INTENT(IN)                         :: y0
REAL(DP), INTENT(IN)                         :: z0
REAL(DP), INTENT(IN)                         :: fact

INTEGER :: ind, ir, ix, iy, iz
REAL(DP) :: deltar, deltarinv, rcon, rrel
REAL(DP) :: rr, rx, ry, rz, rx2, ry2, rz2, x1, y1, z1

deltar = dx/1733D0

ind = 0
deltarinv = 1D0/deltar

DO iz=minz,maxz
  z1=(iz-nzsh)*dz
  rz2=z1-z0
  rz2=rz2*rz2
  DO iy=miny,maxy
    y1=(iy-nysh)*dy
    ry2=y1-y0
    ry2=ry2*ry2+rz2
    DO ix=minx,maxx
      x1=(ix-nxsh)*dx
      rx=x1-x0
      rr=SQRT(rx*rx+ry2)
      
!                      rr=max(rr,SMALL)
      
      rrel = rr*deltarinv
      ir = INT(rrel)
      rrel = fact*(rrel-ir)
      rcon = fact-rrel
      ir = ir+1
      
      ind = ind + 1
      
      field(ind) = field(ind) + rcon*table(ir)+rrel*table(ir+1)
      
    END DO
  END DO
END DO

RETURN
END SUBROUTINE addtabtofield
!------------------------------------------------------------


!------------------------------------------------------------

SUBROUTINE addtabtofieldonsubgrid(field,table,x0,y0,z0, fact,nsgsize)
!------------------------------------------------------------
USE params
IMPLICIT NONE
!     same as above routine but with subgrids


REAL(DP), INTENT(IN OUT)                         :: field(kdfull2)
REAL(DP), INTENT(IN)                         :: table(kfermi)
REAL(DP), INTENT(IN)                         :: x0
REAL(DP), INTENT(IN)                         :: y0
REAL(DP), INTENT(IN)                         :: z0
REAL(DP), INTENT(IN)                         :: fact
INTEGER, INTENT(IN)                      :: nsgsize

INTEGER :: ind, ir, ix, iy, iz
REAL(DP) :: deltar, deltarinv, rcon, rrel
REAL(DP) :: rr, rx, ry, rz, rx2, ry2, rz2, x1, y1, z1

INTEGER :: getnearestgridpoint
INTEGER :: conv3to1

deltar = dx/1733D0
deltarinv = 1D0/deltar


ind = getnearestgridpoint(x0,y0,z0)


CALL conv1to3(ind)

DO iz=iindtmp(3)-nsgsize,iindtmp(3)+nsgsize
  IF (iz >= minz .AND. iz <= maxz) THEN
    z1=(iz-nzsh)*dz
    rz2=z1-z0
    rz2=rz2*rz2
    
    
    DO iy=iindtmp(2)-nsgsize,iindtmp(2)+nsgsize
      IF (iy >= miny .AND. iy <= maxy) THEN
        y1=(iy-nysh)*dy
        ry2=y1-y0
        ry2=ry2*ry2+rz2
        
        DO ix=iindtmp(1)-nsgsize,iindtmp(1)+nsgsize
          
          IF (ix >= minx .AND. ix <= maxx) THEN
            
            x1=(ix-nxsh)*dx
            rx=x1-x0
            
            rr=SQRT(rx*rx+ry2)
            
!                      rr=max(rr,SMALL) ! avoid zero
            
            rrel = rr*deltarinv
            ir = INT(rrel)
            rrel = fact*(rrel-ir)
            rcon = fact-rrel
            ir = ir+1
            
            ind=conv3to1(ix,iy,iz)
            
            
            
            field(ind) = field(ind) + rcon*table(ir)+rrel*table(ir+1)
            
            
            
          END IF
        END DO
      END IF
    END DO
  END IF
  
END DO


RETURN
END SUBROUTINE addtabtofieldonsubgrid
!------------------------------------------------------------
#endif


!------------------------------------------------------------

SUBROUTINE dintfieldfunc(field,func,xx,yy,zz,param)
!------------------------------------------------------------
USE params
IMPLICIT NONE

!     folds function around xx,yy,zz

!     calculates the integral

!     \int d^3r field(\vec r) \times \nabla_R func(|\vec r-\vec R|)

!     where \vec R = (xx,yy,zz)

REAL(DP), INTENT(IN)                         :: field(kdfull2)
REAL(DP), INTENT(IN)                         :: xx
REAL(DP), INTENT(IN)                         :: yy
REAL(DP), INTENT(IN)                         :: zz
REAL(DP), INTENT(IN)         :: param
! Function name passed as an argument :
INTERFACE 
  REAL(DP) FUNCTION func(r,s)
    USE params
    REAL(DP), INTENT(IN) :: r
    REAL(DP), INTENT(IN) :: s
  END FUNCTION func
END INTERFACE

INTEGER :: ind, ix, iy, iz
REAL(DP) :: sumx, sumy, sumz, radial, rr, rx, ry, rz, x1, y1, z1




ind = 0
sumx = 0D0
sumy = 0D0
sumz = 0D0


DO iz=minz,maxz
  z1=(iz-nzsh)*dz
  rz=zz-z1
  DO iy=miny,maxy
    y1=(iy-nysh)*dy
    ry=yy-y1
    DO ix=minx,maxx
      x1=(ix-nxsh)*dx
      rx=xx-x1
      rr=SQRT(rx*rx+ry*ry+rz*rz)
      rr=MAX(rr,small) ! avoid zero
      ind=ind+1
      
      radial = func(rr,param)/rr
      
      sumx = sumx + field(ind)*radial*rx
      sumy = sumy + field(ind)*radial*ry
      sumz = sumz + field(ind)*radial*rz
      
    END DO
  END DO
END DO



rvectmp(1) = sumx * dvol
rvectmp(2) = sumy * dvol
rvectmp(3) = sumz * dvol


RETURN
END SUBROUTINE dintfieldfunc

!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE dintfieldfunconsubgrid(field,func,xx,yy,zz,param)
!------------------------------------------------------------
USE params
IMPLICIT NONE

!     folds function around xx,yy,zz using a subgrid only

!     calculates the integral

!     \int d^3r field(\vec r) \times \nable_R func(|\vec r-\vec R|)

!     where \vec R = (xx,yy,zz)

REAL(DP), INTENT(IN)                         :: field(kdfull2)
! Function name passed as an argument :
INTERFACE 
  REAL(DP) FUNCTION func(r,s)
    USE params
    REAL(DP), INTENT(IN) :: r
    REAL(DP), INTENT(IN) :: s
  END FUNCTION func
END INTERFACE

REAL(DP), INTENT(IN)                         :: xx
REAL(DP), INTENT(IN)                         :: yy
REAL(DP), INTENT(IN)                         :: zz
REAL(DP), INTENT(IN)         :: param

INTEGER :: ind, ix, iy, iz
REAL(DP) :: sumx, sumy, sumz, radial, rr, rx, ry, rz, x1, y1, z1

INTEGER,EXTERNAL :: getnearestgridpoint
INTEGER,EXTERNAL :: conv3to1


ind = getnearestgridpoint(xx,yy,zz)

CALL conv1to3(ind)


sumx = 0D0
sumy = 0D0
sumz = 0D0


DO iz=iindtmp(3)-nzsg,iindtmp(3)+nzsg
  z1=(iz-nzsh)*dz
  rz=zz-z1
  
  DO iy=iindtmp(2)-nysg,iindtmp(2)+nysg
    y1=(iy-nysh)*dy
    ry=yy-y1
    
    DO ix=iindtmp(1)-nxsg,iindtmp(1)+nxsg
      x1=(ix-nxsh)*dx
      rx=xx-x1
      
      rr=SQRT(rx*rx+ry*ry+rz*rz)
      rr=MAX(rr,small)
      
      ind = conv3to1(ix,iy,iz)
      
      radial = func(rr,param)/rr
      
      sumx = sumx + field(ind)*radial*rx
      sumy = sumy + field(ind)*radial*ry
      sumz = sumz + field(ind)*radial*rz
      
    END DO
  END DO
END DO
rvectmp(1) = sumx * dvol
rvectmp(2) = sumy * dvol
rvectmp(3) = sumz * dvol


RETURN
END SUBROUTINE dintfieldfunconsubgrid




!GB




!------------------------------------------------------------

#if(raregas)
!------------------------------------------------------------

SUBROUTINE makesubgrid(ipart,xpos,ypos,zpos)
!------------------------------------------------------------
USE params
IMPLICIT NONE


INTEGER, INTENT(IN)                  :: ipart
REAL(DP), INTENT(IN)                     :: xpos
REAL(DP), INTENT(IN)                     :: ypos
REAL(DP), INTENT(IN)                     :: zpos

INTEGER,EXTERNAL :: getnearestgridpoint

! get the index of the grid point closest to the exact position

isubgcenter(ipart) = getnearestgridpoint(xpos,ypos,zpos)


END SUBROUTINE makesubgrid
!------------------------------------------------------------
#endif
!------------------------------------------------------------

SUBROUTINE updatesubgrids
!------------------------------------------------------------
USE params
IMPLICIT NONE

! if a particle has moved so far that we have to move the
! subgrid as well, then do so!
#if(raregas)
INTEGER :: i, ind
INTEGER,EXTERNAL :: getnearestgridpoint


DO i=1,nc
  ind = getnearestgridpoint(xc(i),yc(i),zc(i))
  IF (ind /= isubgcenter(i)) THEN
!            call makeSubGrid(i,xc(i),yc(i),zc(i))
    isubgcenter(i) = ind
  END IF
END DO
DO i=1,NE
  ind = getnearestgridpoint(xe(i),ye(i),ze(i))
  IF (ind /= isubgcenter(i+nc)) THEN
!            call makeSubGrid(i,xe(i),ye(i),ze(i))
    isubgcenter(i+nc) = ind
  END IF
END DO
DO i=1,nk
  ind = getnearestgridpoint(xk(i),yk(i),zk(i))
  IF (ind /= isubgcenter(i+NE+nc))  &
      THEN
!            call makeSubGrid(i,xk(i),yk(i),zk(i))
  isubgcenter(i+nc+NE) = ind
END IF
END DO
#endif

!$$$c now we check if the particles are in the box
!$$$      do i=1,nc+ne+nk
!$$$         call conv1To3(isubgcenter(i))
!$$$         iz = iIndTmp(3)
!$$$         iy = iIndTmp(2)
!$$$         ix = iIndTmp(1)
!$$$
!$$$         if (iz.ge.minz+nzsg .and. iz.le.maxz-nzsg .and.
!$$$     &       iy.ge.miny+nysg .and. iy.le.maxy-nysg .and.
!$$$     &       ix.ge.minx+nxsg .and. ix.le.maxx-nxsg) then
!$$$
!$$$            iOutOfBox(i) = 0
!$$$
!$$$         else
!$$$            iOutOfBox(i) = 1
!$$$
!$$$         endif
!$$$
!$$$      enddo
!$$$
!$$$
!$$$      do i=1,nc
!$$$         if (iOutOfBox(i).eq.0 .or. iOutOfBox(i+nc).eq.0) then
!$$$            iOutOfBox(i) = 0
!$$$            iOutOfBox(i+nc) = 0
!$$$         endif
!$$$      enddo
!$$$
!$$$


END SUBROUTINE updatesubgrids
!------------------------------------------------------------
#if(raregas)

!------------------------------------------------------------

SUBROUTINE putfunctosubgrid(qfield,indpart,func)
!------------------------------------------------------------
USE params
IMPLICIT NONE

REAL(DP), INTENT(OUT)                        :: qfield(kdsub)
INTEGER, INTENT(IN)                      :: indpart
! Function name passed as an argument :
INTERFACE 
  REAL(DP) FUNCTION func(r,s)
    USE params
    REAL(DP), INTENT(IN) :: r
    REAL(DP), INTENT(IN) :: s
  END FUNCTION func
END INTERFACE

INTEGER :: ind, ix, iy, iz, ix0, iy0, iz0
REAL(DP) :: r, ss, x1, y1, z1, x0, y0, z0


CALL conv1to3(isubgcenter(indpart))

iz0 = iindtmp(3)
iy0 = iindtmp(2)
ix0 = iindtmp(1)


ind = 0
DO iz=iz0-nzsg,iz0+nzsg
  z1 = (iz-nz)*dz
  DO iy=iy0-nysg,iy0+nysg
    y1 = (iy-ny)*dy
    DO ix=ix0-nxsg,ix0+nxsg
      x1 = (ix-nx)*dx
      
      ind = ind + 1
      
      IF (indpart <= nc) THEN
        x0 = xc(indpart)
        y0 = yc(indpart)
        z0 = zc(indpart)
        ss = sigmac
      ELSE IF (indpart <= nc+NE) THEN
        x0 = xe(indpart)
        y0 = ye(indpart)
        z0 = ze(indpart)
        ss = sigmav
      ELSE
        x0 = xk(indpart)
        y0 = yk(indpart)
        z0 = zk(indpart)
        ss = sigmak
      END IF
      
      r = (x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0)
      r = SQRT(r)
      
      qfield(ind) = func(r,ss)
      
    END DO
  END DO
END DO


END SUBROUTINE putfunctosubgrid
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE putgausstosubgrid(isgcenter,subgfield,x0,y0,z0,s)
!------------------------------------------------------------
USE params
IMPLICIT NONE

INTEGER, INTENT(IN)              :: isgcenter
REAL(DP), INTENT(OUT)            :: subgfield(kdsub)
REAL(DP), INTENT(IN)             :: x0
REAL(DP), INTENT(IN)             :: y0
REAL(DP), INTENT(IN)             :: z0
REAL(DP), INTENT(IN)             :: s


INTEGER :: ind, ix, iy, iz, ix0, iy0, iz0
REAL(DP) :: r2, x1, y1, z1

CALL conv1to3(isgcenter)

iz0 = iindtmp(3)
iy0 = iindtmp(2)
ix0 = iindtmp(1)


ind = 0
DO iz=iz0-nzsg,iz0+nzsg
  z1 = (iz-nz)*dz
  DO iy=iy0-nysg,iy0+nysg
    y1 = (iy-ny)*dy
    DO ix=ix0-nxsg,ix0+nxsg
      x1 = (ix-nx)*dx
      
      ind = ind + 1
      
      r2 = (x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0)
      
      subgfield(ind) = EXP(-r2/(s*s))
      
    END DO
  END DO
END DO



END SUBROUTINE putgausstosubgrid
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE addsubgrids(field1,field2,signum)
!------------------------------------------------------------
USE params
IMPLICIT NONE

REAL(DP), INTENT(IN OUT)                        :: field1(kdsub)
REAL(DP), INTENT(IN)                         :: field2(kdsub)
REAL(DP), INTENT(IN)             :: signum

INTEGER :: ii

DO ii=1,kdsub
  field1(ii) = field1(ii) + signum*field2(ii)
END DO

END SUBROUTINE addsubgrids
!------------------------------------------------------------
#endif

!------------------------------------------------------------

INTEGER FUNCTION conv3to1(k,j,i)
!------------------------------------------------------------
USE params
IMPLICIT NONE

INTEGER, INTENT(IN):: i, j, k

conv3to1=(i-1)*2*nx*2*ny+(j-1)*2*nx+k

RETURN
END FUNCTION conv3to1
!------------------------------------------------------------


!------------------------------------------------------------

INTEGER FUNCTION getnearestgridpoint(xpos, ypos,zpos)
!------------------------------------------------------------
USE params
IMPLICIT NONE
REAL(DP),INTENT(IN)              :: xpos 
REAL(DP),INTENT(IN)              :: ypos
REAL(DP),INTENT(IN)              :: zpos
!------------------------------------------------------------
!  getNearestGridPoint

INTEGER :: ind, ix0, iy0, iz0
REAL(DP) :: x0, y0, z0

!     THE NINT FUNCTION TURNED OUT TO BE EXPENSIVE!!
!      ix0 = nint(xpos/dx) + nx
!      iy0 = nint(ypos/dy) + ny
!      iz0 = nint(zpos/dz) + nz

x0=xpos/dx
ix0=NINT(x0) ! cPW
!xt=x0-ix0   ! cPW

y0=ypos/dy
iy0=NINT(y0) ! cPW
!yt=y0-iy0   ! cPW

z0=zpos/dz
iz0=NINT(z0) ! cPW
!zt=z0-iz0   ! cPW

!IF (xt > 00.5D0) ix0=ix0+1     ! cPW
!IF (yt > 00.5D0) iy0=iy0+1     ! cPW
!IF (zt > 00.5D0) iz0=iz0+1     ! cPW

ix0=ix0+nx
iy0=iy0+ny
iz0=iz0+nz

IF(ix0 > maxx) ix0=maxx ! avoid grid point out of box
IF(ix0 < minx) ix0=minx
IF(iy0 > maxy) iy0=maxy
IF(iy0 < miny) iy0=miny
IF(iz0 > maxz) iz0=maxz
IF(iz0 < minz) iz0=minz

! convert triple indices to single index:

ind =(iz0-1)*2*nx*2*ny+(iy0-1)*2*nx+ix0

getnearestgridpoint = ind

IF (ind > kdfull2 .OR. ind <= 0) STOP 'invalid grid point in getnearestgridpoint'

RETURN
END FUNCTION getnearestgridpoint

!------------------------------------------------------------
INTEGER FUNCTION getnearestgridpoint2(xpos, ypos,zpos)
!------------------------------------------------------------
USE params
IMPLICIT NONE

!------------------------------------------------------------
!  getNearestGridPoint

REAL(DP),INTENT(IN)              :: xpos 
REAL(DP),INTENT(IN)              :: ypos
REAL(DP),INTENT(IN)              :: zpos

INTEGER :: ind, ix0, iy0, iz0
REAL(DP) :: x0, y0, z0

!     THE NINT FUNCTION TURNED OUT TO BE EXPENSIVE!!
!      ix0 = nint(xpos/dx) + nx
!      iy0 = nint(ypos/dy) + ny
!      iz0 = nint(zpos/dz) + nz

x0=xpos/dx
ix0=INT(x0+nxsh)

y0=ypos/dy
iy0=INT(y0+nysh)

z0=zpos/dz
iz0=INT(z0+nzsh)


IF(ix0 > maxx) ix0=maxx ! avoid grid point out of box
IF(ix0 < minx) ix0=minx
IF(iy0 > maxy) iy0=maxy
IF(iy0 < miny) iy0=miny
IF(iz0 > maxz) iz0=maxz
IF(iz0 < minz) iz0=minz

WRITE(*,*) ' ix0,iy0,iz0=',ix0,iy0,iz0

! convert triple indices to single index:

ind =(iz0-1)*nx2*ny2+(iy0-1)*nx2+ix0

getnearestgridpoint2 = ind

IF (ind > kdfull2 .OR. ind <= 0) STOP 'invalid grid point in getnearestgridpoint2'

RETURN
END FUNCTION getnearestgridpoint2
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE conv1to3(ind)
!------------------------------------------------------------
USE params
IMPLICIT NONE
!     single index to triple index:

INTEGER,INTENT(IN) :: ind

INTEGER :: ind1, indx, indy, indz, iupper
REAL(DP) :: rupper

ind1 = ind

rupper = REAL(ind1,DP)/(2*nx*2*ny)

IF (rupper == INT(rupper)) THEN
  iupper = INT(rupper)
ELSE
  iupper = INT(rupper)+1
END IF


indz = iupper

ind1 = ind1 - (indz-1)*2*nx*2*ny

rupper = REAL(ind1,DP)/(2*nx)
IF (rupper == INT(rupper)) THEN
  iupper = INT(rupper)
ELSE
  iupper = INT(rupper)+1
END IF

indy =  iupper
indx = ind1 - (indy-1)*2*nx


iindtmp(1) = indx
iindtmp(2) = indy
iindtmp(3) = indz


!         conv1To3 = 1
!     return
END SUBROUTINE conv1to3

!----------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION getxval(ind)
!------------------------------------------------------------
USE params
IMPLICIT NONE
INTEGER, INTENT(IN) :: ind

CALL conv1to3(ind)

getxval = (iindtmp(1)-nx)*dx
RETURN
END FUNCTION getxval
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION getyval(ind)
!------------------------------------------------------------
USE params
IMPLICIT NONE
INTEGER, INTENT(IN) :: ind

CALL conv1to3(ind)

getyval = (iindtmp(2)-ny)*dy
RETURN
END FUNCTION getyval
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION getzval(ind)
!------------------------------------------------------------
USE params
IMPLICIT NONE
INTEGER, INTENT(IN) :: ind

CALL conv1to3(ind)

getzval = (iindtmp(3)-nz)*dz
RETURN
END FUNCTION getzval
!------------------------------------------------------------

#if(raregas)
!------------------------------------------------------------

SUBROUTINE updatearpos(iions,iels,ikats)
!------------------------------------------------------------
USE params
IMPLICIT NONE

! ONLY FOR THE ARGON AND TRADITIONAL SODIUM CASE!!



INTEGER, INTENT(IN)                      :: iions
INTEGER, INTENT(IN)                      :: iels
INTEGER, INTENT(IN)                      :: ikats

#if(raregas)

INTEGER :: i

IF (iions /= 0) THEN
  DO i=1,nc
    xc(i) = cx(i)
    yc(i) = cy(i)
    zc(i) = cz(i)
  END DO
END IF
IF (iels /= 0) THEN
  DO i=1,NE
    xe(i) = cx(i+nrare)
    ye(i) = cy(i+nrare)
    ze(i) = cz(i+nrare)
  END DO
END IF
IF (ikats /= 0) THEN
  DO i=1,nk
    xk(i) = cx(i+2*nrare)
    yk(i) = cy(i+2*nrare)
    zk(i) = cz(i+2*nrare)
  END DO
END IF
#endif

END SUBROUTINE updatearpos
!------------------------------------------------------------
#endif


!------------------------------------------------------------

SUBROUTINE renormsg(field,isgc,fromv,tov)
!------------------------------------------------------------
USE params
IMPLICIT NONE


REAL(DP), INTENT(IN OUT)                        :: field(kdfull2)
INTEGER, INTENT(IN)                      :: isgc
REAL(DP), INTENT(IN)                         :: fromv
REAL(DP), INTENT(IN)                         :: tov

INTEGER :: ind, ix, iy, iz

INTEGER,EXTERNAL :: conv3to1

ind = isgc

CALL conv1to3(ind)

DO iz=iindtmp(3)-nzsg,iindtmp(3)+nzsg
  DO iy=iindtmp(2)-nysg,iindtmp(2)+nysg
    DO ix=iindtmp(1)-nxsg,iindtmp(1)+nxsg
      
      ind = conv3to1(ix,iy,iz)
      
      field(ind) = field(ind)*tov/fromv
      
    END DO
  END DO
END DO


END SUBROUTINE renormsg
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION foldfunc(field,func,xx,yy,zz,param)
!------------------------------------------------------------
USE params
IMPLICIT NONE

!     folds function around xx,yy,zz

!     calculates the integral

!     \int d^3r field(\vec r) \times func(|\vec r-\vec R|)

!     where \vec R = (xx,yy,zz)

REAL(DP), INTENT(IN)                         :: field(kdfull2)
REAL(DP), INTENT(IN)                         :: xx
REAL(DP), INTENT(IN)                         :: yy
REAL(DP), INTENT(IN)                         :: zz
REAL(DP), INTENT(IN)                     :: param
! Function name passed as an argument :
INTERFACE 
  REAL(DP) FUNCTION func(r,s)
    USE params
    REAL(DP), INTENT(IN) :: r
    REAL(DP), INTENT(IN) :: s
  END FUNCTION func
END INTERFACE

INTEGER :: ind, ix, iy, iz
REAL(DP) :: acc, rr, rx, ry, rz, x1, y1, z1


ind = 0
acc = 0D0

DO iz=minz,maxz
  z1=(iz-nzsh)*dz
  rz=z1-zz
  DO iy=miny,maxy
    y1=(iy-nysh)*dy
    ry=y1-yy
    DO ix=minx,maxx
      x1=(ix-nxsh)*dx
      rx=x1-xx
      rr=SQRT(rx*rx+ry*ry+rz*rz)
      rr=MAX(rr,small) ! avoid zero
      ind=ind+1
      
      acc = acc + field(ind)*func(rr,param)
      
    END DO
  END DO
END DO

foldfunc = acc * dvol

RETURN
END FUNCTION foldfunc
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE foldgradfunc(field,func,xx,yy,zz,param)
!------------------------------------------------------------
USE params
IMPLICIT NONE

!     folds function around xx,yy,zz

!     calculates the integral

!     \int d^3r field(\vec r) \times \nabla_R func(|\vec r-\vec R|)

!     where \vec R = (xx,yy,zz)



REAL(DP), INTENT(IN)                         :: field(kdfull2)
REAL(DP), INTENT(IN)                         :: xx
REAL(DP), INTENT(IN)                         :: yy
REAL(DP), INTENT(IN)                         :: zz
REAL(DP), INTENT(IN)         :: param
! Function name passed as an argument :
INTERFACE 
  REAL(DP) FUNCTION func(r,s)
    USE params
    REAL(DP), INTENT(IN) :: r
    REAL(DP), INTENT(IN) :: s
  END FUNCTION func
END INTERFACE

INTEGER :: ind, ix, iy, iz
REAL(DP) :: rder, sumx, sumy, sumz
REAL(DP):: radial, rr, rx, ry, rz, x1, y1, z1

ind = 0
sumx = 0D0
sumy = 0D0
sumz = 0D0
rder = 1.0D-5


DO iz=minz,maxz
  z1=(iz-nzsh)*dz
  rz=zz-z1
  DO iy=miny,maxy
    y1=(iy-nysh)*dy
    ry=yy-y1
    DO ix=minx,maxx
      x1=(ix-nxsh)*dx
      rx=xx-x1
      rr=SQRT(rx*rx+ry*ry+rz*rz)
      rr=MAX(rr,small) ! avoid zero
      ind=ind+1
      
      radial = func(rr+rder,param)-func(rr-rder,param)
      
      
      sumx = sumx + field(ind)*radial*rx/rr
      sumy = sumy + field(ind)*radial*ry/rr
      sumz = sumz + field(ind)*radial*rz/rr
      
    END DO
  END DO
END DO



rvectmp(1) = sumx * dvol/(rder+rder)
rvectmp(2) = sumy * dvol/(rder+rder)
rvectmp(3) = sumz * dvol/(rder+rder)


RETURN
END SUBROUTINE foldgradfunc
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE foldgradfunconsubgrid(field,func,xx,yy,zz,param)
!------------------------------------------------------------
USE params
IMPLICIT NONE

!     folds function around xx,yy,zz using a subgrid only

!     calculates the integral

!     \int d^3r field(\vec r) \times \nable_R func(|\vec r-\vec R|)

!     where \vec R = (xx,yy,zz)


REAL(DP), INTENT(IN)                         :: field(kdfull2)
REAL(DP), INTENT(IN)                         :: xx
REAL(DP), INTENT(IN)                         :: yy
REAL(DP), INTENT(IN)                         :: zz
REAL(DP), INTENT(IN)         :: param
! Function name passed as an argument :
INTERFACE 
  REAL(DP) FUNCTION func(r,s)
    USE params
    REAL(DP), INTENT(IN) :: r
    REAL(DP), INTENT(IN) :: s
  END FUNCTION func
END INTERFACE

INTEGER,EXTERNAL :: getnearestgridpoint
INTEGER,EXTERNAL :: conv3to1

INTEGER :: ind, ix, iy, iz
REAL(DP) :: radial, rr, rder, rx, ry, rz, sumx, sumy, sumz, x1, y1, z1

ind = getnearestgridpoint(xx,yy,zz)

CALL conv1to3(ind)


sumx = 0.0D0
sumy = 0.0D0
sumz = 0.0D0
rder = 1.0D-5


DO iz=iindtmp(3)-nzsg,iindtmp(3)+nzsg
  IF(iz >= minz.AND.iz <= maxz)THEN
    z1=(iz-nzsh)*dz
    rz=zz-z1
    
    DO iy=iindtmp(2)-nysg,iindtmp(2)+nysg
      IF(iy >= miny.AND.iy <= maxy)THEN
        y1=(iy-nysh)*dy
        ry=yy-y1
        
        DO ix=iindtmp(1)-nxsg,iindtmp(1)+nxsg
          IF(ix >= minx.AND.ix <= maxx)THEN
            x1=(ix-nxsh)*dx
            rx=xx-x1
            
            rr=SQRT(rx*rx+ry*ry+rz*rz)
            rr=MAX(rr,small)
            
            ind = conv3to1(ix,iy,iz)
            
            radial = func(rr+rder,param)-func(rr-rder,param)
            
            sumx = sumx + field(ind)*radial*rx/rr
            sumy = sumy + field(ind)*radial*ry/rr
            sumz = sumz + field(ind)*radial*rz/rr
            
          END IF
        END DO
      END IF
    END DO
  END IF
END DO


rvectmp(1) = sumx * dvol/(rder+rder)
rvectmp(2) = sumy * dvol/(rder+rder)
rvectmp(3) = sumz * dvol/(rder+rder)


RETURN
END SUBROUTINE foldgradfunconsubgrid
!------------------------------------------------------------



!------------------------------------------------------------

REAL(DP) FUNCTION dintprod(field1,field2)
!------------------------------------------------------------
USE params
IMPLICIT NONE


REAL(DP), INTENT(IN)                         :: field1(kdfull2)
REAL(DP), INTENT(IN)                         :: field2(kdfull2)


INTEGER :: ind, ix, iy, iz
REAL(DP):: acc

ind = 0
acc = 0.0D0
DO iz=minz,maxz
  DO iy=miny,maxy
    DO ix=minx,maxx
      ind=ind+1
      acc = acc + field1(ind)*field2(ind)
    END DO
  END DO
END DO

dintprod = acc*dvol


RETURN
END FUNCTION dintprod
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION dintprodsubgrid(field1,field2,xx,yy,zz)
!------------------------------------------------------------
USE params
IMPLICIT NONE

!     calculates the integral of the product of two scalar
!     fields not over the whole box, but only over a subgrid
!     which is centered at xx,yy,zz


REAL(DP), INTENT(IN)                         :: field1(kdfull2)
REAL(DP), INTENT(IN)                         :: field2(kdfull2)
REAL(DP), INTENT(IN)                         :: xx
REAL(DP), INTENT(IN)                         :: yy
REAL(DP), INTENT(IN)                         :: zz

INTEGER :: ind, ix, iy, iz
REAL(DP) :: acc

INTEGER,EXTERNAL :: getnearestgridpoint
INTEGER,EXTERNAL :: conv3to1



ind = getnearestgridpoint(xx,yy,zz)

CALL conv1to3(ind)

acc = 0.0D0
DO iz=iindtmp(3)-nzsg,iindtmp(3)+nzsg
  DO iy=iindtmp(2)-nysg,iindtmp(2)+nysg
    DO ix=iindtmp(1)-nxsg,iindtmp(1)+nxsg
      
      ind = conv3to1(ix,iy,iz)
      
      acc = acc + field1(ind)*field2(ind)
      
    END DO
  END DO
END DO

dintprodsubgrid = acc * dvol

RETURN
END FUNCTION dintprodsubgrid
!------------------------------------------------------------



!------------------------------------------------------------

REAL(DP) FUNCTION dintfieldtimesgauss(field,xx,yy,zz,sigm,ccharg)
!------------------------------------------------------------
USE params
IMPLICIT NONE

!     calculates the integral:

!             \int d^3r f(\vec{r})\cdot A exp(-|\vec{r}-\vec{R}|/(2\sigma^2))
!     with A being the proper factor normalising the Gaussian to ccharg

!     ATTENTION: be aware that the factor e2 is NOT included here,
!                but must be attached in the calling routine



!     xx,yy,zz are the center of the Gaussian, sigm its width
!     according to the convention in this code:
!     exp(-r**2/(2.*sigm**2)) !!!



REAL(DP), INTENT(IN)                         :: field(kdfull2)
REAL(DP), INTENT(IN)                         :: xx
REAL(DP), INTENT(IN)                         :: yy
REAL(DP), INTENT(IN)                         :: zz
REAL(DP), INTENT(IN)                         :: sigm
REAL(DP), INTENT(IN)                         :: ccharg

INTEGER :: ind, ix, iy, iz
REAL(DP) :: a, acc, sigm2 
REAL(DP) :: r2, rx, ry, rz, x1, y1, z1

INTEGER,EXTERNAL :: getnearestgridpoint
INTEGER,EXTERNAL :: conv3to1



a = ccharg/(pi**1.5D0*2D0**1.5D0*sigm**3)
sigm2 = sigm*sigm
acc = 0D0

IF (ipseudo == 0) THEN
  
  ind = 0
  
  DO iz=minz,maxz
    z1=(iz-nzsh)*dz
    rz=z1-zz
    DO iy=miny,maxy
      y1=(iy-nysh)*dy
      ry=y1-yy
      DO ix=minx,maxx
        x1=(ix-nxsh)*dx
        rx=x1-xx
        
        r2 = rx*rx+ry*ry+rz*rz
        
        acc = acc + field(ind)*EXP(-r2/2D0/sigm2)
        
      END DO
    END DO
  END DO
  
ELSE ! ipseudo = 1
  
  ind = getnearestgridpoint(xx,yy,zz)
  
  CALL conv1to3(ind)
  
  
  
  
  DO iz=iindtmp(3)-nzsg,iindtmp(3)+nzsg
    z1=(iz-nzsh)*dz
    rz=z1-zz
    
    DO iy=iindtmp(2)-nysg,iindtmp(2)+nysg
      y1=(iy-nysh)*dy
      ry=y1-yy
      
      DO ix=iindtmp(1)-nxsg,iindtmp(1)+nxsg
        x1=(ix-nxsh)*dx
        rx=x1-xx
        
        r2=rx*rx+ry*ry+rz*rz
        
        
        ind = conv3to1(ix,iy,iz)
        
        acc = acc + field(ind)*EXP(-r2/2D0/sigm2)
        
      END DO
    END DO
  END DO
  
END IF


acc = acc * a * dvol

dintfieldtimesgauss = acc

RETURN
END FUNCTION dintfieldtimesgauss
!------------------------------------------------------------

#if(raregas)
!------------------------------------------------------------

INTEGER FUNCTION iptyp(ind)
!------------------------------------------------------------
USE params
IMPLICIT NONE
!     yields the particle type of particle with "long" index ind


!     result means:

!                  1 ---> anion core
!                  2 ---> anion valence shell
!                  3 ---> cation
INTEGER, INTENT(IN) :: ind


IF (ind <= nc) THEN
  iptyp = 1
ELSE IF (ind <= nc+NE) THEN
  iptyp = 2
ELSE
  iptyp = 3
END IF

RETURN
END FUNCTION iptyp
!------------------------------------------------------------

!------------------------------------------------------------

INTEGER FUNCTION ismobile(ind)
!------------------------------------------------------------
USE params
IMPLICIT NONE
INTEGER, INTENT(IN) :: ind

INTEGER :: ireturn
ireturn = 1

IF (ind <= nc) THEN
  IF (imobc(ind) /= 1) ireturn = 0
ELSE IF (ind <= nc+NE) THEN
  IF (imobe(ind-nc) /= 1) ireturn = 0
ELSE
  IF (imobk(ind-nc-NE) /= 1) ireturn = 0
END IF

ismobile = ireturn

RETURN
END FUNCTION ismobile
!------------------------------------------------------------


!------------------------------------------------------------

SUBROUTINE getpositionold(ii)
!------------------------------------------------------------
USE params
IMPLICIT NONE
INTEGER, INTENT(IN) :: ii

IF (ii <= nc) THEN
  rvectmp(1)=xcold(ii)
  rvectmp(2)=ycold(ii)
  rvectmp(3)=zcold(ii)
ELSE IF (ii <= nc+NE) THEN
  rvectmp(1)=xeold(ii)
  rvectmp(2)=yeold(ii)
  rvectmp(3)=zeold(ii)
ELSE
  rvectmp(1)=xkold(ii)
  rvectmp(2)=ykold(ii)
  rvectmp(3)=zkold(ii)
END IF

RETURN
END SUBROUTINE getpositionold
!------------------------------------------------------------
#endif

!------------------------------------------------------------

INTEGER FUNCTION isoutofbox(x,y,z)
!------------------------------------------------------------
USE params
IMPLICIT NONE

!     return values:
!        0 --> particle is inside of the box
!        1 --> particle is inside of the box but close to the
!              border, so that it cannot be described by
!              pseudodensity formalism
!        2 --> particle is far out of the box

REAL(DP),INTENT(IN) :: x,y,z

IF (z < (1-nzsh)*dz+dinmargin .OR. z > (nzsh*dz-dinmargin)  &
    .OR. y < (1-nysh)*dy+dinmargin .OR. y > (nysh*dy-dinmargin)  &
    .OR.(x < (1-nxsh)*dx+dinmargin .OR. x > (nxsh*dx-dinmargin))) THEN
  isoutofbox = 1   ! in the box but close to the border
ELSE
  isoutofbox = 0   ! in the box
END IF

! in the box and not close to the border  
IF (z < (1-nzsh)*dz .OR. z > (nzsh*dz) &
    .OR. y < (1-nysh)*dy .OR. y > (nysh*dy)  &
    .OR.(x < (1-nxsh)*dx .OR. x > (nxsh*dx))) THEN
  isoutofbox = 2
END IF





RETURN
END FUNCTION isoutofbox
!------------------------------------------------------------
