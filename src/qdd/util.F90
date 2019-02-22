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



MODULE util

IMPLICIT NONE

INTERFACE projmoms
  MODULE PROCEDURE r_projmoms, c_projmoms
END INTERFACE projmoms

INTERFACE wfovlp
  MODULE PROCEDURE r_wfovlp, c_wfovlp
END INTERFACE wfovlp

INTERFACE wfnorm
  MODULE PROCEDURE r_wfnorm, c_wfnorm
END INTERFACE wfnorm

INTERFACE project
  MODULE PROCEDURE r_project, c_project
END INTERFACE project

CONTAINS

!------------------------------------------------------------
SUBROUTINE shiftfield(q0,shix,shiy,shiz)
!------------------------------------------------------------
USE params
IMPLICIT NONE

!  Shifts array 'q0' by shix,shiy,shiz along grid points. Shift values
!  shix,shiy,shiz must be positive
!  
!  Input/Output:
!    q0   = real array for one spatial wavefunction to be shifted
!  Input:
!    shix = shift distance along x direction, nearest multiple of 
!           grid spacing is actually taken; must be > 0.
!    shiy = as shix, but for y-direction
!    shiz = as shix, but for y-direction

REAL(DP), INTENT(IN OUT) :: q0(kdfull2)
REAL(DP), INTENT(IN)     :: shix
REAL(DP), INTENT(IN)     :: shiy
REAL(DP), INTENT(IN)     :: shiz

REAL(DP),ALLOCATABLE ::  q1(:)

INTEGER :: ind, ind1, ix, iy, iz , ix1, iy1, iz1, nshift
INTEGER, EXTERNAL :: conv3to1

ALLOCATE(q1(kdfull2))

q1=q0

nshift = nint(shix/dx)
if(nshift < 0) STOP ' SHIX must be >= 0'

IF(nshift > 0) then
  DO ix=minx+nshift,maxx
    DO iy=miny,maxy
      DO iz=minz,maxz
      
        ix1 = ix - nshift
      
        ind = conv3to1(ix,iy,iz)
        ind1 = conv3to1(ix1,iy,iz)
      
        q0(ind) = q1(ind1)
      
      END DO
    END DO
  END DO

  DO ix=minx,minx+nshift-1
    DO iy=miny,maxy
      DO iz=minz,maxz
      
        ind = conv3to1(ix,iy,iz)
      
        q0(ind)=0D0
      
      END DO
    END DO
  END DO
END IF

!::::::::::::::::::::::::::::::


q1=q0

nshift = nint(shiy/dy)
if(nshift < 0) STOP ' SHIY must be >= 0'

IF(nshift > 0) then
  DO iy=miny+nshift,maxy
    DO ix=minx,maxx
      DO iz=minz,maxz
      
        iy1 = iy - nshift
      
        ind = conv3to1(ix,iy,iz)
        ind1 = conv3to1(ix,iy1,iz)
      
        q0(ind) = q1(ind1)
      
      END DO
    END DO
  END DO

  DO iy=miny,miny+nshift-1
    DO ix=minx,maxx
      DO iz=minz,maxz
      
        ind = conv3to1(ix,iy,iz)
      
        q0(ind)=0D0
      
      END DO
    END DO
  END DO
END IF


!::::::::::::::::::::::::::::::

q1=q0

nshift = nint(shiz/dz)
if(nshift < 0) STOP ' SHIZ must be >= 0'

IF(nshift > 0) then
  DO iz=minz+nshift,maxz
    DO ix=minx,maxx
      DO iy=miny,maxy
      
        iz1 = iz - nshift
      
        ind = conv3to1(ix,iy,iz)
        ind1 = conv3to1(ix,iy,iz1)
      
        q0(ind) = q1(ind1)
      
      END DO
    END DO
  END DO

  DO iz=minz,minz+nshift-1
    DO ix=minx,maxx
      DO iy=miny,maxy
      
        ind = conv3to1(ix,iy,iz)
      
        q0(ind)=0D0
      
      END DO
    END DO
  END DO
END IF

DEALLOCATE(q1)


RETURN
END SUBROUTINE shiftfield
!------------------------------------------------------------

!------------------------------------------------------------
CHARACTER*9 FUNCTION inttostring(inumber)
!------------------------------------------------------------

! translates integer number to corresponding character variable

IMPLICIT NONE
INTEGER, INTENT(IN)                      :: inumber

SELECT CASE(inumber)
  CASE(0:999999999)
    WRITE(inttostring,'(i9)') inumber
  CASE DEFAULT
    STOP 'ERROR in intToString : input argument outside range 0 to 999999999'
END SELECT

RETURN
END FUNCTION inttostring
!------------------------------------------------------------

!------------------------------------------------------------
SUBROUTINE pm3dcut(iunit,iax1,iax2,value,field)
!------------------------------------------------------------

! Prepares 3D plotting with 'pm3d' in 'gnuplot.
!
! Input:
!   iunit      = write unit number-
!   iax1,iax2  = choices of the two axes (from 1..3) along which is plotted.
!   value      = coordinate on third axis (complement of iax1,iax2) along
!                which field is plotted.
!   field      = 3D array from which 2D cut is retrieved.

USE params
IMPLICIT NONE


INTEGER, INTENT(IN)   :: iunit
INTEGER, INTENT(IN)   :: iax1
INTEGER, INTENT(IN)   :: iax2
REAL(DP), INTENT(IN)  :: value
REAL(DP), INTENT(IN)  :: field(kdfull2)

INTEGER :: ind, ix, iy, iz, x1, y1, z1

IF (iax1 == 1 .AND. iax2 == 2) THEN ! xy-plane
  
  ind=0
  DO iz=1,maxz
    z1=(iz-nzsh)*dz
    
    DO iy=1,maxy
      y1=(iy-nysh)*dy
      DO ix=1,maxx
        x1=(ix-nxsh)*dx
        
        ind=ind+1
        
        IF (z1 == value) THEN
          WRITE(iunit,'(3f10.3,1e15.6)') x1,y1,z1,field(ind)
!               write(6,'(3f10.3,1e15.6)') x1,y1,z1,field(ind)
        END IF
      END DO
      IF (z1 == value) THEN
        WRITE(iunit,*)
!                  write(6,*)
      END IF
    END DO
   
  END DO
  
  
ELSE IF (iax1 == 1 .AND. iax2 == 3) THEN ! xz-plane
  
  ind=0
  DO iz=1,maxz
    z1=(iz-nzsh)*dz
    
    DO iy=1,maxy
      y1=(iy-nysh)*dy
      
      DO ix=1,maxx
        x1=(ix-nxsh)*dx
        
        ind=ind+1
        IF (y1 == value) THEN
          WRITE(iunit,'(3f10.3,1e15.6)') x1,y1,z1,field(ind)
        END IF
        
      END DO
      IF (y1 == value) THEN
        WRITE(iunit,*)
      END IF
      
    END DO
    
  END DO
  
  
  
ELSE IF (iax1 == 2 .AND. iax2 == 3) THEN ! yz-plane
  
  ind=0
  DO iz=1,maxz
    z1=(iz-nzsh)*dz
    
    DO iy=1,maxy
      y1=(iy-nysh)*dy
      
      
      DO ix=1,maxx
        x1=(ix-nxsh)*dx
        
        ind=ind+1
        IF (x1 == value) THEN
          WRITE(iunit,'(3f10.3,1e15.6)') x1,y1,z1,field(ind)
        END IF
        
      END DO
      
    END DO
    WRITE(iunit,*)
  END DO
  
END IF



RETURN
END SUBROUTINE pm3dcut
!------------------------------------------------------------

!------------------------------------------------------------
COMPLEX(DP) FUNCTION orbitaloverlap(q1,q2)
!------------------------------------------------------------
USE params

IMPLICIT NONE

!  Overlap <q1|q2> = \int\d^3r q_1^* q_2 between the two complex
!  wavefunctions q1 and q2.


COMPLEX(DP), INTENT(IN) :: q1(kdfull2)
COMPLEX(DP), INTENT(IN) :: q2(kdfull2)

INTEGER :: ind
REAL(DP) :: sumr, sumi
COMPLEX(DP) :: csum
!WRITE(6,*) 'Entering orbitalOverlap'

sumr=0D0
sumi=0D0

DO ind=1,kdfull2
  sumr=sumr+REAL(q1(ind),DP)*REAL(q2(ind),DP) + AIMAG(q1(ind))*AIMAG(q2(ind))
  sumi=sumi+REAL(q1(ind),DP)*AIMAG(q2(ind))- REAL(q2(ind),DP)*AIMAG(q1(ind))
END DO


csum=CMPLX(sumr*dx*dy*dz,sumi*dx*dy*dz,DP)

orbitaloverlap=csum

RETURN
END FUNCTION orbitaloverlap
!------------------------------------------------------------

!------------------------------------------------------------
REAL(DP) FUNCTION realoverlap(q1,q2)
!------------------------------------------------------------
USE params
IMPLICIT NONE

!     returns \int\d^3r q_1^* q_2

COMPLEX(DP), INTENT(IN)                  :: q1(kdfull2)
COMPLEX(DP), INTENT(IN)                  :: q2(kdfull2)

INTEGER :: ind
REAL(DP) :: sumr

sumr=0D0

DO ind=1,kdfull2
  sumr=sumr+REAL(q1(ind),DP)*REAL(q2(ind),DP) + AIMAG(q1(ind))*AIMAG(q2(ind))
END DO

realoverlap=sumr*dvol

RETURN
END FUNCTION realoverlap
!------------------------------------------------------------

!-------------------------------------------------------------
REAL(DP) FUNCTION realovsubgrid(q1,q2,ion)

!     Real part of overlap between wavefunctions 'q1' and 'q2'
!     accumulated on the subgrid of ion 'ion' defined by
!     the non-local PsP.

USE params
IMPLICIT NONE

COMPLEX(DP), INTENT(IN)                  :: q1(kdfull2)
COMPLEX(DP), INTENT(IN)                  :: q2(kdfull2)
INTEGER, INTENT(IN)                  :: ion

INTEGER :: i, ii
REAL(DP) :: sumr
!-------------------------------------------------------------


sumr=0D0
DO i=1,ifin(ion)
  ii = icount(i,ion)
  sumr = REAL(q2(ii),DP)*REAL(q1(ii),DP) +AIMAG(q2(ii))*AIMAG(q1(ii)) + sumr
END DO

realovsubgrid=sumr*dvol

RETURN
END FUNCTION realovsubgrid
!------------------------------------------------------------

!------------------------------------------------------------
SUBROUTINE gettemperature(iflag)
!------------------------------------------------------------

! Compute ionic kinetic energy and print corresponding ionic temperature 
! on unit 175.

USE params
IMPLICIT NONE

INTEGER,INTENT(IN) :: iflag

INTEGER ::  i
REAL(DP) :: ek, ekcm, pcmx, pcmy, pcmz, rkt, rm, rmcm, tt
#if(raregas)
REAL(DP) :: rmc, rmk
#endif

IF (iflag == 4) THEN
  
  pcmx = 0D0
  pcmy = 0D0
  pcmz = 0D0
  
  DO i=1,nion
    pcmx=pcmx+cpx(i)
    pcmy=pcmy+cpy(i)
    pcmz=pcmz+cpz(i)
  END DO
  
  rmcm = nion * amu(np(1))*1836D0*0.5D0
  rm = amu(np(1))*1836D0*0.5D0
  ekcm= pcmx**2+pcmy**2+pcmz**2
  ekcm = ekcm/2D0/rmcm
  
  ek=0D0
  DO i=1,nion
    ek = ek + (cpx(i)**2+cpy(i)**2+cpz(i)**2)/2D0/rm
  END DO
  
  rkt = 2D0/(3D0*nion-3D0)*(ek - ekcm) ! thermal energy kT in Ry
  
! Temperature in Kelvin
  
  tt = rkt/(6.507D-6)
  
  WRITE(151,'(2f12.4)') tfs,tt

#if(raregas)
ELSE IF (iflag <= 3) THEN
  
  pcmx = 0D0
  pcmy = 0D0
  pcmz = 0D0
  
  DO i=1,nc
    pcmx=pcmx+pxc(i)
    pcmy=pcmy+pyc(i)
    pcmz=pcmz+pzc(i)
  END DO
  DO i=1,nk
    pcmx=pcmx+pxk(i)
    pcmy=pcmy+pyk(i)
    pcmz=pcmz+pzk(i)
  END DO
  
  
  
  rmcm = (nc * mion+nk*mkat)*1836*ame
  rmc = mion*1836*ame
  rmk = mkat*1836*ame
  
  ekcm= pcmx**2+pcmy**2+pcmz**2
  ekcm = ekcm/2D0/rmcm
  
  ek=0D0
  DO i=1,nc
    ek = ek + (pxc(i)**2+pyc(i)**2+pzc(i)**2)/2./rmc
  END DO
  DO i=1,nk
    ek = ek + (pxk(i)**2+pyk(i)**2+pzk(i)**2)/2./rmk
  END DO
  
  rkt = 2D0/(3D0*(nc+nk)-3D0)*(ek - ekcm) ! thermal energy kT in Ry
  
! Temperature in Kelvin
  
  tt = rkt/(6.507D-6)
  
  WRITE(157,'(2f12.4)') tfs,tt
#endif  
  
END IF

RETURN
END SUBROUTINE gettemperature
!------------------------------------------------------------

SUBROUTINE densbelow(field)
!------------------------------------------------------------

! Integrate 'field' ober x-y plane, accumulate along z axis and print
! on unit 650.

USE params
IMPLICIT NONE

REAL(DP), INTENT(IN) :: field(kdfull2)

INTEGER :: ind, ix, iy, iz
REAL(DP) :: acc, x1, y1, z1
REAL(DP),DIMENSION(:),ALLOCATABLE :: vs(:)

ALLOCATE(vs(maxz+1))
ind = 0

acc = 0D0

DO iz=minz,maxz
  z1=(iz-nzsh)*dz
  DO iy=miny,maxy
    y1=(iy-nysh)*dy
    DO ix=minx,maxx
      x1=(ix-nxsh)*dx
      ind = ind + 1
      acc = acc + field(ind)
    END DO
  END DO
  vs(iz) = acc * dx*dy*dz
  WRITE(650,'(2f10.5,1e20.9)') tfs, z1, vs(iz)
END DO

DEALLOCATE(vs)

RETURN
END SUBROUTINE densbelow
!------------------------------------------------------------

!-----priCM-------------------------------------------------------
SUBROUTINE pricm(rho)
!------------------------------------------------------------
USE params
IMPLICIT NONE

!     calculates center of density rho
!     returns result on rVecTmp(1:3) via module params

REAL(DP), INTENT(IN)  :: rho(kdfull2)

INTEGER :: ind, ix, iy, iz
REAL(DP) :: acc, xcm, ycm, zcm, x1, y1, z1

!------------------------------------------------------------

xcm = 0D0
ycm = 0D0
zcm = 0D0
acc = 0D0
!      xxcm = 0D0
!      yycm = 0D0
!      zzcm = 0D0


ind=0
DO iz=minz,maxz
  z1=(iz-nzsh)*dz
  DO iy=miny,maxy
    y1=(iy-nysh)*dy
    DO ix=minx,maxx
      x1=(ix-nxsh)*dx
!old               rr=sqrt(rx*rx+ry*ry+rz*rz)
!old               rr=max(rr,SMALL) ! avoid zero
      ind=ind+1
      
      xcm = xcm + rho(ind)*x1
      ycm = ycm + rho(ind)*y1
      zcm = zcm + rho(ind)*z1
!               xxcm = xxcm + rho(ind)*x1*rho(ind)*x1
!               yycm = yycm + rho(ind)*y1*rho(ind)*y1
!               zzcm = zzcm + rho(ind)*z1*rho(ind)*z1
      acc = acc + rho(ind)
    END DO
  END DO
END DO

rvectmp(1) = xcm/acc
rvectmp(2) = ycm/acc
rvectmp(3) = zcm/acc

END SUBROUTINE pricm
!------------------------------------------------------------

!-----priCM_state-------------------------------------------------------

SUBROUTINE pricm_state(psir)
!------------------------------------------------------------
USE params
IMPLICIT NONE

!     calculates center of density for states in wavefunctions 'psir'
!     and prints results on standard output

REAL(DP), INTENT(IN)  :: psir(kdfull2,kstate)

INTEGER :: i1, i2, i3, ind, nbr
REAL(DP) :: xx, yy, zz

!------------------------------------------------------------

DO nbr=1,nstate
  xx = 0D0
  yy = 0D0
  zz = 0D0
  ind = 0
  DO i3=minz,maxz
    DO i2=miny,maxy
      DO i1=minx,maxx
        ind = ind+1
        xx = xx + psir(ind,nbr)*psir(ind,nbr)*(i1-nxsh)*dx
        yy = yy + psir(ind,nbr)*psir(ind,nbr)*(i2-nysh)*dy
        zz = zz + psir(ind,nbr)*psir(ind,nbr)*(i3-nzsh)*dz
      END DO
    END DO
  END DO
  xx = xx*dvol
  yy = yy*dvol
  zz = zz*dvol
  WRITE(6,'(a,i5,a,i5,a,3(1pg13.5))') 'node: myn=',  &
      myn,', state: n=',nbr,'<x>,<y>,<z>=',xx,yy,zz
END DO

RETURN
END SUBROUTINE pricm_state


!------------------------------------------------------------

SUBROUTINE rotatevec(x,y,z,iax,alpha)
!------------------------------------------------------------
! rotation of vector '(x,y,z)' about axis 'iax' with angle 'alpha'
! returning the result on vector 'rvectmp' via module 'params'..
USE params
IMPLICIT NONE

REAL(DP), INTENT(IN)::x
REAL(DP), INTENT(IN)::y
REAL(DP), INTENT(IN)::z
INTEGER,INTENT(IN)::iax
REAL(DP),INTENT(in)::alpha


IF (iax == 1) THEN ! rotate along x-axis
  rvectmp(1)=x
  rvectmp(2)=COS(alpha)*y-SIN(alpha)*z     
  rvectmp(3)=SIN(alpha)*y+COS(alpha)*z     
ELSE IF (iax == 2) THEN
  rvectmp(1)=COS(alpha)*x+SIN(alpha)*z
  rvectmp(2)=y
  rvectmp(3)=-SIN(alpha)*x+COS(alpha)*z
ELSE IF (iax == 3) THEN
  rvectmp(1)=COS(alpha)*x-SIN(alpha)*y     
  rvectmp(2)=SIN(alpha)*x+COS(alpha)*y     
  rvectmp(3)=z
ELSE
  STOP 'Error in subroutine rotateVec'
END IF

RETURN
END SUBROUTINE rotatevec
!------------------------------------------------------------


!------------------------------------------------------------

SUBROUTINE rotatevec3D(vecin,vecout,vecalpha)
!------------------------------------------------------------
! rotation of 3D vector 'vecin' by rotation vector 'vecalpha'
! returning result on 'vecout'.
! Angle to be given in radian.
!
USE params, ONLY: DP
IMPLICIT NONE

REAL(DP), INTENT(IN) :: vecin(3)
REAL(DP), INTENT(IN) :: vecalpha(3)
REAL(DP), INTENT(OUT) :: vecout(3)

REAL(DP) :: absalpha,sa,c,c2a2

! prepare auxiliary variables

absalpha = SQRT(SUM(vecalpha**2))
c = COS(absalpha)
sa = SIN(absalpha)/absalpha
c2a2 = (c-1D0)/(absalpha*absalpha)

! compute transformation

vecout(1) = (c-c2a2*vecalpha(1)**2)*vecin(1) &
           -(sa*vecalpha(3)+c2a2*vecalpha(1)*vecalpha(2))*vecin(2) &
           +(sa*vecalpha(2)-c2a2*vecalpha(1)*vecalpha(3))*vecin(3) 

vecout(2) = (sa*vecalpha(3)-c2a2*vecalpha(1)*vecalpha(2))*vecin(1) &
           +(c-c2a2*vecalpha(2)**2)*vecin(2) &
           -(sa*vecalpha(1)+c2a2*vecalpha(2)*vecalpha(3))*vecin(3) 
           

vecout(3) =-(sa*vecalpha(2)+c2a2*vecalpha(1)*vecalpha(3))*vecin(1) &
           +(sa*vecalpha(1)-c2a2*vecalpha(2)*vecalpha(3))*vecin(2) &
           +(c-c2a2*vecalpha(3)**2)*vecin(3) 
           


RETURN
END SUBROUTINE rotatevec3D
!------------------------------------------------------------



!------------------------------------------------------------

SUBROUTINE swaparr(arr1,arr2,n)
!------------------------------------------------------------
USE params
IMPLICIT NONE

!     swaps arrays arr1 and arr2



REAL(DP), INTENT(IN OUT)                     :: arr1(kdfull2)
REAL(DP), INTENT(IN OUT)                     :: arr2(kdfull2)
INTEGER, INTENT(IN)                      :: n

INTEGER :: i
REAL(DP),ALLOCATABLE ::  arrt(:)

ALLOCATE(arrt(kdfull2))
DO i=1,n
  arrt(i)=arr1(i)
  arr1(i)=arr2(i)
  arr2(i)=arrt(i)
END DO
DEALLOCATE(arrt)

RETURN
END SUBROUTINE swaparr
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE cparr(arr1,arr2,n)
!------------------------------------------------------------
USE params
IMPLICIT NONE

!     copies array arr1 to arr2

REAL(DP), INTENT(IN)                         :: arr1(kdfull2)
REAL(DP), INTENT(OUT)                        :: arr2(kdfull2)
INTEGER, INTENT(IN)                      :: n

INTEGER :: i

DO i=1,n
  arr2(i)=arr1(i)
END DO

RETURN
END SUBROUTINE cparr
!------------------------------------------------------------



!------------------------------------------------------------

SUBROUTINE makexbsfile(ionlymob,iwithsh)
!------------------------------------------------------------

! Prepares input for plotting ionic configuration with the chemical
! grafical freeware 'xbs'.

USE params
IMPLICIT NONE

INTEGER,INTENT(IN) :: ionlymob
INTEGER,INTENT(IN) :: iwithsh

#if(raregas)
INTEGER :: i
REAL(DP) :: xx, yy, zz
#endif

OPEN(789,STATUS='unknown',FILE='xbs.Na_MgO')

#if(raregas)
DO i=1,nc
  IF (ionlymob == 0 .OR. imobc(i) == 1) THEN
    WRITE(789,'(a,3f17.7)') 'atom O  ',xc(i),yc(i),zc(i)
  END IF
END DO

IF (iwithsh /= 0) THEN
  DO i=1,NE
    IF (ionlymob == 0 .OR. imobe(i) == 1) THEN
      xx = xc(i) +150*(xe(i)-xc(i))
      yy = yc(i) +150*(ye(i)-yc(i))
      zz = zc(i) +150*(ze(i)-zc(i))
      WRITE(789,'(a,3f17.7)') 'atom Osh  ',xx,yy,zz
    END IF
  END DO
END IF

DO i=1,nk
  IF (ionlymob == 0 .OR. imobk(i) == 1) THEN
    WRITE(789,'(a,3f17.7)') 'atom Mg  ',xk(i),yk(i),zk(i)
  END IF
END DO
#endif



WRITE(789,*) 'spec   O   0.6  0.7 0.3 0.3'
WRITE(789,*) 'spec   Mg  0.6  1.0 0.0 0.0'
WRITE(789,*) 'spec   Osh  1.4  0.0 0.0 1.0'
WRITE(789,*) 'tmat  1.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0'
WRITE(789,*) 'dist    51.629'
WRITE(789,*) 'inc      5.000'
WRITE(789,*) 'scale   15.671'
WRITE(789,*) 'rfac 1.00'
WRITE(789,*)'bfac 1.00'
WRITE(789,*) 'pos    0.000    0.000'
WRITE(789,*) 'switches 1 0 1 0 0 1 1 0 0'

CLOSE(789)

RETURN
END SUBROUTINE makexbsfile
!------------------------------------------------------------


SUBROUTINE printfieldx(iunit,field,y,z)
!------------------------------------------------------------

! Print 'field' along x-axis for grid points y and z.

USE params
IMPLICIT NONE


INTEGER, INTENT(IN)     :: iunit
REAL(DP), INTENT(IN)    :: field(kdfull2)
REAL(DP), INTENT(IN)    :: y
REAL(DP), INTENT(IN)    :: z

INTEGER :: ind, ix
REAL(DP) :: x 
INTEGER,EXTERNAL :: getnearestgridpoint

DO ix=minx,maxx
  x = (ix-nxsh)*dx
  
  ind = getnearestgridpoint(x,y,z)
  
  WRITE(iunit,'(3f12.4,2e17.7)') x,y,z,field(ind)
  
END DO

RETURN
END SUBROUTINE printfieldx
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE printfieldy(iunit,field,x,z)
!------------------------------------------------------------

! Print 'field' along y-axis for grid points x and z.

USE params
IMPLICIT NONE


INTEGER, INTENT(IN)     :: iunit
REAL(DP), INTENT(IN)    :: field(kdfull2)
REAL(DP), INTENT(IN)    :: x
REAL(DP), INTENT(IN)    :: z

INTEGER :: ind, iy
REAL(DP) :: y 
INTEGER,EXTERNAL :: getnearestgridpoint

DO iy=miny,maxy
  y = (iy-nysh)*dy
  
  ind = getnearestgridpoint(x,y,z)
  
  WRITE(iunit,'(3f12.4,2e17.7)') x,y,z,field(ind)
  
END DO


RETURN
END SUBROUTINE printfieldy
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE printfieldz(iunit,field,x,y)
!------------------------------------------------------------

! Print 'field' along z-axis for grid points x and y.

USE params
IMPLICIT NONE


INTEGER, INTENT(IN)     :: iunit
REAL(DP), INTENT(IN)    :: field(kdfull2)
REAL(DP), INTENT(IN)    :: x
REAL(DP), INTENT(IN)    :: y


INTEGER :: ind, iz
REAL(DP) :: z 
INTEGER,EXTERNAL :: getnearestgridpoint

DO iz=minz,maxz
  z = (iz-nzsh)*dz
  
  ind = getnearestgridpoint(x,y,z)
  
  WRITE(iunit,'(3f12.4,2e17.7)') x,y,z,field(ind)
  
END DO


RETURN
END SUBROUTINE printfieldz
!------------------------------------------------------------



!------------------------------------------------------------

SUBROUTINE printfield(iunit,field,comment)
!------------------------------------------------------------

! Print full 3D field with associated coordinates

USE params
IMPLICIT NONE


INTEGER, INTENT(IN)                      :: iunit
REAL(DP), INTENT(IN)                     :: field(kdfull2)
CHARACTER (LEN=*), INTENT(IN)            :: comment


LOGICAL :: topenunit
INTEGER :: ind, ix, iy, iz
REAL(DP) :: x1, y1, z1

INQUIRE(iunit,OPENED=topenunit)
IF(.NOT.topenunit) THEN
  OPEN(iunit,FILE=comment//'.'//outnam)
  WRITE(iunit,'(a)') '#    x   y  z   '//comment
END IF

ind=0
DO iz=1,maxz
  z1=(iz-nzsh)*dz
  DO iy=1,maxy
    y1=(iy-nysh)*dy
    DO ix=1,maxx
      x1=(ix-nxsh)*dx
      ind=ind+1
      WRITE(iunit,'(3f12.4,1e17.7)') x1,y1,z1,field(ind)
    END DO
  END DO
END DO

IF(.NOT.topenunit)  CLOSE(iunit)

RETURN
END SUBROUTINE printfield
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE printfield2(iunit,field1,field2)
!------------------------------------------------------------

! Print two full 3D fields with associated coordinates

USE params
IMPLICIT NONE

INTEGER, INTENT(IN)                  :: iunit
REAL(DP), INTENT(IN)                     :: field1(kdfull2)
REAL(DP), INTENT(IN)                     :: field2(kdfull2)

INTEGER :: ind, ix, iy, iz
REAL(DP) :: ft1, ft2, x1, y1, z1

ind=0
DO iz=1,maxz
  z1=(iz-nzsh)*dz
  DO iy=1,maxy
    y1=(iy-nysh)*dy
    DO ix=1,maxx
      x1=(ix-nxsh)*dx
      
      ind=ind+1
      
      ft1=field1(ind)
      ft2=field2(ind)
      
      IF (ABS(ft1) < 1D-90) ft1=0D0
      IF (ABS(ft2) < 1D-90) ft2=0D0
      
      WRITE(iunit,'(3f12.4,2e17.7)') x1,y1,z1,ft1,ft2
      
    END DO
  END DO
END DO

RETURN
END SUBROUTINE printfield2
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE printforces(iflag,iterat)
!------------------------------------------------------------

! Print forces on ions

USE params
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: iflag
INTEGER, INTENT(IN)                      :: iterat

INTEGER :: i

IF (iterat > 0) THEN
  
  WRITE(6,*)'CHECKING FORCES !!',iflag
  
  IF(iflag == 4)THEN
    DO i=1,nion
      WRITE(6,'(3e17.8)') fx(i),fy(i),fz(i)
    END DO
#if(raregas)
  ELSE IF(iflag == 1)THEN
    DO i=1,nc
      WRITE(6,'(3e17.8)') fxc(i),fyc(i),fzc(i)
    END DO
  ELSE IF(iflag == 2)THEN
    DO i=1,NE
      WRITE(6,'(3e17.8)') fxe(i),fye(i),fze(i)
    END DO
  ELSE IF(iflag == 3)THEN
    DO i=1,nk
      WRITE(6,'(3e17.8)') fxk(i),fyk(i),fzk(i)
    END DO
#endif
  END IF
  
END IF

RETURN
END SUBROUTINE printforces


!------------------------------------------------------------

SUBROUTINE getcm(iflag,iflagc,iflagk)
!------------------------------------------------------------

!    Calculates the center of mass and stores it in the vector 
!    rVecTmp(1:3) which is communicated via module 'params'.

USE params
IMPLICIT NONE
INTEGER,INTENT(IN)::iflag
INTEGER,INTENT(IN)::iflagc
INTEGER,INTENT(IN)::iflagk

INTEGER :: i
REAL(DP) :: summ, sumx, sumy, sumz

summ = 0D0
sumx = 0D0
sumy = 0D0
sumz = 0D0

!     the common scaling factor 1836.0*ame can be skipped here!

IF (iflag /= 0) THEN
  DO i=1,nion
    summ = summ + amu(np(i))
    sumx = sumx + amu(np(i))*cx(i)
    sumy = sumy + amu(np(i))*cy(i)
    sumz = sumz + amu(np(i))*cz(i)
  END DO
END IF

#if(raregas)
IF (iflagc /= 0) THEN
  DO i=1,nc
    summ = summ + mion + me
    sumx = sumx + mion*xc(i) + me*xe(i)
    sumy = sumy + mion*yc(i) + me*ye(i)
    sumz = sumz + mion*zc(i) + me*ze(i)
  END DO
END IF

IF (iflagk /= 0) THEN
  DO i=1,nk
    summ = summ + mkat
    sumx = sumx + mkat*xk(i)
    sumy = sumy + mkat*yk(i)
    sumz = sumz + mkat*zk(i)
  END DO
END IF
#endif

rvectmp(1) = sumx/summ
rvectmp(2) = sumy/summ
rvectmp(3) = sumz/summ

RETURN
END SUBROUTINE getcm
!------------------------------------------------------------


!-----emoms-------------------------------------------------------emoms

!  Various multipole moments for the given electronic distribution 
!  'rho'. Result is array of moments 'qe' communicated via module 'params'.

SUBROUTINE emoms(rho)
USE params
IMPLICIT NONE


REAL(DP), INTENT(IN)  :: rho(2*kdfull2)

INTEGER :: ind, ix, iy, iz, k
REAL(DP) :: updens, dodens
REAL(DP) :: sux, sdx, suy, sdy, suz, sdz, s, s1, s2
REAL(DP) :: x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4
!----------------------------------------------------------------------

sux=0D0
sdx=0D0
suy=0D0
sdy=0D0
suz=0D0
sdz=0D0
s1=0D0
s2=0D0

nrmom=35
IF(nrmom > kmom) STOP ' too many moments in EMOMS'

DO k=1,nrmom
  qe(k)=0D0
END DO


!     switch for calculating moments relative to center of mass (1)
!     or center of box (0)

rvectmp = 0D0
IF(iemomsrel == 1 .AND. nion2 > 0) CALL getcm(1,0,0)

ind=0
DO iz=minz,maxz
  z1=(iz-nzsh)*dz
  z1=z1-rvectmp(3)
  
  z2=z1*z1
  z3=z2*z1
  z4=z2*z2
  DO iy=miny,maxy
    y1=(iy-nysh)*dy
    y1=y1-rvectmp(2)
    
    y2=y1*y1
    y3=y2*y1
    y4=y2*y2
    DO ix=minx,maxx
      ind=ind+1
!test               if((ix.ne.nx2).and.(iy.ne.ny2).and.(iz.ne.nz2)) then
      IF((ix <= nx2).AND.(iy <= ny2).AND.(iz <= nz2)) THEN
        x1=(ix-nxsh)*dx
        x1=x1-rvectmp(1)
        
        x2=x1*x1
        x3=x2*x1
        x4=x2*x2
        s=rho(ind)
!                                                     monopole
        qe(1)=qe(1)+s
!                                                     dipole
        qe(2)=qe(2)+s*x1
        qe(3)=qe(3)+s*y1
        qe(4)=qe(4)+s*z1
!                                                     quadrupole
        qe(5)=qe(5)+s*x2
        qe(6)=qe(6)+s*y2
        qe(7)=qe(7)+s*z2
        qe(8)=qe(8)+s*x1*y1
        qe(9)=qe(9)+s*z1*x1
        qe(10)=qe(10)+s*z1*y1
!                                                     octupole
        qe(11)=qe(11)+s*x3
        qe(12)=qe(12)+s*x2*y1
        qe(13)=qe(13)+s*z1*x2
        qe(14)=qe(14)+s*x1*y2
        qe(15)=qe(15)+s*x1*y1*z1
        qe(16)=qe(16)+s*z2*x1
        qe(17)=qe(17)+s*y3
        qe(18)=qe(18)+s*z1*y2
        qe(19)=qe(19)+s*z2*y1
        qe(20)=qe(20)+s*z3
!                                                     hexadecupole
        qe(21)=qe(21)+s*x4
        qe(22)=qe(22)+s*x3*y1
        qe(23)=qe(23)+s*x3*z1
        qe(24)=qe(24)+s*x2*y2
        qe(25)=qe(25)+s*z1*x2*y1
        qe(26)=qe(26)+s*x2*z2
        qe(27)=qe(27)+s*x1*z3
        qe(28)=qe(28)+s*z1*x1*y2
        qe(29)=qe(29)+s*z2*x1*y1
        qe(30)=qe(30)+s*z3*x1
        qe(31)=qe(31)+s*y4
        qe(32)=qe(32)+s*z1*y3
        qe(33)=qe(33)+s*z2*y2
        qe(34)=qe(34)+s*z3*y1
        qe(35)=qe(35)+s*z4
!                                                    spin dipole
        updens = 0.5D0*(rho(ind)*rho(ind+nxyz)+rho(ind))
        sux = sux+updens*x1
        suy = suy+updens*y1
        suz = suz+updens*z1
        s1 = s1+updens
        dodens = -0.5D0*(rho(ind)*rho(ind+nxyz)-rho(ind))
        sdx = sdx+dodens*x1
        sdy = sdy+dodens*y1
        sdz = sdz+dodens*z1
        s2 = s2+dodens !integral over downdensity
      END IF
    END DO
  END DO
END DO

s1=s1*dvol
s2=s2*dvol

!      endif
!      write(6,'(a,2f12.7)') 'total number of spinup/down=',s1,s2

se(1)=sux/(s1)-sdx/(s2)   !spindipole-x
se(2)=suy/(s1)-sdy/(s2)   !spindipole-y
se(3)=suz/(s1)-sdz/(s2)   !spindipole-z
se(4)=s1                  ! total number spin-up
se(5)=s2                  ! total number spin-down

DO k=1,nrmom
  qe(k)=qe(k)*dvol
END DO

DO k=2,nrmom
  qe(k)=qe(k)/qe(1)      !normalization
END DO
!      write(6,'(a,6(f9.3))') 'monop.,dipole=',qe(1),qe(2),qe(3),qe(4)
!      write(6,'(a,6(f10.4))')
!     &   'up/down,mon.,dip,=',s1,s2,qe(1),qe(2),qe(3),qe(4)
!      write(7,'(a,f8.4,a,f7.2,a,3f9.2)')
!     &  'In emoms: t=',tfs,' moments: monop.=',qe(1),
!     &  ' dip.: ',qe(2),qe(3),qe(4)
!      write(6,'(a,f8.4,a,f7.2,a,3f9.2)')
!     &  'In emoms: t=',tfs,' moments: monop.=',qe(1),
!     &  ' dip.: ',qe(2),qe(3),qe(4)

DO k=1,3
  se(k)=se(k)*dvol
END DO


RETURN
END SUBROUTINE emoms


!       *********************

SUBROUTINE laserp(vlaser,rho)

!       **********************

!      computes a laser potential corresponding to an electric
!      field polarized in unit direction given by ux,uy,uz

!      laser potential is :

!      vlaser (x,y,z) = e_0 (x*ux + y*uy + z*uz) f(t) sin(omega*t)

!      caution vlaser is stored as a 1-d vector (as any other grid field)

!      the energy absorbed from the laser field is accumulated
!      on the variable 'elaser' (stored in module 'params')

!      photon pulsation is          omega
!      electric field intensity is  e_0

!      the laser pulse is characterized by
!          its length      deltat
!          and peak time   tpeak
!          and begins at   tnode

!      pulse profile is given by f(t) (t is time)
!      2 options are presently implemented :

!      itft = 1 --> ramp laser pulse, sine switching on/off up to/from f(t) = 1

!          tstart = tnode  + tpeak
!          tend   = tstart + deltat
!          tvend  = tend   + tpeak

!          0.     < t < tnode    f(t) = 0.
!          tnode  < t < tstart   f(t) = sin(.5 *PI * (t-tnode) / tpeak)
!          tstart < t < tend     f(t) = 1.
!          tend   < t < tvend    f(t) = sin(.5 *PI* (t - tend) / tpeak)
!          tvend  < t < ...      f(t) = 0.

!      itft = 2 --> Gaussian laser pulse
!          tmax = tnode + tpeak
!          0.     < t < tnode    f(t) = 0.
!          tnode  < t < ...      f(t) = exp(-((t-tmax)/deltat)**2)

!      itft = 3 --> cos**2 envelope
!          tnode   = start of pulse
!          tmax    = tnode + 2*tpeak    = end of pulse
!          2*tpeak = peak length
!          deltat  = obsolete, must be set to zero
!
!      this pulse allows a second cos**2 pulse (probe) switched by
!      start time tstart2 and length 2*tpeak2.


!      itft = 4 --> cos**4 envelope
!          tnode   = start of pulse
!          tmax    = tnode + 2*tpeak    = end of pulse
!          2*tpeak = peak length
!          deltat  = obsolete, must be set to zero


! The routine accumulates the integrated pulse profile in 'fpulseinteg1/2'.
! Trapezoidal integration is used and the last pulse value is
! saved in 'foft1/2old'.

USE params
IMPLICIT NONE

REAL(DP), INTENT(OUT)  :: vlaser(kdfull2)
REAL(DP), INTENT(IN)   :: rho(2*kdfull2)

INTEGER :: ind, ix, iy, iz
REAL(DP) :: ascal,acc1,acc2, ppower
REAL(DP) :: foft,foft1,foft2
REAL(DP) :: tstart,tend,tmax,tpulse2,tvend,snorm,scal1,scal2
REAL(DP) :: ex, ey, ez, x1, y1, z1
IF (ABS(e0) <= 1D-20) THEN
  ind = 0
  DO  iz=minz,maxz
    DO  iy=miny,maxy
      DO  ix=minx,maxx
        ind  = ind + 1
        vlaser(ind) = 0.0D0
      END DO
    END DO
  END DO
  RETURN
END IF


IF(ilas == 0) THEN     ! switch 'ilas' is maintained in 'params.F90'
  ilas=1
  ascal=e1x*e2x+e1y*e2y+e1z*e2z
!  sx=e1x+e2x*COS(phi)
!  sy=e1y+e2y*COS(phi)
!  sz=e1z+e2z*COS(phi)
!  snorm=sx*sx+sy*sy+sz*sz
  snorm = e1x*e1x+e1y*e1y+e1z*e1z+e2x*e2x+e2y*e2y+e2z*e2z
  snorm=SQRT(snorm)
  e1x=e1x/snorm
  e1y=e1y/snorm
  e1z=e1z/snorm
  e2x=e2x/snorm
  e2y=e2y/snorm
  e2z=e2z/snorm
  elaser = 0D0
  acc1old = 0D0
  acc2old = 0D0
  foft1old = 0D0
  foft2old = 0D0
  fpulseinteg1 = 0D0
  fpulseinteg2 = 0D0
  timeold = 0D0
END IF

!     prepare time profile

foft2 = 0D0                ! standard case is: no second pulse

IF(itft == 1) THEN
  tstart = tnode  + tpeak
  tend   = tstart + deltat
  tvend  = tend   + tpeak
  
  IF(tfs <= tnode.OR.tfs >= tvend) foft = 0D0
!mb         if(tfs.gt.tnode.and.tfs.le.tstart)
!mb     &        foft = sin (.5 * pi * tfs / tpeak)
!mb it should rather be:
  IF(tfs > tnode.AND.tfs <= tstart) foft = SIN(0.5D0*pi*(tfs-tnode)/tpeak)
  IF(tfs > tstart.AND.tfs <= tend) foft = 1D0
  IF(tfs > tend.AND.tfs < tvend) foft = SIN (0.5D0 * pi * (tvend-tfs) / tpeak)
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
    foft = COS((-0.5D0+(tfs-tnode)/(2D0*tpeak+deltat))*pi)**2
  END IF

  IF(ABS(e0_2) > small) THEN
    tpulse2 = tpeak2+tpeak2
    IF(tfs >= tstart2 .AND. tfs <= tstart2+tpulse2) THEN
      foft2 = COS((-0.5D0+(tfs-tstart2)/tpulse2)*pi)**2
    END IF
  END IF

END IF

IF(itft == 4) THEN
  tvend = tnode + 2*tpeak + deltat
  IF(tfs <= tnode.OR.tfs >= tvend) THEN
    foft = 0D0
  ELSE
    foft = COS((-0.5D0+(tfs-tnode)/(2*tpeak+deltat))*pi)**4
  END IF
END IF

IF(itft <= 0.OR.itft >= 5) THEN
  STOP ' this pulse profile not yet implemented'
END IF
power=e0*e0*foft*foft

!     implement time profile in space


foft1 = COS(omega*tfs/0.0484D0)*foft
foft2 = COS(omega2*tfs/0.0484D0+phase2)*foft2
fpulseinteg1 = fpulseinteg1 + (foft1+foft1old)*0.5D0*(tfs/0.0484D0-timeold)
fpulseinteg2 = fpulseinteg2 + (foft2+foft2old)*0.5D0*(tfs/0.0484D0-timeold)

      write(7,'(a,2f9.3,2(1pg12.4))') &
         ' tfs,tRy,foft1,foft2=',tfs,tfs/0.048D0,foft1,foft2
ind = 0
acc1 = 0D0
acc2 = 0D0
DO  iz=minz,maxz
  z1=(iz - nzsh)*dz
  DO  iy=miny,maxy
    y1=(iy - nysh)*dy
    DO  ix=minx,maxx
      x1=(ix - nxsh)*dx
      ind  = ind + 1
      scal1 = x1 * e1x + y1 * e1y + z1 * e1z
      scal2 = x1 * e2x + y1 * e2y + z1 * e2z
      vlaser(ind) = - e0 * (scal1*foft1+scal2*foft2)
      acc1 = acc1 + e0*scal1*rho(ind)
      acc2 = acc2 + e0*scal2*rho(ind)
    END DO
  END DO
END DO
acc1 = acc1*dvol
acc2 = acc2*dvol
elaser = elaser + (acc1-acc1old)*(foft1+foft1old)*0.5D0 & 
                + (acc2-acc2old)*(foft2+foft2old)*0.5D0 
!!test:
! WRITE(*,*) ' laser:',acc1,acc1old,foft1,foft1old, &
!                   (acc1-acc1old)*(foft1+foft1old)*0.5D0 & 
!                 + (acc2-acc2old)*(foft2+foft2old)*0.5D0,elaser
acc1old = acc1
acc2old = acc2
foft1old = foft1
foft2old = foft2 
timeold = tfs/0.0484D0

ex=e1x*foft1+e2x*foft2
ey=e1y*foft1+e2y*foft2
ez=e1z*foft1+e2z*foft2

ppower = e0*e0*(ex*ex+ey*ey+ez*ez)


WRITE(38,'(1f15.5,7(1pg14.5))') tfs,ex,ey,ez,ppower,elaser, &
                                fpulseinteg1,fpulseinteg2



RETURN
END SUBROUTINE laserp


!-----projectp-----------------------------------------------------

SUBROUTINE projectp(Vproj)

! Calculates the potential 'Vproj' from the point charge projectile
! on the valence electrons.

USE params
IMPLICIT NONE

REAL(DP), INTENT(OUT) :: Vproj(kdfull2)

INTEGER :: ind, ix, iy, iz
REAL(DP) :: x1, y1, z1
!-----------------------------------------------------------------

    
IF (ABS(projcharge) < 1D-3) THEN
   ind = 0
   DO  iz=minz,maxz
      DO  iy=miny,maxy
         DO  ix=minx,maxx
            ind  = ind + 1
            Vproj(ind) = 0D0
         ENDDO
      ENDDO
   ENDDO
ELSE
   ind = 0
   DO  iz=minz,maxz
      z1=(iz - nzsh)*dz
      DO  iy=miny,maxy
         y1=(iy - nysh)*dy
         DO  ix=minx,maxx
            x1=(ix - nxsh)*dx
            ind  = ind + 1
            Vproj(ind) = - projcharge*e2/sqrt &
                 ((x1-(projvelx*tfs+projinix))**2 + &
                 (y1-(projvely*tfs+projiniy))**2 + &
                 (z1-(projvelz*tfs+projiniz))**2)
         ENDDO
      ENDDO
   ENDDO
ENDIF
 
RETURN
END SUBROUTINE projectp

!     **************************

REAL(DP) FUNCTION c_wfnorm(psi)

!  Norm of wavefunction 'psi'.     COMPLEX psi version

USE params
IMPLICIT NONE
COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2)

INTEGER :: ii
REAL(DP) :: acc
acc = 0D0
DO ii=1,nxyz
  acc   = REAL(psi(ii),DP)*REAL(psi(ii),DP) +AIMAG(psi(ii))*AIMAG(psi(ii)) + acc
END DO
c_wfnorm = acc*dvol
RETURN
END FUNCTION c_wfnorm

!     **************************

REAL(DP) FUNCTION r_wfnorm(psi1)

!  Norm of wavefunction 'psi'.     REAL psi version

USE params
IMPLICIT NONE
REAL(DP), INTENT(IN)                         :: psi1(kdfull2)

INTEGER :: ii 
REAL(DP):: acc

acc = 0D0
DO ii=1,nxyz
  acc = (psi1(ii))*psi1(ii)  + acc
END DO
r_wfnorm = acc*dvol
RETURN
END FUNCTION r_wfnorm

!     **************************

COMPLEX(DP) FUNCTION c_wfovlp(psi1,psi2)

!     *************************

!      Overlap   <psi1|psi2>
!         COMPLEX version

USE params
IMPLICIT NONE

COMPLEX(DP), INTENT(IN)      :: psi1(kdfull2)
COMPLEX(DP), INTENT(IN)      :: psi2(kdfull2)

INTEGER :: ii
COMPLEX(DP) :: csum

csum = 0D0
DO ii=1,nxyz
  csum = CONJG(psi1(ii))*psi2(ii)  + csum
END DO
c_wfovlp = csum*dvol
RETURN
END FUNCTION c_wfovlp
!     **************************

REAL(DP) FUNCTION r_wfovlp(psi1,psi2)

!     *************************

!      Overlap   <psi1|psi2>
!         REAL version

USE params
IMPLICIT NONE

REAL(DP), INTENT(IN)                     :: psi1(kdfull2)
REAL(DP), INTENT(IN)                     :: psi2(kdfull2)

INTEGER :: ii
REAL(DP) :: acc

acc = 0D0
DO ii=1,nxyz
  acc = (psi1(ii))*psi2(ii)  + acc
END DO
r_wfovlp = acc*dvol

RETURN
END FUNCTION r_wfovlp


SUBROUTINE rhointxy(rho,itime)

!       **********************

!      provides data file to analyse time-evolution of density.
!      density 'rho' is integrated over x and y .
!      the integrated density versus z is written on file 29
!      for later analysis.
!      at present a print modulus of 10 is hardwired.

USE params
IMPLICIT NONE

REAL(DP), INTENT(IN)   :: rho(2*kdfull2)
INTEGER, INTENT(IN)    :: itime

INTEGER :: ind, ix, iy, iz
REAL(DP) ::acc
REAL(DP) :: rhoint(minz:maxz)
LOGICAL :: tfirst=.true.
!DATA tfirst/.true./
!      data modrho/10/
IF(modrho == 0) THEN
  STOP 'modrho must be larger than zero'
END IF
IF(MOD(itime,modrho) /= 0) RETURN
time=itime*dt1
tfs=time*0.0484D0
IF(tfirst) THEN
  OPEN(UNIT=28,FORM='unformatted',FILE='rhointxy')
  WRITE(28) minz,maxz,nzsh,dz,centfz,dt1,itime,modrho
  tfirst = .false.
  WRITE(6,*) ' output for integrated rho initialised'
END IF

ind = 0
DO  iz=minz,maxz
  acc   =  0D0
  DO  iy=miny,maxy
    DO  ix=minx,maxx
      ind  = ind + 1
      acc   = rho(ind) + acc
    END DO
  END DO
  rhoint(iz) = acc*dx*dy
END DO
WRITE(28) tfs,rhoint

RETURN
END SUBROUTINE rhointxy

!       *********************

SUBROUTINE rhointyz(rho,itime)

!       **********************

!      provides data file to analyse time-evolution of density.
!      density 'rho' is integrated over x and y .
!      the integrated density versus z is written on file 29
!      for later analysis.
!      at present a print modulus of 10 is hardwired.

USE params
IMPLICIT NONE

REAL(DP), INTENT(IN)   :: rho(2*kdfull2)
INTEGER, INTENT(IN)    :: itime



LOGICAL :: tfirst
INTEGER :: ind, ix, iy, iz
REAL(DP) :: acc
REAL(DP) :: rhoint(minx:maxx)


DATA tfirst/.true./
!      data modrho/10/

IF(MOD(itime,modrho) /= 0) RETURN
time=itime*dt1
tfs=time*0.0484D0
IF(tfirst) THEN
  OPEN(UNIT=29,FORM='unformatted',FILE='rhointyz')
  WRITE(29) minx,maxx,nxsh,dx,centfx,dt1,itime,modrho
  tfirst = .false.
  WRITE(6,*) ' output for integrated rho initialised'
END IF

ind = 0
DO  ix=minx,maxx
  acc   =  0D0
  DO  iz=minz,maxz
    DO  iy=miny,maxy
      ind  =(iz-1)*nxyf+(iy-1)*nyf+ix
      acc   = rho(ind) + acc
    END DO
  END DO
  rhoint(ix) = acc*dz*dy
END DO
WRITE(29) tfs,rhoint

RETURN
END SUBROUTINE rhointyz

!       *********************

SUBROUTINE rhointxz(rho,itime)

!       **********************

!      provides data file to analyse time-evolution of density.
!      density 'rho' is integrated over x and y .
!      the integrated density versus z is written on file 29
!      for later analysis.
!      at present a print modulus of 10 is hardwired.

USE params
IMPLICIT NONE

REAL(DP), INTENT(IN) :: rho(2*kdfull2)
INTEGER, INTENT(IN)  :: itime

LOGICAL :: tfirst
INTEGER :: ind, ix, iy, iz
REAl(DP) :: acc
REAL(DP) :: rhoint(miny:maxy)

DATA tfirst/.true./
!      data modrho/10/

IF(MOD(itime,modrho) /= 0) RETURN
time=itime*dt1
tfs=time*0.0484D0
IF(tfirst) THEN
  OPEN(UNIT=27,FORM='unformatted',FILE='rhointxz')
  WRITE(27) miny,maxy,nysh,dy,centfy,dt1,itime,modrho
  tfirst = .false.
  WRITE(6,*) ' output for integrated rho initialised'
END IF

ind = 0
DO  iy=miny,maxy
  acc   =  0D0
  DO  iz=minz,maxz
    DO  ix=minx,maxx
      ind  =(iz-1)*nxyf+(iy-1)*nyf+ix
      acc   = rho(ind) + acc
    END DO
  END DO
  rhoint(iy) = acc*dx*dz
END DO
WRITE(27) tfs,rhoint

RETURN
END SUBROUTINE rhointxz


!       **************************

SUBROUTINE prifld(field,comment)

! Print 'field' along x-axis for y=0 and z=0.

USE params
IMPLICIT NONE

REAL(DP), INTENT(IN)             :: field(kdfull2)
CHARACTER (LEN=*), INTENT(IN)   :: comment

LOGICAL :: tnopri
INTEGER :: ind, jx, jy, jz

DATA tnopri/.false./

!c       print field along x-axis

IF(tnopri) RETURN
WRITE(6,'(/2a)') 'field along x:',comment
ind = 0
DO jz=minz,maxz
  DO jy=miny,maxy
    DO jx=minx,maxx
      ind    = 1 + ind
      IF(jz == nzsh .AND. jy == nysh)  &
           WRITE(6,'(1x,f6.2,g13.5)') (jx-nxsh)*dx,field(ind)
    END DO
  END DO
END DO

RETURN
END SUBROUTINE prifld

SUBROUTINE prifldz(field,comment)

! Print 'field' along z-axis for x=0 and y=0.

USE params
IMPLICIT NONE

REAL(DP), INTENT(IN)             :: field(kdfull2)
CHARACTER (LEN=*), INTENT(IN)   :: comment

LOGICAL :: tnopri
INTEGER :: ind, jx, jy, jz
DATA tnopri/.false./

!c       print field along z-axis

IF(tnopri) RETURN
WRITE(6,'(/2a)') 'field along z:',comment
ind = 0
DO jz=minz,maxz
  DO jy=miny,maxy
    DO jx=minx,maxx
      ind    = 1 + ind
      IF(jx == nxsh .AND. jy == nysh)  &
          WRITE(6,'(1x,f6.2,g13.5)') (jz-nzsh)*dx,field(ind)
    END DO
  END DO
END DO

RETURN
END SUBROUTINE prifldz



!       **************************

SUBROUTINE prifld2(iunit,field,comment)

! Print 'field' along y-axis for z=0 and y=0.
! 'iunit' is the file unti at which output is written.

USE params
IMPLICIT NONE

INTEGER, INTENT(IN)              :: iunit
REAL(DP), INTENT(IN)             :: field(kdfull2)
CHARACTER (LEN=*), INTENT(IN)   :: comment

LOGICAL :: tnopri
INTEGER :: ind, jx, jy, jz
DATA tnopri/.false./

!c       print field along x-axis

IF(tnopri) RETURN
WRITE(iunit,'(/2a)') 'field along x:',comment
ind = 0
DO jz=minz,maxz
  DO jy=miny,maxy
    DO jx=minx,maxx
      ind    = 1 + ind
      IF(jz == nzsh .AND. jy == nysh)  &
          WRITE(iunit,'(1x,f6.2,g13.5)') (jx-nxsh)*dx,field(ind)
    END DO
  END DO
END DO

RETURN
END SUBROUTINE prifld2




!       **************************

SUBROUTINE mtv_fld(field,i)

! Prepares plotting of 'field' with MTV software.

USE params
IMPLICIT NONE

INTEGER, INTENT(IN)              :: i  ! 1=static, 2=dynamic (only change is filenaming)
REAL(DP), INTENT(IN)             :: field(kdfull2)

INTEGER :: ind, jx, jy, num


!     print 2d-field

num = nz
ind = nxyf*num
IF(i == 1) THEN
  OPEN(70,FILE='denstat.'//outnam,STATUS='unknown')
ELSE IF(i == 2) THEN
  OPEN(70,FILE='densdyn.'//outnam,STATUS='unknown')
END IF
WRITE(70,*) '$ data=contour'
WRITE(70,'(a,i4,a,f7.2,a,f7.2)')  &
    '% nx =',nx2,' xmin=',(minx-nxsh)*dx,' xmax=',(maxx-nxsh)*dx
WRITE(70,'(a,i4,a,f7.2,a,f7.2)')  &
    '% ny =',ny2,' ymin=',(miny-nysh)*dy,' ymax=',(maxy-nysh)*dy
WRITE(70,*) '% contstyle=2'
WRITE(70,*) '% nsteps=50'
WRITE(70,*) '% interp=2'
WRITE(70,*) '% xlabel="x"'
WRITE(70,*) '% ylabel="y"'
WRITE(70,*) '% toplabel="mg2"'
WRITE(70,*) '% equalscale=false'
WRITE(70,*) '% fitpage=false'
WRITE(70,*) '% xyratio=1.0'

DO jy=miny,maxy
  DO jx=minx,maxx
    ind = ind + 1
    WRITE(70,'(g13.5)') field(ind)
  END DO
  WRITE(70,*) ' '
END DO

WRITE(70,*) '$ end'
CLOSE(70)

RETURN
END SUBROUTINE mtv_fld

!     ****************************

SUBROUTINE fmtv_fld(psi,field,i)

! Prepares plotting of 'psi**2' with MTV software.

USE params
IMPLICIT NONE
CHARACTER (LEN=6) :: ext
CHARACTER (LEN=8) :: ext1

REAL(DP), INTENT(IN)        :: field(kdfull2)
COMPLEX(DP), INTENT(IN)     :: psi(kdfull2,kstate)
INTEGER, INTENT(IN)         :: i

INTEGER :: ind, ishift, j, jx, jy, jz, k, nb
REAL(DP) :: rhoz, xmin1, ymin1, zmin1, xmax1, ymax1, zmax1
REAL(DP), ALLOCATABLE :: rho(:)

ALLOCATE(rho(2*kdfull2))

!     print integrated 3d-field (dimension to integrate on depends of values of i3dx, i3dz)

xmin1=-nx*dx
ymin1=-ny*dy
zmin1=-nz*dz
xmax1=nx*dx
ymax1=ny*dy
zmax1=nz*dz
WRITE(ext,'(a,i5.5)') '.',i

IF(i3dz == 1) THEN
  OPEN(70,FILE='zdensdyn'//ext//'.'//outnam,STATUS='unknown')
  WRITE(70,'(a)') '$ data=contour'
  WRITE(70,'(a,i2,a,f7.3,a,f5.2)') '% nx = ',nx2,  &
      ' xmin=',xmin1,' xmax=',xmax1
  WRITE(70,'(a,i2,a,f7.3,a,f5.2)') '% ny = ',ny2,  &
      ' ymin=',ymin1,' ymax=',ymax1
  WRITE(70,'(a)') '% contstyle=2'
  WRITE(70,'(a)') '% nsteps=50'
  WRITE(70,'(a)') '% interp=2'
  WRITE(70,'(a)') '% xlabel="x"'
  WRITE(70,'(a)') '% ylabel="y"'
  WRITE(70,'(a,a3,a)') '% toplabel="',outnam,'"'
  WRITE(70,'(a)') '% equalscale=false'
  WRITE(70,'(a)') '% fitpage=false'
  WRITE(70,'(a)') '% xyratio=1.0'
  ind = 0
  DO jy=miny,maxy
    DO jx=minx,maxx
      ind = ind + 1
      j = 0
      rhoz = 0D0
      DO jz=minz,maxz
        j=j+1
        k = (j-1)*nxyf + ind
        rhoz = rhoz + field(k)
      END DO
      WRITE(70,'(g13.5)') rhoz
    END DO
    WRITE(70,*) ' '
  END DO
  WRITE(70,'(a)') '$ end'
  CLOSE(70)
END IF

IF(i3dx == 1) THEN
  OPEN(73,FILE='xdensdyn'//ext//'.'//outnam,STATUS='unknown')
  WRITE(73,'(a)') '$ data=contour'
  WRITE(73,'(a,i2,a,f7.3,a,f5.2)') '% nx = ',ny2,  &
      ' xmin=',ymin1,' xmax=',ymax1
  WRITE(73,'(a,i2,a,f7.3,a,f5.2)') '% ny = ',nz2,  &
      ' ymin=',zmin1,' ymax=',zmax1
  WRITE(73,'(a)') '% contstyle=2'
  WRITE(73,'(a)') '% nsteps=50'
  WRITE(73,'(a)') '% interp=2'
  WRITE(73,'(a)') '% xlabel="y"'
  WRITE(73,'(a)') '% ylabel="z"'
  WRITE(73,'(a,a3,a)') '% toplabel="',outnam,'"'
  WRITE(73,'(a)') '% equalscale=false'
  WRITE(73,'(a)') '% fitpage=false'
  WRITE(73,'(a)') '% xyratio=1.0'
  
  ind = 0
  DO jz=minz,maxz
    DO jy=miny,maxy
      rhoz=0D0
      DO jx=minx,maxx
        ind=ind+1
        rhoz=rhoz+field(ind)
      END DO
      WRITE(73,'(g13.5)') rhoz
    END DO
    WRITE(73,*) ' '
  END DO
  WRITE(73,'(a)') '$ end'
  CLOSE(73)
END IF

IF(i3dstate == 1) THEN
  DO nb=1,nstate
    ishift = (ispin(nrel2abs(nb))-1)*nxyz
    DO ind=1,nxyz
      rho(ind+ishift)= occup(nb)*(CONJG(psi(ind,nb)))*psi(ind,nb)
    END DO
    WRITE(ext1,'(a1,i1,a1,i5.5)') '.',nb,'.',i
    OPEN(71,FILE='zdenspsi'//ext1//'.'//outnam,STATUS='unknown')
    WRITE(71,'(a)') '$ data=contour'
    WRITE(71,'(a,i2,a,f7.3,a,f5.2)') '% nx = ',nx2,  &
        ' xmin=',xmin1,' xmax=',xmax1
    WRITE(71,'(a,i2,a,f7.3,a,f5.2)') '% ny = ',ny2,  &
        ' ymin=',ymin1,' ymax=',ymax1
    WRITE(71,'(a)') '% contstyle=2'
    WRITE(71,'(a)') '% nsteps=50'
    WRITE(71,'(a)') '% interp=2'
    WRITE(71,'(a)') '% xlabel="x"'
    WRITE(71,'(a)') '% ylabel="y"'
    WRITE(71,'(a,a3,a,i1,a)') '% toplabel="',outnam,'_',nb,'"'
    WRITE(71,'(a)') '% equalscale=false'
    WRITE(71,'(a)') '% fitpage=false'
    WRITE(71,'(a)') '% xyratio=1.0'
    ind=0
    DO jy=miny,maxy
      DO jx=minx,maxx
        ind=ind+1
        j=0
        rhoz=0D0
        DO jz=minz,maxz
          j=j+1
          k=(j-1)*nxyf+ind
          rhoz=rhoz+rho(k)
        END DO
        WRITE(71,'(g13.5)') rhoz
      END DO
      WRITE(71,*) ' '
    END DO
    WRITE(71,'(a)') '$ end'
    CLOSE(71)
    
    OPEN(72,FILE='xdenspsi'//ext1//'.'//outnam,STATUS='unknown')
    WRITE(72,'(a)') '$ data=contour'
    WRITE(72,'(a,i2,a,f7.3,a,f5.2)') '% nx = ',ny2,  &
        ' xmin=',ymin1,' xmax=',ymax1
    WRITE(72,'(a,i2,a,f7.3,a,f5.2)') '% ny = ',nz2,  &
        ' ymin=',zmin1,' ymax=',zmax1
    WRITE(72,'(a)') '% contstyle=2'
    WRITE(72,'(a)') '% nsteps=50'
    WRITE(72,'(a)') '% interp=2'
    WRITE(72,'(a)') '% xlabel="y"'
    WRITE(72,'(a)') '% ylabel="z"'
    WRITE(72,'(a,a3,a,i1,a)') '% toplabel="',outnam,'_',nb,'"'
    WRITE(72,'(a)') '% equalscale=false'
    WRITE(72,'(a)') '% fitpage=false'
    WRITE(72,'(a)') '% xyratio=1.0'
    ind=0
    DO jz=minz,maxz
      DO jy=miny,maxy
        rhoz=0D0
        DO jx=minx,maxx
          ind=ind+1
          rhoz=rhoz+rho(ind)
        END DO
        WRITE(72,'(g13.5)') rhoz
      END DO
      WRITE(72,*) ' '
    END DO
    WRITE(72,'(a)') '$ end'
    CLOSE(72)
  END DO
END IF

DEALLOCATE(rho)


RETURN
END SUBROUTINE fmtv_fld
! **************************************************

SUBROUTINE mtv_3dfld(rho,iframe)

! Prepares 3D plotting of 'rho' with MTV software.

USE params

IMPLICIT NONE

INTEGER, INTENT(IN) :: iframe
REAL(DP),INTENT(IN) :: rho(2*kdfull2)

CHARACTER (LEN=5) :: ext
INTEGER :: ind, ix, iy, iz
!REAL(DP) :: sumr


!sumr=0
WRITE(ext,'(a,i4.4)') '.',iframe
!       OPEN(unit=80,STATUS='UNKNOWN',FILE='xdens'//ext//'.'//outnam)
!       OPEN(unit=81,STATUS='UNKNOWN',FILE='ydens'//ext//'.'//outnam)
OPEN(UNIT=82,STATUS='UNKNOWN',FILE='dens'//ext//'.'//outnam)
!       OPEN(unit=83,STATUS='UNKNOWN',FILE='xydens'//ext//'.'//outnam)
ind=0
DO iz=minz,maxz
!  sumr=0
  DO iy=miny,maxy
    DO ix=minx,maxx
      ind=ind+1
!                if (ix.eq.24) then
!                   write(80,*)Rho(ind)
!                endif
!                if (iy.eq.24) then
!                   write(81,*)Rho(ind)
!                endif
!                if (iz.eq.24) then
      WRITE(82,*)rho(ind)
!                endif
!                sumr=sumr+rho(ind)
    END DO
  END DO
!         write(83,*)sumr
END DO
!       close(80)
!       close(81)
CLOSE(82)
!       close(83)

RETURN
END SUBROUTINE mtv_3dfld




!     ******************************

SUBROUTINE orthogwfr(wfr,isp,q0)

!     ******************************

!     Orthogonalizes real wavefunction 'wfr' with spin label 'isp'
!     on the set of occupied states in 'q0'.
!     Occupation numbers must be 0 or 1, i.e. 'temp=0.0' required.

USE params
IMPLICIT NONE

REAL(DP), INTENT(IN)         :: q0(kdfull2,kstate)
REAL(DP), INTENT(IN OUT)     :: wfr(kdfull2)
INTEGER, INTENT(IN)          :: isp

INTEGER :: i,nbe
REAL(DP) :: overlap
!*********************************************************

!     orthogonalization

IF(temp /= 0D0) STOP ' ORTHOGWFR requires temperature 0'
DO nbe=1,nstate
  IF(ispin(nbe) == isp .AND. occup(nbe) > 0.999D0) THEN
    overlap = wfovlp(wfr,q0(:,nbe))
    DO  i=1,nxyz
      wfr(i)=wfr(i)-overlap*q0(i,nbe)
    END DO
  END IF
END DO

RETURN
END SUBROUTINE orthogwfr







!------------------------------------------------------------
REAL(DP) FUNCTION ran1(idum)
!------------------------------------------------------------
USE params, ONLY: DP
IMPLICIT NONE

INTEGER, INTENT(IN OUT)::idum
INTEGER :: ia,im,iq,ir,ntab,ndiv
REAL(DP) :: am,eps,rnmx
PARAMETER (ia=16807,im=2147483647,am=1D0/im,iq=127773,ir=2836,  &
    ntab=32,ndiv=1+(im-1)/ntab,eps=1.2D-7,rnmx=1D0-eps)

!     Minimal random number generator of Park and Miller with Bays-Durham
!     shuffle and added safeguards. Returns a uniform random deviate
!     between 0.0 and 1.0 (exclusive of the endpoints). Call with idum
!     a negative integer to initialize; thereafter, do not alter
!     idum between successive deviates in a sequence. RNMX should
!     approximate the largest floating value that is less than 1.


INTEGER :: j,k,iv(ntab),iy
SAVE iv,iy
DATA iv /ntab*0/,iy /0/

IF (idum <= 0 .OR. iy == 0) THEN
  idum=MAX(-idum,1)
  DO j=ntab+8,1,-1
    k=idum/iq
    idum=ia*(idum-k*iq)-ir*k
    IF (idum < 0) idum=idum+im
    IF (j <= ntab) iv(j)=idum
  END DO
  iy=iv(1)
END IF
k=idum/iq
idum=ia*(idum-k*iq)-ir*k
IF (idum < 0) idum=idum+im
j=1+iy/ndiv
iy=iv(j)
iv(j)=idum
ran1=MIN(am*iy,rnmx)

RETURN
END FUNCTION ran1
!------------------------------------------------------------


!------------------------------------------------------------
REAL(DP) FUNCTION gasdev(idum)

! Random Gaussian distribution

USE params, ONLY: DP
IMPLICIT NONE

INTEGER,INTENT(IN OUT) :: idum

INTEGER :: iset
REAL(DP) :: fac,gset,rsq,v1,v2
SAVE iset,gset

DATA iset/0/
IF (iset == 0) THEN
  rsq=2D0  ! some value higher than 1.0 to start to enter the loop
  DO WHILE(rsq >= 1D0 .OR. rsq == 0D0)
    v1=2D0*ran1(idum)-1D0
    v2=2D0*ran1(idum)-1D0
    rsq=v1**2+v2**2
  END DO
  fac=SQRT(-2D0*LOG(rsq)/rsq)
  gset=v1*fac
  gasdev=v2*fac
  iset=1
ELSE
  gasdev=gset
  iset=0
END IF

RETURN
END FUNCTION gasdev
!------------------------------------------------------------


!------------------------------------------------------------
SUBROUTINE givetemperature(pxt,pyt,pzt,nteil,temperat,masse, ipflag)
!------------------------------------------------------------

! Produces thermal distribution of ionic velocities.

USE params
IMPLICIT NONE

REAL(DP),INTENT(OUT) :: pxt(:)
REAL(DP),INTENT(OUT) :: pyt(:)
REAL(DP),INTENT(OUT) :: pzt(:)
INTEGER,INTENT(IN) :: nteil
REAL(DP),INTENT(IN) :: temperat
REAL(DP),INTENT(IN)  :: masse
INTEGER,INTENT(IN) :: ipflag

!     temper given in Kelvin, convert to kT in Ryd now

INTEGER :: i,irand,ntmob
REAL(DP) :: ekhaben,eksoll,fac,sumx,sumy,sumz,sc,temper

WRITE(6,*) 'ENTERING GIVETEMPERATURE'

temper=temperat*6.507D-6

fac=SQRT(temper/masse)

sumx=0D0
sumy=0D0
sumz=0D0
irand=-4

DO i=1,nteil
  pxt(i)=gasdev(irand)*fac
  pyt(i)=gasdev(irand)*fac
  pzt(i)=gasdev(irand)*fac
  sumx = sumx + pxt(i)
  sumy = sumy + pyt(i)
  sumz = sumz + pzt(i)
END DO

DO i=1,nteil
  pxt(i)=pxt(i)-sumx/nteil
  pyt(i)=pyt(i)-sumy/nteil
  pzt(i)=pzt(i)-sumz/nteil
END DO


ekhaben = 0D0


ntmob=0

IF (ipflag == 4) THEN
  IF (imob /= 0) ntmob=nion
#if(raregas)
ELSE IF (ipflag == 1) THEN
  DO i=1,nc
    IF (imobc(i) /= 0) ntmob=ntmob+1
  END DO
ELSE IF (ipflag == 2) THEN
  DO i=1,NE
    IF (imobe(i) /= 0) ntmob=ntmob+1
  END DO
ELSE IF (ipflag == 3) THEN
  DO i=1,nk
    IF (imobk(i) /= 0) ntmob=ntmob+1
  END DO
#endif
ELSE
  STOP 'Error in routine giveTemperature'
END IF

eksoll = (3D0*ntmob-3D0)/2D0*temper

!      write(6,*) 'eksoll: ',eksoll

DO i=1,nteil
  IF (ipflag == 4 .AND. imob /= 0) THEN
    ekhaben=ekhaben+ pxt(i)**2+pyt(i)**2+pzt(i)**2
#if(raregas)
  ELSE IF (ipflag == 1 .AND. imobc(i) /= 0) THEN
    ekhaben=ekhaben+ pxt(i)**2+pyt(i)**2+pzt(i)**2
  ELSE IF (ipflag == 2 .AND. imobe(i) /= 0) THEN
    ekhaben=ekhaben+ pxt(i)**2+pyt(i)**2+pzt(i)**2
  ELSE IF (ipflag == 3 .AND. imobk(i) /= 0) THEN
    ekhaben=ekhaben+ pxt(i)**2+pyt(i)**2+pzt(i)**2
#endif
  END IF
END DO

ekhaben = ekhaben /2D0/masse

!      write(6,*) 'ekhaben: ',ekhaben

sc = SQRT(eksoll/ekhaben)

DO i=1,nteil
  pxt(i)=pxt(i)*sc
  pyt(i)=pyt(i)*sc
  pzt(i)=pzt(i)*sc
END DO

ekhaben = 0D0

DO i=1,nteil
  IF (ipflag == 4 .AND. imob /= 0) THEN
    ekhaben=ekhaben+ pxt(i)**2+pyt(i)**2+pzt(i)**2
#if(raregas)
  ELSE IF (ipflag == 1 .AND. imobc(i) /= 0) THEN
    ekhaben=ekhaben+ pxt(i)**2+pyt(i)**2+pzt(i)**2
  ELSE IF (ipflag == 2 .AND. imobe(i) /= 0) THEN
    ekhaben=ekhaben+ pxt(i)**2+pyt(i)**2+pzt(i)**2
  ELSE IF (ipflag == 3 .AND. imobk(i) /= 0) THEN
    ekhaben=ekhaben+ pxt(i)**2+pyt(i)**2+pzt(i)**2
#endif
  END IF
END DO

ekhaben = ekhaben /2D0/masse
!      write(6,*) 'ekhaben nachher: ',ekhaben

RETURN
END SUBROUTINE givetemperature
!------------------------------------------------------------


!-----safeopen-----------------------------------------------

SUBROUTINE safeopen(iunit,iact,icond,filename)

! Opens file unit 'iunit' and asks for status before.

USE params
IMPLICIT NONE


INTEGER, INTENT(IN)           :: iunit,iact,icond
CHARACTER (LEN=*), INTENT(IN) :: filename


LOGICAL :: topenf

IF(icond > 0)  THEN
  IF(MOD(iact,icond) == 0) THEN
    INQUIRE(iunit,OPENED=topenf)
    IF(.NOT.topenf) OPEN(iunit,POSITION='append',  &
        FILE=filename//'.'//outnam)
  END IF
END IF

RETURN
END SUBROUTINE safeopen
!     ******************************


!     ******************************
SUBROUTINE cleanfile(iunit)

! Closes file unit 'iunit' if it was open.

IMPLICIT NONE

INTEGER,INTENT(IN) :: iunit
LOGICAL :: topen

INQUIRE(iunit,OPENED=topen)
IF(topen) CLOSE(iunit)

RETURN
END SUBROUTINE cleanfile
!     ******************************


#if(parano)
!     ******************************

SUBROUTINE r_project(qin,qout,ispact,q0)

!     ******************************
!     projects all occupied states 'q0' out of 'qin'.
!      q0     = set of s.p. wavefunctions (real)
!      qin    = wavefunction from which 'q0' are to be removed (real)
!      qout   = resulting wavefunction
!      ispact = spin of 'qin'

USE params
IMPLICIT NONE 

REAL(DP), INTENT(IN)   :: q0(kdfull2,kstate)
REAL(DP), INTENT(IN)   :: qin(kdfull2)
REAL(DP), INTENT(OUT)  :: qout(kdfull2)
INTEGER, INTENT(IN)    :: ispact

INTEGER :: i, nbe
REAL(DP) :: ovl
!*********************************************************

DO  i=1,nxyz
  qout(i)=qin(i)
END DO

DO nbe=1,nstate
  IF(ispin(nbe) == ispact .AND.occup(nbe) > 0.99999999D0) THEN
    ovl=0D0
    DO i=1,nxyz
      ovl=ovl+q0(i,nbe)*qout(i)
    END DO
    ovl=ovl*dvol
    DO  i=1,nxyz
      qout(i)=qout(i)-ovl*q0(i,nbe)
    END DO
  END IF
END DO

RETURN
END SUBROUTINE r_project
!     ******************************

SUBROUTINE c_project(qin,qout,ispact,q0)

!     ******************************
!     projects all occupied states 'q0' out of 'qin'.
!      q0     = set of s.p. wavefunctions (complex)
!      qin    = wavefunction from which 'q0' are to be removed (complex)
!      qout   = resulting wavefunction
!      ispact = spin of 'qin'

USE params
IMPLICIT NONE

COMPLEX(DP), INTENT(IN)   :: q0(kdfull2,kstate)
COMPLEX(DP), INTENT(IN)   :: qin(kdfull2)
COMPLEX(DP), INTENT(OUT)  :: qout(kdfull2)
INTEGER, INTENT(IN)    :: ispact

INTEGER :: nbe
COMPLEX(DP) :: ovl

!*********************************************************

qout=qin

DO nbe=1,nstate
  IF(ispin(nbe) == ispact .AND.occup(nbe) > 0.99999999D0) THEN
    ovl=dvol*SUM(CONJG(q0(:,nbe))*qout)
    qout(:)=qout(:)-q0(:,nbe)*ovl
  END IF
END DO

RETURN
END SUBROUTINE c_project
#endif

#if(parayes)
!     ******************************

SUBROUTINE r_project(qin,qout,ispact,q0)

!     ******************************
!     projects all occupied states 'q0' out of 'qin'.
!      q0     = set of s.p. wavefunctions (real)
!      qin   = wavefunction from which 'q0' are to be removed (real)
!      ispact = spin of 'qin'

USE params
IMPLICIT NONE

REAL(DP), INTENT(IN)   :: q0(kdfull2,kstate)
REAL(DP), INTENT(IN)   :: qin(kdfull2)
REAL(DP), INTENT(OUT)  :: qout(kdfull2)
INTEGER, INTENT(IN)    :: ispact

!*********************************************************
qout=0D0  !to avoid warning at compile time
STOP ' subroutine PROJECT not yet for parallel'

RETURN
END SUBROUTINE r_project
!     ******************************

SUBROUTINE c_project(qin,qout,ispact,q0)                                   ! cPW

!     ******************************
!     projects all occupied states 'q0' out of 'qin'.
!      q0     = set of s.p. wavefunctions (real)
!      qin    = wavefunction from which 'q0' are to be removed (real)
!      qout   = resulting wavefunction
!      ispact = spin of 'qin'

USE params
IMPLICIT NONE

COMPLEX(DP), INTENT(IN)   :: q0(kdfull2,kstate)
COMPLEX(DP), INTENT(IN)   :: qin(kdfull2)
COMPLEX(DP), INTENT(OUT)  :: qout(kdfull2)
INTEGER, INTENT(IN)    :: ispact

COMPLEX(DP) :: ovl

!*********************************************************
qout=0D0  !to avoid warning at compile time
STOP ' subroutine CPROJECT not yet for parallel'

RETURN
END SUBROUTINE c_project
#endif




!-----pair  ------------ part of the Hartree-Fock package --------------

SUBROUTINE pair(e,gw,ph,nz,nmax,gp,eferm,delta,partnm,  &
    iter,ipair,eps,iab,kstate)
USE params, ONLY: DP,one,zero
IMPLICIT NONE

!     * * * * * * * * * *     pairing iteration      * * * * * * * * * *

!     determines the Fermi-energy 'eferm' and the gap 'delta'
!     such that the "gap-equation" and the "particle number condition"
!     are fulfilled simultaneously. this is achieved by iteration
!     with a gradient scheme (newton's tangential formula)
!     where the gradient is evaluated by finite differences.

!     there are two schemes possible:
!     1. the "constant force" approach, where the pairing force is given
!                                       in 'gp'.
!     2. the "constant gap" approach,   where the gap is kept fixed and
!                                       is given in 'delta'.
!     there is only 'eferm' to be determined in the "constant gap"
!     approach. the gap-equation is used at the end to determine an
!     according force 'gp'.

!     input via list:
!      e       = array(1:nmax) of single particle energies.
!      gw      = array(1:nmax) of occupation probabilities,
!                input is used as initial guess for the iteration.
!      ph      = array(1:nmax) of the multiplicities of the single
!                particle levels.
!      nz      = desired particle number
!      nmax    = maximum number of single particle states
!                (needs to be smaller than the parameter 'kstate'
!                 and smaller than 46, for dimensioning reasons).
!      gp      = size of pairing-force, is input in case of
!                "constant force" approach.
!      eferm   = Fermi-energy, input is used as initial guess;
!                if eferm=0.0 is input then a suitable initial guess
!                is computed from input 'gw'.
!      delta   = pairing gap, is input in case of "constant gap" approach
!      iter    = number of pairing iterations,
!                abs(iter) is the number of iterations
!      ipair   = switch for pairing type
!                 ipair=1 -->  "constant force" approach, using
!                              pairing force 'gp'
!                 ipair=0 -->  "constant gap" approach, using
!                              given gap 'delta'
!                 ipair=4 -->  thermal occupation with temperature
!                              given by 'delta'
!      eps     = end condition for pairing iteration;
!                iteration terminates if error in particle number
!                is smaller than 'eps' and if 'eferm' and/or 'delta'
!                are smaller than 'eps'.
!      iab     = print level; significant values are -2, -1, 0 and 1;
!                larger value prints more.

!     output via list:
!      gw      = array(1..nmax) of the occupation probabilities
!                finally achieved.
!      gp      = size of pairing-force, is output in case of
!                "constant gap" approach.
!      delta   = pairing gap, is output in case of
!                "constant force" approach.
!      eferm   = Fermi-energy.
!      sum     = particle number, from adding 'gw' with degeneracy 'ph'.

!     some auxiliary fields are:
!      ph      = array(1:46) of the multiplicities of the single
!                particle levels. array is ordered according to
!                the harmonic oscillator scheme, with rising shell
!                number and lowering angular momentum within each shell.
!      dparn   = array(1:3) of variation in the particle number as
!                consequence of 'deferm' and 'ddelta'.
!      dgap    = array(1:3) of variation in the gap equation as
!                consequence of 'deferm' and 'ddelta'.
!      gr      = double array of variation of particle number and of
!                gap equation, used for construction of the tangential
!                step

REAL(DP),INTENT(IN) :: e(kstate)
REAL(DP),INTENT(IN OUT) :: gw(kstate)
REAL(DP),INTENT(IN OUT) :: ph(kstate)
INTEGER,INTENT(IN) :: nz
INTEGER,INTENT(IN) :: nmax
REAL(DP),INTENT(IN OUT) :: gp
REAL(DP),INTENT(IN OUT) :: eferm
REAL(DP),INTENT(IN OUT) :: delta
REAL(DP),INTENT(IN OUT) :: partnm
INTEGER,INTENT(IN) :: iter
INTEGER,INTENT(IN) :: ipair
REAL(DP),INTENT(IN) ::eps
INTEGER,INTENT(IN) :: iab
INTEGER,INTENT(IN) :: kstate
    
LOGICAL :: converged= .false.   ! will be true when convergence is reached
INTEGER :: i, isum, it, itmaxp, itsum, k, ka
REAL(DP):: dampin, deferm, ddelta, delti, delact, delold, del2
REAL(DP):: elam, equasi, efold, efstep, parnum, gapeq, gphalf, grsav
REAL(DP):: emax, emin
REAl(Dp):: eferm1, eferm2, eferm3, parn1, parn2, parn3
REAL(DP):: efmax, dlmax, det, detmin, defmin, ddlmin, dlstep
REAL(DP):: xmaxlw,x10lw
REAL(DP):: dparn(3), dgap(3),gr(2,2)


DATA gr,dparn,dgap/10*0D0/

!     further numerical variables are:
!      deferm  = step in the Fermi-energy, to evaluate
!                numerically the derivatives.
!      ddelta  = step in the gap, to evaluate numerically the
!                derivatives; used only in the "constant force" approach
!      dampin  = damping factor for the tangential step
!      efmax,dlmax,detmin,defmin,ddlmin  = limitations for the step

!     note: these parameters are adjusted for an application with
!           the nuclear shell model. it is assumed that the single
!           particle energies are given in units of [mev] and display
!           a reasonable nuclear spectrum.

DATA dampin,deferm,ddelta/0.9D0,1D0,0.1D0/
DATA efmax,dlmax,detmin,defmin,ddlmin/3D0,.6D0,.1D-3,.1D-3,.1D-3/
!clust      data dampin,deferm,ddelta/0.9,0.01,0.1/
!clust      data efmax,dlmax,detmin,defmin,ddlmin/0.03,.6,.1e-3,.1e-3,.1

DATA xmaxlw,x10lw/1.0D-30,1.0D-10/

!-----------------------------------------------------------------------

IF(iab > 1) THEN
  
  WRITE(7,'(/a,i4,3(a,f9.4)/a,10(/1x,10f8.3))')  &
      ' *** pair entered with: iter=',iter,'  delta=',delta,  &
      '  eferm=',eferm,'  gp=',gp,  &
      '                        single particle energies:', (e(i),i=1,nmax)
  WRITE(6,'(/a,i4,3(a,f9.4)/a,10(/1x,10f8.3))')  &
      ' *** pair entered with: iter=',iter,'  delta=',delta,  &
      '  eferm=',eferm,'  gp=',gp,  &
      '                        single particle energies:', (e(i),i=1,nmax)
  WRITE(6,*) 'kstate',kstate
  WRITE(7,*) 'kstate',kstate
  WRITE(6,*) 'nmax',nmax
  WRITE(7,*) 'nmax',nmax
  WRITE(6,'(a,500f6.2)') 'gw',(gw(i),i=1,nmax)
  WRITE(6,'(a,500f6.2)') 'ph',(ph(i),i=1,nmax)
  WRITE(6,*) 'nz',nz
  WRITE(7,*) 'gw',(gw(i),i=1,nmax)
  WRITE(7,*) 'ph',(ph(i),i=1,nmax)
  WRITE(7,*) 'nz',nz
  WRITE(6,*) partnm,ipair,eps,iab
  WRITE(7,*) partnm,ipair,eps,iab
END IF

!      prepare Fermi-energy 'eferm', if not yet initialised.

IF(eferm == zero) THEN
  isum   = 0
  DO it=1,nmax
    itsum  = it
    isum   = INT(ph(it)+.000001D0) +isum
    IF(isum+isum-nz >= 0) EXIT
  END DO
  IF(it > nmax) STOP ' not enough states to reach particle number'
  eferm  = e(itsum)+deferm
  IF(iab >= 1) WRITE(7,'(1x,a,i5,(a,g11.4),a,i5)')  &
      'initialistion: isum=',isum,'  eferm=',eferm,'  nz=',nz
  
END IF

IF(ipair == 1) THEN
  
!-----------------------------------------------------------------------
  
!     case of constant pairing-force, given by 'gp'.
!      gap 'delta', Fermi-energy 'eferm', and occupations 'gw'
!      are computed.
  
!-----------------------------------------------------------------------
  
  delti  = delta
  IF(ABS(delti) < x10lw) delti = one
  DO it=1,ABS(iter)
    
!       gap-equation & particle-number at 'eferm' & 'delti' (for k=1)
!       and at 'eferm' varied (k=2) and at 'delta' varied (k=3).
    
    DO ka=1,3
!                                          actual Fermi-energy and gap
      k      = 4-ka
      IF(k == 1) THEN
        elam   = eferm
        delact = delti
      ELSE IF(k == 2) THEN
        elam   = eferm+deferm
        delact = delti
      ELSE
        elam   = eferm
        delact = delti+ddelta
      END IF
!                                         compute actual quasiparticle e
!                                         and occupation-weights,
!                                         accumulate gap-equation and
!                                         particle number
      parnum = zero
      gapeq  = zero
      DO i=1,nmax
        equasi = SQRT((e(i)-elam)*(e(i)-elam)+delact*delact)
        gw(i)  = 0.5D0-0.5D0*(e(i)-elam)/equasi
        gapeq  = ph(i)/equasi   + gapeq
        parnum = gw(i)*ph(i)    + parnum
      END DO
!                                         store actual part.number and
!                                         gap-equation
      dgap(k)  = 0.25D0*gp*gapeq-one          !?? extra 1/2 because of w
      dparn(k) = parnum-nz
    END DO
    
!       construction of the gradient-matrix and its inverse
!        ( determinant 'det' is limited in value to avoid singularities
    
    gr(1,1) = ( dgap(2)- dgap(1))/deferm
    gr(1,2) = (dparn(2)-dparn(1))/deferm
    gr(2,1) = ( dgap(3)- dgap(1))/ddelta
    gr(2,2) = (dparn(3)-dparn(1))/ddelta
    det     = gr(1,1)*gr(2,2)-gr(1,2)*gr(2,1)
    det     = SIGN(MAX(ABS(det),detmin),det)
    IF(ABS(det) < xmaxlw)  det    = SIGN(xmaxlw,det)
    
    IF(iab >= 1) THEN
      WRITE(7,'(a,i4,3(a,g11.3)/a,3g11.3,a,3g11.3/a,4g11.3)')  &
          ' iteration at it=',it,' :  delti=',delti,'  eferm=',eferm,  &
          '  det=',det, ' dgap=',dgap,' ,  dparn=',dparn,  &
          '   gr=',gr
      WRITE(7,'(a,10f7.4,10(/6x,10f7.4))') '   gw=',gw
    END IF
    
    grsav   =  gr(1,1)
    gr(1,1) =  gr(2,2)/det
    gr(1,2) = -gr(1,2)/det
    gr(2,1) = -gr(2,1)/det
    gr(2,2) = grsav/det
    
!       tangential step  (including some step limiters)
    
    efold  = eferm
    delold = delti
    efstep = dampin*(gr(1,1)*dgap(1)+gr(2,1)*dparn(1))
    dlstep = dampin*(gr(1,2)*dgap(1)+gr(2,2)*dparn(1))
    eferm  = eferm-SIGN(MIN(ABS(efstep),efmax),efstep)
    delti  = ABS(delti-SIGN(MIN(ABS(dlstep),dlmax),dlstep) )
    deferm = efold-eferm
    ddelta = delold-delti
    IF(ABS(deferm) < defmin) deferm = defmin
    IF(ABS(ddelta) < ddlmin) ddelta = ddlmin
    
    IF(iab == 0) WRITE(7,'(a,i4,2(a,g11.3))')  &
        ' pairing iteration it=',it,'  :  delti=',delti, '  eferm=',eferm
    
    IF(ABS(dgap(1)) < eps .AND. ABS(dparn(1)) < eps) THEN
      converged=.true.
      EXIT
    ENDIF
  END DO

  delta  = delti
  
ELSE

!-----------------------------------------------------------------------
  
!      case of constant paring-gap, given by 'delta'.
!      using now the safer secant step.
  
!-----------------------------------------------------------------------
  
!      determine upper and lower bound for secant step
  
  emin = +1.0D30
  emax = -1.0D30
  DO i=1,nmax
    emin = MIN(e(i),emin)
    emax = MAX(e(i),emax)
  END DO
  itmaxp = ABS(iter)
  IF(itmaxp < 0) THEN
    eferm1 = eferm
    eferm2 = eferm
    eferm3 = eferm
    parn1  = parnm(e,gw,ph,nmax,delta,eferm,ipair)
    parn2  = parn1
    parn3  = parn1
    IF(parn1-nz > zero) THEN
      deferm = (eferm-emin)/itmaxp
      IF(iab >= 1) WRITE(7,'(a,3(1pg12.4))')  &
          ' search brackets: parn2,eferm,emin=',parn2,eferm,emin
      DO it=1,itmaxp
        eferm1 = eferm1-deferm
        parn1  = parnm(e,gw,ph,nmax,delta,eferm1,ipair)
        IF(iab >= 1) WRITE(7,'(a,i5,3(1pg12.4))')  &
            ' search brackets: it,parn1,eferm1=',it,parn1,eferm1
        IF(parn1-nz < zero) EXIT
        eferm2 = eferm1
      END DO
      IF(it>itmaxp) STOP ' particle number not embraced in pair'
    ELSE IF(parn1-nz < zero) THEN
      deferm = (emax-eferm)/itmaxp
      IF(iab >= 1) WRITE(7,'(a,3(1pg12.4))')  &
          ' search brackets: parn1,eferm,emax=',parn1,eferm,emax
      DO it=1,itmaxp
        eferm2 = eferm2+deferm
        parn2  = parnm(e,gw,ph,nmax,delta,eferm2,ipair)
        IF(iab >= 1) WRITE(7,'(a,i5,3(1pg12.4))')  &
            ' search brackets: it,parn2,eferm2=',it,parn2,eferm2
        IF(parn2-nz > zero) EXIT
        eferm1 = eferm2
      END DO
      IF(it>itmaxp) STOP ' particle number not embraced in pair'
    ELSE
      IF(iab >= 1) WRITE(7,'(a)') ' solution immediately found '
      converged=.true.   ! jump to end
    END IF


    IF(.NOT. converged) THEN
!   Perform the secant search
      IF(iab >= 1) WRITE(7,'(a)') ' brackets found '
     
      DO it=1,itmaxp
        eferm3 = eferm1 + (nz-parn1)*(eferm2-eferm1)/(parn2-parn1)
        parn3  = parnm(e,gw,ph,nmax,delta,eferm3,ipair)
        IF(iab >= 0) WRITE(7,'(a,i4,2(a,g11.3))')  &
            ' secant iteration it=',it,'  :  eferm=',eferm3, '  dparn(1)=',parn3-nz
        IF(ABS(parn3-nz) < eps) THEN
          converged=.true.
          EXIT
        ELSE IF(parn3-nz > zero) THEN
          eferm2 = eferm3
          parn2  = parn3
        ELSE
          eferm1 = eferm3
          parn1  = parn3
        END IF
      END DO
    END IF
  END IF
  

  IF(.NOT. converged) THEN
!      the secant step did not succeed. use bisection
    eferm1 = emin-4*delta
    eferm2 = emax+4*delta
    DO it=1,itmaxp
      eferm3 = 0.5D0*(eferm2+eferm1)
      parn3  = parnm(e,gw,ph,nmax,delta,eferm3,ipair)
      IF(iab >= 0) WRITE(7,'(a,i4,2(a,g11.3))')  &
          ' bisection it=',it,'  :  eferm=',eferm3, '  dparn(1)=',parn3-nz
      IF(ABS(parn3-nz) < eps) THEN
        converged=.true.
        EXIT
      ELSE IF(parn3-nz < zero) THEN
        eferm1 = eferm3
        parn1  = parn3
      ELSE IF(parn3-nz > zero) THEN
        eferm2 = eferm3
        parn2  = parn3
      END IF
    END DO
  END IF
  
  eferm = eferm3
  dparn(1) = parn3 - nz
  dgap(1)  = zero
  partnm   = parn3
  
!      determine force 'gp' according to given 'delta' and the
!      equilibrium reached.
  
  IF(ipair == 0) THEN
    gapeq  = zero
    DO i=1,nmax
      gapeq = ph(i)/SQRT((e(i)-eferm)*(e(i)-eferm)+delta*delta) + gapeq
    END DO
    gp     = 4D0/gapeq                       !?? account for wate=2.
  END IF
  
!     end of big switch between cases
  
END IF

!-----------------------------------------------------------------------

!     finally the gap-equation and particle-number are checked and print

!-----------------------------------------------------------------------


!     in case of incomplete convergence print a warning

IF(iab >= -1 .AND. .NOT. converged) WRITE(7,'(2(a,g11.3),2(a,i3))')  &
    ' ---> no convergence in pair! err(gap-eq)=',dgap(1),  &
    '  err(part.nr)=',dparn(1),'  nz=',nz,'  nmzx=',nmax

!     check gap equation

IF(iab >= 0) THEN
  IF(ipair /= 4) THEN
    gphalf = gp*0.5D0
    del2   = MAX(xmaxlw,delta*delta)
    gapeq  = zero
    partnm = zero
    DO i=1,nmax
      partnm = gw(i)*ph(i) + partnm
      gapeq  = gphalf*ph(i)/SQRT(del2+(e(i)-eferm)**2) + gapeq
    END DO
  ELSE
    gapeq = zero
  END IF
  
  WRITE(7,'(a,i4,2(a,g13.5))')  &
      ' pair finished with ',it,' iterations: gap-eq.=',gapeq, '  part.nr=',partnm
END IF

RETURN
END SUBROUTINE pair
!-----parnm------------ part of the Hartree-Fock package ---------------

REAL(DP) FUNCTION parnm(e,gw,ph,nmax,delta,elam,ipair)

!     the include file 'params.inc' communicates the precision
!     transfering a line 'implicit double precision (a-h,o-z)' or not.
!     furthermore it carries the common block /constn/ from which
!     the pseudo-constant one=1.0 is taken here.


!     the function computes the particle number in a pairing or
!     temperature distribution. the parameters are:
!       e      = array of s.p. energies
!       gw     = array of occupation weights (is workspace and output!)
!       ph     = array of degeneracy factors
!       nmax   = number of states in the above three arrays
!       delta  = pairing gap or temperature
!       elam   = Fermi energy
!       ipair  = switch to pairing (0) or temperature (4).

USE params
IMPLICIT NONE

REAL(DP),INTENT(IN) :: e(ksttot)
REAL(DP),INTENT(OUT) :: gw(ksttot)
REAL(DP),INTENT(IN) :: ph(ksttot)
INTEGER,INTENT(IN) :: nmax
REAL(DP),INTENT(IN) :: delta
REAL(DP),INTENT(IN) :: elam
INTEGER,INTENT(IN) :: ipair

INTEGER :: i
REAL(DP) :: acc, equasi

!----------------------------------------------------------------------

acc    = zero
IF(ipair == 0) THEN
  DO i=1,nmax
    equasi = SQRT( (e(i)-elam)*(e(i)-elam) + delta*delta )
    gw(i)  = 0.5D0-0.5D0*(e(i)-elam)/equasi
    acc    = gw(i)*ph(i)    + acc
  END DO
ELSE IF(ipair == 4) THEN            ! case of temperature
  DO i=1,nmax
    equasi = (e(i)-elam)/delta
    IF(equasi > one*50D0) THEN
      gw(i) = 0D0
    ELSE IF(equasi < -one*50D0) THEN
      gw(i)  = 1D0
    ELSE
      gw(i)  = 1D0/(1D0+EXP(equasi))
    END IF
    acc    = gw(i)*ph(i)    + acc
  END DO
END IF

parnm = acc

RETURN
END FUNCTION parnm




!-----omega_mieplasmon--------------------------------------------------

REAL(DP) FUNCTION omega_mieplasmon(rho)

!     Estimates the Mie plasmon frequency
!     (? seems to be correct only for Na with soft local PsP ?)

USE params
IMPLICIT NONE

REAL(DP), INTENT(IN) :: rho(kdfull2)


INTEGER :: i,ii, ix, iy, iz
REAL(DP) :: psrho1, psrho2, acc, omegam
REAL(DP) :: rx, ry, rz, r2, x1, y1, z1

REAL(DP),DIMENSION(:),ALLOCATABLE :: rhops

REAL(DP),PARAMETER:: prho1_data=-0.46073D0, prho2_data=0.13287D0



!------------------------------------------------------------

ALLOCATE(rhops(kdfull2))

IF(nion2 == 2) RETURN

IF(ipsptyp == 0) THEN
  
!     total electronic charge from given density 'rho'
  
  apnum=0D0
  DO i=1,nxyz
    apnum=apnum+rho(i)
  END DO
  apnum=apnum*dvol
  
!     positive background charge distribution for different cases:
!     jellium or local Gaussian pseudo-densities
  
  IF(nion2 == 0) THEN
    DO i=1,nxyz
      rhops(i)=rhojel(i)
    END DO
  ELSE
    DO i=1,nxyz
      rhops(i)=0D0
    END DO
    
    DO ii=1,nion
      i=0
      DO iz=minz,maxz
        z1=(iz-nzsh)*dz
        rz=z1-cz(ii)
        DO iy=miny,maxy
          y1=(iy-nysh)*dy
          ry=y1-cy(ii)
          DO ix=minx,maxx
            x1=(ix-nxsh)*dx
            i=i+1
            rx=x1-cx(ii)
            r2=rx*rx+ry*ry+rz*rz
            psrho1=prho1_data*EXP(-r2/(2.0D0*sgm1(11)*sgm1(11)))
            psrho2=prho2_data*EXP(-r2/(2.0D0*sgm2(11)*sgm2(11)))
            rhops(i)=rhops(i)+psrho1+psrho2
          END DO
        END DO
      END DO
    END DO
  END IF
  
  rhomix=0D0
  acc = 0D0
  DO i=1,nxyz
    acc = acc+rhops(i)
    rhomix=rhomix+rho(i)*rhops(i)
  END DO
  rhopss = acc*dvol
  rhomix=rhomix*dvol
  
!     16*pi = 4*pi * e^2 * hbar^2/4m
  
  omegam=SQRT(16.0D0*pi/3.0D0*rhomix/apnum)
ELSE
  omegam = 0D0                 ! Mie plasmon not installed for Goedecker
  apnum = 0D0
  acc = 0D0
  rhomix = 0D0
END IF

DEALLOCATE(rhops)

omega_mieplasmon = omegam

RETURN
END FUNCTION omega_mieplasmon
!-----view3D------------------------------------------------------

SUBROUTINE view3d()

!     Prepares file for viewing ionic structure with 'xbs'.

USE params
IMPLICIT NONE
CHARACTER (LEN=3) :: ext
INTEGER :: ion
!-----------------------------------------------------------------

ext=outnam

OPEN(39, FILE=ext//'.bs', STATUS='UNKNOWN')
DO ion=1,nion
  WRITE(39,'(a4,a6,a1,a6,3f10.4)') 'atom',' ','H',' ', cx(ion),cy(ion),cz(ion)
END DO
WRITE(39,*) ' '
WRITE(39,'(a4,a6,a1,2(a6,a3))') 'spec',' ','H',' ', '0.5',' ','1.0'
WRITE(39,*) ' '
WRITE(39,'(a5,a5,a1,a6,a1,4(a6,a3))') 'bonds',' ','H',' ',  &
    'H',' ','0.0',' ','9.0',' ','0.1',' ','1.0'
WRITE(39,*) ' '
WRITE(39,'(a4,9(a2,a3))') 'tmat',' ','1.0',' ','0.0',' ','0.0',  &
    ' ','0.0',' ','1.0',' ','0.0',' ','0.0',' ','0.0',' ','1.0'
WRITE(39,'(a4,a6,a4)') 'dist',' ','12.0'
WRITE(39,'(a3,a8,a3)') 'inc',' ','5.0'
WRITE(39,'(a5,a5,a4)') 'scale',' ','20.0'
WRITE(39,'(a4,a7,a3)') 'rfac',' ','1.0'
WRITE(39,'(a4,a7,a3)') 'bfac',' ','1.0'
WRITE(39,'(a3,2(a6,a3))') 'pos',' ','0.0',' ','0.0'
WRITE(39,'(a8,9(a2,a1))') 'switches',' ','1',' ','0',' ','1',' ','0',  &
    ' ','0',' ','1',' ','1',' ','0',' ','0'
CLOSE(39)

RETURN
END SUBROUTINE view3d



!-----timer------------------------------------------------------

SUBROUTINE timer(iaction)

!     Initializes timer and retrieves CPU times
!       iaction  = 1  for initial call
!                  2  for subsequent call
!     The relative time is taken between step 2 and 1.

USE params
IMPLICIT NONE

INTEGER, INTENT(IN) :: iaction

INTEGER :: i 
REAL(DP),ALLOCATABLE,SAVE :: cputimenew(:),cputimeold(:)

LOGICAL,SAVE :: tfirst=.true.

!-----------------------------------------------------------------

IF(tfirst) THEN
  ALLOCATE(cputimenew(0:knode),cputimeold(0:knode))
!  WRITE(*,*) ' cputime allocated. knode=',knode
!  WRITE(*,*) cputimeold
  tfirst=.FALSE.
ELSE
!  WRITE(*,*) ' next call to timer. cputimeold=',cputimeold
END IF

IF(iaction == 1) THEN
  CALL cpu_time(cputimeold(0))
ELSE IF(iaction == 2) THEN
  CALL cpu_time(cputimenew(0))
  IF(myn == 0) WRITE(6,'(a,1pg20.10)') 'time for step in node 0=',  &
      cputimenew(0)-cputimeold(0)
  DO i=0,knode
    cputimeold(i) = cputimenew(i)
  END DO
ELSE
  STOP 'TIMER: this IACTION is not valid'
END IF

RETURN
END SUBROUTINE timer


!-----stimer------------------------------------------------------

SUBROUTINE stimer(iaction)

!     Initializes timer and retrieves CPU times
!       iaction  = 1  for initial call
!                  2  for subsequent call
!     The relative time is taken between step 2 and 1.

USE params
IMPLICIT NONE

INTEGER, INTENT(IN) :: iaction

INTEGER, SAVE :: itimeold=0,itimenew=0

!-----------------------------------------------------------------


IF(iaction == 1) THEN
  CALL system_clock(itimeold)
ELSE IF(iaction == 2) THEN
  CALL system_clock(itimenew)
  IF(myn == 0) WRITE(6,'(a,1pg13.5)') 'sys.time for step in node 0=',  &
      (itimenew-itimeold)*1D-4
  itimeold = itimenew
ELSE
  STOP 'SYS.TIMER: this IACTION is not valid'
END IF

RETURN
END SUBROUTINE stimer



!---probab----------------------------------------------------------

SUBROUTINE probab(psitmp)

! Calculates the probabilities to find a certain charge stae from
! the distribution of loss of norms of the s.p. wavefunctions
! given set of wavefunctions in input array 'psitmp'.

USE params
IMPLICIT NONE
COMPLEX(DP), INTENT(IN)                :: psitmp(kdfull2,kstate)

#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
INTEGER :: ic, k, n, nact, nb, nod2
REAL(DP)                               :: prob(ksttot+1),pold(ksttot+1)
REAL(DP)                               :: rtmpuse(nstate)
REAL(DP)                               :: rtmpuse_all(ksttot)
#else
INTEGER :: k
REAL(DP)                               :: prob(nstate+1),pold(nstate+1)
REAL(DP)                               :: rtmpuse(nstate)
#endif

INTEGER :: i,it
COMPLEX(DP)                            :: cscal
!-----------------------------------------------------------------

DO i=1,nstate
   cscal=orbitaloverlap(psitmp(:,i),psitmp(:,i))
   rtmpuse(i)=SQRT(REAL(cscal,DP)**2 + AIMAG(cscal)**2)
   prob(i)=0D0
ENDDO

#if(parayes)
CALL  mpi_comm_rank(mpi_comm_world,myn,mpi_ierror)
IF(myn/=0) THEN
  prob(ksttot+1)=0D0
  CALL mpi_send(rtmpuse,nstate,mpi_double_precision,0,myn,mpi_comm_world,mpi_ierror)
  DO n=1,nstate
    nact = nrel2abs_other(n,0) 
!    rtmpuse_all(nact) = rtmuse(n)
  ENDDO
ELSE
  DO nod2=0,knode-1
    IF(nod2 > 0) CALL mpi_recv(rtmpuse,nstate,mpi_double_precision,  &
          nod2,mpi_any_tag,mpi_comm_world,is,ic)
    DO nb=1,nstate_node(nod2)
      nact = nrel2abs_other(nb,nod2) 
      rtmpuse_all(nact) = rtmpuse(nb)
    END DO
  END DO
ENDIF
#else
prob(nstate+1)=0D0
#endif

#if(parayes)
IF(myn==0) THEN
   prob(1)=rtmpuse_all(1)
   prob(2)=1D0-rtmpuse_all(1)
   DO i=2,ksttot
      DO k=1,i
      	 pold(k)=prob(k)
      ENDDO
   
      DO k=1,i
         prob(k)=pold(k)*rtmpuse_all(i)
      ENDDO
   
      DO k=2,i+1
         prob(k)=prob(k)+pold(k-1)*(1D0-rtmpuse_all(i))
      ENDDO
   ENDDO

   CALL safeopen(808,it,jnorms,'pproba')
   WRITE(808,'(500f12.8)') tfs,prob(1:(ksttot+1))    
   CALL flush(808)
ENDIF
#else
prob(1)=rtmpuse(1)
prob(2)=1D0-rtmpuse(1)
           
DO i=2,nstate
   DO k=1,i
      pold(k)=prob(k)
   ENDDO
   
   DO k=1,i
      prob(k)=pold(k)*rtmpuse(i)
   ENDDO
   
   DO k=2,i+1
      prob(k)=prob(k)+pold(k-1)*(1D0-rtmpuse(i))
   ENDDO
ENDDO

CALL safeopen(808,it,jnorms,'pproba')
WRITE(808,'(500f12.8)') tfs,prob(1:(nstate+1))      
CALL flush(808)
#endif

RETURN
END SUBROUTINE probab

!-----------------------------------------------------------------
 
!-----phstate---------------------------------------------------------

SUBROUTINE phstate(psi)

!  Picks one particular 1ph state from the given wavefunctions 
!  'psi' and mixes this to a rotated ph configuration. The state
!  is rotated by 'phangle'.
!  The particle state is 'npstate' and the hole state 'nhstate'.
!  The old hole wavefunction is stored and returned on 'oldhole'.
!  The parameters 'phangle', 'npstate' and 'nhstate' are communicated
!  through module PARAMS.
!  'oldhole' is also communicated through PARAMS and allocated here.


USE params
IMPLICIT NONE

COMPLEX(DP), INTENT(IN OUT)   :: psi(kdfull2,kstate)

REAL(DP) :: phangpi, ca, sa, cb, sb

  IF(npstate > nstate) STOP ' PHSTATE: particle state out of range'
  IF(nhstate > nstate) STOP ' PHSTATE: hole state out of range'
  IF(occup(npstate) > 0.5D0) STOP 'PHSTATE: particle state already occupied'
  IF(occup(nhstate) < 0.5D0) STOP 'PHSTATE: hole state not occupied'

  ALLOCATE(oldhole(kdfull2))
  ALLOCATE(newhole(kdfull2))

  oldhole = psi(:,nhstate)
  phangpi = phangle/180D0*PI
  ca = cos(phangpi)
  sa = sin(phangpi)
  cb = cos(phphase/180D0*PI)
  sb = sin(phphase/180D0*PI)
  
  psi(:,nhstate) = (ca*cb+CMPLX(0D0,sa*sb,DP))*psi(:,nhstate) &
                  +(sa*cb-CMPLX(0D0,ca*sb,DP))*psi(:,npstate)
  psi(:,npstate) = -(sa*cb+CMPLX(0D0,ca*sb,DP))*oldhole &
                   +(ca*cb-CMPLX(0D0,sa*sb,DP))*psi(:,npstate)

  newhole = psi(:,nhstate)


END SUBROUTINE phstate


 
!-----phoverl---------------------------------------------------------

SUBROUTINE phoverl(psi)

!  Computes overlap of original Slater state with dynamically rotated
!  1ph state. This routine makes sense only after having called
!  'phstate'.


USE params
IMPLICIT NONE

COMPLEX(DP), INTENT(IN)   :: psi(kdfull2,kstate)
LOGICAL,SAVE :: tfirst=.true.
COMPLEX(DP) :: covo,covn

IF(TFIRST) THEN
  OPEN(381,file='phoverlap.'//outnam)
  WRITE(381,'(a/a,2f10.2,2i5/a)') &
  '# overlap of ph-rotated state with orinigal hole state:', &
  '# phangle,phphase,npstate,nhstate=',phangle,phphase,npstate,nhstate, &
  '#  time (fs)   real(ov)   imag(ov)   abs(ov)    phase(ov)'
  tfirst = .false.
END IF

covo = SUM(CONJG(oldhole)*psi(:,nhstate))*dvol
covn = SUM(CONJG(newhole)*psi(:,nhstate))*dvol
WRITE(381,'(f10.3,8(1pg13.5))') tfs, &
  covo,SQRT(REAL(covo,DP)**2+AIMAG(covo)**2),atan2(AIMAG(covo),REAl(covo,DP)), &
  covn,SQRT(REAL(covn,DP)**2+AIMAG(covn)**2),atan2(AIMAG(covn),REAl(covn,DP))
CALL FLUSH(381)



RETURN

END SUBROUTINE phoverl

 
!-----stateoverl---------------------------------------------------------

SUBROUTINE stateoverl(psi1,psi2)

!  Computes overlap of two Slater states with possibly
!  different sets of s.p. states.
!  Both states must have the same sequence of occupation numbers.
!  The set of s.p. states is given by 'psi1' and 'psi2'.
!  The occupation numbers 'occup' are communicated via 'params'.
!  All states should have occupation 1.

USE params
IMPLICIT NONE

COMPLEX(DP), INTENT(IN)   :: psi1(kdfull2,kstate)
COMPLEX(DP), INTENT(IN)   :: psi2(kdfull2,kstate)

INTEGER :: nbe,nb2
COMPLEX(DP), ALLOCATABLE :: spoverlaps(:,:)

LOGICAL,SAVE :: tfirst=.true.
LOGICAL,PARAMETER :: ttest=.false.

IF(TFIRST) THEN
  OPEN(382,file='stateoverlap.'//outnam)
  WRITE(382,'(a)') &
  '# overlap of static state with dynamically propagated one:', &
  '#  time (fs)   overlap'
  DO nbe=1,nstate
    IF(occup(nbe) < 0.9999999999999D0) &
      STOP ' STATEOVERL requires full occupation'
  END DO
  tfirst = .false.
END IF

ALLOCATE(spoverlaps(nstate,nstate))

DO nbe=1,nstate
  DO nb2=1,nstate
    IF(ispin(nrel2abs(nbe)) == ispin(nrel2abs(nb2))) THEN
      spoverlaps(nb2,nbe) = SUM(CONJG(psi1(:,nbe))*psi2(:,nb2))*dvol
    ELSE
      spoverlaps(nb2,nbe) = 0D0
    END IF
  END DO
END DO


WRITE(382,'(f10.4,2(1pg13.5))') time,determinant(spoverlaps,nstate,nstate)
IF(ttest) THEN
  WRITE(382,'(a)') ' matrix up:'
  DO nbe=1,nstate/2
    WRITE(382,'(16(1pg13.5))') (spoverlaps(nb2,nbe),nb2=1,nstate/2)
  END DO
  WRITE(382,'(a)') ' matrix down:'
  DO nbe=nstate/2+1,nstate
    WRITE(382,'(16(1pg13.5))') (spoverlaps(nb2,nbe),nb2=nstate/2+1,nstate)
  END DO
END IF
CALL FLUSH(382)

DEALLOCATE(spoverlaps)

RETURN

END SUBROUTINE stateoverl


!-----projmoms-------------------------------------------------------projmoms
! REAL version
!-----------------------------------------------------------------------
SUBROUTINE r_projmoms(rho,psi)

! Multipole moments relative to c.m. on 'qetarget' and relative to 
! projectile coordinate 'cz' on 'qeproj'.

USE params
IMPLICIT NONE

REAL(DP), INTENT(IN) :: rho(2*kdfull2)
REAL(DP), INTENT(IN) :: psi(kdfull2,kstate)
#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

INTEGER :: ind, ix, iy, iz, k
REAL(DP) :: sproj,starget
REAL(DP) :: x1, y1, z1, x1t, y1t, z1t, x1p, y1p, z1p 

#if(parano)
INTEGER :: ik, ikk
#else
INTEGER :: kk,nbe,nbee
REAL(DP) :: sprojec
#endif
!----------------------------------------------------------------------
nrmom=35
IF(nrmom > kmom) STOP ' too many moments in projmoms'

DO k=1,nrmom
  qetarget(k)=0D0
  qeproj(k)=0D0
END DO

!     switch for calculating moments relative to center of mass (1)
!     or center of box (0)
rvectmp = 0D0
IF(iemomsrel == 1 .AND. nion2 > 0) CALL getcm(1,0,0)


ind=0
DO iz=minz,maxz
  z1=(iz-nzsh)*dz
  z1t=z1-rvectmp(3)
  z1p=z1-cz(nproj)
  DO iy=miny,maxy
    y1=(iy-nysh)*dy
    y1t=y1-rvectmp(2)
    y1p=y1-cy(nproj)
    DO ix=minx,maxx
      ind=ind+1
      IF((ix <= nx2).AND.(iy <= ny2).AND.(iz <= nz2)) THEN
        x1=(ix-nxsh)*dx
        x1t=x1-rvectmp(1)
        x1p=x1-cx(nproj)
        sproj=0D0
#if(parano)
        DO ik=1,nproj_states
          ikk=proj_states(ik)
          sproj=sproj+psi(ind,ikk)*psi(ind,ikk) 
        END DO
#else
        sprojec=0D0 
        DO nbe=1,nstate
          nbee=nrel2abs(nbe)
          DO kk=1,nproj_states
            IF (nbee == proj_states(kk)) THEN
              sprojec=sprojec+psi(ind,nbe)*psi(ind,nbe) 
            END IF
          END DO
        END DO
        CALL mpi_barrier(mpi_comm_world, mpi_ierror)
        CALL mpi_allreduce(sprojec,sproj,1,mpi_double_precision,  &
                mpi_sum,mpi_comm_world,mpi_ierror)
        CALL mpi_barrier(mpi_comm_world, mpi_ierror)
#endif
        starget=rho(ind)-sproj
!                                                     monopole
        qetarget(1)=qetarget(1)+starget
        qeproj(1)=qeproj(1)+sproj
!                                                     dipole
        qetarget(2)=qetarget(2)+starget*x1t
        qetarget(3)=qetarget(3)+starget*y1t
        qetarget(4)=qetarget(4)+starget*z1t

        qeproj(2)=qeproj(2)+sproj*x1p
        qeproj(3)=qeproj(3)+sproj*y1p
        qeproj(4)=qeproj(4)+sproj*z1p
      END IF
    END DO
  END DO
END DO

DO k=1,nrmom
  qetarget(k)=qetarget(k)*dvol
  qeproj(k)=qeproj(k)*dvol
END DO

DO k=2,nrmom
  qetarget(k)=qetarget(k)/qetarget(1)      !normalization
  qeproj(k)=qeproj(k)/qeproj(1)      !normalization
END DO

RETURN
END SUBROUTINE r_projmoms

!-----------------------------------------------------------------------
! COMPLEX version
!-----------------------------------------------------------------------
SUBROUTINE c_projmoms(rho,psi)
USE params
IMPLICIT NONE

REAL(DP), INTENT(IN)    :: rho(2*kdfull2)
COMPLEX(DP), INTENT(IN) :: psi(kdfull2,kstate)
#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

INTEGER :: ind, ix, iy, iz, k
REAL(DP) :: sproj,starget
REAL(DP) :: x1, y1, z1, x1t, y1t, z1t, x1p, y1p, z1p 

#if(parano)
INTEGER :: ik, ikk
#else
INTEGER :: kk, nbe, nbee
REAL(DP) :: sprojec
#endif
!----------------------------------------------------------------------
nrmom=35
IF(nrmom > kmom) STOP ' too many moments in projmoms'

DO k=1,nrmom
  qetarget(k)=0D0
  qeproj(k)=0D0
END DO

!     switch for calculating moments relative to center of mass (1)
!     or center of box (0)
rvectmp = 0D0
IF(iemomsrel == 1 .AND. nion2 > 0) CALL getcm(1,0,0)

ind=0
DO iz=minz,maxz
  z1=(iz-nzsh)*dz
  z1t=z1-rvectmp(3)
  z1p=z1-cz(nproj)
  DO iy=miny,maxy
    y1=(iy-nysh)*dy
    y1t=y1-rvectmp(2)
    y1p=y1-cy(nproj)
    DO ix=minx,maxx
      ind=ind+1
      IF((ix <= nx2).AND.(iy <= ny2).AND.(iz <= nz2)) THEN
        x1=(ix-nxsh)*dx
        x1t=x1-rvectmp(1)
        x1p=x1-cx(nproj)
        sproj=0D0
#if(parano)
        DO ik=1,nproj_states
          ikk=proj_states(ik)
          sproj=sproj+REAL(CONJG(psi(ind,ikk))*psi(ind,ikk),DP)
        END DO
#else
        sprojec=0D0 
        DO nbe=1,nstate
          nbee=nrel2abs(nbe)
          DO kk=1,nproj_states
            IF (nbee == proj_states(kk)) THEN
              sprojec=sprojec+REAL(CONJG(psi(ind,nbe))*psi(ind,nbe),DP)
            END IF
          END DO
        END DO
        CALL mpi_barrier(mpi_comm_world, mpi_ierror)
        CALL mpi_allreduce(sprojec,sproj,1,mpi_double_precision,  &
          mpi_sum,mpi_comm_world,mpi_ierror)
        CALL mpi_barrier(mpi_comm_world, mpi_ierror)
#endif
        starget=rho(ind)-sproj
!                                                     monopole
        qetarget(1)=qetarget(1)+starget
        qeproj(1)=qeproj(1)+sproj
!                                                     dipole
        qetarget(2)=qetarget(2)+starget*x1t
        qetarget(3)=qetarget(3)+starget*y1t
        qetarget(4)=qetarget(4)+starget*z1t

        qeproj(2)=qeproj(2)+sproj*x1p
        qeproj(3)=qeproj(3)+sproj*y1p
        qeproj(4)=qeproj(4)+sproj*z1p
      END IF
    END DO
  END DO
END DO

DO k=1,nrmom
  qetarget(k)=qetarget(k)*dvol
  qeproj(k)=qeproj(k)*dvol
END DO

DO k=2,nrmom
  qetarget(k)=qetarget(k)/qetarget(1)      !normalization
  qeproj(k)=qeproj(k)/qeproj(1)      !normalization
END DO

RETURN
END SUBROUTINE c_projmoms




!------------------------------------------------------------

SUBROUTINE shiftsmall(q0,shix,shiy,shiz)
!------------------------------------------------------------
USE params
USE kinetic
IMPLICIT NONE

! Shifts complex wavefunction array 'q0' by 'shix', 'shiy', 'shiz' in
! x-, y-, and z-direction.
! The shift is produced by an exponential in Fourier space.


COMPLEX(DP), INTENT(IN OUT)                     :: q0(kdfull2)
REAL(DP), INTENT(IN)                         :: shix
REAL(DP), INTENT(IN)                         :: shiy
REAL(DP), INTENT(IN)                         :: shiz

INTEGER :: ind
REAL(DP) :: dkx, dky, dkz
COMPLEX(DP),ALLOCATABLE ::  q1(:)


IF(.NOT.ALLOCATED(akv)) STOP "SHIFTSMALL requires FFT"
ALLOCATE(q1(kdfull2))

CALL fftf(q0,q1)


dkx=pi/(dx*REAL(nx,DP))
dky=pi/(dy*REAL(ny,DP))
dkz=pi/(dz*REAL(nz,DP))

DO ind=1,kdfull2
  q1(ind) = q1(ind)*EXP((shix*akx(ind)+shiy*aky(ind)+shiz*akz(ind)))
ENDDO

CALL fftback(q1,q0)

DEALLOCATE(q1)

RETURN
END SUBROUTINE shiftsmall
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE dipole_qp(q0,shix,shiy,shiz,bostx,bosty,bostz)
!------------------------------------------------------------
USE params
USE kinetic
IMPLICIT NONE

! Initialization of a dipole deformation in space and momentum.
! Shifts complex wavefunction array 'q0' by 'shix', 'shiy', 'shiz' in
! x-, y-, and z-direction and boosts it by 'bostx', 'bosty' and 'bostz'.


INTEGER :: ind, ix, iy, iz, nb
REAL(DP) :: arg, x1, y1, z1
COMPLEX(DP), INTENT(IN OUT)   :: q0(kdfull2,kstate)
REAL(DP), INTENT(IN)          :: shix,shiy,shiz,bostx,bosty,bostz

LOGICAL,PARAMETER :: ttest=.true.
REAL(DP),ALLOCATABLE :: rhotest(:)

IF(ttest) THEN
  ALLOCATE(rhotest(kdfull2))
  CALL calcrho(rhotest,q0)
  WRITE(*,*) ' test DIPOLE_QP'
  CALL prifld(rhotest,'rho before ')
END IF

DO nb=1,nstate

!   shift

  CALL  shiftsmall(q0(1,nb),shix,shiy,shiz)

!   boost

  ind = 0
  DO iz=minz,maxz
    z1=(iz-nzsh)*dz
    DO iy=miny,maxy
      y1=(iy-nysh)*dy
      DO ix=minx,maxx
        x1=(ix-nxsh)*dx
        ind = ind + 1
        arg = x1*bostx+y1*bosty+z1*bostz
        q0(ind,nb)= CMPLX(COS(arg),SIN(arg),DP)*q0(ind,nb)
      END DO
    END DO
  END DO

END DO

IF(ttest) THEN
  CALL calcrho(rhotest,q0)
  CALL prifld(rhotest,'rho after  ')
  DEALLOCATE(rhotest)
END IF


RETURN


END SUBROUTINE dipole_qp

SUBROUTINE testcurrent(wfin,iteract)

! ******************************


!  Compute currents and prints them on standard output.
!  

USE params
USE kinetic
IMPLICIT NONE

COMPLEX(DP), INTENT(IN) :: wfin(kdfull2)
INTEGER, INTENT(IN)     :: iteract

REAL(DP), ALLOCATABLE :: curr0(:,:)


ALLOCATE(curr0(kdfull2,3))
CALL calc_current(curr0,wfin)
 IF(iteract==0) THEN
  CALL prifld(curr0(:,1),'j_1')
  CALL prifld(curr0(:,2),'j_2')
  CALL prifld(curr0(:,3),'j_3')
 END IF
WRITE(*,*) 'step nr.:',iteract,', integrated current:', sum(abs(curr0))
DEALLOCATE(curr0)



END SUBROUTINE testcurrent

! ******************************

SUBROUTINE testgradient(wfin)

! ******************************


!  Tests the routines for gradients in x, y, and z direction.
!  The test-field 'wfin' is a complex field.
!  

USE params
USE kinetic
IMPLICIT NONE

COMPLEX(DP), INTENT(IN OUT)                      :: wfin(kdfull2)

COMPLEX(DP), ALLOCATABLE :: wftest(:)

REAL(DP) :: ekintestx,ekintesty,ekintestz,ekintot

IF(.NOT.ALLOCATED(akv)) STOP "TESTGRADIENT requires FFT"
ALLOCATE(wftest(kdfull2))

CALL xgradient_rspace(wfin,wftest)
CALL xgradient_rspace(wftest,wftest)
ekintestx = dvol*SUM(wfin*wftest)

CALL ygradient_rspace(wfin,wftest)
CALL ygradient_rspace(wftest,wftest)
ekintesty = dvol*SUM(wfin*wftest)

CALL zgradient_rspace(wfin,wftest)
CALL zgradient_rspace(wftest,wftest)
ekintestz = dvol*SUM(wfin*wftest)

CALL fftf(wfin,wftest)
wftest = akv*wftest
CALL fftback(wftest,wftest)
ekintot = dvol*SUM(wfin*wftest)

WRITE(6,'(a,5(1pg13.5))') ' test gradients:', &
  ekintestx,ekintesty,ekintestz,ekintestx+ekintesty+ekintestz, &
  ekintot
WRITE(7,'(a,5(1pg13.5))') ' test gradients:', &
  ekintestx,ekintesty,ekintestz,ekintestx+ekintesty+ekintestz, &
  ekintot


DEALLOCATE(wftest)

END SUBROUTINE testgradient


COMPLEX(DP) FUNCTION determinant(a,n,np)
! Determinant of a complex matrix 'a'. This function is merely an interface for cludcmp.
USE params, ONLY:DP
IMPLICIT NONE

COMPLEX(DP), INTENT(IN OUT)         :: a(np,np)
INTEGER, INTENT(IN)                :: n
INTEGER, INTENT(IN)                :: np

REAL(DP)     :: d
COMPLEX(DP)  :: det
INTEGER :: indx(n)
INTEGER :: ierror

  CALL cludcmp(a,n,np,indx,d,det,ierror)
  determinant = det

RETURN
END FUNCTION determinant


SUBROUTINE cludcmp(a,n,np,indx,d,det,ierror)
 
!     Lower upper decomposition for a complex matrix:                                     
!     a        matrix to be decomposed, at the end decomposed matrix                       
!     n        actual dimension of matrix                                                  
!     np       physical dimension of matrix                                               
!     indx     array keeping the permutations from pivoting                               
!     d        sign  of number of row interchanges                                        
!     det      determinant of 'a'    
USE params, ONLY:DP
IMPLICIT NONE

COMPLEX(DP), INTENT(IN OUT)         :: a(np,np)
INTEGER, INTENT(IN)                 :: n
INTEGER, INTENT(IN)                 :: np
INTEGER, INTENT(OUT)                :: indx(n)
REAL(DP), INTENT(OUT)               :: d
COMPLEX(DP), INTENT(OUT)            :: det
INTEGER, INTENT(OUT)                :: ierror

INTEGER, PARAMETER :: nmax=500
DOUBLE PRECISION, PARAMETER :: tiny=1.0D-20
INTEGER :: i,imax,j,k
COMPLEX(DP) :: csum,cdum
REAL(DP) :: dum,aamax,vv(nmax)

IF(np > nmax) STOP 'LUDCMP: too large matrix'
ierror = 0

d=1.d0
DO  i=1,n
  aamax=0.d0
  DO  j=1,n
    IF (ABS(a(i,j)) > aamax) aamax=ABS(a(i,j))
  END DO
  IF (aamax == CMPLX(0D0,0D0,DP)) THEN
    STOP 'singular matrix in cludcmp. Press ENTER to continue.'
!~     WRITE(6,*) 'singular matrix in cludcmp. Press ENTER to continue.'
!~     READ(*,*)
    STOP 'singular matrix in cludcmp. Press ENTER to continue.'
    ierror=99
  ENDIF
  vv(i)=1.d0/aamax
END DO
DO  j=1,n
  DO  i=1,j-1
    csum=a(i,j)
    DO  k=1,i-1
      csum=csum-a(i,k)*a(k,j)
    END DO
    a(i,j)=csum
  END DO
  aamax=0.d0
  DO  i=j,n
    csum=a(i,j)
    DO  k=1,j-1
      csum=csum-a(i,k)*a(k,j)
    END DO
    a(i,j)=csum
    dum=vv(i)*ABS(csum)
    IF (dum >= aamax) THEN
      imax=i
      aamax=dum
    END IF
  END DO
  IF (j /= imax)THEN
    DO  k=1,n
      cdum=a(imax,k)
      a(imax,k)=a(j,k)
      a(j,k)=cdum
    END DO
    d=-d
    vv(imax)=vv(j)
  END IF
  indx(j)=imax
  IF(a(j,j) == CMPLX(0D0,0D0,DP)) a(j,j)=tiny
  IF(j /= n)THEN
    cdum=1.d0/a(j,j)
    DO  i=j+1,n
      a(i,j)=a(i,j)*cdum
    END DO
  END IF
END DO

!calculate determinant
det=CMPLX(d,0D0,DP)
DO i=1,n
  det=det*a(i,i)
ENDDO
RETURN
END SUBROUTINE cludcmp

END MODULE util
