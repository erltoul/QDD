! Package containing absorbing boundary conditions and
! subsequent analysis of observables from electron emission.


! Editorial remarks:
!  The treatment of spherical absorbing bounds in 'spherAbso'
!  uses a full 3D mask array while cartesian absorbing bounds
!  in 'abso' act merely in boundary stripes. The handling of
!  3D arrays, although somewhat more space consuming, is much
!  simpler. The code should be modified to use only these
!  3D array. The distinction between spherical and cartesian
!  absorbing bound appears then only in the initialization
!  of the mask array.

!  One may also consider to perform the cumulation of total
!  absorbed density within the absorbing loop.



#include"define.h"



!-----absbc----------------------------------------------------------

SUBROUTINE absbc(psi,rho,it)

!     apply absorbing bounds, optionally accumulate absorbed density

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN OUT)                     :: psi(kdfull2,kstate)
REAL(DP), INTENT(IN)                         :: rho(2*kdfull2)
INTEGER, INTENT(IN OUT)                  :: it

REAL(DP),DIMENSION(:),ALLOCATABLE :: w3


LOGICAL :: firstcall=.true.

!------------------------------------------------------------

! WRITE(*,*) ' ABSBC first line'

IF(nabsorb <= 0) RETURN


!     optional initialization of mask

IF (firstcall) THEN
  IF(ispherabso /= 0) THEN
    CALL init_spherabso()
!    WRITE(*,*) ' ABSBC spherical initialization'
  ELSE
    CALL init_abso()
!    WRITE(*,*) ' ABSBC cartesian initialization'
  END IF
  firstcall = .false.
!  WRITE(*,*) ' ABSBC: first call'
END IF



!     save old density

IF(myn == 0) THEN
  ALLOCATE(w3(kdfull2))
!  WRITE(*,*) ' ABSBC: after allocate'
  w3 = rho(1:kdfull2)
!  WRITE(*,*) ' ABSBC: after w3=rho'
END IF

!     apply mask function (and accumulate absorption per state)

IF(jescmaskorb /=0) THEN
  DO nbe=1,nstate
!  WRITE(*,*) ' ABSBC-escmascorb: nbe=',nbe
    DO ind=1,kdfull2
      IF(tgridabso(ind)) THEN
        rhoabsoorb(ind,nbe) = rhoabsoorb(ind,nbe) + spherloss(ind)*(  &
            REAL(psi(ind,nbe))*REAL(psi(ind,nbe))  &
            +imag(psi(ind,nbe))*imag(psi(ind,nbe))  )
      END IF
    END DO
!  WRITE(*,*) ' ABSBC-escmaskorb: nbe=',nbe
  END DO
END IF

DO nbe=1,nstate                                            ! cPW diesen Teil verschoben
!  WRITE(*,*) ' ABSBC: nbe=',nbe
  DO ind=1,kdfull2
    IF(tgridabso(ind)) THEN
      psi(ind,nbe)=sphermask(ind)*psi(ind,nbe)
    END IF
  END DO
!  WRITE(*,*) ' ABSBC: nbe=',nbe
END DO

!         call abso(psi)

!     accumulate absorbed density

CALL calcrho(rho,psi)
IF(myn == 0) THEN
  DO ind=1,kdfull2
    rhoabso(ind)=rhoabso(ind)+w3(ind)-rho(ind)
  END DO
  DEALLOCATE(w3)
END IF



RETURN
END SUBROUTINE absbc


!-----init_abs_accum----------------------------------------------------------

SUBROUTINE init_abs_accum()

!     Initializes accumulator for absorbed densities with zero.

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

!--------------------------------------------------------------------

IF(nabsorb <= 0) RETURN

DO ind=1,kdfull2
  rhoabso(ind)=0D0
END DO


IF(jescmaskorb /=0)  rhoabsoorb(1:kdfull2,1:nstate)=0D0

RETURN
END SUBROUTINE init_abs_accum


!-----init_absbc-----------------------------------------------------

SUBROUTINE init_absbc(rho)

!     initializes geometry parameters for absorbing bounds

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)


REAL(DP), INTENT(IN)                         :: rho(2*kdfull2)


IF (iangabso == 1) THEN ! origin is at center of density
  
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
        codx = codx + rho(ind)*x1
        cody = cody + rho(ind)*y1
        codz = codz + rho(ind)*z1
        
      END DO
    END DO
  END DO
  
  codx = codx * dvol
  cody = cody * dvol
  codz = codz * dvol
  
  xango=codx
  yango=cody
  zango=codz
  
ELSE IF (iangabso == 2) THEN ! origin is at center of cluster mass
  
  xango=0D0
  yango=0D0
  zango=0D0
  
  DO i=1,nion
    xango=xango+cx(i)
    yango=yango+cy(i)
    zango=zango+cz(i)
  END DO
  
  xango=xango/nion
  yango=yango/nion
  zango=zango/nion

ELSE
  
  xango=0D0
  yango=0D0
  zango=0D0
  
END IF
! //

IF (iangabso /= 0 .AND. ipes /= 0) THEN
! prepare arrays for photo-electron-spectra
! get distance from polar origin to edge of box
  
  rmin = 1D10
  
  IF (ABS(xango-(minx-nxsh)*dx)-dx*(nabsorb+1) < rmin)  &
      rmin=ABS(xango-(minx-nxsh)*dx)
  IF (ABS(xango-(maxx-nxsh)*dx-dx*(nabsorb+1)) < rmin)  &
      rmin=ABS(xango-(maxx-nxsh)*dx)
  IF (ABS(yango-(miny-nysh)*dy-dy*(nabsorb+1)) < rmin)  &
      rmin=ABS(yango-(miny-nysh)*dy)
  IF (ABS(yango-(maxy-nysh)*dy-dy*(nabsorb+1)) < rmin)  &
      rmin=ABS(yango-(maxy-nysh)*dy)
  IF (ABS(zango-(minz-nzsh)*dz-dz*(nabsorb+1)) < rmin)  &
      rmin=ABS(zango-(minz-nzsh)*dz)
  IF (ABS(zango-(maxz-nzsh)*dz-dz*(nabsorb+1)) < rmin)  &
      rmin=ABS(zango-(maxz-nzsh)*dz)
  
  
! now make sphere of gridpoints with radius rmin
! these grippoints shall serve as measure points for the
! outgoing waves
  
  jj=0
  
  DO itheta=1,nangtheta
    
    tt = (angthetah-angthetal)/(nangtheta-1) * (itheta-1) + angthetal
    
    
    DO iphi=1,nangphi
      
      jj=jj+1
      
      pp = (angphih-angphil)/(nangphi-1) * (iphi-1) + angphil
      
      xn = rmin*COS(pp)*SIN(tt)
      yn = rmin*SIN(pp)*SIN(tt)
      zn = rmin*COS(tt)
      
      indicesmp(jj) = getnearestgridpoint(xn,yn,zn)
      
    END DO
  END DO
  
END IF

RETURN
END SUBROUTINE init_absbc


!------------------------------------------------------------

SUBROUTINE init_spherabso()

!     Initializes mask function for spherical boundary conditions

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

!------------------------------------------------------------

WRITE(6,*) 'x,y,zango', xango,yango,zango

bcrad = nabsorb*dx
dmin2 = 1D10

!IF (ABS((minz-nzsh)*dz-zango) < dmin2) dmin2=ABS((minz-nzsh)*dz-zango)
IF (ABS((maxz-nzsh)*dz-zango) < dmin2) dmin2=ABS((maxz-nzsh)*dz-zango)
!IF (ABS((miny-nysh)*dy-yango) < dmin2) dmin2=ABS((miny-nysh)*dy-yango)
IF (ABS((maxy-nysh)*dy-yango) < dmin2) dmin2=ABS((maxy-nysh)*dy-yango)
IF (ABS((maxx-nxsh)*dx-xango) < dmin2) dmin2=ABS((maxx-nxsh)*dx-xango)
!IF (ABS((minx-nxsh)*dx-xango) < dmin2) dmin2=ABS((minx-nxsh)*dx-xango)

WRITE(6,*) 'Setting spherical absorbing mask...'
WRITE(6,*) 'Distance to edge of box: '  ,dmin2

DMIN1 = dmin2 - bcrad

IF (DMIN1 < 0D0) STOP 'Error in abso: dmin1<0'

dmin22=dmin2**2
dmin12=DMIN1**2

!         sum = 0D0
ind = 0
DO iz=minz,maxz
  z1=(iz-nzsh)*dz
  rz=z1-zango
  DO iy=miny,maxy
    y1=(iy-nysh)*dy
    ry=y1-yango
    DO ix=minx,maxx
      x1=(ix-nxsh)*dx
      rx=x1-xango
      dist2 = rx*rx+ry*ry+rz*rz
      
      ind=ind+1
      IF (dist2 <= dmin12) THEN
        sphermask(ind)=1.0
        tgridabso(ind)=.false.
      ELSE IF (dist2 > dmin22) THEN
        sphermask(ind)=0D0
        tgridabso(ind)=.true.
      ELSE
        dist = MAX(small,SQRT(dist2))
        cosact = COS((dist-DMIN1)*0.5D0*pi/bcrad)
        IF(cosact > 0D0) THEN
          sphermask(ind)= cosact**powabso
        ELSE
          sphermask(ind)= 0D0
        END IF
        tgridabso(ind)=.true.
      END IF
      spherloss(ind) = 1D0-sphermask(ind)**2
!               sum = sum + spherloss(ind)
    END DO
  END DO
END DO
!      write(6,'(a,1pg13. 5)') ' summed loss=',sum*dvol

CALL prifld(sphermask,' mask ')

DO i=1,nmps
  WRITE(*,*) ' pt.,mask=',imps(i),sphermask(imps(i)),tgridabso(imps(i))
END DO

RETURN
END SUBROUTINE init_spherabso


!-----abso------------------------------------------------------------

SUBROUTINE init_abso()

!     Initializes mask for rectangular absorbing boundaries conditions

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
REAL(DP) :: xmask(nx2),ymask(ny2),zmask(nz2)

LOGICAL,PARAMETER :: wflag=.true.

!-------------------------------------------------------------------------

!     prepare mask functions in each direction separately

bcrad = nabsorb*dx

DO iz = 1,nz2
  zmask(iz) = 1.0D0
END DO
DO iz = 1,nabsorb
  zmask(iz) = COS(pi*0.5D0*(nabsorb+1.0-iz)/nabsorb)**powabso
  zmask(nz2+1-iz) = zmask(iz)
END DO
DO iy = 1,ny2
  ymask(iy) = 1.0D0
END DO
DO iy = 1,nabsorb
  ymask(iy) = COS(pi*0.5D0*(nabsorb+1.0-iy)/nabsorb)**powabso
  ymask(ny2+1-iy) = ymask(iy)
END DO
DO ix = 1,nx2
  xmask(ix) = 1.0D0
END DO
DO ix = 1,nabsorb
  xmask(ix) = COS(pi*0.5D0*(nabsorb+1.0D0-ix)/nabsorb)**powabso
  xmask(nx2+1-ix) = xmask(ix)
END DO
IF(wflag) THEN
  WRITE(6,'(a)') ' ZMASK:'
  WRITE(6,'(1x,5(1pg12.4))') zmask
  WRITE(6,'(a)') ' YMASK:'
  WRITE(6,'(1x,5(1pg12.4))') ymask
  WRITE(6,'(a)') ' XMASK:'
  WRITE(6,'(1x,5(1pg12.4))') xmask
END IF


!     compose to one mask function on all grid points

ind = 0
DO iz=minz,maxz
  DO iy=miny,maxy
    DO ix=minx,maxx
      ind=ind+1
      sphermask(ind)=xmask(ix)*ymask(iy)*zmask(iz)
      tgridabso(ind)= sphermask(ind) < 0.999999999999D0
      spherloss(ind) = 1D0-sphermask(ind)**2
    END DO
  END DO
END DO

RETURN
END SUBROUTINE init_abso

!------------------------------------------------------------

SUBROUTINE initmeasurepoints
!------------------------------------------------------------
USE params
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER :: getnearestgridpoint

dmin2 = 1D10

!IF (ABS((minz-nzsh)*dz-zango) < dmin2) dmin2=ABS((minz-nzsh)*dz-zango)
IF (ABS((maxz-nzsh)*dz-zango) < dmin2) dmin2=ABS((maxz-nzsh)*dz-zango)
!IF (ABS((miny-nysh)*dy-yango) < dmin2) dmin2=ABS((miny-nysh)*dy-yango)
IF (ABS((maxy-nysh)*dy-yango) < dmin2) dmin2=ABS((maxy-nysh)*dy-yango)
IF (ABS((maxx-nxsh)*dx-xango) < dmin2) dmin2=ABS((maxx-nxsh)*dx-xango)
!IF (ABS((minx-nxsh)*dx-xango) < dmin2) dmin2=ABS((minx-nxsh)*dx-xango)

!r=dmin2-bcrad-dx
bcrad = nabsorb*dx
r=dmin2-bcrad
WRITE(*,*) ' analyzing point at r,ir=',r,nint(r/dx)

nmps=0

DO ith=1,nmptheta+1
  DO iph=1,nmpphi
    
    p=(iph-1)*2D0*pi/nmpphi
    t=(ith-1)*pi/nmptheta
    
    nmps=nmps+1

    IF(nmps > maxmps) STOP ' array for analyzing points exhausted'
    
    x=r*COS(p)*SIN(t)
    y=r*SIN(p)*SIN(t)
    z=r*COS(t)
    
    
    imps(nmps) = getnearestgridpoint(x,y,z)
    DO ik=1,nmps-1                              
      IF(imps(nmps) == imps(ik)) THEN           
       nmps = nmps - 1                          
      END IF                                    
    END DO                                      
    
  END DO
END DO

!WRITE(*,*) ' mask in INITMEASUREPOINTS'
!CALL prifld(sphermask,' mask ')

DO ik=1,nmps                                    
 impsact = imps(nmps)
 impsx = mod(impsact-1,nx2)+1
 impsy = mod(impsact/nx2,ny2)+1
 impsz = mod(impsact/(nx2*ny2),nz2)+1
 WRITE(*,*) ' analyzing pts: nmps,imps=',nmps,impsact
 WRITE(*,*) ' nx,ny,nz=',impsx,impsy,impsz
 WRITE(*,*) ' x,y,z=',x,y,z
 WRITE(*,*) ' theta,phi=',t,p
END DO                                          



RETURN
END SUBROUTINE initmeasurepoints
!------------------------------------------------------------


!-----escmask---------------------------------------------------------

SUBROUTINE escmask(it)

!     print collected information on escaping electrons


USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: it
CHARACTER (LEN=1) :: inttostring1
CHARACTER (LEN=2) :: inttostring2
CHARACTER (LEN=3) :: inttostring3

CHARACTER (LEN=1) :: str1
CHARACTER (LEN=2) :: str2
CHARACTER (LEN=3) :: str3
CHARACTER (LEN=4) :: str4
CHARACTER (LEN=5) :: str5
CHARACTER (LEN=6) :: str6
CHARACTER (LEN=7) :: str7
CHARACTER (LEN=8) :: str8
CHARACTER (LEN=9) :: str9

!--------------------------------------------------------------------


IF (jescmaskorb /=0 .AND. MOD(it,jescmaskorb) == 0) THEN
  DO nbe=1,nstate
    nbeabs = nrel2abs(nbe)
    IF (nbeabs < 10) THEN
      str1=inttostring1(nbeabs)
      OPEN(588,STATUS='unknown', FILE='pescmaskOrb.'//str1//'.'//outnam)
    ELSE IF (nbeabs < 100) THEN
      str2=inttostring2(nbeabs)
      OPEN(588,STATUS='unknown', FILE='pescmaskOrb.'//str2//'.'//outnam)
    ELSE IF (nbeabs < 1000) THEN
      str3=inttostring3(nbeabs)
      OPEN(588,STATUS='unknown', FILE='pescmaskOrb.'//str3//'.'//outnam)
    ELSE
      STOP 'ERROR: Too many states for jescmaskOrb'
    END IF
    CALL printfield(588,rhoabsoorb(1,nbe),'x')
    CLOSE(588)
  END DO
END IF

IF(myn == 0) THEN
  IF (jescmask .NE. 0 .AND. MOD(it,jescmask) == 0) THEN
    IF (itof == 0) THEN
      OPEN(589,STATUS='unknown',FILE='pescmask.'//outnam)
    ELSE
      
      it0=it
      ndig=1
      
      DO i=1,10
        IF (it0/10. >= 1.) THEN
          ndig=ndig+1
          it0=INT(it0/10)
        END IF
      END DO
      
      IF (ndig == 1) THEN
        WRITE (str1, '(I1)') it
        OPEN(589,STATUS='unknown', FILE='pescmask.'//str1//'.'//outnam)
      ELSE IF (ndig == 2) THEN
        WRITE (str2, '(I2)') it
        OPEN(589,STATUS='unknown', FILE='pescmask.'//str2//'.'//outnam)
      ELSE IF (ndig == 3) THEN
        WRITE (str3,'(I3)') it
        OPEN(589,STATUS='unknown', FILE='pescmask.'//str3//'.'//outnam)
      ELSE IF (ndig == 4) THEN
        WRITE (str4, '(I4)') it
        OPEN(589,STATUS='unknown', FILE='pescmask.'//str4//'.'//outnam)
      ELSE IF (ndig == 5) THEN
        WRITE (str5, '(I5)') it
        OPEN(589,STATUS='unknown', FILE='pescmask.'//str5//'.'//outnam)
      ELSE IF (ndig == 6) THEN
        WRITE (str6, '(I6)') it
        OPEN(589,STATUS='unknown', FILE='pescmask.'//str6//'.'//outnam)
      ELSE IF (ndig == 7) THEN
        WRITE (str7, '(I7)') it
        OPEN(589,STATUS='unknown', FILE='pescmask.'//str7//'.'//outnam)
      ELSE IF (ndig == 8) THEN
        WRITE (str8, '(I8)') it
        OPEN(589,STATUS='unknown', FILE='pescmask.'//str8//'.'//outnam)
      ELSE IF (ndig == 9) THEN
        WRITE (str9, '(I9)') it
        OPEN(589,STATUS='unknown', FILE='pescmask.'//str9//'.'//outnam)
      ELSE
        STOP '::too many time steps::'
      END IF
      
    END IF
    CALL printfield(589,rhoabso,'x')
    CLOSE(589)
  END IF


  IF (it == 0) THEN
    OPEN(589,STATUS='unknown',FILE='pescmask.0.'//outnam)
    CALL printfield(589,rhoabso,'x')
    CLOSE(589)
  END IF

END IF

RETURN
END SUBROUTINE escmask


!------------------------------------------------------------

SUBROUTINE angabso(it)
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)



INTEGER, INTENT(IN OUT)                  :: it
REAL(DP) :: absocone(maxnang)

DO jj=1,maxnang
  absocone(jj)=0D0
END DO

jj=0

DO itheta=1,nangtheta
  
  tt = (angthetah-angthetal)/(nangtheta-1) * (itheta-1) + angthetal
  
  DO iphi=1,nangphi
    
    jj=jj+1
    
    pp = (angphih-angphil)/(nangphi-1) * (iphi-1) + angphil
    
    xn = COS(pp)*SIN(tt)
    yn = SIN(pp)*SIN(tt)
    zn = COS(tt)
    
    ind=0
    DO iz=minz,maxz
      z1=(iz-nzsh)*dz
      DO iy=miny,maxy
        y1=(iy-nysh)*dy
        DO ix=minx,maxx
          x1=(ix-nxsh)*dx
          
          ind=ind+1
          
          tau = xn*(x1-xango)+ yn*(y1-yango)+  &
              zn*(z1-zango)
          
          IF (tau >= 0D0) THEN
            
            x = xango + tau*xn
            y = yango + tau*yn
            z = zango + tau*zn
            
            dist2 = (x-x1)*(x-x1)+(y-y1)*(y-y1)+ (z-z1)*(z-z1)
            
            dia2 = TAN(delomega/2)
            dia2 = dia2*dia2
            dia2 = dia2*( (x-xango)*(x-xango)+  &
                (y-yango)*(y-yango) + (z-zango)*(z-zango))
            
            IF (dist2 <= dia2) THEN ! Hit! Grid point is in cone
              
              absocone(jj)=absocone(jj)+rhoabso(ind)
              
            END IF
            
          END IF
          
        END DO
      END DO
    END DO
    
    WRITE(47,'(1f14.5,2f12.3,1e17.5)') tfs,tt/pi*180.0,  &
        pp/pi*180.0,absocone(jj)
    
  END DO
  
END DO


RETURN
END SUBROUTINE angabso
!------------------------------------------------------------



!-----nesacpe------------------------------------------------

SUBROUTINE nescape(it,rho)

!     Compute and print total number of escaped electrons

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)


INTEGER, INTENT(IN)                      :: it
REAL(DP), INTENT(IN)                         :: rho(kdfull2)


!------------------------------------------------------------


tinfs=it*dt1*0.0484/2D0/ame
!apnum=0D0
!DO i=1,nxyz
!  apnum=apnum+rho(i)
!END DO
!apnum=apnum*dvol
apnum = dvol*SUM(rho)
absosum = dvol*SUM(rhoabso)

CALL safeopen(23,it,jesc,'pescel')
!      open(23,position='append',file='pescel.'//outnam)
WRITE(23,'(f13.5,3(1pg13.5))') tinfs,1.0-apnum/nclust,nclust-apnum,absosum
CALL flush(23)

RETURN
END SUBROUTINE nescape

!------------------------------------------------------------

SUBROUTINE evalmp(iunit,q0)
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)


INTEGER, INTENT(IN) :: iunit
COMPLEX(DP), INTENT(IN) :: q0(kdfull2,kstate)

REAL(DP) :: rp(2*maxmps)
!REAL(DP) :: cp(10000)
COMPLEX(DP) :: q0phase
LOGICAL :: topenf
#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
LOGICAL,PARAMETER :: ttestpar=.FALSE.
#endif

IF(myn==0) THEN
  INQUIRE(803,OPENED=topenf)
  IF(.NOT.topenf) &
    OPEN(803,POSITION='append',FORM='unformatted',FILE='pMP.'//outnam)
END IF


#if(parano)
DO nbe=1,nstate
  j=0
  DO i=1,nmps
    ii = imps(i)
    CALL conv1to3(ii)
    x1 = (iindtmp(1)- nxsh)*dx
    y1 = (iindtmp(2)- nysh)*dy
    z1 = (iindtmp(3)- nzsh)*dz
    scal1 = x1 * e1x + y1 * e1y + z1 * e1z
    scal2 = x1 * e2x + y1 * e2y + z1 * e2z
    q0phase = q0(ii,nbe)*  &
             EXP(CMPLX(0D0,-e0*(scal1*fpulseinteg1+scal2*fpulseinteg2)))
    j=j+1
    rp(j)=REAL(q0phase)
    j=j+1
    rp(j)=imag(q0phase)
  END DO
!  WRITE(iunit,'(1f14.5,1000e18.8)') tfs,rp(1:2*nmps)                       ! cPW: maxmps = 500
  WRITE(iunit) tfs,rp(1:2*nmps)                                             ! cPW
END DO
#endif

#if(parayes)
DO nba=1,nstate_all
  mynact = nhome(nba)
  IF(ttestpar) WRITE(*,*) ' before barrier: nba,node=',nba,mynact,myn
  CALL mpi_barrier (mpi_comm_world, mpi_ierror)
  IF(ttestpar) WRITE(*,*) ' after barrier: nba,node=',nba,mynact,myn
  IF(mynact == myn) THEN
    nbe = nabs2rel(nba) 
    j=0
    DO i=1,nmps
      ii = imps(i)
      CALL conv1to3(ii)
      x1 = (iindtmp(1)- nxsh)*dx
      y1 = (iindtmp(2)- nysh)*dy
      z1 = (iindtmp(3)- nzsh)*dz
      scal1 = x1 * e1x + y1 * e1y + z1 * e1z
      scal2 = x1 * e2x + y1 * e2y + z1 * e2z
      q0phase = q0(ii,nbe)*  &
               EXP(CMPLX(0D0,-e0*(scal1*fpulseinteg1+scal2*fpulseinteg2)))
      j=j+1
      rp(j)=REAL(q0phase)
      j=j+1
      rp(j)=imag(q0phase)
    END DO
    IF(myn /= 0) THEN
      CALL mpi_send(rp,2*nmps,mpi_double_precision,0,nba,mpi_comm_world,ic)
      IF(ttestpar) WRITE(*,*) ' sent: nba,node=',nba,myn
    END IF
  END IF
  IF(myn == 0) THEN
    IF(mynact /= 0) THEN
       CALL mpi_recv(rp,2*nmps,mpi_double_precision,mynact,nba,  &
                     mpi_comm_world,is,ic)
       IF(ttestpar) WRITE(*,*) ' received: nba,from node=',nba,mynact
    END IF
!    WRITE(iunit,'(1f14.5,1000e18.8)') tfs,rp(1:2*nmps)                     ! cPW: maxmps = 500
    WRITE(iunit) tfs,rp(1:2*nmps)                                           ! cPW
  END IF
END DO
IF(ttestpar) WRITE(*,*) ' all done. before last barrier'
CALL mpi_barrier (mpi_comm_world, mpi_ierror)
#endif

IF(myn==0) CLOSE(iunit)


RETURN
END SUBROUTINE evalmp
!------------------------------------------------------------


