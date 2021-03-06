#include "define.h"
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!                 The Erlangen - Toulouse

!            SHORT RANGE INTERACTION PACKAGE

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!------------------------------------------------------------

SUBROUTINE getshortforce(ityp1,ityp2,ind1,ind2,rho,iflag,iflag2)
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
!     distributing short force routine
!     ind1 and ind2 are of the short type!


INTEGER, INTENT(IN OUT)                  :: ityp1
INTEGER, INTENT(IN OUT)                  :: ityp2
INTEGER, INTENT(IN OUT)                  :: ind1
INTEGER, INTENT(IN OUT)                  :: ind2
REAL(DP), INTENT(IN)                         :: rho(kdfull2*2)
INTEGER, INTENT(IN OUT)                  :: iflag
INTEGER, INTENT(IN)                      :: iflag2

INTEGER :: getnearestgridpoint
INTEGER :: conv3to1


rder = 1.0D-5


#if(raregas)

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
IF (isrtyp(ityp1,ityp2) == 1) THEN ! MgO-case
  
  signum = 1D0
  
! sort by particle type
!$$$         if (ityp1.gt.ityp2) then
!$$$            itmp = ityp1
!$$$            ityp1 = ityp2
!$$$            ityp2 = itmp
!$$$            itmp = ind1
!$$$            ind1 = ind2
!$$$            ind2 = itmp
!$$$            signum = -1.0
!$$$c            stop 'warning: signum is -1, CHECK!'
!$$$         endif
  
! determine parameters
  
  
  IF (ityp1 <= 3 .AND. ityp2 <= 3) THEN ! GSM-GSM is of Born-Mayer-type
    
    
    
    itmp1 = iconvshorttolong(ityp1,ind1)
    itmp2 = iconvshorttolong(ityp2,ind2)
    
    
    rr = getdistance2(itmp1,itmp2)
    rr = SQRT(MAX(rr,1D-14))
    
    
    
    
    CALL setbornparas(ityp1,ityp2) ! paras are stored in rVecTmp
    
    
    
    
    radfor = fbornmayer(rr,rvectmp(1),rvectmp(2),rvectmp(3))
    
    
    
    CALL getrelvec(itmp1,itmp2,rr) ! direction stored in rVecTmp
    
    
    
    forx = radfor * rvectmp(1) * signum
    fory = radfor * rvectmp(2) * signum
    forz = radfor * rvectmp(3) * signum
    
    
    
    
    IF (iflag2 == 0) THEN ! as usual
      CALL addforce(itmp1,-forx,-fory,-forz)
      CALL addforce(itmp2,forx,fory,forz)
    ELSE IF (iflag2 == 1) THEN ! ibh second part
      CALL addforce(itmp1,-forx,-fory,-forz)
!               call addForce(itmp2,-forx,-fory,-forz)
    ELSE
      CALL addforce(itmp1,forx,fory,forz)
    END IF
    
    
  ELSE IF ((ityp1 <= 3 .AND. ityp2 == 4) .OR.  &
        (ityp1 == 4 .AND.ityp2 <= 3))  THEN  !GSM-Na or Na-Na force
    
    IF (ityp1 == 1) THEN  !
      
      xr = xc(ind1)-cx(ind2)
      yr = yc(ind1)-cy(ind2)
      zr = zc(ind1)-cz(ind2)
      
      rr = MAX(SQRT(xr*xr+yr*yr+zr*zr),small)
      
      radfor = fbornmayermod(rr,bcn,ecn,fcn,sigcn,ccncud ,ccnd,dcn)
      
      
      forx = radfor * xr / rr * signum
      fory = radfor * yr / rr * signum
      forz = radfor * zr / rr * signum
      
      fxc(ind1)=fxc(ind1)-forx
      fyc(ind1)=fyc(ind1)-fory
      fzc(ind1)=fzc(ind1)-forz
      
      
      fx(ind2)=fx(ind2)+forx
      fy(ind2)=fy(ind2)+fory
      fz(ind2)=fz(ind2)+forz
    ELSE IF (ityp2 == 1) THEN
      xr = cx(ind1)-xc(ind2)
      yr = cy(ind1)-yc(ind2)
      zr = cz(ind1)-zc(ind2)
      
      rr = MAX(SQRT(xr*xr+yr*yr+zr*zr),small)
      
      radfor = fbornmayermod(rr,bcn,ecn,fcn,sigcn,ccncud ,ccnd,dcn)
      
      
      forx = radfor * xr / rr * signum
      fory = radfor * yr / rr * signum
      forz = radfor * zr / rr * signum
      
      fx(ind1)=fx(ind1)-forx
      fy(ind1)=fy(ind1)-fory
      fz(ind1)=fz(ind1)-forz
      
      
      fxc(ind2)=fxc(ind2)+forx
      fyc(ind2)=fyc(ind2)+fory
      fzc(ind2)=fzc(ind2)+forz
      
      
      
      
    ELSE IF (ityp1 == 3) THEN
      
      xr = xk(ind1)-cx(ind2)
      yr = yk(ind1)-cy(ind2)
      zr = zk(ind1)-cz(ind2)
      
      rr = MAX(SQRT(xr*xr+yr*yr+zr*zr),small)
      
!               if(BkN2.eq.0D0)then
      
      radfor = fbornmayermod(rr,bkn,ekn,fkn,sigkn,ckncud, cknd,dkn)
      
!               else
!                radfor = funkderFermi(r,BkN,1D0/sigkN,FkN)+
!     &                    funkderFermi(r,BkN2,1D0/sigkN2,FkN2)
!               endif
      
      forx = radfor * xr / rr * signum
      fory = radfor * yr / rr * signum
      forz = radfor * zr / rr * signum
      
      fxk(ind1)=fxk(ind1)-forx
      fyk(ind1)=fyk(ind1)-fory
      fzk(ind1)=fzk(ind1)-forz
      
      fx(ind2)=fx(ind2)+forx
      fy(ind2)=fy(ind2)+fory
      fz(ind2)=fz(ind2)+forz
      
    ELSE IF (ityp2 == 3) THEN
      
      xr = cx(ind1)-xk(ind2)
      yr = cy(ind1)-yk(ind2)
      zr = cz(ind1)-zk(ind2)
      
      rr = MAX(SQRT(xr*xr+yr*yr+zr*zr),small)
      
!               if(BkN2.eq.0D0)then
      
      radfor = fbornmayermod(rr,bkn,ekn,fkn,sigkn,ckncud, cknd,dkn)
      
!               else
!                radfor = funkderFermi(r,BkN,1D0/sigkN,FkN)+
!     &                    funkderFermi(r,BkN2,1D0/sigkN2,FkN2)
!               endif
      
      
      forx = radfor * xr / rr * signum
      fory = radfor * yr / rr * signum
      forz = radfor * zr / rr * signum
      
      fx(ind1)=fx(ind1)-forx
      fy(ind1)=fy(ind1)-fory
      fz(ind1)=fz(ind1)-forz
      
      fxk(ind2)=fxk(ind2)+forx
      fyk(ind2)=fyk(ind2)+fory
      fzk(ind2)=fzk(ind2)+forz
      
      
    ELSE                ! Na-Na
! no short range interaction as in MgO case
    END IF
    
  ELSE IF (ityp2 == 5) THEN ! GSM-DFT force
    
    IF (ityp1 <= 3) THEN
      
! to be implemented by fitting the Muenchen data
      
      
    ELSE
! nothing to do
    END IF
    
    
  ELSE
    STOP 'Error in MgO case of getShortForce'
  END IF
  
  
  
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ELSE IF (isrtyp(ityp1,ityp2) == 2) THEN ! Ar-case
  
  signum = 1D0
  
! sort by particle type
  IF (ityp1 > ityp2) THEN
    itmp = ityp1
    ityp1 = ityp2
    ityp2 = itmp
    itmp = ind1
    ind1 = ind2
    ind2 = itmp
    signum = -1D0
    STOP 'warning: signum is -1'
  END IF
  
  IF (ityp1 == 1) THEN ! 1=c
    IF (ityp2 == 1) THEN ! 2=c
!    Ar case: Core-Core Interaction
      
      
      xr = xc(ind1)-xc(ind2)
      yr = yc(ind1)-yc(ind2)
      zr = zc(ind1)-zc(ind2)
      
      dist2 = xr*xr+yr*yr+zr*zr
      dist = SQRT(dist2)
      
!               radfor = VArArShort(dist+rder)-VArArShort(dist-rder)
      IF(iararlj == 0) THEN
        radfor = v_ar_ar(dist+rder)-v_ar_ar(dist-rder)
      ELSE
        radfor = vararlj(dist+rder)-vararlj(dist-rder)
      END IF
      
      
      radfor = radfor /2./rder*signum
      
      forx = radfor*xr/dist
      fory = radfor*yr/dist
      forz = radfor*zr/dist
      
      fxc(ind1) = fxc(ind1) - forx
      fyc(ind1) = fyc(ind1) - fory
      fzc(ind1) = fzc(ind1) - forz
      
      fxc(ind2) = fxc(ind2) + forx
      fyc(ind2) = fyc(ind2) + fory
      fzc(ind2) = fzc(ind2) + forz
      
    ELSE IF (ityp2 == 4) THEN ! 2=Na
!     Ar case: Core-Na interaction
      
      
      xr = xc(ind1) - cx(ind2)
      yr = yc(ind1) - cy(ind2)
      zr = zc(ind1) - cz(ind2)
      
!               write(6,*) 'xryrzr:',xr,yr,zr
      
      dist2 = xr*xr+yr*yr+zr*zr
      dist = SQRT(dist2)
      
      radfor = varnashort(dist+rder)-varnashort(dist-rder)
      
      
      radfor = radfor/(rder+rder)*signum
      
      forx = radfor*xr/dist
      fory = radfor*yr/dist
      forz = radfor*zr/dist
      
      
      fxc(ind1) = fxc(ind1) - forx
      fyc(ind1) = fyc(ind1) - fory
      fzc(ind1) = fzc(ind1) - forz
      
      fx(ind2) = fx(ind2) + forx
      fy(ind2) = fy(ind2) + fory
      fz(ind2) = fz(ind2) + forz
      
    ELSE IF (ityp2 == 5) THEN !2=DFT
!     Ar case: Core-electron interaction
      
      ind=0
      sumx = 0D0
      sumy = 0D0
      sumz = 0D0
      
      IF (ipseudo == 0) THEN
        
        DO iz=minz,maxz
          z1=(iz-nzsh)*dz
          zr=z1-zc(ind1)
          
          DO iy=miny,maxy
            y1=(iy-nysh)*dy
            yr=y1-yc(ind1)
            
            DO ix=minx,maxx
              x1=(ix-nxsh)*dx
              xr=x1-xc(ind1)
              
              ind=ind+1
              
              r2=xr*xr + yr*yr + zr*zr
              rr=MAX(SQRT(r2),small)
              
              chpddr = rho(ind)*(varelcore(rr+rder)  &
                  - varelcore(rr-rder))/(rr*2.*rder)
              
! attention again: rho is positive!!!
              
              sumx = sumx - chpddr*xr
              sumy = sumy - chpddr*yr
              sumz = sumz - chpddr*zr
              
              
              
            END DO
          END DO
        END DO
        
      ELSE ! ipseudo =1
        
        ind = getnearestgridpoint(xc(ind1),yc(ind1),zc(ind1))
        
        CALL conv1to3(ind)
        
        nsgsize = 10 ! just for the argon case!
        
        DO iz=iindtmp(3)-nsgsize,iindtmp(3)+nsgsize
          IF (iz >= minz .AND. iz <= maxz) THEN
            z1=(iz-nzsh)*dz
            zr=zc(ind1)-z1
            
            DO iy=iindtmp(2)-nsgsize,iindtmp(2)+nsgsize
              IF (iy >= miny .AND. iy <= maxy) THEN
                
                y1=(iy-nysh)*dy
                yr=yc(ind1)-y1
                
                DO ix=iindtmp(1)-nsgsize,iindtmp(1)+nsgsize
                  IF (ix >= minx .AND. ix <= maxx) THEN
                    
                    x1=(ix-nxsh)*dx
                    xr=xc(ind1)-x1
                    
                    ind = conv3to1(ix,iy,iz)
                    
                    r2=xr*xr + yr*yr + zr*zr
                    rr=MAX(SQRT(r2),small)
                    
                    
                    chpddr = rho(ind)*(varelcore(rr+rder)  &
                        - varelcore(rr-rder))/(rr*2.*rder)
                    
                    
!                        chpddr = rho(ind)*(V_Ar_el_core(rr+rder)
!     &                        - V_Ar_el_core(rr-rder))/(rr*2.*rder)
! attention again: rho is positive!!!
                    
                    sumx = sumx + chpddr*xr
                    sumy = sumy + chpddr*yr
                    sumz = sumz + chpddr*zr
                    
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
        
        
        
      END IF
      
      
      fxc(ind1) = fxc(ind1) + sumx*dvol
      fyc(ind1) = fyc(ind1) + sumy*dvol
      fzc(ind1) = fzc(ind1) + sumz*dvol
      
      
    END IF
  END IF
  
ELSE  !other cases
  
! OTHER CASES TO BE IMPLEMENTED HERE
  
END IF
#endif

RETURN
END SUBROUTINE getshortforce
!------------------------------------------------------------




!------------------------------------------------------------

SUBROUTINE initisrtyp
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
!     Initializes the short range interaction type matrix
!     isrtyp(i,j)

!     i,j each range from 1..5, meaning:

!     i  --->    1   --->  GSM core
!                2   --->  GSM valence shell
!                3   --->  GSM kation
!                4   --->  Na core
!                5   --->  DFT electron

!     The value of the matrix elements can be:

!                0   --->  no short range interaction
!                1   --->  Born-Mayer type
!                2   --->  Argon case

!     The matrix must be symmetric, of course.


!$$$      open(156,status='old',file='for005sr')
!$$$        do i=1,5
!$$$         read(156,*) isrtyp(i,1),isrtyp(i,2),isrtyp(i,3),isrtyp(i,4),
!$$$     &           isrtyp(i,5)
!$$$        enddo
!$$$      close(156)
!     OBSOLETE by namelists!!


!     check consistency
#if(raregas)
IF (isrtypall /= 0) THEN
  WRITE(6,*) 'Setting short range parameters...'
  DO i=1,5
    DO j=1,5
      isrtyp(i,j)=isrtypall
    END DO
  END DO
  
END IF


DO i=1,5
  DO j=1,5
    IF (isrtyp(i,j) /= isrtyp(j,i)) THEN
      STOP 'Short range interaction type matrix not symmetric'
    END IF
  END DO
END DO

!     check if it is necessary to handle function parameters by input file


DO i=1,5
  DO j=i,5
    IF (isrtyp(i,j) == 3) THEN
      IF (fermia == 0D0) STOP 'no parameters for fermi function given'
    END IF
  END DO
END DO
#endif

RETURN
END SUBROUTINE initisrtyp
!------------------------------------------------------------

#if(raregas)
!------------------------------------------------------------

SUBROUTINE getrelvec(ind1,ind2,r)
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

CALL getparas(ind1)
x1 = rvectmp(1)
y1 = rvectmp(2)
z1 = rvectmp(3)

CALL getparas(ind2)

IF (r > 0) THEN
  rvectmp(1) = (x1 - rvectmp(1))/r
  rvectmp(2) = (y1 - rvectmp(2))/r
  rvectmp(3) = (z1 - rvectmp(3))/r
ELSE
  r = MAX(SQRT(getdistance2(ind1,ind2)),small)
  rvectmp(1) = (x1 - rvectmp(1))/r
  rvectmp(2) = (y1 - rvectmp(2))/r
  rvectmp(3) = (z1 - rvectmp(3))/r
END IF

RETURN
END SUBROUTINE getrelvec
!------------------------------------------------------------

!-----VArArShort------------------------------------------------
!     only short range part of Ar-Ar interaction

FUNCTION vararshort(r)

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

!     Ar-Ar potential


REAL(DP), INTENT(IN OUT)                     :: r
DATA  alpha_ar,beta_ar, c6_ar, c8_ar,  a_ar  &
    /  1.7301, 1.7966,55.465,3672.9,794.21/


rabs = MAX(ABS(r),small)
core = a_ar*EXP(-alpha_ar*rabs)/rabs -2.0/(1D0+EXP((beta_ar/rabs)**2))  &
    *(c6_ar/rabs**6+c8_ar/rabs**8)

vararshort = e2*core

RETURN
END FUNCTION vararshort
!------------------------------------------------------------

!-----VArArLJ------------------------------------------------
!     only short range part of Ar-Ar interaction

FUNCTION vararlj(r)

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

!     Ar-Ar potential


REAL(DP), INTENT(IN OUT)                     :: r
DATA  epslj,sigmalj /  0.0007647,6.42503/


rabs = MAX(ABS(r),small)

vararlj = 4.0*epslj*((sigmalj/rabs)**12.0- (sigmalj/rabs)**6.0)

RETURN
END FUNCTION vararlj
!------------------------------------------------------------


!------------------------------------------------------------

FUNCTION vbornmayer(r,a,s,c)
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)


vbornmayer = a*EXP(-r/s) - c/r**6


RETURN
END FUNCTION vbornmayer
!------------------------------------------------------------

!------------------------------------------------------------

FUNCTION fbornmayer(r,a,s,c)
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

rder=1D-5

!      ftmp = VBornMayer(r+rder,A,s,C) - VBornMayer(r-rder,A,s,C)
!      ftmp = ftmp/(2.0*rder)

fbornmayer = -a/s*EXP(-r/s) + 6*c*(erf(r)/r)**7

!      FBornMayer = ftmp

RETURN
END FUNCTION fbornmayer
!------------------------------------------------------------

!------------------------------------------------------------

FUNCTION fbornmayer46(r,a,s,c4,c6)
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)



fbornmayer46 = -a/s*EXP(-r/s) +4*c4/r**5 + 6*c6/r**7

!      FBornMayer = ftmp

RETURN
END FUNCTION fbornmayer46
!------------------------------------------------------------

!------------------------------------------------------------

FUNCTION fbornmayer2345678d(r,a,e,f,s,s2,c2,c3,c4,c5,c6,c7,c8, c9,c10,cd,d)
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

rder=1D-5

acc=0D0

sum2= -a/s*EXP((r-f)/s)/(e+EXP((r-f)/s))**2

rrp=rr+rder
rrm=rr-rder

IF (c2 /= 0D0) acc=acc+c2*((erf(rrp/s2)/rrp)**2 -(erf(rrm/s2)/rrm)**2)
IF (c3 /= 0D0) acc=acc+c3*((erf(rrp/s2)/rrp)**3 -(erf(rrm/s2)/rrm)**3)
IF (c4 /= 0D0) acc=acc+c4*((erf(rrp/s2)/rrp)**4 -(erf(rrm/s2)/rrm)**4)
IF (c5 /= 0D0) acc=acc+c5*((erf(rrp/s2)/rrp)**5 -(erf(rrm/s2)/rrm)**5)
IF (c6 /= 0D0) acc=acc+c6*((erf(rrp/s2)/rrp)**6 -(erf(rrm/s2)/rrm)**6)
IF (c7 /= 0D0) acc=acc+c7*((erf(rrp/s2)/rrp)**7 -(erf(rrm/s2)/rrm)**7)
IF (c8 /= 0D0) acc=acc+c8*((erf(rrp/s2)/rrp)**8 -(erf(rrm/s2)/rrm)**8)
IF (c9 /= 0D0) acc=acc+c9*((erf(rrp/s2)/rrp)**9 -(erf(rrm/s2)/rrm)**9)
IF (c10 /= 0D0) acc=acc+c10*((erf(rrp/s2)/rrp)**10 -(erf(rrm/s2)/rrm)**10)
IF (cd /= 0D0) acc=acc+cd*((erf(rrp/s2)/rrp)**d -(erf(rrm/s2)/rrm)**d)


acc = acc/2./rder

!$$$      if (C2.ne.0.0) acc=acc+2*C2/r**3
!$$$      if (C3.ne.0.0) acc=acc+3*C3/r**4
!$$$      if (C4.ne.0.0) acc=acc+4*C4/r**5
!$$$      if (C5.ne.0.0) acc=acc+5*C5/r**6
!$$$      if (C6.ne.0.0) acc=acc+6*C6/r**7
!$$$      if (C7.ne.0.0) acc=acc+7*C7/r**8
!$$$      if (C8.ne.0.0) acc=acc+8*C8/r**9
!$$$      if (C9.ne.0.0) acc=acc+9*C9/r**10
!$$$      if (C10.ne.0.0) acc=acc+10*C10/r**11
!$$$      if (CD.ne.0.0) acc=acc+D*CD/r**(D+1)


fbornmayer2345678d = acc+sum2

RETURN
END FUNCTION fbornmayer2345678d
!------------------------------------------------------------

!------------------------------------------------------------

FUNCTION fbornmayermod(r,a,e,f,s,s2,cd,d)
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

rder=1D-5


acc=0D0

sum2= -a/s*EXP((r-f)/s)/(e+EXP((r-f)/s))**2

rrp=r+rder
rrm=r-rder

IF (cd /= 0D0) THEN
  acc=acc-cd*((erf(rrp/s2)/rrp)**d -(erf(rrm/s2)/rrm)**d)
END IF


acc = acc/2./rder

fbornmayermod = acc+sum2


RETURN
END FUNCTION fbornmayermod
!------------------------------------------------------------

!------------------------------------------------------------

FUNCTION varnashort(r)
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

!     Ar-Na potential

!     parameters for repulsive core  (energies in Ha, rescale to Ry by e2 !)

!      data vrep/214.14/          ! strength of repulsive core
!      data rcor/1.6906/          ! inverse radius of repulsive core
!      data vdip2/308.98/         ! 2*strength of higher dipole part
!      data rdip2/2.20/           ! radius of cutoff in higher dipole part

!----------------------------------------------------------------------

!      rabs = abs(r)

!      V_Ar_Na =
!     &         e2*(vrep*exp(-rcor*rabs)/rabs
!     &             - vdip2/((1.0d0+exp(rdip2/rabs)**2)*rabs**6)
!     &             + ch(18)*ch(11)*V_soft(r,sgm1(18))
!     &            )


REAL(DP), INTENT(IN)                         :: r
DATA a1/334.85D0/
DATA a2/52.5D0/
DATA a3/1383.0D0/
DATA b1/-1.7624D0/
DATA b2/1.815D0/

varnashort=a1*EXP(b1*r)/r-2.0D0/(1.d0+EXP(b2/r)) *(a2/r**6+a3/r**8)
varnashort=e2*varnashort


RETURN
END FUNCTION varnashort


!------------------------------------------------------------

SUBROUTINE getshortenergy(ityp1,ityp2,ind1,ind2)
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

! sort by particle type
IF (ityp1 > ityp2) THEN
  itmp = ityp1
  ityp1 = ityp2
  ityp2 = itmp
  itmp = ind1
  ind1 = ind2
  ind2 = itmp
  signum = -1D0
!            stop 'warning: signum is -1'
END IF


sumion=0D0

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
IF (isrtyp(ityp1,ityp2) == 1) THEN ! MgO-case
  
  
  
  IF (ityp2 <= 3) THEN ! GSM-GSM case
    itmp1 = iconvshorttolong(ityp1,ind1)
    itmp2 = iconvshorttolong(ityp2,ind2)
    rr = getdistance2(itmp1,itmp2)
    rr = SQRT(MAX(rr,1D-14))

    CALL setbornparas(ityp1,ityp2) ! paras are stored in rVecTmp
    
    sumion = sumion + vbornmayer(rr,rvectmp(1),rvectmp(2), rvectmp(3))
    
  ELSE IF (ityp2 == 4) THEN ! GSM-Na or Na-Na
    
    IF (ityp1 <= 3) THEN
      
      
      
      IF (ityp1 == 1) THEN
        rr = (xc(ind1)-cx(ind2))**2+(yc(ind1)-cy(ind2))**2+  &
            (zc(ind1)-cz(ind2))**2
        rr = MAX(small,SQRT(rr))
        sumion = sumion + bcn/(ecn+EXP((rr-fcn)/sigcn))
        IF (ccn2 /= 0D0) sumion=sumion-  &
            ccn2*(erf(rr/ccncu2)/rr)**2*EXP(-rr/ccncul2)
        IF (ccn3 /= 0D0) sumion=sumion-  &
            ccn3*(erf(rr/ccncu3)/rr)**3*EXP(-rr/ccncul3)
        IF (ccn4 /= 0D0) sumion=sumion-  &
            ccn4*(erf(rr/ccncu4)/rr)**4*EXP(-rr/ccncul4)
        IF (ccn5 /= 0D0) sumion=sumion-  &
            ccn5*(erf(rr/ccncu5)/rr)**5*EXP(-rr/ccncul5)
        IF (ccn6 /= 0D0) sumion=sumion-  &
            ccn6*(erf(rr/ccncu6)/rr)**6*EXP(-rr/ccncul6)
        IF (ccn7 /= 0D0) sumion=sumion-  &
            ccn7*(erf(rr/ccncu7)/rr)**7*EXP(-rr/ccncul7)
        IF (ccn8 /= 0D0) sumion=sumion-  &
            ccn8*(erf(rr/ccncu8)/rr)**8*EXP(-rr/ccncul8)
        IF (ccn9 /= 0D0) sumion=sumion-  &
            ccn9*(erf(rr/ccncu9)/rr)**9*EXP(-rr/ccncul9)
        IF (ccn10 /= 0D0) sumion=sumion-  &
            ccn10*(erf(rr/ccncu10)/rr)**10*EXP(-rr/ccncul10)
        IF (ccnd /= 0D0) sumion=sumion-  &
            ccnd*(erf(rr/ccncud)/rr)**dcn*EXP(-rr/ccnculd)
        
      ELSE IF(ityp1 == 2) THEN
        rr = (xe(ind1)-cx(ind2))**2+(ye(ind1)-cy(ind2))**2+  &
            (ze(ind1)-cz(ind2))**2
      ELSE
        rr = (xk(ind1)-cx(ind2))**2+(yk(ind1)-cy(ind2))**2+  &
            (zk(ind1)-cz(ind2))**2
        rr = MAX(small,SQRT(rr))
!               write(6,*) 'DEBUG: ',rr,ind1
        sumion = sumion + bkn/(ekn+EXP((rr-fkn)/sigkn))
        IF (bkn2 /= 0D0) sumion=sumion+ bkn2/(ekn2+EXP((rr-fkn2)/sigkn2))
        IF (ckn2 /= 0D0) sumion=sumion-  &
            ckn2*(erf(rr/ckncu2)/rr)**2*EXP(-rr/ckncul2)
        IF (ckn3 /= 0D0) sumion=sumion-  &
            ckn3*(erf(rr/ckncu3)/rr)**3*EXP(-rr/ckncul3)
        IF (ckn4 /= 0D0) sumion=sumion-  &
            ckn4*(erf(rr/ckncu4)/rr)**4*EXP(-rr/ckncul4)
        IF (ckn5 /= 0D0) sumion=sumion-  &
            ckn5*(erf(rr/ckncu5)/rr)**5*EXP(-rr/ckncul5)
        IF (ckn6 /= 0D0) sumion=sumion-  &
            ckn6*(erf(rr/ckncu6)/rr)**6*EXP(-rr/ckncul6)
        IF (ckn7 /= 0D0) sumion=sumion-  &
            ckn7*(erf(rr/ckncu7)/rr)**7*EXP(-rr/ckncul7)
        IF (ckn8 /= 0D0) sumion=sumion-  &
            ckn8*(erf(rr/ckncu8)/rr)**8*EXP(-rr/ckncul8)
        IF (ckn9 /= 0D0) sumion=sumion-  &
            ckn9*(erf(rr/ckncu9)/rr)**9*EXP(-rr/ckncul9)
        IF (ckn10 /= 0D0) sumion=sumion-  &
            ckn10*(erf(rr/ckncu10)/rr)**10*EXP(-rr/ckncul10)
        IF (cknd /= 0D0) sumion=sumion-  &
            cknd*(erf(rr/ckncud)/rr)**dkn*EXP(-rr/cknculd)
      END IF ! if ityp1.eq.1
      
      srenergy=srenergy+sumion
      
      
    ELSE ! Na-Na
! no short range interaction in this case
    END IF
    
    
  ELSE IF (ityp2 == 5) THEN ! DFT
    
    IF (ityp1 <= 3) THEN ! DFT-GSM
      
      IF (ityp1 == 1) THEN ! core - DFT
        
      END IF
      
      IF (ityp1 == 3) THEN ! kation - DFT
        
      END IF
      
      
    END IF
    
    
  END IF
  
  
  
  
  
  
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ELSE IF (isrtyp(ityp1,ityp2) == 2) THEN ! Ar case
  
  IF (ityp1 == 1 .AND. ityp2 == 4) THEN
    
    dist2 = (cx(ind2)-xc(ind1))**2+(cy(ind2)-yc(ind1))**2+  &
        (cz(ind2)-zc(ind1))**2
    dist = SQRT(dist2)
    dist = MAX(dist,small)
    
    sumion = varnashort(dist)
    
  ELSE IF (ityp1 == 1  .AND. ityp2 == 1) THEN
    
    
    dist2 =  (xc(ind1)-xc(ind2))**2+(yc(ind1)-yc(ind2))**2+  &
        (zc(ind1)-zc(ind2))**2
    dist = SQRT(dist2)
    dist = MAX(dist,small)
    
    IF(iararlj == 0) THEN
      sumion = vararshort(dist)
    ELSE
      sumion = vararlj(dist)
    END IF
    
  END IF
  
  
ELSE ! other cases
  
! to be implemented
  
END IF

rscaltmp = sumion

RETURN
END SUBROUTINE getshortenergy
!------------------------------------------------------------

!************************************************************

SUBROUTINE addshortrepulsivepot(field,iswitch)
!************************************************************
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
!     gives each GSM particle its short range repulsive (Pauli)
!     potential, i.e. short range interaction GSM-el
!     right now only the argon case is implemented


REAL(DP), INTENT(IN OUT)                     :: field(kdfull2)
INTEGER, INTENT(IN)                      :: iswitch


EXTERNAL varelcore,vfermi,funkfermi

delinv = 1733D0/dx




!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
IF (isrtyp(1,5) == 2) THEN ! The Argon Case; c-DFT
  
  DO is=1,nc
    
    IF (imobc(is) == iswitch) THEN
      
      CALL addtabtofield(field,varelcore0,xc(is),yc(is),zc(is),1D0)
!old          call addFuncToField(field,vArElCore,xc(is),yc(is),zc(is),1D0)
      
    END IF
    
  END DO !nc-loop
  
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ELSE IF (isrtyp(1,5) == 1) THEN ! the MgO-case; c-DFT
  
  
  DO is=1,nc
    
    IF (imobc(is) == iswitch) THEN
      
      CALL addfunctofield3(field,funkfermi,xc(is),yc(is),zc(is),1D0,  &
          fermiac,fermibc,fermicc)
      
      CALL addfunctofield3(field,funkfermi,xc(is),yc(is),zc(is),1D0,  &
          fermia2c,fermib2c,fermic2c)
      
    END IF
    
  END DO ! nc-loop
  
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ELSE IF (isrtyp(1,5) == 3) THEN ! the general fermi function; c-DFT
  DO is=1,nc
    
    IF (imobc(is) == iswitch) THEN
      
      CALL addfunctofield(field,vfermi,xc(is),yc(is),zc(is),1D0)
      
    END IF
    
  END DO
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ELSE ! c-DFT
  
! To be implemented for other cases
  
END IF
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::




!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
IF (isrtyp(2,5) == 1) THEN ! the MgO-case, e-DFT
  
! nothing to do in this model
  
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ELSE IF (isrtyp(2,5) == 2) THEN ! the Ar-case, e-DFT
  
! nothing to do in this model
  
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ELSE ! e-DFT
  
! to be implemented for other cases
  
END IF
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
IF (isrtyp(3,5) == 1) THEN ! the MgO-case, k-DFT
  
  DO is=1,nk
    
    IF (imobk(is) == iswitch) THEN
      
      CALL addfunctofield3(field,funkfermi,xk(is),yk(is),zk(is),1D0,  &
          fermiak,fermibk,fermick)
      
      CALL addfunctofield3(field,funkfermi,xk(is),yk(is),zk(is),1D0,  &
          fermia2k,fermib2k,fermic2k)
      
      
    END IF
    
  END DO ! nk-loop
  
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ELSE IF (isrtyp(3,5) == 2) THEN ! the Ar-case, k-DFT
  
! there are no kations in this model
  
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ELSE IF (isrtyp(3,5) == 3) THEN ! k-DFT by fermi function
  DO is=1,nk
    
    IF (imobk(is) == iswitch) THEN
      CALL addfunctofield(field,vfermi,xk(is),yk(is),zk(is),1D0)
    END IF
  END DO
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ELSE ! k-DFT
  
! to be implemented for other cases
  
END IF
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


RETURN
END SUBROUTINE addshortrepulsivepot
!------------------------------------------------------------

!************************************************************

SUBROUTINE addshortrepulsivepotonsubgrid(field,iswitch)
!************************************************************
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
!     gives each GSM particle its short range repulsive (Pauli)
!     potential, i.e. short range interaction GSM-el
!     right now only the argon case is implemented

!     nsgsize is the length of the subgrid in one direction
!     the total length of the subgrid box is 2*nsgsize+1


REAL(DP), INTENT(IN OUT)                     :: field(kdfull2)
INTEGER, INTENT(IN)                      :: iswitch
INTEGER :: getnearestgridpoint
INTEGER :: conv3to1


EXTERNAL varelcore,vfermi,funkfermi,funkpower



delinv = 1733D0/dx

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
IF (isrtyp(1,5) == 1) THEN ! MgO-case ! c-DFT
  
  DO is=1,nc
    
    IF (imobc(is) == iswitch) THEN
      
      CALL addfunctofieldonsubgrid3(field,funkfermi,  &
          xc(is),yc(is),zc(is),1D0, fermiac,fermibc,fermicc,15)
      
      CALL addfunctofieldonsubgrid3(field,funkfermi,  &
          xc(is),yc(is),zc(is),1D0, fermia2c,fermib2c,fermic2c,15)
      
      
      IF (ccel6 /= 0) THEN
        
        CALL addfunctofieldonsubgrid1(field,funkpower,  &
            xc(is),yc(is),zc(is),ccel6,6)
        
      ELSE IF (ccel8 /= 0) THEN
        
        CALL addfunctofieldonsubgrid1(field,funkpower,  &
            xc(is),yc(is),zc(is),ccel8,8)
        
      ELSE IF (ccel10 /= 0) THEN
        
        CALL addfunctofieldonsubgrid1(field,funkpower,  &
            xc(is),yc(is),zc(is),ccel10,10)
        
      END IF
      
      
    END IF
    
    
    
    
    
  END DO ! nc-loop
  
  
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ELSE IF (isrtyp(1,5) == 2) THEN ! The Argon Case, c-DFT
  
  DO is=1,nc
    
    IF(imobc(is) == iswitch) THEN  ! choose only cores which are
! mobile or unmobile
      
!old          call addFuncToFieldOnSubgrid(field,vArElCore,xc(is),
      CALL addtabtofieldonsubgrid(field,varelcore0,xc(is),  &
          yc(is),zc(is),1D0,10)
      
    END IF
    
  END DO ! nc-loop
  
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ELSE IF (isrtyp(1,5) == 3) THEN ! c-DFT by Fermi function
  
  CALL addfunctofieldonsubgrid(field,vfermi,xc(is), yc(is),zc(is),1D0,10)
  
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ELSE ! c-DFT
! To be implemented for other case
  
END IF
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
IF (isrtyp(2,5) == 1) THEN ! The MgO Case, e-DFT
! nothing to do in this model
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ELSE IF (isrtyp(2,5) == 2) THEN ! The Ar Case, e-DFT
! nothing to do in this model
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ELSE ! e-DFT
! to be implemented for other cases
END IF
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
IF (isrtyp(3,5) == 1) THEN ! The MgO Case, k-DFT
  
  DO is=1,nk
    
    IF (imobk(is) == iswitch) THEN
      CALL addfunctofieldonsubgrid3(field,funkfermi,xk(is)  &
          ,yk(is),zk(is),1D0, fermiak,fermibk,fermick,15)
      
      CALL addfunctofieldonsubgrid3(field,funkfermi,xk(is)  &
          ,yk(is),zk(is),1D0, fermia2k,fermib2k,fermic2k,15)
      
      
      IF (ckel6 /= 0) THEN
        
        CALL addfunctofieldonsubgrid1(field,funkpower,  &
            xk(is),yk(is),zk(is),ckel6,6,15)
        
      ELSE IF (ckel8 /= 0) THEN
        
        CALL addfunctofieldonsubgrid1(field,funkpower,  &
            xk(is),yk(is),zk(is),ckel8,8,15)
        
      ELSE IF (ckel10 /= 0) THEN
        
        CALL addfunctofieldonsubgrid1(field,funkpower,  &
            xk(is),yk(is),zk(is),ckel10,10,15)
        
      END IF
      
      
      
    END IF
    
!          call addFuncToFieldOnSubgrid(potion,vArElCore,xk(is),
!     &       yk(is),zk(is),1D0,10)
  END DO
  
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ELSE IF (isrtyp(3,5) == 2) THEN ! The Ar Case, k-DFT
! nothing to do in this model, no kations present
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ELSE IF (isrtyp(3,5) == 3) THEN ! k-DFT, by Fermi function
  
  CALL addfunctofieldonsubgrid(field,vfermi,xk(is), yk(is),zk(is),1D0,10)
  
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ELSE ! k-DFT
! to be implemented for other cases
END IF
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



!      stop

RETURN
END SUBROUTINE addshortrepulsivepotonsubgrid
!------------------------------------------------------------


!------------------------------------------------------------

SUBROUTINE setbornparas(ityp1,ityp2)
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

!      write(6,*) ityp1,ityp2

IF (ityp1 == 1) THEN ! 1:c
  IF (ityp2 == 1) THEN   ! 2:c
    rvectmp(1) = bcc
    rvectmp(2) = sigcc
    rvectmp(3) = ccc
  ELSE IF (ityp2 == 2) THEN ! 2:e
    rvectmp(1) = bcv
    rvectmp(2) = sigcv
    rvectmp(3) = ccv
  ELSE IF (ityp2 == 3) THEN ! 2:k
    rvectmp(1) = bck
    rvectmp(2) = sigck
    rvectmp(3) = cck
  END IF
ELSE IF(ityp1 == 2) THEN   ! 1:e
  IF (ityp2 == 2) THEN ! 2:e
    rvectmp(1) = bvv
    rvectmp(2) = sigvv
    rvectmp(3) = cvv
  ELSE IF (ityp2 == 3) THEN    ! 2:k
    rvectmp(1) = bkv
    rvectmp(2) = sigkv
    rvectmp(3) = ckv
  END IF
ELSE IF (ityp1 == 3) THEN  ! 1:k
  rvectmp(1) = bkk
  rvectmp(2) = sigkk
  rvectmp(3) = ckk
END IF

RETURN
END SUBROUTINE setbornparas
!------------------------------------------------------------
#endif

!------------------------------------------------------------

FUNCTION funkpower(r,p)
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

funkpower=erf(r/2D0)/(r**p)

RETURN
END FUNCTION funkpower
!------------------------------------------------------------


!------------------------------------------------------------

FUNCTION funkfermi(r,a,b,c)
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

funkfermi=a/(1D0+EXP(b*(r-c)))

RETURN
END FUNCTION funkfermi
!------------------------------------------------------------


!------------------------------------------------------------

FUNCTION funkderfermi(r,a,b,c)
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

funkderfermi=-a*b*EXP(b*(r-c))/(1+EXP(b*(r-c)))**2D0

RETURN
END FUNCTION funkderfermi
!------------------------------------------------------------
