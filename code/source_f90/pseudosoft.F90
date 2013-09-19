#include"define.h"
 
!     ***********************

SUBROUTINE calcpseudo()

!     ***********************
USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
!      dimension rho(2*kdfull2)
!REAL(DP) :: fieldtmp(kdfull2)


IF(nion2 == 2) THEN
  call  pseudo_external()
  RETURN
ELSE IF(nion2 == 0) THEN
  RETURN
END IF


IF(ipsptyp == 0) THEN
  IF(idielec == 0) THEN
    CALL pseudosoft()
  ELSE
    CALL pseudosoft_dielec()
  END IF
ELSE IF(ipsptyp >= 1) THEN
  CALL pseudogoed()
ELSE
  STOP ' CALCPSEUDO: this type of PsP not yet implemented'
END IF

#if(raregas)
IF(idielec == 1) THEN
  DO ind=1,kdfull2
    rfieldtmp(ind)=potion(ind)
  END DO
  IF(ipsptyp == 0) THEN
    CALL pseudosoft2()
  ELSE IF(ipsptyp >= 1) THEN
    STOP ' Goedecker PsP and dielectric layer incompatible'
  END IF
  DO ind=1,kdfull2
    CALL conv1to3(ind)
    IF(iindtmp(1) > nint(xdielec/dx)+nx) THEN
      potion(ind)=rfieldtmp(ind)
    END IF
  END DO
END IF
#endif

RETURN
END SUBROUTINE calcpseudo


!------pseudo_external----------------------------------------------

SUBROUTINE pseudo_external()


!     In this routine we read the PsP from a file


USE params
IMPLICIT REAL(DP) (A-H,O-Z)

  OPEN(48,FORM='unformatted',FILE='potion.dat')
  READ(48) potion
  CLOSE(48)

RETURN
END

!------pseudosoft----------------------------------------------

SUBROUTINE pseudosoft()


!     In this routine we calculate ONLY the PsP from
!     the cluster cores.
!     Potentials from substrate ions are included by a call
!     to a separate subroutine. The case of dielectric layer
!     is dealt with in a different routine.

!--------------------------------------------------------------
!     ATTENTION: the definition of Gaussians
!                used in the code:
!                sgm() is the width of a Gaussian defined like
!                      exp(-r**2/(2*sgm**2))

!--------------------------------------------------------------

USE params
USE kinetic
#if(netlib_fft|fftw_cpu)
USE coulsolv, ONLY: falr
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

!    size of subgrid in units of mesh size
!INTEGER,SAVE :: nxsgp=7,nysgp=7,nzsgp=7
!INTEGER,PARAMETER :: nsubgridpsp=7

!      dimension rho(2*kdfull2)
REAL(DP),DIMENSION(:),ALLOCATABLE :: pseudorho,potsave,potshort
!DIMENSION pseudorho(kdfull2)
!DIMENSION potsave(kdfull2)
!DIMENSION potshort(kdfull2)
REAL(DP) :: ri(3)
INTEGER :: conv3to1
INTEGER :: getnearestgridpoint
EXTERNAL v_soft

!EQUIVALENCE (pseudorho,w1)
!EQUIVALENCE (potsave,w2)
!EQUIVALENCE (potshort,w3)

!--------------------------------------------------------------


IF (ifreezekspot == 1 .AND. tfs > 0) RETURN

ALLOCATE(pseudorho(kdfull2))
ALLOCATE(potsave(kdfull2))
ALLOCATE(potshort(kdfull2))

DO ind=1,nxyz
  potion(ind)=0D0
  potsave(ind)=0D0
  potshort(ind)=0D0
  IF(ipseudo == 1) pseudorho(ind) = 0D0
END DO


IF(ipseudo == 1)THEN
  
!     traditional PsP of the Na cores by pseudodensities:
  

  DO ist=1,nion
    
!    radpsp = sgm1(np(ist))/0.8493218
!    nzsgp = nsubgridpsp*NINT(radpsp/dz)
!    nysgp = nsubgridpsp*NINT(radpsp/dy)
!    nxsgp = nsubgridpsp*NINT(radpsp/dx)
    
    IF (isoutofbox(cx(ist),cy(ist),cz(ist)) == 0) THEN
      
      ind = getnearestgridpoint(cx(ist),cy(ist),cz(ist))
      
      CALL conv1to3(ind)
!      WRITE(*,*) ' PsD: ion,ind,iindtmp=',ion,ind,iindtmp
      
      cfac1 = chg1(np(ist))/(pi**1.5*2D0**1.5*sgm1(np(ist))**3.)
      cfac2 = chg2(np(ist))/(pi**1.5*2D0**1.5*sgm2(np(ist))**3.)
      exfac1 = 1D0/(2D0*sgm1(np(ist))**2D0)
      exfac2 = 1D0/(2D0*sgm2(np(ist))**2D0)
      DO iz=iindtmp(3)-nzsg,iindtmp(3)+nzsg
        rz=(iz-nzsh)*dz-cz(ist)
        DO iy=iindtmp(2)-nysg,iindtmp(2)+nysg
          ry=(iy-nysh)*dy-cy(ist)
          DO ix=iindtmp(1)-nxsg,iindtmp(1)+nxsg
            rx=(ix-nxsh)*dx-cx(ist)
            rr=rx*rx+ry*ry+rz*rz
            
            ind = conv3to1(ix,iy,iz)
            
            pseudorho(ind) = pseudorho(ind)  &
                + cfac1*EXP(-rr*exfac1) + cfac2*EXP(-rr*exfac2)
            
          END DO
        END DO
      END DO
!  CALL prifld(pseudorho,'pseudo-dens')
      
      
    ELSE ! isOutOfBox not equal 0
      
      CALL addfunctofield1(potsave,v_soft,cx(ist),cy(ist),  &
          cz(ist),chg1(np(ist))*e2,sgm1(np(ist))*sq2)
      CALL addfunctofield1(potsave,v_soft,cx(ist),cy(ist),  &
          cz(ist),chg2(np(ist))*e2,sgm2(np(ist))*sq2)
      
      
    END IF ! isOutOfBox
    
  END DO ! sodium core loop
  
  
#if(raregas)  
  IF(isurf /= 0) CALL pseudosoft_substrate(pseudorho,potsave)
#endif
  
  
#if(gridfft)
  CALL falr(pseudorho,potion,nx2,ny2,nz2,kdfull2)
#endif
#if(findiff|numerov)
  CALL solv_fft(pseudorho,potion,dx,dy,dz)
#endif
  
  
#if(raregas)
  IF (isurf /= 0 .AND. ivdw /= 2) THEN
    CALL addshortrepulsivepotonsubgrid(potshort,1)
  ELSE IF (isurf /= 0 .AND. ivdw == 2) THEN
    CALL addshortrepulsivepot(potshort,1)
  END IF
#endif
  
  
ELSE ! ipseudo=0
  
!        pseudo-potentials directly
  
  DO ist=1,nion
    
    CALL addfunctofield1(potion,v_soft,cx(ist),cy(ist),cz(ist),  &
        chg1(np(ist))*e2,sgm1(np(ist))*sq2)
    CALL addfunctofield1(potion,v_soft,cx(ist),cy(ist),cz(ist),  &
        chg2(np(ist))*e2,sgm2(np(ist))*sq2)
    
  END DO
  
#if(raregas)
  IF (isurf /= 0) THEN
    CALL addgsmpot(potion,1)
    CALL addshortrepulsivepot(potshort,1)
  END IF
#endif
  
END IF ! ipseudo

#if(raregas)
IF(nc > 0 .AND. ivdw == 1) CALL getvdwpot
#endif

!     now add contribution from fixed ions and out-of-box ions

DO ind=1,kdfull2
  potion(ind) = potion(ind) + potfixedion(ind) + potsave(ind)+ potshort(ind)
#if(raregas)
  rfieldtmp(ind)=potshort(ind)
#endif
END DO

DEALLOCATE(pseudorho)
DEALLOCATE(potsave)
DEALLOCATE(potshort)


RETURN
END SUBROUTINE pseudosoft

!GB


!GB

!     ***********************

FUNCTION dvsdr(r,sigma)

!     ***********************

!     derivation of V_soft
USE params, ONLY: DP,pi
IMPLICIT REAL(DP) (A-H,O-Z)


REAL(DP), INTENT(OUT)                        :: r
REAL(DP), INTENT(IN)                         :: sigma
!DATA pi/3.141592653589793/

r = ABS(r)

fac = 2D0/(SQRT(pi)*sigma*r)

dvsdr = fac*EXP(-r*r/(sigma*sigma))-v_soft(r,sigma)/r

RETURN
END FUNCTION dvsdr

!-------------------------------------------------------------------------


!     ***********************

FUNCTION d2vsdr2(r,sigma)
USE params, ONLY: DP,pi
IMPLICIT REAL(DP) (A-H,O-Z)

!     ***********************

!     second derivation of V_soft



REAL(DP), INTENT(OUT)                        :: r
REAL(DP), INTENT(IN)                         :: sigma
!DATA pi/3.141592653589793/

r = ABS(r)

fac = -2D0/(SQRT(pi)*sigma)

d2vsdr2 = fac*EXP(-r*r/(sigma*sigma))*(r**(-2D0)+2*sigma**(-2D0))  &
    + r**(-2D0)*v_soft(r,sigma) - dvsdr(r,sigma)/r

RETURN
END FUNCTION d2vsdr2
!GB







!-----V_soft------------------------------------------------------------

FUNCTION v_soft(r,sigma)
USE params, ONLY: DP
IMPLICIT REAL(DP) (A-H,O-Z)

!     soft Coulomb potential from Gaussian density,
!     uses error function from Chebyshev approximation.
!       r     =  distance at which potential is computed
!       sigma =  width parameter of underlying Gaussian


REAL(DP), INTENT(IN)                         :: r
REAL(DP), INTENT(IN)                         :: sigma
REAL(DP), PARAMETER :: pi=3.141592653589793D0

!------------------------------------------------------------------------
rabs = ABS(r)

!                       the error function  (good for 10**-7 precision)

z=rabs/sigma


!     use Coulomb cut-off for V_soft
!     relative error is smaller than 1e-20

IF (z > 6D0) THEN
  v_soft = 1D0/rabs
  RETURN
END IF

IF(z <= 1D-1)THEN             ! use Taylor expansion for z < 0.1
  v_soft=(2D0-2D0*z**2/3D0+z**4/5D0-z**6/21D0) /(SQRT(pi)*sigma)
ELSE
  t=1D0/(1D0+0.5D0*z)
  erfcc=t*EXP(-z*z-1.26551223D0+t*(1.00002368D0+t*(.37409196D0+t*  &
      (.09678418D0+t*(-.18628806D0+t*(.27886807D0+t*(-1.13520398D0+t*  &
      (1.48851587D0+t*(-.82215223D0+t*.17087277D0)))))))))
  IF (r < 0D0) erfcc=2D0-erfcc
  f=1D0-erfcc
  
!     final composition
  v_soft = f/rabs
END IF

RETURN
END FUNCTION v_soft
!-------------------------------------------------------------------------



!old c
!old c
!old c       ******************************
!old c
!old         function f(x)
!old c
!old c       ******************************
!old c
!old       z=abs(x)
!old       t=1./(1.+0.5*z)
!old       erfcc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+t*
!old      *(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+t*
!old      *(1.48851587+t*(-.82215223+t*.17087277)))))))))
!old       if (x.lt.0.) erfcc=2.-erfcc
!old       f=1.0-erfcc
!old       return
!old       end


!-----V_ion_ion------------------------------------------------------------


FUNCTION v_ion_ion(dist,n1,n2)

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

!     effective ion-ion potential
!       r      =  distance at which potential is computed
!       n1     =  index of the first ion (1<=n1<=nion)
!       n2     =  index of the second ion (1<=n2<=nion)
!       (see array np(*) in 'init.F')


!------------------------------------------------------------------------

npmin = MIN(np(n1),np(n2))
npmax = MAX(np(n1),np(n2))
#if(raregas)
IF(npmin == -18) THEN
  IF(npmax == -18) THEN
    v_ion_ion = ch(-18)**2*e2*v_soft(dist,2D0*sgm1(-18))
!g          V_ion_ion = e2*ch(-18)**2*V_soft(dist,sgm1(-18)*SQ2)
  ELSE IF(npmax == 18) THEN
    IF(ABS(n1-n2) == nrare)THEN
      v_ion_ion = 0.5*c_dipmod*dist*dist
    ELSE
      v_ion_ion = ch(18)*ch(-18)*e2*v_soft(dist,2D0*sgm1(18))
!g              V_ion_ion = e2*ch(18)*ch(-18)*V_soft(dist,sgm1(18)*SQ2)
    END IF
  ELSE
!test          V_ion_ion =-V_Ar_Na(dist)
!g          V_ion_ion = ch(-18)*ch(npmax)*e2*V_soft(dist,sgm1(18))
    v_ion_ion = e2*ch(-18)*ch(npmax)*v_soft(dist,sgm1(18)*sq2)
  END IF
ELSE IF(npmin /= 18) THEN
  IF(npmax == 18) THEN
!g          V_ion_ion = ch(18)*ch(npmin)*e2*V_soft(dist,sgm1(18))
    v_ion_ion = e2*ch(18)*ch(npmin)*v_soft(dist,sgm1(18)*sq2) + v_ar_na(dist)
  ELSE
    v_ion_ion=e2*ch(npmin)*ch(npmax)/dist
  END IF
ELSE
  v_ion_ion = v_ar_ar(dist)
END IF
#else
v_ion_ion=e2*ch(npmin)*ch(npmax)/dist
#endif


RETURN
END FUNCTION v_ion_ion


!------------------------------------------------------------

FUNCTION dv_softdr(r,s)
!------------------------------------------------------------
! returns the derivative of erf(r/s)/r by finite differences

!      double precision r,s


USE params, ONLY: DP
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: r
REAL(DP), INTENT(IN OUT)                     :: s
DOUBLE PRECISION :: rder

rder = 1.0E-5

!         ftemp = erf((r+rder)/s)/(r+rder)
!         ftemp = ftemp - erf((r-rder)/s)/(r-rder)
!         ftemp = ftemp/2./rder

ftemp = v_soft(r+rder,s)
ftemp = ftemp - v_soft(r-rder,s)
ftemp = ftemp/(rder+rder)

dv_softdr = ftemp

RETURN
END FUNCTION dv_softdr
!------------------------------------------------------------

!------------------------------------------------------------

FUNCTION v_coulomb(r)
!------------------------------------------------------------
USE params, ONLY: DP
IMPLICIT REAL(DP) (A-H,O-Z)


r = MAX(small,r)

v_coulomb = 1/r

RETURN
END FUNCTION v_coulomb
!------------------------------------------------------------
