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

 
!     ***********************

SUBROUTINE calcpseudo()

!     ***********************
!     Call the appropriate routine for the calculation of pseudopotentials (PsP)

USE params
USE kinetic
IMPLICIT NONE
!      dimension rho(2*kdfull2)
#if(raregas)
INTEGER :: ind
#endif

! Choice of ionic background :
IF(nion2 == 2) THEN         ! read background from a file (potion.dat)
  CALL  pseudo_external()
  RETURN
ELSE IF(nion2 == 0) THEN    ! Jellium (Homogeneous Electron Gas)
  RETURN                    
END IF
! Or else, background from ionic PsP :

SELECT CASE(ipsptyp)
  CASE(0)   ! soft local PsP
    IF(idielec == 0) THEN
      CALL pseudosoft()         ! Psp, nothing special
#if(raregas)
    ELSE
      CALL pseudosoft_dielec()  ! PsP with dielectric support
#endif
    END IF
    
  CASE(1:4)  ! Goedecker PsP (soft or local)
    CALL pseudogoed()       
    
  CASE DEFAULT
    STOP ' CALCPSEUDO: this type of PsP not yet implemented'
END SELECT


#if(raregas)
IF(idielec == 1) THEN
  ! Dielectric layer
  DO ind=1,kdfull2
    rfieldtmp(ind)=potion(ind) ! potion will be erased in pseudosoft2
  END DO
  IF(ipsptyp == 0) THEN
    CALL pseudosoft2()    !!! Not clear what happens here : seems to be similar to pseudosoft_dielec, but outside the layer (and slightly different formula)
  ELSE IF(ipsptyp >= 1) THEN
    STOP ' Goedecker PsP and dielectric layer incompatible'
  END IF
  DO ind=1,kdfull2
    CALL conv1to3(ind)
    IF(iindtmp(1) > nint(xdielec/dx)+nx) THEN 
      potion(ind)=rfieldtmp(ind)  ! Outside the layer, take PsP as it was first computed
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
IMPLICIT NONE

  OPEN(48,FORM='unformatted',FILE='potion.dat')
  READ(48) potion
  CLOSE(48)

RETURN
END SUBROUTINE pseudo_external
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
USE coulsolv, ONLY:solv_poisson_f,solv_poisson_e,tcoulfalr

IMPLICIT NONE

!    size of subgrid in units of mesh size

INTEGER :: ind, ist, ix, iy, iz
REAL(DP) :: cfac1, cfac2, exfac1, exfac2, rr, rx, ry, rz

REAL(DP),DIMENSION(:),ALLOCATABLE :: pseudorho,potsave,potshort

INTEGER,EXTERNAL :: isoutofbox
INTEGER,EXTERNAL :: conv3to1
INTEGER,EXTERNAL :: getnearestgridpoint
REAL(DP),EXTERNAL :: v_soft
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
      
      cfac1 = chg1(np(ist))/(pi**1.5D0*2D0**1.5D0*sgm1(np(ist))**3.D0)
      cfac2 = chg2(np(ist))/(pi**1.5D0*2D0**1.5D0*sgm2(np(ist))**3.D0)
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

IF(tcoulfalr) THEN
  WRITE(*,*) 'in PSEUDOSOFT FALR switch'
  CALL solv_poisson_f(pseudorho,potion,kdfull2)
  WRITE(*,*) 'in PSEUDOSOFT after FALR switch'
ELSE
  WRITE(*,*) 'in PSEUDOSOFT COULEX switch'
  CALL solv_poisson_e(pseudorho,potion,kdfull2)
  WRITE(*,*) 'in PSEUDOSOFT after COULEX switch'
END IF
!  CALL solv_poisson(pseudorho,potion,kdfull2)
  
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


!     ***********************

REAL(DP) FUNCTION dvsdr(r,sigma)

!     ***********************

!     derivation of V_soft
USE params, ONLY: DP,pi
IMPLICIT NONE


REAL(DP), INTENT(IN)  :: r
REAL(DP), INTENT(IN)   :: sigma

REAL(DP), EXTERNAL :: v_soft
REAL(DP):: fac, rabs

rabs = ABS(r)

fac = 2D0/(SQRT(pi)*sigma*rabs)

dvsdr = fac*EXP(-rabs*rabs/(sigma*sigma))-v_soft(rabs,sigma)/rabs

RETURN
END FUNCTION dvsdr

!-------------------------------------------------------------------------


!     ***********************

REAL(DP) FUNCTION d2vsdr2(r,sigma)
USE params, ONLY: DP,pi
IMPLICIT NONE

!     ***********************

!     second derivation of V_soft

REAL(DP), INTENT(IN OUT)                        :: r
REAL(DP), INTENT(IN)                         :: sigma

REAl(DP) :: fac
REAl(DP), EXTERNAL :: dvsdr
REAl(DP), EXTERNAL :: v_soft
r = ABS(r)

fac = -2D0/(SQRT(pi)*sigma)

d2vsdr2 = fac*EXP(-r*r/(sigma*sigma))*(r**(-2D0)+2*sigma**(-2D0))  &
    + r**(-2D0)*v_soft(r,sigma) - dvsdr(r,sigma)/r

RETURN
END FUNCTION d2vsdr2








!-----V_soft------------------------------------------------------------

REAL(DP) FUNCTION v_soft(r,sigma)
USE params, ONLY: DP,PI
IMPLICIT NONE

!     soft Coulomb potential from Gaussian density,
!     uses error function from Chebyshev approximation.
!       r     =  distance at which potential is computed
!       sigma =  width parameter of underlying Gaussian

REAL(DP), INTENT(IN)                         :: r
REAL(DP), INTENT(IN)                         :: sigma

REAl(DP) :: rabs
!------------------------------------------------------------------------
rabs = ABS(r)
v_soft = erf(rabs/sigma)/rabs

RETURN
END FUNCTION v_soft
!-------------------------------------------------------------------------

!old        Archaic pre-f95 way to approach erf(z),  simple precision
!old            Left here as a witness of oldest times. 
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


REAL(DP) FUNCTION v_ion_ion(dist,n1,n2)

USE params
USE kinetic
IMPLICIT NONE

!     effective ion-ion potential
!       dist   =  distance at which potential is computed
!       n1     =  index of the first ion (1<=n1<=nion)
!       n2     =  index of the second ion (1<=n2<=nion)
!       (see array np(*) in 'init.F')

REAL(DP),INTENT(IN) :: dist
INTEGER,INTENT(IN)  :: n1
INTEGER,INTENT(IN)  :: n2

INTEGER :: npmin, npmax
#if(raregas)
REAL(DP),EXTERNAL:: v_soft
REAL(DP),EXTERNAL:: v_ar_ar
REAL(DP),EXTERNAL:: v_ar_na
#endif
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
      v_ion_ion = 0.5D0*c_dipmod*dist*dist
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

REAL(DP) FUNCTION dv_softdr(r,s)
!------------------------------------------------------------
! returns the derivative of erf(r/s)/r by finite differences

USE params, ONLY: DP
IMPLICIT NONE

REAL(DP), INTENT(IN)                     :: r
REAL(DP), INTENT(IN)                     :: s

REAL(DP):: ftemp, rder
REAL(DP),EXTERNAL::v_soft
rder = 1.0D-5

!         ftemp = erf((r+rder)/s)/(r+rder)
!         ftemp = ftemp - erf((r-rder)/s)/(r-rder)
!         ftemp = ftemp/2./rder

ftemp = v_soft(r+rder,s)
ftemp = ftemp - v_soft(r-rder,s)
ftemp = ftemp/(rder+rder)

dv_softdr = ftemp

RETURN
END FUNCTION dv_softdr
