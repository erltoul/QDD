#include "define.h"

#if(raregas)
!------------------------------------------------------------

SUBROUTINE testimage(testpot,testrho)
!------------------------------------------------------------
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)



REAL(DP), INTENT(IN OUT)                     :: testpot(kdfull2)
REAL(DP), INTENT(IN OUT)                     :: testrho(kdfull2)
REAL(DP) :: rho(2*kdfull2),tphi(kdfull2)

EXTERNAL gauss



!      do ii=1,kdfull2/2
!         call  conv1To3(ii)

!         if (iIndTmp(1).lt.nx)
!            ind=conv3To1(iIndTmp(1)+2*(nx-iIndTmp(1)),iIndTmp(2),iIndTmp(3))

!         rho(ii)=rho(ii)-(epsilon-1.)/(epsilon+1.)*rho(ind)



!      enddo








!         ss=1.4
!         fac=1./((PI)**(1.5D0)*(ss)**3.)

!      call addFuncToFieldOnSubgrid1(rho,gauss,25.6,12.8,0.,fac,ss,20)
!      call addFuncToFieldOnSubgrid1(rho,gauss,25.6,-12.8,0.,fac,ss,20)
!         fact=-(epsilon-1)/((PI)**(1.5D0)*((ss)**3.)*(epsilon+1))

!      call addFuncToFieldOnSubgrid1(rho,gauss,-25.6,12.8,0.,fact,ss,20)
!      call addFuncToFieldOnSubgrid1(rho,gauss,-25.6,-12.8,0.,fact,ss,20)
!      write(6,'(a,i)')'KZBOX=',kzbox


CALL printfield(912,testrho,'ptrho')
CALL printfield(913,testpot,'ptpotion')


STOP







RETURN
END SUBROUTINE testimage
!------------------------------------------------------------




!------------------------------------------------------------

SUBROUTINE addimage(rh,ihalfspace)
!------------------------------------------------------------
!     ihalfspace=1: calculates image charges for potential x>xDielec
!     ihalfspace=0: calculates image charges for potential x<xDielec

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)


REAL(DP), INTENT(OUT)                        :: rh(kdfull2)
INTEGER, INTENT(IN)                      :: ihalfspace

INTEGER :: conv3to1



!          ss=1.4
!          fac=1./((PI)**(1.5D0)*(ss)**3.)

!       call addFuncToFieldOnSubgrid1(rho,gauss,0.,0.,0.,fac,ss,20)

IF(ihalfspace == 1) THEN
  DO ii=1,kdfull2
    CALL conv1to3(ii)
    IF (iindtmp(1) < nint(xdielec/dx)+nxsh) THEN
      rh(ii)=rh(ii)*2D0*epsdi/(epsdi+1D0)
      k=-iindtmp(1)+2*(nint(xdielec/dx)+nxsh)
      IF(k <= maxx) THEN
        ind=conv3to1(k,iindtmp(2),iindtmp(3))
        rh(ii)=rh(ii)-(epsdi-1.)/(epsdi+1.)*rh(ind)
!                 write(6,'(4i)') iIndTmp(1),k,ii,ind
      END IF
    END IF
  END DO
ELSE IF(ihalfspace == 0) THEN
  DO ii=1,kdfull2
    CALL conv1to3(ii)
    IF (iindtmp(1) > nint(xdielec/dx)+nxsh) THEN
      rh(ii)=rh(ii)*2D0/(epsdi+1D0)
      k=-iindtmp(1)+2*(nint(xdielec/dx)+nxsh)
      IF(k >= minx) THEN
        ind=conv3to1(k,iindtmp(2),iindtmp(3))
        rh(ii)=rh(ii)+(epsdi-1.)/(epsdi+1.)*rh(ind)
      END IF
    END IF
  END DO
END IF


RETURN
END SUBROUTINE addimage

!------------------------------------------------------------
!------------------------------------------------------------



!-----pseudosoft_dielec-----------------------------------------

SUBROUTINE pseudosoft_dielec()

!     This is the twin brother of 'pseudosoft' which here includes
!     the image potentials in a dielectric layer.
!     The routine compute the pseud-potentials from the cluster ions.

!--------------------------------------------------------------

USE params
!USE kinetic
USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)
!      dimension rho(2*kdfull2)
REAL(DP),ALLOCATABLE :: pseudorho(:)
REAL(DP),ALLOCATABLE :: potsave(:)
REAL(DP),ALLOCATABLE :: potshort(:)
REAL(DP) :: ri(3)
INTEGER :: conv3to1
INTEGER :: getnearestgridpoint

EXTERNAL v_soft

!EQUIVALENCE (pseudorho,w1)
!EQUIVALENCE (potsave,w2)
!EQUIVALENCE (potshort,w3)

!------------------------------------------------------------------

IF(idielec == 0) RETURN
IF (ifreezekspot == 1 .AND. tfs > 0D0) RETURN

ALLOCATE(pseudorho(kdfull2),potsave(kdfull2),potshort(kdfull2))

DO ind=1,nxyz
  potion(ind)=0D0
  potsave(ind)=0D0
  potshort(ind)=0D0
  IF(ipseudo == 1) pseudorho(ind) = 0D0
END DO

IF(ipseudo == 1)THEN
  
!     traditional PsP of the Na cores by pseudodensities:
  
  DO ist=1,nion
    
    IF (isoutofbox(cx(ist),cy(ist),cz(ist)) == 0) THEN
      
      ind = getnearestgridpoint(cx(ist),cy(ist),cz(ist))
      
      CALL conv1to3(ind)
      
      DO iz=iindtmp(3)-nzsg,iindtmp(3)+nzsg
        z1=(iz-nzsh)*dz
        rz=z1-cz(ist)
        DO iy=iindtmp(2)-nysg,iindtmp(2)+nysg
          y1=(iy-nysh)*dy
          ry=y1-cy(ist)
          DO ix=iindtmp(1)-nxsg,iindtmp(1)+nxsg
            x1=(ix-nxsh)*dx
            rx=x1-cx(ist)
            rr=rx*rx+ry*ry+rz*rz
            
            ind = conv3to1(ix,iy,iz)
            
            IF(idielec == 0) THEN
              
              pseudorho(ind) = pseudorho(ind)  &
                  + chg1(np(ist))*EXP(-rr/(2D0*sgm1(np(ist))**2D0))  &
                  /(pi**1.5D0*2D0**1.5D0*sgm1(np(ist))**3.)  &
                  + chg2(np(ist))*EXP(-rr/(2D0*sgm2(np(ist))**2D0))  &
                  /(pi**1.5D0*2D0**1.5D0*sgm2(np(ist))**3.)
              
            ELSE
              IF(cx(ist) < xdielec) THEN
                
                pseudorho(ind) = pseudorho(ind)  &
                    + 2D0*epsdi/(epsdi+1D0)*chg1(np(ist))*  &
                    EXP(-rr/(2D0*sgm1(np(ist))**2D0))  &
                    /(pi**1.5D0*2D0**1.5D0*sgm1(np(ist))**3.)  &
                    + 2D0*epsdi/(epsdi+1D0)*chg2(np(ist))*  &
                    EXP(-rr/(2D0*sgm2(np(ist))**2D0))  &
                    /(pi**1.5D0*2D0**1.5D0*sgm2(np(ist))**3.)
                
              ELSE
                pseudorho(ind) = pseudorho(ind)  &
                    + chg1(np(ist))*EXP(-rr/(2D0*sgm1(np(ist))**2D0))  &
                    /(pi**1.5D0*2D0**1.5D0*sgm1(np(ist))**3.)  &
                    + chg2(np(ist))*EXP(-rr/(2D0*sgm2(np(ist))**2D0))  &
                    /(pi**1.5D0*2D0**1.5D0*sgm2(np(ist))**3.)
              END IF
            END IF
            
          END DO
        END DO
      END DO
      
      IF (idielec /= 0) THEN
        
!     image densities of the Na cores
        
        IF(cx(ist) > xdielec) THEN
          
          IF (isoutofbox(-cx(ist)+2D0*xdielec,cy(ist),cz(ist))  == 0) THEN
            
!     image charge in the box
            
            DO iz=iindtmp(3)-nzsg,iindtmp(3)+nzsg
              z1=(iz-nzsh)*dz
              rz=z1-cz(ist)
              DO iy=iindtmp(2)-nysg,iindtmp(2)+nysg
                y1=(iy-nysh)*dy
                ry=y1-cy(ist)
                DO ix=-iindtmp(1)+2*(nint(xdielec/dx)+nxsh)-nxsg,  &
                      -iindtmp(1)+2*(nint(xdielec/dx)+nxsh)+nxsg
                  x1=(ix-nxsh)*dx
                  rx=x1+cx(ist)-2D0*xdielec
                  rr=rx*rx+ry*ry+rz*rz
                  
                  ind = conv3to1(ix,iy,iz)
                  
                  
                  pseudorho(ind) = pseudorho(ind) - (epsdi-1)/(epsdi+1)*  &
                      (chg1(np(ist))*EXP(-rr/(2D0*sgm1(np(ist))**2D0))  &
                      /(pi**1.5D0*2D0**1.5D0*sgm1(np(ist))**3.)  &
                      + chg2(np(ist))*EXP(-rr/(2D0*sgm2(np(ist))**2D0))  &
                      /(pi**1.5D0*2D0**1.5D0*sgm2(np(ist))**3.))
                  
                  
                END DO
              END DO
            END DO
            
          ELSE ! isOutOfBox(-cx(is)+2D0*xDielec....ne.0)
            
            
            CALL addfunctofield1(potsave,v_soft,-cx(ist)+2D0*xdielec,  &
                cy(ist),cz(ist), -chg1(np(ist))*e2*(epsdi-1D0)/(epsdi+1D0),  &
                sgm1(np(ist))*sq2)
            CALL addfunctofield1(potsave,v_soft,-cx(ist)+2D0*xdielec,  &
                cy(ist),cz(ist), -chg2(np(ist))*e2*(epsdi-1D0)/(epsdi+1D0),  &
                sgm2(np(ist))*sq2)
            
            
          END IF ! isOutOfBox(-cx(is)+2D0*xDielec....ne.0)
          
        END IF ! cx .gt. xDielec
        
      END IF ! iDielec
      
    ELSE ! isOutOfBox not equal 0
      
      IF(idielec == 0) THEN
        
        CALL addfunctofield1(potsave,v_soft,cx(ist),cy(ist),  &
            cz(ist),chg1(np(ist))*e2,sgm1(np(ist))*sq2)
        CALL addfunctofield1(potsave,v_soft,cx(ist),cy(ist),  &
            cz(ist),chg2(np(ist))*e2,sgm2(np(ist))*sq2)
        
      ELSE
        
        IF(cx(ist) < xdielec) THEN
          
          CALL addfunctofield1(potsave,v_soft,cx(ist),cy(ist),  &
              cz(ist),2D0*epsdi/(epsdi+1D0)*chg1(np(ist))*e2, sgm1(np(ist))*sq2)
          CALL addfunctofield1(potsave,v_soft,cx(ist),cy(ist),  &
              cz(ist),2D0*epsdi/(epsdi+1D0)*chg2(np(ist))*e2, sgm2(np(ist))*sq2)
        ELSE
          
          CALL addfunctofield1(potsave,v_soft,cx(ist),cy(ist),  &
              cz(ist),chg1(np(ist))*e2,sgm1(np(ist))*sq2)
          CALL addfunctofield1(potsave,v_soft,cx(ist),cy(ist),  &
              cz(ist),chg2(np(ist))*e2,sgm2(np(ist))*sq2)
          
          CALL addfunctofield1(potsave, v_soft,-cx(ist)+2D0*xdielec,  &
              cy(ist),cz(ist), -chg1(np(ist))*e2*(epsdi-1D0)/(epsdi+1D0),  &
              sgm1(np(ist))*sq2)
          CALL addfunctofield1(potsave, v_soft,-cx(ist)+2D0*xdielec,  &
              cy(ist),cz(ist), -chg2(np(ist))*e2*(epsdi-1D0)/(epsdi+1D0),  &
              sgm2(np(ist))*sq2)
        END IF
        
        
      END IF ! iDielec
      
    END IF ! isOutOfBox
    
  END DO ! sodium core loop
  
  
  IF(isurf /= 0) CALL pseudosoft_substrate(pseudorho,potsave)
  
  
#if(gridfft)
  CALL falr(pseudorho,potion,nx2,ny2,nz2,kdfull2)
#endif
#if(findiff|numerov)
  CALL solv_fft(pseudorho,potion,dx,dy,dz)
#endif
  
  
  IF (isurf /= 0 .AND. ivdw /= 2) THEN
    CALL addshortrepulsivepotonsubgrid(potshort,1)
  ELSE IF (isurf /= 0 .AND. ivdw == 2) THEN
    CALL addshortrepulsivepot(potshort,1)
  END IF
  
  
ELSE ! ipseudo=0
  
  DO ist=1,nion
    
    IF(idielec == 0) THEN
      
      
      CALL addfunctofield1(potion,v_soft,cx(ist),cy(ist),cz(ist),  &
          chg1(np(ist))*e2,sgm1(np(ist))*sq2)
      CALL addfunctofield1(potion,v_soft,cx(ist),cy(ist),cz(ist),  &
          chg2(np(ist))*e2,sgm2(np(ist))*sq2)
      
      
    ELSE
      
      IF(cx(ist) < xdielec) THEN
        
        CALL addfunctofield1(potion,v_soft,cx(ist),cy(ist),cz(ist),  &
            2D0*epsdi/(epsdi+1D0)*chg1(np(ist))*e2,sgm1(np(ist))*sq2)
        CALL addfunctofield1(potion,v_soft,cx(ist),cy(ist),cz(ist),  &
            2D0*epsdi/(epsdi+1D0)*chg2(np(ist))*e2,sgm2(np(ist))*sq2)
        
      ELSE
        
        CALL addfunctofield1(potion,v_soft,cx(ist),cy(ist),cz(ist),  &
            chg1(np(ist))*e2,sgm1(np(ist))*sq2)
        CALL addfunctofield1(potion,v_soft,cx(ist),cy(ist),cz(ist),  &
            chg2(np(ist))*e2,sgm2(np(ist))*sq2)
        
        CALL addfunctofield1(potion,v_soft,-cx(ist)+2D0*xdielec,  &
            cy(ist),cz(ist),-chg1(np(ist))*e2*(epsdi-1D0)/  &
            (epsdi+1D0),sgm1(np(ist))*sq2)
        CALL addfunctofield1(potion,v_soft,-cx(ist)+2D0*xdielec,  &
            cy(ist),cz(ist),-chg2(np(ist))*e2*(epsdi-1D0)/  &
            (epsdi+1D0),sgm2(np(ist))*sq2)
        
      END IF ! cx.lt.xDielec
    END IF ! iDielec.eq.0
  END DO
  
  
  IF (isurf /= 0) THEN
    CALL addgsmpot(potion,1)
    CALL addshortrepulsivepot(potshort,1)
  END IF
  
END IF ! ipseudo

IF(nc > 0 .AND. ivdw == 1) CALL getvdwpot

!     now add contribution from fixed ions and out-of-box ions

DO ind=1,kdfull2
  potion(ind) = potion(ind) + potfixedion(ind) + potsave(ind)+ potshort(ind)
  rfieldtmp(ind)=potshort(ind)
END DO

DEALLOCATE(pseudorho,potsave,potshort)


RETURN
END SUBROUTINE pseudosoft_dielec


!------------------------------------------------------------

SUBROUTINE pseudosoft2()

!------------------------------------------------------------
USE params
!USE kineticx
USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)
!      dimension rho(2*kdfull2)
REAL(DP),ALLOCATABLE :: pseudorho(:)
REAL(DP) :: potsave(kdfull2)
REAL(DP) :: potshort(kdfull2)
REAL(DP) :: ri(3)
INTEGER :: conv3to1
INTEGER :: getnearestgridpoint

EXTERNAL v_soft
ALLOCATE(pseudorho(kdfull2))
!DO ind=1,kdfull2
  potion=0D0
  potsave=0D0
  potshort=0D0
  IF(ipseudo == 1) pseudorho = 0D0
!END DO


IF(ipseudo == 1) THEN
  
!     traditional PsP of the Na cores by pseudodensities:
!     calculates potion for in the dielectric medium
  
  DO is=1,nion
    
    IF (isoutofbox(cx(is),cy(is),cz(is)) == 0) THEN
      
      ind = getnearestgridpoint(cx(is),cy(is),cz(is))
      
      CALL conv1to3(ind)
      
      DO iz=iindtmp(3)-nzsg,iindtmp(3)+nzsg
        z1=(iz-nzsh)*dz
        rz=z1-cz(is)
        DO iy=iindtmp(2)-nysg,iindtmp(2)+nysg
          y1=(iy-nysh)*dy
          ry=y1-cy(is)
          DO ix=iindtmp(1)-nxsg,iindtmp(1)+nxsg
            x1=(ix-nxsh)*dx
            rx=x1-cx(is)
            rr=rx*rx+ry*ry+rz*rz
            
            ind = conv3to1(ix,iy,iz)
            
            IF(idielec == 0) THEN
              
              pseudorho(ind) = pseudorho(ind)  &
                  + chg1(np(is))*EXP(-rr/(2D0*sgm1(np(is))**2D0))  &
                  /(pi**1.5D0*2D0**1.5D0*sgm1(np(is))**3.)  &
                  + chg2(np(is))*EXP(-rr/(2D0*sgm2(np(is))**2D0))  &
                  /(pi**1.5D0*2D0**1.5D0*sgm2(np(is))**3.)
              
            ELSE
              IF(cx(is) > xdielec) THEN
                
                pseudorho(ind) = pseudorho(ind)  &
                    + 2D0/(epsdi+1D0)*chg1(np(is))*  &
                    EXP(-rr/(2D0*sgm1(np(is))**2D0))  &
                    /(pi**1.5D0*2D0**1.5D0*sgm1(np(is))**3.)  &
                    + 2D0/(epsdi+1D0)*chg2(np(is))*  &
                    EXP(-rr/(2D0*sgm2(np(is))**2D0))  &
                    /(pi**1.5D0*2D0**1.5D0*sgm2(np(is))**3.)
                
              ELSE
                pseudorho(ind) = pseudorho(ind)  &
                    + chg1(np(is))*EXP(-rr/(2D0*sgm1(np(is))**2D0))  &
                    /(pi**1.5D0*2D0**1.5D0*sgm1(np(is))**3.)  &
                    + chg2(np(is))*EXP(-rr/(2D0*sgm2(np(is))**2D0))  &
                    /(pi**1.5D0*2D0**1.5D0*sgm2(np(is))**3.)
              END IF
            END IF
            
            
          END DO
        END DO
      END DO
      
      
      
      IF (idielec /= 0) THEN
        
!     image densities of the Na cores
        
        IF(cx(is) < xdielec) THEN
          
          IF (isoutofbox(-cx(is)+2D0*xdielec,cy(is),cz(is)) == 0) THEN
            
!     image charge in the box
            
            DO iz=iindtmp(3)-nzsg,iindtmp(3)+nzsg
              z1=(iz-nzsh)*dz
              rz=z1-cz(is)
              DO iy=iindtmp(2)-nysg,iindtmp(2)+nysg
                y1=(iy-nysh)*dy
                ry=y1-cy(is)
                DO ix=-iindtmp(1)+2*(nint(xdielec/dx)+nxsh)-nxsg,  &
                      -iindtmp(1)+2*(nint(xdielec/dx)+nxsh)+nxsg
                  x1=(ix-nxsh)*dx
                  rx=x1+cx(is)-2D0*xdielec
                  rr=rx*rx+ry*ry+rz*rz
                  
                  ind = conv3to1(ix,iy,iz)
                  
                  
                  pseudorho(ind) = pseudorho(ind) +  &
                      (epsdi-1D0)/(epsdi+1D0)*  &
                      (chg1(np(is))*EXP(-rr/(2D0*sgm1(np(is))**2D0))  &
                      /(pi**1.5D0*2D0**1.5D0*sgm1(np(is))**3.)  &
                      + chg2(np(is))*EXP(-rr/(2D0*sgm2(np(is))**2D0))  &
                      /(pi**1.5D0*2D0**1.5D0*sgm2(np(is))**3.))
                  
                  
                END DO
              END DO
            END DO
            
          ELSE ! isOutOfBox(-cx(is)+2D0*xDielec....ne.0)
            
            
            CALL addfunctofield1(potsave,v_soft,-cx(is)+2D0*xdielec,  &
                cy(is),cz(is), chg1(np(is))*e2*(epsdi-1D0)/(epsdi+1D0),  &
                sgm1(np(is))*sq2)
            CALL addfunctofield1(potsave,v_soft,-cx(is)+2D0*xdielec,  &
                cy(is),cz(is), chg2(np(is))*e2*(epsdi-1)/(epsdi+1),  &
                sgm2(np(is))*sq2)
            
            
          END IF ! isOutOfBox(-cx(is)+2D0*xDielec....ne.0)
          
        END IF ! cx .lt. xDielec
        
      END IF ! iDielec
      
    ELSE  ! isOutOfBox not equal 0
      
      IF(idielec == 0) THEN
        
        CALL addfunctofield1(potsave,v_soft,cx(is),cy(is),  &
            cz(is),chg1(np(is))*e2,sgm1(np(is))*sq2)
        CALL addfunctofield1(potsave,v_soft,cx(is),cy(is),  &
            cz(is),chg2(np(is))*e2,sgm2(np(is))*sq2)
        
      ELSE
        
        IF(cx(is) > xdielec) THEN
          
          CALL addfunctofield1(potsave,v_soft,cx(is),cy(is),  &
              cz(is),2D0/(epsdi+1D0)*chg1(np(is))*e2, sgm1(np(is))*sq2)
          CALL addfunctofield1(potsave,v_soft,cx(is),cy(is),  &
              cz(is),2D0/(epsdi+1D0)*chg2(np(is))*e2, sgm2(np(is))*sq2)
        ELSE
          
          CALL addfunctofield1(potsave,v_soft,cx(is),cy(is),  &
              cz(is),chg1(np(is))*e2,sgm1(np(is))*sq2)
          CALL addfunctofield1(potsave,v_soft,cx(is),cy(is),  &
              cz(is),chg2(np(is))*e2,sgm2(np(is))*sq2)
          
          CALL addfunctofield1(potsave, v_soft,-cx(is)+2D0*xdielec,  &
              cy(is),cz(is), chg1(np(is))*e2*(epsdi-1D0)/(epsdi+1D0),  &
              sgm1(np(is))*sq2)
          CALL addfunctofield1(potsave, v_soft,-cx(is)+2D0*xdielec,  &
              cy(is),cz(is), chg2(np(is))*e2*(epsdi-1D0)/(epsdi+1D0),  &
              sgm2(np(is))*sq2)
        END IF
        
        
      END IF ! iDielec
      
    END IF ! isOutOfBox(cx,cy,cz)
    
  END DO ! sodium core loop
  
!mb         if (isurf.ne.0) call addGSMDensities(pseudorho)
  
  IF (isurf /= 0) THEN
    DO i=1,nc+NE+nk
      CALL getparas(i)
      
      IF (ioutofbox(i) == 0 .AND. imobtmp == 1) THEN
        
        IF(idielec == 0) THEN
          
          CALL addgsmdensity(pseudorho, rvectmp(1),rvectmp(2),  &
              rvectmp(3),sigtmp,chgtmp,0)
          
        ELSE
          
          IF(rvectmp(1) > xdielec) THEN
            
            CALL addgsmdensity(pseudorho, rvectmp(1),rvectmp(2),rvectmp(3),  &
                sigtmp,2D0/(epsdi+1D0)*chgtmp,0)
            
          ELSE
            CALL addgsmdensity(pseudorho, rvectmp(1),rvectmp(2),  &
                rvectmp(3),sigtmp,chgtmp,0)
            
            IF (isoutofbox(-rvectmp(1)+2D0*xdielec,  &
                  rvectmp(2),rvectmp(3)) == 0) THEN
              CALL addgsmdensity(pseudorho,  &
                  -rvectmp(1)+2D0*xdielec,rvectmp(2),rvectmp(3),  &
                  sigtmp,chgtmp*(epsdi-1)/(epsdi+1),0)
            ELSE ! image charge is outside
              CALL addfunctofield1(potsave,v_soft, -rvectmp(1)+2D0*xdielec,  &
                  rvectmp(2),rvectmp(3), e2*chgtmp*(epsdi-1)/(epsdi+1),sigtmp*sq2)
            END IF ! isOutOfBox
          END IF ! rVecTmp(1)
          
        END IF ! iDielec
        
        
      ELSE
        
        IF (imobtmp == 1) THEN ! particle is outside but mobile
          
          IF(idielec == 0) THEN
            
            CALL addfunctofield1(potsave,v_soft,rvectmp(1),  &
                rvectmp(2),rvectmp(3),e2*chgtmp,sigtmp*sq2)
            
          ELSE
            
            IF(rvectmp(1) > xdielec) THEN
              CALL addfunctofield1(potsave,v_soft,  &
                  rvectmp(1),rvectmp(2),rvectmp(3),  &
                  2D0/(epsdi+1D0)*e2*chgtmp,sigtmp*sq2)
            ELSE
              CALL addfunctofield1(potsave,v_soft,  &
                  rvectmp(1),rvectmp(2),rvectmp(3), e2*chgtmp,sigtmp*sq2)
              CALL addfunctofield1(potsave,v_soft, -rvectmp(1)+2D0*xdielec,  &
                  rvectmp(2),rvectmp(3), e2*chgtmp*(epsdi-1)/(epsdi+1),sigtmp*sq2)
            END IF
            
          END IF ! iDielec
          
          
!             write(6,*) 'Warning: particle outside the box, but mobile'
!             write(6,'(3i,3f12.4)')
!     &       i,imobTmp,iOutOfBox(i),rVecTmp(1),rVecTmp(2),rVecTmp(3)
        END IF
        
      END IF
    END DO
  END IF
  
  
#if(gridfft)
  CALL falr(pseudorho,potion,nx2,ny2,nz2,kdfull2)
#endif
#if(findiff|numerov)
  CALL solv_fft(pseudorho,potion,dx,dy,dz)
#endif
  
  IF (isurf /= 0 .AND. ivdw /= 2) THEN
    CALL addshortrepulsivepotonsubgrid(potshort,1)
  ELSE IF (isurf /= 0 .AND. ivdw == 2) THEN
    CALL addshortrepulsivepot(potshort,1)
  END IF
  
ELSE ! ipseudo=0
  
  DO is=1,nion
    
    IF(idielec == 0) THEN
      
      CALL addfunctofield1(potion,v_soft,cx(is),cy(is),cz(is),  &
          chg1(np(is))*e2,sgm1(np(is))*sq2)
      CALL addfunctofield1(potion,v_soft,cx(is),cy(is),cz(is),  &
          chg2(np(is))*e2,sgm2(np(is))*sq2)
      
      
    ELSE
      
      IF(cx(is) > xdielec) THEN
        
        CALL addfunctofield1(potion,v_soft,cx(is),cy(is),cz(is),  &
            2D0/(epsdi+1D0)*chg1(np(is))*e2, sgm1(np(is))*sq2)
        CALL addfunctofield1(potion,v_soft,cx(is),cy(is),cz(is),  &
            2D0/(epsdi+1D0)*chg2(np(is))*e2, sgm2(np(is))*sq2)
      ELSE
        
        CALL addfunctofield1(potion,v_soft,cx(is),cy(is),cz(is),  &
            chg1(np(is))*e2,sgm1(np(is))*sq2)
        CALL addfunctofield1(potion,v_soft,cx(is),cy(is),cz(is),  &
            chg2(np(is))*e2,sgm2(np(is))*sq2)
        
        CALL addfunctofield1(potion,v_soft,-cx(is)+2D0*xdielec,  &
            cy(is),cz(is),chg1(np(is))*e2*(epsdi-1D0)/  &
            (epsdi+1D0),sgm1(np(is))*sq2)
        CALL addfunctofield1(potion,v_soft,-cx(is)+2D0*xdielec,  &
            cy(is),cz(is),chg2(np(is))*e2*(epsdi-1D0)/  &
            (epsdi+1D0),sgm2(np(is))*sq2)
        
      END IF ! cx.gt.xDielec
    END IF ! iDielec.eq.0
    
    
  END DO
  
  
  IF (isurf /= 0) THEN
    CALL addgsmpot2(potion,1)
    CALL addshortrepulsivepot(potshort,1)
  END IF
  
END IF ! ipseudo

IF(nc > 0 .AND. ivdw == 1) CALL getvdwpot

!     now add contribution from fixed ions and out-of-box ions

DO ind=1,kdfull2
  potion(ind) = potion(ind) + potfixedion(ind) + potsave(ind)+ potshort(ind)
  rfieldtmp(ind)=potshort(ind)
  
!cccccccc testBF
!         call conv1To3(ind)
!         if(iIndTmp(1).le.int(xDielec/dx)+nx) then
!            potion(ind)=potion(ind)-2D0
!         endif
!ccccccccccccccccc
  
END DO

!      call testImage(potion,pseudorho)

DEALLOCATE(pseudorho)

RETURN
END SUBROUTINE pseudosoft2


!------------------------------------------------------------

!************************************************************

SUBROUTINE addgsmpot2(field,iswitch)
!************************************************************
USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
!   add electrostatic potentials from GSM particles to potion(kdfull2)
!     iswitch = 0 --> calculate potential of fixed ions
!               1 --> calculate potential of mobile ions



REAL(DP), INTENT(OUT)                        :: field(kdfull2)
INTEGER, INTENT(IN)                      :: iswitch


!      if (nclust.eq.0) return

IF(idielec == 0) THEN
  
  DO is=1,nc
    
    
    
    IF (imobc(is) == iswitch) THEN  ! if .eq. 0 then treat is as static potential
      
      sigc = sigmac
      
      chgtmp = chgc(is)
      
      ind = 0
      DO iz=minz,maxz
        z1=(iz-nzsh)*dz
        rz=z1-zc(is)
        DO iy=miny,maxy
          y1=(iy-nysh)*dy
          ry=y1-yc(is)
          DO ix=minx,maxx
            x1=(ix-nxsh)*dx
            rx=x1-xc(is)
            rr=SQRT(rx*rx+ry*ry+rz*rz)
            rr=MAX(rr,small) ! avoid zero
            ind=ind+1
            
            field(ind)=field(ind) + chgtmp*e2*v_soft(rr,sigc*sq2)
            
            
          END DO
        END DO
      END DO
      
    END IF
    
    
  END DO
  
  DO is=1,NE
    
    
    
    IF (imobe(is) == iswitch) THEN
      
      sigv = sigmav
      chgtmp = chge(is)
      
      
      ind = 0
      DO iz=minz,maxz
        z1=(iz-nzsh)*dz
        rz=z1-ze(is)
        DO iy=miny,maxy
          y1=(iy-nysh)*dy
          ry=y1-ye(is)
          DO ix=minx,maxx
            x1=(ix-nxsh)*dx
            rx=x1-xe(is)
            rr=SQRT(rx*rx+ry*ry+rz*rz)
            rr=MAX(rr,small) ! avoid zero
            ind=ind+1
            
            field(ind)=field(ind) + chgtmp*e2*v_soft(rr,sigv*sq2)
            
            
          END DO
        END DO
      END DO
      
    END IF
    
  END DO
  
  DO is=1,nk
    
    
    
    IF (imobk(is) == iswitch) THEN
      
      sigk = sigmak
      chgtmp = chgk(is)
      
      ind = 0
      DO iz=minz,maxz
        z1=(iz-nzsh)*dz
        rz=z1-zk(is)
        DO iy=miny,maxy
          y1=(iy-nysh)*dy
          ry=y1-yk(is)
          DO ix=minx,maxx
            x1=(ix-nxsh)*dx
            rx=x1-xk(is)
            rr=SQRT(rx*rx+ry*ry+rz*rz)
            rr=MAX(rr,small) ! avoid zero
            ind=ind+1
            
            field(ind)=field(ind) + chgtmp*e2*v_soft(rr,sigk*sq2)
            
            
          END DO
        END DO
      END DO
      
    END IF
    
  END DO
  
END IF ! iDielec.eq.0


IF(idielec /= 0) THEN ! image potential
  
  DO is=1,nc
    
    
    
    IF (imobc(is) == iswitch) THEN  ! if .eq. 0 then treat is as static potential
      sigc = sigmac
      
      chgtmp = chgc(is)
      
      IF(xc(is) < xdielec) THEN
        
        
        ind = 0
        DO iz=minz,maxz
          z1=(iz-nzsh)*dz
          rz=z1-zc(is)
          DO iy=miny,maxy
            y1=(iy-nysh)*dy
            ry=y1-yc(is)
            DO ix=minx,maxx
              x1=(ix-nxsh)*dx
              rx=x1-xc(is)
              rr=SQRT(rx*rx+ry*ry+rz*rz)
              rr=MAX(rr,small) ! avoid zero
              ind=ind+1
              
              field(ind)=field(ind) + chgtmp*e2*v_soft(rr,sigc*sq2)
              
              
            END DO
          END DO
        END DO
        
        ind = 0
        DO iz=minz,maxz
          z1=(iz-nzsh)*dz
          rz=z1-zc(is)
          DO iy=miny,maxy
            y1=(iy-nysh)*dy
            ry=y1-yc(is)
            DO ix=minx,maxx
              x1=(ix-nxsh)*dx
              rx=x1+xc(is)-2D0*xdielec
              rr=SQRT(rx*rx+ry*ry+rz*rz)
              rr=MAX(rr,small) ! avoid zero
              ind=ind+1
              
              field(ind)=field(ind) + (epsdi-1D0)/(epsdi+1D0)*  &
                  chgtmp*e2*v_soft(rr,sigc*sq2)
              
              
            END DO
          END DO
        END DO
        
      ELSE
        
        ind = 0
        DO iz=minz,maxz
          z1=(iz-nzsh)*dz
          rz=z1-zc(is)
          DO iy=miny,maxy
            y1=(iy-nysh)*dy
            ry=y1-yc(is)
            DO ix=minx,maxx
              x1=(ix-nxsh)*dx
              rx=x1-xc(is)
              rr=SQRT(rx*rx+ry*ry+rz*rz)
              rr=MAX(rr,small) ! avoid zero
              ind=ind+1
              
              field(ind)=field(ind) + 2D0/(epsdi+1D0)*  &
                  chgtmp*e2*v_soft(rr,sigc*sq2)
              
              
            END DO
          END DO
        END DO
        
        
      END IF
    END IF
    
  END DO
  
  DO is=1,NE
    
    
    
    IF (imobe(is) == iswitch) THEN
      
      sigv = sigmav
      chgtmp = chge(is)
      
      IF(xe(is) < xdielec) THEN
        
        
        ind = 0
        DO iz=minz,maxz
          z1=(iz-nzsh)*dz
          rz=z1-ze(is)
          DO iy=miny,maxy
            y1=(iy-nysh)*dy
            ry=y1-ye(is)
            DO ix=minx,maxx
              x1=(ix-nxsh)*dx
              rx=x1-xe(is)
              rr=SQRT(rx*rx+ry*ry+rz*rz)
              rr=MAX(rr,small) ! avoid zero
              ind=ind+1
              
              field(ind)=field(ind) + chgtmp*e2*v_soft(rr,sigv*sq2)
              
              
            END DO
          END DO
        END DO
        
        ind = 0
        DO iz=minz,maxz
          z1=(iz-nzsh)*dz
          rz=z1-ze(is)
          DO iy=miny,maxy
            y1=(iy-nysh)*dy
            ry=y1-ye(is)
            DO ix=minx,maxx
              x1=(ix-nxsh)*dx
              rx=x1+xe(is)-2D0*xdielec
              rr=SQRT(rx*rx+ry*ry+rz*rz)
              rr=MAX(rr,small) ! avoid zero
              ind=ind+1
              
              field(ind)=field(ind) + (epsdi-1D0)/(epsdi+1D0)*  &
                  chgtmp*e2*v_soft(rr,sigv*sq2)
              
              
            END DO
          END DO
        END DO
        
      ELSE ! xe(is).gt.xDielec
        
        ind = 0
        DO iz=minz,maxz
          z1=(iz-nzsh)*dz
          rz=z1-ze(is)
          DO iy=miny,maxy
            y1=(iy-nysh)*dy
            ry=y1-ye(is)
            DO ix=minx,maxx
              x1=(ix-nxsh)*dx
              rx=x1-xe(is)
              rr=SQRT(rx*rx+ry*ry+rz*rz)
              rr=MAX(rr,small) ! avoid zero
              ind=ind+1
              
              field(ind)=field(ind) + 2D0/(epsdi+1D0)*  &
                  chgtmp*e2*v_soft(rr,sigv*sq2)
              
              
            END DO
          END DO
        END DO
        
        
      END IF ! xe(is).lt.xDielec
    END IF ! imobe(is).eq.iswitch
    
  END DO
  
  DO is=1,nk
    
    
    
    IF (imobk(is) == iswitch) THEN
      
      sigk = sigmak
      chgtmp = chgk(is)
      
      IF(xk(is) < xdielec) THEN
        
        ind = 0
        DO iz=minz,maxz
          z1=(iz-nzsh)*dz
          rz=z1-zk(is)
          DO iy=miny,maxy
            y1=(iy-nysh)*dy
            ry=y1-yk(is)
            DO ix=minx,maxx
              x1=(ix-nxsh)*dx
              rx=x1-xk(is)
              rr=SQRT(rx*rx+ry*ry+rz*rz)
              rr=MAX(rr,small) ! avoid zero
              ind=ind+1
              
              field(ind)=field(ind) + chgtmp*e2*v_soft(rr,sigk*sq2)
              
              
            END DO
          END DO
        END DO
        
        ind = 0
        DO iz=minz,maxz
          z1=(iz-nzsh)*dz
          rz=z1-zk(is)
          DO iy=miny,maxy
            y1=(iy-nysh)*dy
            ry=y1-yk(is)
            DO ix=minx,maxx
              x1=(ix-nxsh)*dx
              rx=x1+xk(is)-2D0*xdielec
              rr=SQRT(rx*rx+ry*ry+rz*rz)
              rr=MAX(rr,small) ! avoid zero
              ind=ind+1
              
              field(ind)=field(ind) + (epsdi-1D0)/(epsdi+1D0)*  &
                  chgtmp*e2*v_soft(rr,sigk*sq2)
              
              
            END DO
          END DO
        END DO
        
      ELSE ! xk(is).gt.xDielec
        
        ind = 0
        DO iz=minz,maxz
          z1=(iz-nzsh)*dz
          rz=z1-zk(is)
          DO iy=miny,maxy
            y1=(iy-nysh)*dy
            ry=y1-yk(is)
            DO ix=minx,maxx
              x1=(ix-nxsh)*dx
              rx=x1-xk(is)
              rr=SQRT(rx*rx+ry*ry+rz*rz)
              rr=MAX(rr,small) ! avoid zero
              ind=ind+1
              
              field(ind)=field(ind) + 2D0/(epsdi+1D0)*  &
                  chgtmp*e2*v_soft(rr,sigk*sq2)
              
              
            END DO
          END DO
        END DO
        
      END IF ! xk(is).gt.xDielec
    END IF ! imobk(is).eq.iswitch
    
  END DO
  
END IF

!      write(6,*) 'addGSMP: ',field(1)

RETURN
END SUBROUTINE addgsmpot2

!##############################################################

!-----energ_Dielec-----------------------------------------------

SUBROUTINE energ_dielec(rho)

!     Computes Coulomb energy of dielectric background and
!     returns result via common.

USE params
!USE kinetic
USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)


REAL(DP), ALLOCATABLE :: rhotmp(:)
REAL(DP), ALLOCATABLE :: chpcoulimage(:)
!EQUIVALENCE(w1,chpcoulimage)
!EQUIVALENCE(w2,rhotmp)                 ! occupies also w3

!------------------------------------------------------------------

IF(idielec /= 1) THEN
  ecrhoimage=0D0
  RETURN
END IF


!     check workspace

ALLOCATE(rhotmp(2*kdfull2))
ALLOCATE(chpcoulimage(kdfull2))

!     compute dielectric background



DO ii=1,2*kdfull2
  rhotmp(ii)=rho(ii)
END DO

CALL addimage(rho,1)

DO ii=1,2*kdfull2
  CALL conv1to3(ii)
  IF(iindtmp(1) > nint(xdielec/dx)+nxsh) rho(ii)=rho(ii)-rhotmp(ii)
END DO

#if(gridfft)
CALL falr(rho,chpcoulimage,nx2,ny2,nz2,kdfull2)
#endif
#if(findiff|numerov)
CALL solv_fft(rho,chpcoulimage,dx,dy,dz)
#endif

DO ii=1,2*kdfull2
  rho(ii)=rhotmp(ii)
END DO
DO ii=1,kdfull2
  rfieldtmp(ii)=chpcoulimage(ii)
END DO

CALL addimage(rho,0)

DO ii=1,2*kdfull2
  CALL conv1to3(ii)
  IF(iindtmp(1) < nint(xdielec/dx)+nxsh) rho(ii)=rho(ii)-rhotmp(ii)
END DO

#if(gridfft)
CALL falr(rho,chpcoulimage,nx2,ny2,nz2,kdfull2)
#endif
#if(findiff|numerov)
CALL solv_fft(rho,chpcoulimage,dx,dy,dz)
#endif
!test       call prifld(chpcoul,'Dielec.Coul.')


DO ii=1,kdfull2
  CALL conv1to3(ii)
  IF(iindtmp(1) > nint(xdielec/dx)+nxsh)THEN
    chpcoulimage(ii)=rfieldtmp(ii)
  END IF
END DO
DO ii=1,2*kdfull2
  rho(ii)=rhotmp(ii)
END DO



DO ind=1,nxyz
  ecrhoimage=ecrhoimage+rho(ind)*chpcoulimage(ind)
END DO
ecrhoimage=ecrhoimage*dvol/2D0

DEALLOCATE(rhotmp)
DEALLOCATE(chpcoulimage)


RETURN
END SUBROUTINE energ_dielec

#else
!-----energ_Dielec-----------------------------------------------

SUBROUTINE energ_dielec(rho)

!     Computes Coulomb energy of dielectric background and
!     returns result via common.

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)


REAL(DP), INTENT(IN)                     :: rho(2*kdfull2)


!------------------------------------------------------------------

ecrhoimage=0D0
RETURN
END SUBROUTINE energ_dielec
!------------------------------------------------------------

SUBROUTINE addimage(rh,ihalfspace)
!------------------------------------------------------------
!     ihalfspace=1: calculates image charges for potential x>xDielec
!     ihalfspace=0: calculates image charges for potential x<xDielec

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)


REAL(DP), INTENT(IN)                     :: rh(kdfull2)
INTEGER, INTENT(IN OUT)                  :: ihalfspace


STOP ' code not compiled for image potential '

RETURN

END SUBROUTINE addimage
!-----pseudosoft_dielec-----------------------------------------

SUBROUTINE pseudosoft_dielec()
STOP ' code not compiled for dielectric layer'

RETURN
END SUBROUTINE pseudosoft_dielec
!------------------------------------------------------------

SUBROUTINE pseudosoft2()
STOP ' code not compiled for dielectric layer'

RETURN
END SUBROUTINE pseudosoft2

#endif
