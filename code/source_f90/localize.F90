#include"define.h"
 
#if(gridfft)
! ---localize-----------------------------------------------------------

#ifdef REALSWITCH
SUBROUTINE localizer(rho,psi)
#else
SUBROUTINE localize(rho,psi)
#endif

!  computes localization criterion of Becke et al
!  dynamical version - not up to date !!!!!!!!!!!!

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN)                     :: rho(2*kdfull2)
#ifdef REALSWITCH
REAL(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)
REAL(DP) :: p(kdfull2)
#else
COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)
COMPLEX(DP) :: p(kdfull2)
#endif
COMPLEX(DP) :: q2(kdfull2)
!REAL(DP),DIMENSION(:),ALLOCATABLE :: akk
COMPLEX(DP) :: test
REAL(DP) :: average_localization(2)
REAL(DP),DIMENSION(:),ALLOCATABLE :: ajt
REAL(DP),DIMENSION(:),ALLOCATABLE :: drho,arho
REAL(DP),DIMENSION(:),ALLOCATABLE :: tau


LOGICAL,PARAMETER :: tupdate=.true.

IF(.NOT.tupdate) STOP ' LOCALIZE not up to date '

!ALLOCATE(akk(kdfull2))
ALLOCATE(ajt(kdfull2))
ALLOCATE(drho(kdfull2))
ALLOCATE(arho(kdfull2))
ALLOCATE(tau(2*kdfull2))

!   init derivative

dkx=pi/(dx*REAL(nx))
dky=pi/(dy*REAL(ny))
dkz=pi/(dz*REAL(nz))



!   scaling factor for Thomas-Fermi energy
tf_fac = 0.6D0*(6.0D0*pi**2)**0.666666666667D0

!   compute the currents, gradient of density and kinetic density

average_localization = 0D0
DO is=2,1,-1
  sign = (3-2*is)*1D0
  
  DO i=1,kdfull2
    tau(i)=0D0
  END DO

  DO idirection=1,3

!         prepare momentum-space factor for actual direction
!    IF(idirection.EQ.1) THEN
!      ind=0
!      DO i3=1,nz2
!        DO i2=1,ny2
!          DO i1=1,nx2
!            IF(i1 >= (nx+1)) THEN
!              zkx=(i1-nx2-1)*dkx
!            ELSE
!              zkx=(i1-1)*dkx
!            END IF
!            ind=ind+1
!            akk(ind)=-zkx
!          END DO
!        END DO
!      END DO
!    ELSEIF(idirection.EQ.2) THEN
!      ind=0
!      DO i3=1,nz2
!        DO i2=1,ny2
!          IF(i2 >= (ny+1)) THEN
!            zky=(i2-ny2-1)*dky
!          ELSE
!            zky=(i2-1)*dky
!          END IF
!          DO i1=1,nx2
!            ind=ind+1
!            akk(ind)=-zky
!          END DO
!        END DO
!      END DO
!    ELSEIF(idirection.EQ.3) THEN
!      ind=0
!      DO i3=1,nz2
!        IF(i3 >= (nz+1)) THEN
!          zkz=(i3-nz2-1)*dkz
!        ELSE
!          zkz=(i3-1)*dkz
!        END IF
!        DO i2=1,ny2
!          DO i1=1,nx2
!            ind=ind+1
!            akk(ind)=-zkz
!          END DO
!        END DO
!      END DO
!    END IF
  
    DO i=1,kdfull2
      ajt(i)=0D0
      drho(i)=0D0
      arho(i)=0D0
    END DO

    DO nb=1,nstate
      IF(ispin(nb).NE.is) CYCLE
      ocfac = occup(nb)*2D0
!    WRITE(*,*) ' nb,is,ocfac=',nb,is,ocfac


#ifdef REALSWITCH
#if(netlib_fft|fftw_cpu)
      CALL rftf(psi(1,nb),q2)
#endif
#if(fftw_gpu)
      CALL rftf(psi(1,nb),q2,ffta,gpu_ffta)
#endif
#else
#if(netlib_fft|fftw_cpu)
      CALL fftf(psi(1,nb),q2)
#endif
#if(fftw_gpu)
      CALL fftf(psi(1,nb),q2,ffta,gpu_ffta)
#endif
#endif
      IF(idirection.EQ.1) THEN
        DO ind=1,kdfull2
#if(netlib_fft|fftw_cpu)
          q2(ind)=q2(ind)*akx(ind)
#endif
#if(fftw_gpu)
          CALL multiply_ak2(gpu_ffta,gpu_akxfft,kdfull2)
#endif
        END DO
      ELSEIF(idirection.EQ.2) THEN
        DO ind=1,kdfull2
#if(netlib_fft|fftw_cpu)
          q2(ind)=q2(ind)*aky(ind)
#endif
#if(fftw_gpu)
          CALL multiply_ak2(gpu_ffta,gpu_akyfft,kdfull2)
#endif
        END DO
      ELSEIF(idirection.EQ.3) THEN
        DO ind=1,kdfull2
#if(netlib_fft|fftw_cpu)
          q2(ind)=q2(ind)*akz(ind)
#endif
#if(fftw_gpu)
          CALL multiply_ak2(gpu_ffta,gpu_akzfft,kdfull2)
#endif
        END DO
      ENDIF
#ifdef REALSWITCH
#if(netlib_fft|fftw_cpu)
      CALL rfftback(q2,p)
#endif
#if(fftw_gpu)
      CALL rfftback(q2,p,ffta,gpu_ffta)
#endif
#else

#if(netlib_fft|fftw_cpu)
      CALL fftback(q2,p)
#endif
#if(fftw_gpu)
      CALL fftback(q2,p,ffta,gpu_ffta)
#endif

#endif
      DO ind=1,kdfull2
#ifdef REALSWITCH
        drho(ind)=drho(ind)+ ocfac*(psi(ind,nb)*p(ind))
        tau(ind)=tau(ind)+ occup(nb)*(p(ind)*p(ind))
        arho(ind)=arho(ind)+ occup(nb)*(psi(ind,nb)*psi(ind,nb))
#else
        drho(ind)=drho(ind)+ ocfac*REAL(CONJG(psi(ind,nb))*p(ind),DP)
        ajt(ind)=ajt(ind)- occup(nb)*AIMAG(CONJG(psi(ind,nb))*p(ind))
        tau(ind)=tau(ind)+ occup(nb)*REAL(CONJG(p(ind))*p(ind),DP)
        arho(ind)=arho(ind)+ occup(nb)*REAL(CONJG(psi(ind,nb))*psi(ind,nb),DP)
#endif
      END DO


  
    END DO  !loop over states

!         accumulate contribution from actual direction on 'tau'
    DO ind=1,kdfull2
!      rp = max(rho(ind)*(0.5D0+sign*rho(ind+kdfull2)),1D-30)
       rp =MAX(arho(ind),1D-90)
      tau(ind)=tau(ind) - (ajt(ind)**2+ 0.25D0*drho(ind)**2)/(rp)
!      arho(ind) = tau(ind)*rp- (ajt(ind)**2+ 0.25D0*drho(ind)**2)
    END DO

!     midx = (minx+maxx)/2
!     midy = (miny+maxy)/2
!     midz = (minz+maxz)/2
!     WRITE(6,'(a,i3)')  '#  test localization, direction=',idirection
!     indadd  = (midz-1)*nxyf+(midy-1)*nyf
!     WRITE(6,'(3(1pg13.5))') ((ix-nxsh)*dx,  &
!         arho(indadd+ix),tau(indadd+ix),ix=minx,maxx)

  
  END DO  ! loop over directions


!     compose localization measure, result stored on 'tau'

  sumpart = 0D0
  DO ind=1,kdfull2
!    rp = max(rho(ind)*(0.5D0+sign*rho(ind+kdfull2)),1D-30)
    rp =MAX(arho(ind),1D-30)
    tau(ind) = 1D0/(1D0+(tau(ind)/(tf_fac*rp**1.666666666667D0))**2)
    IF(tau(ind).LT.1D-99) tau(ind)=0D0
    average_localization(is) = average_localization(is) + tau(ind)*rp
    sumpart = sumpart + rp
  END DO
  average_localization(is) = average_localization(is)/sumpart


  IF(is==2) THEN
    DO i=1,kdfull2
      tau(i+kdfull2)=tau(i)    ! lower spin in upper storage
    END DO
  END IF

END DO  ! loop over spins

!DEALLOCATE(akk)
DEALLOCATE(ajt)
DEALLOCATE(drho)
DEALLOCATE(arho)

!     rescale to final criterion


!     preliminary print along axes


midx = (minx+maxx)/2
midy = (miny+maxy)/2
midz = (minz+maxz)/2

#ifdef REALSWITCH
OPEN(33,STATUS='unknown',FILE='pelfstat.'//outnam)
!CALL safeopen(33,0,1,'pelfstat')
WRITE(33,'(/a)')  '# localization after static iteration'
#else
OPEN(33,STATUS='unknown',FILE='pelf.'//outnam)
!CALL safeopen(33,it,jelf,'pelf')
WRITE(33,'(/a,1f13.5)')  '# localization at time=',tfs
#endif
WRITE(33,'(a,2(1pg13.5))') '# average localizations:',average_localization
WRITE(33,'(a)') &
  '#   coordinate   localiz.-up   localiz.-down  rho-up  rho-down'

WRITE(33,'(/a)')  '#   x    localization '
indadd  = (midz-1)*nxyf+(midy-1)*nyf
WRITE(33,'(5(1pg13.5))') ((ix-nxsh)*dx,  &
    tau(indadd+ix),tau(indadd+ix+kdfull2), &
    rho(indadd+ix)*(0.5D0+rho(indadd+ix+kdfull2)), &
    rho(indadd+ix)*(0.5D0-rho(indadd+ix+kdfull2)), &
    ix=minx,maxx)

WRITE(33,'(/,/a)')  '#   y    localization '
indadd  = (midz-1)*nxyf+midx
WRITE(33,'(5(1pg13.5))') ((iy-nysh)*dy, &
   tau(indadd+(iy-1)*nyf),tau(indadd+(iy-1)*nyf+kdfull2), &
   rho(indadd+(iy-1)*nyf)*(0.5D0+rho(indadd+(iy-1)*nyf+kdfull2)), &
   rho(indadd+(iy-1)*nyf)*(0.5D0-rho(indadd+(iy-1)*nyf+kdfull2)), &
   iy=miny,maxy)

WRITE(33,'(/,/a)')  '#   z    localization '
indadd  = (midy-1)*nyf+midx
WRITE(33,'(5(1pg13.5))') ((iz-nzsh)*dz,  &
   tau(indadd+(iz-1)*nxyf),tau(indadd+(iz-1)*nxyf+kdfull2), &
   rho(indadd+(iz-1)*nxyf)*(0.5D0+rho(indadd+(iz-1)*nxyf+kdfull2)), &
   rho(indadd+(iz-1)*nxyf)*(0.5D0-rho(indadd+(iz-1)*nxyf+kdfull2)), &
   iz=minz,maxz)
WRITE(33,'(/)')
CALL flush(33)
CLOSE(33)

#ifdef REALSWITCH
OPEN(33,STATUS='unknown',FILE='pelfstat2Dxy.'//outnam)
!CALL safeopen(33,0,1,'pelfstat2Dxy')
WRITE(33,'(/a)')  '# localization after static iteration'
#else
CALL safeopen(33,it,jelf,'pelf2Dxy')
WRITE(33,'(/a,1f13.5)')  '# localization at time=',tfs
#endif
ind=0
DO iz=1,maxz
  z1=(iz-nzsh)*dz
  DO iy=1,maxy
    y1=(iy-nysh)*dy
    DO ix=1,maxx
      x1=(ix-nxsh)*dx
      ind=ind+1
      temp=tau(ind)
      IF (z1 == 0D0) WRITE(33,'(3f10.3,3(1pg15.6))') &
            x1,y1,z1,tau(ind),tau(ind+kdfull2),rho(ind)
    END DO
    IF (z1 == 0D0) WRITE(33,'(1x)')
  END DO
END DO
CALL flush(33)
CLOSE(33)

#ifdef REALSWITCH
OPEN(33,STATUS='unknown',FILE='pelfstat2Dyz.'//outnam)
!CALL safeopen(33,0,1,'pelfstat2Dyz')
WRITE(33,'(/a)')  '# localization after static iteration'
#else
CALL safeopen(33,it,jelf,'pelf2Dyz')
WRITE(33,'(/a,1f13.5)')  '# localization at time=',tfs
#endif
ind=0
DO iz=1,maxz
  z1=(iz-nzsh)*dz
  DO iy=1,maxy
    y1=(iy-nysh)*dy
    DO ix=1,maxx
      x1=(ix-nxsh)*dx
      ind=ind+1
      temp=tau(ind)
      IF (x1 == 0D0) WRITE(33,'(3f10.3,3(1pg15.6))') &
            x1,y1,z1,tau(ind),tau(ind+kdfull2),rho(ind)
    END DO
    IF (x1 == 0D0) WRITE(33,'(1x)')
  END DO
  WRITE(33,'(1x)')
END DO
CALL flush(33)
CLOSE(33)

#ifdef REALSWITCH
OPEN(33,STATUS='unknown',FILE='pelfstat2Dxz.'//outnam)
!CALL safeopen(33,0,1,'pelfstat2Dxz')
WRITE(33,'(/a)')  '# localization after static iteration'
#else
CALL safeopen(33,it,jelf,'pelf2Dxz')
WRITE(33,'(/a,1f13.5)')  '# localization at time=',tfs
#endif
ind=0
DO iz=1,maxz
  z1=(iz-nzsh)*dz
  DO iy=1,maxy
    y1=(iy-nysh)*dy
    DO ix=1,maxx
      x1=(ix-nxsh)*dx
      ind=ind+1
      temp=tau(ind)
      IF (y1 == 0D0) WRITE(33,'(3f10.3,3(1pg15.6))') &
            x1,y1,z1,tau(ind),tau(ind+kdfull2),rho(ind)
    END DO
    IF (y1 == 0D0) WRITE(33,'(1x)')
  END DO
END DO

WRITE(33,'(1x)')
WRITE(33,'(1x)')

CALL flush(33)
CLOSE(33)


#ifdef REALSWITCH
OPEN(33,STATUS='unknown',FILE='pelfstat3D.'//outnam)
!CALL safeopen(33,0,1,'pelfstat3D')
WRITE(33,'(/a)')  '# localization after static iteration'
!  no 3D info for dynamic case <--> huge files !!
!CALL safeopen(33,it,jelf,'pelf3D')
!#else
!WRITE(33,'(/a,1f13.5)')  '# localization at time=',tfs
#endif
DO ix=minx,maxx
  DO iy=minz,maxz
    DO iz=minz,maxz
      indadd  = (iy-1)*nyf+ix+(iz-1)*nxyf
      WRITE(33,'(1x,3f10.2,3g13.5)')  &
          (ix-nxsh)*dx,(iy-nysh)*dy,(iz-nzsh)*dz, &
          tau(indadd), tau(indadd+kdfull2), rho(indadd)
    END DO
    WRITE(33,'(1x)') 
  END DO
END DO
CALL flush(33)
CLOSE(33)

DEALLOCATE(tau)

RETURN
#ifdef REALSWITCH
END SUBROUTINE localizer
#else
END SUBROUTINE localize
#endif




#else
!----localize-----------------------------------------------------------

SUBROUTINE localize(rho,psi)

!  dummy

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)


!-----------------------------------------------------------------------

RETURN

END SUBROUTINE localize
!----localizer-----------------------------------------------------------

SUBROUTINE localizer(rho,psi)

!  dummy

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: rho(2*kdfull2)
COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)


!-----------------------------------------------------------------------

RETURN

END SUBROUTINE localizer
#endif



