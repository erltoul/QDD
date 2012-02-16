#include "define.h"
 

!       ******************************

SUBROUTINE calclocal(rho,aloc)

!       ******************************

!       computes local part of hamiltonian
!       'time' <= 0  signal static iteration i.e. without laser field


USE params
!USE kinetic
USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN)                         :: rho(2*kdfull2)
REAL(DP), INTENT(OUT)                        :: aloc(2*kdfull2)

REAL(DP),DIMENSION(:),ALLOCATABLE :: rhon,chpdft,vlaser,homoe,Vproj


IF (ifreezekspot == 1 .AND. tfs > 0D0) RETURN

!     check workspace

!IF(usew1) STOP ' in SSTEP: workspace W1 already active '
!      if(usew4) stop ' in SSTEP: workspace W3 already active '
!usew1 = .true.
!      usew4 = .true.
ALLOCATE(rhon(kdfull2))
ALLOCATE(chpdft(2*kdfull2))
!ALLOCATE(vlaser(kdfull2))



!       first, the netto charge density




!test      write(6,'(a)') 'CALCLOCAL:'

DO ind=1,nxyz
  IF(nion2 == 0) THEN
    rhon(ind)=rho(ind)-rhojel(ind)                 !jellium
  ELSE IF(nion2 /= 0) THEN
    rhon(ind)=rho(ind)                             !pseudo
  END IF
END DO

#if(raregas)
IF(idielec /= 0) THEN
  CALL addimage(rhon,1) ! writes rhon+rhonimage over rhon
END IF
#endif

!     the coulombic part
!     warning : counet inserts the esquar factor



#if(gridfft)
IF (nion2 == 0) CALL falr(rhon,chpcoul,nx2,ny2,nz2,kdfull2)
#endif
#if(findiff|numerov)
IF (nion2 == 0) CALL solv_fft(rhon,chpcoul,dx,dy,dz)
#endif
!usew1 = .false.
!test      call prifld(rhon,'Coulomb dens.')
!test      call prifld(chpcoul,'Coulomb pot.')

!     the lda part

IF(ifsicp /= 5) THEN
  CALL calc_lda(rho,chpdft)
!test        call prifld(chpcoul,'xc potential')
ELSE
  chpdft = 0D0
END IF




!     the laser part


IF(tfs > 0D0) THEN
  ALLOCATE(vlaser(kdfull2))
  CALL laserp(vlaser,rho)
  ALLOCATE(Vproj(kdfull2))
  CALL projectp(Vproj)
END IF


!       the sum



DO ind=1,nxyz
  IF(nion2 /= 0) THEN
    AINT=chpcoul(ind)-potion(ind)
!            aint=-potion(ind)
  ELSE IF(nion2 == 0) THEN
    AINT=chpcoul(ind)
  END IF
  IF(tfs > 0D0) AINT = AINT + vlaser(ind) + Vproj(ind)
  aloc(ind)=chpdft(ind)+AINT
  aloc(ind+nxyz)=chpdft(ind+nxyz)+AINT

#if(raregas)  
  IF (nc > 0 .AND. ivdw == 1) THEN
    aloc(ind) = aloc(ind) + potvdw(ind)
    aloc(ind+nxyz) = aloc(ind+nxyz) + potvdw(ind)
  END IF
#endif

END DO


IF(tfs > 0D0)  THEN
   DEALLOCATE(vlaser)
   DEALLOCATE(Vproj)
ENDIF


!      optionally static dipole potential

IF(tdipolxyz) THEN
  
  ind=0
  DO jz=1,nz2
    addz = (jz-nz)*dz*dpolz
    DO jy=1,ny2
      addy = (jy-ny)*dy*dpoly
      DO jx=1,nx2
        addx = (jx-nx)*dx*dpolx
        add = addx+addy+addz
        ind = ind + 1
        aloc(ind) = aloc(ind) + add
        aloc(ind+nxyz) = aloc(ind+nxyz) + add
      END DO
    END DO
  END DO
END IF

!usew1 = .false.
!      usew4 = .false.
DEALLOCATE(rhon)
DEALLOCATE(chpdft)
!DEALLOCATE(vlaser)


IF (izforcecorr == 1) CALL zeroforce(aloc,rho)


RETURN
END SUBROUTINE calclocal
!#if(gunnar)

!     ******************************

SUBROUTINE calc_lda_gunnar(rho,chpdft)

!     ******************************


!    computes the lsda potential and the rear.energ.
!    with the gunnarsson & lundqvist functional 1976


USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
REAL(DP), INTENT(IN)                         :: rho(2*kdfull2)
REAL(DP), INTENT(OUT)                        :: chpdft(2*kdfull2)


#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif
REAL(DP),DIMENSION(:),ALLOCATABLE :: p1,p2

!     the gunnarson lundqvist parameters

!      data small/1.0e-20/
!#if(exonly)
!DATA cppar,rppar,expar /0D0,11.400D0,0.916D0/
!DATA cfpar,rfpar,exfpar/0D0,15.900D0,1.154D0/
!#else
!DATA cppar,rppar,expar /0.0666D0,11.400D0,0.916D0/
!DATA cfpar,rfpar,exfpar/0.0406D0,15.900D0,1.154D0/
!#endif
REAL(DP) :: cppar=0.0666D0,rppar=11.400D0,expar=0.916D0
REAL(DP) :: cfpar=0.0406D0,rfpar=15.900D0,exfpar=1.154D0

!     check workspace

#if(parayes)
ALLOCATE(p1(kdfull2))
ALLOCATE(p2(kdfull2))
#endif

!     optionally override correlation functional

IF(idenfunc==3) THEN
  cppar = 0D0
  cfpar = 0D0
END IF

!        nxyz=nx2*ny2*nz2
trd4pi = 3D0/(4D0*3.141592653589793D0)
onetrd = (1D0/3D0)
enrear = 0D0
expot  = 4.0*expar/3.0
exfpot = 4.0*exfpar/3.0
e2fac  = e2          ! *0.5
ii     = 0
ec=0D0
ec1=0D0
#if(parayes)
CALL  mpi_comm_rank(mpi_comm_world,myn,icode)
#endif
#if(parano)
myn=0
#endif
nn=myn+1
n=nn-1
ntranche=nxyz/knode
mini = (n)*ntranche+1
maxi =( nn)*ntranche

DO ii=mini,maxi
  rp     = rho(ii)+1D-20
  rspinv = (MAX(small,rp)/trd4pi)**onetrd
  xi     = rho(ii+nxyz)+1D-20
  xip    = 1D0+xi
  IF(xip == 0D0) THEN
    xip3 = 0D0
  ELSE
    xip3   =  (1+xi)**(1D0/3D0)
  END IF
  xim    = 1D0-xi
  IF(xim == 0D0) THEN
    xim3   = 0D0
  ELSE
    xim3   =  (1-xi)**(1D0/3D0)
  END IF
  fofxi  = -1.9236826D0*(xip*xip3+xim*xim3-2D0)
  fpofxi = -2.5648817D0*(xip3-xim3)
  rspip  = rppar*rspinv
  rlogp  = LOG(1D0+rspip)
  upotp  = (-expot*rspinv -cppar*rlogp)
  rspp   = 1D0/rspip
  rspp2  = rspp*rspp
  excp   = -expar*rspinv -cppar*((1D0+rspp2*rspp)*rlogp+0.5D0*rspp-rspp2-onetrd)
  rspif  = rfpar*rspinv
  rlogf  = LOG(1D0+rspif)
  upotf  = (-exfpot*rspinv -cfpar*rlogf)
  rspf   = 1D0/rspif
  rspf2  = rspf*rspf
  excf   = -exfpar*rspinv  &
      -cfpar*((1D0+rspf2*rspf)*rlogf+0.5D0*rspf-rspf2-onetrd)
  ubase  = (1D0+fofxi)*upotp-fofxi*upotf
  ebase  = (excp-excf)*fpofxi
  chpdft(ii) = 0.5D0*(ubase+ebase*xim)*e2fac
!old        chpdfu(ii) = 0.5D0*(ubase-ebase*xip)*e2fac
  chpdft(ii+nxyz) = 0.5D0*(ubase-ebase*xip)*e2fac
  exx1   = rp*(excp+(excp-excf)*fofxi)*0.5D0*e2fac
  exx    = exx1-0.25*rp*(chpdft(ii)*xip+chpdft(ii+nxyz)*xim)
  ec1=exx1+ec1
  ec=exx+ec
END DO
e=ec                    ! this is the rearrangement energy for e_tot
e1=ec1                  ! this is the full functional

!k    the v_xc for spinup/ is in the lower/upper block of chpdft

#if(parayes)
CALL pi_allgather(chpdft,mini,ntranche,p1,kdfull)
CALL pi_allgather(chpdft(nxyz+1),mini,ntranche,p2,kdfull)
CALL pi_allreduce(ec,e,1,mpi_double_precision, mpi_sum,mpi_comm_world,ic)
CALL pi_allreduce(ec1,e1,1,mpi_double_precision, mpi_sum,mpi_comm_world,ic)
DO i=1,nxyz
  chpdft(i)=p1(i)
!mb if number of down spin electrons is 0, then print 0
!        if (neldw.gt.0) then
  chpdft(i+nxyz)=p2(i)
!        else
!           chpdft(i+nxyz)=0.0
!        endif
END DO
#endif

!old#if(parano)
!old      do i=1,nxyz
!old        chpdft(i+nxyz)=chpdfu(i)
!old      enddo
!old#endif
IF(nion == 1) THEN
  enrear=0D0
  enrea1=0D0
ELSE
  enrear=e*dvol     !rearrangement energy for e_tot
  enrea1=e1*dvol    !total xc-energy
END IF

!old      usew1 = .false.
#if(parayes)
DEALLOCATE(p1)
DEALLOCATE(p2)
#endif

RETURN
END SUBROUTINE calc_lda_gunnar
!#endif


!#if(pw92)

!     ******************************

SUBROUTINE calc_lda_pw92(rho,chpdft)

!     ******************************


!    computes the lsda potential and the rear.energ.
!    with the perdew-wang functional
!    attention: the rearrangement energy for e-tot has to be computed
!               the same way as for gunnarsson & lundqvist
!               the mean-field part is o.k.!!!!!

USE params
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN)                         :: rho(2*kdfull2)
REAL(DP), INTENT(OUT)                        :: chpdft(2*kdfull2)

!        parameter (pi=3.141592654)

! pade approximant to perdew-wang 92 density functional


!!!!!      icount = 0


DATA a0  /0.458165293D0 /
DATA da0 /0.119086804D0 /

DATA a1  /2.2170586D0 /
DATA da1 /0.615740256D0 /

DATA a2  /0.740555173D0 /
DATA da2 /0.157420151D0 /

DATA a3  /0.019682278D0 /
DATA da3 /0.003532336D0 /

DATA b1  /1.000000000D0 /
DATA db1 /0.000000000D0 /

DATA b2  /4.504130959D0 /
DATA db2 /0.2361297D0 /

DATA b3  /1.1106363D0 /
DATA db3 /0.205200460D0 /

DATA b4  /0.023592917D0/
DATA db4 /0.004200005D0 /

enrear = 0D0
IF(directenergy) THEN
  enerpw = 0D0
END IF
ec=0D0
DO ii=1,nxyz
  rp     = MAX(rho(ii),1D-16)
  xi     = rho(ii+nxyz)
  
  t1 = xi*rp
  t2 = rp
  t3 = 1D0/t2
  t4 = xi
  t6 = (1D0+t4)**(1D0/3D0)
  t7 = t6**2
  t8 = t7**2
  t10 = (1D0-t4)**(1D0/3D0)
  t11 = t10**2
  t12 = t11**2
  t13 = t8+t12-2
  t15 = 2D0**(1D0/3D0)
  t17 = 1D0/(2D0*t15-2D0)
  t22 = 3D0**(1D0/3D0)
  t23 = (a1+da1*t13*t17)*t22
  t24 = 4D0**(1D0/3D0)
  t25 = t24**2
  t26 = 1/pi
  t28 = (t26*t3)**(1D0/3D0)
  t29 = t25*t28
  t34 = t22**2
  t35 = (a2+da2*t13*t17)*t34
  t36 = t28**2
  t37 = t24*t36
  t42 = (a3+da3*t13*t17)*t26
  t44 = a0+da0*t13*t17+t23*t29/4.0+t35*t37/4+3D0/4D0*t42*t3
  t48 = (b1+db1*t13*t17)*t22
  t53 = (b2+db2*t13*t17)*t34
  t58 = (b3+db3*t13*t17)*t26
  t63 = (b4+db4*t13*t17)*t22
  t64 = t36**2
  t = t48*t29/4+t53*t37/4+3D0/4D0*t58*t3+3D0/16D0*t63*t25*t64
  t68 = 1/t
  t70 = t2**2
  t71 = 1/t70
  t72 = t1*t71
  t77 = 4D0/3D0*t6*(t3-t72)+4D0/3D0*t10*(-t3+t72)
  t82 = t22*t25
  t83 = t82*t28
  t88 = 1/t36*t26*t71
  t93 = t34*t24*t36
  t98 = 1/t28*t26*t71
  t102 = t17*t26*t3
  t109 = t**2
  t135 = t44*t68+t2*(da0*t77*t17+da1*t77*t17*t83/4-t23*t25*t88/12+da2  &
      *t77*t17*t93/4-t35*t24*t98/6+3D0/4D0*da3*t77*t102-3D0/4D0*t42  &
      *t71)*t68-t2*t44/t109*(db1*t77*t17*t83/4-t48*t25*t88/12+db2*t77*t17  &
      *t93/4-t53*t24*t98/6+3D0/4D0*db3*t77*t102-3D0/4D0*t58*t71+3D0  &
      /16D0*db4*t77*t17*t82*t64-t63*t25*t28*t26*t71/4)
  
  
  
  
  
  chpdft(ii)      = -t135  * e2
  
  
!      if (neldw.gt.0) then
  
  
  
  
  
  t77 = 4D0/3D0*t6*(-t3-t72)+4D0/3D0*t10*(t3+t72)
  t82 = t22*t25
  t83 = t82*t28
  t88 = 1/t36*t26*t71
  t93 = t34*t24*t36
  t98 = 1/t28*t26*t71
  t102 = t17*t26*t3
  t109 = t**2
  t135 = t44*t68+t2*(da0*t77*t17+da1*t77*t17*t83/4-t23*t25*t88/12+da2  &
      *t77*t17*t93/4-t35*t24*t98/6+3D0/4D0*da3*t77*t102-3D0/4D0*t42  &
      *t71)*t68-t2*t44/t109*(db1*t77*t17*t83/4-t48*t25*t88/12+db2*t77*t17 &
      *t93/4-t53*t24*t98/6+3D0/4D0*db3*t77*t102-3D0/4D0*t58*t71+3D0  &
      /16D0*db4*t77*t17*t82*t64-t63*t25*t28*t26*t71/4)
  
  
  
  
  
  chpdft(ii+nxyz) = -t135  * e2
  
!      else
!         chpdft(ii+nxyz) = 0.0
!      endif
  
  
  
  
  t1=rp
  t4 = xi
  t6 = (1+t4)**(1D0/3D0)
  t7 = t6**2
  t8 = t7**2
  t10 = (1-t4)**(1D0/3D0)
  t11 = t10**2
  t12 = t11**2
  t13 = t8+t12-2
  t15 = 2**(1D0/3D0)
  t17 = 1/(2*t15-2)
  t22 = 3**(1D0/3D0)
  t24 = 4**(1D0/3D0)
  t25 = t24**2
  t26 = 1/pi
  t28 = (t26*t3)**(1D0/3D0)
  t29 = t25*t28
  t34 = t22**2
  t36 = t28**2
  t37 = t24*t36
  t65 = t36**2
  t70 = t1*(a0+da0*t13*t17+(a1+da1*t13*t17)*t22*t29/4+(a2+da2*t13*t17  &
      )*t34*t37/4+3D0/4D0*(a3+da3*t13*t17)*t26*t3)/((b1+db1*t13*t17)*  &
      t22*t29/4+(b2+db2*t13*t17)*t34*t37/4+3D0/4D0*(b3+db3*t13*t17)*t26 &
      *t3+3D0/16D0*(b4+db4*t13*t17)*t22*t25*t65)
  
!    xc energy-density is now:   -e2*t70/rp
  
!    next step is to compose rearrangement energy
  
  t1 = xi*rp
  t2 = rp
  t3=ABS((t1+t2)/2.0)   !  *e2
  t4=ABS((t1-t2)/2.0)   !  *e2
  t5= chpdft(ii)*t3+chpdft(ii+nxyz)*t4
  
  IF(directenergy) THEN
    enerpw = -t70*e2 + enerpw
  END IF
  
  
  ec = (-t70*e2 - 0.5D0*t5) + ec
!old        ec=-t70/2.0*e2+ec
  
  
  
END DO

enrear=ec*dvol

IF(directenergy) THEN
  enerpw = enerpw*dvol
END IF


RETURN
END SUBROUTINE calc_lda_pw92
!#endif




