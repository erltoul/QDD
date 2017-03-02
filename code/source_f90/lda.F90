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

#include "define.h"

!       ******************************

SUBROUTINE calclocal(rho,aloc)

!       ******************************

!       computes local part of Hamiltonian
!       'time' <= 0  signal static iteration i.e. without laser field


USE params
USE util, ONLY:laserp,projectp
#if(netlib_fft|fftw_cpu)
USE coulsolv
#endif
IMPLICIT NONE

REAL(DP), INTENT(IN)                         :: rho(2*kdfull2)
REAL(DP), INTENT(OUT)                        :: aloc(2*kdfull2)

INTEGER :: ind, jx, jy, jz 
REAL(DP) :: add, addx, addy, addz
REAL(DP),DIMENSION(:),ALLOCATABLE :: rhon,chpdft,vlaser,Vproj

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

!     the Coulombic part
!     warning : counet inserts the esquar factor



#if(gridfft)
IF (nion2 == 0) CALL falr(rhon,chpcoul,kdfull2)
#endif
#if(findiff|numerov)
IF (nion2 == 0) CALL solv_fft(rhon,chpcoul,dx,dy,dz)
#endif
!usew1 = .false.
!test      call prifld(rhon,'Coulomb dens.')
!test      call prifld(chpcoul,'Coulomb pot.')

!chpcoul=0D0

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
    add=chpcoul(ind)-potion(ind)
  ELSE
    add=chpcoul(ind)
  END IF
  IF(tfs > 0D0) THEN
     add = temp + vlaser(ind) + Vproj(ind)
  END IF

  aloc(ind)=chpdft(ind)+add
  aloc(ind+nxyz)=chpdft(ind+nxyz)+add

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
!    with the Gunnarsson & Lundqvist functional 1976


USE params
IMPLICIT NONE
REAL(DP), INTENT(IN)                         :: rho(2*kdfull2)
REAL(DP), INTENT(OUT)                        :: chpdft(2*kdfull2)


#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
INTEGER :: i 
REAL(DP),DIMENSION(:),ALLOCATABLE :: p1,p2
#endif

INTEGER :: ii, mini, maxi, ntranche
REAL(DP) :: e, e1, e2fac, ebase, ec, ec1, enrea1, excf, excp, expot, exfpot , exx, exx1 
REAl(DP) :: fofxi, fpofxi, rlogf, rlogp, rp, rspf, rspif, rspinv, rspip, rspf2, rspp, rspp2
REAl(Dp) :: ubase, upotf, upotp
REAL(DP) :: xi, xip, xip3, xim, xim3
REAL(DP),PARAMETER ::  trd4pi = 3D0/(4D0*pi), onetrd = (1D0/3D0)
!     the Gunnarsson Lundqvist parameters

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

enrear = 0D0
expot  = 4D0*expar/3D0
exfpot = 4D0*exfpar/3D0
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
ntranche=nxyz/knode
mini = (myn)*ntranche+1
maxi =(myn+1)*ntranche

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
  exx    = exx1-0.25D0*rp*(chpdft(ii)*xip+chpdft(ii+nxyz)*xim)
  ec1=exx1+ec1
  ec=exx+ec
END DO
e=ec                    ! this is the rearrangement energy for e_tot
e1=ec1                  ! this is the full functional

!k    the v_xc for spinup/ is in the lower/upper block of chpdft

#if(parayes)
CALL pi_allgather(chpdft,mini,ntranche,p1,kdfull2)
CALL pi_allgather(chpdft(nxyz+1),mini,ntranche,p2,kdfull2)
CALL pi_allreduce(ec,e,1,mpi_double_precision, mpi_sum,mpi_comm_world,icode)
CALL pi_allreduce(ec1,e1,1,mpi_double_precision, mpi_sum,mpi_comm_world,icode)
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
!    with the Perdew-Wang functional
!    attention: the rearrangement energy for e-tot has to be computed
!               the same way as for Gunnarsson & Lundqvist
!               the mean-field part is o.k.!!!!!

USE params
IMPLICIT NONE

REAL(DP), INTENT(IN)                         :: rho(2*kdfull2)
REAL(DP), INTENT(OUT)                        :: chpdft(2*kdfull2)
!~ REAL(DP), SAVE                               :: et
INTEGER :: mysize

#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
INTEGER :: nod
REAL(DP),DIMENSION(:),ALLOCATABLE :: rhonod1,chpdftnod1,rhonod2,chpdftnod2
REAL(DP) :: e, ep

#endif
INTEGER :: ii
REAL(DP) :: ec, rp, xi
REAL(DP) :: a0, a1, a2, a3, da0, da1, da2, da3, da4
REAL(DP) :: b0, b1, b2, b3, b4, db0, db1, db2, db3, db4
REAl(DP) :: t, t1, t2, t3, t4, t5, t6, t7, t8, t10, t11, t12, t13, t15, t17, &
          & t22, t23, t24, t25, t26, t28, t29, t34, t35, t36, t37, t42, t44, t48, &
          & t53, t58, t63, t64, t65, t68, t70, t71, t72, t77, t82,  t83, t88, t93, t98, t102, t109, t135
!        parameter (pi=3.141592654)

! Pade' approximant to Perdew-Wang 92 density functional


!!!!!      icount = 0

#if(!lda_gpu)
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

!write(6,*)rho(1),rho(1+nxyz)
!write(6,*)chpdft(1)

!CALL cpu_time(time_start)

#if(parayes)
CALL mpi_comm_rank(mpi_comm_world,nod,icode)
mysize=lengnod(nod+1)

ALLOCATE(rhonod1(mysize))
ALLOCATE(rhonod2(mysize))
CALL pi_scatterv(rho,nxyz,rhonod1,mysize,icode)
CALL pi_scatterv(rho(nxyz+1),nxyz,rhonod2,mysize,icode)

ALLOCATE(chpdftnod1(mysize))
ALLOCATE(chpdftnod2(mysize))
#else
mysize=nxyz
#endif

!!!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nxyz,rho,chpdft) SCHEDULE(STATIC) REDUCTION(+: ec,enerpw)
!DO ii=1,nxyz
DO ii=1,mysize
#if(parayes)
  rp     = MAX(rhonod1(ii),1D-16)
  xi     = rhonod2(ii)
#else
  rp     = MAX(rho(ii),1D-16)
  xi     = rho(ii+nxyz)
#endif
  
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
  t44 = a0+da0*t13*t17+t23*t29/4.0D0+t35*t37/4D0+3D0/4D0*t42*t3
  t48 = (b1+db1*t13*t17)*t22
  t53 = (b2+db2*t13*t17)*t34
  t58 = (b3+db3*t13*t17)*t26
  t63 = (b4+db4*t13*t17)*t22
  t64 = t36**2
  t = t48*t29/4D0+t53*t37/4D0+3D0/4D0*t58*t3+3D0/16D0*t63*t25*t64
  t68 = 1D0/t
  t70 = t2**2
  t71 = 1D0/t70
  t72 = t1*t71
  t77 = 4D0/3D0*t6*(t3-t72)+4D0/3D0*t10*(-t3+t72)
  t82 = t22*t25
  t83 = t82*t28
  t88 = 1D0/t36*t26*t71
  t93 = t34*t24*t36
  t98 = 1D0/t28*t26*t71
  t102 = t17*t26*t3
  t109 = t**2
  t135 = t44*t68+t2*(da0*t77*t17+da1*t77*t17*t83/4D0-t23*t25*t88/12D0+da2  &
      *t77*t17*t93/4D0-t35*t24*t98/6D0+3D0/4D0*da3*t77*t102-3D0/4D0*t42  &
      *t71)*t68-t2*t44/t109*(db1*t77*t17*t83/4-t48*t25*t88/12D0+db2*t77*t17  &
      *t93/4D0-t53*t24*t98/6D0+3D0/4D0*db3*t77*t102-3D0/4D0*t58*t71+3D0  &
      /16D0*db4*t77*t17*t82*t64-t63*t25*t28*t26*t71/4D0)
  
  
  
  
  
#if(parayes)  
  chpdftnod1(ii)      = -t135  * e2
#else
  chpdft(ii)      = -t135  * e2
#endif
  
  
!      if (neldw.gt.0) then
  
  
  
  
  
  t77 = 4D0/3D0*t6*(-t3-t72)+4D0/3D0*t10*(t3+t72)
  t82 = t22*t25
  t83 = t82*t28
  t88 = 1D0/t36*t26*t71
  t93 = t34*t24*t36
  t98 = 1D0/t28*t26*t71
  t102 = t17*t26*t3
  t109 = t**2
  t135 = t44*t68+t2*(da0*t77*t17+da1*t77*t17*t83/4D0-t23*t25*t88/12+da2  &
      *t77*t17*t93/4D0-t35*t24*t98/6D0+3D0/4D0*da3*t77*t102-3D0/4D0*t42  &
      *t71)*t68-t2*t44/t109*(db1*t77*t17*t83/4D0-t48*t25*t88/12D0+db2*t77*t17 &
      *t93/4D0-t53*t24*t98/6D0+3D0/4D0*db3*t77*t102-3D0/4D0*t58*t71+3D0  &
      /16D0*db4*t77*t17*t82*t64-t63*t25*t28*t26*t71/4D0)
  
  
  
  
  
#if(parayes)  
  chpdftnod2(ii)      = -t135  * e2
#else
  chpdft(ii+nxyz) = -t135  * e2
#endif

!      else
!         chpdft(ii+nxyz) = 0.0
!      endif
  
  
  
  
  t1=rp
  t4 = xi
  t6 = (1D0+t4)**(1D0/3D0)
  t7 = t6**2
  t8 = t7**2
  t10 = (1D0-t4)**(1D0/3D0)
  t11 = t10**2
  t12 = t11**2
  t13 = t8+t12-2D0
  t15 = 2D0**(1D0/3D0)
  t17 = 1D0/(2D0*t15-2D0)
  t22 = 3D0**(1D0/3D0)
  t24 = 4D0**(1D0/3D0)
  t25 = t24**2
  t26 = 1D0/pi
  t28 = (t26*t3)**(1D0/3D0)
  t29 = t25*t28
  t34 = t22**2
  t36 = t28**2
  t37 = t24*t36
  t65 = t36**2
  t70 = t1*(a0+da0*t13*t17+(a1+da1*t13*t17)*t22*t29/4D0+(a2+da2*t13*t17  &
      )*t34*t37/4D0+3D0/4D0*(a3+da3*t13*t17)*t26*t3)/((b1+db1*t13*t17)*  &
      t22*t29/4D0+(b2+db2*t13*t17)*t34*t37/4D0+3D0/4D0*(b3+db3*t13*t17)*t26 &
      *t3+3D0/16D0*(b4+db4*t13*t17)*t22*t25*t65)
  
!    xc energy-density is now:   -e2*t70/rp
  
!    next step is to compose rearrangement energy
  
  t1 = xi*rp
  t2 = rp
  t3=ABS((t1+t2)/2D0)   !  *e2
  t4=ABS((t1-t2)/2D0)   !  *e2
#if(parayes)
  t5= chpdftnod1(ii)*t3+chpdftnod2(ii)*t4
#else
 t5= chpdft(ii)*t3+chpdft(ii+nxyz)*t4
#endif
  
  IF(directenergy) THEN
    enerpw = -t70*e2 + enerpw
  END IF
  
  
  ec = (-t70*e2 - 0.5D0*t5) + ec
!old        ec=-t70/2.0*e2+ec

END DO
!!!!$OMP END PARALLEL DO
#if(parayes)
CALL pi_allgatherv(chpdftnod1,mysize,chpdft,nxyz,icode)
CALL pi_allgatherv(chpdftnod2,mysize,chpdft(nxyz+1),nxyz,icode)
DEALLOCATE(chpdftnod1,chpdftnod2)
DEALLOCATE(rhonod1,rhonod2)
CALL mpi_allreduce(ec,e,1,mpi_double_precision,mpi_sum,mpi_comm_world,icode)
IF(directenergy) &
     CALL mpi_allreduce(enerpw,ep,1,mpi_double_precision,&
     mpi_sum,mpi_comm_world,icode)
ec=e
enerpw=ep
#endif

#else
! lda_gpu
enrear = 0D0
ec=0D0

IF(directenergy) THEN
  enerpw = 0D0
  CALL calc_lda_enerpw_gpu(rho,chpdft,ec,enerpw)
ELSE
  CALL calc_lda_gpu(rho,chpdft,ec)
END IF
#endif

enrear=ec*dvol
!  CALL cpu_time(time_end)
!  write (6,*)ec
!  et=et+time_end-time_start
!  write(6,*)"Time lda:",et
!  STOP


IF(directenergy) THEN
  enerpw = enerpw*dvol
END IF


RETURN
END SUBROUTINE calc_lda_pw92
!#endif




