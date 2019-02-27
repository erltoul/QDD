
!-----savings----------------------------------------------------

SUBROUTINE savings(psi,it)

!   Check status of computing resources and optionally save wavefunctions.
!
!   Input:
!      psi    = set of complex s.p. wavefunctions
!      it     = time step in calling routine

USE params
IMPLICIT NONE

COMPLEX(DP), INTENT(IN)   :: psi(kdfull2,kstate)
INTEGER, INTENT(IN)       :: it

REAL(DP) :: time_act, time_elapse



!     if walltime has almost expired, save all relevant data to
!     prepare for restart and stop:

IF(trequest > 0D0) THEN
 IF(myn == 0) THEN                                                        ! cPW
  CALL cpu_time(time_act)
  time_elapse=time_act-time_absinit
  WRITE(6,'(a,5(1pg13.5))') &
   ' cpu_time: trequest,timefrac,t_elapse=',trequest,timefrac,time_elapse
  CALL FLUSH(6)
 END IF
#if(parayes)
 CALL pi_scatter(time_elapse)
#endif
 IF (time_elapse > trequest*timefrac) THEN
  CALL SAVE(psi,it,outnam)
  IF(myn == 0) THEN
    OPEN(660,STATUS='unknown',FILE='progstatus')
    WRITE(660,*) '0'
    CLOSE(660)
  END IF
#if(parayes)
  CALL mpi_barrier (mpi_comm_world, mpi_ierror)
  !CALL mpi_finalize (mpi_ierror)
#endif
  STOP ' finish at TREQUEST'
 END IF
END IF                                                                    ! cPW
!     if program has received user message to shut down program, make it so!
IF (ishutdown /= 0) THEN
  CALL SAVE(psi,it,outnam)
  STOP 'Terminated by user message.'
END IF

RETURN
END SUBROUTINE savings



SUBROUTINE escmask(it)

! Print collected information on escaping electrons.
!
! Input:
!   it   = nr. of time step in calling routine.

USE params
USE util, ONLY: inttostring,printfield
IMPLICIT NONE

INTEGER, INTENT(IN)  :: it

INTEGER :: nbe, nbeabs


IF (jescmaskorb /=0 .AND. MOD(it,jescmaskorb) == 0) THEN
  DO nbe=1,nstate
    nbeabs = nrel2abs(nbe)
    IF (nbeabs < 1000) THEN
      OPEN(588,STATUS='unknown', FILE='pescmaskOrb.'//trim(adjustl(inttostring(nbeabs)))//'.'//outnam)
    ELSE
      STOP 'ERROR: Too many states for jescmaskOrb'
    END IF
    CALL printfield(588,rhoabsoorb(1,nbe),'x')
    CLOSE(588)
  END DO
END IF

IF(myn == 0) THEN
  IF (jescmask .NE. 0 .AND. MOD(it,jescmask) == 0) THEN
    IF(it < 1000000000) THEN
      OPEN(589,STATUS='unknown', FILE='pescmask.'//trim(adjustl(inttostring(it)))//'.'//outnam)
    ELSE
      STOP '::too many time steps::'
    END IF
    CALL printfield(589,rhoabso,'x')
    CLOSE(589)
  END IF


  IF (it == 0) THEN
    OPEN(589,STATUS='unknown',FILE='pescmask.0.'//outnam)    ! Why do it again ?? 
    CALL printfield(589,rhoabso,'x')
    CLOSE(589)
  END IF

END IF

RETURN
END SUBROUTINE escmask


!-----init_scattel-------------------------------------------------

SUBROUTINE init_scattel(psi)

!  Adds one electron in scattering state as Gaussian wavepacket
!  with a certain velocity given by 'scatterelectronv...'.
!  The wavefunctions bookkeeping fields are extended accordingly.


USE params
IMPLICIT NONE

COMPLEX(DP), INTENT(IN OUT)   :: psi(kdfull2,kstate)


INTEGER :: i, ind, ix, iy, iz
REAL(DP) :: arg, fac, fr, pnorm, rr, x1, y1, z1
!------------------------------------------------------------------

pnorm = scatterelectronvxn**2+scatterelectronvyn**2+ scatterelectronvzn**2
pnorm = SQRT(pnorm)

IF (pnorm == 0D0) STOP 'Momentum of Scatt. Electron vanishes!'

scatterelectronvxn = scatterelectronvxn/pnorm
scatterelectronvyn = scatterelectronvyn/pnorm
scatterelectronvzn = scatterelectronvzn/pnorm

pnorm = SQRT(scatterelectronenergy*2D0*ame)

scatterelectronvxn = scatterelectronvxn*pnorm
scatterelectronvyn = scatterelectronvyn*pnorm
scatterelectronvzn = scatterelectronvzn*pnorm

nstate = nstate + 1
IF(nstate > kstate) STOP ' insufficient KSTATE in INIT_SCATTEL'
occup(nstate)=1D0
nrel2abs(nstate)=nstate
nabs2rel(nstate)=nstate
ispin(nstate)=1

! test actual number of active states, only occcupied states allowed
DO i=1,nstate
  IF(occup(i)<0.5D0) STOP "only occupied states allowed in case of attachement"
END DO

fac = 1D0/SQRT(pi**1.5D0*scatterelectronw**3)

ind = 0

DO iz=minz,maxz
  z1=(iz-nzsh)*dz
  DO iy=miny,maxy
    y1=(iy-nysh)*dy
    DO ix=minx,maxx
      x1=(ix-nxsh)*dx
      ind = ind + 1
      arg = (x1-scatterelectronx)*scatterelectronvxn+  &
          (y1-scatterelectrony)*scatterelectronvyn+  &
          (z1-scatterelectronz)*scatterelectronvzn
      rr = (x1 - scatterelectronx)**2 &
          + (y1 - scatterelectrony)**2 &
          + (z1 - scatterelectronz)**2
      fr = fac*EXP(-rr/2D0/scatterelectronw**2)
      psi(ind,nstate) = CMPLX(fr,0D0,DP)*CMPLX(COS(arg),SIN(arg),DP)
    END DO
  END DO
END DO

RETURN
END SUBROUTINE init_scattel


!-----init_projwf-------------------------------------------------

SUBROUTINE init_projwf(psi)

!  Boosts the electronic wavefunctions of the projectile 
!  by the same velocity, whose norm is given by 'eproj' and
!  the direction by 'vpx', 'vpy',and 'vpz'.
!  This routine is similar to 'init_velocel', but applied
!  only to projectile wavefunctions, distinguished by an
!  entry in the array 'proj_states'.

USE params
IMPLICIT NONE

COMPLEX(DP), INTENT(OUT)                     :: psi(kdfull2,kstate)

INTEGER :: ind, ix, iy, iz, kk, nbe, nbee
REAL(DP) :: arg, rnorm, v0, x1, y1, z1
COMPLEX(DP) :: cfac


!------------------------------------------------------------------
! lionel : np(nion) ==> np(nproj)

v0 = SQRT(2D0*eproj/(amu(np(nproj))*1836.0D0*ame))
rnorm = vpx**2 + vpy**2+ vpz**2
rnorm = SQRT(rnorm)

IF (rnorm == 0) STOP 'Velocity vector not normalizable'

vpx = vpx/rnorm*v0*ame
vpy = vpy/rnorm*v0*ame
vpz = vpz/rnorm*v0*ame

IF(taccel>0D0) RETURN               ! boost is done adiabatically

IF(init_lcao.NE.1) THEN
  WRITE(*,*) ' instantaneous acceleration only for INIT_LCAO==1'
  STOP  ' instantaneous acceleration only for INIT_LCAO==1'
END IF
IF(nproj_states==0) THEN
  WRITE(*,*) ' CAUTION : atomic projectile without wf to boost'
  WRITE(*,*) ' if there are electrons on the projectile, please use'
  WRITE(*,*) ' nproj_states in GLOBAL, proj_states and nproj in DYNAMIC'
ELSE
  WRITE(*,*)'Input states of the projectile',proj_states(:)
  ind = 0
  DO iz=minz,maxz
     z1=(iz-nzsh)*dz
     DO iy=miny,maxy
        y1=(iy-nysh)*dy
        DO ix=minx,maxx
           x1=(ix-nxsh)*dx
           ind = ind + 1
           arg = x1*vpx+y1*vpy+z1*vpz
           cfac = CMPLX(COS(arg),SIN(arg),DP)
           DO nbe=1,nstate
              nbee=nrel2abs(nbe)
              DO kk=1,nproj_states
                 IF (nbee == proj_states(kk)) psi(ind,nbe)=cfac*psi(ind,nbe)
              END DO
           END DO
        END DO
     END DO
  END DO
END IF

RETURN
END SUBROUTINE init_projwf


!-----projmoms-------------------------------------------------------projmoms
! REAL version
!-----------------------------------------------------------------------
SUBROUTINE r_projmoms(rho,psi)

! Multipole moments relative to c.m. on 'qetarget' and relative to 
! projectile coordinate 'cz' on 'qeproj'.

USE params
USE util, ONLY:getcm
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
USE util, ONLY:getcm
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

!  print integrated 3d-field (dimension to integrate on depends of values of i3dx, i3dz)

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




SUBROUTINE hpsi_boostinv(qact,aloc,current,rho,nbe)

!     Action of boost-invariant Hamiltonian on one s.p. wavefunction:
!       qact     = wavefunction on which H acts and resulting w.f.
!       aloc     = local potential for the actual spin component
!       current  = average local momentum in x,y, and z-direction
!       nbe      = number of state
!     The routine requires that 'current' has been accumulated before.


USE params
USE kinetic
IMPLICIT NONE

COMPLEX(DP), INTENT(IN OUT)  :: qact(kdfull2)
REAL(DP), INTENT(IN)         :: aloc(2*kdfull2)
REAL(DP), INTENT(IN)         :: current(kdfull2,3)
REAL(DP), INTENT(IN)         :: rho(2*kdfull2)
INTEGER, INTENT(IN)          :: nbe

!                                   workspaces
COMPLEX(DP),ALLOCATABLE :: q1(:),q2(:)

COMPLEX(DP) :: cf



!----------------------------------------------------------------------


! check availability of FFT
IF(.NOT.ALLOCATED(akv)) STOP "HPSI_BOOSTINVARIANT requires FFT"
ALLOCATE(q1(kdfull2),q2(kdfull2))

q1=qact
CALL hpsi(qact,aloc,nbe,0)

CALL xgradient_rspace(q1,q2)
qact = qact - h2m*current(:,1)*q2 
q2 = current(:,1)*q1
CALL xgradient_rspace(q2,q2)
qact = qact - h2m*q2 


CALL ygradient_rspace(q1,q2)
qact = qact - h2m*current(:,2)*q2 
q2 = current(:,2)*q1
CALL ygradient_rspace(q2,q2)
qact = qact - h2m*q2 


CALL zgradient_rspace(q1,q2)
qact = qact - eye*h2m*current(:,3)*q2 
q2 = current(:,3)*q1
CALL zgradient_rspace(q2,q2)
qact = qact - eye* h2m*q2 

qact = qact + h2m* &
 (current(:,1)**2+current(:,2)**2+current(:,3)**2)

WRITE(*,*) ' HPSI_BOOSTINV: E_coll=',dvol*h2m*SUM(rho(:)* &
 (current(:,1)**2+current(:,2)**2+current(:,3)**2))


RETURN
END SUBROUTINE hpsi_boostinv
