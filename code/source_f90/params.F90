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

MODULE params
IMPLICIT NONE
SAVE
INTEGER,PARAMETER :: DP=KIND(1D0)  ! precision  setting

!
! general settings for dimensions, part.numbers etc
!

! number of nodes (=1 for serial version)
INTEGER :: knode=2
! max. nr. electron states per node
!fix! INTEGER,PARAMETER :: kstate=20
INTEGER :: kstate=0
! max. total nr. electron  states
INTEGER :: ksttot
INTEGER,PRIVATE :: ksttot2

!  settings ad definitions for openmp parallel computing
#if(paropenmp)
INTEGER :: numthr = 4  ! actual number of threads in openmp
INTEGER,EXTERNAL :: OMP_GET_MAX_THREADS, OMP_GET_NUM_PROCS, OMP_NUM_THREADS
INTEGER,EXTERNAL :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
EXTERNAL :: OMP_SET_NUM_THREADS
#else
INTEGER,PARAMETER :: numthr = 1  ! actual number of threads in openmp
#endif
INTEGER :: nthr                ! max number of threads -- 1

! maximum number of ions
!fix! INTEGER,PARAMETER :: ng=8
INTEGER :: ng


INTEGER,PARAMETER :: maxnang=100    ! max. angular bins for PAD
INTEGER,PARAMETER :: maxmps=180    ! max. nr. analyzing points for PES
! parameter(drtab=1e-4)


! physical units and constants (here Rydberg units)
REAL(DP),PARAMETER :: e2=2.0D0
REAL(DP),PARAMETER :: hbar=1.0D0
REAL(DP),PARAMETER :: ame=0.5D0
REAL(DP),PARAMETER :: h2m=hbar*hbar/2.0/ame  !  h2m =  hbar**2/2m

! frequently used mathematical constants
REAL(DP),PARAMETER :: zero=0.0D0
REAL(DP),PARAMETER :: half=0.5D0
REAL(DP),PARAMETER :: one=1.0D0
REAL(DP),PARAMETER :: PI=3.141592653589793D0
REAL(DP),PARAMETER :: fourpi=4.0D0*PI
REAL(DP),PARAMETER :: small=1.0D-20
REAL(DP),PARAMETER :: sq2=1.4142135623730950D0
REAL(DP),PARAMETER :: sqrtpi=1.772453850905516D0
COMPLEX(DP),PARAMETER :: eye=(0D0,1D0)


!     maximum sizes of the box in x,y, and z
INTEGER :: kxbox=0,kybox=0,kzbox=0
! deduced grid parameters
INTEGER :: kxmax,kymax,kzmax
INTEGER :: nx2,ny2,nz2
!INTEGER,PARAMETER :: nxyzf=nx2*ny2*nz2                 !?
INTEGER :: nxyz,nyf,nxyf
INTEGER :: kdfull2
!INTEGER,PARAMETER :: kfbox=kdfull2                    ! ??
INTEGER :: nx,ny,nz
!INTEGER,PARAMETER :: nx1=nx+1,ny1=ny+1,nz1=nz+1,nzi=nz1,nzr=nz+nz  !?
!INTEGER,PARAMETER :: nxy1=nx1*ny1    !?
!INTEGER,PARAMETER :: ksmax=nx+1,kdfull=(nx+1)*(ny+1)*(nz+1)    !?
INTEGER :: knodem     ! =knode-1
#if(gridfft)
! bounds of loops
INTEGER :: minx=1,maxx
INTEGER :: miny=1,maxy
INTEGER :: minz=1,maxz
! bounds of esc. el. ???
INTEGER :: nbnx=2,nbxx
INTEGER :: nbny=2,nbxy
INTEGER :: nbnz=2,nbxz
! mid-point for x,y,z-values
INTEGER :: nxsh,nysh,nzsh
#endif
#if(findiff|numerov)
! bounds of loops
INTEGER :: minx,maxx
INTEGER :: miny,maxy
INTEGER :: minz,maxz
! bounds of esc. el.
INTEGER :: nbnx,nbxx
INTEGER :: nbny,nbxy
INTEGER :: nbnz,nbxz
! offset for x,y,z-values   ???
INTEGER :: nxsh=0,nysh=0,nzsh=0
#endif


REAL(DP) :: dx=0D0,dy=0D0,dz=0D0,dvol                  !  mesh spacing, volume
REAL(DP),ALLOCATABLE :: xval(:),yval(:),zval(:)        !  grid coordinates
REAL(DP),ALLOCATABLE :: xt2(:),yt2(:),zt2(:)           !  coordinates**2


REAL(DP),ALLOCATABLE :: enonlo(:)               !  non-local s.p.energy
REAL(DP),ALLOCATABLE :: amoy(:)                 !  single particle energies
REAL(DP),ALLOCATABLE :: spvariance(:)           !  s.p. energy variances
REAL(DP),ALLOCATABLE :: spvariancep(:)          !  s.p. en.varian. projected
REAL(DP),ALLOCATABLE :: spvariancebi(:)         !  s.p. en.varian. boost-inv.
REAL(DP),ALLOCATABLE :: spenergybi(:)           !  s.p. energy boost-inv.
REAL(DP),ALLOCATABLE :: spnorm(:)               !  norm of s.p. wf
REAL(DP),ALLOCATABLE :: occup(:)                !  occupation weight
INTEGER :: nstate=1,nspdw=0                     !  Nr. of states
INTEGER :: nstate_all                           !  total Nr. of states (parayes)
INTEGER,ALLOCATABLE :: nq(:,:)                  !  nodes of init. states
INTEGER,ALLOCATABLE :: ipar(:,:)                !  xyz-parities
INTEGER,ALLOCATABLE :: ispin(:)                 !  spin of states
INTEGER,ALLOCATABLE :: nrel2abs(:),nabs2rel(:)  !  pointer to wfs
INTEGER,ALLOCATABLE :: nrel2abs_other(:,:)      !  pointer to wfs
INTEGER,ALLOCATABLE :: nhome(:)                 !  home node of wf


#if(parayes)
INTEGER,ALLOCATABLE ::  nstate_node(:)   ! number of active states in a node
INTEGER,ALLOCATABLE ::  nstart_node(:)   ! offset address for counting in a node
INTEGER,ALLOCATABLE ::  ispin_node(:,:)
#endif


!     basic parameters

LOGICAL :: directenergy=.false.            ! to compute energy directly
INTEGER :: numspin=2                       ! number of spin components
REAL(DP) :: omeg,eferm                     !  initial Gaussians
REAL(DP) :: dInMargin=6D0
INTEGER :: iforce=0
INTEGER :: ipseudo=1          ! switch on/off subgrids and pseudo-densities
REAL(DP) :: xfac,yfac,zfac                       !  initial. auxiliary
REAL(DP) :: epswf=0.2D0,e0dmp=2D0,epsorc=1D-8        !  convergence
REAL(DP) :: b2occ=0.4D0,gamocc=10D0,deocc=0D0,osfac=1D0 !  h.o.initalization
REAL(DP) :: temp=0D0,occmix=0.5D0,epsoro=1D-6        !  electron temperature
INTEGER :: isurf=0
REAL(DP) :: endcon=1D-5,radjel=4D0,surjel=1D0     !  jellium params
INTEGER :: itback=200                            !  jellium params
REAL(DP) :: betatj,gammtj,bet4tj                 !      "      "
REAL(DP) :: bbeta=0D0,gamma=0D0,beta4=0D0        !      "      "
REAL(DP) :: falph,fbeta,fhexe                    !  jellium auxiliary
REAL(DP) :: dpolx=0D0,dpoly=0D0,dpolz=0D0       !  static dipole potential
LOGICAL :: tdipolxyz                     ! switch to dipole fields
INTEGER :: iexcit=0,irotat=0,ispidi=0          !  dynam. initialisation
REAL(DP) :: phirot=0D0
INTEGER :: nclust=1,nion=1,nion2=1
!LOGICAL ::  xmnotallocated=.true. ! xm array allocated or not (to fix bug with gfortran on Debian64)
REAL(DP) :: charge             !  Nr. el. & ions
REAL(DP) :: scaleClust=1D0,scaleClustx=1D0,scaleClusty=1D0
REAL(DP) :: scaleClustz=1D0
INTEGER :: iemomsRel=1   ! 1 = relative to c.m.,  0 = relative to box
REAL(DP) :: shiftClustx=0D0,shiftClusty=0D0,shiftClustz=0D0
REAL(DP) :: rotClustx=0D0,rotClusty=0D0,rotClustz=0D0
INTEGER :: imob=1   ! fix (=0) or unfix (=1) ions
INTEGER :: ifixcmion
REAL(DP) :: shiftWFx=0D0,shiftWFy=0D0,shiftWFz=0D0
INTEGER :: ispinsep=0
REAL(DP) :: comx,comy,comz,dion(3),qtion(3,3)
REAL(DP) :: apnum,rmsion,dmdistion,rmsel,rhopss,rhomix
REAL(DP) :: codx,cody,codz,del(3),qtel(3,3)
REAL(DP) :: time_absinit
REAL(DP) :: phangle=0D0,phphase=0D0         ! angle and phase of ph rotation

INTEGER :: nhstate,npstate
INTEGER :: idebug=0,ifreezekspot=0
INTEGER :: itof=0,jescmask=0,ishiftcmtoorigin=0,jescmaskorb=0
INTEGER :: ishutdown=0
! INTEGER :: icheckformessages=1,jcheckformessages=50
INTEGER :: jplotdensitydiff=0,jplotdensity2d=0,jplotdensitydiff2d=0
INTEGER :: nmptheta=2,nmpphi=1,nmps,jmp=0,imps(maxmps)
INTEGER :: jovlp=100000000,jnorms=50
INTEGER :: iscatterelectron=0,jcharges=0,jattach=0

INTEGER,ALLOCATABLE :: ispin_target(:)
REAL(DP),ALLOCATABLE :: occ_target(:)
REAL(DP),ALLOCATABLE :: spe_target(:)
COMPLEX(DP),ALLOCATABLE :: psi_target(:,:)
INTEGER,ALLOCATABLE :: match(:,:)
REAL(DP) :: totintegprob
REAL(DP) :: aver_estar,emin_target,emax_target
INTEGER  :: nstate_target,nmatch

INTEGER :: iaddcluster=0,iswforce=0,iplotorbitals=0, ievaluate=0
REAL(DP) :: ekin0pp=0D0,vxn0=0D0,vyn0=0D0,vzn0=-1D0
REAL(DP) :: eproj=0D0,vpx=0D0,vpy=0D0,vpz=-1D0,taccel=0D0
INTEGER :: nproj=1,nproj_states=0
INTEGER,ALLOCATABLE :: proj_states(:)
REAL(DP) :: trequest=0D0,timefrac=0.98D0
REAL(DP) :: rheatclust
REAL(DP) :: igeneratesurffile
!REAL(DP),ALLOCATABLE :: rnormsinit(kstate)
REAL(DP) :: ehom0=0D0,ehomx=0D0,ehomy=0D0,ehomz=1D0
INTEGER :: ihome=0
REAL(DP) :: scatterelectronenergy=0D0,scatterelectronw=1D0
REAL(DP) :: scatterelectronvxn=0D0,scatterelectronvyn=0D0,scatterelectronvzn=1D0
REAL(DP) :: scatterelectronx=0D0,scatterelectrony=0D0,scatterelectronz=0D0
REAL(DP) :: reference_energy=0D0
REAL(DP) :: drcharges=5D0


CHARACTER (LEN=13) :: outnam
INTEGER :: iflocaliz=0                           ! evaluate localization
INTEGER :: myn                                 ! nr. of actual node
INTEGER :: ifls,ismax=1000,itmax=1000,istinf=10,ipasinf=1
INTEGER :: isitmax=0         ! number of imaginary-time steps (afterburn)
INTEGER :: idyniter=0        ! number iterations to start dynamic E0DMP 
INTEGER :: iffastpropag=1,ifexpevol=0
INTEGER :: irest=0,istat=0, isave=0,idenspl=0
INTEGER :: i3dz=0,i3dx=0,i3dstate=0,istream=0,modrho=999999
INTEGER :: jpos=-9999,jvel=-9999,jener=10,jesc=-9999,jforce=0,jposcm=0,jgeomion=0
INTEGER :: jinfo=10,jdip=-9999,jdiporb=0,jquad=0,jang=0,jspdp=0,jenergy=10
INTEGER :: jgeomel=0,jangabso=0,jelf=0,jstinf=10,jstboostinv=0
INTEGER :: jstateoverlap=0
INTEGER :: nabsorb=0,ifsicp=2,ifredmas=0,ionmdtyp=0,icooltyp=0
INTEGER :: init_lcao=0,ipsptyp=0,ivdw=0,idenfunc=1
INTEGER :: izforcecorr=-1,icooltimes=0, ntref=0
INTEGER :: jheatmod=0         ! modulus for re-heating the system
INTEGER :: ifrhoint_time=0,iangmo=0,ifspemoms=0,iftransme=0
INTEGER :: iterat,itersicp6
LOGICAL :: tstinf
REAL(DP) :: rheattemp=0D0        ! re-heat temperature
REAL(DP) :: tempion=0D0,dt1=0D0
REAL(DP) :: centfx=0D0,centfy=0D0,centfz=0D0
REAL(DP) :: shiftinix=0D0,shiftiniy=0D0,shiftiniz=0D0

REAL(DP) :: bcol1=0D0,bcol23=0D0,dbcol=0.1D0,betacol=0.99D0,chgcol=0D0
INTEGER :: ntheta=0,nphi=0
#if(parayes)
INTEGER :: ifhamdiag=0
#endif
#if(parano)
INTEGER :: ifhamdiag=20
#endif



!     spatial fields as densities and potentials
REAL(DP),ALLOCATABLE :: rhojel(:)                !  jellium density
REAL(DP),ALLOCATABLE :: potion(:)                !  pseudopotentials
REAL(DP),ALLOCATABLE :: potFixedIon(:)           !  potential from frozen ions
REAL(DP),ALLOCATABLE :: chpcoul(:)               !  Coulomb potential



!     fields and variables for analysis of electron emission

REAL(DP),ALLOCATABLE :: rhoabso(:)           !  storage for absorbed density
REAL(DP),ALLOCATABLE :: spherMask(:)         !  mask for spherical absorbing bounds
REAL(DP),ALLOCATABLE :: spherloss(:)         !  loss factor for spher. abs. bounds
REAL(DP) :: bcrad,powabso=0.125D0           ! width & power of abs. bounds
INTEGER :: ispherAbso=1               !  switch to spherical abs. bounds
INTEGER :: iangabso=0,nangtheta=1,nangphi=1
INTEGER :: ipes=0, indicesmp(maxnang*maxnang)
REAL(DP) :: angthetah=PI,angthetal=0D0
REAL(DP) :: angphil=0D0,angphih=2*PI,delomega,xango,yango,zango
LOGICAL,ALLOCATABLE :: tgridabso(:)       !  array tagging absorbing points
REAL(DP),ALLOCATABLE :: rhoabsoorb(:,:)




!     common fields for the spatial moments

!      common /moment/ ql(kmom),
INTEGER,PARAMETER :: kmom=35
INTEGER :: nrmom
REAL(DP) :: qe(kmom),qeproj(kmom),qetarget(kmom),se(5),ajx,ajy,ajz
REAL(DP),ALLOCATABLE :: qeorb_all(:,:)
!COMMON /moment/ qe,se,ajx,ajy,ajz,nrmom

! storage for the case of 1ph rotation (see 'phangle')
COMPLEX(DP),ALLOCATABLE :: oldhole(:),newhole(:)

! storage for base wavefunctions in case of dynamics with exact exchange 
! or in case of computing overlaps with initial state
COMPLEX(DP),ALLOCATABLE :: psisavex(:,:)

!     the energy transmitted from calc-lda to info etc

REAL(DP) :: enrear,ecback,ecrho,ecorr,dt12,sgaus,ekion,energy
REAL(DP) :: energ2,enerpw,encoulsp,entrop,epot,espnb,esh1
REAL(DP) :: etot,ekionold,qold2,qold3,qold4
REAL(DP) :: ekmat=0D0,engg,enii,enig,ecrhoimage=0D0
REAL(DP),ALLOCATABLE :: ekinsp(:),evarsp(:),evarsp2(:),epotsp(:)
INTEGER :: jekion,iquery4
#if(twostsic)  
REAL(DP),ALLOCATABLE :: hmatrix(:,:)
COMPLEX(DP),ALLOCATABLE :: expmatrix(:,:)
REAL(DP) :: symcon
INTEGER :: ndims(2)
#endif


! dynamic variables of ionic motion

INTEGER,PARAMETER :: nxsg=7,nysg=7,nzsg=7      ! size of subgrids
INTEGER :: modionstep=1                        ! modulus for ion step
INTEGER :: inewforce
INTEGER :: mzforce=0,myforce=0,mxforce=0       ! symmetrized forces
INTEGER :: nrare=0,nfix=0                    !  Nr. of raregas atoms
INTEGER,ALLOCATABLE :: nfixed(:)                    !  Nr. of fixed ions
INTEGER :: idielec=0                     !  switch to dielectricum
LOGICAL,ALLOCATABLE :: tblock(:)
REAL(DP),ALLOCATABLE :: cx(:),cy(:),cz(:) ,cpx(:),cpy(:),cpz(:) 
REAL(DP),ALLOCATABLE :: dcx(:),dcy(:),dcz(:) ,dcpx(:),dcpy(:),dcpz(:) 
REAL(DP),ALLOCATABLE :: fx(:),fy(:),fz(:), flx(:),fly(:),flz(:)
REAL(DP),ALLOCATABLE :: fprojx(:),fprojy(:),fprojz(:)

!                                      book keeping for LCGO initialization
REAL(DP),ALLOCATABLE :: radini(:)
INTEGER,ALLOCATABLE :: initord(:,:),ipol(:)
INTEGER,ALLOCATABLE :: nmaxst(:)



!     parameters for simulated annealing
REAL(DP),ALLOCATABLE :: eloc(:),enoloc(:),eion(:,:),eiinew(:)
REAL(DP) :: cptemp,delpos,ERR,binerg, errtot,erfac1,trfac1,prfac1
REAL(DP) :: trfac2,prfac2,errsim,eiontot,enoloctot,facann
REAL(DP) :: eloctot,errks0,errks,sumvar,sumvar2
INTEGER :: ionsin,nrun,nloop1,nloop2,loop1,nyes,ncon, ncsim,iknow
INTEGER :: isize_seed
INTEGER,ALLOCATABLE :: rand_seed(:)



!     parameters for external excitations by laser or projectile
INTEGER :: itft=3
REAL(DP) :: tnode=0D0,deltat=0D0,tpeak=0D0,omega=0D0,e0=0D0,time,tfs=0D0  
REAL(DP) :: e1x=1D0,e1y=0D0,e1z=0D0,phi=0D0
REAL(DP) :: e2x=0D0,e2y=0D0,e2z=0D0,phase2=0D0,omega2=0D0,e0_2=0D0
REAL(DP) :: tstart2=0D0,tpeak2
REAL(DP) :: fl(6),power
REAL(DP) :: elaser
INTEGER  :: ijel
REAL(DP) :: acc1old,acc2old,foft1old,foft2old,timeold
INTEGER  :: ilas=0
REAL(DP) :: fpulseinteg1        ! integrated pulse for gauge trasnf
REAL(DP) :: fpulseinteg2        ! integrated pulse for gauge trasnf.
REAL(DP) :: projcharge=0D0                   ! projectile charge
REAL(DP) :: projvelx=0D0,projvely=0D0,projvelz=0D0   ! projectile velocity
REAL(DP) :: projinix=0D0,projiniy=0D0,projiniz=0D0   ! initial projectile position
                   ! impact parameter = min(projinix,projiniy,projiniz)



! workspace for communication
REAL(DP) :: rvectmp2(3),rvectmp(3)
INTEGER :: iindtmp(3)

! pointer for energy-density functional
PROCEDURE(),POINTER :: calc_lda

#if(parayes)
INTEGER,ALLOCATABLE :: lengnod(:),displ(:)
#endif

#if(fftw_gpu)
INTEGER :: num_gpus !total number of gpus on the node
INTEGER :: mygpu !number of the actual gpu used by the node
#endif

!                          these includes should be shifted to own modules
#if(raregas)
#include "surf.F90"
#endif

#include "pseudo.F90"

CONTAINS

SUBROUTINE init_baseparams()

nthr=numthr-1  
! deduced grid parameters
 kxmax=kxbox/2+1;kymax=kybox/2+1;kzmax=kzbox/2+1
 nx2=kxbox;ny2=kybox;nz2=kzbox
 nxyz=nx2*ny2*nz2;nyf=nx2;nxyf=nx2*ny2
 kdfull2=kxbox*kybox*kzbox
 nx=nx2/2;ny=ny2/2;nz=nz2/2
! deduce nr. of states per node, readjust total nr. of states
! note: the input variable 'kstate' means total nr. of states
!       and is copied first to the correct variable 'ksttot'.
 ksttot = kstate
 kstate = ksttot/knode
 ksttot2 = knode*kstate    ! trial value
 if(ksttot2 < ksttot) kstate=1+kstate
 ksttot = knode*kstate    ! final value
 knodem=knode-1
#if(gridfft)
! bounds of loops
 minx=1;maxx=nx2
 miny=1;maxy=ny2
 minz=1;maxz=nz2
! bounds of esc. el. ???
 nbnx=2;nbxx=nx2-1
 nbny=2;nbxy=ny2-1
 nbnz=2;nbxz=nz2-1
! mid-point for x,y,z-values
 nxsh=nx2/2;nysh=ny2/2;nzsh=nz2/2
#endif
#if(findiff|numerov)
! bounds of loops
 minx=-nx;maxx=nx
 miny=-ny;maxy=ny
 minz=-nz;maxz=nz
! bounds of esc. el.
 nbnx=-nx+1;nbxx=nx-1
 nbny=-ny+1;nbxy=ny-1
 nbnz=-nz+1;nbxz=nz-1
! offset for x,y,z-values   ???
 nxsh=0;nysh=0;nzsh=0
#endif

! max. nr. of ions
ng=nion


END SUBROUTINE init_baseparams


SUBROUTINE init_fields()


ALLOCATE(nq(3,ksttot))                         !  nodes of init. states
ALLOCATE(ipar(3,ksttot))                       !  xyz-parities
ALLOCATE(ispin(ksttot))                        !  spin of states
ALLOCATE(nrel2abs(kstate),nabs2rel(ksttot))    !  pointer to wfs
ALLOCATE(nrel2abs_other(kstate,0:knode-1))     !  pointer to wfs
ALLOCATE(nhome(ksttot))                        !  home node of wf
#if(parayes)
ALLOCATE(nstate_node(0:knodem))                !  book-keeping parallele
ALLOCATE(nstart_node(0:knodem))                !          "
ALLOCATE(ispin_node(ksttot,0:knodem))
#endif


ALLOCATE(amoy(kstate))                        !  single particle energies
amoy=0D0
ALLOCATE(enonlo(kstate))                        !  single particle energies
enonlo=0D0
ALLOCATE(spvariance(kstate))                  !  s.p. energy variances
ALLOCATE(spvariancep(kstate))                 !  s.p. energy variances
ALLOCATE(spvariancebi(kstate))                !  s.p. energy variances
ALLOCATE(qeorb_all(ksttot,11))                !  s.p. dipole moments
ALLOCATE(spenergybi(kstate))
ALLOCATE(spnorm(kstate))                      !  norm of s.p. wf
ALLOCATE(occup(kstate))                       !  occupation weight

ALLOCATE(rhojel(kdfull2))                !  jellium density
ALLOCATE(potion(kdfull2))                !  pseudopotentials
ALLOCATE(potFixedIon(kdfull2))           !  potential from frozen ions
ALLOCATE(chpcoul(kdfull2))               !  Coulomb potential

IF(nabsorb > 0) THEN
   ALLOCATE(rhoabso(kdfull2))           !  storage for absorbed density
   ALLOCATE(spherMask(kdfull2))         !  mask for spherical absorbing bounds
   ALLOCATE(spherloss(kdfull2))         !  loss factor for spher. abs. bounds
   ALLOCATE(tgridabso(kdfull2))
END IF

ALLOCATE(ekinsp(kstate), evarsp(kstate),evarsp2(kstate),epotsp(kstate))
ekinsp=0D0
evarsp=0D0
evarsp2=0D0
epotsp=0D0

#if(twostsic)  
ALLOCATE(hmatrix(kstate,kstate))
#endif

!                                       fields for PsP projectors
IF(ipsptyp /=0) THEN
  ALLOCATE(ifin(0:ng),icount(knl,0:ng))
  ALLOCATE(p0_1(knl,0:ng),p0_2(knl,0:ng),p1_1(knl,0:ng),p1_1x(knl,0:ng))
  ALLOCATE(p1_1y(knl,0:ng),p1_1z(knl,0:ng) )
  ALLOCATE(p0_3(knl,0:ng),p1_2(knl,0:ng),p1_2x(knl,0:ng) )
  ALLOCATE(p1_2y(knl,0:ng),p1_2z(knl,0:ng) )
  ALLOCATE(p2_1(knl,0:ng),p2_xy(knl,0:ng),p2_xz(knl,0:ng) )
  ALLOCATE(p2_yz(knl,0:ng),p2_xy2(knl,0:ng),p2_z2(knl,0:ng))
  ALLOCATE(tnonloc(0:ng))
END IF

!                                      fields for kionic variables
ALLOCATE(tblock(0:ng))
ALLOCATE(nfixed(ng),np(0:ng)) 
ALLOCATE(cx(ng),cy(ng),cz(ng) ,cpx(ng),cpy(ng),cpz(ng))
ALLOCATE(dcx(ng),dcy(ng),dcz(ng) ,dcpx(ng),dcpy(ng),dcpz(ng))
ALLOCATE(fx(ng),fy(ng),fz(ng), flx(ng),fly(ng),flz(ng))
ALLOCATE(fprojx(ng),fprojy(ng),fprojz(ng))

!                                      book keeping for LCGO initialization
ALLOCATE(radini(1:ng),ipol(1:ng))
ALLOCATE( initord(1:3,1:ng),nmaxst(1:ng))

IF(icooltyp == 3) ALLOCATE(eloc(0:ng),enoloc(0:ng),eion(ng,ng),eiinew(ng))


#if(raregas)
IF(ivdw == 1) ALLOCATE(potvdw(kdfull2),frho(ng,3))
ALLOCATE(rfieldtmp(kdfull2))
ALLOCATE(potstat(kdfull2),potesfixed(kdfull2))
ALLOCATE(potesmob(kdfull2), phim(kdfull2),phimv(kdfull2),phimd(kdfull2))
ALLOCATE(fxt(ng),fyt(ng),fzt(ng))
maxpar = MAX(ng,ngpar) ! max of any used particle
#endif

RETURN

END SUBROUTINE init_fields

#if(raregas)

SUBROUTINE init_raregas()
iprifixed=0
ipotfixed=0
ifmdshort=1
ifadiadip=0
enerinfty=0D0
chgc0=6.119D0
chge0=-6.119D0
chgk0=2D0

iararlj=1
scaledist=1D0
distlayers=4D0
disttolerance=0.5D0
nunflayc=20
nunflaye=20
nunflayk=20
runfrowc=20D0
runfrowe=20D0
runfrowk=20D0
fermia=0D0  ! serves in addition as a switch: if isrtyp(i,j)=3 for any i,j then fermia must be given explicitly
fermib=1D0
fermic=1D0
fermiac=0D0
fermibc=0D0
fermicc=0D0
fermiae=0D0
fermibe=0D0
fermice=0D0
fermiak=0D0
fermibk=0D0
fermick=0D0
fermia2c=0D0
fermib2c=0D0
fermic2c=0D0
fermia2e=0D0
fermib2e=0D0
fermic2e=0D0
fermia2k=0D0
fermib2k=0D0
fermic2k=0D0
ccel6=0D0
ccel8=0D0
ccel10=0D0
ckel6=0D0
ckel8=0D0
ckel10=0D0



sigkv=1D0
bkv=0D0
ckv=0D0
sigvv=1D0
bvv=0D0
cvv=0D0
sigkk=1D0
bkk=0D0
ckk=0D0
sigcv=1D0
bcv=0D0
ccv=0D0
sigck=1D0
bck=0D0
cck=0D0
sigcc=1D0
bcc=0D0
ccc=0D0
sigcn=1D0
bcn=0D0
ecn=0D0
fcn=0D0
ccn=0D0
ccn2=0D0
ccn3=0D0
ccn4=0D0
ccn5=0D0
ccn6=0D0
ccn7=0D0
ccn8=0D0
ccn9=0D0
ccn10=0D0
ccnd=0D0
dcn=1D0
ccncu2=1.0D-10
ccncu3=1.0D-10
ccncu4=1.0D-10
ccncu5=1.0D-10
ccncu6=1.0D-10
ccncu7=1.0D-10
ccncu8=1.0D-10
ccncu9=1.0D-10
ccncu10=1.0D-10
ccncud=1.0D-10
ccncul2=1.0D10
ccncul3=1.0D10
ccncul4=1.0D10
ccncul5=1.0D10
ccncul6=1.0D10
ccncul7=1.0D10
ccncul8=1.0D10
ccncul9=1.0D10
ccncul10=1.0D10
ccnculd=1.0D10

sigkn=1D0
bkn=0D0
ekn=0D0
fkn=0D0
sigkn2=1D0
bkn2=0D0
ekn2=0D0
fkn2=0D0
ckn=0D0
ckn2=0D0
ckn3=0D0
ckn4=0D0
ckn5=0D0
ckn6=0D0
ckn7=0D0
ckn8=0D0
ckn9=0D0
ckn10=0D0
cknd=0D0
dkn=1.0D0
ckncu2=1.0D-10
ckncu3=1.0D-10
ckncu4=1.0D-10
ckncu5=1.0D-10
ckncu6=1.0D-10
ckncu7=1.0D-10
ckncu8=1.0D-10
ckncu9=1.0D-10
ckncu10=1.0D-10
ckncud=1.0D-10
ckncul2=1.0D10
ckncul3=1.0D10
ckncul4=1.0D10
ckncul5=1.0D10
ckncul6=1.0D10
ckncul7=1.0D10
ckncul8=1.0D10
ckncul9=1.0D10
ckncul10=1.0D10
cknculd=1.0D10


END SUBROUTINE init_raregas
#endif


END MODULE params
