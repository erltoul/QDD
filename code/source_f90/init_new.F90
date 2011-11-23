#include"define.h"
 
!------------------------------------------------------------

SUBROUTINE initnamelists

!     Sets defaults for input variables
!     and then reads input variables through namelist.

!------------------------------------------------------------
USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

CHARACTER (LEN=80) :: title
!      character*10 fname(0:3)
!      data fname/'for005.011',
!     &           'for005.100','for005.101','for005.110'/

PARAMETER (kparall=11)
CHARACTER (LEN=10) :: fname(0:kparall)
!      data fname/'for005.001','for005.010','for005.011',
!     &           'for005.100','for005.101','for005.110'/
DATA fname/'for005.001','for005.010','for005.011',  &
    'for005.100','for005.101','for005.110',  &
    'forjel.001','forjel.010','forjel.011',  &
    'forjel.100','forjel.101','forjel.110'/

NAMELIST /global/   nclust,nion,charge,nspdw,nion2,nc,nk,  &
    temp,occmix,isurf,b2occ,gamocc,deocc,osfac,  &
    init_lcao,nxin,nyin,nzin,dx,dy,dz,  &
    radjel,surjel,bbeta,gamma,beta4,endcon,itback,  &
    epswf,e0dmp,epsoro,iomax,dpolx,dpoly,dpolz,  &
    bcol1,bcol2,dbcol,ntheta,nphi,betacol,chgcol, iappendtounit,scaleclust,  &
    scaleclustx,scaleclusty,scaleclustz, shiftclustx,shiftclusty,shiftclustz,  &
    isavelastwf,iuselastwf,isavelastion,  &
    rotclustx,rotclusty,rotclustz,imob,  &
    istopcont,idebug,ishiftcmtoorigin,iaddcluster,  &
    iswforce,iplotorbitals,ievaluate,ehom0,ihome,  &
    ehomx,ehomy,ehomz,shiftwfx,shiftwfy,shiftwfz, ispinsep


!*************************************************************

!k use of parameters for dynamics:

!  use of discrete ions via pseudopotentials:      nion2 not=0
!  use of homogenous background via jellium:       nion2=0

!  iforce=0
!  polarized clusters (e.g. polar. isomer of na_12)   iforce=1

!  dipoleshift or spindipoleshift      iexcit=0:
!    shift of electronic density by centfx/y/z       ispidi=0
!    shift of spinup densitiy by +0.5*centfx/y/z     ispidi=1
!             spindown densitiy by -0.5*centfx/y/z   ispidi=1

!  scissor mode excitation             iexcit=1:
!    chose axis of rotation by 'irotat'
!    chose angle of rotation by 'phirot'

!*************************************************************

NAMELIST /dynamic/ nabsorb,  &
    iemomsrel,ifsicp,ionmdtyp,ifredmas,icooltyp,ipsptyp,  &
    ipseudo,ismax,itmax,isave,istinf,ipasinf,dt1,irest,  &
    centfx,centfy,centfz,ispidi,iforce,iexcit,  &
    irotat,phirot,i3dz,i3dx,i3dstate,istream,iflocaliz,  &
    modrho,jpos,jvel,jener,jesc,jforce,istat,jgeomion,  &
    jdip,jquad,jang,jangabso,jspdp,jinfo,jenergy,ivdw,  &
    jposcm,mxforce,myforce,mzforce,jgeomel,jelf,jstinf, &
    tempion,idenspl,ekmat,nfix,  &
    itft,tnode,deltat,tpeak,omega,e0,  &
    projcharge,projvelo,impactb,projini,pimpactb, e1x,e1y,e1z,e2x,e2y,e2z,phi,  &
    izforcecorr, iclassicmd, dinmargin,  &
    iextfield,ntref,iangabso,ipes,nangtheta,nangphi,  &
    delomega,angthetal,angthetah,angphil,angphih,  &
    ifreezekspot,powabso,ispherabso,ifixcmion,  &
    ekin0pp,vxn0,vyn0,vzn0,jescmask,itof,jescmaskorb,  &
    trequest,timefrac,icheckformessages,  &
    jcheckformessages,nmptheta,nmpphi,jmp,jovlp,  &
    jnorms,jplotdensitydiff,jplotdensitydiff2d,  &
    jplotdensity2d,jcharges,drcharges, iscatterelectron,scatterelectronenergy,  &
    scatterelectronvxn,scatterelectronvyn, scatterelectronvzn,scatterelectronx,  &
    scatterelectrony,scatterelectronz,scatterelectronw


NAMELIST /surface/  ipotstatic,iprifixed,surftemp,ipotfixed,  &
    sigmac,sigmav,sigmak,isystem,jsavesurf, chgc0,chge0,chgk0,ifreezedipoles,  &
    cspr,mion,me,mkat,isrtyp,isrtypall,  &
    shiftx,shifty,shiftz,scaledist,scaledistx,scaledisty,  &
    scaledistz,iprintonlyifmob, iaxis,distlayers,disttolerance,  &
    nunflayc,nunflaye,nunflayk, runfrowc,runfrowe,runfrowk,  &
    fermia,fermib,fermic, fermiac,fermibc,fermicc,  &
    fermiae,fermibe,fermice, fermiak,fermibk,fermick,  &
    fermia2c,fermib2c,fermic2c, fermia2e,fermib2e,fermic2e,  &
    fermia2k,fermib2k,fermic2k, sigkv,bkv,ckv,  &
    sigvv,bvv,cvv, sigkk,bkk,ckk,  &
    sigcv,bcv,ccv, sigck,bck,cck,  &
    sigcc,bcc,ccc, sigcn,bcn,ccn,  &
    sigkn,bkn,ckn, sigkn2,bkn2,  &
    ccn2,ccn3,ccn4,ccn4,ccn5,ccn6,ccn7,ccn8,ccnd,dcn,  &
    ccn9,ccn10,ccncu2,ccncu3,ccncu4,ccncu5,ccncu6,ccncu7,  &
    ccncu8,ccncu9,ccncu10,ccncud,  &
    ccncul2,ccncul3,ccncul4,ccncul5,ccncul6,ccncul7,  &
    ccncul8,ccncul9,ccncul10,ccnculd,  &
    ckn2,ckn3,ckn4,ckn4,ckn5,ckn6,ckn7,ckn8,cknd,dkn,  &
    ckn9,ckn10,ckncu2,ckncu3,ckncu4,ckncu5,ckncu6,ckncu7,  &
    ckncu8,ckncu9,ckncu10,ckncud,  &
    ckncul2,ckncul3,ckncul4,ckncul5,ckncul6,ckncul7,  &
    ckncul8,ckncul9,ckncul10,cknculd, ccel6,ccel8,ccel10,  &
    ckel6,ckel8,ckel10, ecn,fcn,ekn,fkn,ekn2,fkn2,  &
    iunfixall,ifixall,unfixclateralrad, unfixelateralrad,unfixklateralrad,  &
    unfixclateralrady,unfixelateralrady, unfixklateralrady,  &
    unfixclateralradx,unfixelateralradx, unfixklateralradx,  &
    fixcbelow,fixebelow,fixkbelow, fixcbelowy,fixebelowy,fixkbelowy,  &
    fixcbelowx,fixebelowx,fixkbelowx, unfixcrad,unfixerad,unfixkrad,  &
    iusecell,iuselast,ibh,cbh,iforcecl2co, enerinfty,epsdi,idielec,xdielec,iararlj




WRITE(6,*) 'Reading for005.// ...'


!     initialize the variables with default values


nclust = 1
nion = 1
charge = 0D0
nspdw = 0
idebug=0
temp = 0D0
occmix = 0D0
isurf = 0
b2occ = 0.4
gamocc = 0D0
deocc = 0.2
osfac = 1D0
init_lcao = 0
nxin = nx2/2
nyin = ny2/2
nzin = nz2/2
dx = 0D0
dy = 0D0
dz = 0D0
radjel = 4.0
surjel = 1D0
bbeta = 0D0
gamma = 0D0
beta4 = 0D0
endcon = 1.0D-5
itback = 200
epswf = 0.2
e0dmp = 0.9
epsoro = 1.0D-6
dpolx = 0D0
dpoly = 0D0
dpolz = 0D0
bcol1 = 0D0
bcol2 = 0D0
dbcol= 0.1
ntheta = 0
nphi = 0
betacol = 0.99
chgcol = 0
iappendtounit = 0
icomefromabso=0
isavelastwf=0
isavelastion=0
iuselastwf=0
iaddcluster=0
iswforce=0
iplotorbitals=0
ievaluate=0
iemomsrel=1      ! 1 = relat. to center of mass of cluster
! 0 = relat. to origin of box


ehom0=0D0
ihome=0
ehomx=0D0
ehomy=0D0
ehomz=1D0


nion2 = 1
nc=0
nk=0
tfs=0D0
powabso=0.125D0            ! power of mask in absorbing bounds
ispinsep=0                 ! initialize some spin asymmetry
ispherabso=0
scaleclust=1D0             ! scales inter-ion distances of
scaleclustx=1D0            ! cluster ions; x,y,z variables allow
scaleclusty=1D0            ! irregular distortion
scaleclustz=1D0
shiftclustx=0D0            ! shifts Na ions in space as compared
shiftclusty=0D0            !   to for005ion file
shiftclustz=0D0            !
ishiftcmtoorigin=0
rotclustx=0D0
rotclusty=0D0
rotclustz=0D0
imob=1            ! fix (=0) or unfix (=1) Na ions
istopcont=0
shiftwfx=0D0
shiftwfy=0D0
shiftwfz=0D0
nabsorb=0         ! size of mask for adsorbing b.c.
ifsicp=2          ! type of SIC
ionmdtyp=1        ! switch on/off ionic propagation
ifredmas=0        ! scale down ionic masses
icooltyp=0        ! type of cooling dynamics
icooltimes=0      ! not part of namelist, but initialized here, nevertheless
jheatmod=0         ! modulus for re-heating the system
rheattemp=0D0       ! re-heat temperature
ipsptyp=0          ! traditional (0) or Goedecker (1) PsP
ipseudo=1          ! swich on/off subgrids and pseudo-densities
dinmargin=6.
ismax=4000
itmax=0
isave=itmax
istinf=10
ipasinf=1
dt1=0.1
irest=0             ! restart: index of last iteration before crash
centfx=0D0        ! give boost in x-direction
centfy=0D0        ! give boost in y-direction
centfz=0D0        ! give boost in z-direction
ispidi=0
iforce=0
iexcit=0
irotat=0
phirot=0D0
i3dz=0
i3dx=0
i3dstate=0
istream=0
iflocaliz = 0
modrho=50
jpos=10
jposcm=10
jvel=10
jener=10
jesc=10
jelf=0
jstinf=0
jforce=10
jgeomion=10
jgeomel=0
jescmask=1000
jescmaskorb=1000
itof=0
istat=0
jdip=10
jquad=0
jangabso=10
jang=0
jspdp=0
jinfo=10
jenergy=10
ivdw=0
mxforce=0
myforce=0
mzforce=0
tempion=0D0
idenspl=0
ekmat=0D0
nfix=0
itft=4
tnode=1.e-5
deltat=100D0
tpeak=50D0
omega=0.34
e0=0D0
e1x=1D0
e1y=0D0
e1z=0D0
e2x=0D0                    ! default is: no second pulse
e2y=0D0
e2z=0D0
projcharge=0D0
projvelo=0D0
impactb=0D0
pimpactb=0D0
projini=0D0
phi=0D0
izforcecorr=-1
iextfield=0
ntref=1E9
iangabso=0
ipes=0
nangtheta=1
angthetal=0D0
angthetah=pi
nangphi=1
angphil=0D0
angphih=2*pi
delomega=(angthetah-angthetal)/4./nangtheta
delomega=MIN(delomega,(angphih-angphil)/4./nangphi)
ifreezekspot=0
ekin0pp=0D0
vxn0=0D0
vyn0=0D0
vzn0=-1.
trequest=1.0D6              ! requested walltime in minutes, will be converted to
! seconds later
timefrac=0.98            ! percentage of trequest, after which program shall shut down
! and save

icheckformessages=1
jcheckformessages=50
ishutdown=0
nmptheta=2
nmpphi=1
jmp=0
jovlp=1E9
jnorms=1E9
jplotdensitydiff=0
jplotdensitydiff2d=0
jplotdensity2d=0
iscatterelectron=0
scatterelectronenergy=0D0
scatterelectronvxn=0D0
scatterelectronvyn=0D0
scatterelectronvzn=1D0
scatterelectronw=1.5
scatterelectronx=0D0
scatterelectrony=0D0
scatterelectronz=nzsh*dz-4.*scatterelectronw
jcharges=0
drcharges=5.

ifreezedipoles=0
ipotstatic=0
iprintonlyifmob=0
iprifixed=0
ipotfixed=0
idielec=0
xdielec=0D0
iforcecl2co=0     ! every iforcecl2co time step the force from cluster to Arcores is given in the FILE force_clust2cores.outnam
isystem=2
jsavesurf=1E9
iunfixall=0
iusecell=0
iuselast=0
ibh=0
cbh=0D0
unfixclateralrad=-1.
unfixelateralrad=-1.
unfixklateralrad=-1.
unfixclateralrady=-1.
unfixelateralrady=-1.
unfixklateralrady=-1.
unfixclateralradx=-1.
unfixelateralradx=-1.
unfixklateralradx=-1.
unfixcrad=-1.
unfixerad=-1.
unfixkrad=-1.
fixcbelowy=-1E6
fixebelowy=-1E6
fixkbelowy=-1E6
fixcbelowx=-1E6
fixebelowx=-1E6
fixkbelowx=-1E6
fixcbelow=-1E6
fixebelow=-1E6
fixkbelow=-1E6
surftemp=0D0
enerinfty=0D0
epsdi = 1D0
sigmac=1.43/sq2
sigmav=1.43/sq2
sigmak=1D0
chgc0=6.119
chge0=-6.119
chgk0=2.0
cspr=6.119**2*e2/11.08
mion=39.95
me=4.38/1836.0
mkat=10D0
DO i=1,5
  DO j=1,5
    isrtyp(i,j) = 0
  END DO
END DO
isrtypall=0
iararlj=1
shiftx=0D0
shifty=0D0
shiftz=0D0
scaledist=1D0
iaxis=0
distlayers=4.0
disttolerance=0.5
nunflayc=20
nunflaye=20
nunflayk=20
runfrowc=20D0
runfrowe=20D0
runfrowk=20D0
fermia=0D0  ! serves in addition as a switch: if isrtyp(i,j)=3 for any i,j then fermia must be given explicitely
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
dcn=1.
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

sigkn=1.
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
dkn=1.
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



!     open input files


#if(simpara)
IF(myn >= 0 .AND. myn <= kparall) THEN
  WRITE(*,*) ' open for005: myn=',myn,fname(myn)
  iu = 7+myn
  WRITE(iu,*) ' open for005: myn=',myn,fname(myn)
  OPEN(UNIT=5,STATUS='old',FORM='formatted',FILE=fname(myn))
  WRITE(6,*) ' for005 opened: myn=',myn
  WRITE(iu,*) ' for005 opened: myn=',myn
#else
  IF(myn == 0)THEN
    OPEN(UNIT=5,STATUS='old',FORM='formatted',FILE='for005')
    iu=7
#endif
    WRITE(*,*) ' enter title (=qualifier) for that run:'
    READ(5,*)  title
    outnam = title(1:13)
    IF(outnam == '   ')  STOP " no title given "
    WRITE(*,*) ' title is now: '//outnam
    WRITE(iu,*) ' title is now: '//outnam
    CLOSE(5)
    
    OPEN(5,STATUS='old',FORM='formatted',FILE='for005.'//outnam)
    
    
    
!     read input from namelists
    
    READ(5,global,END=99999)
    WRITE(*,*) ' GLOBAL read'
    IF(dx*dy*dz == 0D0) STOP ' DX & DY & DZ must be given explicitely'
    READ(5,dynamic,END=99999)
    WRITE(*,*) ' DYNAMIC read'
    READ(5,surface,END=99999)
    WRITE(*,*) ' SURFACE read'
    
    
    99999    CLOSE(5)
#if(simpara)
  ELSE IF(myn > 5) THEN
    STOP ' this node not  active'
#endif
  END IF
  
  
#if(parayes)
  CALL comm_inputparams()
#endif
  
  
  
!     initialize and read parameters for simulated annealing:
  
  IF(icooltyp == 3) THEN
    ionsin    = 3
    iknow     = 0
    facann    = 1.5
    nrun      = 100
    nloop1    = 100
    nloop2    = 200
    cptemp    = 0.008
    delpos    = 0.2
    ERR       = 6.0D-5
    ifall     = 121466
    trfac2    = 0.03
    prfac2    = 0.05
    errsim    = 0.001
    ncsim     = 10
    errtot    = 0.002
    ncon      = 8
    erfac1    = 0.1
    trfac1    = 0.08
    prfac1    = 0.05
    errks0    = 5.0D-3
    IF(myn == 0)THEN
      
      OPEN(UNIT=44,STATUS='old',FORM='formatted', FILE='ionen-in.'//outnam)
      WRITE(6,'(/a)') 'PARAMETERS FOR METROPOLIS:'
      READ(44,*) ionsin,iknow,facann
      WRITE(6,'(a,2i2,f7.3)') ' ionsin,iknow,facann: ',ionsin, iknow,facann
      READ(44,*) nrun,nloop1,nloop2
      WRITE(6,'(a,3i4)') ' nrun,nloop1,nloop2: ',nrun,nloop1 ,nloop2
      READ(44,*) cptemp,delpos,ERR,errks0
      WRITE(6,'(a,4f12.6)') ' cptemp,delpos,err,errks0: ',cptemp  &
          ,delpos,ERR,errks0
      READ(44,*) ifall
      WRITE(6,'(a,i10)') ' ifall: ',ifall
      READ(44,*) trfac2,prfac2,errsim,ncsim
      WRITE(6,'(a,3f8.3,i3)') ' trfac2,prfac2,errsim,ncsim: '  &
          ,trfac2,prfac2,errsim,ncsim
      READ(44,*) errtot,ncon
      WRITE(6,'(a,f8.3,i3)') ' errtot,ncon: ' ,errtot,ncon
      READ(44,*) erfac1,trfac1,prfac1
      WRITE(6,'(a,3f8.3)') ' erfac1,trfac1,prfac1: ' ,erfac1,trfac1,prfac1
      WRITE(6,*) ' '
      WRITE(6,*) ' '
      CLOSE(44)
      
    END IF
    
    
#if(parayes)
    CALL comm_simann()
#endif
  END IF
  
!      write(*,*) ' INIT: nion2=',nion2
  RETURN
  
END SUBROUTINE initnamelists



!-----changePerio--------------------------------------------

SUBROUTINE changeperio

!     Reads pseudo potential parameters from namelist

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

!      namelist /perio/ ch(-92:92),amu(-92:92),
NAMELIST /perio/ ch,amu, cc1,cc2,crloc,  &
    h0_11g,h0_22g,h0_33g, h1_11g,h1_22g,h2_11g,  &
    dr1,dr2,prho1,prho2, r0g,r1g,r2g,radiong,  &
    nrow,tblock
!     &                crs,chs,
!x     &                cr0,cr1,cr2,
!     &                sgm1,sgm2,chg1,chg2,

OPEN(5,STATUS='old',FORM='formatted',FILE='for005.'//outnam)
READ(5,perio,END=99999)
99999 CONTINUE
CLOSE(5)
!test         write(iu,*) ' changeperio done. myn=',myn

END SUBROUTINE changeperio


!-----iparams--------------------------------------------------

SUBROUTINE iparams(psi)

!     check consistency of input paramaters, do some initializations

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(OUT)                     :: psi(kdfull2,kstate)

CHARACTER (LEN=80) :: title

!--------------------------------------------------------------


n = myn


!     some initializations

NE = nc
trequest=trequest*60D0
phi=phi*2*pi

!usew1 = .false.
!usew2 = .false.
!usew3 = .false.
! ??      nn=0
! ??      n=24

DO j=1,kstate
  DO i=1,kdfull2
    psi(i,j) = CMPLX(0D0,0D0,DP)
  END DO
END DO


WRITE(6,*) 'i3dz,i3dx,i3dstate,idenspl=' ,i3dz,i3dx,i3dstate,idenspl
qold2=0.01
qold3=0.01
qold4=0.01
jekion=0
iquery4=0
WRITE(6,*) 'jekion=',jekion


!     ceck consistency of input options


DO i=1,5
  DO j=1,5
    IF (isrtyp(i,j) == 2 .AND. iswforce == 1)  &
        STOP 'New force routine not implemented for Ar. Use iswForce=0!'
  END DO
END DO



IF(myn == 0)THEN
  IF(irest > itmax) STOP ' IREST > ITMAX is nonsense'
END IF

IF(nion2 == 0 .AND. ifsicp >= 3)  &
    STOP 'Jellium not compatible with KLI or Slater-SIC'
IF(nion2 == 2 .AND. ifsicp >= 3)  &
    STOP 'external PsP not compatible with KLI or Slater-SIC'

IF(nc+NE+nk.gt.ngpar) STOP ' not enough NGPAR for substrate'

#if(hamdiag&parayes)
IF(MOD(kstate,2) == 1) STOP ' KSTATE must be even'
#endif

#if(hamdiag&parayes)
STOP ' step with H diagonalization only on serial code'
#endif

#if(hamdiag&(findiff|numerov))
STOP ' step with H diagonalization not yet for fin.diff.'
#endif


#if(directenergy)
IF(ifsicp /= 3 .AND. ifsicp /= 4 .AND. ifsicp /= 0 .AND. ifsicp /= 6)  &
    STOP ' directenergy=1 only for Slater and KLI'
#endif

#if(!fullspin)
IF(ifsicp >= 3) STOP 'IFSICP>2 requires fullspin code'
#endif

#if(exchange&fullsic)
STOP ' do not compile EXCHANGE and FULLSIC simultaneously'
#endif

#if(twostsic)
#if(fullsic)
STOP ' TWOSTSIC and FULLSIC cannot run simultaneously'
#endif
#if(locsic)
STOP ' TWOSTSIC and LOCSIC cannot run simultaneously'
#endif
#if(parayes)
STOP ' TWOSTSIC cannot yet run in parallel code'
#endif
#endif


#if(expevol)
IF(ionmdtyp /= 0)  &
    STOP ' exponential evolution not with ionic motion'    !  why?
#endif

#if(fullsic)
IF(ifsicp == 6 .AND. itmax > 0)  &
    STOP ' full SIC not yet adapetd to dynamic case'
#endif

!#if(fastpropag)
!      if(nabsorb.gt.0)
!     & stop ' absorbing bound incompatible with FASTPROPG'
!      if(ionmdtyp.gt.0)
!     & stop ' ionic dynamics incompatible with FASTPROPG'
!#endif

!#if(selpara)
!      if(idielec.ne.0 .or. isurf.ne.0)
!     &  stop 'option SELPARA not compatible with substrates'
!#endif

IF(nion2 == 2) THEN
  IF(ipseudo /= 0) STOP ' IPSEUDO=0 needed for external psp (NION2=2)'
  IF(iexcit /= 0) STOP ' IEXCIR=0 needed for external psp (NION2=2)'
  IF(ipsptyp /= 0) STOP ' IPSPTYP=0 needed for external psp (NION2=2)'
  IF(imob /= 0) STOP ' fixed ions needed for external psp (IMOB=0)'
END if


RETURN
END SUBROUTINE iparams


!-----ocoption--------------------------------------------------

SUBROUTINE ocoption(iu)

!     output of compiled and of read-in options

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

!----------------------------------------------------------------

IF(iu == 8) OPEN(8,STATUS='unknown',FORM='formatted',  &
    FILE='poptions.'//outnam)

WRITE(iu,'(a)') 'the following options are used:'
WRITE(iu,*)
#if(gunnar)
WRITE(iu,'(a)') 'gl 76 exchange-correlation functional'
#endif
#if(pw92)
WRITE(iu,'(a)') 'perdew-wang 92 exchange-correlation functional'
#endif
! parallel or serial:
#if(parayes)
WRITE(iu,*) 'parallel code: number of nodes=',knode
WRITE(iu,*) 'max. number of wf per node',kstate
#endif
#if(parano)
WRITE(iu,*) 'serial code: max. number of wf=',kstate
#endif
#if(fullspin)
WRITE(iu,'(a)') ' spin code '
#else
WRITE(iu,'(a)') ' no-spin code '
#endif
#if(hamdiag)
WRITE(iu,'(a)') ' static step: Hamiltonian in subspace diagonalized'
#else
WRITE(iu,'(a)') ' static step: Hamiltonian in subspace not diagonalized'
#endif
! fft:
#if(gridfft)
WRITE(iu,'(a)') 'fourier propagation'
WRITE(iu,'(a)') 'using netlib ffts'
#endif
#if(coufou)
WRITE(iu,'(a)') 'falr coulomb solver'
#endif
#if(coudoub)
WRITE(iu,'(a)') 'exact coulomb solver'
#endif
#if(fastpropag)
WRITE(iu,'(a)') 'faster TV propagator'
#endif
#if(expevol)
WRITE(iu,'(a)') 'time step by exponential evolution'
#endif

#if(homo)
WRITE(iu,'(a)') 'hom. electric field switched on'
#endif

#if(adiadip)
WRITE(iu,'(a)') 'substrate: adiabatic dipoles'
#else
WRITE(iu,'(a)') 'substrate: dynamic dipoles'
#endif
#if(vdw)
WRITE(iu,'(a)') 'substrate: VdW activated'
#endif


!       output compiled parameters

WRITE(iu,*)'_______________________________________________'
WRITE(iu,*)
WRITE(iu,*)'compiled parameters'
WRITE(iu,*)
WRITE(iu,*) 'maximum compiled sizes x y z',kxbox,kybox,kzbox
WRITE(iu,*) 'maximum number of ions',ng
WRITE(iu,*)
WRITE(iu,*) 'units'
WRITE(iu,*)
WRITE(iu,*) 'h bar',hbar
WRITE(iu,*) 'e2',e2
WRITE(iu,*) 'm electron',ame
WRITE(iu,*)'_______________________________________________'
WRITE(iu,*)


! pseudos
IF(ipsptyp == 0) THEN
  WRITE(iu,'(a)') 'soft local pseudopotentials (errf)'
ELSE IF(ipsptyp == 1) THEN
  WRITE(iu,'(a)') 'full goedecker pseudopotentials'
ELSE IF(ipsptyp == 2) THEN
  WRITE(iu,'(a)') 'local goedecker pseudopotentials'
ELSE
  STOP ' this type IPSPTYP not yet implemented'
END IF
IF(init_lcao == 1) THEN
  WRITE(iu,'(a)') ' initialize wavefunctions with LCGO'
ELSE
  WRITE(iu,'(a)') ' initialize wavefunctions with h.o.'
END IF
! abs bounds
IF(nabsorb > 0) THEN
  IF(ispherabso == 1) THEN
    WRITE(iu,'(a,i4,1pg13.5)')  &
        ' spherical absorbing b.c.: nabsorb,powabso=', nabsorb,powabso
  ELSE
    WRITE(iu,'(a,i4,1pg13.5)')  &
        ' cartesian absorbing b.c.: nabsorb,powabso=', nabsorb,powabso
  END IF
ELSE
  WRITE(iu,'(a)') ' no absorbing b.c.'
END IF
! MD type
IF(ionmdtyp == 0) THEN
  WRITE(iu,'(a)') 'no molecular dynamics '
ELSE IF(ionmdtyp == 1) THEN
  WRITE(iu,'(a)') 'molecular dynamics: leap-frog step'
ELSE IF(ionmdtyp == 2) THEN
  WRITE(iu,'(a)') 'molecular dynamics: Runge-Kutta step'
ELSE
  STOP ' this option for MD not valid '
END IF
IF(ifixcmion == 1) WRITE(iu,*) ' fix ionic c.m. during ionic motion'
IF(ifredmas == 1) WRITE(iu,'(a)') ' reduced ionic mass = half proton mass '
IF(icooltyp == 1) THEN
  WRITE(iu,'(a)') 'cooling with pseudo dynamics '
ELSE IF(icooltyp == 2) THEN
  WRITE(iu,'(a)') 'cooling with steepest descent '
ELSE IF(icooltyp == 3) THEN
  WRITE(iu,'(a)') 'cooling with Monte Carlo '
END IF

IF(ifsicp == 0) THEN
  WRITE(iu,'(a)') 'no sic'
ELSE IF(ifsicp == 1)  THEN
  WRITE(iu,'(a)') 'sic activated: GAM'
ELSE IF(ifsicp == 2)  THEN
  WRITE(iu,'(a)') 'sic activated: ADSIC'
ELSE IF(ifsicp == 3)  THEN
  WRITE(iu,'(a)') 'sic activated: Slater'
ELSE IF(ifsicp == 4)  THEN
  WRITE(iu,'(a)') 'sic activated: KLI'
#if(exchange)
ELSE IF(ifsicp == 5)  THEN
  WRITE(iu,'(a)') 'sic activated: exact exchange'
#else
ELSE IF(ifsicp == 5)  THEN
  STOP ' code not compiled for exact exchange'
#endif
#if(fullsic)
ELSE IF(ifsicp == 6)  THEN
  WRITE(iu,'(a)') 'sic activated: full SIC'
#else
ELSE IF(ifsicp == 6)  THEN
  STOP ' code not compiled for full SIC'
#endif
#if(twostsic)
ELSE IF(ifsicp == 7)  THEN
  WRITE(iu,'(a)') 'sic activated: GSLat'
ELSE IF(ifsicp == 8)  THEN
  WRITE(iu,'(a)') 'sic activated: full SIC'
#else
ELSE IF(ifsicp == 7)  THEN
  STOP ' code not compiled for GSlat'
ELSE IF(ifsicp == 8)  THEN
  STOP ' code not compiled for double-set SIC'
#endif
ELSE
  WRITE(iu,'(a)') 'this version of SIC not available'
END IF

#if(!expevol)
IF(itmax > 0 .AND. ifsicp >= 6) STOP  &
    ' full TDSIC requires exponential evolution'
#endif


!     dynamical options

WRITE(iu,'(a,3i6)') ' kxbox,kybox,kzbox=',kxbox,kybox,kzbox
#if(raregas)
WRITE(iu,'(a,4i6)') ' nelect,nion,nrare,nstate=',nclust,nion,nrare,nstate
#else
WRITE(iu,'(a,3i6)') ' nelect,nion,nstate=',nclust,nion,nstate
#endif
WRITE(iu,'(a,4i3,f7.2)') ' ispidi,iforce,iexcit,irotat,phirot=',  &
    ispidi,iforce,iexcit,irotat,phirot
WRITE(iu,'(a,3f8.2)') ' boost: centfx,centfy,centfz=',centfx,centfy,centfz
WRITE(iu,'(a,2i3,f7.3)') ' irest,istat,dt1=',irest,istat,dt1
WRITE(iu,'(a/a,i3,5f8.3/a,7f8.3)') ' laser:',  &
    ' itft,tnode,deltat,tpeak,omega,e0=', itft,tnode,deltat,tpeak,omega,e0,  &
    ' e1x,e1y,e1z,e2x,e2y,e2z,phi=', e1x,e1y,e1z,e2x,e2y,e2z,phi

IF (isurf /= 0) THEN
  WRITE(iu,*) '*************************************************'
  WRITE(iu,*) 'SURFACE/MATRIX PRESENT:'
  WRITE(iu,'(a,i4,a,i4)') 'GSM particles: ',nc,'  kations: ',nk
  WRITE(iu,*) '*************************************************'
END IF


WRITE(iu,*) 'CODE VERSION: ',iversion

IF(iu == 8) CLOSE(8)

RETURN
END SUBROUTINE ocoption

!-----init_output--------------------------------------------------

SUBROUTINE init_output()

!     write headers of output files

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)
CHARACTER (LEN=2) :: num

!------------------------------------------------------------------

IF(myn < 9) THEN
  WRITE(num,'(i1)') myn
  maxnum = 1
ELSE
  WRITE(num,'(i2)') myn
  maxnum=2
END IF
OPEN(UNIT=7,STATUS='unknown', FILE='for006.'//num(1:maxnum)//outnam)


!     output grid properties

IF (myn == 0) THEN
  WRITE(6,*) '::::::::::::::::::::::::::::::::::::::::'
  WRITE(6,*) 'Dimensions of the box:'
  WRITE(6,'(a,3i5,i8)')  &
      ' KXBOX,KYBOX,KZBOX,KDFULL2=',kxbox,kybox,kzbox,kdfull2
  WRITE(6,'(a,f10.2,a,f10.2)') 'x = ',(1-nxsh)*dx, ' .. ',nxsh*dx
  WRITE(6,'(a,f10.2,a,f10.2)') 'y = ',(1-nysh)*dy, ' .. ',nysh*dy
  WRITE(6,'(a,f10.2,a,f10.2)') 'z = ',(1-nzsh)*dz, ' .. ',nzsh*dz
  WRITE(6,*) '::::::::::::::::::::::::::::::::::::::::'
  WRITE(6,*) 'Dimensions of the inner margin:'
  WRITE(6,'(a,f10.2,a,f10.2)') 'x = ',(1-nxsh)*dx+dinmargin, ' .. ',  &
      nxsh*dx-dinmargin
  WRITE(6,'(a,f10.2,a,f10.2)') 'y = ',(1-nysh)*dy+dinmargin, ' .. ',  &
      nysh*dy-dinmargin
  WRITE(6,'(a,f10.2,a,f10.2)') 'z = ',(1-nzsh)*dz+dinmargin, ' .. ',  &
      nzsh*dz-dinmargin
  WRITE(6,*) '::::::::::::::::::::::::::::::::::::::::'
END IF


!     deduced grid parameters

WRITE(6,'(a,6i5,i8,g12.4)') ' MAIN: nx,ny,nz,nx1,ny1,nzi,nxyz,dvol=',  &
    nx,ny,nz,nx1,ny1,nzi,nxyz,dvol




!       output static parameters (traditional sodium case)

WRITE(7,'(5x,a,i2,a,i2/a/)')  &
    '# electr.: ',nclust,'# ions: ',nion,'=========='
WRITE(7,'(a,i3)') 'number of wave-fkts nstate = ',nstate
WRITE(7,'(a,f4.2)') 'initialisation : osfac=',osfac
WRITE(7,'(a,10(/t2,3(3i4,5x),3i4))')  &
    'quantumnumbers :  ',((nq(i,j),i=1,3),j=1,nstate)
IF(temp > 0D0) THEN
  WRITE(7,'(a,f10.5,a,f4.2)') 'temperature kt= ',temp,'Ry,   occmix=',occmix
  WRITE(7,'(a,1g13.5,a,i4)') 'epsoro=',epsoro,'  iomax=',iomax
END IF
WRITE(7,'(a,i3)') 'number of wave-fkts nclust = ',nclust
WRITE(7,'(3(a,i4,3x))') 'nx=',nx,'ny=',ny,'nz=',nz
WRITE(7,'(3(a,f4.2,3x))') 'dx=',dx,'dy=',dy,'dz=',dz
WRITE(7,'(a)') 'damped gradient step :'
WRITE(7,'(2(a,1pg12.4,3x))') 'epswf=',epswf,'e0dmp=',e0dmp
WRITE(7,'(a)') 'jellium background :'
WRITE(7,'(2(a,1pg12.4,3x))') 'radjel=',radjel,'surjel=',surjel

WRITE(7,'(a,a)')  'set of fixed deformation parameters ',  &
    '(in hill-wheeler-coordinates) :'
WRITE(7,'(a,f4.2,3x,a,f5.0)') 'bbeta=',bbeta,'gamma=',gamma
WRITE(7,'(//)')


RETURN


!  301 format('/t2 * * * * * parameter   ksmax   too small  * * * * *')
!  302 format('/t2 * * * * * parameter   kdfull  too small  * * * * *')
!  304 format('/t2 * * * * * parameter   ktax    too small  * * * * *')
!  305 format('/t2 * * * * * parameter   ktay    too small  * * * * *')
!  306 format('/t2 * * * * * parameter   ktaz    too small  * * * * *')
!  307 format('/t2 * * * * * parameter   kstsq   too small  * * * * *')
!  308 format('/t2 * * * * * parameter   kstate  too small  * * * * *')
!  309 format('/t2 * * * * * parameter   krun    too small  * * * * *')
!  310 format('/t2 * * * * * parameter   kxmax   too small  * * * * *')
!  311 format('/t2 * * * * * parameter   kymax   too small  * * * * *')
!  312 format('/t2 * * * * * parameter   kzmax   too small  * * * * *')
!  320 format('/t2 --> kstate greater initialisation of quantum numbers')


WRITE(7,*)'_______________________________________________'
WRITE(7,*)
WRITE(7,*)'        dynamic parameters have been read :'
WRITE(7,*)'_______________________________________________'
WRITE(7,*)

IF(iforce == 0) THEN
  WRITE(7,*)
  WRITE(7,*) 'unpolarized clusters (e.g. groundstate of na_12)   '
END IF
IF(iforce == 1) THEN
  WRITE(7,*)
  WRITE(7,*) 'polarized clusters (e.g. polar. isomer of na_12)   '
END IF
IF(nion2 == 0) THEN
!        write(6,*)
!       write(6,*) 'jellium code'
  WRITE(7,*)
  WRITE(7,*) 'nions2=0 : jellium code'
  ecorr=0D0
ELSE IF(nion2 == 0) THEN
  WRITE(7,*) 'external local PsP from file'
ELSE
  WRITE(7,*) 'ionic code - number of ions',nion,nion2
END IF
IF(irest == 0) THEN
  WRITE(7,*) 'number of static iterations',ismax
ELSE
  WRITE(7,*) 'restart option activated - restart at',irest
END IF
WRITE(7,*) 'saving wfs every ',isave,'iterations'
WRITE(7,*) 'number of dynamic iterations',itmax
WRITE(7,*) 'printing out observables every',ipasinf,'iterations'
WRITE(7,*) 'timestep=',dt1,'units'


IF(iexcit == 0) THEN
  WRITE(7,*) 'dipoleshift or spindipoleshift'
  IF(ispidi == 0) THEN
    WRITE(7,*) 'shift of electronic density by', centfx,centfy,centfz
  ELSE
    IF(ispidi == 1) THEN
      WRITE(7,*) 'shift of spinup density by'
      WRITE(7,*) 0.5*centfx,0.5*centfy,0.5*centfz
    ELSE IF(ispidi == -1) THEN
      WRITE(7,*) 'shift of spindown density by'
      WRITE(7,*) -0.5*centfx,-0.5*centfy,-0.5*centfz
    END IF
  END IF
ELSE
  IF(iexcit == 1) THEN
    WRITE(7,*) 'scissor mode excitation'
    WRITE(7,*) 'axis of rotation',irotat
    WRITE(7,*) 'angle of rotation',phirot
  END IF
END IF
IF(tempion /= 0D0) THEN
  WRITE(7,*) 'initial thermalization of the ions at t=',tempion
ELSE
  WRITE(7,*) 'zero initial ionic kinetic energy'
END IF

! laser

IF(e0 /= 0D0) THEN
  WRITE(7,*) 'laser pulse active -  parameters : '
  IF(itft == 1) THEN
    WRITE(7,*) 'ramp laser pulse, sine switching'
    WRITE(7,*) ' on/off up to/from f(t) = 1'
  END IF
  IF(itft == 2) THEN
    WRITE(7,*) 'gaussian laser pulse'
  END IF
  WRITE(7,*) 'length of the pulse',deltat
  WRITE(7,*) 'peak time',tpeak
  WRITE(7,*) 'begins at',tnode
  WRITE(7,*) 'pulsation',omega,'units'
  WRITE(7,*) 'field strength',e0
  IF(phi == 0D0) WRITE(7,*) 'linear polarization'
  IF(phi == 0.25) WRITE(7,*) 'circular right polarization'
  IF(phi == -0.25) WRITE(7,*) 'circular left polarization'
  ascal=e1x*e2x+e1y*e2y+e1z*e2z
  IF(ascal /= 0D000) WRITE(7,*) 'warning : e1 e2 non orthogon.'
  WRITE(7,*)  e1x,e1y,e1z,e2x,e2y,e2z,phi
ELSE
  WRITE(7,*) 'laser pulse switched off'
END IF
WRITE(7,*)'_______________________________________________'
WRITE(7,*)


RETURN
END SUBROUTINE init_output


!-----iperio-----------------------------------------------------

SUBROUTINE iperio

!     initializes the "periodic table" :
!      masses of implemented elements, charges
!      and parameters for the pseudo-potentials

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)
!      dimension dr1(-ng:ng),dr2(-ng:ng) ! bugBF
!      dimension prho1(-ng:ng),prho2(-ng:ng) ! bugBF
REAL(DP) :: dr1(-92:92),dr2(-92:92)
REAL(DP) :: prho1(-92:92),prho2(-92:92)

!----------------------------------------------------------------


WRITE(6,*) 'Entering iperio'


DO iel=-92,92
  amu(iel)=0D0   ! masses as a function of atomic number
END DO

!     soft PsP (error functions)

IF(ipsptyp == 0) THEN
  
!       hydrogen
  
  amu(1)=1D0
  ch(1)=1D0
  dr1(1)=0.3
  dr2(1)=0.45
  prho1(1)=6.94527
  prho2(1)=-0.92056
  sgm1(1)=dr1(1)*0.8493218
  sgm2(1)=dr2(1)*0.8493218
  chg1(1)=((sgm1(1)*SQRT(pi*2.0))**3)*prho1(1)*ch(1)
  chg2(1)=((sgm2(1)*SQRT(pi*2.0))**3)*prho2(1)*ch(1)
  WRITE(6,*) 'ch(1)=',ch(1)
  WRITE(6,*) 'total charge of Psphydr=',chg1(1)+chg2(1)
  
!       sodium
  
  amu(11)=23.0
  ch(11)=1D0  ! charge of pseudopotential
  dr1(11)=0.8018
  dr2(11)=1.3693
  prho1(11)=-0.46073
  prho2(11)=0.13287
  sgm1(11)=dr1(11)*0.8493218
  sgm2(11)=dr2(11)*0.8493218
  chg1(11)=((sgm1(11)*SQRT(pi*2.0))**3)*prho1(11)*ch(11)
  chg2(11)=((sgm2(11)*SQRT(pi*2.0))**3)*prho2(11)*ch(11)
  WRITE(6,*) 'ch(11)=',ch(11)
  WRITE(6,*) 'total charge of Pspsodi=',chg1(11)+chg2(11)
  
!       magnesium
  
  amu(12)=24.312
  ch(12)=2.0
  dr1(12)=0.5
  dr2(12)=1D0
  prho1(12)=-1.6582
  prho2(12)=0.31091
  sgm1(12)=dr1(12)*0.8493218
  sgm2(12)=dr2(12)*0.8493218
  chg1(12)=((sgm1(12)*SQRT(pi*2.0))**3)*prho1(12)*ch(12)
  chg2(12)=((sgm2(12)*SQRT(pi*2.0))**3)*prho2(12)*ch(12)
  WRITE(6,*) 'ch(12)=',ch(12)
  WRITE(6,*) 'total charge of Pspmagn=',chg1(12)+chg2(12)
  
!       Argon , positive core part
  
  amu(18)   = 39.95
  ch(18)    = 6.119  ! charge of pseudopotential
  sgm1(18)  = 1.43/sq2
  dr1(18)   = sgm1(18)/0.8493218
  sgm2(18)  = 1
  dr2(18)   = 1
  chg1(18)  = ch(18)
  prho1(18) = chg1(18)/((sgm1(18)*SQRT(pi*2.0))**3)
  chg2(18)  = 0D0
  prho2(18) = 0D0
  WRITE(6,*) 'ch(18)=',ch(18)
  WRITE(6,*) 'total charge of PspAr_core=',chg1(18)+chg2(18)
  
!       Argon , negative valence cloud (inert)
  
  amu(-18)  = 4.38/1836.0
  ch(-18)   = -6.119  ! charge of pseudopotential
  sgm1(-18) = 1.43/sq2
  dr1(-18)  = sgm1(-18)/0.8493218
  sgm2(-18) = 1
  dr2(-18)  = 1
  chg1(-18) = ch(-18)
  prho1(-18)= chg1(-18)/((sgm1(-18)*SQRT(pi*2.0))**3)
  chg2(-18) = 0D0
  prho2(-18)= 0D0
!  c_dipmod  = 3.378  !  ch(-18)**2*e2/11.08
  WRITE(6,*) 'ch(-18)=',ch(-18)
  WRITE(6,*) 'total charge of PspAr_el=',chg1(-18)+chg2(-18)
  
!       potassium
  
  amu(19)=39.0
  ch(19)=1D0  ! charge of pseudopotential
  dr1(19)=0.9973
  dr2(19)=1.9957
  prho1(19)=-0.13771
  prho2(19)=0.030223
  sgm1(19)=dr1(19)*0.8493218
  sgm2(19)=dr2(19)*0.8493218
  chg1(19)=((sgm1(19)*SQRT(pi*2.0))**3)*prho1(19)*ch(19)
  chg2(19)=((sgm2(19)*SQRT(pi*2.0))**3)*prho2(19)*ch(19)
  WRITE(6,*) 'ch(19)=',ch(19)
  WRITE(6,*) 'total charge of Psppotassium=',chg1(19)+chg2(19)
  
!       cesium
  
  amu(58)=140D0
  ch(58)=1D0  ! charge of pseudopotential
  dr1(58)=1.12
  dr2(58)= 2.552
  prho1(58)=-0.62923E-01
  prho2(58)=0.11554E-01
  sgm1(58)=dr1(58)*0.8493218
  sgm2(58)=dr2(58)*0.8493218
  chg1(58)=((sgm1(58)*SQRT(pi*2.0))**3)*prho1(58)*ch(58)
  chg2(58)=((sgm2(58)*SQRT(pi*2.0))**3)*prho2(58)*ch(58)
  WRITE(6,*) 'ch(58)=',ch(58)
  WRITE(6,*) 'total charge of Psppotassium=',chg1(58)+chg2(58)
  
!     standard Goedecker PsP (local part initialized here)
  
ELSE IF(ipsptyp == 1) THEN
  
!       hydrogen
  
  amu(1)   = 1D0
  ch(1)    = 1D0
  cc1(1)   =-4.180237D0
  cc2(1)   = 0.725075D0
  crloc(1) = 0.2D0
  nrow(1)  = 1
  
!       helium
  
  amu(2)   = 4.0D0
  ch(2)    = 2.0D0
  cc1(2)   =-9.112023D0
  cc2(2)   = 1.698368D0
  crloc(2) = 0.2D0
  nrow(2)  = 1
  
!       bor
  
  amu(5)   = 10.81D0
  ch(5)    = 3.0D0
  cc1(5)   =-5.578642D0
  cc2(5)   = 0.804251D0
  crloc(5) = 0.43393D0
  nrow(5)  = 2
  r0g(5)=0.373843D0
  r1g(5)=0.360393D0
  radiong(5)=2.5D0
  h0_11g(5)=6.233928D0
  h1_11g(5)=0D0
  
  
!       carbon
  
  amu(6)   = 12.0D0
  ch(6)    = 4.0D0
  cc1(6)   =-8.513771D0
  cc2(6)   = 1.228432D0
  crloc(6) = 0.348830D0
  nrow(6)  = 2
  r0g(6)=0.304533D0
  r1g(6)=0.232677D0
  radiong(6)=2.0
  h0_11g(6)=9.522842D0
  h1_11g(6)=0D0
  
!       nitrogen
  
  amu(7)   = 14.0D0
  ch(7)    = 5.0D0
  cc1(7)   =-12.2046419D0
  cc2(7)   = 1.7558249D0
  crloc(7) = 0.2889046D0
  nrow(7)  = 2
  r0g(7)=0.256605D0
  r1g(7)=0.270134D0
  radiong(7)=2.0D0
  h0_11g(7)=13.552433D0
  h1_11g(7)=0D0
  
!       oxygen
  
  amu(8)   = 16.0D0
  ch(8)    = 6.0D0
  cc1(8)   =-16.4822284D0
  cc2(8)   = 2.3701353D0
  crloc(8) = 0.2477535D0
  nrow(8)  = 2
  r0g(8)=0.221786D0
  r1g(8)=0.256829D0
  radiong(8)=1.5D0
  h0_11g(8)=18.266917D0
  h1_11g(8)=0D0
  
!       neon   (row 2)
  
  amu(10)  = 20.2D0
  ch(10)   = 8.0D0
  cc1(10)  = -27.692852D0
  cc2(10)  = 4.005906D0
  crloc(10)= 0.19D0
  nrow(10) = 2
  r0g(10) = 0.174268D0
  r1g(10) = 0.214913D0
  h0_11g(10)= 28.506098D0
  h0_22g(10)= -1.076245D0
  h1_11g(10)= -0.000090D0
  
!       sodium
  
  amu(11)  = 23.0D0
  ch(11)   = 1D0
  cc1(11)  = -1.238867D0    !nbf2.25  ! -1.238867
  cc2(11)  = 0D0
  crloc(11)= 0.885509D0     !nbf0.8   ! 0.885509
  nrow(11) = 3
  r0g(11) = 0.661104D0
  r1g(11) = 0.857119D0
  h0_11g(11)= 1.847271D0
  h0_22g(11)= 0.582004D0
  h1_11g(11)= 0.471133D0
  
!       magnesium
  
  amu(12)  = 24.312D0
  ch(12)   = 2.0D0
  cc1(12)  =-2.864297D0
  cc2(12)  = 0D0
  crloc(12)=0.651812D0
  nrow(12) = 3
  r0g(12)=0.556478D0
  r1g(12)=0.677569D0
  h0_11g(12)=2.970957D0
  h0_22g(12)=1.329941D0
  h1_11g(12)=1.049881D0
  
!       aluminium
  
  amu(13)  = 26.98D0
  ch(13)   = 3.0D0
  cc1(13)  =-8.491315D0
  cc2(13)  = 0D0
  crloc(13)= 0.45D0
  nrow(13) = 3
  r0g(13)=0.460104D0
  r1g(13)=0.536744D0
  h0_11g(13)=5.088340D0
  h0_22g(13)=2.679700D0
  h1_11g(13)=2.193438D0
  
!       silicium
  
  amu(14)  = 28.09D0
  ch(14)   = 4.0D0
  cc1(14)  =-7.336103D0
  cc2(14)  = 0D0
  crloc(14)= 0.44D0
  nrow(14) = 3
  
!       phosphor
  
  amu(15)  = 30.974D0
  ch(15)   = 5.0D0
  cc1(15)  =-6.65422D0
  cc2(15)  = 0D0
  crloc(15)= 0.43D0
  nrow(15) = 3
  r0g(15)=0.389803D0
  r1g(15)=0.440796D0
  h0_11g(15)=6.842136D0
  h0_22g(15)=3.856693D0
  h1_11g(15)=3.282606D0
  
!       argon   (row 3)
  
  amu(18)  = 40
  ch(18)   = 8.0D0
  cc1(18)  = -7.1D0 !nbf9.883 !  9.236   !  -7.1
  cc2(18)  = 0D0     !  -2.5265 !  0.0
  crloc(18)= 0.4D0     !  0.6     !  0.4
  nrow(18) = 3
  r0g(18) = 0.317381D0 ! 0.4
  r1g(18) = 0.351619D0 ! 0.4
  h0_11g(18)= 10.249487D0 ! 2.4247           ! 0.2808  ! 10.249487
  h0_22g(18)= 5.602516D0  ! 0.0              ! 0.0     ! 5.602516
  h1_11g(18)= 4.978801D0  ! 0.0              ! 0.0     ! 4.978801
  
  
!       calcium
  
  amu(20)  = 40.08D0
  ch(20)   = 2.0D0
  cc1(20)  = 0D0
  cc2(20)  = 0D0
  crloc(20)= 0.8D0
  nrow(20) = 4
  r0g(20)=0.669737D0
  r1g(20)=0.946474D0
  r2g(20)=0.526550D0
  h0_11g(20)=1.645014D0
  h0_22g(20)=1.523491D0
  h0_33g(20)=0.295996D0
  h1_11g(20)=0.585479D0
  h1_22g(20)=0.126329D0
  h2_11g(20)=-3.0323321D0
  
!       Cu
  
  amu(29)  = 63.546D0
  ch(29)   = 1D0
  cc1(29)  = 0D0
  cc2(29)  = 0D0
  crloc(29)= 0.58D0
  nrow(29) = 4
  r0g(29)=0.843283D0
  r1g(29)=1.089543D0
  r2g(29)=1.291602D0
  h0_11g(29)=0.975787D0
  h0_22g(29)=-0.822070D0
  h0_33g(29)=-0.133237D0
  h1_11g(29)=0.024580D0
  h1_22g(29)=-0.249001D0
  h2_11g(29)=-0.065292D0
  
!       Ag
  
  amu(47)  = 107.8682D0
  ch(47)   = 1D0
  cc1(47)  = -2.376061D0
  cc2(47)  = 0D0
  crloc(47)= 0.65D0
  nrow(47) = 4
  r0g(47)=1.012705D0
  r1g(47)=1.235842D0
  r2g(47)=1.016159D0
  h0_11g(47)=0.897931D0
  h0_22g(47)=-0.7483230D0
  h0_33g(47)=0.029787D0
  h1_11g(47)=0.130081D0
  h1_22g(47)=-0.277495D0
  h2_11g(47)=-0.038842D0
  
!     purely local Goedecker PsP
  
ELSE IF(ipsptyp == 2) THEN
  
!       hydrogen Gianocci PsP (approx. local goedecker)
  
  amu(1)   = 1D0
  ch(1)    = 1D0
  nrow(1)  = 1
  
!       Na (approximate local pseudo built on Geodecker local part)
  
  amu(11)  = 23.0
  ch(11)   = 1D0
  cc1(11)  = 1.5 ! 1.55 ! 2.5 ! 3.0   !  3.1
  cc2(11)  = 0D0
  crloc(11)= 0.9
  nrow(11) = 3
!        write(6,*)'params pseudo, c1,c2',cc1(18),cc2(18),crloc(18)
  
!       Ar (approximate local pseudo built on Geodecker local part)
  
  amu(18)  = 39.95
  ch(18)   = 8.0
!        cc1(18)  = 60.37    ! 13.375   !  2.5034  ! 2.482
!        cc2(18)  = -8.972   ! -3.4762  ! -1.4664  ! -1.4526
!        crloc(18)= 0.4      ! 0.6      ! 0.9      ! 0.9
  
  cc1(18)  =   2.482   ! 2.5034  !
  cc2(18)  =  -1.4526  ! -1.4664  !
  crloc(18)=  0.9
  
  
!        cc1(18)  = 4.9558   !  9.9117
!        cc2(18)  = -2.0318  ! -4.0637
!        crloc(18)= 0.8
!09        cc1(18)  = 3.85
!09        cc2(18)  = -2.44
!09        crloc(18)= 0.9
  nrow(18) = 3
!        write(6,*)'params pseudo, c1,c2',cc1(18),cc2(18),crloc(18)
  
ELSE
  STOP ' IPERIO: this type PsP not yet implemented '
END IF

!     reset ionic masses if asked for


WRITE(6,*) 'resettig ionic masses...'
IF(ifredmas == 1) THEN
  WRITE(6,*) 'in progress...'
  DO i=-92,92
!           amu(i)=0.0484
    amu(i)=0.5
  END DO
!         do i=-92,92
!            amu(i) = amu(i)*0.5/amu(11)
!mb            amu(i)=amu(i)*0.5/amu(18)
!         enddo
!mb         write(6,*)'Argon core mass is reduced to half of hydrogen'
  WRITE(6,*)'Na mass is reduced to half of hydrogen'
  WRITE(6,*) 'reduced mass is a half of hydrogen'
  scpx=0D0
  scpy=0D0
  scpz=0D0
END IF

WRITE(6,*) 'Mass of na: ',amu(11)

WRITE(6,*) 'END OF IPERIO.'

RETURN
END SUBROUTINE iperio


!-----initions-------------------------------------------------

SUBROUTINE initions()

!     read ionic positions and initialize ion-related parameters

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)
CHARACTER (LEN=3) :: orderxyz
REAL(DP) :: vecin(3),vecout(3),vecalpha(3)

IF(nion2 == 0) THEN
  ecorr = 0D0
  RETURN       ! this is the case of jellium
END IF

n = myn

WRITE (6,*) 'Entering initions()'



!     readings and basic transformation are only done for node "0"

#if(parayes)
IF(myn == 0)THEN
#endif
  
  
  OPEN(UNIT=9,STATUS='old',FORM='formatted', FILE='for005ion.'//outnam)
  
  
  WRITE(6,*) 'Reading positions...'
  DO ion=1,nion
    
    IF(init_lcao == 1) THEN
      READ(9,*) cx(ion),cy(ion),cz(ion),np(ion),orderxyz, radini(ion)
!                          translate ordering of atomic states
      IF(orderxyz == 'xyz') THEN
        initord(1,ion) = 1
        initord(2,ion) = 2
        initord(3,ion) = 3
      ELSE IF(orderxyz == 'xzy') THEN
        initord(1,ion) = 1
        initord(2,ion) = 3
        initord(3,ion) = 2
      ELSE IF(orderxyz == 'yxz') THEN
        initord(1,ion) = 2
        initord(2,ion) = 1
        initord(3,ion) = 3
      ELSE IF(orderxyz == 'yzx') THEN
        initord(1,ion) = 2
        initord(2,ion) = 3
        initord(3,ion) = 1
      ELSE IF(orderxyz == 'zxy') THEN
        initord(1,ion) = 3
        initord(2,ion) = 1
        initord(3,ion) = 2
      ELSE IF(orderxyz == 'zyx') THEN
        initord(1,ion) = 3
        initord(2,ion) = 2
        initord(3,ion) = 1
      ELSE
        STOP ' this ordering of XYZ not provided '
      END IF
    ELSE
      READ(9,*) cx(ion),cy(ion),cz(ion),np(ion)
    END IF
    
    
!       initial kinetic momenta=0
    
    cpx(ion)=0D0
    cpy(ion)=0D0
    cpz(ion)=0D0
    
  END DO
  CLOSE(UNIT=9)
  
  
  
!       re-scale cluster if desired
  
  
  IF (ABS(scaleclust-1D0) > small .OR. ABS(scaleclustx-1D0) > small .OR.  &
        ABS(scaleclusty-1D0) > small .OR. ABS(scaleclustz-1D0) > small) THEN
    
    WRITE(6,*) 'Rescaling...'
    
    
    IF (scaleclust /= 1D0) THEN
      scaleclustx=scaleclust
      scaleclusty=scaleclust
      scaleclustz=scaleclust
    END IF
    
    
    xcm = 0D0
    ycm = 0D0
    zcm = 0D0
    DO i=1,nc
      xcm = xcm + cx(i)
      ycm = ycm + cy(i)
      zcm = zcm + cz(i)
    END DO
    xcm = xcm/nion
    ycm = ycm/nion
    zcm = zcm/nion
    
    DO i=1,nion
      cx(i) = (cx(i)-xcm)*scaleclustx + xcm
      cy(i) = (cy(i)-ycm)*scaleclusty + ycm
      cz(i) = (cz(i)-zcm)*scaleclustz + zcm
    END DO
    
  END IF
  
!       re-scaling done.
  
  
  
  
!       shift cluster in space if desired
  
  
  IF (shiftclustx /= 0D0 .OR. shiftclusty /= 0D0 .OR.  &
        shiftclustz /= 0D0) THEN
    
    WRITE(6,*) 'Shifting cluster ions...'
    
    DO i=1,nion
      cx(i)=cx(i)+shiftclustx
      cy(i)=cy(i)+shiftclusty
      cz(i)=cz(i)+shiftclustz
    END DO
    
  END IF
  
  
  IF (ishiftcmtoorigin == 1) THEN
    CALL getcm(1,0,0,0)
    rvectmp2(1)=rvectmp(1)
    rvectmp2(2)=rvectmp(2)
    rvectmp2(3)=rvectmp(3)
    DO ii=1,nion
      cx(ii)=cx(ii)-rvectmp(1)
      cy(ii)=cy(ii)-rvectmp(2)
      cz(ii)=cz(ii)-rvectmp(3)
    END DO
  END IF
  
!       shift done
  
  
  
!   Rotate cluster along rotation vector 
!    (/rotclustx,rotclusty,rotclustz/).
!   Rotation vector is given initially in degree and converted
!   internally to radian.

  IF(abs(rotclustx)+abs(rotclusty)+abs(rotclustz)>0D0) THEN
    
    vecalpha(1)=rotclustx/180D0*pi
    vecalpha(2)=rotclusty/180D0*pi
    vecalpha(3)=rotclustz/180D0*pi
    
! temporary shift to CM
    ctmpx=SUM(cx(1:nion))/nion
    ctmpy=SUM(cy(1:nion))/nion
    ctmpz=SUM(cz(1:nion))/nion
    cx=cx-ctmpx
    cy=cy-ctmpy
    cz=cz-ctmpz
    
! apply rotation   
    DO i=1,nion

      vecin(1) = cx(ion)
      vecin(2) = cy(ion)
      vecin(3) = cz(ion)

      CALL rotatevec3D(vecin,vecout,vecalpha)

      cx(ion) = vecout(1) 
      cy(ion) = vecout(2) 
      cz(ion) = vecout(3) 
      
    END DO
    
!                   shift center of mass back to original position
    cx=cx+ctmpx
    cy=cy+ctmpy
    cz=cz+ctmpz
    
  END IF
!     rotation completed
  
  
  DO ion=1,nion
    IF(ipsptyp == 0) THEN
      tblock(ion) = .NOT.(ABS(np(ion)) > 99)
      np(ion) = MOD(np(ion),100)
      WRITE(6,'(a,i4,a,3g12.4,i4,1x,l1)')  &
          ' ion nr.',ion,':  x,y,z,type,block=',  &
          cx(ion),cy(ion),cz(ion),np(ion),tblock(ion)
      
      IF(np(ion) /= 11 .AND. np(ion) /= 12 .AND. np(ion) /= 1  &
            .AND. np(ion) /= 18 .AND. np(ion) /= -18  &
            .AND. np(ion) /= 19 .AND. np(ion) /= 58) THEN
        WRITE(6,'(a,1i6,a,1i6)') 'np(',ion,')=',np(ion)
        STOP 'element not provided with erf pseudo potential'
      END IF
    END IF
    IF(np(ion) > 92) STOP 'element out of range'
    IF(amu(np(ion)) == 0D0) STOP 'unknown elem. found'
    
    WRITE(7,*)'np(',ion,')=',np(ion)
    WRITE(6,*)'np(',ion,')=',np(ion)
    WRITE(7,*)'amu(',np(ion),')=',amu(np(ion))
    WRITE(6,*)'amu(',np(ion),')=',amu(np(ion))
    DO iunit=6,7
      WRITE(iunit,*) ' PsP parameters for ion:'
      WRITE(iunit,'(a,1pg14.6)') 'ch=',ch(np(ion)),  &
          'amu=',amu(np(ion)), 'cc1=',cc1(np(ion)),  &
          'cc2=',cc2(np(ion)), 'crloc=',crloc(np(ion)),  &
          'crs=',crs(np(ion)), 'chs=',chs(np(ion)),  &
          'h0_11g=',h0_11g(np(ion)), 'h0_22g=',h0_22g(np(ion)),  &
          'h0_33g=',h0_33g(np(ion)), 'h1_11g=',h1_11g(np(ion)),  &
          'h1_22g=',h1_22g(np(ion)), 'h2_11g=',h2_11g(np(ion)),  &
          'sgm1=',sgm1(np(ion)), 'sgm2=',sgm2(np(ion)),  &
          'chg1=',chg1(np(ion)), 'chg2=',chg2(np(ion)),  &
          'nrow=',nrow(np(ion)), 'r0g=',r0g(np(ion)),  &
          'r1g=',r1g(np(ion)), 'r2g=',r2g(np(ion)),  &
          'radiong=',radiong(np(ion)), 'tblock=',tblock(np(ion))
    END DO
    
    IF(ipsptyp == 1) THEN
      dgrid = crloc(np(ion))/0.8493218
      WRITE(7,*) ' optimal grid spacing=',dgrid
      WRITE(6,*) ' optimal grid spacing=',dgrid
    END IF
  END DO
  
  
  IF (tempion > 0D0 .AND.imob /= 0) THEN
    CALL givetemperature(cpx,cpy,cpz,nion,tempion, amu(np(1))*1836*ame,4)
  END IF
  
#if(parayes)
END IF                                             ! myn=0
#endif

!     Part for node "0" finished. Now distribute to all nodes.

#if(parayes)
CALL comm_ionconfig()
#endif


!        check consistency of ions  (obsolete?)


!     check consistency 'nrow' <--> nr. of projectors


IF(ipsptyp == 1 .AND. nion > 1) THEN
!old     &   (ionmdtyp.ne.0 .or. icooltyp.ne.0)) then
  DO ion=1,nion
!old            if(np(ion).ne.np(1))
!old     &           stop ' all ions must be in the same row of elements'
    npact = np(ion)
    IF(ABS(h0_11g(npact))+ABS(h1_11g(npact)) > 0D0 .AND. nrow(npact) < 2)  &
        STOP "NROW must be at least 2 for proper projectors"
    IF(ABS(h0_22g(npact)) > 0D0 .AND. nrow(npact) < 3)  &
        STOP "NROW must be at least 3 for proper projectors"
    IF(ABS(h0_33g(npact))+ABS(h1_22g(npact)) +ABS(h2_11g(npact)) > 0D0  &
        .AND. nrow(npact) < 4)  &
        STOP "NROW must be at least 4 for proper projectors"
  END DO
END IF


!         if(nrare.gt.0 .and. nfix.eq.1)then
!            call fixion
!c if there is a bulk to simule
!         elseif(nfix.eq.0 .or. nrare.eq.0)then
!            do i=1,nion
!               nfixed(i)=0
!            enddo
!         endif




!       np(0) for use in monte-carlo


IF(icooltyp == 3)  CALL cenmass()
np(0) = np(1)
dt12=dt1




!       initialize peudopotential background

!g         call calcpseudo(rho)

IF(ipsptyp == 1) THEN
  DO ion=1,nion
    CALL calc_proj(cx(ion),cy(ion),cz(ion),ion)
  END DO
END IF


!     background energy


WRITE(6,*) 'Calculating ionic energy ...'

sumion = energ_ions()     !  ???
WRITE(7,'(a,f7.4)') 'sumion=',sumion
WRITE(6,'(a,f7.4)') 'sumion=',sumion
WRITE(7,*)
WRITE(6,*)
!k
ecorr=sumion           !analytical case   ??
!k


!     initialize common velocity

IF (ekin0pp > 0) THEN
  v0 = SQRT(2.*ekin0pp/(amu(np(nion))*1836.0*ame))
  rnorm = vxn0**2 + vyn0**2+ vzn0**2
  rnorm = SQRT(rnorm)
  IF (rnorm == 0) STOP 'Velocity vector not normalizable'
  vxn02 = vxn0/rnorm*v0*ame
  vyn02 = vyn0/rnorm*v0*ame
  vzn02 = vzn0/rnorm*v0*ame
  xm = amu(np(nion))*1836.0
  DO i=1,nion
    cpx(i)=vxn02*xm
    cpy(i)=vyn02*xm
    cpz(i)=vzn02*xm
  END DO
  WRITE(*,*) ' initial CP:',cpx(1),cpy(1),cpz(1)
END IF ! initial kinetic energy for ions


!     short protocol


IF(myn == 0) THEN
  WRITE(7,'(a)')'initial ionic positions and velocities:'
  WRITE(6,'(a)')'initial ionic positions and velocities:'
  DO ion=1,nion
    WRITE(7,'(a,i2,a,6f9.4)') 'ion',ion, '=',cx(ion),cy(ion),cz(ion),  &
        cpx(ion),cpy(ion),cpz(ion)
    WRITE(6,'(a,i2,a,6f9.4)') 'ion',ion, '=',cx(ion),cy(ion),cz(ion),  &
        cpx(ion),cpy(ion),cpz(ion)
  END DO
END IF




WRITE(6,*) 'Initions Done.'

RETURN
END SUBROUTINE initions


!-----init_jellium-------------------------------------------------

SUBROUTINE init_jellium()

!     initialize jellium background


USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

!------------------------------------------------------------------

IF(nion2 /= 0) RETURN              ! case of detailed ions

!     factors for transformation from hill-wheeler coordinates

gamarg = gamma * pi / 180D0
a20fac = COS(gamarg)
a22fac = SIN(gamarg) / SQRT(2.0)

!     transformation from hill-wheeler coordinates

falph = bbeta * a20fac
fbeta = bbeta * a22fac
fhexe = beta4
WRITE(7,'(/a,f8.3,a,f8.1,a,f8.3)')  &
    ' electronic input: beta=',bbeta,'  gamma=',gamma, '  beta4=',beta4



xclust=nion*1D0
CALL jelbak(xclust,falph,fbeta,fhexe,srms,0)

sum1=0D0
DO ind=1,nxyz
  sum1=sum1+rhojel(ind)
END DO
sum1=sum1*dvol
WRITE(6,*) 'jellium-norm=',sum1
WRITE(7,*) 'jellium-norm=',sum1

RETURN
END SUBROUTINE init_jellium

!-----initwf--------------------------------------------------------


SUBROUTINE initwf(psir)

!     master routine to initialize the electronic wavefunctions

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)
#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

REAL(DP), INTENT(IN OUT)                     :: psir(kdfull2,kstate)
REAL(DP) :: rho(kdfull2*2)

!----------------------------------------------------------------------


n=myn

!      h2m=hbar*hbar/2.0/ame
IF(myn == 0)THEN
  IF(nspdw > nclust/2) STOP 'nspdw must be less or equal nclust/2'
END IF


!     prepare book-keeping for wavefunctions

omeg=0.25*h2m
!      if(nclust.gt.0)
!     &    call ininqb(nclust,deocc,b2occ,gamocc*pi/180D0)
CALL ininqb(nclust,deocc,b2occ,gamocc*pi/180D0)


!     initialize H.O. wavefunctions

IF(init_lcao /= 1) THEN
#if(parano)
  IF (iuselastwf == 0) THEN
    CALL initho(psir)
  ELSE
!             stop
  END IF
#endif
#if(parayes)
    CALL initho(psir)
#endif
  
  
!       remove symmetry restrictions
  
!  IF(init_lcao /= 1) CALL rem_sym(psir)
  CALL rem_sym(psir)
  WRITE(7,'(a)')'after ordering:'
  WRITE(7,'(a,2i3)') 'nstate,nclust',nstate,nclust
  WRITE(6,'(a)')'after ordering:'
  WRITE(6,'(a,2i3)') 'nstate,nclust',nstate,nclust
  IF(ispinsep /= 0) CALL spinsep(psir)
  
  WRITE(7,'(a)') 'after oscillator initialization:'
  DO i=1,nstate
    nr=nrel2abs(i)
    IF (myn == 0) THEN
      WRITE(7,'(a,i3,a,f12.4,a,i2,a,3i3)')'wf ',i,  &
          ' occ',occup(i),' sp',ispin(nr), ' knots:',nq(1,nr),nq(2,nr),nq(3,nr)
    END IF
  END DO
  
ELSE
  
!       optionally LCGO initialization
  
  CALL genermowf(psir)
  WRITE(6,'(a)') 'after LCAO initialization:'
  WRITE(7,'(a)') 'after LCAO initialization:'
  DO nbr=1,nstate
    en=0D0
    DO i=1,kdfull2
      en=en+psir(i,nbr)*psir(i,nbr)
    END DO
    en=en*dx*dy*dz
    WRITE(6,*) 'node=',myn,', state=',nbr,  &
        ', en=',en,', occ=',occup(nbr),', spin=',ispin(nbr)
    WRITE(7,*) 'node=',myn,', state=',nbr,  &
        ', en=',en,', occ=',occup(nbr),', spin=',ispin(nbr)
  END DO
END IF




!     optional shift of wavefunctions

#if(simpara)
WRITE(7,*) ' SHIFTFIELD overridden'
#else
IF(ABS(shiftwfx)+ABS(shiftwfy)+ABS(shiftwfz) > 1D-20) THEN
  DO nbe=1,nstate
  
    DO ii=1,kdfull2
      rfieldtmp(ii)=psir(ii,nbe)
    END DO
  
    CALL shiftfield(rfieldtmp,shiftwfx,shiftwfy,shiftwfz)
  
    DO ii=1,kdfull2
      psir(ii,nbe)=rfieldtmp(ii)
    END DO
  
  END DO
END IF
#endif


!     final ortho-normalization to clean up

#if(parayes)
CALL  mpi_comm_rank(mpi_comm_world,myn,icode)
WRITE(*,*) ' wfs initialized: myn=',myn
CALL mpi_barrier (mpi_comm_world, mpi_ierror)
WRITE(6,*) 'myn=',myn,' before SCHMID'
#endif
CALL schmidt(psir)
#if(parayes)
WRITE(6,*) 'myn=',myn,' after SCHMID'
CALL mpi_barrier (mpi_comm_world, mpi_ierror)
#endif
RETURN
END SUBROUTINE initwf


!-----ininqb-------------------------------------------------------

SUBROUTINE ininqb(nelect,deoccin,betain,gamin)
USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)
#if(parayes)
INCLUDE 'mpif.h'

INTEGER, INTENT(IN)                      :: nelect
REAL(DP), INTENT(IN)                         :: deoccin
REAL(DP), INTENT(IN)                         :: betain
REAL(DP), INTENT(IN OUT)                     :: gamin
INTEGER :: is(mpi_status_size)
#endif

!     initialization of book-keeping arrays of states nq, ispin.

!     estimates fermi energy from particle number and considers
!     all states up "fermi-energy plus deoccin" in the ordering
!     of a deformed harmonic oscillator for given deformation.
!     the input parameters are:
!       nelect  = number of electrons
!       deoccin   = number of osc. shells above fermi shell
!       betain  = quadrupole deformation of jellium background
!       gamin   = triaxiality angle of jellium background
!       temp    = (via common) temperature, temp=0 cuts to occupied only

!     the output goes on the occuoation fields nq and ispin on
!     common /option/ .

INTEGER, PARAMETER :: kmxsav=kdfull/3
REAL(DP) :: esp(ksttot)           ! storage for s.p. energies
REAL :: efacto                ! factor to get energies from h.o.
REAL :: efermi                ! estimate for fermi shell
REAL :: q20fac                ! sqrt(5/16pi)
REAL :: cosfac,sinfac         ! weightes deduced from 'gamin'
!     real      xfac,yfac,zfac        ! effective osc. energies in x,y,z
REAL :: speact                ! actual s.p. energy in loop
INTEGER :: noscmx                ! maximum oscillator number
INTEGER :: n                     ! nr. of state
INTEGER :: noscx,noscy,noscz     ! osc. nr. in each direction

REAL(DP) :: ph(2*ksttot)             ! degeneracy of wavefunction, for
REAL(DP) :: occu(ksttot)
#if(fullspin)
DATA ph/ksttot*1D0,ksttot*1D0/
#else
DATA ph/ksttot*2.0/
#endif
LOGICAL :: tocc
DATA tocc/.false./

!-----------------------------------------------------------------------

!     prepare initial parameters
!     estimate of Fermi energy relies on spherical oscillator shells.
!     the shell label N is determined by

!     nelect,spin = N*(N+1)*(N+2)/6

!     where 'nelect,spin' is the nr. of electrons for given spin.
!     this relation is resolved approximately for N.

q20fac = SQRT(5.0/(16.0*pi))
#if(fullspin)
nelup  = nelect-nspdw
neldw  = nspdw
!test
WRITE(6,*) 'nelup',nelup
WRITE(6,*) 'neldw',nspdw
!test
#else
IF(MOD(nelect,2) == 1)  &
    STOP ' nr. of electrons must be even for spin degeneracy'
nelup  = nclust/2
#endif
efacto = 0.25/(1D0*nelect)**0.3333333
efrmup = (6.0*nelup)**0.3333333
efrmup = efrmup/(1D0-1D0/(efrmup*efrmup))**0.3333333-1.5
ecutup = efrmup+deoccin
nomxup = ecutup*(1D0+2.0*q20fac*betain)+0.5
ecutup = efacto*ecutup
#if(fullspin)
IF(neldw > 0) THEN
  efrmdw = (6.0*neldw)**0.3333333
  efrmdw = efrmdw/(1D0-1D0/(efrmdw*efrmdw))**0.3333333-1.5
ELSE
  efrmdw = -0.00001
END IF
ecutdw = efrmdw+deoccin
nomxdw = ecutdw*(1D0+2.0*q20fac*betain)+0.5
ecutdw = efacto*ecutdw
#endif
cosfac = q20fac*COS(gamin)
sinfac = q20fac*SQRT(2.0)*SIN(gamin)
xfac   = efacto/(1D0-betain*(cosfac-sinfac))
yfac   = efacto/(1D0-betain*(cosfac+sinfac))
zfac   = efacto/(1D0+2.0*betain*cosfac)
!      write(*,*) ' efacto,efrmup,ecutup,deoccin=',
!     &   efacto,efrmup,ecutup,deoccin

!     loop over all possible osc. triplets, determine s.p. energy
!     and keep if below ecut

!     select spin up states
!     distribute preliminary occupation numbers

WRITE(6,'(a,3f8.3)') ' zfac,yfac,xfac=',zfac,yfac,xfac
WRITE(6,'(a)') ' is  n  nx   ny   nz   speact'
nmaxdiff = (ksttot - nelect)/2
n = 0
lb1: DO noscz=0,nomxup
     DO noscy=0,nomxup
     DO noscx=0,nomxup
      speact = noscz*zfac+noscy*yfac+noscx*xfac
!        write(*,*) ' noscx,noscy,noscz,speact=',noscx,noscy,noscz,speact
      IF(speact < ecutup) THEN
        n     = 1+n
!          write(*,*) 'n=',n
!$$$          if(n.gt.ksttot)
!$$$     &      stop ' deoccin or part.number to large for given ksttot'
        nq(1,n)  = noscx
        nq(2,n)  = noscy
        nq(3,n)  = noscz
        ispin(n) = 1
        esp(n)   = speact
!test
!          write(6,*) 'n,nelup,ksttot',n,nelup,ksttot
!test
        IF(n <= nelup) THEN
          occu(n) = 1D0
        ELSE
          occu(n) = 0D0
        END IF
        
        IF(n > ksttot)  &
            STOP ' deoccin or part.number to large for given ksttot'
        
!          write(6,*)
!         write(6,*) 'speact',speact,'n',n,'isp',ispin(n)
        IF(n > (nelup+nmaxdiff)) THEN
          WRITE(6,*) 'n=',n,' greater than nelup+dnmax=',nelup+nmaxdiff
          n=n-1
          EXIT lb1
        END IF
!        WRITE(6,'(5i5,1f10.3)') ispin(n),n,noscx,noscy,noscz,speact
        WRITE(6,'(2i5,3f10.3)') ispin(n),n,occu(n),speact,ecutup
      END IF
     END DO
     END DO
     END DO lb1
!19   continue
!test
WRITE(6,*) 'n(nelup)',n
!test
IF(n < nelup) STOP ' not enough states to reach (spin-up) particle number'
#if(fullspin)

!     select spin down states
!     distribute preliminary occupation numbers

mspindw = 0
lb2: DO noscz=0,nomxdw
     DO noscy=0,nomxdw
     DO noscx=0,nomxdw
      speact = noscz*zfac+noscy*yfac+noscx*xfac
      IF(nelup /= neldw) speact = speact*1.01    ! enhance for spin asymmetry
      IF(speact <= ecutdw) THEN
        n     = 1+n
!test
!        WRITE(6,*) 'n,neldw,ksttot',n,neldw,ksttot
!test
!GB          if(n.gt.ksttot)
        nq(1,n)  = noscx
        nq(2,n)  = noscy
        nq(3,n)  = noscz
        ispin(n) = 2
        mspindw = mspindw + 1
        IF(mspindw <= neldw) THEN
          occu(n) = 1D0
        ELSE
          occu(n) = 0D0
        END IF
        esp(n)   = efacto*speact
        
        IF(n > ksttot) STOP ' deocc or part.number to large for given ksttot'
        
        IF(mspindw > (neldw+nmaxdiff)) THEN
           WRITE(6,*) 'mspindw=',mspindw,&
                ' greater than neldw+dnmax=',neldw+nmaxdiff
!          WRITE(6,*) 'n greater that neldw+dnmax'
           mspindw=mspindw-1
          EXIT lb2
        ENDIF
!        WRITE(6,'(5i5,1f10.3)') ispin(n),n,noscx,noscy,noscz,speact
        WRITE(6,'(3i5,3f10.3)') ispin(n),n,mspindw,occu(n),speact,ecutdw
        WRITE(7,*)
        WRITE(7,*) 'speact',speact,'n',n,'isp',ispin(n)
        WRITE(7,*) 'no',noscx,noscy,noscz,zfac,yfac,xfac
      END IF
     END DO
     END DO
     END DO lb2
!29   CONTINUE
!test
WRITE(6,*) 'mspindw',mspindw
WRITE(6,*) 'neldw',neldw
!test
IF(mspindw < neldw)  &
    STOP ' not enough states to reach spin-down particle number'
#endif

nstate=n
WRITE(7,*) 'nstate ininqb',nstate
!      if(iforce.eq.1) ispin(nstate) = ispin(nstate-1)
IF(iforce == 1) STOP "option IFORCE obsolete"

!    preliminary occupation numbers from estimated s.p. energies

!      eferm  = 0.0         ! restart search of eferm from scratch
!      epstmp = 1.e-8
!      nnumb = 0
!      nph   = ph(1)
!      do i=1,nstate
!        occu(i)=1.0
!        nnumb = nph + nnumb
!        write(7,'(t2, 3i2, tr5, f9.3, tr3, i2)')
!     &      (nq(k,i), k=1,3), occu(i), ispin(i)
!      enddo
IF(tocc .AND. nstate > nelect)  &
    CALL pair(esp,occu,ph,nclust,nstate,gp,eferm,temp,partnm,  &
    90,4,epstmp,-1,ksttot)

!     reorder states in reverse order of 'occu' (largest first)

DO i=1,nstate
  DO j=i+1,nstate
    IF(occu(j) > occu(i)) THEN
      sav      = occu(i)
      occu(i) = occu(j)
      occu(j) = sav
      sav    = esp(i)
      esp(i) = esp(j)
      esp(j) = sav
      isav    = ispin(i)
      ispin(i) = ispin(j)
      ispin(j) = isav
      DO k=1,3
        isav    = nq(k,i)
        nq(k,i) = nq(k,j)
        nq(k,j) = isav
      END DO
    END IF
  END DO
END DO
DO i=1,nstate
  occu(i) = ph(i)*occu(i)
  WRITE(7,'(t2, 3i2, tr5, e9.3, tr5, f9.3, tr3, i2)')  &
      (nq(k,i), k=1,3), occu(i), esp(i), ispin(i)
END DO
!      occu(nstate)=0.0


!-----------------------------------------------------------------------

!  initial parameters:


!    parallel version : knode=number of nodes
!    in the rest we mean absolute = relative to all the wf
!                        relative = relative to the wf on our node

!    nstate=number of wavefunctions on our node
!       (in the nonparallel case : all the wavefunctions !)
!    nrel2abs=array  to get the absolute index of the wf
!       (not the one on the node !)
!    nabs2rel=array  to get the relative index of the wf
!       ( on the node !)
!    nhome=home of the wf (absolute !)
!    ispin=array of spins as a function of relative index:
!          up and down-wfs are alternated
!    occup=array of occupation as a function of  relative index
!          -> always 1.0 in the present code

#if(parano)
DO i=1,nstate
  nrel2abs(i)=i
  nabs2rel(i)=i
  nhome(i)=0
  occup(i)=occu(i)
END DO
nstate_all = nstate

!     if uneven number : better enforce by hand a zero occupation
!     number for the last wf, with an even nstate !

#endif
#if(parayes)
CALL  mpi_comm_rank(mpi_comm_world,myn,icode)
!old      nsttest=0
!old      nloc=0
!old      do ipr=0,knode-1
!old         nprov=0
!old         do ifull=1,kstate
!old            nsttest=nsttest+1
!old            if(nsttest.le.nstate) nprov=ifull
!old            if(myn.eq.ipr)  nloc=nprov
!old            write(30+myn,*)
!old     &           'ipr',ipr,'ifull',ifull,'nprov',nprov,'nloc',nloc
!old         enddo
!old      enddo
!old      nstate=nloc
!old      write(6,*) 'changing  on  proc',myn,'nstate to',nstate
!old      write(7,*) 'changing  on  proc',myn,'nstate to',nstate
!old      do i=1,kstate
!old        nrel2abs(i)=i+myn*kstate
!old      enddo
!old      do i=1,ksttot
!old        nabs2rel(i)=mod(i-1,kstate)+1
!old        nhome(i)=(i-1)/kstate
!old      enddo
!old      do i=1,nstate
!old        occup(i)=occu(nrel2abs(i))
!old      enddo
!old      write(7,'(8i1)') (ispin(nrel2abs(i)),i=1,nstate)
!old      write(7,'(8i1)')   (nhome(i),i=1,ksttot)

!     compute equidistribution of wfs on nodes

nstate_all = nstate
nstpernode = nstate_all/knode
nodeplus   = nstate_all-nstpernode*knode
nstaccum = 0
DO nod=0,knode-1
  IF(nod < nodeplus) THEN
    nstate_node(nod) = nstpernode+1
  ELSE
    nstate_node(nod) = nstpernode
  END IF
  nstart_node(nod) = nstaccum
  nstaccum = nstaccum+nstate_node(nod)
END DO
nstate = nstate_node(myn)
WRITE(6,'(a,i3,a,2i5)') 'myn=',myn,  &
    ': nstaccum,nstate_all=',nstaccum,nstate_all
WRITE(6,'(a,100i3)') ' nstate_node=',(nstate_node(nod),nod=0,knode-1)
WRITE(6,'(a,100i3)') ' nstart_node=',(nstart_node(nod),nod=0,knode-1)
IF(nstaccum /= nstate_all) STOP '  number of states does not match '

DO i=1,kstate
  nrel2abs(i)=i+nstart_node(myn)
END DO
DO i=1,ksttot
  DO nod=knode-1,0,-1
    IF(i > nstart_node(nod)) THEN
      nabs2rel(i)=i-nstart_node(nod)
      nhome(i)=nod
      ispin_node(nabs2rel(i),nod) = ispin(i)
      GO TO 999
    END IF
  END DO
  999    CONTINUE
END DO
DO i=1,nstate
  occup(i)=occu(nrel2abs(i))
END DO
WRITE(7,'(a,i3,a,80i1)') 'myn=',myn,': ispin',(ispin(nrel2abs(i)),i=1,nstate)
WRITE(6,'(a,i3,a,80i1)') 'myn=',myn,': ispin',(ispin(nrel2abs(i)),i=1,nstate)
WRITE(7,'(a,i3,a,80i1)') 'myn=',myn,': nhome',(nhome(i),i=1,ksttot)
WRITE(6,'(a,i3,a,80i1)') 'myn=',myn,': nhome',(nhome(i),i=1,ksttot)
#endif
RETURN
END SUBROUTINE ininqb


!-----jelbak-----------------------------------------------------jelbak

SUBROUTINE jelbak(partn,alphel,betael,hexel,sqr,iturn)

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

!     computes the jellium background such that the radius is
!     given by 'radjel'*part.numb.**1/3, the surface by
!     'surjel', and the deformation is adjusted to the
!     electron deformation given on 'alphel' and 'betael'

!     y40 is expressed in terms of y20_effective to turn it automaticall
!     into the principal axes.

!     this version allows three-dimensional rotation of the
!     jellium background if iturn is one. this happens as an excitation
!     mechanism shortly before the dynamics is started.




REAL(DP), INTENT(IN)                         :: partn
REAL(DP), INTENT(IN OUT)                     :: alphel
REAL(DP), INTENT(IN OUT)                     :: betael
REAL(DP), INTENT(IN)                         :: hexel
REAL(DP), INTENT(OUT)                        :: sqr
INTEGER, INTENT(IN)                      :: iturn
DATA astep/0.6/

REAL(DP) :: vecin(3),vecout(3),vecalpha(3)

!----------------------------------------------------------------------

!     evaluate cases according to 'irotat'

IF(iturn == 1) THEN
!  case for scissor modes
  IF(irotat == 1) THEN
    vecalpha(1)  = phirot*PI/180D0
    vecalpha(2)  = 0D0
    vecalpha(3)  = 0D0
  ELSE IF(irotat == 2) THEN
    vecalpha(1)  = 0D0
    vecalpha(2)  = phirot*PI/180D0
    vecalpha(3)  = 0D0
  ELSE IF(irotat == 3) THEN
    vecalpha(1)  = 0D0
    vecalpha(2)  = 0D0
    vecalpha(3)  = phirot*PI/180D0
  ELSE IF(irotat == 4) THEN
    vecalpha(1)  = phirot*PI/180D0
    vecalpha(2)  = vecalpha(1)
    vecalpha(3)  = vecalpha(1)
  END IF
  
  WRITE(6,'(/a,3f8.3)')  &
      ' jellium background rotated by anglex,angley,anglez=', anglex,angley,anglez
  WRITE(7,'(/a,3f8.3)')  &
      ' jellium background rotated by anglex,angley,anglez=', anglex,angley,anglez
  
ELSE IF(iturn == 0) THEN
!  case for static iteration
  anglex  = 0D0
  angley  = 0D0
  anglez  = 0D0
END IF

!      itmax = itback

!     dimensionless deformation parameters
!      zero = 0.0
onetrd = 1D0 / 3.0
d4pi   = 1D0/(4.0*pi)
alpfac = 4.0*pi/5.0
q20fac = SQRT(5.0/(16.0*pi))
q22fac = SQRT(15.0/(8.0*pi))
q40fac = SQRT(9.0/(4.0*pi))
!----------------------------------------------------------------------

!     set initial values for iteration of deformation
alphbk = alphel
betabk = betael
hexabk = hexel


!     effective angle 'gamma' to determine principle axis for y40

bet2ef = (alphbk*alphbk+2.0*betabk*betabk)
IF(bet2ef /= 0D0) THEN
  bet2ef  = 1D0/SQRT(bet2ef)
  y20fac  = alphbk*bet2ef
  y20obs  = y20fac*q20fac
  y22fac  = betabk*bet2ef
  y22obs  = y22fac*q22fac
  q40red = 7.0*SQRT(pi/9.0)          ! factor for y40 from y20^2
  q20red = SQRT(5.0/pi)/7.0          ! cofactor on y20 in y40 from
  q40obs = 8.0*q40red/q40fac
END IF

radius = radjel * partn**onetrd
rhoc   = 3.0*d4pi / radjel**3

argum  = radius / surjel
IF(argum < 38.) THEN
  rho0 = rhoc * (1D0 + EXP(-argum))
ELSE
  rho0 = rhoc
END IF

!     iterate density to give the same deformation as the electrons
alpha  = alphel
beta   = betael


DO iter=1,itback
  
!       compute the distribution and accumulate its moments
  qoct  = 0D0
  sqhe   = 0D0
  sqt   = 0D0
  sqr   = 0D0
  dpm   = 0D0
  sqq   = 0D0
  sqn   = 0D0
  ii    = 0
  thetax = 0D0
  thetay = 0D0
  thetaz = 0D0
  theta0 = 0D0
  
  ii=0
  DO i3=minz,maxz
    zact = (i3-nzsh)*dz
    DO i2=miny,maxy
      yact = (i2-nysh)*dy
      DO i1=minx,maxx
        xact = (i1-nxsh)*dx
        ii    = ii+1
        IF(iturn == 1) THEN
!          CALL rotxyz(xact,yact,zact,x,y,z,anglex,angley,anglez)
          vecin(1) = xact
          vecin(2) = yact
          vecin(3) = zact
          CALL rotatevec3D(vecin,vecout,vecalpha)
          x = vecout(1)
          y = vecout(2)
          z = vecout(3)
        ELSE IF(iturn == 0) THEN
          x=xact
          y=yact
          z=zact
        END IF
        xx=x*x
        yy=y*y
        zz=z*z
        rr    = xx + yy + zz
        r     = SQRT(rr)
        IF(rr /= zero) THEN
          y20   = q20fac*(zz+zz-xx-yy)/rr
          y22   = q22fac*(xx-yy)/rr
          
          IF(bet2ef == 0D0) THEN
            y40   = q40fac*(8.0*zz*zz+3.0*xx*xx+ 3.0*yy*yy-24.0*zz*xx-  &
                24.0*zz*yy+6.0*xx*yy)/ (8.0*rr*rr)
          ELSE
            y2eff = y20fac*y20+y22fac*y22
            y40   = q40red*(y2eff*y2eff-q20red*y2eff-d4pi)
          END IF
        ELSE
          y20   = 0D0
          y22   = 0D0
          y40   = 0D0
        END IF
        reff = radius*(1D0+alphbk*y20+betabk*y22 +hexabk*y40  &
            -(alphbk*alphbk+2.0*betabk*betabk+hexabk*hexabk)*d4pi)
        argum = (r-reff)/surjel
        IF(argum > +38.0) THEN
          rhojel(ii) = 0D0
        ELSE IF(argum < -38.0) THEN
          rhojel(ii) = rho0
        ELSE
          rhojel(ii) = rho0/(1D0+EXP(argum))
        END IF
        rh    = rhojel(ii)
        sqn   = sqn   + rh
        sqr   = sqr   + rr*rh
        y20   = (zz+zz-xx-yy)
        sqq   = sqq   + y20*rh
        y22   = (xx-yy)
        sqt   = sqt   + y22*rh
        y2eff  = y20obs*y20+y22obs*y22
        y40    = q40obs*(y2eff*y2eff-q20red*y2eff*rr-d4pi*rr*rr)
        sqhe   = sqhe + y40*rh
        thetax = thetax  + (yy+zz)*rh
        thetay = thetay  + (xx+zz)*rh
        thetaz = thetaz  + (xx+yy)*rh
      END DO
    END DO
  END DO
  
  volel  = dx*dy*dz
  sqn    = volel*sqn
  dpm    = volel*dpm
  sqr    = SQRT(volel*sqr/sqn)
  sqq    = volel*sqq
  qoct   = volel*qoct
  sqt    = volel*sqt
  sqhe   = volel*sqhe
  theta0 = (thetax+thetay+thetaz)/3.0
  thetax = thetax/theta0
  thetay = thetay/theta0
  thetaz = thetaz/theta0
  alpha  = alpfac*q20fac*sqq/(sqn*sqr*sqr)
  beta   = alpfac*0.5*q22fac*sqt/(sqn*sqr*sqr)
  deralp = 1D0
  derbet = 1D0
  delalp = -astep*(alpha-alphel)/deralp
  delbet = -astep*(beta -betael)/derbet
  alphbk = alphbk+delalp
  betabk = betabk+delbet
  precis = (alpha-alphel)**2+(beta-betael)**2+(partn-sqn)**2
  
  IF((precis < endcon).AND.(ABS(partn-sqn) < endcon)) GO TO 199
  
  radius = radius * (partn/sqn)**onetrd
  argum  = radius / surjel
  IF(argum < 38.) THEN
    rho0 = rhoc * (1D0 + EXP(-argum))
  ELSE
    rho0 = rhoc
  END IF
!      write(7,1000)    iter,precis,sqn,dpm,sqr,sqq,sqt,
!     &        sqhe,alpha,beta,(alphbk-delalp),(betabk-delbet),alphel,
!     &            betael,deralp,derbet,delalp,delbet
!      write(7,*) 'thetax=',thetax
  
END DO

WRITE(7,'(a/a,g11.3)') ' ---> background deformation did not converge!',  &
    '      residual error in"alpha"=',precis
STOP ' no convergence in "jelbak"'

199  CONTINUE

WRITE(7,'(a,i4,a,f12.7,a,2(/3(a,f12.7)),5(/2(a,f12.7)))')  &
    ' iteration nr.',iter,' precision=',precis,':',  &
    ' sqn   =',sqn,  ' dpm   =',dpm,  ' sqr   =',sqr,  &
    ' sqq   =',sqq,  ' sqt   =',sqt,  ' sqhe  =',sqhe,  &
    ' alpha =',alpha,' beta  =',beta,  &
    ' alphbk=',alphbk-delalp,' betabk=',betabk-delbet,  &
    ' alphel=',alphael,' betael=',betael, ' deralp=',deralp,' derbet=',derbet,  &
    ' delalp=',delalp,' delbet=',delbet
WRITE(7,*) 'thetax=',thetax
WRITE(7,*) 'thetay=',thetay
WRITE(7,*) 'thetaz=',thetaz
beta2j = SQRT(alpha*alpha+2.0*beta*beta)
gammaj = ATAN(1.4142136*beta/alpha)*180D0/pi
WRITE(7,*) 'effective jellium deformations: beta2j=',beta2j,  &
    '  gammaj=',gammaj


RETURN
END SUBROUTINE jelbak

!-----initho------------------------------------------------------initho

SUBROUTINE initho(psir)

!     initailizes harmonic oscillator wavefunctions

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)


REAL(DP), INTENT(OUT)                        :: psir(kdfull2,kstate)
REAL(DP) :: valx(nx1),valy(ny1),valz(nz1)
REAL(DP), ALLOCATABLE :: phix(:)

!EQUIVALENCE (phix(1),w1(1))

!----------------------------------------------------------------------

!     check workspace

ALLOCATE(phix(kdfull))
!IF(usew1) STOP ' in SSTEP: workspace W1 already active '
!usew1 = .true.

an  = REAL(2*nclust)
IF(temp > 0D0) THEN
  homx = omeg*an**(-0.33333333)*xfac
  homy = omeg*an**(-0.33333333)*yfac
  homz = omeg*an**(-0.33333333)*zfac
  bxx  = (2.0*h2m/homx)**3
  bxy  = (2.0*h2m/homy)**3
  bxz  = (2.0*h2m/homz)**3
  bxx  = osfac*(bxx**0.16666666)
  bxy  = osfac*(bxy**0.16666666)
  bxz  = osfac*(bxz**0.16666666)
ELSE
  hom  = omeg*an**(-0.33333333)
  bk1  = (2.0*h2m/hom)**3
  bk1   = osfac*(bk1**0.166666666)
END IF


DO nb=1,nstate
  
!       number of knots
  
  
!       nq is  relative to the proc
  
!       this way we implicitely parallelize the computation of phix
  
  inx=nq(1,nrel2abs(nb))
  iny=nq(2,nrel2abs(nb))
  inz=nq(3,nrel2abs(nb))
  
  DO  i=1,3
    ipar(i,nb)=1-2*MOD(nq(i,nrel2abs(nb)),2)
  END DO
  
  IF(temp > 0D0) THEN
    CALL clust(inx,bxx,0D0,xval,valx,nx1)
    CALL clust(iny,bxy,0D0,yval,valy,ny1)
    CALL clust(inz,bxz,0D0,zval,valz,nzi)
  ELSE IF(temp <= 0D0) THEN
    CALL clust(inx,bk1,0D0,xval,valx,nx1)
    CALL clust(iny,bk1,0D0,yval,valy,ny1)
    CALL clust(inz,bk1,0D0,zval,valz,nz1)
  END IF
  
!       composition of the factorised wave-function
!       occupies only upper box (1/8 part), but is returned on 'psir'
  
  ii=0
  DO iz=1,nzi
    vz=valz(iz)
    DO iy=1,ny1
      vy=valy(iy)
      DO ix=1,nx1
        vx=valx(ix)
        ii=ii+1
        psir(ii,nb)=vx*vy*vz
      END DO
    END DO
  END DO
  
  
END DO

DEALLOCATE(phix)
!usew1 = .false.

RETURN
END SUBROUTINE initho


!-----clust-------------------------------------------------------

SUBROUTINE clust(in,b,z,x,val,n1)


!     initializes one part (x-,y- or z-direction) of a factorised
!     wave-function
!     input:  in    number of knots
!             b     inverse width of harmonic oscillator solution
!             z     location of the nucleus resp. cluster (see initho)
!             x     value of grid-point
!             n1    number of grid-points
!     output: val   part of the factorised wave-function

!     the function is
!                sqrt(2)^n*h(sqrt(2)*x)
!     where h(x) is the hermite polynomial as explained in rottmann, p.1

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)


INTEGER, INTENT(IN OUT)                  :: in
REAL(DP), INTENT(IN)                         :: b
REAL(DP), INTENT(IN OUT)                     :: z
REAL(DP), INTENT(IN OUT)                     :: x(ksmax)
REAL(DP), INTENT(OUT)                        :: val(ksmax)
INTEGER, INTENT(IN)                      :: n1


!-----------------------------------------------------------------------

fak = 1
DO i=2,in
  fak = fak * i
END DO
coef=1./SQRT(b*1.772454*(2.**in)*fak)


DO j=1,n1
  
  tau=(x(j)-z)/b
  argum = -0.5*tau*tau
  IF(argum < -38D0) THEN
    val(j) = 0D0
  ELSE
    v=coef*EXP(argum)
    IF(in == 0) THEN
      val(j)=v
    ELSE IF(in == 1) THEN
      val(j)=v*2*tau
    ELSE IF(in == 2) THEN
      val(j)=v*(4.*tau*tau-2.)
    ELSE IF(in == 3) THEN
      val(j)=v*(8.*tau**3-12.*tau)
    ELSE IF(in == 4) THEN
      val(j)=v*(16.*tau**4-48.*tau**2+12.)
    ELSE IF(in == 5) THEN
      val(j)=v*(32.*tau**5-160D0*tau**3+120D0*tau)
    ELSE IF(in == 6) THEN
      val(j)=v*8.0*(8.*tau**6-60D0*tau**4+90D0*tau**2-15.0)
    ELSE IF(in == 7) THEN
      val(j)=v*16.0*(8.*tau**7-84.*tau**5+210D0*tau**3-105.0*tau)
    ELSE
      WRITE(6,'(a,1pg12.4)') ' wrong radial quantum number in clust: in=',in
      STOP 'wrong radial quantum number in clust'
    END IF
  END IF
  
END DO

RETURN
END SUBROUTINE clust

!-----rotxyz-----------------------------------------------------------

SUBROUTINE rotxyz(xin,yin,zin,xout,yout,zout,anglex,angley,anglez)

!     rotation in three dimensions in three separate steps,
!      first, a rotation by an angle 'anglex' about the x-axis,
!      second, a rotation by an angle 'angley' about the y-axis, and
!      third, a rotation by an angle 'anglez' about the z-axis.
!     (xin,yin,zin) is the input vector and
!     (xout,yout,zout) the result.

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

!----------------------------------------------------------------------

!     rotation about x-axis, intermediate result on x1,y1,z1

x1     = xin
cphi   = COS(anglex)
sphi   = SIN(anglex)
y1     = yin*cphi+zin*sphi
z1     = -yin*sphi+zin*cphi

!     rotation about y-axis, intermediate result on x2,y2,z2

y2     = y1
cphi   = COS(angley)
sphi   = SIN(angley)
z2     = z1*cphi+x1*sphi
x2     = -z1*sphi+x1*cphi

!     rotation about z-axis, final result on xout,yout,zout

zout   = z2
cphi   = COS(anglez)
sphi   = SIN(anglez)
xout   = x2*cphi+y2*sphi
yout   = -x2*sphi+y2*cphi

RETURN
END SUBROUTINE rotxyz


!-----rem_sym----------------------------------------------------------

SUBROUTINE rem_sym(psi)

!     takes symmetric cube (upper 1/8) stored in real wavefunction 'psi'
!     and unfolds to full box using given parities


USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: psi(kdfull2,kstate)

REAL(DP), ALLOCATABLE :: phix(:)
!EQUIVALENCE (phix,w1)

!       dimension occu(kstate),ipa(3,kstate),ispi(kstate)

!     check workspace

ALLOCATE(phix(kdfull))
!IF(usew1) THEN
!  STOP ' in REM_SYM: workspace W1 already in use'
!ELSE
!  usew1 = .true.
!END IF

!      nxyz=nx1*ny1*nzi
DO nbe=1,nstate
  ii    = 0
  si=0
  DO i3=1,nzi
    DO i2=1,ny1
      DO i1=1,nx1
        ii    = ii+1
        y=psi(ii,nbe)
        rh=y*y
        IF((i1 /= 1).AND.(i1 /= nx1)) rh=rh+rh
        IF((i2 /= 1).AND.(i2 /= ny1)) rh=rh+rh
        IF((i3 /= 1).AND.(i3 /= nzi)) rh=rh+rh
        si   = si+rh
      END DO
    END DO
  END DO
  si=si
  WRITE(7,'(a,i3,a,f9.4)') 'level',nbe,'sum',si
!        do i=1,kdfull2
!          psi(i,nbe)=0.01
!        enddo
END DO

!     centering the density rhod on nx1 ny1 nzi
!     taking symmetry into account

nxyc=nx1*ny1
!      nxyf=nx2*ny2
!      nyf=nx2
!      nx=nx1-1
!      ny=ny1-1
!      nz=nzi-1

!     loop on states

DO nb=1,nstate
  
!     wf parity vs x y z plans
  
  xp=ipar(1,nb)
  yp=ipar(2,nb)
  zp=ipar(3,nb)
  ind = 0
  DO i3=1,nzi
    DO i2=1,ny1
      DO i1=1,nx1
        ind = 1 + ind
        phix(ind) = psi(ind,nb)
      END DO
    END DO
  END DO
  DO i3=1,nzi
    DO i2=1,ny1
      DO i1=1,nx1
        
!         init
        
        indsym=((i3-1)*nxyc+(i2-1)*nx1+i1)
!          pjel=poty(indsym)
        psp=phix(indsym)
        
!         direct part
        
        i1f=i1+nx-1
        i2f=i2+ny-1
        i3f=i3+nz-1
        ind=(i3f-1)*nxyf+(i2f-1)*nyf+i1f
        psi(ind,nb)=psp
!          if((i1.eq.1).and.(i2.eq.1).and.(i3.eq.1)) goto 10
        IF((i1 /= 1).OR.(i2 /= 1).OR.(i3 /= 1)) THEN
          
!           symmetry
          
          IF(i1 /= nx1) THEN
            i1f=nx-i1+1
            i2f=i2+ny-1
            i3f=i3+nz-1
            ind=(i3f-1)*nxyf+(i2f-1)*nyf+i1f
            psi(ind,nb)=psp*xp
          END IF
          
!           symmetry
          
          IF(i2 /= ny1) THEN
            i1f=i1+nx-1
            i2f=ny-i2+1
            i3f=i3+nz-1
            ind=(i3f-1)*nxyf+(i2f-1)*nyf+i1f
            psi(ind,nb)=psp*yp
          END IF
          
!           symmetry
          
          IF(i3 /= nzi) THEN
            i1f=i1+nx-1
            i2f=i2+ny-1
            i3f=nz-i3+1
            ind=(i3f-1)*nxyf+(i2f-1)*nyf+i1f
            psi(ind,nb)=psp*zp
          END IF
          
!           symmetry
          
          IF(i1 /= nx1) THEN
            IF(i2 /= ny1) THEN
              i1f=nx-i1+1
              i2f=ny-i2+1
              i3f=i3+nz-1
              ind=(i3f-1)*nxyf+(i2f-1)*nyf+i1f
              psi(ind,nb)=psp*xp*yp
            END IF
          END IF
          
!           symmetry
          
          IF(i2 /= ny1) THEN
            IF(i3 /= nzi) THEN
              i1f=i1+nx-1
              i2f=ny-i2+1
              i3f=nz-i3+1
              ind=(i3f-1)*nxyf+(i2f-1)*nyf+i1f
              psi(ind,nb)=psp*yp*zp
            END IF
          END IF
          
!           symmetry
          
          IF(i1 /= nx1) THEN
            IF(i3 /= nzi) THEN
              i1f=nx-i1+1
              i2f=i2+ny-1
              i3f=nz-i3+1
              ind=(i3f-1)*nxyf+(i2f-1)*nyf+i1f
              psi(ind,nb)=psp*xp*zp
            END IF
          END IF
          
!           symmetry
          
          IF(i1 /= nx1) THEN
            IF(i2 /= ny1) THEN
              IF(i3 /= nzi) THEN
                i1f=nx-i1+1
                i2f=ny-i2+1
                i3f=nz-i3+1
                ind=(i3f-1)*nxyf+(i2f-1)*nyf+i1f
                psi(ind,nb)=psp*xp*yp*zp
              END IF
            END IF
          END IF
        END IF
! 10       continue
      END DO
    END DO
  END DO
  
END DO

!     release workspace
!usew1 = .false.
DEALLOCATE(phix)

RETURN
END SUBROUTINE rem_sym

!-----spinsep-----------------------------------------------------------

SUBROUTINE spinsep(psir)

!     to induce local spin current by shifting spin-up and down in
!     different directions

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)


REAL(DP), INTENT(OUT)                        :: psir(kdfull2,kstate)


!----------------------------------------------------------------------

#if(fullspin)
DO nb=1,nstate
  sgeps = (3-2*ispin(nb))*0.01
  ii  = 0
  DO iz=minz,maxz
    z1=(iz-nzsh)*dz
    DO iy=miny,maxy
      y1=(iy-nysh)*dy
      DO ix=minx,maxx
        x1=(ix-nxsh)*dx
        ii = ii+1
        psir(ii,nb) = psir(ii,nb)*(1D0+sgeps*(x1+y1+z1))
      END DO
    END DO
  END DO
END DO
#endif
RETURN
END SUBROUTINE spinsep

!-----fixion------------------------------------------------------------

SUBROUTINE fixion

!     to fixe two layers along y and z of an argon cluster (only cores)
!     to simule the bulk
!     the input file must have already the positions of the bulk

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER :: icy(nion),icz(nion)
INTEGER :: nfixedyu(nion),nfixedyd(nion)
INTEGER :: nfixedzu(nion),nfixedzd(nion)

! two layers along y and z are fixed
! first layer
ind = 0
nfixedyu(1) = 1
nfixedyd(1) = 1
nfixedzu(1) = 1
nfixedzd(1) = 1
icy(1) = cy(1)
IF(ABS(cy(1)-icy(1)) > 0.5D0 .AND. cy(1) < 0D0) icy(1) = icy(1) - 1
IF(ABS(cy(1)-icy(1)) > 0.5D0 .AND. cy(1) > 0D0) icy(1) = icy(1) + 1
icz(1) = cz(1)
IF(ABS(cz(1)-icz(1)) > 0.5D0 .AND. cz(1) < 0D0) icz(1) = icz(1) - 1
IF(ABS(cz(1)-icz(1)) > 0.5D0 .AND. cz(1) > 0D0) icz(1) = icz(1) + 1
DO ion=2,nrare
  icy(ion) = cy(ion)
  IF(ABS(cy(ion)-icy(ion)) > 0.5D0 .AND. cy(ion) < 0D0) icy(ion) = icy(ion) - 1
  IF(ABS(cy(ion)-icy(ion)) > 0.5D0 .AND. cy(ion) > 0D0) icy(ion) = icy(ion) + 1
  DO ind1=1,ion-1
    IF(nfixedyu(ind1) > 0)THEN
      IF(icy(ion) > icy(ind1))THEN
        nfixedyu(ion) = ion
        nfixedyu(ind1) = 0
      ELSE IF(icy(ion) == icy(ind1))THEN
        nfixedyu(ion) = ion
      ELSE
        nfixedyu(ion) = 0
      END IF
    END IF
    IF(nfixedyd(ind1) > 0)THEN
      IF(icy(ion) < icy(ind1))THEN
        nfixedyd(ion) = ion
        nfixedyd(ind1) = 0
      ELSE IF(icy(ion) == icy(ind1))THEN
        nfixedyd(ion) = ion
      ELSE
        nfixedyd(ion) = 0
      END IF
    END IF
  END DO
END DO
DO ion=2,nrare
  icz(ion) = cz(ion)
  IF(ABS(cz(ion)-icz(ion)) > 0.5D0 .AND. cz(ion) < 0D0) icz(ion) = icz(ion) - 1
  IF(ABS(cz(ion)-icz(ion)) > 0.5D0 .AND. cz(ion) > 0D0) icz(ion) = icz(ion) + 1
  DO ind1=1,ion-1
    IF(nfixedzu(ind1) > 0)THEN
      IF(icz(ion) > icz(ind1))THEN
        nfixedzu(ion) = ion
        nfixedzu(ind1) = 0
      ELSE IF(icz(ion) == icz(ind1))THEN
        nfixedzu(ion) = ion
      ELSE
        nfixedzu(ion) = 0
      END IF
    END IF
    IF(nfixedzd(ind1) > 0)THEN
      IF(icz(ion) < icz(ind1))THEN
        nfixedzd(ion) = ion
        nfixedzd(ind1) = 0
      ELSE IF(icz(ion) == icz(ind1))THEN
        nfixedzd(ion) = ion
      ELSE
        nfixedzd(ion) = 0
      END IF
    END IF
  END DO
END DO
DO ion=1,nion
  IF(ion <= nrare)THEN
    IF(nfixedyu(ion) > 0 .OR. nfixedyd(ion) > 0 .OR. nfixedzu(ion) > 0  &
          .OR. nfixedzd(ion) > 0)THEN
      nfixed(ion) = ion
      ind = ind +1
    ELSE
      nfixed(ion) = 0
    END IF
  ELSE
    nfixed(ion) = 0
  END IF
END DO

! second layer
ind2 = 1
DO WHILE(nfixed(ind2) > 0)
  ind2 = ind2 + 1
END DO
nfixedyu(ind2) = ind2
nfixedyd(ind2) = ind2
nfixedzu(ind2) = ind2
nfixedzd(ind2) = ind2
DO ion=ind2+1,nrare
  IF(ion /= nfixed(ion))THEN
    DO ind1=ind2,ion-1
      IF(ind1 /= nfixed(ind1))THEN
        IF(nfixedyu(ind1) > 0)THEN
          IF(icy(ion) > icy(ind1))THEN
            nfixedyu(ion) = ion
            nfixedyu(ind1) = 0
          ELSE IF(icy(ion) == icy(ind1))THEN
            nfixedyu(ion) = ion
          ELSE
            nfixedyu(ion) = 0
          END IF
        END IF
        IF(nfixedyd(ind1) > 0)THEN
          IF(icy(ion) < icy(ind1))THEN
            nfixedyd(ion) = ion
            nfixedyd(ind1) = 0
          ELSE IF(icy(ion) == icy(ind1))THEN
            nfixedyd(ion) = ion
          ELSE
            nfixedyd(ion) = 0
          END IF
        END IF
      END IF
    END DO
  END IF
END DO
DO ion=ind2+1,nrare
  IF(ion /= nfixed(ion))THEN
    DO ind1=ind2,ion-1
      IF(ind1 /= nfixed(ind1))THEN
        IF(nfixedzu(ind1) > 0)THEN
          IF(icz(ion) > icz(ind1))THEN
            nfixedzu(ion) = ion
            nfixedzu(ind1) = 0
          ELSE IF(icz(ion) == icz(ind1))THEN
            nfixedzu(ion) = ion
          ELSE
            nfixedzu(ion) = 0
          END IF
        END IF
        IF(nfixedzd(ind1) > 0)THEN
          IF(icz(ion) < icz(ind1))THEN
            nfixedzd(ion) = ion
            nfixedzd(ind1) = 0
          ELSE IF(icz(ion) == icz(ind1))THEN
            nfixedzd(ion) = ion
          ELSE
            nfixedzd(ion) = 0
          END IF
        END IF
      END IF
    END DO
  END IF
END DO
DO ion=1,nion
  IF(nfixed(ion) == 0)THEN
    IF(ion <= nrare)THEN
      IF(nfixedyu(ion) > 0 .OR. nfixedyd(ion) > 0 .OR. nfixedzu(ion) > 0  &
            .OR. nfixedzd(ion) > 0)THEN
        nfixed(ion) = ion
        ind = ind +1
      ELSE
        nfixed(ion) = 0
      END IF
    ELSE
      nfixed(ion) = 0
    END IF
  END IF
END DO
WRITE(6,*)'number of atoms fixed:',ind

RETURN
END SUBROUTINE fixion




!-----checkoptions-------------------------------------------------


SUBROUTINE checkoptions()

!     to check consistency of compile-time options

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

!------------------------------------------------------------------

#if(parayes)
#if(findiff|numerov)
STOP ' parallele computing not active with finite differences'
#endif
#if(exchange)
STOP ' exact exchange not compatible with parallele code'
#endif
!#if(escmaskorb)
!STOP " parayes and escmaskorb incompatible"
!#endif
IF(jmp /=0) STOP ' evaluation of PES not possible in parallele code'
!IF(jescmaskorb /=0) &
! STOP ' evaluation orbital PAD not possible in parallele code'
IF(istat /= 0) STOP ' static restart not possible in parallele code'
#endif

#if(kli && exchange)
STOP ' exact exchange and KLI not compatible '
#endif
#if(kli && symmcond)
STOP ' propagated symm.cond. and KLI not compatible '
#endif
#if(kli && !directenergy)
STOP ' KLI requires directenergy=1 in "define.h" '
#endif
#if(directenergy && !pw92)
STOP ' directenergy=1 requires Perdew&Wang functional '
#endif
#if(directenergy && exchange)
STOP ' directenergy=1 not yet prepared for exact exchange '
#endif
#if(vdw && !raregas)
STOP " set raregas=1 when using VdW"
#endif


RETURN
END SUBROUTINE checkoptions


!-----init_lastwf---------------------------------------------------

SUBROUTINE init_lastewf(psir,psi)

!     overwrites initialization by old wavefunction
!     saved in from previous run.

!     !! This routine requires reconsideration.
!        It is not yet suited for parallel computing, etc.

USE params
USE kinetic
USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP), INTENT(IN OUT)                     :: psir(kdfull2,kstate)
COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)

REAL(DP) :: rho(kdfull2)   ! ???

!------------------------------------------------------------------


IF (istopcont == 3) THEN
  iuselastwf=0
  isavelastwf=1
  iuselast=-1
END IF

IF (istopcont == 4) THEN
  iuselast=1
  iusecell=0
  iuselastwf=1
  isavelastwf=0
END IF

IF (iuselastwf /= 0) THEN
  
  IF (myn == 0) THEN
    OPEN(709,STATUS='old',FILE='for005wf.init')
    READ(709,*) nstate,ngrpoints
    IF (ngrpoints /= kdfull2) STOP 'in init: wrong dimension of box'
    
    DO ist=1,nstate
      DO ind=1,kdfull2
        READ(709,*) psir(ind,ist)
      END DO
    END DO
    
  END IF ! myn
  
  
  CALL calcrhor(rho,psir)
  
  CLOSE(709)
  
END IF


IF (istopcont == 1) THEN
  IF (NE > 0) THEN
    WRITE(*,*) ' ADJUSTDIP from INIT'
    CALL adjustdip(rho)
  END IF
  CALL infor(psir,rho)
  IF (myn == 0) THEN
    OPEN(508,POSITION='append', FILE='energies.res')
    WRITE(508,'(f10.2,1e20.11)') cz(1),energy-enerinfty
    CLOSE(508)
  END IF
  STOP
END IF



IF (istopcont == 4) THEN
  
  CALL tinit(psir,psi)
  CALL calcrho(rho,psi)
!         call adjustdip(rho)
  CALL calcpseudo()
  
#if(gridfft)
  CALL falr(rho,chpcoul,nx2,ny2,nz2,kdfull2)
#endif
#if(findiff|numerov)
  CALL solv_fft(rho,chpcoul,dx,dy,dz)
#endif
  
  CALL calclocal(rho,aloc)
  CALL calc_sic(rho,aloc,psi)
!         call getSPenergiesr(psi,rho,aloc,akv,it)
  CALL info(psi,rho,aloc,akv,0)
  IF (myn == 0) THEN
    OPEN(508,POSITION='append', FILE='energies.res')
    WRITE(508,'(f10.2,1e20.11)') cz(1),energy-enerinfty
    CLOSE(508)
  END IF
  STOP  ' termination invoked by iStopCont=4'
  
END IF


RETURN
END SUBROUTINE init_lastewf



!-----init_homfield-------------------------------------------------

SUBROUTINE init_homfield()

!     initialize a homogeneous electrocal field and
!     adds it to the background field 'potFixedIon'.

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

!------------------------------------------------------------------

sc=ehomx**2+ehomy**2+ehomz**2
sc=SQRT(sc)
ehomx=ehomx/sc
ehomy=ehomy/sc
ehomz=ehomz/sc
ind=0
DO iz=minz,maxz
  z1=(iz-nzsh)*dz
  DO iy=miny,maxy
    y1=(iy-nysh)*dy
    DO ix=minx,maxx
      x1=(ix-nxsh)*dx
      ind = ind + 1
      vhom = e2*(x1*ehomx+y1*ehomy+z1*ehomz)*ehom0
      potfixedion(ind)=potfixedion(ind)+vhom
    END DO
  END DO
END DO

RETURN
END SUBROUTINE init_homfield


!-----init_grid-----------------------------------------------------

SUBROUTINE init_grid(ak,akv)

!     initialize Coulomb solver, kinetic energy and other
!     grid properties

USE params
USE kinetic
USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(OUT)                     :: ak(kdfull2)
REAL(DP), INTENT(OUT)                        :: akv(kdfull2)



!------------------------------------------------------------------



!     check and initialize

IF(myn == 0)THEN
  
  WRITE(6,*) 'nxin,nyin,nzin:', nxin,nyin,nzin
  WRITE(7,*) 'nxin,nyin,nzin:', nxin,nyin,nzin
  IF(nxin /= nx2/2) STOP ' input NX does not matched compiled grid'
  IF(nyin /= ny2/2) STOP ' input NY does not matched compiled grid'
  IF(nzin /= nz2/2) STOP ' input NZ does not matched compiled grid'
  IF(nabsorb > MIN(nxin,nyin,nzin)/2) STOP " NABSO too large"
END IF
!old      if(kxmax.lt.nx1) then
!tab         stop ' error in parameter'
!tab      elseif(kymax.lt.ny1) then
!tab         stop ' error in parameter'
!tab      elseif(kzmax.lt.nzi) then
!tab         stop ' error in parameter'
!tab      endif

dvol=dx*dy*dz

xx1=-dx
DO ix=1,nx1
  xx1=xx1+dx
  xval(ix)=xx1
  xt2(ix)=xx1*xx1
END DO

xy1=-dy
DO iy=1,ny1
  xy1=xy1+dy
  yval(iy)=xy1
  yt2(iy)=xy1*xy1
END DO

xz1=-dz
DO iz=1,nzi
  xz1=xz1+dz
  zval(iz)=xz1
  zt2(iz)=xz1*xz1
END DO


!     init coulomb solver

#if(findiff|numerov)
CALL d3sinfinit (dx,dy,dz)
#else
#if(coufou || coudoub)
CALL intcoul(dx,dy,dz,nx2,ny2,nz2)
#endif
#endif

!     init kinetic energy array

#if(findiff)
CALL inv3p_ini(dt1)
#endif
#if(numerov)
CALL inv5p_ini(dt1)
#endif
#if(gridfft)
dkx=pi/(dx*REAL(nx))
dky=pi/(dy*REAL(ny))
dkz=pi/(dz*REAL(nz))

WRITE(7,*) 'h bar squared over two m electron',h2m
WRITE(7,*) ' testprint EYE=',eye
WRITE(7,*) ' dkx,dky,dkz=',dkx,dky,dkz

!     prepare k**2 and kinetic propagation factor in 3D momentum space

ind=0
DO i3=1,nz2
  IF(i3 >= (nz+1)) THEN
    zkz=(i3-nz2-1)*dkz
  ELSE
    zkz=(i3-1)*dkz
  END IF
  DO i2=1,ny2
    IF(i2 >= (ny+1)) THEN
      zky=(i2-ny2-1)*dky
    ELSE
      zky=(i2-1)*dky
    END IF
    DO i1=1,nx2
      IF(i1 >= (nx+1)) THEN
        zkx=(i1-nx2-1)*dkx
      ELSE
        zkx=(i1-1)*dkx
      END IF
      ind=ind+1
      ak(ind)=EXP(-eye*dt1*(zkx**2+zky**2+zkz**2)*h2m)
      akv(ind)=(zkx**2+zky**2+zkz**2)*h2m
    END DO
  END DO
END DO
!MB
!      do ind=1,kdfull2
!         akvhelp(ind)=akv(ind)
!      enddo
!MB/

!     prepare kinetic propagation factors in 1D momentum spaces

DO i3=1,nz2
  IF(i3 >= (nz+1)) THEN
    zkz=(i3-nz2-1)*dkz
  ELSE
    zkz=(i3-1)*dkz
  END IF
  akpropz(i3)=EXP(-eye*dt1*zkz**2*h2m)
END DO

DO i2=1,ny2
  IF(i2 >= (ny+1)) THEN
    zky=(i2-ny2-1)*dky
  ELSE
    zky=(i2-1)*dky
  END IF
  akpropy(i2)=EXP(-eye*dt1*zky**2*h2m)
END DO

DO i1=1,nx2
  IF(i1 >= (nx+1)) THEN
    zkx=(i1-nx2-1)*dkx
  ELSE
    zkx=(i1-1)*dkx
  END IF
  akpropx(i1)=EXP(-eye*dt1*zkx**2*h2m)
END DO

#endif

RETURN
END SUBROUTINE init_grid





!-----ithion--------------------------------------------------------

SUBROUTINE ithion

!       initial thermalization of the ions (only executed if needed)
!       we assign 0.5*kb*t per degree of freedom

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

!---------------------------------------------------------------------


bk  = 1D0          ! use temperature in units of Ry
ekin=0.5D0*bk*tempion
WRITE(7,'(a,f9.4)')  'wanted temp',tempion
WRITE (7,'(a,f9.4)') 'wanted energ  per degree of freedom',ekin
xm=0.5D0*1836.0*amu(np(1))*ame
WRITE(7,*) 'warning : ithion not yet able to treat unhomogeneous systems'

!     attention to masses !

ekin=ekin/3/nion*(3*nion-6.0)
vmoy=SQRT(3.0*2.0*ekin/xm)
WRITE(7,'(a,f9.4)') 'corresponding speed',vmoy
pmoy=xm*vmoy

!     we choose random directions of the speeds

CALL spheric(pmoy,cpx,cpy,cpz,nion)
CALL conslw(cx,cy,cz,cpx,cpy,cpz,nion)

!     we rescale  the speeds

DO ion=1,nion
  ek=cpx(ion)*cpx(ion)
  ek=ek+cpy(ion)*cpy(ion)
  ek=ek+cpz(ion)*cpz(ion)
  pn=SQRT(ek)
  cpx(ion)=cpx(ion)/pn*pmoy
  cpy(ion)=cpy(ion)/pn*pmoy
  cpz(ion)=cpz(ion)/pn*pmoy
END DO
CALL conslw(cx,cy,cz,cpx,cpy,cpz,nion)
ekion=0D0
DO ion=1,nion
  ek=0D0
  ek=ek+cpx(ion)*cpx(ion)
  ek=ek+cpy(ion)*cpy(ion)
  ek=ek+cpz(ion)*cpz(ion)
  xm=0.5D0*1836.0*amu(np(1))*ame
  ek=ek/2.0/xm
  WRITE(7,'(a,i2,a,f9.4)') 'ion',ion,'ek=',ek
  ekion=ekion+ek
  WRITE(7,'(a,i2,3f9.4)') 'ion',ion,cpx(ion),cpy(ion),cpz(ion)
END DO
WRITE(7,'(a,f9.4)') 'kin.energ after renormalization',ekion
ekion=ekion/(3*nion-6.0)
WRITE(7,'(a,f9.4)') 'av.kin.energ per net degree of freedom',ekion
ekion=ekion*2.0/bk
WRITE(7,'(a,f9.4)') 'corresponding temperature',ekion



RETURN
END SUBROUTINE ithion


!-----transf ----------------------------------------------------------

SUBROUTINE transf

!     transform ions to the main axis of inertia

USE params
USE kinetic
!USE coulsolv
IMPLICIT REAL(DP) (A-H,O-Z)

REAL(DP) :: trafo(3,3),cenmas(3),tiner(3,3),dminer(3)
REAL(DP) :: pos(ng,3),help(3)
CHARACTER (LEN=3) :: ext

!-------------------------------------------------------------------------


ext=outnam

OPEN(9, FILE='center.dat', STATUS='UNKNOWN')
OPEN(10, FILE=ext//'.pdb', STATUS='UNKNOWN')

DO ion=1,nion
  pos(ion,1) = cx(ion)
  pos(ion,2) = cy(ion)
  pos(ion,3) = cz(ion)
END DO


!   compute center of mass

cenmas(1) = zero
cenmas(2) = zero
cenmas(3) = zero
gamu      = zero

DO ion = 1,nion
  DO icoo = 1,3
    cenmas(icoo) = cenmas(icoo) + amu(np(ion))*pos(ion,icoo)
  END DO
  gamu = gamu + amu(np(ion))
END DO

DO icoo = 1,3
  cenmas(icoo) = cenmas(icoo)/gamu
END DO

!  transform to center of mass coordinates

DO ion = 1,nion
  DO icoo = 1,3
    pos(ion,icoo) = pos(ion,icoo) - cenmas(icoo)
  END DO
END DO

WRITE(9,*) 'center of mass coordinates:'
WRITE(9,*) ' '

DO ion=1,nion
  cx(ion) = pos(ion,1)
  cy(ion) = pos(ion,2)
  cz(ion) = pos(ion,3)
  WRITE(9,'(3f10.4)') cx(ion),cy(ion),cz(ion)
END DO
WRITE(9,*) ' '

!  compute tensor of inertia (tiner)
!  cenmas = origin = fix point

DO i1=1,3
  DO i2=1,3
    tiner(i1,i2) = zero
  END DO
END DO

DO ion = 1,nion
  tiner(1,1) = tiner(1,1) + pos(ion,2)*pos(ion,2) + pos(ion,3)*pos(ion,3)
  tiner(2,2) = tiner(2,2) + pos(ion,3)*pos(ion,3) + pos(ion,1)*pos(ion,1)
  tiner(3,3) = tiner(3,3) + pos(ion,1)*pos(ion,1) + pos(ion,2)*pos(ion,2)
  tiner(1,2) = tiner(1,2) - pos(ion,1)*pos(ion,2)
  tiner(1,3) = tiner(1,3) - pos(ion,1)*pos(ion,3)
  tiner(2,3) = tiner(2,3) - pos(ion,3)*pos(ion,2)
END DO

tiner(2,1) = tiner(1,2)
tiner(3,1) = tiner(1,3)
tiner(3,2) = tiner(2,3)

!  diagonalize tensor of inertia

CALL jacobi(tiner,3,3,dminer,trafo,idummy)

!  transform to system of inertia
!  trafo = matrix with normalized eigenvectors in columns,
!  1/trafo = trafo (t) transforms to new coordinates

WRITE(9,'(5x,a,/,3g15.6,/)') 'Moments of inertia (x,y,z) :',  &
    dminer(1),dminer(2),dminer(3)

DO ion = 1,nion
  
  DO i=1,3
    help(i) = 0D0
    DO  j=1,3
      help(i) = help(i) + trafo(j,i)*pos(ion,j)
    END DO
  END DO
  
  DO i=1,3
    pos(ion,i) = help(i)
  END DO
  
END DO

iswap = 1

IF(iswap == 1) THEN
  
!   find and fix z-axis = axis with outstanding moment of inertia
  
  aviner = 1./3. * (dminer(1) + dminer(2) + dminer(3))
  delmx  = dminer(3) - aviner
  iax    = 3
!   check wether x or y are outstanding
  delx   = dminer(1) - aviner
  IF (ABS(delx) > ABS(delmx)) THEN
    iax    = 1
    delmx = delx
  END IF
  dely   = dminer(2) - aviner
  IF (ABS(dely) > ABS(delmx)) THEN
    iax    = 2
    delmx = dely
  END IF
!  change coordinates if necessary
  IF (iax /= 3) THEN
    DO ion = 1,nion
      zpos   = pos(ion,3)
      pos(ion,3)   = pos(ion,iax)
      pos(ion,iax) = zpos
    END DO
    WRITE(9,'(5x,a,i5,/)') 'Neue Achse : ', iax
  END IF
  IF (delmx > zero) THEN
    WRITE(9,'(5x,a,/)') ' Form: Linse'
  ELSE
    WRITE(9,'(5x,a,/)') ' Form: Zigarre'
  END IF
  
END IF    !iswap

WRITE(9,*) 'after transformation on the main axis of inertia:'
WRITE(9,*) ' '
WRITE(10,'(a10,a3)') 'COMPND    ',ext

DO ion=1,nion
  cx(ion) = pos(ion,1)
  cy(ion) = pos(ion,2)
  cz(ion) = pos(ion,3)
  WRITE(9,'(3f10.4)') cx(ion),cy(ion),cz(ion)
  WRITE(10,'(a6,a4,i1,a1,a2,a11,a1,a6,f6.3,a2, f6.3,a2,f6.3,a2,a4,a2,a4)')  &
      'HETATM',' ',ion,' ','Na',' ','1',' ', cx(ion),' ',cy(ion),' ',cz(ion),  &
      ' ','1.00',' ','0.00'
END DO

WRITE(10,'(a3)') 'END'
WRITE(9,*) ' '

DO ion1=1,nion
  k=0
  DO ion2=1,nion
    IF(ion1 /= ion2 .AND. ion1 < ion2) THEN
      dist3=SQRT((cx(ion1)-cx(ion2))**2+(cy(ion1)-cy(ion2))**2  &
          +(cz(ion1)-cz(ion2))**2)
      IF(k == 0) THEN
        WRITE(9,'(a16,i2,a4,i2,a2,f8.3)')  &
            'distance between',ion1,' and',ion2,' =',dist3
      ELSE
        WRITE(9,'(a22,i2,a2,f8.3)') '                   and',ion2,' =',dist3
      END IF
      k=k+1
    END IF
  END DO
END DO

CLOSE(9)
CLOSE(10)
RETURN
END SUBROUTINE transf

!---------------------------------------------------------------------

!     adapted from numerical recipes (p. 456):
!     compute eigenvalues (d) and eigenvectors of a symmetrical matrix
!     a(n x n)and returns the matrix v(n x n), whose columns contain
!     the normalized eigenvectors of a

!---------------------------------------------------------------------

SUBROUTINE jacobi(a,n,np,d,v,nrot)
USE params, ONLY: DP

REAL(DP), INTENT(IN OUT)                     :: a(np,np)
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN OUT)                  :: np
REAL(DP), INTENT(OUT)                        :: d(np)
REAL(DP), INTENT(OUT)                        :: v(np,np)
INTEGER, INTENT(OUT)                     :: nrot


INTEGER, PARAMETER :: nmax=500
INTEGER :: i,ip,iq,j
REAL :: c,g,h,s,sm,t,tau,theta,tresh,b(nmax),z(nmax)

DO ip=1,n
  DO iq=1,n
    v(ip,iq)=0D0
  END DO
  v(ip,ip)=1D0
END DO

DO ip=1,n
  b(ip)=a(ip,ip)
  d(ip)=b(ip)
  z(ip)=0D0
END DO

nrot=0
DO i=1,50
  sm=0D0
  DO ip=1,n-1
    DO iq=ip+1,n
      sm=sm+ABS(a(ip,iq))
    END DO
  END DO
  IF(sm == 0D0)RETURN
  IF(i < 4)THEN
    tresh=0.2D0*sm/n**2
  ELSE
    tresh=0D0
  END IF
  DO ip=1,n-1
    DO iq=ip+1,n
      g=1D2*ABS(a(ip,iq))
      IF((i > 4).AND.(ABS(d(ip))+g == ABS(d(ip)))  &
            .AND.(ABS(d(iq))+g == ABS(d(iq))))  THEN
        a(ip,iq)=0D0
      ELSE IF(ABS(a(ip,iq)) > tresh) THEN
        h=d(iq)-d(ip)
        IF(ABS(h)+g == ABS(h))THEN
          t=a(ip,iq)/h
        ELSE
          theta=0.5D0*h/a(ip,iq)
          t=1./(ABS(theta)+SQRT(1.+theta**2))
          IF(theta < 0D0)t=-t
        END IF
        c=1./SQRT(1+t**2)
        s=t*c
        tau=s/(1.+c)
        h=t*a(ip,iq)
        z(ip)=z(ip)-h
        z(iq)=z(iq)+h
        d(ip)=d(ip)-h
        d(iq)=d(iq)+h
        a(ip,iq)=0D0
        DO j=1,ip-1
          g=a(j,ip)
          h=a(j,iq)
          a(j,ip)=g-s*(h+g*tau)
          a(j,iq)=h+s*(g-h*tau)
        END DO
        DO j=ip+1,iq-1
          g=a(ip,j)
          h=a(j,iq)
          a(ip,j)=g-s*(h+g*tau)
          a(j,iq)=h+s*(g-h*tau)
        END DO
        DO j=iq+1,n
          g=a(ip,j)
          h=a(iq,j)
          a(ip,j)=g-s*(h+g*tau)
          a(iq,j)=h+s*(g-h*tau)
        END DO
        DO j=1,n
          g=v(j,ip)
          h=v(j,iq)
          v(j,ip)=g-s*(h+g*tau)
          v(j,iq)=h+s*(g-h*tau)
        END DO
        nrot=nrot+1
      END IF
    END DO
  END DO
  DO ip=1,n
    b(ip)=b(ip)+z(ip)
    d(ip)=b(ip)
    z(ip)=0D0
  END DO
END DO

STOP 'too many iterations in jacobi'

RETURN
END SUBROUTINE jacobi

