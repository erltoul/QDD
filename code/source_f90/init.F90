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

#include"define.h"
 
!------------------------------------------------------------

SUBROUTINE initnamelists

!     Sets defaults for input variables
!     and then reads input variables through namelist.

!------------------------------------------------------------
USE params
USE kinetic

IMPLICIT NONE

CHARACTER (LEN=80) :: title
!      character*10 fname(0:3)
!      data fname/'for005.011',
!     &           'for005.100','for005.101','for005.110'/

INTEGER,PARAMETER :: kparall=11
INTEGER:: iu, nnx2
REAL(DP)::dx2
#if(simpara)
CHARACTER (LEN=10) :: fname(0:kparall)
!      data fname/'for005.001','for005.010','for005.011',
!     &           'for005.100','for005.101','for005.110'/
DATA fname/'for005.001','for005.010','for005.011',  &
    'for005.100','for005.101','for005.110',  &
    'forjel.001','forjel.010','forjel.011',  &
    'forjel.100','forjel.101','forjel.110'/
#endif
#if(raregas)
NAMELIST /global/   nclust,nion,nspdw,nion2,nc,nk,numspin,  &
    temp,occmix,isurf,b2occ,gamocc,deocc,osfac,  &
    init_lcao,kstate,kxbox,kybox,kzbox,dx,dy,dz,  &
    radjel,surjel,bbeta,gamma,beta4,endcon,itback,  &
    epswf,e0dmp,epsoro,  dpolx,dpoly,dpolz,  &
    scaleclust,scaleclustx,scaleclusty,scaleclustz, &
    shiftclustx,shiftclusty,shiftclustz,  &
    rotclustx,rotclusty,rotclustz,imob,iswitch_interpol,  &
    idebug,ishiftcmtoorigin,iaddcluster,  &
    iswforce,iplotorbitals,ievaluate,ehom0,ihome,  &
    ehomx,ehomy,ehomz,shiftwfx,shiftwfy,shiftwfz, ispinsep, &
    nproj_states
#else
NAMELIST /global/   nclust,nion,nspdw,nion2,numspin,  &
    temp,occmix,isurf,b2occ,gamocc,deocc,osfac,  &
    init_lcao,kstate,kxbox,kybox,kzbox,dx,dy,dz,  &
    radjel,surjel,bbeta,gamma,beta4,endcon,itback,  &
    epswf,e0dmp,epsoro,  dpolx,dpoly,dpolz,  &
    scaleclust,scaleclustx,scaleclusty,scaleclustz, &
    shiftclustx,shiftclusty,shiftclustz,  &
    rotclustx,rotclusty,rotclustz,imob,iswitch_interpol,  &
    idebug,ishiftcmtoorigin,iaddcluster,  &
    iswforce,iplotorbitals,ievaluate,ehom0,ihome,  &
    ehomx,ehomy,ehomz,shiftwfx,shiftwfy,shiftwfz, ispinsep, &
    nproj_states
#endif


!*************************************************************

!k use of parameters for dynamics:

!  use of discrete ions via pseudopotentials:      nion2 not=0
!  use of homogeneous background via jellium:       nion2=0

!  iforce=0
!  polarized clusters (e.g. polar. isomer of na_12)   iforce=1

!  dipoleboost or spindipoleboost      iexcit=0:
!    boost of electronic density by centfx/y/z       ispidi=0
!    boost of spinup density by +0.5*centfx/y/z     ispidi=1
!             spindown density by -0.5*centfx/y/z   ispidi=1
!  dipoleshift (iexcit=0) by shiftinix/y/z


!  scissor mode excitation             iexcit=1:
!    chose axis of rotation by 'irotat'
!    chose angle of rotation by 'phirot'

!*************************************************************

NAMELIST /dynamic/ directenergy,nabsorb,idenfunc,  &
    iemomsrel,ifsicp,ionmdtyp,ifredmas,modionstep,icooltyp,ipsptyp,  &
    ipseudo,ismax,itmax,isitmax,isave,istinf,ipasinf,dt1,irest,  &
    centfx,centfy,centfz, shiftinix,shiftiniy,shiftiniz, &
    ispidi,iforce,iexcit,iangmo,  &
    irotat,phirot,i3dz,i3dx,i3dstate,istream,iflocaliz,  &
    idyniter,ifrhoint_time,ifhamdiag,iffastpropag, &
    modrho,jpos,jvel,jener,jesc,jforce,istat,jgeomion,  &
    jdip,jdiporb,jquad,jang,jangabso,jspdp,jinfo,jenergy,ivdw,  &
    jposcm,mxforce,myforce,mzforce,jgeomel,jelf,jstinf, &
    jstboostinv,ifspemoms,iftransme,ifexpevol,ifcnevol, &
    tempion,idenspl,ekmat,nfix,  &
    itft,tnode,deltat,tpeak,omega,e0,  &
    projcharge,projvelx,projvely,projvelz, &
    projinix,projiniy,projiniz, &
    e1x,e1y,e1z,e2x,e2y,e2z,phi,  &
    phase2,omega2,e0_2,tstart2,tpeak2, &
    izforcecorr, dinmargin,  &
    ntref,iangabso,ipes,nangtheta,nangphi,  &
    delomega,angthetal,angthetah,angphil,angphih,  &
    ifreezekspot,powabso,ispherabso,ifixcmion,  &
    ekin0pp,vxn0,vyn0,vzn0,jescmask,itof,jescmaskorb,  &
    eproj,nproj,nproj_states,vpx,vpy,vpz,taccel,   &
    trequest,timefrac,  &
    nmptheta,nmpphi,jmp,jovlp,  &
    jnorms,jplotdensitydiff,jplotdensitydiff2d,  &
    jplotdensity2d,jcharges,drcharges, &
    iscatterelectron,jattach,scatterelectronenergy,  &
    scatterelectronvxn,scatterelectronvyn, &
    scatterelectronvzn,scatterelectronx,  &
    scatterelectrony,scatterelectronz,scatterelectronw,jattach, &
    phangle,phphase,nhstate,npstate, &
    jstateoverlap


NAMELIST /surface/  &
#if(raregas)
    surftemp,ipotfixed,ifmdshort,ifadiadip,  &
    sigmac,sigmav,sigmak,isystem,jsavesurf, chgc0,chge0,chgk0, &
    cspr,mion,me,mkat,isrtyp,isrtypall,  &
    shiftx,shifty,shiftz,scaledist, &
    iprintonlyifmob, iaxis,  &
    fermia,fermib,fermic, fermiac,fermibc,fermicc,  &
    fermiae,fermibe,fermice, fermiak,fermibk,fermick,  &
    fermia2c,fermib2c,fermic2c, fermia2e,fermib2e,fermic2e,  &
    fermia2k,fermib2k,fermic2k, sigkv,bkv,ckv,  &
    sigvv,bvv,cvv, sigkk,bkk,ckk,  &
    sigcv,bcv,ccv, sigck,bck,cck,  &
    sigcc,bcc,ccc, sigcn,bcn,ccn,  &
    sigkn,bkn,ckn, sigkn2,bkn2,  &
    ccn2,ccn3,ccn4,ccn5,ccn6,ccn7,ccn8,ccnd,dcn,  &
    ccn9,ccn10,ccncu2,ccncu3,ccncu4,ccncu5,ccncu6,ccncu7,  &
    ccncu8,ccncu9,ccncu10,ccncud,  &
    ccncul2,ccncul3,ccncul4,ccncul5,ccncul6,ccncul7,  &
    ccncul8,ccncul9,ccncul10,ccnculd,  &
    ckn2,ckn3,ckn4,ckn5,ckn6,ckn7,ckn8,cknd,dkn,  &
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
    iusecell,iuselast,ibh,cbh,iforcecl2co, &
    enerinfty,epsdi,idielec,xdielec,iararlj
#else
    epsdi,idielec,xdielec
#endif


WRITE(6,*) 'Reading for005.// ...'


!     initialize the variables with default values



delomega=(angthetah-angthetal)/4D0/nangtheta
delomega=MIN(delomega,(angphih-angphil)/4D0/nangphi)

scatterelectronz=nzsh*dz-4D0*scatterelectronw

#if(raregas)
CALL init_raregas()
#endif


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
    WRITE(6,*) ' title is now: '//outnam
    WRITE(iu,*) ' title is now: '//outnam
    CLOSE(5)
write(6,*) outnam, "**************************************"
    
    OPEN(5,STATUS='old',FORM='formatted',FILE='for005.'//outnam)
    
!     read input from namelists

    READ(5,global,END=99999)
    WRITE(*,*) ' GLOBAL read'
    IF(dx*dy*dz == 0D0) STOP ' DX & DY & DZ must be given explicitely'
    IF(dx<0D0) THEN
      WRITE(6,*) 'negative DX'
      WRITE(7,*) 'negative DX'
      OPEN(119, FILE='dx', STATUS='UNKNOWN')
      READ (119,*) dx2
      CLOSE(119)
      OPEN(120, FILE='nx', STATUS='UNKNOWN')
      READ(120,*) nnx2
      CLOSE(120)
      IF(dx2<0D0) THEN
        WRITE(6,*) 'negative DX in file DX'
        WRITE(6,*) 'DX,NX are computed then program restart'
        WRITE(7,*) 'negative DX in file DX'
        WRITE(7,*) 'DX,NX are computed then program restart'
      ELSE
        dx=dx2
        dy=dx2
        dz=dx2
        kxbox=nnx2
        kybox=nnx2
        kzbox=nnx2
        WRITE(6,*) 'DXDYDZ changed to',DX,'from file DX'
        WRITE(7,*) 'DXDYDZ changed to',DX,'from file DX'
        WRITE(6,*) 'KXBOX KYBOX KZBOX changed to',nnx2,'from file NX'
        WRITE(7,*) 'KXBOX KYBOX KZBOX changed to',nnx2,'from file NX'
      ENDIF
    ENDIF
            
    IF(nproj_states>0)THEN
      ALLOCATE(proj_states(nproj_states))
      proj_states(:)=0
    END IF
    READ(5,dynamic,END=99999)
    WRITE(*,*) ' DYNAMIC read'
#if(raregas)    
    READ(5,surface,END=99999)
    WRITE(*,*) ' SURFACE read'
#endif
    99999    CLOSE(5)

    IF(jdip<0) STOP "you must specify JDIP in namelist DYNAMIC"
    IF(jesc<0) STOP "you must specify JESC in namelist DYNAMIC"
    IF(nclust<=0) THEN ! total charge of the cluster is read
      WRITE(*,*) 'zero or negative nclust : global charge of',nclust
    ENDIF


#if(simpara)
  ELSE IF(myn > 5) THEN
    STOP ' this node not  active'
#endif
  END IF
  
  IF (ionmdtyp==0) THEN
     jpos=0
     jvel=0
  ELSE
    IF(jpos<0) STOP "you must specify JPOS in namelist DYNAMIC"
    IF(jvel<0) STOP "you must specify JVEL in namelist DYNAMIC"
  END IF
! adapt input parameters if necessary
  IF(nion2 == 0) iemomsrel=0    ! relat. to center of box for jellium

  
#if(parayes|paraworld)
  CALL comm_inputparams()
#endif
  
  
  tdipolxyz = dpolx*dpolx+dpoly*dpoly+dpolz*dpolz .GT. 0D0
  
!      write(*,*) ' INIT: nion2=',nion2
  RETURN
  
END SUBROUTINE initnamelists



!-----changePerio--------------------------------------------

SUBROUTINE changeperio

!     Reads pseudo potential parameters from namelist

USE params
IMPLICIT NONE

INTEGER :: info
!      namelist /perio/ ch(-99:99),amu(-99:99),
NAMELIST /perio/ ch,amu, cc1,cc2,crloc,  &
    h0_11g,h0_12g,h0_22g,h0_33g, h1_11g,h1_22g,h2_11g,  &
    dr1,dr2,prho1,prho2, r0g,r1g,r2g,radiong

OPEN(5,STATUS='old',FORM='formatted',FILE='for005.'//outnam)
READ(5,perio, IOSTAT=info)
IF(IS_IOSTAT_END(info)) WRITE(6,*) 'End of file reached while reading for perio'
CLOSE(5)
!test         write(iu,*) ' changeperio done. myn=',myn

#if(parayes)
CALL comm_periodic()
#endif

END SUBROUTINE changeperio


!-----iparams--------------------------------------------------

SUBROUTINE iparams()

!     check consistency of input parameters, do some initializations

USE params
IMPLICIT NONE

INTERFACE 
  SUBROUTINE calc_lda_gunnar(rho,chpdft)
  USE params, ONLY: DP,kdfull2
  REAL(DP), INTENT(IN)                         :: rho(2*kdfull2)
  REAL(DP), INTENT(OUT)                        :: chpdft(2*kdfull2)
  END SUBROUTINE calc_lda_gunnar
END INTERFACE
INTERFACE 
  SUBROUTINE calc_lda_pw92(rho,chpdft)
  USE params, ONLY: DP,kdfull2
  REAL(DP), INTENT(IN)                         :: rho(2*kdfull2)
  REAL(DP), INTENT(OUT)                        :: chpdft(2*kdfull2)
  END SUBROUTINE calc_lda_pw92
END INTERFACE

#if(raregas)
INTEGER :: i,j
#endif
!--------------------------------------------------------------

!     some initializations
#if(raregas)
NE = nc
#endif
trequest=trequest*60D0
phi=phi*pi/180D0              ! convert input 'phi' from degree




WRITE(6,*) 'i3dz,i3dx,i3dstate,idenspl=' ,i3dz,i3dx,i3dstate,idenspl
qold2=0.01D0
qold3=0.01D0
qold4=0.01D0
jekion=0
iquery4=0
WRITE(6,*) 'jekion=',jekion


!     check consistency of input options

IF(numspin.NE.2 .AND. iftransme ==1) STOP ' IFTRANSME needs full spin'


#if(raregas)
DO i=1,5
  DO j=1,5
    IF (isrtyp(i,j) == 2 .AND. iswforce == 1)  &
        STOP 'New force routine not implemented for Ar. Use iswForce=0!'
  END DO
END DO
#endif

IF(nion > ng) STOP " more ions than dimensioned. enhance NG"

IF(myn == 0)THEN
  IF(irest > itmax) STOP ' IREST > ITMAX is nonsense'
END IF

IF(nion2 == 0 .AND. ifsicp >= 3)  &
    STOP 'Jellium not compatible with KLI or Slater-SIC'
IF(nion2 == 2 .AND. ifsicp >= 3)  &
    STOP 'external PsP not compatible with KLI or Slater-SIC'

#if(raregas)
IF(nc+NE+nk.gt.ngpar) STOP ' not enough NGPAR for substrate'
#endif

!#if(hamdiag&parayes)
!IF(MOD(kstate,2) == 1) STOP ' KSTATE must be even'
!#endif

#if(parayes)
IF(ifhamdiag == 1) STOP ' step with H diagonalization only on serial code'
!IF(jmp /=0) STOP ' evaluation of PES not possible in parallele code'     cPW
!IF(istat /= 0) STOP ' static restart not possible in parallele code'
#endif

#if(findiff|numerov)
!IF(ifhamdiag>0) STOP ' step with H diagonalization not yet for fin.diff.'
#endif


IF(directenergy .AND.  &
   ifsicp /= 3 .AND. ifsicp /= 4 .AND. ifsicp /= 0 .AND. ifsicp /= 8)  &
    STOP ' directenergy=1 only for Slater, KLI, or 2st-SIC'

IF(numspin.NE.2 .AND. ifsicp >= 3) STOP 'IFSICP>2 requires fullspin code'

#if(cmplxsic)
#if(!twostsic) 
  STOP 'CMPLXSIC requires TWOSTSIC=1'
#endif
#endif

#if(twostsic)
#if(locsic)
STOP ' TWOSTSIC and LOCSIC cannot run simultaneously'
#endif
#if(parayes)
STOP ' TWOSTSIC cannot yet run in MPI parallel code'
#endif
IF(ifsicp==8 .AND. .NOT.directenergy) &
   STOP 'full SIC (IFSICP=8) requires DIRECTENERGY=.TRUE.'
#endif

IF(ifsicp.NE.8 .AND. isitmax.NE.0) STOP 'ISITMAX.NE.0 only for SIC'

IF(isitmax>0 .AND. ifexpevol== 0) &
    STOP ' imaginary-time step only for exponential evolution'

IF(ifexpevol == 1 .AND. ionmdtyp /= 0)  &
    STOP ' exponential evolution not with ionic motion'    !  why?

if(nabsorb == 0 .AND. jesc .NE. 0) &
    STOP ' JESC must be zero for NABSORB=0'

!#if(selpara)
!      if(idielec.ne.0 .or. isurf.ne.0)
!     &  stop 'option SELPARA not compatible with substrates'
!#endif

#if(parayes)
!  IF(jstinf /= 0) &
!   STOP ' print of s.p.energies & variances not for parallel version'
  IF(jgeomel /= 0) &
   STOP ' print electronic geometry not for parallel version'
#endif

IF(nion2 == 2) THEN
  IF(ipseudo /= 0) STOP ' IPSEUDO=0 needed for external psp (NION2=2)'
  IF(iexcit /= 0) STOP ' IEXCIT=0 needed for external psp (NION2=2)'
  IF(ipsptyp /= 0) STOP ' IPSPTYP=0 needed for external psp (NION2=2)'
  IF(imob /= 0) STOP ' fixed ions needed for external psp (IMOB=0)'
END if

IF(ABS(phangle) > small .AND. istat /= 1) &
  STOP ' ph rotation only with ISTAT=1'

! set the pointer for the energy-density functional
IF(idenfunc==1) THEN
  calc_lda => calc_lda_pw92
ELSE IF(idenfunc==2 .OR. idenfunc==3) THEN
  calc_lda => calc_lda_gunnar
ELSE
  STOP ' invalid value for density-functional selector IDENFUNC'
END IF


RETURN
END SUBROUTINE iparams


!-----ocoption--------------------------------------------------

SUBROUTINE ocoption(iu)

!     output of compiled and of read-in options

USE params
IMPLICIT NONE
INTEGER :: iu
!----------------------------------------------------------------

IF(iu == 8) OPEN(8,STATUS='unknown',FORM='formatted',  &
    FILE='poptions.'//outnam)

WRITE(iu,'(a)') 'the following options are used:'
WRITE(iu,*)
IF(idenfunc==1) THEN
  WRITE(iu,'(a)') 'Perdew-Wang 92 exchange-correlation functional'
ELSE IF(idenfunc==2) THEN
  WRITE(iu,'(a)') 'gl 76 exchange-correlation functional'
ELSE IF(idenfunc==3) THEN
  WRITE(iu,'(a)') 'pure exchange functional in LDA'
ELSE
  STOP ' invalid value for density-functional selector IDENFUNC'
END IF

! parallel or serial:
#if(parayes)
WRITE(iu,*) 'parallel code: number of nodes=',knode
WRITE(iu,*) 'max. number of wf per node',kstate
#endif
#if(parano)
WRITE(iu,*) 'serial code: max. number of wf=',kstate
#endif
IF(numspin==2) THEN
  WRITE(iu,'(a)') ' spin code '
ELSE IF(numspin==1) THEN
  WRITE(iu,'(a)') ' no-spin code '
ELSE
  STOP "invalid uinput parameter NUMSPIN"
END IF
IF(ifhamdiag>0) THEN
  WRITE(iu,'(a)') ' static step: Hamiltonian in subspace diagonalized'
ELSE
  WRITE(iu,'(a)') ' static step: Hamiltonian in subspace not diagonalized'
END IF
! fft:
#if(gridfft)
WRITE(iu,'(a)') 'Fourier propagation'
#if(netlib_fft)
WRITE(iu,'(a)') 'using netlib ffts'
#endif
#if(fftw_cpu)
WRITE(iu,'(a)') 'using fftw@cpu'
#endif
#endif
#if(coufou)
WRITE(iu,'(a)') 'falr coulomb solver'
#endif
#if(coudoub)
WRITE(iu,'(a)') 'exact coulomb solver'
#endif
IF(iffastpropag == 1) WRITE(iu,'(a)') 'faster TV propagator'

IF(ifexpevol == 1) THEN
  WRITE(iu,'(a)') 'time step by exponential evolution'
ELSE
  WRITE(iu,'(a)') 'time step by T-V splitting'
END IF
IF(tdipolxyz) WRITE(iu,'(a)') 'hom. electric field switched on'

#if(raregas)
IF(ifadiadip == 1) THEN
  WRITE(iu,'(a)') 'substrate: adiabatic dipoles'
ELSE
  WRITE(iu,'(a)') 'substrate: dynamic dipoles'
END IF
IF(ivdw ==1) WRITE(iu,'(a)') 'substrate: VdW activated'
WRITE(iu,'(a,i2)') 'contribution short-range V to energy: yes/no=',ifmdshort
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
  WRITE(iu,'(a)') 'full Goedecker pseudopotentials'
ELSE IF(ipsptyp == 2) THEN
  WRITE(iu,'(a)') 'local Goedecker pseudopotentials'
ELSE IF(ipsptyp == 3) THEN
  WRITE(iu,'(a)') 'full Goedecker pseudopotentials form file goed.asci'
ELSE IF(ipsptyp == 4) THEN
  WRITE(iu,'(a)') 'full semicore Goedecker pseudopotentials form file goed.asci'
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
  WRITE(iu,'(a)') 'molecular dynamics: velocity Verlet'
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

IF(directenergy) THEN
  WRITE(iu,'(a)') ' energy computed directly'
ELSE
  WRITE(iu,'(a)') ' energy computed using s.p. energies'
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
ELSE IF(ifsicp == 5) THEN
  WRITE(iu,'(a)') 'sic activated: exact exchange'
IF(ifsicp == 6)  STOP ' IFSICP=6 presently not provided'
#if(twostsic)
ELSE IF(ifsicp == 7)  THEN
  WRITE(iu,'(a)') 'sic activated: localized SIC'
#if(parayes)
  STOP " IFSICP=7 not possible in parallel code"
#endif
ELSE IF(ifsicp == 8)  THEN
  WRITE(iu,'(a)') 'sic activated: double-set SIC'
#else
ELSE IF(ifsicp == 7)  THEN
  STOP ' code not compiled for localized SIC'
ELSE IF(ifsicp == 8)  THEN
  STOP ' code not compiled for double-set SIC'
#endif
ELSE
  WRITE(iu,'(a)') 'this version of SIC not available'
END IF

!IF(.NOT.(ifexpevol==1) .AND. itmax > 0 .AND. ifsicp >= 6) STOP  &
!    ' full TDSIC requires exponential evolution'


!     dynamical options


WRITE(iu,'(a,3i6)') ' kxbox,kybox,kzbox=',kxbox,kybox,kzbox
IF(nclust<=0) then
#if(raregas)
  WRITE(iu,'(a,4i6)') ' charge,nion,nrare,nstate=',nclust,nion,nrare,nstate
#else
  WRITE(iu,'(a,3i6)') ' charge,nion,nstate=',nclust,nion,kstate
#endif
ELSE
#if(raregas)
  WRITE(iu,'(a,4i6)') ' nelect,nion,nrare,nstate=',nclust,nion,nrare,nstate
#else
  WRITE(iu,'(a,3i6)') ' nelect,nion,nstate=',nclust,nion,kstate
#endif
ENDIF
WRITE(iu,'(a,4i3,f7.2)') ' ispidi,iforce,iexcit,irotat,phirot=',  &
    ispidi,iforce,iexcit,irotat,phirot
WRITE(iu,'(a,3f8.2)') ' boost: centfx,centfy,centfz=',centfx,centfy,centfz
WRITE(iu,'(a,3f8.2)') ' shift: shiftinix,y,z=',shiftinix,shiftiniy,shiftiniz
WRITE(iu,'(a,2i3,f7.3)') ' irest,istat,dt1=',irest,istat,dt1
WRITE(iu,'(a/a,i3,5f8.3/a,7f8.3)') ' laser:',  &
    ' itft,tnode,deltat,tpeak,omega,e0=', itft,tnode,deltat,tpeak,omega,e0,  &
    ' e1x,e1y,e1z,e2x,e2y,e2z,phi=', e1x,e1y,e1z,e2x,e2y,e2z,phi

#if(raregas)
IF (isurf /= 0) THEN
  WRITE(iu,*) '*************************************************'
  WRITE(iu,*) 'SURFACE/MATRIX PRESENT:'
  WRITE(iu,'(a,i4,a,i4)') 'GSM particles: ',nc,'  cations: ',nk
  WRITE(iu,*) '*************************************************'
END IF
#endif


WRITE(iu,*) 'CODE VERSION: ', IVERSION

IF(iu == 8) CLOSE(8)

RETURN
END SUBROUTINE ocoption

!-----init_output--------------------------------------------------

SUBROUTINE init_output()

!     write headers of output files

USE params
IMPLICIT NONE
INTEGER :: i,j,maxnum
REAL(DP) :: ascal
CHARACTER (LEN=3) :: num

!------------------------------------------------------------------

IF(myn < 10) THEN
  WRITE(num,'(i1)') myn
  maxnum = 1
ELSE IF(myn < 100 .AND. myn > 9) THEN
  WRITE(num,'(i2)') myn
  maxnum=2
ELSE
  WRITE(num,'(i3)') myn
  maxnum=3
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

WRITE(6,'(a,3i5,i8,g12.4)') ' INIT: nx,ny,nz,nxyz,dvol=',  &
    nx,ny,nz,nxyz,dvol




!       output static parameters (traditional sodium case)

IF(nclust>0) THEN
  WRITE(7,'(5x,a,i2,a,i2/a/)')  &
  '# electr.: ',nclust,'# ions: ',nion,'=========='
ELSE
  WRITE(7,'(5x,a,i2,a,i2/a/)')  &
  '# charge.: ',nclust,'# ions: ',nion,'=========='
       nstate=nion+nclust
  IF(ipsptyp<1  .AND. ipseudo<1) THEN
    nclust=nstate
    nspdw=nstate/2
  ENDIF
ENDIF
WRITE(7,'(a,i3)') 'number of alkali wave-fkts nstate = ',nstate

WRITE(7,'(a,f4.2)') 'initialisation : osfac=',osfac
WRITE(7,'(a,10(/t2,3(3i4,5x),3i4))')  &
    'quantumnumbers :  ',((nq(i,j),i=1,3),j=1,nstate)
IF(temp > 0D0) THEN
  WRITE(7,'(a,f10.5,a,f4.2)') 'temperature kt= ',temp,'Ry,   occmix=',occmix
  WRITE(7,'(a,1g13.5)') 'epsoro=',epsoro
  WRITE(6,'(a,f10.5,a,f4.2)') 'temperature kt= ',temp,'Ry,   occmix=',occmix
  WRITE(6,'(a,1g13.5)') 'epsoro=',epsoro
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
      WRITE(7,*) 0.5D0*centfx,0.5D0*centfy,0.5D0*centfz
    ELSE IF(ispidi == -1) THEN
      WRITE(7,*) 'shift of spindown density by'
      WRITE(7,*) -0.5D0*centfx,-0.5D0*centfy,-0.5D0*centfz
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
    WRITE(7,*) 'Gaussian laser pulse'
  END IF
  WRITE(7,*) 'length of the pulse',deltat
  WRITE(7,*) 'peak time',tpeak
  WRITE(7,*) 'begins at',tnode
  WRITE(7,*) 'pulsation',omega,'units'
  WRITE(7,*) 'field strength',e0
  IF(phi == 0D0) WRITE(7,*) 'linear polarization'
  IF(phi == 0.25D0) WRITE(7,*) 'circular right polarization'
  IF(phi == -0.25D0) WRITE(7,*) 'circular left polarization'
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
IMPLICIT NONE
!      dimension dr1(-ng:ng),dr2(-ng:ng) ! bugBF
!      dimension prho1(-ng:ng),prho2(-ng:ng) ! bugBF
INTEGER :: i,iel,natom,icountt,inew,iskip
INTEGER :: l,nactual,nval
REAL(DP) :: amfac,c1,c2,c3,c4,crr,h11,h22,h33,rloc
REAL(DP) :: scpx,scpy,scpz
character (len=3) ::  naml
character (len=2)  ::  namc,symb(99)
character(len=80) ::  a

!----------------------------------------------------------------


WRITE(6,*) 'Entering iperio'


DO iel=-99,99
  amu(iel)= 0D0   ! masses as a function of atomic number
END DO

!     soft PsP (error functions)

IF(ipsptyp == 0) THEN
  
!       hydrogen
  
  amu(1)=1D0
  ch(1)=1D0
  dr1(1)=0.3D0
  dr2(1)=0.45D0
  prho1(1)=6.94527D0
  prho2(1)=-0.92056D0
  sgm1(1)=dr1(1)*0.8493218D0
  sgm2(1)=dr2(1)*0.8493218D0
  chg1(1)=((sgm1(1)*SQRT(pi*2D0))**3)*prho1(1)*ch(1)
  chg2(1)=((sgm2(1)*SQRT(pi*2D0))**3)*prho2(1)*ch(1)
  WRITE(6,*) 'ch(1)=',ch(1)
  WRITE(6,*) 'total charge of Psphydr=',chg1(1)+chg2(1)
  
!       sodium
  
  amu(11)=23.0D0
  ch(11)=1D0  ! charge of pseudopotential
  dr1(11)=0.8018D0
  dr2(11)=1.3693D0
  prho1(11)=-0.46073D0
  prho2(11)=0.13287D0
  sgm1(11)=dr1(11)*0.8493218D0
  sgm2(11)=dr2(11)*0.8493218D0
  chg1(11)=((sgm1(11)*SQRT(pi*2D0))**3)*prho1(11)*ch(11)
  chg2(11)=((sgm2(11)*SQRT(pi*2D0))**3)*prho2(11)*ch(11)
  WRITE(6,*) 'ch(11)=',ch(11)
  WRITE(6,*) 'total charge of Pspsodi=',chg1(11)+chg2(11)
  
!       magnesium
  
  amu(12)=24.312D0
  ch(12)=2D0
  dr1(12)=0.5D0
  dr2(12)=1D0
  prho1(12)=-1.6582D0
  prho2(12)=0.31091D0
  sgm1(12)=dr1(12)*0.8493218D0
  sgm2(12)=dr2(12)*0.8493218D0
  chg1(12)=((sgm1(12)*SQRT(pi*2D0))**3)*prho1(12)*ch(12)
  chg2(12)=((sgm2(12)*SQRT(pi*2D0))**3)*prho2(12)*ch(12)
  WRITE(6,*) 'ch(12)=',ch(12)
  WRITE(6,*) 'total charge of Pspmagn=',chg1(12)+chg2(12)
  
!       Argon , positive core part
  
  amu(18)   = 39.95D0
  ch(18)    = 6.119D0  ! charge of pseudopotential
  sgm1(18)  = 1.43D0/sq2
  dr1(18)   = sgm1(18)/0.8493218D0
  sgm2(18)  = 1D0
  dr2(18)   = 1D0
  chg1(18)  = ch(18)
  prho1(18) = chg1(18)/((sgm1(18)*SQRT(pi*2D0))**3)
  chg2(18)  = 0D0
  prho2(18) = 0D0
  WRITE(6,*) 'ch(18)=',ch(18)
  WRITE(6,*) 'total charge of PspAr_core=',chg1(18)+chg2(18)
  
!       Argon , negative valence cloud (inert)
  
  amu(-18)  = 4.38D0/1836.0D0
  ch(-18)   = -6.119D0  ! charge of pseudopotential
  sgm1(-18) = 1.43D0/sq2
  dr1(-18)  = sgm1(-18)/0.8493218D0
  sgm2(-18) = 1D0
  dr2(-18)  = 1D0
  chg1(-18) = ch(-18)
  prho1(-18)= chg1(-18)/((sgm1(-18)*SQRT(pi*2D0))**3)
  chg2(-18) = 0D0
  prho2(-18)= 0D0
!  c_dipmod  = 3.378  !  ch(-18)**2*e2/11.08
  WRITE(6,*) 'ch(-18)=',ch(-18)
  WRITE(6,*) 'total charge of PspAr_el=',chg1(-18)+chg2(-18)
  
!       potassium
  
  amu(19)=39.0D0
  ch(19)=1D0  ! charge of pseudopotential
  dr1(19)=0.9973D0
  dr2(19)=1.9957D0
  prho1(19)=-0.13771D0
  prho2(19)=0.030223D0
  sgm1(19)=dr1(19)*0.8493218D0
  sgm2(19)=dr2(19)*0.8493218D0
  chg1(19)=((sgm1(19)*SQRT(pi*2D0))**3)*prho1(19)*ch(19)
  chg2(19)=((sgm2(19)*SQRT(pi*2D0))**3)*prho2(19)*ch(19)
  WRITE(6,*) 'ch(19)=',ch(19)
  WRITE(6,*) 'total charge of Psppotassium=',chg1(19)+chg2(19)
  
!       cesium
  
  amu(58)=140D0
  ch(58)=1D0  ! charge of pseudopotential
  dr1(58)=1.12D0
  dr2(58)= 2.552D0
  prho1(58)=-0.62923D-01
  prho2(58)=0.11554D-01
  sgm1(58)=dr1(58)*0.8493218D0
  sgm2(58)=dr2(58)*0.8493218D0
  chg1(58)=((sgm1(58)*SQRT(pi*2D0))**3)*prho1(58)*ch(58)
  chg2(58)=((sgm2(58)*SQRT(pi*2D0))**3)*prho2(58)*ch(58)
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
!  nrow(1)  = 1
  
!       helium
  
  amu(2)   = 4.0D0
  ch(2)    = 2.0D0
  cc1(2)   =-9.112023D0
  cc2(2)   = 1.698368D0
  crloc(2) = 0.2D0
!  nrow(2)  = 1

  
!       bor
  
  amu(5)   = 10.81D0
  ch(5)    = 3.0D0
  cc1(5)   =-5.578642D0
  cc2(5)   = 0.804251D0
  crloc(5) = 0.43393D0
!  nrow(5)  = 2
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
!  nrow(6)  = 2
  r0g(6)=0.304533D0
  r1g(6)=0.232677D0
  radiong(6)=2D0
  h0_11g(6)=9.522842D0
  h1_11g(6)=0D0
  
!       nitrogen
  
  amu(7)   = 14.0D0
  ch(7)    = 5.0D0
  cc1(7)   =-12.2046419D0
  cc2(7)   = 1.7558249D0
  crloc(7) = 0.2889046D0
!  nrow(7)  = 2
  r0g(7)=0.256605D0
  r1g(7)=0.270134D0
  radiong(7)=2.0D0
  h0_11g(7)=13.552433D0
  h1_11g(7)=0D0
  
!       oxygen
  
  amu(8)   = 16.0D0
  ch(8)    = 6.0D0
  cc1(8)   =-16.4822284D0    ! -16.580318
  cc2(8)   = 2.3701353D0     ! 2.395701
  crloc(8) = 0.2477535D0     ! 0.247621
!  nrow(8)  = 2
  r0g(8)=0.2222028D0      !  0.221786D0
  r1g(8)=0.256829D0
  radiong(8)=1.5D0
  h0_11g(8)=18.19996387D0   ! 18.266917D0
  h1_11g(8)=0D0
  
!       fluor   (row 2)
  
  amu(9)  = 18D0
  ch(9)   = 7.0D0
 cc1(9)   =-21.307361D0
 cc2(9)   = 3.072869D0
 crloc(9) = 0.218525D0
 r0g(9)=0.195567D0
 r1g(9)=0.2D0
 radiong(9)=1.3D0
 h0_11g(9)=23.58494D0
 h1_11g(9) =0D0
  
!       neon   (row 2)
  
  amu(10)  = 20.2D0
  ch(10)   = 8.0D0
  cc1(10)  = -27.692852D0
  cc2(10)  = 4.005906D0
  crloc(10)= 0.19D0
!  nrow(10) = 2
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
!  nrow(11) = 3
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
!  nrow(12) = 3
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
!  nrow(13) = 3
  r0g(13)=0.460104D0
  r1g(13)=0.536744D0
  h0_11g(13)=5.088340D0
  h0_22g(13)=2.679700D0
  h1_11g(13)=2.193438D0
  
!       silicon
  
  amu(14)  = 28.09D0
  ch(14)   = 4.0D0
  cc1(14)  =-7.336103D0
  cc2(14)  = 0D0
  crloc(14)= 0.44D0
!  nrow(14) = 3
  
!       phosphor
  
  amu(15)  = 30.974D0
  ch(15)   = 5.0D0
  cc1(15)  =-6.65422D0
  cc2(15)  = 0D0
  crloc(15)= 0.43D0
!  nrow(15) = 3
  r0g(15)=0.389803D0
  r1g(15)=0.440796D0
  h0_11g(15)=6.842136D0
  h0_22g(15)=3.856693D0
  h1_11g(15)=3.282606D0
  
!       argon   (row 3)
  
  amu(18)  = 40D0
  ch(18)   = 8.0D0
  cc1(18)  = -7.1D0 !nbf9.883 !  9.236   !  -7.1
  cc2(18)  = 0D0     !  -2.5265 !  0.0
  crloc(18)= 0.4D0     !  0.6     !  0.4
!  nrow(18) = 3
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
!  nrow(20) = 4
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
!  nrow(29) = 4
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
!  nrow(47) = 4
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
  
!       hydrogen Gianocci PsP (approx. local Goedecker)
  
  amu(1)   = 1D0
  ch(1)    = 1D0
!  nrow(1)  = 1
  
!       Na (approximate local pseudo built on Goedecker local part)
  
  amu(11)  = 23.0D0
  ch(11)   = 1D0
  cc1(11)  = 1.5D0 ! 1.55 ! 2.5 ! 3.0   !  3.1
  cc2(11)  = 0D0
  crloc(11)= 0.9D0
!  nrow(11) = 3
!        write(6,*)'params pseudo, c1,c2',cc1(18),cc2(18),crloc(18)
  
!       Ar (approximate local pseudo built on Goedecker local part)
  
  amu(18)  = 39.95D0
  ch(18)   = 8D0
!        cc1(18)  = 60.37    ! 13.375   !  2.5034  ! 2.482
!        cc2(18)  = -8.972   ! -3.4762  ! -1.4664  ! -1.4526
!        crloc(18)= 0.4      ! 0.6      ! 0.9      ! 0.9
  
  cc1(18)  =   2.482D0   ! 2.5034  !
  cc2(18)  =  -1.4526D0  ! -1.4664  !
  crloc(18)=  0.9D0
  
  
!        cc1(18)  = 4.9558   !  9.9117
!        cc2(18)  = -2.0318  ! -4.0637
!        crloc(18)= 0.8
!09        cc1(18)  = 3.85
!09        cc2(18)  = -2.44
!09        crloc(18)= 0.9
!  nrow(18) = 3
!        write(6,*)'params pseudo, c1,c2',cc1(18),cc2(18),crloc(18)
ELSE IF((ipsptyp == 3).OR.(ipsptyp ==4)) THEN
!Goedecker read from two files, with valence electrons only
  do iel=1,92
    amu(iel)=0D0   ! masses as a function of atomic number
  enddo

!       read atomic number and masses of the elements

  OPEN(UNIT=96,STATUS='OLD',FORM='FORMATTED',FILE='periodic')
  do i=1,86
     read(96,*) natom,symb(natom),amu(natom)
  enddo
  close(unit=96)
!write(6,*) 'periodic table read from file periodic'

!       read pseudopotential parameters

  OPEN(UNIT=96,STATUS='OLD',FORM='FORMATTED',FILE='goed.asci')

!       beginning of reading loop

  icountt=-1
  iskip=-1

!       start of loop

open(unit=19,status='scratch')
!~ 50	continue
icountt=0
ch(0)=0D0
crloc(0)=0D0
cc1(0)=0D0
cc2(0)=0D0
r0g(0)=0D0
r1g(0)=0D0
r2g(0)=0D0

do ! infinite loop
!we read a new line
  read(96,'(a)',end=1000) a     ! exits loop at EOF
  if(a(1:1).ne.' ') then
!we have a new element
    inew=0
    icountt=icountt+1
    
    ch(icountt)=0D0
    crloc(icountt)=0D0
    cc1(icountt)=0D0
    cc2(icountt)=0D0
    r0g(icountt)=0D0
    r1g(icountt)=0D0
    r2g(icountt)=0D0
    l=0
    rewind 19
    write(19,'(a)') a   ! use a temporary file to re-read the line
    rewind 19
    read(19,'(a3,i2,5(f10.6))') naml,nval,rloc,c1,c2,c3,c4
    if(naml(3:3).eq.'1') then  ! new element is a semi-core
!   write(6,*) 'semi-core'
!   write(6,103) naml(1:3),' ',nval,rloc,c1,c2,c3,c4
!103   format(2A,i2,5(f10.6))
      isave=icountt-1 ! semi-core not counted as en element: remember previous element 
      icountt=99 
      iskip=99

!we skip semicore pseudos excepted if asked 

      IF(ipsptyp == 4) THEN ! Do not skip
        ! previous element will be replaced by the new semi-core version. 
        icountt=isave  
        iskip=icountt
      ENDIF
    else
      iskip=icountt
      rewind 19
      read(19,'(a2,i2,5(f10.6))') namc,nval,rloc,c1,c2,c3,c4
  !    write(6,103) namc(1:2),' ',nval,rloc,c1,c2,c3,c4
    endif
    ! stores properties of new element with index 'icountt'
    ch(icountt)=nval
    crloc(icountt)=rloc
    cc1(icountt)=c1
    cc2(icountt)=c2
  !       we search for the element in the periodic table
    nactual=0
    do iel=1,86
      if(naml(1:2).eq.symb(iel)) then
!        write(6,*) naml(1:2),' found in periodic table, Z=',iel,icountt
        nactual=iel
      endif
    enddo    
!    if(nval.eq.nactual) write(6,*) 'all electron pseudopotential'
    if(nactual.eq.0) stop 'not found in periodic table'
    if(icountt.ne.99) then ! it is not a semicore
      if(nactual.ne.icountt) then  ! the element is not on it usual location in the periodic table
        write(6,*)  'nactual',nactual,'icount',icountt
        stop 'lost in pseudop. coef table'
      endif
    endif
  else
!       we have nonlocal coefficients in this line
    rewind 19
    write(19,'(a)') a(6:80)
    rewind 19
    read(19,'(4(f10.6))') crr,h11,h22,h33
!    write(6,'(A,4(f10.6))') 'nonlocal reread',crr,h11,h22,h33
    if(iskip.ne.99) then
!    if(crr.eq.0D0) l=l-1
      if(l.eq.0) then
        r0g(icountt)=crr
        h0_11g(icountt)=h11
        h0_22g(icountt)=h22
        h0_33g(icountt)=h33
      endif
      if(l.eq.1) then
        r1g(icountt)=crr
        h1_11g(icountt)=h11
        h1_22g(icountt)=h22
      endif
      if(l.eq.2) then
        r2g(icountt)=crr
        h2_11g(icountt)=h11
      endif
      l=l+1  ! nonlocal coefficients line counter
    endif
  endif

!       go back to read a new line
  if(icountt.eq.99) icountt=isave  ! if it was a skiped Semi-core, continue from the previous element
end do !next line !

1000 write(6,*) 'all parameters have been read in file goed.asci'

close(unit=19)
ipsptyp =1
ELSE
  STOP ' IPERIO: this type PsP not yet implemented '
END IF

!     reset ionic masses if asked for


WRITE(6,*) 'resettig ionic masses...'
IF(ifredmas == 1) THEN
!  WRITE(6,*) 'in progress...'
  IF (ipsptyp == 0) THEN
     DO i=-92,92
        !           amu(i)=0.0484
        amu(i)=0.5D0
     END DO
  ELSE
     amfac = 1D0/amu(6)
     DO i=-92,92
        amu(i) = amu(i)*amfac
        amu(i)=0.5D0
     ENDDO
  ENDIF
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
USE util, ONLY:getcm,givetemperature,rotatevec3D
IMPLICIT NONE
CHARACTER (LEN=3) :: orderxyz
REAL(DP) :: dd, distmax, xcm, ycm, zcm, ctmpx, ctmpy, ctmpz
REAL(DP) :: optis, dgrid, sumion, v0, rnorm, tempv, xm
REAL(DP) :: vxn02,vyn02,vzn02
REAL(DP) :: vecin(3),vecout(3),vecalpha(3),totvalec=0D0
INTEGER :: i,ii,inx,inxg,ion,iunit,n,nxopti

INTEGER ::  igrid(7)
data igrid /32,48,64,72,96,128,160/

REAL(DP),EXTERNAL :: energ_ions

IF(nion2 == 0) THEN
  ecorr = 0D0
  RETURN       ! this is the case of jellium
END IF

WRITE (6,*) 'Entering initions()'



!     readings and basic transformation are only done for node "0"
! why ? (fc)

#if(parayes)
!IF(myn == 0)THEN
#endif
  
  
  OPEN(UNIT=9,STATUS='old',FORM='formatted', FILE='for005ion.'//outnam)
  
  
  WRITE(6,*) 'Reading positions...'
  distmax=0D0
  DO ion=1,nion
    
    IF(init_lcao == 1) THEN
      READ(9,*) cx(ion),cy(ion),cz(ion),np(ion),orderxyz, radini(ion)&
           ,ipol(ion)
        dd=sqrt(cx(ion)*cx(ion)+cy(ion)*cy(ion)+cz(ion)*cz(ion))
        if(dd.gt.distmax) distmax=dd

        totvalec=totvalec+ch(np(ion))
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
        totvalec=totvalec+ch(np(ion))
    END IF
    WRITE(6,*) ' ion,initord=',ion,initord(:,ion)    
    
!       initial kinetic momenta=0
    
    cpx(ion)=0D0
    cpy(ion)=0D0
    cpz(ion)=0D0
    
  END DO
  CLOSE(UNIT=9)
  
  WRITE(6,*) 'total number of valence electrons',totvalec
  IF(nclust<=0) THEN
    nstate=totvalec+nclust
    nclust=nstate
    nspdw=nclust/2
    WRITE(6,*) 'nstate set to',nstate
    WRITE(6,*) 'kstate*knode is',kstate*knode
    WRITE(6,*) 'nclust set to',nclust
    WRITE(6,*) 'nspdw set to',nspdw
    IF(nstate>kstate*knode) STOP 'kstate too small - increase' 
  END IF

  
  
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
    DO i=1,nion
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
    CALL getcm(1,0,0)  !  c.m. now on 'rvectmp(1:3)'
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
    
    !WRITE(*,*) ' rotate ions by anglex,y,z=',rotclustx,rotclusty,rotclustz
    vecalpha(1)=rotclustx/180D0*pi
    vecalpha(2)=rotclusty/180D0*pi
    vecalpha(3)=rotclustz/180D0*pi
    !WRITE(*,*) ' rotate ions by anglex,y,z=',vecalpha
    
! temporary shift to CM
    ctmpx=SUM(cx(1:nion))/nion
    ctmpy=SUM(cy(1:nion))/nion
    ctmpz=SUM(cz(1:nion))/nion
    cx=cx-ctmpx
    cy=cy-ctmpy
    cz=cz-ctmpz
    
! apply rotation   
    DO ion=1,nion

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

    !WRITE(*,*) 'ionic positions now:'
    !DO ion=1,nion
    !  WRITE(*,*) cx(ion),cy(ion),cz(ion)
    !END DO
    
  END IF
!     rotation completed
  
  tnonlocany = .false.  
  optis=1D10
  DO ion=1,nion
      tblock(ion) = .FALSE.
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
    WRITE(6,'(a,i4,a,3(g12.4),i4,1x,2(1pg13.5))')  &
      ' ion nr.',ion,':  x,y,z,type,params=',  &
        cx(ion),cy(ion),cz(ion),np(ion),amu(np(ion)),crloc(np(ion))
    IF(np(ion) > 92) STOP 'element out of range'
    IF(amu(np(ion)) == 0D0) STOP 'unknown elem. found'


    ! set flag for non-local PsP

    IF(ipsptyp==1) THEN
      tnonloc(ion) = (ABS(h2_11g(np(ion))) + ABS(h1_22g(np(ion))) &
                    + ABS(h0_33g(np(ion))) + ABS(h1_11g(np(ion))) &
                    + ABS(h0_22g(np(ion))) + ABS(h0_11g(np(ion)))) &
                    .GT. small
      IF(tnonloc(ion)) tnonlocany = .true.
    END IF

    
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
          'h0_11g=',h0_11g(np(ion)), 'h0_12g=',h0_12g(np(ion)), &
          'h0_22g=',h0_22g(np(ion)),  &
          'h0_33g=',h0_33g(np(ion)), 'h1_11g=',h1_11g(np(ion)),  &
          'h1_22g=',h1_22g(np(ion)), 'h2_11g=',h2_11g(np(ion)),  &
          'sgm1=',sgm1(np(ion)), 'sgm2=',sgm2(np(ion)),  &
          'chg1=',chg1(np(ion)), 'chg2=',chg2(np(ion)),  &
          'r0g=',r0g(np(ion)),  &
          'r1g=',r1g(np(ion)), 'r2g=',r2g(np(ion)),  &
          'radiong=',radiong(np(ion))
      WRITE(iunit,*) 'tblock=',tblock(ion)
      IF(ipsptyp == 1) WRITE(iunit,*) 'tnonloc=',tnonloc(ion)
    END DO
    
    IF(ipsptyp == 1) THEN
      dgrid = crloc(np(ion))/0.8493218D0
      WRITE(7,*) ' optimal grid spacing for this element=',dgrid
      WRITE(6,*) ' optimal grid spacing for this element=',dgrid
      if(dgrid.lt.optis) optis=dgrid
    END IF
  END DO
      WRITE(7,*) 'dx=',dx,' absolute optimal grid spacing =',optis
      WRITE(6,*) 'dx=',dx,' absolute optimal grid spacing =',optis
           OPEN(119, FILE='dx', STATUS='UNKNOWN')
           write(119,*) optis
           close(119)

           nxopti=(distmax+1D0)/optis
           nxopti=(nxopti/4)*8
           inx=0
           inxg=0
           do while(inx<=nxopti) 
               inxg=inxg+1
               if(inxg.lt.7) then
                inx=igrid(inxg)
               else
                inx=nxopti
              endif
           ENDdo
      WRITE(6,*) 'distmax of ions',distmax,'recommended nx ny nz',inx
           OPEN(120, FILE='nx', STATUS='UNKNOWN')
           write(120,*) inx
           close(120)
  
  IF (tempion > 0D0 .AND.imob /= 0) THEN
    CALL givetemperature(cpx,cpy,cpz,nion,tempion, amu(np(1))*1836.0D0*ame,4)
  END IF
  
#if(parayes)
!END IF                                             ! myn=0
#endif

!     Part for node "0" finished. Now distribute to all nodes.

#if(parayes)
CALL comm_ionconfig()
#endif


!        check consistency of ions  (obsolete?)


!     check consistency 'nrow' <--> nr. of projectors



!       np(0) for use in Monte-Carlo


IF(icooltyp == 3)  CALL cenmass()
np(0) = np(1)
dt12=dt1




!       initialize pseudopotential background

!g         call calcpseudo(rho)

IF(ipsptyp == 1) THEN
  DO ion=1,nion
    CALL calc_proj(cx(ion),cy(ion),cz(ion),cx(ion),cy(ion),cz(ion),ion)
    IF (iswitch_interpol==1) CALL calc_projFine(cx(ion),cy(ion),cz(ion),cx(ion),cy(ion),cz(ion),ion)
  END DO
  IF (iswitch_interpol==1) CALL mergetabs
END IF


!     background energy


WRITE(6,*) 'Calculating ionic energy ...'

sumion = energ_ions()     !  ???
WRITE(7,'(a,f17.4)') 'sumion=',sumion
WRITE(6,'(a,f17.4)') 'sumion=',sumion
WRITE(7,*)
WRITE(6,*)
!k
ecorr=sumion           !analytical case   ??
!k


!     initialize common velocity

IF (ekin0pp > 0D0) THEN
  v0 = SQRT(2D0*ekin0pp/(amu(np(nion))*1836.0D0*ame))
  rnorm = vxn0**2 + vyn0**2+ vzn0**2
  rnorm = SQRT(rnorm)
  IF (rnorm == 0) STOP 'Velocity vector not normalizable'
  vxn02 = vxn0/rnorm*v0*ame
  vyn02 = vyn0/rnorm*v0*ame
  vzn02 = vzn0/rnorm*v0*ame
  xm = amu(np(nion))*1836.0D0
  DO i=1,nion
    cpx(i)=vxn02*xm
    cpy(i)=vyn02*xm
    cpz(i)=vzn02*xm
  END DO
  WRITE(*,*) ' EKIN0PP: initial CP:',cpx(1),cpy(1),cpz(1)
END IF ! initial kinetic energy for ions


!     Initialization of the projectile velocity
!     In the case of an atom, the w.f. are boosted accordingly in tinit
! lionel : np(nion)=> np(nproj)
IF (eproj > 0D0 .AND. taccel<1D-5) THEN
  v0 = SQRT(2D0*eproj/(amu(np(nproj))*1836.0D0*ame))
  rnorm = vpx**2 + vpy**2+ vpz**2
  rnorm = SQRT(rnorm)
  IF (rnorm == 0) STOP 'Velocity vector not normalizable'
  tempv=v0*ame*amu(np(nproj))*1836.0D0/rnorm
  cpx(nproj) = vpx*tempv
  cpy(nproj) = vpy*tempv
  cpz(nproj) = vpz*tempv
  WRITE(*,*) ' projectile momenta initialized: eproj,v0=',eproj,v0
ENDIF

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
IMPLICIT NONE

INTEGER :: ind
REAL(DP):: a20fac,a22fac,gamarg,xclust
REAL(DP):: srms,sum1
!------------------------------------------------------------------

IF(nion2 /= 0) RETURN              ! case of detailed ions

!     factors for transformation from hill-wheeler coordinates

gamarg = gamma * pi / 180D0
a20fac = COS(gamarg)
a22fac = SIN(gamarg) / SQRT(2D0)

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
USE util, ONLY:shiftfield
IMPLICIT NONE
#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

REAL(DP), INTENT(IN OUT)                     :: psir(kdfull2,kstate)
INTEGER ::i,nr,nbe,nbr
REAL(DP)::en
REAL(DP), ALLOCATABLE :: rfieldaux(:)

!----------------------------------------------------------------------

!      h2m=hbar*hbar/2.0/ame
IF(myn == 0)THEN
  IF(nspdw > nclust/2) STOP 'nspdw must be less or equal nclust/2'
END IF


!     prepare book-keeping for wavefunctions

omeg=0.25D0*h2m
!      if(nclust.gt.0)
!     &    call ininqb(nclust,deocc,b2occ,gamocc*pi/180D0)

IF(ifhamdiag>0 .AND. nstate>nclust) THEN
  WRITE(6,'(2a,2i5)')  ' IFHAMDIAG>0 only allowed for NSTATE=NCLUST', &
   ' Presently: nstate,nclust=',nstate,nclust
  WRITE(7,'(2a,2i5)')  ' IFHAMDIAG>0 only allowed for NSTATE=NCLUST', &
   ' Presently: nstate,nclust=',nstate,nclust
  STOP  ' IFHAMDIAG>0 only allowed for NSTATE=NCLUST'
END IF
WRITE(*,*) 'before ininqb'
CALL ininqb(nclust,deocc,b2occ,gamocc*pi/180D0)
WRITE(*,*) 'after ininqb'


!     initialize H.O. wavefunctions

IF(init_lcao /= 1) THEN
   CALL initho(psir)   
   CALL initho(psir)

   
!       remove symmetry restrictions
  
!  IF(init_lcao /= 1) CALL rem_sym(psir)
!  CALL rem_sym(psir)
  WRITE(7,'(a)')'after ordering:'
  WRITE(7,'(a,2i3)') 'nstate,nclust',nstate,nclust
  WRITE(6,'(a)')'after ordering:'
  WRITE(6,'(a,2i3)') 'nstate,nclust',nstate,nclust
  IF(ispinsep /= 0) CALL spinsep(psir)
  
  WRITE(7,'(a)') 'after oscillator initialization:'
  DO i=1,nstate
    nr=nrel2abs(i)
!    IF (myn == 0) THEN
      WRITE(7,'(a,i3,a,f12.4,a,i2,a,3i3)')'wf ',i,  &
          ' occ',occup(i),' sp',ispin(nr), ' knots:',nq(1,nr),nq(2,nr),nq(3,nr)
!    END IF
  END DO

ELSE

!       optionally LCGO initialization
  
!  CALL ininodes()
  occup=1D0
#if(parano) 
  nstate_all=nclust
  nstate=nclust
#endif
  CALL genermowf(psir,nmaxst)
  WRITE(6,'(a)') 'after LCAO initialization:'
  WRITE(7,'(a)') 'after LCAO initialization:'
  DO nbr=1,nstate
!    en=0D0
!    DO i=1,kdfull2
!      en=en+psir(i,nbr)*psir(i,nbr)
!    END DO
!    en=en*dx*dy*dz
    en=SUM(psir(:,nbr)**2)*dx*dy*dz
    WRITE(6,'(2(a,i5),3(a,1pg13.5))') 'node=',myn,', state=',nbr,  &
        ', en=',en,', occ=',occup(nbr),', spin=',ispin(nbr)
    WRITE(7,'(2(a,i5),3(a,1pg13.5))') 'node=',myn,', state=',nbr,  &
        ', en=',en,', occ=',occup(nbr),', spin=',ispin(nbr)
  END DO
END IF




!     optional shift of wavefunctions

#if(simpara)
WRITE(7,*) ' SHIFTFIELD overridden'
#else
IF(ABS(shiftwfx)+ABS(shiftwfy)+ABS(shiftwfz) > 1D-20) THEN
  DO nbe=1,nstate
    CALL shiftfield(psir(:,nbe),shiftwfx,shiftwfy,shiftwfz)
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
USE util, ONLY:pair
IMPLICIT NONE
#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

INTEGER, INTENT(IN)                      :: nelect
REAL(DP), INTENT(IN)                         :: deoccin
REAL(DP), INTENT(IN)                         :: betain
REAL(DP), INTENT(IN)                     :: gamin

!     initialization of book-keeping arrays of states nq, ispin.

!     estimates Fermi energy from particle number and considers
!     all states up "Fermi-energy plus deoccin" in the ordering
!     of a deformed harmonic oscillator for given deformation.
!     the input parameters are:
!       nelect  = number of electrons
!       deoccin   = number of osc. shells above Fermi shell
!       betain  = quadrupole deformation of jellium background
!       gamin   = triaxiality angle of jellium background
!       temp    = (via common) temperature, temp=0 cuts to occupied only

!     the output goes on the occupation fields nq and ispin on
!     common /option/ .

!INTEGER, PARAMETER :: kmxsav=kdfull/3
REAL(DP) :: esp(ksttot)           ! storage for s.p. energies
REAL(DP) :: efacto                ! factor to get energies from h.o.
!~ REAL(DP) :: efermi                ! estimate for fermi shell
REAL(DP) :: q20fac                ! sqrt(5/16pi)
REAL(DP) :: cosfac,sinfac         ! weightes deduced from 'gamin'
!     real      xfac,yfac,zfac        ! effective osc. energies in x,y,z
REAL(DP) :: speact                ! actual s.p. energy in loop
!~ INTEGER :: noscmx                ! maximum oscillator number
INTEGER :: n                     ! nr. of state
INTEGER :: noscx,noscy,noscz     ! osc. nr. in each direction
INTEGER :: nomxup, nomxdw, nmaxdiff, mspindw
INTEGER :: i,isav,j,k

#if(parayes)
INTEGER :: nstpernode, nstaccum, nod, nodeplus
#endif

REAL(DP),ALLOCATABLE :: ph(:)             ! degeneracy of wavefunction, for
REAL(DP) :: occu(ksttot)
REAL(DP) :: ecutdw, ecutup, efrmup, efrmdw,nelup,neldw
REAL(DP) :: gp,partnm,epstmp,sav
LOGICAL :: tocc
DATA tocc/.false./

REAL(DP),PARAMETER :: third=1D0/3D0

!-----------------------------------------------------------------------

ALLOCATE(ph(2*ksttot))
IF(numspin==2) THEN
  ph = 1D0
ELSE IF(numspin==1) THEN
  ph=0D0
  ph(1:ksttot) = 2D0
ELSE
  STOP "invalid value for input parameter NUMSPIN"
END IF


!     prepare initial parameters
!     estimate of Fermi energy relies on spherical oscillator shells.
!     the shell label N is determined by

!     nelect,spin = N*(N+1)*(N+2)/6

!     where 'nelect,spin' is the nr. of electrons for given spin.
!     this relation is resolved approximately for N.

q20fac = SQRT(5D0/(16D0*pi))
IF(numspin==2) THEN
  nelup  = nelect-nspdw
  neldw  = nspdw
  !test
  WRITE(6,*) 'nelup',nelup
  WRITE(6,*) 'neldw',nspdw
  !test
ELSE
  IF(MOD(nelect,2) == 1)  &
    STOP ' nr. of electrons must be even for spin degeneracy'
  nelup  = nclust/2
END IF
efacto = 0.25D0/(1D0*nelect)**third
efrmup = (6D0*nelup)**third
efrmup = efrmup/(1D0-1D0/(efrmup*efrmup))**third-1.5D0
ecutup = efrmup+deoccin
nomxup = ecutup*(1D0+2D0*q20fac*betain)+0.5D0
ecutup = efacto*ecutup
IF(numspin==2) THEN
  IF(neldw > 0) THEN
    efrmdw = (6D0*neldw)**third
    efrmdw = efrmdw/(1D0-1D0/(efrmdw*efrmdw))**third-1.5D0
  ELSE
    efrmdw = -0.00001D0
  END IF
  ecutdw = efrmdw+deoccin
  nomxdw = ecutdw*(1D0+2.0D0*q20fac*betain)+0.5D0
  ecutdw = efacto*ecutdw
END IF
cosfac = q20fac*COS(gamin)
sinfac = q20fac*SQRT(2D0)*SIN(gamin)
xfac   = efacto/(1D0-betain*(cosfac-sinfac))
yfac   = efacto/(1D0-betain*(cosfac+sinfac))
zfac   = efacto/(1D0+2D0*betain*cosfac)
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
!old?        IF(n >= (nelup+nmaxdiff)) THEN
        IF(n > (nelup+nmaxdiff)) THEN
          WRITE(6,*) 'n=',n,' greater than nelup+dnmax=',nelup+nmaxdiff
          n=n-1                              ! new from MD
          EXIT lb1
        END IF
        WRITE(6,'(a,5i5,3f10.3)') 'init states: n,...=',&
          n,ispin(n),noscx,noscy,noscz,occu(n),speact,ecutup
      END IF
     END DO
     END DO
     END DO lb1
!19   continue
!test
WRITE(6,*) 'n(nelup)',n
!test
IF(n < nelup) STOP ' not enough states to reach (spin-up) particle number'
IF(numspin==2) THEN

!     select spin down states
!     distribute preliminary occupation numbers

  mspindw = 0
  lb2: DO noscz=0,nomxdw
       DO noscy=0,nomxdw
       DO noscx=0,nomxdw
        speact = noscz*zfac+noscy*yfac+noscx*xfac
        IF(nelup /= neldw) speact = speact*1D01    ! enhance for spin asymmetry
        IF(speact <= ecutdw) THEN
          n     = 1+n
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
        
!old?        IF(mspindw >= (neldw+nmaxdiff)) THEN
          IF(mspindw > (neldw+nmaxdiff)) THEN
             WRITE(6,*) 'mspindw=',mspindw,&
                ' greater than neldw+dnmax=',neldw+nmaxdiff
            mspindw=mspindw-1                                   ! new by MD
            EXIT lb2
          ENDIF
          WRITE(6,'(a,3i5,3f10.3)') 'check spin down:', &
                 ispin(n),n,mspindw,occu(n),speact,ecutdw
          WRITE(7,*)
          WRITE(7,*) 'speact',speact,'n',n,'isp',ispin(n)
          WRITE(7,*) 'no',noscx,noscy,noscz,zfac,yfac,xfac
        END IF
       END DO
       END DO
       END DO lb2
!test
WRITE(6,*) 'mspindw',mspindw
WRITE(6,*) 'neldw',neldw
!test
  IF(mspindw < neldw)  &
      STOP ' not enough states to reach spin-down particle number'
END IF

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

DEALLOCATE(ph)


!
!  initialize book-keeping fields assoxiating states with nodes

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
IF (init_lcao==1) nstate_all = nclust !lionel
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
  DO nod=0,knode-1
    nrel2abs_other(i,nod) = i+nstart_node(nod)
  END DO
END DO
DO i=1,ksttot
  DO nod=knode-1,0,-1
    IF(i > nstart_node(nod)) THEN
      nabs2rel(i)=i-nstart_node(nod)
      nhome(i)=nod
      ispin_node(nabs2rel(i),nod) = ispin(i)
      EXIT
    END IF
  END DO
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
USE util, ONLY:rotatevec3D
IMPLICIT NONE

!     computes the jellium background such that the radius is
!     given by 'radjel'*part.numb.**1/3, the surface by
!     'surjel', and the deformation is adjusted to the
!     electron deformation given on 'alphel' and 'betael'

!     y40 is expressed in terms of y20_effective to turn it automatically
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

INTEGER :: i1,i2,i3,ii,iter
REAL(DP) :: vecin(3),vecout(3),vecalpha(3),vecalp(3)
REAL(DP) :: anglex,angley,anglez,alpha,beta,thetax,thetay,thetaz,theta0
REAL(DP) :: alphael,beta2j,gamma2j,gammaj
REAL(DP) :: onetrd,d4pi,alpfac,q20fac,q22fac,q40fac
REAL(DP) :: alphbk,betabk,hexabk,bet2EF
REAL(DP) :: deralp,derbet,delalp,delbet
REAL(DP) :: y2eff,y20,y22,y40,y20fac,y20obs,y22fac,y22obs,q40obs,q20red,q40red
REAL(DP) :: radius,argum,rho0,rhoc,precis,volel
REAL(DP) :: dpm,qoct,sqhe,sqt,sqq,sqn
REAL(DP) :: astep
REAL(DP) :: r,reff,rh,rr
REAL(DP) :: x,y,z,xac,yac,zac,xact,yact,zact,xx,yy,zz
DATA astep/0.6D0/
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
      ' jellium background rotated by anglex,angley,anglez=', vecalpha*180D0/PI    
  WRITE(7,'(/a,3f8.3)')  &
      ' jellium background rotated by anglex,angley,anglez=', vecalpha*180D0/PI    
  
ELSE IF(iturn == 0) THEN
!  case for static iteration
  anglex  = 0D0
  angley  = 0D0
  anglez  = 0D0
END IF

!      itmax = itback

!     dimensionless deformation parameters
!      zero = 0.0
onetrd = 1D0 / 3D0
d4pi   = 1D0/(4D0*pi)
alpfac = 4D0*pi/5D0
q20fac = SQRT(5D0/(16D0*pi))
q22fac = SQRT(15D0/(8D0*pi))
q40fac = SQRT(9D0/(4D0*pi))
!----------------------------------------------------------------------

!     set initial values for iteration of deformation
alphbk = alphel
betabk = betael
hexabk = hexel


!     effective angle 'gamma' to determine principle axis for y40

bet2ef = (alphbk*alphbk+2D0*betabk*betabk)
IF(bet2ef /= 0D0) THEN
  bet2ef  = 1D0/SQRT(bet2ef)
  y20fac  = alphbk*bet2ef
  y20obs  = y20fac*q20fac
  y22fac  = betabk*bet2ef
  y22obs  = y22fac*q22fac
  q40red = 7D0*SQRT(pi/9D0)          ! factor for y40 from y20^2
  q20red = SQRT(5D0/pi)/7D0          ! cofactor on y20 in y40 from
  q40obs = 8D0*q40red/q40fac
END IF

radius = radjel * partn**onetrd
rhoc   = 3D0*d4pi / radjel**3

argum  = radius / surjel
IF(argum < 38D0) THEN
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
        IF(abs(rotclustx)+abs(rotclusty)+abs(rotclustz)>0D0) THEN    
         vecalp(1)=rotclustx/180D0*pi                                
         vecalp(2)=rotclusty/180D0*pi                                
         vecalp(3)=rotclustz/180D0*pi                                
         vecin(1) = xact                                             
         vecin(2) = yact                                             
         vecin(3) = zact                                             
         CALL rotatevec3D(vecin,vecout,vecalp)                       
         xac = vecout(1)                                             
         yac = vecout(2)                                             
         zac = vecout(3)                                             
        ELSE                                                         
         xac = xact                                                  
         yac = yact                                                  
         zac = zact                                                  
        END IF                                                       
        IF(iturn == 1) THEN
!          CALL rotxyz(xact,yact,zact,x,y,z,anglex,angley,anglez)
          vecin(1) = xac    
          vecin(2) = yac    
          vecin(3) = zac    
          CALL rotatevec3D(vecin,vecout,vecalpha)
          x = vecout(1)
          y = vecout(2)
          z = vecout(3)
        ELSE IF(iturn == 0) THEN
          x=xac             
          y=yac             
          z=zac             
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
            y40   = q40fac*(8D0*zz*zz+3D0*xx*xx+ 3D0*yy*yy-24D0*zz*xx-  &
                24D0*zz*yy+6D0*xx*yy)/ (8D0*rr*rr)
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
            -(alphbk*alphbk+2D0*betabk*betabk+hexabk*hexabk)*d4pi)
        argum = (r-reff)/surjel
        IF(argum > +38D0) THEN
          rhojel(ii) = 0D0
        ELSE IF(argum < -38D0) THEN
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
  theta0 = (thetax+thetay+thetaz)/3D0
  thetax = thetax/theta0
  thetay = thetay/theta0
  thetaz = thetaz/theta0
  alpha  = alpfac*q20fac*sqq/(sqn*sqr*sqr)
  beta   = alpfac*0.5D0*q22fac*sqt/(sqn*sqr*sqr)
  deralp = 1D0
  derbet = 1D0
  delalp = -astep*(alpha-alphel)/deralp
  delbet = -astep*(beta -betael)/derbet
  alphbk = alphbk+delalp
  betabk = betabk+delbet
  precis = (alpha-alphel)**2+(beta-betael)**2+(partn-sqn)**2
  
  IF((precis < endcon).AND.(ABS(partn-sqn) < endcon)) EXIT
  
  radius = radius * (partn/sqn)**onetrd
  argum  = radius / surjel
  IF(argum < 38D0) THEN
    rho0 = rhoc * (1D0 + EXP(-argum))
  ELSE
    rho0 = rhoc
  END IF
!      write(7,1000)    iter,precis,sqn,dpm,sqr,sqq,sqt,
!     &        sqhe,alpha,beta,(alphbk-delalp),(betabk-delbet),alphel,
!     &            betael,deralp,derbet,delalp,delbet
!      write(7,*) 'thetax=',thetax
  
END DO

IF(iter== itback+1 )THEN
  WRITE(7,'(a/a,g11.3)') ' ---> background deformation did not converge!',  &
      '      residual error in"alpha"=',precis
  STOP ' no convergence in "jelbak"'
END IF

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
beta2j = SQRT(alpha*alpha+2D0*beta*beta)
gammaj = ATAN(1.4142136D0*beta/alpha)*180D0/pi
WRITE(7,*) 'effective jellium deformations: beta2j=',beta2j,  &
    '  gammaj=',gammaj


RETURN
END SUBROUTINE jelbak

!-----initho------------------------------------------------------initho

SUBROUTINE initho(psir)

!     initializes harmonic oscillator wavefunctions

USE params
IMPLICIT NONE

#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif


REAL(DP), INTENT(OUT)                        :: psir(kdfull2,kstate)
REAL(DP) :: valx(nx2),valy(ny2),valz(nz2)
!~ REAL(DP), ALLOCATABLE :: phix(:)
INTEGER :: i,ii,inx,iny,inz,ix,iy,iz,nb
REAL(DP) :: an,bk1,bxx,bxy,bxz,hom,homx,homy,homz
REAL(DP) :: vx,vy,vz

REAL(DP),PARAMETER :: third=1D0/3D0
REAL(DP),PARAMETER :: sixth=1D0/6D0

!EQUIVALENCE (phix(1),w1(1))

!----------------------------------------------------------------------

!     check workspace


an  = REAL(2*nclust,DP)
IF(temp > 0D0) THEN
  homx = omeg*an**(-third)*xfac
  homy = omeg*an**(-third)*yfac
  homz = omeg*an**(-third)*zfac
  bxx  = (2D0*h2m/homx)**3
  bxy  = (2D0*h2m/homy)**3
  bxz  = (2D0*h2m/homz)**3
  bxx  = osfac*(bxx**sixth)
  bxy  = osfac*(bxy**sixth)
  bxz  = osfac*(bxz**sixth)
ELSE
  hom  = omeg*an**(-third)
  bk1  = (2D0*h2m/hom)**3
  bk1   = osfac*(bk1**sixth)
END IF


DO nb=1,nstate
  
!       number of knots
  
  
!       nq is  relative to the proc
  
!       this way we implicitly parallelize the computation of phix
  
  inx=nq(1,nrel2abs(nb))
  iny=nq(2,nrel2abs(nb))
  inz=nq(3,nrel2abs(nb))
  
  DO  i=1,3
    ipar(i,nb)=1-2*MOD(nq(i,nrel2abs(nb)),2)
  END DO
  
  IF(temp > 0D0) THEN
    CALL clust(inx,bxx,0D0,xval,valx,nx2)
    CALL clust(iny,bxy,0D0,yval,valy,ny2)
    CALL clust(inz,bxz,0D0,zval,valz,nz2)
  ELSE IF(temp <= 0D0) THEN
    CALL clust(inx,bk1,0D0,xval,valx,nx2)
    CALL clust(iny,bk1,0D0,yval,valy,ny2)
    CALL clust(inz,bk1,0D0,zval,valz,nz2)
  END IF
  
!       composition of the factorised wave-function
!       occupies only upper box (1/8 part), but is returned on 'psir'
 
  
  ii=0
  DO iz=1,nz2
    vz=valz(iz)
    DO iy=1,ny2
      vy=valy(iy)
      DO ix=1,nx2
        vx=valx(ix)
        ii=ii+1
        psir(ii,nb)=vx*vy*vz
      END DO
    END DO
  END DO

  
END DO



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
!     where h(x) is the Hermite polynomial as explained in Rottmann, p.1

USE params
IMPLICIT NONE


INTEGER, INTENT(IN OUT)                  :: in
REAL(DP), INTENT(IN)                     :: b
REAL(DP), INTENT(IN)                     :: z
REAL(DP), INTENT(IN)                     :: x(*)
REAL(DP), INTENT(OUT)                        :: val(*)
INTEGER, INTENT(IN)                      :: n1

INTEGER :: i,j
REAL(DP) :: argum, coef, fak, tau, v
!-----------------------------------------------------------------------

fak = 1
DO i=2,in
  fak = fak * i
END DO
coef=1D0/SQRT(b*1.772454D0*(2D0**in)*fak)


DO j=1,n1
  
  tau=(x(j)-z)/b
  argum = -0.5D0*tau*tau
  IF(argum < -38D0) THEN
    val(j) = 0D0
  ELSE
    v=coef*EXP(argum)
    SELECT CASE(in)
      CASE(0)
        val(j)=v
      CASE(1)
        val(j)=v*2*tau
      CASE(2)
        val(j)=v*(4D0*tau*tau-2D0)
      CASE(3)
        val(j)=v*(8D0*tau**3-12D0*tau)
      CASE(4)
        val(j)=v*(16D0*tau**4-48D0*tau**2+12D0)
      CASE(5)
        val(j)=v*(32D0*tau**5-160D0*tau**3+120D0*tau)
      CASE(6)
        val(j)=v*8D0*(8D0*tau**6-60D0*tau**4+90D0*tau**2-15D0)
      CASE(7)
        val(j)=v*16D0*(8D0*tau**7-84D0*tau**5+210D0*tau**3-105D0*tau)
      CASE(8)
        val(j)=v*tau**8              ! ortho-normalize later
      CASE DEFAULT
        WRITE(6,'(a,1pg12.4)') ' wrong radial quantum number in clust: in=',in
        STOP 'wrong radial quantum number in clust'
    END SELECT
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
IMPLICIT NONE

REAL(DP),INTENT(IN) :: xin
REAL(DP),INTENT(IN) :: yin
REAL(DP),INTENT(IN) :: zin
REAL(DP),INTENT(OUT) :: xout
REAL(DP),INTENT(OUT) :: yout
REAL(DP),INTENT(OUT) :: zout
REAL(DP),INTENT(IN) ::anglex
REAL(DP),INTENT(IN) ::angley
REAL(DP),INTENT(IN) ::anglez

REAL(DP) :: x1,y1,z1,x2,y2,z2
REAL(DP) :: cphi,sphi
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



!-----spinsep-----------------------------------------------------------

SUBROUTINE spinsep(psir)

!     to induce local spin current by shifting spin-up and down in
!     different directions

USE params
IMPLICIT NONE


REAL(DP), INTENT(OUT)                        :: psir(kdfull2,kstate)

INTEGER :: ii,ix,iy,iz,nb
REAL(DP) :: sgeps,x1,y1,z1

!----------------------------------------------------------------------

IF(numspin==2) THEN
  WRITE(*,*) ' SPINSEP invoked'
  DO nb=1,nstate
    sgeps = (3-2*ispin(nb))*0.01D0
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
END IF
RETURN
END SUBROUTINE spinsep

!-----fixion------------------------------------------------------------

SUBROUTINE fixion

!     to fix two layers along y and z of an argon cluster (only cores)
!     to simulate the bulk
!     the input file must have already the positions of the bulk

USE params
IMPLICIT NONE

INTEGER :: ind, ind1, ind2, ion
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
IMPLICIT NONE

!------------------------------------------------------------------

#if(parayes)
#if(findiff|numerov)
STOP ' parallele computing not active with finite differences'
#endif
IF(ifsicp==5) STOP ' exact exchange not compatible with parallele code'
#endif

#if(symmcond)
IF(ifsicp==4) STOP ' propagated symm.cond. and KLI not compatible'
#endif
IF(.NOT.directenergy .AND. ifsicp==4) STOP " KLI requires directenergy=.true."
#if(!pw92)
IF(directenergy) STOP ' directenergy=.true. requires Perdew&Wang functional '
#endif
IF(directenergy .AND. ifsicp==5) &
   STOP ' directenergy=.true. not yet prepared for exact exchange '
#if(!raregas)
IF(ivdw /=0) STOP " set raregas=1 when using VdW"
#endif

#if(twostsic)
!IF(ifhamdiag>0 .AND. ifsicp==8) &
!  STOP "Hamiltonian diagonalization presently not compatible with full SIC"
IF(ifhamdiag>0 .AND. ifsicp==7) &
  STOP "Hamiltonian diagonalization presently not compatible with localized SIC"
#endif

RETURN
END SUBROUTINE checkoptions


!-----init_homfield-------------------------------------------------

SUBROUTINE init_homfield()

!     initialize a homogeneous electrical field and
!     adds it to the background field 'potFixedIon'.

USE params
IMPLICIT NONE

INTEGER :: ind,ix,iy,iz
REAL(DP) :: sc,vhom,x1,y1,z1
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

SUBROUTINE init_grid()

!     initialize Coulomb solver, kinetic energy and other
!     grid properties

USE params
USE kinetic
#if(netlib_fft|fftw_cpu)
USE coulsolv
#endif
IMPLICIT NONE

#if(paraworld||parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

INTEGER :: nrank=0

INTEGER::coeff,level

INTEGER :: ix, iy, iz
REAL(DP) :: dxfine, dyfine, dzfine, dvolfine
REAL(DP) :: x1,y1,z1
!------------------------------------------------------------------
#if(paraworld)

CALL  mpi_comm_rank(mpi_comm_world,nrank,icode)
#endif

level=nrank

IF(level>=1) THEN
  coeff = 2**level
  dx=dx*coeff
  dy=dy*coeff
  dz=dz*coeff
  kxbox=kxbox/coeff
  kybox=kybox/coeff
  kzbox=kzbox/coeff
  WRITE(6,*) 'level',level,dx,kxbox
!    write(outnam,*) level 

  nx2=kxbox
  ny2=kybox
  nz2=kzbox
  kdfull2=nx2*ny2*nz2

  kxmax=kxbox/2+1;kymax=kybox/2+1;kzmax=kzbox/2+1
  nx2=kxbox;ny2=kybox;nz2=kzbox;nx2fine=2*nx2-1;ny2fine=2*ny2-1
  nxyz=nx2*ny2*nz2;nyf=nx2;nxyf=nx2*ny2;nxyfine=nx2fine*ny2fine
  kdfull2=kxbox*kybox*kzbox
  kdfull2fine=(2*kxbox-1)*(2*kybox-1)*(2*kzbox-1)
  kdfull2fine=(2*kxbox)*(2*kybox)*(2*kzbox)
  nx=nx2/2;ny=ny2/2;nz=nz2/2
  nxfine=kxbox;nyfine=kybox;nzfine=kzbox
  END IF
!#if(gridfft)
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

!#endif
#if(tfindiff|tnumerov)
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





 

ALLOCATE(xval(nx2),yval(ny2),zval(nz2))        !  grid coordinates
ALLOCATE(xt2(nx2),yt2(ny2),zt2(nz2))           !  coordinates**2

!     check and initialize

IF(myn == 0)THEN
 
  WRITE(6,*) 'kxbox,kybox,kzbox:', kxbox,kybox,kzbox
  WRITE(7,*) 'kxbox,kybox,kzbox:', kxbox,kybox,kzbox
  IF(kxbox*kybox*kzbox == 0) &
     STOP ' you must specify the box sizes KXBOX,KYBOX,KZBOX'
  IF(kstate == 0) STOP ' you must specify the maximum nr. of states KSTATE'
  IF(nabsorb > MIN(kxbox,kybox,kzbox)/4) STOP " NABSO too large"
END IF
!old      if(kxmax.lt.nx1) then
!tab         stop ' error in parameter'
!tab      elseif(kymax.lt.ny1) then
!tab         stop ' error in parameter'
!tab      elseif(kzmax.lt.nzi) then
!tab         stop ' error in parameter'
!tab      endif

dvol=dx*dy*dz
dxfine=dx/2D0
dyfine=dy/2D0
dzfine=dz/2D0
dvolfine=dvol/8


DO ix=1,nx2
  x1=(ix-nxsh)*dx
  xval(ix)=x1
  xt2(ix)=x1*x1
END DO

DO iy=1,ny2
  y1=(iy-nysh)*dy
  yval(iy)=y1
  yt2(iy)=y1*y1
END DO

DO iz=1,nz2
  z1=(iz-nzsh)*dz
  zval(iz)=z1
  zt2(iz)=z1*z1
END DO

!     init kinetic energy array

#if(findiff)
CALL inv3p_ini(dt1)
#endif
#if(numerov)
CALL inv5p_ini(dt1)
#endif
#if(gridfft)
CALL init_grid_fft(dx,dy,dz,nx2,ny2,nz2,dt1,h2m)
#endif

!     init Coulomb solver

!!$#if(findiff|numerov)
!!$CALL d3sinfinit (dx,dy,dz)
!!$#else
#if(gridfft)
#if(coufou || coudoub)
CALL init_coul(dx,dy,dz,nx2,ny2,nz2)
#endif
#endif


RETURN
END SUBROUTINE init_grid





!-----ithion--------------------------------------------------------

SUBROUTINE ithion

!       initial thermalization of the ions (only executed if needed)
!       we assign 0.5*kb*t per degree of freedom

USE params
IMPLICIT NONE

INTEGER :: ion
REAL(DP) :: bk, ek, ekin, pmoy, pn, vmoy, xm
!---------------------------------------------------------------------


bk  = 1D0          ! use temperature in units of Ry
ekin=0.5D0*bk*tempion
WRITE(7,'(a,f9.4)')  'wanted temp',tempion
WRITE (7,'(a,f9.4)') 'wanted energ  per degree of freedom',ekin
xm=0.5D0*1836.0D0*amu(np(1))*ame
WRITE(7,*) 'warning : ithion not yet able to treat unhomogeneous systems'

!     attention to masses !

ekin=ekin/3/nion*(3*nion-6D0)
vmoy=SQRT(3D0*2D0*ekin/xm)
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
  xm=0.5D0*1836D0*amu(np(1))*ame
  ek=ek/2D0/xm
  WRITE(7,'(a,i2,a,f9.4)') 'ion',ion,'ek=',ek
  ekion=ekion+ek
  WRITE(7,'(a,i2,3f9.4)') 'ion',ion,cpx(ion),cpy(ion),cpz(ion)
END DO
WRITE(7,'(a,f9.4)') 'kin.energ after renormalization',ekion
ekion=ekion/(3*nion-6D0)
WRITE(7,'(a,f9.4)') 'av.kin.energ per net degree of freedom',ekion
ekion=ekion*2D0/bk
WRITE(7,'(a,f9.4)') 'corresponding temperature',ekion



RETURN
END SUBROUTINE ithion


!-----transf ----------------------------------------------------------

SUBROUTINE transf

!     transform ions to the main axis of inertia

USE params
IMPLICIT NONE

INTEGER :: i, i1, i2, iax, idummy, ion, ion1, ion2, icoo, iswap
INTEGER :: j, k
REAL(DP) :: aviner, delmx, delx, dely, dist3, gamu, zpos
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
  
  aviner = 1D0/3D0 * (dminer(1) + dminer(2) + dminer(3))
  delmx  = dminer(3) - aviner
  iax    = 3
!   check whether x or y are outstanding
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
IMPLICIT NONE

REAL(DP), INTENT(IN OUT)                     :: a(np,np)
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: np
REAL(DP), INTENT(OUT)                        :: d(np)
REAL(DP), INTENT(OUT)                        :: v(np,np)
INTEGER, INTENT(OUT)                     :: nrot


INTEGER, PARAMETER :: nmax=500
INTEGER :: i,ip,iq,j
REAL(DP) :: c,g,h,s,sm,t,tau,theta,tresh,b(nmax),z(nmax)

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
          t=1D0/(ABS(theta)+SQRT(1D0+theta**2))
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
