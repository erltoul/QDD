#include"define.h"
 
!-----init_parallele-------------------------------------------------

SUBROUTINE init_parallele()

!     Initializes parallele computing determines actual 'myn'.
!     The actual node 'myn' is communicated via common.

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
#if(parayes||simpara)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
#endif

!------------------------------------------------------------------

myn=0
#if(parayes||simpara)
WRITE(*,*) ' before mpi_init'
CALL mpi_init(icode)
WRITE(*,*) ' before mpi_comm_size'
CALL  mpi_comm_size(mpi_comm_world,nprocs,icode)
WRITE(*,*) nprocs,knode
#if(parayes)
!IF(nprocs /= knode) STOP
knode = nprocs
IF(knode /= kstate) STOP "knode must be equal to kstate !"
#else
knode = 1
nprocs= 1
#endif
WRITE(*,*) ' before mpi_comm_rank'
CALL  mpi_comm_rank(mpi_comm_world,n,icode)
nn=n                           ! ?
nb=n+1                           ! ?
myn=n
#else             
knode = 1         
nprocs= 1         
#endif

#if(paropenmp)
CALL OMP_SET_NUM_THREADS(numthr)
WRITE(*,*) ' init. OMP:  Nr. threads=',numthr,OMP_GET_NUM_THREADS(),OMP_GET_MAX_THREADS()
nthr = OMP_GET_MAX_THREADS()-1
#else
nthr = 0
#endif
WRITE(*,*) ' INIT_PARALLELE: nthr=',nthr



RETURN
END SUBROUTINE init_parallele


#if(parayes)
!-----comm_inputparams-------------------------------------------------

SUBROUTINE comm_inputparams()

!     Communicates input parameters from node "0" to other nodes.

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)

!------------------------------------------------------------------



CALL  mpi_comm_rank(mpi_comm_world,myn,icode)

CALL i_scatter(kxbox)
CALL i_scatter(kybox)
CALL i_scatter(kzbox)
CALL i_scatter(kstate)


CALL pi_scatter(scaleclustx)
CALL pi_scatter(scaleclusty)
CALL pi_scatter(scaleclustz)
CALL pi_scatter(shiftclustx)
CALL pi_scatter(shiftclusty)
CALL pi_scatter(shiftclustz)
CALL pi_scatter(rotclustx)
CALL pi_scatter(rotclusty)
CALL pi_scatter(rotclustz)
CALL pi_scatter(ehom0)
CALL pi_scatter(ehomx)
CALL pi_scatter(ehomy)
CALL pi_scatter(ehomz)
CALL pi_scatter(shiftwfx)
CALL pi_scatter(shiftwfy)
CALL pi_scatter(shiftwfz)

CALL i_scatter(ispinsep)
CALL i_scatter(idebug)
CALL i_scatter(ishiftcmtoorigin)
CALL i_scatter(iaddcluster)
CALL i_scatter(iswforce)
CALL i_scatter(iplotorbitals)
CALL i_scatter(ievaluate)
CALL i_scatter(ihome)


CALL i_scatter(iemomsrel)
CALL i_scatter(ifspemoms)
CALL i_scatter(iffastpropag)
CALL i_scatter(iflocaliz)
CALL i_scatter(jgeomion)
CALL i_scatter(jangabso)
CALL i_scatter(jgeomel)
CALL i_scatter(jelf)
CALL i_scatter(jstinf)
CALL i_scatter(impactb)
CALL i_scatter(izforcecorr)
CALL i_scatter(iclassicmd)
CALL i_scatter(ntref)
CALL i_scatter(iangabso)
CALL i_scatter(ipes)
CALL i_scatter(nangtheta)
CALL i_scatter(nangphi)
CALL i_scatter(ifreezekspot)
CALL i_scatter(ispherabso)
CALL i_scatter(ifixcmion)
CALL i_scatter(jescmask)
CALL i_scatter(itof)
CALL i_scatter(jescmaskorb)
! CALL i_scatter(icheckformessages)
! CALL i_scatter(jcheckformessages)
CALL i_scatter(nmptheta)
CALL i_scatter(nmpphi)
CALL i_scatter(jmp)
CALL i_scatter(jovlp)
CALL i_scatter(jnorms)
CALL i_scatter(jplotdensitydiff)
CALL i_scatter(jplotdensitydiff2d)
CALL i_scatter(jplotdensity2d)
CALL i_scatter(jcharges)
CALL i_scatter(iscatterelectron)


CALL pi_scatter(drcharges)
CALL pi_scatter(tempion)
CALL pi_scatter(projcharge)
CALL pi_scatter(projvelo)
CALL pi_scatter(projini)
CALL pi_scatter(pimpactb)
CALL pi_scatter(dinmargin)
CALL pi_scatter(delomega)
CALL pi_scatter(angthetal)
CALL pi_scatter(angthetah)
CALL pi_scatter(angphil)
CALL pi_scatter(angphih)
CALL pi_scatter(powabso)
CALL pi_scatter(ekin0pp)
CALL pi_scatter(vxn0)
CALL pi_scatter(vyn0)
CALL pi_scatter(vzn0)
CALL pi_scatter(trequest)
CALL pi_scatter(timefrac)
CALL pi_scatter(scatterelectronenergy)
CALL pi_scatter(scatterelectronvxn)
CALL pi_scatter(scatterelectronvyn)
CALL pi_scatter(scatterelectronvzn)
CALL pi_scatter(scatterelectronx)
CALL pi_scatter(scatterelectrony)
CALL pi_scatter(scatterelectronz)
CALL pi_scatter(scatterelectronw)



IF(myn == 0 .AND. knode /= 1)THEN
  DO nod=1,knode-1
    

    CALL mpi_send(outnam,13,mpi_character,nod,1, mpi_comm_world,ic)
    
    CALL mpi_send(nion2,1,mpi_integer,nod,1,mpi_comm_world,ic)
    

    CALL mpi_send(nabsorb,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(ifsicp,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(ionmdtyp,1,mpi_integer,nod,1, mpi_comm_world,ic)
    CALL mpi_send(ifredmas,1,mpi_integer,nod,1, mpi_comm_world,ic)
    CALL mpi_send(icooltyp,1,mpi_integer,nod,1, mpi_comm_world,ic)
    CALL mpi_send(ipsptyp,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(ipseudo,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(imob,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(numspin,1,mpi_integer,nod,1,mpi_comm_world,ic)
    
    CALL mpi_send(ismax,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(itmax,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(isave,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(istinf,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(jstinf,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(ipasinf,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(dt1,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(irest,1,mpi_integer,nod,1,mpi_comm_world,ic)
    
    CALL mpi_send(centfx,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(centfy,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(centfz,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    
    CALL mpi_send(ispidi,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(iforce,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(iexcit,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(irotat,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(phirot,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    
    CALL mpi_send(i3dx,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(i3dz,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(i3dstate,1,mpi_integer,nod,1, mpi_comm_world,ic)
    CALL mpi_send(istream,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(modrho,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    
    CALL mpi_send(jpos,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(jposcm,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(jvel,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(jener,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(jesc,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(jforce,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(istat,1,mpi_integer,nod,1,mpi_comm_world,ic)
    
    CALL mpi_send(jdip,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(jdiporb,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(jquad,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(jang,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(jspdp,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(jinfo,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(jenergy,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(ivdw,1,mpi_integer,nod,1,mpi_comm_world,ic)
    
    CALL mpi_send(mxforce,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(myforce,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(mzforce,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    
    CALL mpi_send(tempion2,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(idenspl,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(ekmat,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(nfix,1,mpi_integer,nod,1,mpi_comm_world,ic)
    
    CALL mpi_send(itft,1,mpi_integer,nod,1,mpi_comm_world,ic)
    CALL mpi_send(tnode,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(deltat,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(tpeak,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(omega,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(e0,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    
    CALL mpi_send(e1x,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(e1y,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(e1z,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(e2x,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(e2y,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(e2z,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(phi,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    

    CALL mpi_send(nclust,1,mpi_integer,nod,2,mpi_comm_world,ic)
    
    CALL mpi_send(nion,1,mpi_integer,nod,2,mpi_comm_world,ic)
    CALL mpi_send(nc,1,mpi_integer,nod,2,mpi_comm_world,ic)
    CALL mpi_send(NE,1,mpi_integer,nod,2,mpi_comm_world,ic)
    CALL mpi_send(nk,1,mpi_integer,nod,2,mpi_comm_world,ic)
    
    CALL mpi_send(scaleclust,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    
    CALL mpi_send(charge,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(nspdw,1,mpi_integer,nod,2,mpi_comm_world,ic)
    CALL mpi_send(temp,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(occmix,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(isurf,1,mpi_integer,nod,2,mpi_comm_world,ic)
    CALL mpi_send(b2occ,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(gamocc,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(deocc,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(osfac,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(init_lcao,1,mpi_integer,nod,2, mpi_comm_world,ic)
    
    CALL mpi_send(dx,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(dy,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(dz,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    

    CALL mpi_send(radjel,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(surjel,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(bbeta,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(gamma,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(beta4,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(endcon,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(itback,1,mpi_integer,nod,2,mpi_comm_world,ic)
    CALL mpi_send(epswf,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(e0dmp,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(epsoro,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(dpolx,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(dpoly,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(dpolz,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(bcol1,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(bcol2,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(dbcol,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(ntheta,1,mpi_integer,nod,2,mpi_comm_world,ic)
    CALL mpi_send(nphi,1,mpi_integer,nod,2,mpi_comm_world,ic)
    CALL mpi_send(betacol,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(chgcol,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    
    CALL mpi_send(ipotstatic,1,mpi_integer,nod,2, mpi_comm_world,ic)
    CALL mpi_send(surftemp,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(sigmac,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(sigmav,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(sigmak,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(chgc0,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(chge0,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(chgk0,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(cspr,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(mion,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(me,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(mkat,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    
#if(raregas)
    DO i=1,5
      DO j=1,5
        CALL mpi_send(isrtyp(i,j),1,mpi_integer,nod,2, mpi_comm_world,ic)
      END DO
    END DO
#endif
    
    CALL mpi_send(shiftx,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(shifty,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(shiftz,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    
    CALL mpi_send(iaxis,1,mpi_integer,nod,2, mpi_comm_world,ic)
    
    CALL mpi_send(distlayers,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(disttolerance,1,mpi_double_precision,nod,2,  &
        mpi_comm_world,ic)
    
    CALL mpi_send(nunflayc,1,mpi_integer,nod,2, mpi_comm_world,ic)
    CALL mpi_send(nunflaye,1,mpi_integer,nod,2, mpi_comm_world,ic)
    CALL mpi_send(nunflayk,1,mpi_integer,nod,2, mpi_comm_world,ic)
    CALL mpi_send(runfrowc,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(runfrowe,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    CALL mpi_send(runfrowk,1,mpi_double_precision,nod,2, mpi_comm_world,ic)
    
  END DO
  
ELSE IF(myn /= 0 .AND. knode /= 1)THEN
  
  CALL mpi_recv(outnam,13,mpi_character,0,mpi_any_tag, mpi_comm_world,is,ic)
  
  CALL mpi_recv(nion2,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  
  CALL mpi_recv(nabsorb,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(ifsicp,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(ionmdtyp,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(ifredmas,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(icooltyp,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(ipsptyp,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(ipseudo,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(imob,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(numspin,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  
  CALL mpi_recv(ismax,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(itmax,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(isave,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(istinf,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(jstinf,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(ipasinf,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(dt1,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(irest,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  
  CALL mpi_recv(centfx,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(centfy,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(centfz,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  
  CALL mpi_recv(ispidi,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(iforce,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(iexcit,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(irotat,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(phirot,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  
  CALL mpi_recv(i3dx,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(i3dz,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(i3dstate,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(istream,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(modrho,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  
  CALL mpi_recv(jpos,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(jposcm,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(jvel,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(jener,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(jesc,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(jforce,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(istat,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  
  CALL mpi_recv(jdip,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(jdiporb,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(jquad,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(jang,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(jspdp,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(jinfo,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(jenergy,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(ivdw,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  
  CALL mpi_recv(mxforce,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(myforce,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(mzforce,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  
  CALL mpi_recv(tempion2,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(idenspl,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(ekmat,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(nfix,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  
  CALL mpi_recv(itft,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(tnode,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(deltat,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(tpeak,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(omega,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  
  CALL mpi_recv(e0,1,mpi_double_precision,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(e1x,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(e1y,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(e1z,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(e2x,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(e2y,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(e2z,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(phi,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  
  CALL mpi_recv(nclust,1,mpi_integer,0,2,mpi_comm_world,is,ic)
  
  CALL mpi_recv(nion,1,mpi_integer,0,2,mpi_comm_world,is,ic)
  CALL mpi_recv(nc,1,mpi_integer,0,2,mpi_comm_world,is,ic)
  CALL mpi_recv(NE,1,mpi_integer,0,2,mpi_comm_world,is,ic)
  CALL mpi_recv(nk,1,mpi_integer,0,2,mpi_comm_world,is,ic)
  
  CALL mpi_recv(scaleclust,1,mpi_double_precision,0,2, mpi_comm_world,is,ic)
  
  CALL mpi_recv(charge,1,mpi_double_precision,0,2, mpi_comm_world,is,ic)
  CALL mpi_recv(nspdw,1,mpi_integer,0,2,mpi_comm_world,is,ic)
  CALL mpi_recv(temp,1,mpi_double_precision,0,2, mpi_comm_world,is,ic)
  CALL mpi_recv(occmix,1,mpi_double_precision,0,2, mpi_comm_world,is,ic)
  CALL mpi_recv(isurf,1,mpi_integer,0,2,mpi_comm_world,is,ic)
  CALL mpi_recv(b2occ,1,mpi_double_precision,0,2, mpi_comm_world,is,ic)
  CALL mpi_recv(gamocc,1,mpi_double_precision,0,2, mpi_comm_world,is,ic)
  CALL mpi_recv(deocc,1,mpi_double_precision,0,2, mpi_comm_world,is,ic)
  CALL mpi_recv(osfac,1,mpi_double_precision,0,2, mpi_comm_world,is,ic)
  CALL mpi_recv(init_lcao,1,mpi_integer,0,2, mpi_comm_world,is,ic)
  
  CALL mpi_recv(dx,1,mpi_double_precision,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(dy,1,mpi_double_precision,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(dz,1,mpi_double_precision,0,mpi_any_tag, mpi_comm_world,is,ic)
  
  CALL mpi_recv(radjel,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(surjel,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(bbeta,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(gamma,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(beta4,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(endcon,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(itback,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(epswf,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(e0dmp,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(epsoro,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(dpolx,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(dpoly,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(dpolz,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(bcol1,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(bcol2,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(dbcol,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(ntheta,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(nphi,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(betacol,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(chgcol,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  
  CALL mpi_recv(ipotstatic,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(surftemp,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(sigmac,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(sigmav,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(sigmak,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(chgc0,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(chge0,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(chgk0,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(cspr,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(mion,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(me,1,mpi_double_precision,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(mkat,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  
#if(raregas)
  DO i=1,5
    DO j=1,5
      CALL mpi_recv(isrtyp(i,j),1,mpi_integer,0,mpi_any_tag,  &
          mpi_comm_world,is,ic)
    END DO
  END DO
#endif
  
  CALL mpi_recv(shiftx,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(shifty,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(shiftz,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  
  CALL mpi_recv(iaxis,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  
  CALL mpi_recv(distlayers,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(disttolerance,1,mpi_double_precision,0,  &
      mpi_any_tag,mpi_comm_world,is,ic)
  
  CALL mpi_recv(nunflayc,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(nunflaye,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(nunflayk,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(runfrowc,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(runfrowe,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(runfrowk,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  
END IF

RETURN
END SUBROUTINE comm_inputparams

#endif


#if(parayes)
!-----comm_ionconfig-------------------------------------------------

SUBROUTINE comm_ionconfig()

!     Communicates ion configuration from node "0" to other nodes.

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)

!------------------------------------------------------------------

n = myn
IF(n == 0 .AND. knode /= 1)THEN
  
  DO nod=1,knode-1
    CALL mpi_send(cx(1),nion,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(cy(1),nion,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(cz(1),nion,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(np(1),nion,mpi_integer,nod,1, mpi_comm_world,ic)
    CALL mpi_send(cpx(1),nion,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(cpy(1),nion,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(cpz(1),nion,mpi_double_precision,nod,1, mpi_comm_world,ic)
    IF(init_lcao == 1) THEN
      CALL mpi_send(radini(1),nion,mpi_double_precision,  &
          nod,1,mpi_comm_world,ic)
      DO ion=1,nion
        CALL mpi_send(initord(1,nion),3,mpi_integer,nod,1, mpi_comm_world,ic)
      END DO
    END IF
    CALL mpi_send(tblock,ng+1,mpi_logical,nod,1, mpi_comm_world,ic)
    IF(ipsptyp == 1) THEN
      CALL mpi_send(tnonloc,ng+1,mpi_logical,nod,1, mpi_comm_world,ic)
      CALL mpi_send(tnonlocany,1,mpi_logical,nod,1, mpi_comm_world,ic)
    END IF
  END DO
  
ELSE IF(n /= 0 .AND. knode /= 1)THEN
  
  CALL mpi_recv(cx(1),nion,mpi_double_precision,0,  &
      mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(cy(1),nion,mpi_double_precision,0,  &
      mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(cz(1),nion,mpi_double_precision,0,  &
      mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(np(1),nion,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(cpx(1),nion,mpi_double_precision,0,  &
      mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(cpy(1),nion,mpi_double_precision,0,  &
      mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(cpz(1),nion,mpi_double_precision,0,  &
      mpi_any_tag,mpi_comm_world,is,ic)
  IF(init_lcao == 1) THEN
    CALL mpi_recv(radini(1),nion,mpi_double_precision,  &
        0,mpi_any_tag,mpi_comm_world,is,ic)
    DO ion=1,nion
      CALL mpi_recv(initord(1,nion),3,mpi_integer,0,  &
          mpi_any_tag,mpi_comm_world,is,ic)
    END DO
  END IF
  CALL mpi_recv(tblock,ng+1,mpi_logical,0,  &
        mpi_any_tag,mpi_comm_world,is,ic)
  IF(ipsptyp == 1) THEN
    CALL mpi_recv(tnonloc,ng+1,mpi_logical,0,mpi_any_tag,mpi_comm_world,is,ic)
    CALL mpi_recv(tnonlocany,1,mpi_logical,0,mpi_any_tag,mpi_comm_world,is,ic)
  END IF
  
END IF

RETURN
END SUBROUTINE comm_ionconfig
!-----comm_periodic-------------------------------------------------

SUBROUTINE comm_periodic()

!     Communicates periodic table from node "0" to other nodes.

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)

!------------------------------------------------------------------

n = myn
IF(n == 0 .AND. knode /= 1)THEN
  
  DO nod=1,knode-1
    CALL mpi_send(ch,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(amu,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(cc1,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(cc2,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(crloc,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(crs,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(chs,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(chg1,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(chg2,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(sgm1,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(sgm2,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(r0g,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(r1g,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(r2g,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(h0_11g,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(h0_22g,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(h0_33g,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(h1_11g,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(h1_22g,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(h2_11g,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(h0_12g,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(radiong,185,mpi_double_precision,nod,1, mpi_comm_world,ic)
  END DO
  
ELSE IF(n /= 0 .AND. knode /= 1)THEN
  
  CALL mpi_recv(ch,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(amu,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(cc1,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(cc2,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(crloc,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(crs,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(chs,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(chg1,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(chg2,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(sgm1,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(sgm2,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(r0g,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(r1g,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(r2g,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(h0_11g,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(h0_22g,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(h0_33g,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(h1_11g,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(h1_22g,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(h2_11g,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(h0_12g,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(radiong,185,mpi_double_precision,0,mpi_any_tag,mpi_comm_world,is,ic)
  
END IF

RETURN
END SUBROUTINE comm_periodic
#endif



#if(parayes)
!-----comm_siman-------------------------------------------------

SUBROUTINE comm_simann()

!     Communicates parameters for simulated annealing  (??)

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)

!------------------------------------------------------------------

n = myn

IF(n == 0 .AND. knode /= 1)THEN
  
  DO nod=1,knode-1
    CALL mpi_send(ionsin,1,mpi_integer,nod,1, mpi_comm_world,ic)
    CALL mpi_send(iknow,1,mpi_integer,nod,1, mpi_comm_world,ic)
    CALL mpi_send(facann,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(nrun,1,mpi_integer,nod,1, mpi_comm_world,ic)
    CALL mpi_send(nloop1,1,mpi_integer,nod,1, mpi_comm_world,ic)
    CALL mpi_send(nloop2,1,mpi_integer,nod,1, mpi_comm_world,ic)
    CALL mpi_send(cptemp,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(delpos,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(ERR,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(errks0,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(trfac2,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(prfac2,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(errsim,1,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(ncsim,1,mpi_integer,nod,1, mpi_comm_world,ic)
    CALL mpi_send(errtot,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(ncon,1,mpi_integer,nod,1, mpi_comm_world,ic)
    CALL mpi_send(erfac1,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(trfac1,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(prfac1,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
    CALL mpi_send(rand_seed,isize_seed,mpi_integer,nod,1, mpi_comm_world,ic)
  END DO
  
ELSE IF(n /= 0 .AND. knode /= 1)THEN
  
  CALL mpi_recv(ionsin,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(iknow,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(facann,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(nrun,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(nloop1,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(nloop2,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(cptemp,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(delpos,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(ERR,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(errks0,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(trfac2,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(prfac2,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(errsim,1,1,mpi_double_precision,0,  &
      mpi_any_tag,mpi_comm_world,is,ic)
  CALL mpi_recv(ncsim,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(errtot,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(ncon,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
  CALL mpi_recv(erfac1,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(trfac1,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(prfac1,1,mpi_double_precision,0,mpi_any_tag,  &
      mpi_comm_world,is,ic)
  CALL mpi_recv(rand_seed,isize_seed,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
END IF

RETURN
END SUBROUTINE comm_simann
#endif

#if(parayes)


!       *****************************

SUBROUTINE pi_allreduce(rin,rout,kdf,id,id1,id2,ic)

!       *****************************

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
INCLUDE 'mpif.h'


REAL, INTENT(IN OUT)                     :: rin(kdf)
REAL, INTENT(IN OUT)                     :: rout(kdf)
INTEGER, INTENT(IN OUT)                  :: kdf
INTEGER, INTENT(IN OUT)                  :: id
INTEGER, INTENT(IN OUT)                  :: id1
INTEGER, INTENT(IN OUT)                  :: id2
INTEGER, INTENT(IN OUT)                  :: ic

CALL mpi_allreduce(rin,rout,kdf,mpi_double_precision,  &
    mpi_sum,mpi_comm_world,ic)
RETURN
END SUBROUTINE pi_allreduce

!       *****************************

SUBROUTINE pi_allgather(rin,mint,nt,rout,kdf)

!       *****************************

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
INCLUDE 'mpif.h'


REAL, INTENT(IN OUT)                     :: rin(kdf)
INTEGER, INTENT(IN OUT)                  :: mint
INTEGER, INTENT(IN OUT)                  :: nt
REAL, INTENT(IN OUT)                     :: rout(kdf)
INTEGER, INTENT(IN OUT)                  :: kdf

CALL mpi_allgather(rin(mint),nt,mpi_double_precision,  &
    rout,nt,mpi_double_precision,mpi_comm_world,ic)
RETURN
END SUBROUTINE pi_allgather

!       *****************************

SUBROUTINE i_allreduce(rin,rout,kdf,id,id1,id2,ic)

!       *****************************

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
INCLUDE 'mpif.h'


REAL, INTENT(OUT)                        :: rin(kdf)
REAL, INTENT(IN OUT)                     :: rout(kdf)
INTEGER, INTENT(IN OUT)                  :: kdf
INTEGER, INTENT(IN OUT)                  :: id
INTEGER, INTENT(IN OUT)                  :: id1
INTEGER, INTENT(IN OUT)                  :: id2
INTEGER, INTENT(IN OUT)                  :: ic

INTEGER :: is(mpi_status_size)

CALL  mpi_comm_size(mpi_comm_world,nprocs,icode)
IF(nprocs /= kstate) STOP
CALL  mpi_comm_rank(mpi_comm_world,n,icode)
l=LOG(kstate*1D0)/LOG(2D0)-1
DO it=0,l
  j=(n/(2D0**FLOAT(it)))
  IF(MOD(j,2) == 0) THEN
    nt=n+2D0**FLOAT(it)
  ELSE
    nt=n-2D0**FLOAT(it)
  END IF
  IF(MOD(j,2) == 0) THEN
    CALL mpi_send(rin,mx,mpi_double_precision, nt,0,mpi_comm_world,icode)
    CALL mpi_recv(rout,mx,mpi_double_precision, nt,0,mpi_comm_world,is,icode)
  ELSE
    CALL mpi_recv(rout,mx,mpi_double_precision, nt,0,mpi_comm_world,is,icode)
    CALL mpi_send(rin,mx,mpi_double_precision, nt,0,mpi_comm_world,icode)
  END IF
  DO i=1,mx
!dir$     unroll 8
    rin(i)=rin(i)+rout(i)
  END DO
END DO
DO i=1,mx
  rout(i)=rin(i)
END DO
RETURN
END SUBROUTINE i_allreduce

!       *****************************

SUBROUTINE i_allgather(rin,mint,ntr,rout,kdf)

!       *****************************

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
INCLUDE 'mpif.h'

REAL, INTENT(IN)                         :: rin(kdf)
INTEGER, INTENT(IN)                      :: mint
INTEGER, INTENT(IN)                      :: ntr
REAL, INTENT(OUT)                        :: rout(kdf)
INTEGER, INTENT(IN OUT)                  :: kdf
INTEGER :: is(mpi_status_size)



CALL  mpi_comm_rank(mpi_comm_world,n,icode)
nn=n+1
l=LOG(kstate*1D0)/LOG(2D0)-1
nax=ntr
nin=mint
DO i=nin,nin+nax
!     dir$     unroll 8
  rout(i)=rin(i)
END DO
DO it=0,l
  j=(n/(2D0**FLOAT(it)))
  IF(MOD(j,2) == 0) THEN
    nt=n+2D0**FLOAT(it)
  ELSE
    nt=n-2D0**FLOAT(it)
  END IF
  IF(MOD(j,2) == 0) THEN
    CALL mpi_send(nax,1,mpi_integer,nt,1,mpi_comm_world,icode)
    CALL mpi_send(nin,1,mpi_integer,nt,2,mpi_comm_world,icode)
    CALL mpi_send(rout(nin),nax,mpi_double_precision,  &
        nt,3,mpi_comm_world,icode)
    CALL mpi_recv(maxr,1,mpi_integer, nt,4,mpi_comm_world,is,icode)
    CALL mpi_recv(minr,1,mpi_integer, nt,5,mpi_comm_world,is,icode)
    CALL mpi_recv(rout(minr),maxr,mpi_double_precision,  &
        nt,6,mpi_comm_world,is,icode)
  ELSE
    CALL mpi_recv(maxr,1,mpi_integer, nt,1,mpi_comm_world,is,icode)
    CALL mpi_recv(minr,1,mpi_integer, nt,2,mpi_comm_world,is,icode)
    CALL mpi_recv(rout(minr),maxr,mpi_double_precision,  &
        nt,3,mpi_comm_world,is,icode)
    CALL mpi_send(nax,1,mpi_integer,nt,4,mpi_comm_world,icode)
    CALL mpi_send(nin,1,mpi_integer,nt,5,mpi_comm_world,icode)
    CALL mpi_send(rout(nin),nax,mpi_double_precision,  &
        nt,6,mpi_comm_world,icode)
    nin=nin-maxr
  END IF
  nax=nax+maxr
  
END DO

RETURN
END SUBROUTINE i_allgather
!       *****************************


SUBROUTINE pi_scatter(buf)

!       *****************************

USE params
!IMPLICIT REAL(DP) (A-H,O-Z)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)


REAL,INTENT(IN OUT) :: buf

IF(myn == 0 .AND. knode /= 1)THEN
  DO nod=1,knode-1
    CALL mpi_send(buf,1,mpi_double_precision,nod,1, mpi_comm_world,ic)
  END DO
ELSE IF(myn /= 0 .AND. knode /= 1)THEN
  CALL mpi_recv(buf,1,mpi_double_precision,0,mpi_any_tag, mpi_comm_world,is,ic)
ELSE
  STOP ' PI_SCATTER1 only for KNODE > 1'
END IF

END SUBROUTINE pi_scatter


SUBROUTINE i_scatter(buf)

!       *****************************

USE params
!IMPLICIT REAL(DP) (A-H,O-Z)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)


INTEGER,INTENT(IN OUT) :: buf

IF(myn == 0 .AND. knode /= 1)THEN
  DO nod=1,knode-1
    CALL mpi_send(buf,1,mpi_integer,nod,1, mpi_comm_world,ic)
  END DO
ELSE IF(myn /= 0 .AND. knode /= 1)THEN
  CALL mpi_recv(buf,1,mpi_integer,0,mpi_any_tag, mpi_comm_world,is,ic)
ELSE
  STOP ' PI_SCATTER1 only for KNODE > 1'
END IF

END SUBROUTINE i_scatter


#endif




#if(parayes)
!-----prispe_parallele-------------------------------------------------

SUBROUTINE prispe_parallele(iunit,it)

!     Collects single particle properties and writes it on
!     unit 'iunit'

USE params
!USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN) :: iunit

INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
REAL(DP) :: prisav(kstate,5)
REAL(DP) :: prisav_all(ksttot,5)
INTEGER :: iprisav(kstate,2)

LOGICAL,PARAMETER :: ttest=.FALSE.

!------------------------------------------------------------------

!                                   collect into aux fields
DO nb=1,nstate
  prisav(nb,1) = occup(nb)
  prisav(nb,2) = ekinsp(nb)
  prisav(nb,3) = amoy(nb)
  prisav(nb,4) = evarsp(nb)
  prisav(nb,5) = spnorm(nb)
  iprisav(nb,1)  = nrel2abs(nb)
  iprisav(nb,2)  = 3-2*ispin(nrel2abs(nb))
END DO

!                                   communicate and print
prisav_all = 0D0
IF(myn /= 0) THEN
  nod = myn
  CALL mpi_send(prisav,5*kstate,mpi_double_precision, 0,nod,mpi_comm_world,ic)
  CALL mpi_send(iprisav,2*kstate,mpi_integer, 0,nod,mpi_comm_world,ic)
  IF(ttest) WRITE(*,*) ' INFO: sent at node:',myn
  DO n=1,nstate
    nact = nrel2abs_other(n,0) 
    prisav_all(nact,:) = prisav(n,:)
  ENDDO
ELSE
  DO nod2=0,knode-1
    IF(nod2 > 0) THEN
      CALL mpi_recv(prisav,5*kstate,mpi_double_precision,  &
          nod2,mpi_any_tag,mpi_comm_world,is,ic)
      CALL mpi_recv(iprisav,2*kstate,mpi_integer,  &
          nod2,mpi_any_tag,mpi_comm_world,is,ic)
      IF(ttest) WRITE(*,*) ' INFO: recv from  node=',nod2
    END IF
    DO nb=1,nstate_node(nod2)
      nact = nrel2abs_other(nb,nod2) 
      prisav_all(nact,:) = prisav(nb,:)
      IF(numspin==2) THEN
        WRITE(iunit,'(a,i3,a,i3,3f9.5,1pg12.4)') 'level:',iprisav(nb,1),  &
          '  spin,occup,ekin,esp,variance =', iprisav(nb,2),(prisav(nb,k),k=1,4)
      ELSE
        WRITE(iunit,'(a,i3,a,3f9.5,1pg12.4)')  &
          'level:',nrel2abs(nb), &
          '  occup,ekin,esp,variance=', (prisav(nb,k),k=1,4)
      END IF
    END DO
  END DO

  IF(jstinf > 0 .AND. MOD(it,jstinf)==0) then
    CALL safeopen(77,it,jstinf,'pspenergies')
    WRITE(77,'(1f15.6,500f12.6)') tfs,(prisav_all(nb,3),nb=1,nstate_all)
    CALL flush(77)
    CALL safeopen(76,it,jstinf,'pspvariances')
    WRITE(76,'(1f15.6,500f12.6)') tfs,(prisav_all(nb,4),nb=1,nstate_all)
    CALL flush(76)
  END IF

  IF(jesc > 0 .AND. jnorms>0 .AND. MOD(it,jnorms) == 0) THEN
    WRITE(806,'(1f15.6,500f12.8)') tfs,1D0-prisav_all(1:nstate_all,5)
  END IF


END IF

!WRITE(*,*) 'prispe_parallele done'

RETURN
END SUBROUTINE prispe_parallele

#endif

