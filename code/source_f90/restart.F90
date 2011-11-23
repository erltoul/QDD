#include"define.h"
 
!     **************************

#ifdef REALSWITCH
SUBROUTINE resume(psi,outna)
#else
SUBROUTINE restart2(psi,outna,trealin)
#endif

!     **************************

!     Reads data on wavefunctions, ions, and fields.
!     'resume' is the version for real wavefunctions (static).
!     The variant 'restart2' produces complex wavefunctions,
!     for 'trealin=.false.' from saved complex wavefunctions and
!     for 'trealin=.true.' from real wavefunctions converting them
!     to complex  after reading.
!     All the data are saved in one file called '(r)save' even in the 
!     parallel case.


USE params
USE kinetic
#ifdef REALSWITCH
#if(twostsic)
USE twostr, ONLY: vecsr,ndims
#endif
#endif
IMPLICIT REAL(DP) (A-H,O-Z)


#ifdef REALSWITCH
REAL(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)
LOGICAL,PARAMETER                       :: trealin=.false.
#else
COMPLEX(DP), INTENT(IN OUT)               :: psi(kdfull2,kstate)
LOGICAL, INTENT(IN)                       :: trealin
#endif

CHARACTER (LEN=13), INTENT(IN)            :: outna

REAL(DP), ALLOCATABLE                     :: psiauxr(:)
LOGICAL :: topenf
LOGICAL,PARAMETER :: ttest = .FALSE.

#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)

#ifdef REALSWITCH
REAL(DP),ALLOCATABLE :: psiaux(:)
#else
COMPLEX(DP),ALLOCATABLE :: psiaux(:)
#endif
REAL(DP),DIMENSION(:),ALLOCATABLE :: amoya,epotspa,ekinspa
REAL(DP),DIMENSION(:),ALLOCATABLE :: amoys,epotsps,ekinsps
INTEGER,DIMENSION(:,:),ALLOCATABLE :: nrel2absf


CALL  mpi_comm_rank(mpi_comm_world,myn,icode)

ALLOCATE(psiaux(kdfull2))
ALLOCATE(amoya(ksttot),epotspa(ksttot),ekinspa(ksttot))
ALLOCATE(amoys(kstate),epotsps(kstate),ekinsps(kstate))
ALLOCATE(nrel2absf(kstate,0:knode-1))

DO nod=0,knode-1
  DO i=1,kstate
    nrel2absf(i,nod)=i+nstart_node(nod)
  END DO
END DO

#endif

#if(simpara)
mynact = 0
#else
mynact = myn
#endif





#ifdef REALSWITCH

  IF(mynact==0)  &
    OPEN(UNIT=60,STATUS='old',FORM='unformatted', FILE='rsave.'//outna)  

#else

  IF(mynact==0) THEN
    IF(trealin) THEN
      INQUIRE(60,OPENED=topenf)
      IF(.NOT.topenf) THEN
        OPEN(UNIT=60,STATUS='old',FORM='unformatted',FILE='rsave.'//outna) 
        IF(TTEST) WRITE(*,*) ' rsave opened'
      ELSE
        REWIND(60)
        IF(TTEST) WRITE(*,*) ' unit 60 taken as is'
      END IF
    ELSE
      OPEN(UNIT=60,STATUS='old',FORM='unformatted', FILE='save.'//outna) 
    END IF
  END IF

#endif


IF(mynact==0) THEN
  !  read the iteration where the data has been saved last:
  READ(60) iact,nstate_test,nclust,nion,nspdw
  IF(nstate_test /= nstate_all) &
   STOP ' RESTART: inconsistent nr. of states'
  IF(ttest) WRITE(*,*) ' READ: iact etc at myn=',myn
END IF
#if(parayes)
IF(knode>1) THEN
  CALL mpi_bcast(iact,1,mpi_integer,0,mpi_comm_world,ic)
!  CALL mpi_bcast(nstate,1,mpi_integer,0,mpi_comm_world,ic)
  CALL mpi_bcast(nclust,1,mpi_integer,0,mpi_comm_world,ic)
  CALL mpi_bcast(nion,1,mpi_integer,0,mpi_comm_world,ic)
  CALL mpi_bcast(nspdw,1,mpi_integer,0,mpi_comm_world,ic)
END IF
#endif



#ifdef COMPLEXSWITCH
IF(trealin) THEN
  irest=0
ELSE
  irest=iact
END IF
#endif

!  read wavefunctions:

IF(trealin) ALLOCATE(psiauxr(kdfull2))

#if(parano)
IF(nclust > 0)THEN
  DO nb=1,nstate
    IF(trealin) THEN
      READ(60) occup(nb),psiauxr(1:nxyz)
      psi(1:nxyz,nb) = psiauxr(1:nxyz)
    ELSE
      READ(60) occup(nb),psi(1:nxyz,nb)
    END IF
  END DO
  READ(60) amoy(1:nstate),epotsp(1:nstate),ekinsp(1:nstate)
END IF
#endif
#if(parayes)
CALL mpi_barrier (mpi_comm_world, mpi_ierror)

IF(nclust > 0)THEN
  DO nb=1,nstate_all
    nod = nhome(nb)
    nba = nabs2rel(nb)
    IF(myn == 0) THEN
      IF(trealin) THEN
        READ(60) occupact,psiauxr(1:nxyz)
        psiaux = psiauxr
      ELSE
        READ(60) occupact,psiaux(1:nxyz)
      END IF
        IF(ttest) WRITE(*,*) ' READ: psiaux at myn,nb,occup=',myn,nb,occupact
    END IF
    CALL send_and_receive(occupact,occup(nba),1,0,nod)
#ifdef REALSWITCH
    CALL send_and_receive(psiaux,psi(1,nba),kdfull2,0,nod)
#else
    CALL csend_and_receive(psiaux,psi(1,nba),kdfull2,0,nod)
#endif
    IF(ttest .AND. myn==nod) WRITE(*,*) ' SENT to node,nb=',nod,nb
  END DO

  IF(myn == 0) READ(60) amoya(1:nstate_all),epotspa(1:nstate_all), &
                        ekinspa(1:nstate_all)



  CALL mpi_barrier (mpi_comm_world, mpi_ierror)


  DO nod=0,knode-1
    nstnod = nstate_node(nod)
    DO nb=1,nstnod
      nba = nrel2absf(nb,nod)
      amoys(nb) = amoya(nba)
      epotsps(nb) = epotspa(nba)
      ekinsps(nb) = ekinspa(nba)
    END DO


    CALL send_and_receive(amoys,amoy,nstnod,0,nod)
    CALL send_and_receive(epotsps,epotsp,nstnod,0,nod)
    CALL send_and_receive(ekinsps,ekinsp,nstnod,0,nod)
  END DO
END IF

#endif
IF(trealin) DEALLOCATE(psiauxr)

    


IF(mynact==0) THEN

!  DO i=1,ksttot
!    READ(60) ispin(i),nrel2abs(i),nabs2rel(i),nhome(i)
!  END DO

!  read protonic coordinates and momenta
  IF(nion > 0) THEN
    IF(trealin) THEN
      READ(60) dummy
    ELSE
      READ(60) cx(1:nion),cy(1:nion),cz(1:nion), &
               cpx(1:nion),cpy(1:nion),cpz(1:nion),np(1:nion)
    END IF
  END IF

  IF(ttest) WRITE(*,*) ' ions read in. nion=',nion
!  read substrate coordinates and momenta
#if(raregas)
  IF (isurf /= 0) THEN
    READ(60) imobc(1:nc),xc(1:nc),yc(1:nc),zc(1:nc), &
             pxc(1:nc),pyc(1:nc),pzc(1:nc)
    READ(60) imobe(1:ne),xe(1:ne),ye(1:ne),ze(1:ne), &
             pxe(1:ne),pye(1:ne),pze(1:ne)
    READ(60) imobk(1:nk),xk(1:nk),yk(1:nk),zk(1:nk), &
             pxk(1:nk),pyk(1:nk),pzk(1:nk)
    READ(60) potfixedion(1:kdfull2)
    IF(ttest) WRITE(*,*) ' surface read in. nc,nk,ne=',nc,nk,ne
  END IF
#endif


  IF(nclust > 0) THEN
    READ(60) qe(1:kmom),se(1:3)
    IF(ttest) WRITE(*,*) ' moments read in'
    IF (nabsorb > 0) THEN
      READ(60) rhoabso(1:kdfull2)
      IF(jescmaskorb /=0) THEN
        DO nbe=1,nstate
          READ(60) rhoabsoorb(1:kdfull2,nbe)
        END DO
      END IF
      IF(ttest) WRITE(*,*) ' abso read in'
    END IF
  END IF

END IF

#if(parayes)
IF(knode > 1) THEN
  IF(nclust > 0) THEN
    CALL mpi_bcast(qe(1:kmom),kmom,mpi_double_precision,0,mpi_comm_world,ic)
    CALL mpi_bcast(se(1:3),3,mpi_double_precision,0,mpi_comm_world,ic)
  END IF
  IF (nabsorb > 0) &
    CALL mpi_bcast(rhoabso(1:kdfull2),kdfull2, &
                   mpi_double_precision,0,mpi_comm_world,ic)
END IF
#endif


#ifdef REALSWITCH
!JM
#if(twostsic)
IF(ifsicp >= 6) THEN
  DO i=1,kdim
    DO j=1,kdim
      vecsr(i,j,1) = 1D3
      vecsr(i,j,2) = 2D3
    END DO
  END DO
  ndims(1) = 0
  ndims(2) = 0
  DO iss=1,2
    DO i=1,kdim
      DO j=1,kdim
        READ(60) vecsr(i,j,iss)
      END DO
    END DO
  END DO
  READ(60) ndims(1)
  READ(60) ndims(2)
END IF
#endif
!JM
#endif

IF(trealin) THEN 
  CLOSE(UNIT=60)
ELSE
  CLOSE(UNIT=60,STATUS='keep')
END IF

#if(parayes)
CALL comm_ionconfig()

DEALLOCATE(psiaux)
DEALLOCATE(amoya,epotspa,ekinspa)
DEALLOCATE(amoys,epotsps,ekinsps)
DEALLOCATE(nrel2absf)
#endif


RETURN

#ifdef REALSWITCH
END SUBROUTINE resume
#else
END SUBROUTINE restart2
#endif



!     **************************

#ifdef REALSWITCH
SUBROUTINE RSAVE(psi,isa,outna)
#else
SUBROUTINE SAVE(psi,isa,outna)
#endif

!     **************************

!  writes out the data if mod(iter,isave)=0
!  all the data is saved in the same file called 'save', even in the parallel case

USE params
USE kinetic
#ifdef REALSWITCH
#if(twostsic)
USE twostr, ONLY: vecsr,ndims
#endif
#endif
IMPLICIT REAL(DP) (A-H,O-Z)

#ifdef REALSWITCH
REAL(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)
#else
COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)
#endif

INTEGER, INTENT(IN OUT)                     :: isa
CHARACTER (LEN=13), INTENT(IN OUT)       :: outna
LOGICAL,PARAMETER :: ttest = .FALSE.


#if(parayes)
INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)

#ifdef REALSWITCH
REAL(DP),ALLOCATABLE :: psiaux(:)
#else
COMPLEX(DP),ALLOCATABLE :: psiaux(:)
#endif

REAL(DP) :: occupact
REAL(DP),DIMENSION(:),ALLOCATABLE :: amoya,epotspa,ekinspa
REAL(DP),DIMENSION(:),ALLOCATABLE :: amoys,epotsps,ekinsps
INTEGER,DIMENSION(:,:),ALLOCATABLE :: nrel2absf

CALL  mpi_comm_rank(mpi_comm_world,myn,icode)

ALLOCATE(psiaux(kdfull2))
ALLOCATE(amoya(ksttot),epotspa(ksttot),ekinspa(ksttot))
ALLOCATE(amoys(kstate),epotsps(kstate),ekinsps(kstate))
ALLOCATE(nrel2absf(kstate,0:knode-1))

DO nod=0,knode-1
  DO i=1,kstate
    nrel2absf(i,nod)=i+nstart_node(nod)
  END DO
END DO

CALL mpi_barrier (mpi_comm_world, mpi_ierror)

#endif

#if(simpara)
  mynact = 0
#else
  mynact = myn
#endif

IF(ttest) WRITE(*,*) ' SAVE-BEFORE: myn=',myn
  
#ifdef REALSWITCH

IF(mynact==0) THEN
  IF(isave > 0) THEN
    OPEN(UNIT=60,STATUS='unknown',FORM='unformatted', FILE='rsave.'//outna) 
    WRITE(*,*) ' RSAVE opened'
  ELSE
    OPEN(UNIT=60,STATUS='scratch',FORM='unformatted') 
    WRITE(*,*) ' scratch opened'
  END IF
END IF

#else

IF(mynact==0) &
  OPEN(UNIT=60,STATUS='unknown',FORM='unformatted', FILE='save.'//outna)   

#endif
  
  
!  write iteration at which the data is saved
IF(mynact==0) THEN
  WRITE(60) isa,nstate_all,nclust,nion,nspdw
  IF(TTEST) WRITE(6,*)' SAVE: isa written at myn=',myn
END IF
  
  IF(TTEST) WRITE(6,*)' SAVE: isa,myn=',isa,myn
  
!  write wavefunctions:
#if(parano)
IF(nclust > 0)THEN
  DO nb=1,nstate
    WRITE(60) occup(nb),psi(1:nxyz,nb)
  END DO
  WRITE(60) amoy(1:nstate),epotsp(1:nstate),ekinsp(1:nstate)
END IF
#endif
#if(parayes)

CALL mpi_barrier (mpi_comm_world, mpi_ierror)

IF(nclust > 0)THEN
  DO nb=1,nstate_all
    nod = nhome(nb)
    nba = nabs2rel(nb)
    IF(ttest .AND. myn == nod) WRITE(*,*) ' SENT from node,nb=',nod,nb
    CALL send_and_receive(occup(nba),occupact,1,nod,0)
#ifdef REALSWITCH
    CALL send_and_receive(psi(1,nba),psiaux,kdfull2,nod,0)
#else
    CALL csend_and_receive(psi(1,nba),psiaux,kdfull2,nod,0)
#endif
    IF(myn == 0) THEN
       WRITE(60) occupact,psiaux(1:nxyz)
       IF(ttest) WRITE(*,*) ' RECEIVED and written: node,nb=',nod,nb
    END IF
  END DO
  
  CALL mpi_barrier (mpi_comm_world, mpi_ierror)

  WRITE(*,*) ' before send ekins: myn=',myn
  DO nod=0,knode-1
    nstnod = nstate_node(nod)
    CALL send_and_receive(amoy,amoys,nstnod,nod,0)
    CALL send_and_receive(epotsp,epotsps,nstnod,nod,0)
    CALL send_and_receive(ekinsp,ekinsps,nstnod,nod,0)
    DO nb=1,nstnod
      nba = nrel2absf(nb,nod)
      amoya(nba) = amoys(nb)
      epotspa(nba) = epotsps(nb)
      ekinspa(nba) = ekinsps(nb)
    END DO
  END DO
  IF(myn == 0) WRITE(60) amoya(1:nstate_all),epotspa(1:nstate_all), &
                         ekinspa(1:nstate_all)
  WRITE(*,*) ' after send ekins: myn=',myn

END IF
#endif
    
IF(mynact==0) THEN
!  DO i=1,ksttot
!    WRITE(60) ispin(i),nrel2abs(i),nabs2rel(i),nhome(i)
!  END DO
    
    
    
    
!  write protonic coordinates and momenta:
  IF(nion > 0) &
    WRITE(60) cx(1:nion),cy(1:nion),cz(1:nion), &
             cpx(1:nion),cpy(1:nion),cpz(1:nion),np(1:nion)

!  write substrate coordinates and momenta
#if(raregas)
  IF (isurf /= 0) THEN
    WRITE(60) imobc(1:nc),xc(1:nc),yc(1:nc),zc(1:nc), &
             pxc(1:nc),pyc(1:nc),pzc(1:nc)
    WRITE(60) imobe(1:ne),xe(1:ne),ye(1:ne),ze(1:ne), &
             pxe(1:ne),pye(1:ne),pze(1:ne)
    WRITE(60) imobk(1:nk),xk(1:nk),yk(1:nk),zk(1:nk), &
             pxk(1:nk),pyk(1:nk),pzk(1:nk)
    WRITE(60) potfixedion(1:kdfull2)
  END IF
#endif    
    
!  write dipol moment etc:
    IF(nclust > 0) THEN 
      WRITE(60) qe(1:kmom),se(1:3)
    
      IF (nabsorb > 0) THEN
        WRITE(60) rhoabso(1:kdfull2)
        IF(jescmaskorb /=0) THEN
          DO nbe=1,nstate
            WRITE(60) rhoabsoorb(1:kdfull2,nbe)
          END DO
        END IF
      END IF
    END IF
END IF    
    
#ifdef REALSWITCH
!JM
#if(twostsic)
    IF(ifsicp >= 6) THEN
      DO iss=1,2
        DO i=1,kdim
          DO j=1,kdim
            WRITE(60) vecsr(i,j,iss)
          END DO
        END DO
      END DO
      WRITE(60) ndims(1)
      WRITE(60) ndims(2)
      WRITE(*,*) 'vecsr written'
    END IF
#endif
!JM
#endif

IF(mynact==0 .AND. isave > 0) CLOSE(UNIT=60,STATUS='keep')
    
#if(parayes)
  DEALLOCATE(psiaux)
  DEALLOCATE(amoya,epotspa,ekinspa)
  DEALLOCATE(amoys,epotsps,ekinsps)
  DEALLOCATE(nrel2absf)
#endif
  
  
  IF(jinfo > 0 .AND. MOD(i,jinfo) == 0) THEN
    WRITE(17,'(a,i6,a)') '** data saved at ',isa,' iterations**'
  END IF
  
  RETURN
#ifdef REALSWITCH
END SUBROUTINE RSAVE
#else
END SUBROUTINE SAVE
#endif



#ifdef REALSWITCH

!     **************************

SUBROUTINE restherm()

!     **************************

USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)
OPEN(UNIT=20,STATUS='unknown',FORM='unformatted',FILE='therm')
!     nxyz=nx2*ny2*nz2
READ(20) iact
IF(iact /= irest) THEN
  WRITE(7,*) 'iact=',iact
  STOP 'bad irest in for005dyn !'
END IF
DO ion=1,nion
  READ(20) cx(ion),cy(ion),cz(ion)
  READ(20) cpx(ion),cpy(ion),cpz(ion)
END DO
CLOSE(UNIT=20,STATUS='keep')
RETURN
END SUBROUTINE restherm






!     **************************

SUBROUTINE addcluster(psi,outna)

!     **************************



USE params
USE kinetic
IMPLICIT REAL(DP) (A-H,O-Z)

COMPLEX(DP), INTENT(IN OUT)                  :: psi(kdfull2,kstate)
CHARACTER (LEN=13), INTENT(IN OUT)       :: outna


OPEN(UNIT=60,STATUS='unknown',FORM='unformatted', FILE='save2')


!  read the iteration where the data has been saved last:


READ(60) iact,nstatet,nclustt,niont,nspdwt
DO i=1,nstatet
  READ(60) occup(i+nstate)
END DO


irest=iact


!  read wavefunctions:
IF(nclustt > 0)THEN
  DO nb=1,nstatet
    DO i=1,nxyz
      READ(60) psi(i,nb+nstate)
    END DO
  END DO
END IF

DO i=1,ksttot
  IF (i <= nstatet) THEN
    READ(60) ispin(i+nstate), nrel2abs(i+nstate),nabs2rel(i+nstate)
  ELSE
    READ(60) idum,idum,idum
  END IF
END DO




!  read protonic coordinates and momenta
WRITE(6,*) nion,niont
DO ion=1,niont
  READ(60) cx(ion+nion),cy(ion+nion),cz(ion+nion), np(ion+nion)
  READ(60) cpx(ion+nion),cpy(ion+nion),cpz(ion+nion), np(ion+nion)
END DO


#if(raregas)
IF (isurf /= 0) THEN

  DO i=1,nc
    READ(60) imobc(i)
    READ(60) xc(i),yc(i),zc(i)
    READ(60) pxc(i),pyc(i),pzc(i)
  END DO
  DO i=1,NE
    READ(60) imobe(i)
    READ(60) xe(i),ye(i),ze(i)
    READ(60) pxe(i),pye(i),pze(i)
  END DO
  DO i=1,nk
    READ(60) imobk(i)
    READ(60) xk(i),yk(i),zk(i)
    READ(60) pxk(i),pyk(i),pzk(i)
  END DO
  
  DO i=1,kdfull2
    READ(60) potfixedion(i)
  END DO
  
END IF
#endif


IF(nclustt > 0)THEN
  DO k=1,kmom
    READ(60) qe(k)
  END DO
  DO i=1,3
    READ(60) se(i)
  END DO
END IF


IF (nabsorb > 0) THEN
  DO ii=1,kdfull2
    READ(60) rhoabso(ii)
  END DO
END IF


CLOSE(UNIT=60,STATUS='keep')

nstate=nstate+nstatet
nclust=nclust+nclustt
nion=nion+niont
nspdw=nspdw+nspdwt
nion2=nion


RETURN
END SUBROUTINE addcluster

#endif

#if(parayes)
#ifdef REALSWITCH
SUBROUTINE send_and_receive(instring,outstring,length,in_node,dest_node)

!     Sends 'instring' from odde 'in_node' to
!     'outstring' on node 'out_node'.
!     Both strings are double precsion and have length 'length'

USE params
IMPLICIT REAL(DP) (A-H,O-Z)

INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)

REAL(DP) :: instring(length),outstring(length)
INTEGER :: length,in_node,dest_node


CALL  mpi_comm_rank(mpi_comm_world,myn,icode)

IF(dest_node == in_node .AND. myn == dest_node) THEN
  outstring = instring
  RETURN
END IF

IF(myn == dest_node) &
   CALL mpi_recv(outstring,length,mpi_double_precision,in_node, &
                 mpi_any_tag, mpi_comm_world,is,ic)


IF(myn == in_node)  &
   CALL mpi_send(instring,length,mpi_double_precision,dest_node,1,  &
        mpi_comm_world,ic)


RETURN 

END SUBROUTINE send_and_receive

#else


SUBROUTINE csend_and_receive(instring,outstring,length,in_node,dest_node)

!     Sends 'instring' from odde 'in_node' to
!     'outstring' on node 'out_node'.
!     Both strings are double complex and have length 'length'

USE params
IMPLICIT REAL(DP) (A-H,O-Z)

INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)

COMPLEX(DP) :: instring(length),outstring(length)
INTEGER :: length,in_node,dest_node


CALL  mpi_comm_rank(mpi_comm_world,myn,icode)

IF(dest_node == in_node .AND. myn == dest_node) THEN
  outstring = instring
  RETURN
END IF

IF(myn == dest_node) &
   CALL mpi_recv(outstring,length,mpi_double_complex,in_node, &
                 mpi_any_tag, mpi_comm_world,is,ic)


IF(myn == in_node)  &
   CALL mpi_send(instring,length,mpi_double_complex,dest_node,1,  &
        mpi_comm_world,ic)


RETURN 


END SUBROUTINE csend_and_receive
#endif
#endif
