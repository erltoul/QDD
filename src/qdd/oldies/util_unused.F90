!-----projmoms-------------------------------------------------------projmoms
! REAL version
!-----------------------------------------------------------------------
SUBROUTINE r_projmoms(rho,psi)

! Multipole moments relative to c.m. and relative to ion coordinate 'cz'.
! ??? presently not used ???

USE params
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



!     **************************

REAL(DP) FUNCTION wfnorm_r(psi1,psi2)

!  Overlap of wavefunctions 'psi1' and 'psi2'.     REAL version

USE params
IMPLICIT NONE
COMPLEX(DP), INTENT(IN OUT)                  :: psi1(kdfull2)
COMPLEX(DP), INTENT(IN OUT)                  :: psi2(kdfull2)

INTEGER :: ii
REAL(DP) :: acc

acc = 0D0
DO ii=1,nxyz
  acc = REAL(psi1(ii),DP)*REAL(psi2(ii),DP) + AIMAG(psi1(ii))*AIMAG(psi2(ii)) + acc
END DO
wfnorm_r = acc*dvol
RETURN
END FUNCTION wfnorm_r
!     **************************

REAL(DP) FUNCTION wfnorm_i(psi1,psi2)

!  Overlap of wavefunctions 'psi1' and 'psi2'.     COMPLEX version

USE params
IMPLICIT NONE
COMPLEX(DP), INTENT(IN OUT)                  :: psi1(kdfull2)
COMPLEX(DP), INTENT(IN OUT)                  :: psi2(kdfull2)

INTEGER :: ii 
REAL(DP) :: acc
acc = 0D0
DO ii=1,nxyz
  acc = -REAL(psi1(ii),DP)*AIMAG(psi2(ii)) +AIMAG(psi1(ii))*REAL(psi2(ii),DP) + acc
END DO
wfnorm_i = acc*dvol
RETURN
END FUNCTION wfnorm_i





!**********************   from restart.F90  *******************

!     **************************

SUBROUTINE restherm()

!     **************************

USE params
USE kinetic
IMPLICIT NONE
INTEGER :: iact,ion

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
IMPLICIT NONE

COMPLEX(DP), INTENT(IN OUT)              :: psi(kdfull2,kstate)
CHARACTER (LEN=13), INTENT(IN OUT)       :: outna

INTEGER :: i,iact,idum,ii,ion,k,nb,nclustt,niont,nspdwt,nstatet
OPEN(UNIT=ifile,STATUS='unknown',FORM='unformatted', FILE='save.'//outna)


!  read the iteration where the data has been saved last:


READ(ifile) iact,nstatet,nclustt,niont,nspdwt
DO i=1,nstatet
  READ(ifile) occup(i+nstate)
END DO


irest=iact


!  read wavefunctions:
IF(nclustt > 0)THEN
  DO nb=1,nstatet
    DO i=1,nxyz
      READ(ifile) psi(i,nb+nstate)
    END DO
  END DO
END IF

DO i=1,ksttot
  IF (i <= nstatet) THEN
    READ(ifile) ispin(i+nstate), nrel2abs(i+nstate),nabs2rel(i+nstate)
  ELSE
    READ(ifile) idum,idum,idum
  END IF
END DO




!  read protonic coordinates and momenta
WRITE(6,*) nion,niont
DO ion=1,niont
  READ(ifile) cx(ion+nion),cy(ion+nion),cz(ion+nion), np(ion+nion)
  READ(ifile) cpx(ion+nion),cpy(ion+nion),cpz(ion+nion), np(ion+nion)
END DO


#if(raregas)
IF (isurf /= 0) THEN

  DO i=1,nc
    READ(ifile) imobc(i)
    READ(ifile) xc(i),yc(i),zc(i)
    READ(ifile) pxc(i),pyc(i),pzc(i)
  END DO
  DO i=1,NE
    READ(ifile) imobe(i)
    READ(ifile) xe(i),ye(i),ze(i)
    READ(ifile) pxe(i),pye(i),pze(i)
  END DO
  DO i=1,nk
    READ(ifile) imobk(i)
    READ(ifile) xk(i),yk(i),zk(i)
    READ(ifile) pxk(i),pyk(i),pzk(i)
  END DO
  
  DO i=1,kdfull2
    READ(ifile) potfixedion(i)
  END DO
  
END IF
#endif


IF(nclustt > 0)THEN
  DO k=1,kmom
    READ(ifile) qe(k)
  END DO
  DO i=1,3
    READ(ifile) se(i)
  END DO
END IF


IF (nabsorb > 0) THEN
  DO ii=1,kdfull2
    READ(ifile) rhoabso(ii)
  END DO
END IF


CLOSE(UNIT=ifile,STATUS='keep')

nstate=nstate+nstatet
nclust=nclust+nclustt
nion=nion+niont
nspdw=nspdw+nspdwt
nion2=nion


RETURN
END SUBROUTINE addcluster


