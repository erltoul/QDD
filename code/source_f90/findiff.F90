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

MODULE kinetic
USE params
IMPLICIT NONE

CONTAINS
#if(findiff)
!-----rkin3D_3r------------------------------------------------------------

SUBROUTINE rkin3d(psi,dxpsi)


!      include 'pot3D.inc'

!     computes kinetic energy (- Laplacian) 3 point finite differences
!     and reflecting boundary conditions



REAL(DP), INTENT(IN OUT)                     :: psi(minx:maxx,miny:maxy,minz:maxz)
REAL(DP), INTENT(OUT)                        :: dxpsi(minx:maxx,miny:maxy,minz:maxz)



!-----------------------------------------------------------------------


INTEGER :: ix, iy, iz              ! loop indices
INTEGER :: lengx,lengy,lengz


!-----------------------------------------------------------------------


!     number of mesh points in each direction

lengx  = maxx+maxx+1
lengy  = maxy+maxy+1
lengz  = maxz+maxz+1


!     reset accumulator


DO iz = minz, maxz
  DO iy = miny, maxy
    DO ix = minx, maxx
      dxpsi(ix,iy,iz) = zero
    END DO
  END DO
END DO


!     kinetic energy in x direction


DO iz = minz, maxz
  DO iy = miny, maxy
    CALL rkin1d_3r(psi(minx,iy,iz),dx,lengx,1, dxpsi(minx,iy,iz))
  END DO
END DO



!     kinetic energy in y direction



DO iz = minz, maxz
  DO ix = minx, maxx
    CALL rkin1d_3r(psi(ix,miny,iz),dy,lengy,lengx, dxpsi(ix,miny,iz))
  END DO
END DO


!     kinetic energy in z direction


DO iy = miny, maxy
  DO ix = minx, maxx
    CALL rkin1d_3r(psi(ix,iy,minz),dz,lengz,lengy*lengx, dxpsi(ix,iy,minz))
  END DO
END DO


RETURN
END SUBROUTINE rkin3d
!-----d2_3r1D ----------------------------------------------------------

SUBROUTINE rkin1d_3r(psi,deltax,nmax,inc,dxpsi)

!USE params

!     computes kinetic energy in one dimension and accumulates
!     that on 'dxpsi'
!     uses five point finite differences.

!     psi    = input array ( a 3D array in calling routine)
!     dxpsi  = output containing 2. derivative + input value
!     deltax = mesh spacing
!     nmax   = number of mesh point in given direction
!     inc    = increment which connects neighbouring mesh points


REAL(DP), INTENT(IN)                         :: psi(*)
REAL(DP), INTENT(IN OUT)                     :: deltax
INTEGER, INTENT(IN)                      :: nmax
INTEGER, INTENT(IN)                      :: inc
REAL(DP), INTENT(OUT)                        :: dxpsi(*)

INTEGER :: ninc
INTEGER :: i
REAL(DP):: d2i
!REAL(DP):: ! Wellenfunktion psi
!REAL(DP):: ! Ableitung der Wellenfunktion psi

!-----------------------------------------------------------------------

d2i = -one/(deltax*deltax)
dxpsi(1) = d2i*(psi(inc+1)-2*psi(1)) + dxpsi(1)

ninc = 1
DO i = 2,nmax-1
  ninc=ninc+inc
  dxpsi(ninc) = dxpsi(ninc) +d2i*(psi(ninc+inc)+psi(ninc-inc)-2.0D0*psi(ninc))
END DO

ninc=ninc + inc
dxpsi(ninc) = d2i*(psi(ninc-inc)-2*psi(ninc)) + dxpsi(ninc)

RETURN
END SUBROUTINE rkin1d_3r

!-----ckin3D_3r------------------------------------------------------------

SUBROUTINE ckin3d(psi,dxpsi)

!USE params
!      include 'pot3D.inc'

!     computes kinetic energy (- Laplacian) 3 point finite differences
!     and reflecting boundary conditions



COMPLEX(DP), INTENT(IN OUT)                  :: psi(minx:maxx,miny:maxy,minz:maxz)
COMPLEX(DP), INTENT(OUT)                     :: dxpsi(minx:maxx,miny:maxy,minz:maxz)



!-----------------------------------------------------------------------


INTEGER :: ix, iy, iz              ! loop indices
INTEGER :: lengx,lengy,lengz


!-----------------------------------------------------------------------


!     number of mesh points in each direction

lengx  = maxx+maxx+1
lengy  = maxy+maxy+1
lengz  = maxz+maxz+1


!     reset accumulator


DO iz = minz, maxz
  DO iy = miny, maxy
    DO ix = minx, maxx
      dxpsi(ix,iy,iz) = zero
    END DO
  END DO
END DO


!     kinetic energy in x direction


DO iz = minz, maxz
  DO iy = miny, maxy
    CALL ckin1d_3r(psi(minx,iy,iz),dx,lengx,1, dxpsi(minx,iy,iz))
  END DO
END DO



!     kinetic energy in y direction



DO iz = minz, maxz
  DO ix = minx, maxx
    CALL ckin1d_3r(psi(ix,miny,iz),dy,lengy,lengx, dxpsi(ix,miny,iz))
  END DO
END DO


!     kinetic energy in z direction


DO iy = miny, maxy
  DO ix = minx, maxx
    CALL ckin1d_3r(psi(ix,iy,minz),dz,lengz,lengy*lengx, dxpsi(ix,iy,minz))
  END DO
END DO


RETURN
END SUBROUTINE ckin3d
!-----d2_3r1D ----------------------------------------------------------

SUBROUTINE ckin1d_3r(psi,deltax,nmax,inc,dxpsi)

!USE params

!     computes kinetic energy in one dimension and accumulates
!     that on 'dxpsi'
!     uses five point finite differences.

!     psi    = input array ( a 3D array in calling routine)
!     dxpsi  = output containing 2. derivative + input value
!     deltax = mesh spacing
!     nmax   = number of mesh point in given direction
!     inc    = increment which connects neighbouring mesh points


COMPLEX(DP), INTENT(IN)                      :: psi(*)
REAL(DP), INTENT(IN OUT)                     :: deltax
INTEGER, INTENT(IN)                      :: nmax
INTEGER, INTENT(IN)                      :: inc
COMPLEX(DP), INTENT(OUT)                     :: dxpsi(*)

INTEGER :: ninc
INTEGER :: i
REAL(DP):: d2i
!COMPLEX(DP) :: ! Wellenfunktion psi
!COMPLEX(DP) :: ! Ableitung der Wellenfunktion psi

!-----------------------------------------------------------------------

d2i = -one/(deltax*deltax)
dxpsi(1) = d2i*(psi(inc+1)-2*psi(1)) + dxpsi(1)

ninc = 1
DO i = 2,nmax-1
  ninc=ninc+inc
  dxpsi(ninc) = dxpsi(ninc) +d2i*(psi(ninc+inc)+psi(ninc-inc)-2.0D0*psi(ninc))
END DO

ninc=ninc + inc
dxpsi(ninc) = d2i*(psi(ninc-inc)-2*psi(ninc)) + dxpsi(ninc)

RETURN
END SUBROUTINE ckin1d_3r

!-----d3mixpropag_3----------------------------------------------------

SUBROUTINE d3mixpropag (psi, deltim)

!USE params
!      include 'pot3D.inc'

!     Propagates T Operator for mixed step:

!              1 - 0.5 * i * deltim * T^
!     neupsi = ------------------------  * psi
!              1 + 0.5 * i * deltim * T^

!     in 3 Dimensions with 3- point precision
!     with boundary conditions zero at doth ends (reflecting)

!   . The "Laplacian"
!     operator is here dimensionless with diagonal element -2
!     and off-diagonals 1.


!     list parameters:


COMPLEX(DP), INTENT(IN OUT)                  :: psi(minx:maxx, miny:maxy,minz:maxz)
REAL(DP), INTENT(IN OUT)                     :: deltim

!                      ! old Function
!      complex neupsi (MINX:MAXX, MINY:MAXY, MINZ:MAXZ)
!                      ! new Function

! Timestep
!     variables:

INTEGER :: ix, iy, iz
! loop indices

!-----------------------------------------------------------------------

!      write(6,*) ' enter D3MIXSTEP: norm=',wfnorm(psi)
DO iz = minz, maxz
  DO iy = miny, maxy
    CALL kinprop_1d3(psi(minx,iy,iz),2*maxx+1,1,dx,deltim )
  END DO
END DO            ! x- direction
!      write(6,*) ' after x D3MIXSTEP: norm=',wfnorm(psi)

DO iz = minz, maxz
  DO ix = minx, maxx
    CALL kinprop_1d3(psi(ix,miny,iz),2*maxy+1,2*maxx+1, dy,deltim)
  END DO
END DO            ! y- direction
!      write(6,*) ' after y D3MIXSTEP: norm=',wfnorm(psi)

DO iy = miny, maxy
  DO ix = minx, maxx
    CALL kinprop_1d3(psi(ix,iy,minz),2*maxz+1,  &
        (2*maxx+1)*(2*maxy+1),dz,deltim)
  END DO
END DO            ! z- direction
!      write(6,*) ' after z D3MIXSTEP: norm=',wfnorm(psi)

RETURN
END SUBROUTINE d3mixpropag

!-----inv3p_ini-------------------------------------------------------------

SUBROUTINE inv3p_ini(deltim)

!     initialize Array of fixed values for invers_3r.
!     Only for same delta in all directions, yet

!USE params
!      include 'pot3D.inc'


REAL(DP), INTENT(IN OUT)                     :: deltim
COMPLEX(DP) :: invnum(2*(maxx+maxy+maxz)+3)
!inverse diagonal elements
COMPLEX(DP) :: diag
!diagonalelement of forward-matrix
COMMON /invnum3c/ invnum,diag



COMPLEX(DP) :: fac   !diagonal elements
INTEGER :: n     !loop


!---------------------------------------------------------------------------

WRITE(*,*) 'Inititalisierung der Koeffizienten fuer fin. Differenzen'

IF(dx /= dy .OR. dy /= dz .OR. dz /= dy)  &
    STOP 'only same gridstep in every direction'

fac = CMPLX(-2.0D0,2.0D0*dx*dx/deltim,DP)
invnum(1) = one/fac

DO n = 2, 2*(maxx+maxy+maxz)+3
  invnum (n) = one/(fac - invnum(n-1))
END DO

diag = CMPLX(2.0D0,2.0D0*dx*dx/deltim,DP)

!      write(6,*) ' init: diag,invnum=',diag,invnum

RETURN
END SUBROUTINE inv3p_ini


!-----kinprop_1d3--------------------------------------------------------------

SUBROUTINE kinprop_1d3 (psi, ndiml, inc, deltax, deltim)


!                        1
!     Solve = -------------------------  * Rhs
!             1 + 0.5 * i * deltim * T^


!     T^ is kinetic energy operator gives in 3 point precision
!     (finite differenc)
!     with boundary conditions ZERO at both ends (reflecting)

!USE params


COMPLEX(DP), INTENT(IN OUT)                  :: psi(ndiml)
INTEGER, INTENT(IN)                      :: ndiml
INTEGER, INTENT(IN)                      :: inc
REAL(DP), INTENT(IN OUT)                     :: deltax
REAL(DP), INTENT(IN OUT)                     :: deltim
COMPLEX(DP) :: invnum(2*(maxx+maxy+maxz)+3)
COMPLEX(DP) :: diag
COMMON /invnum3c/ invnum,diag

!INTEGER :: ! array size
!COMPLEX(DP) :: ! wave function to be propagated
COMPLEX(DP) :: solve
COMPLEX(DP) :: psip,psim,psi0
!REAL(DP):: ! mesh size

!INTEGER ::               ! Work array size
INTEGER:: ndimx

COMPLEX(DP) :: reff (ndimx)      ! effective r.h.s.
!c     $                  offdiag,           ! constant off-diagonal el.
!c     $                  offdiag2           ! squared off-diagonal el.

INTEGER :: i, n, ninc            ! loop index

!----------------------------------------------------------------------------

!     direct step (1-i dt/2 H) on psi

ndimx = 2 * (maxx + maxy + maxz) + 3
ninc=1
psi0 = psi(ninc)
ninc=ninc+inc
psip = psi(ninc)
reff(1) = diag*psi0 -  psip

DO i = 2,ndiml-1
  psim = psi0
  psi0 = psip
  ninc = ninc+inc
  psip = psi(ninc)
  reff(i) = diag*psi0 - (psip+psim)
END DO

psim = psi0
psi0 = psip
reff(ndiml) = diag*psi0 -  psim

!     forward loop (eliminate lower off- diagonal elements)

DO n = 2, ndiml
  reff(n) = reff(n) - reff(n-1)*invnum(n-1)
END DO

!     backward loop: eliminate upper off-diagonal and solve

ninc  = (ndiml-1)*inc + 1
psi(ninc) = reff(ndiml)*invnum(ndiml)
solve = psi(ninc)
DO n = ndiml-1, 1, -1
  ninc = ninc - inc
  solve = (reff(n) - solve)*invnum(n)
  psi(ninc) = solve
END DO

RETURN
END SUBROUTINE kinprop_1d3













#endif
#if(numerov)
!-----ckin3D_5r------------------------------------------------------------

SUBROUTINE ckin3d(psi,dxpsi)

!USE params
!      include 'pot3D.inc'

!     computes kinetic energy (- Laplacian) 5 point finite differences
!     and reflecting boundary conditions



COMPLEX(DP), INTENT(IN OUT)                  :: psi(minx:maxx,miny:maxy,minz:maxz)
COMPLEX(DP), INTENT(OUT)                     :: dxpsi(minx:maxx,miny:maxy,minz:maxz)



!-----------------------------------------------------------------------


INTEGER :: ix, iy, iz              ! loop indices
INTEGER :: lengx,lengy,lengz


!-----------------------------------------------------------------------


!     number of mesh points in each direction

lengx  = maxx+maxx+1
lengy  = maxy+maxy+1
lengz  = maxz+maxz+1


!     reset accumulator


DO iz = minz, maxz
  DO iy = miny, maxy
    DO ix = minx, maxx
      dxpsi(ix,iy,iz) = zero
    END DO
  END DO
END DO


!     kinetic energy in x direction


DO iz = minz, maxz
  DO iy = miny, maxy
    CALL ckin1d_5r(psi(minx,iy,iz),dx,lengx,1, dxpsi(minx,iy,iz))
  END DO
END DO



!     kinetic energy in y direction



DO iz = minz, maxz
  DO ix = minx, maxx
    CALL ckin1d_5r(psi(ix,miny,iz),dy,lengy,lengx, dxpsi(ix,miny,iz))
  END DO
END DO


!     kinetic energy in z direction


DO iy = miny, maxy
  DO ix = minx, maxx
    CALL ckin1d_5r(psi(ix,iy,minz),dz,lengz,lengy*lengx, dxpsi(ix,iy,minz))
  END DO
END DO


RETURN
END SUBROUTINE ckin3d
!-----ckin1D_5r ----------------------------------------------------------

SUBROUTINE ckin1d_5r(psi,deltax,nmax,inc,dxpsi)

!USE params

!     computes kinetic energy in one dimension and accumulates
!     that on 'dxpsi'.
!     uses five point finite differences.

!     psi    = input array ( a 3D array in calling routine)
!     dxpsi  = output containing 2. derivative + input value
!     deltax = mesh spacing
!     nmax   = number of mesh point in given direction
!     inc    = increment which connects neighbouring mesh points


COMPLEX(DP), INTENT(IN OUT)                  :: psi(*)
REAL(DP), INTENT(IN OUT)                     :: deltax
INTEGER, INTENT(IN)                      :: nmax
INTEGER, INTENT(IN)                      :: inc
COMPLEX(DP), INTENT(OUT)                     :: dxpsi(*)

INTEGER :: ninc,inc2
INTEGER :: i
REAL(DP):: d2i
!COMPLEX(DP) :: ! Wellenfunktion psi
!COMPLEX(DP) :: ! Ableitung der Wellenfunktion psi

!-----------------------------------------------------------------------

d2i = -one/(deltax*deltax)/12.0D0
inc2  = inc+inc

ninc = 1
dxpsi(ninc) = dxpsi(ninc) +d2i*(-psi(ninc+inc2)+16.0D0*psi(ninc+inc)  &
    -30.0D0*psi(ninc) )

ninc=ninc+inc
dxpsi(ninc) = dxpsi(ninc) +d2i*(-psi(ninc+inc2)+16.0D0*psi(ninc+inc)  &
    -30.0D0*psi(ninc) +16.0D0*psi(ninc-inc)               )


DO i = 3,nmax-2
  ninc=ninc+inc
  dxpsi(ninc) = dxpsi(ninc) +d2i*(-psi(ninc+inc2)+16.0D0*psi(ninc+inc)  &
      -30.0D0*psi(ninc) +16.0D0*psi(ninc-inc)-psi(ninc-inc2))
END DO

ninc=ninc + inc
dxpsi(ninc) = dxpsi(ninc) +d2i*(                16.0D0*psi(ninc+inc)  &
    -30.0D0*psi(ninc) +16.0D0*psi(ninc-inc)-psi(ninc-inc2))

ninc=ninc + inc
dxpsi(ninc) = dxpsi(ninc) +d2i*(  &
    -30.0D0*psi(ninc) +16.0D0*psi(ninc-inc)-psi(ninc-inc2))


RETURN
END SUBROUTINE ckin1d_5r
!-----rkin3D_5r------------------------------------------------------------

SUBROUTINE rkin3d(psi,dxpsi)

!USE params
!      include 'pot3D.inc'

!     computes kinetic energy (- Laplacian) 5 point finite differences
!     and reflecting boundary conditions



REAL(DP), INTENT(IN OUT)                     :: psi(minx:maxx,miny:maxy,minz:maxz)
REAL(DP), INTENT(OUT)                        :: dxpsi(minx:maxx,miny:maxy,minz:maxz)



!-----------------------------------------------------------------------


INTEGER :: ix, iy, iz              ! loop indices
INTEGER :: lengx,lengy,lengz


!-----------------------------------------------------------------------


!     number of mesh points in each direction

lengx  = maxx+maxx+1
lengy  = maxy+maxy+1
lengz  = maxz+maxz+1


!     reset accumulator


DO iz = minz, maxz
  DO iy = miny, maxy
    DO ix = minx, maxx
      dxpsi(ix,iy,iz) = zero
    END DO
  END DO
END DO


!     kinetic energy in x direction


DO iz = minz, maxz
  DO iy = miny, maxy
    CALL rkin1d_5r(psi(minx,iy,iz),dx,lengx,1, dxpsi(minx,iy,iz))
  END DO
END DO



!     kinetic energy in y direction



DO iz = minz, maxz
  DO ix = minx, maxx
    CALL rkin1d_5r(psi(ix,miny,iz),dy,lengy,lengx, dxpsi(ix,miny,iz))
  END DO
END DO


!     kinetic energy in z direction


DO iy = miny, maxy
  DO ix = minx, maxx
    CALL rkin1d_5r(psi(ix,iy,minz),dz,lengz,lengy*lengx, dxpsi(ix,iy,minz))
  END DO
END DO


RETURN
END SUBROUTINE rkin3d
!-----rkin1D_5r ----------------------------------------------------------

SUBROUTINE rkin1d_5r(psi,deltax,nmax,inc,dxpsi)

!USE params

!     computes kinetic energy in one dimension and accumulates
!     that on 'dxpsi'.
!     uses five point finite differences.

!     psi    = input array ( a 3D array in calling routine)
!     dxpsi  = output containing 2. derivative + input value
!     deltax = mesh spacing
!     nmax   = number of mesh point in given direction
!     inc    = increment which connects neighbouring mesh points


REAL(DP), INTENT(IN OUT)                     :: psi(*)
REAL(DP), INTENT(IN OUT)                     :: deltax
INTEGER, INTENT(IN)                      :: nmax
INTEGER, INTENT(IN)                      :: inc
REAL(DP), INTENT(OUT)                        :: dxpsi(*)

INTEGER :: ninc,inc2
INTEGER :: i
REAL(DP):: d2i
!REAL(DP):: ! Wellenfunktion psi
!REAL(DP):: ! Ableitung der Wellenfunktion psi

!-----------------------------------------------------------------------

d2i = -one/(deltax*deltax)/12.0D0
inc2  = inc+inc

ninc = 1
dxpsi(ninc) = dxpsi(ninc) +d2i*(-psi(ninc+inc2)+16.0D0*psi(ninc+inc)  &
    -30.0D0*psi(ninc) )

ninc=ninc+inc
dxpsi(ninc) = dxpsi(ninc) +d2i*(-psi(ninc+inc2)+16.0D0*psi(ninc+inc)  &
    -30.0D0*psi(ninc) +16.0D0*psi(ninc-inc)               )


DO i = 3,nmax-2
  ninc=ninc+inc
  dxpsi(ninc) = dxpsi(ninc) +d2i*(-psi(ninc+inc2)+16.0D0*psi(ninc+inc)  &
      -30.0D0*psi(ninc) +16.0D0*psi(ninc-inc)-psi(ninc-inc2))
END DO

ninc=ninc + inc
dxpsi(ninc) = dxpsi(ninc) +d2i*(                16.0D0*psi(ninc+inc)  &
    -30.0D0*psi(ninc) +16.0D0*psi(ninc-inc)-psi(ninc-inc2))

ninc=ninc + inc
dxpsi(ninc) = dxpsi(ninc) +d2i*(  &
    -30.0D0*psi(ninc) +16.0D0*psi(ninc-inc)-psi(ninc-inc2))


RETURN
END SUBROUTINE rkin1d_5r

!-----d3mixpropag_5----------------------------------------------------

SUBROUTINE d3mixpropag (psi, neupsi, deltim)

!USE params
!      include 'pot3D.inc'

!     Propagates T Operator for mixed step:

!              1 - 0.5 * i * deltim * T^
!     neupsi = ------------------------  * psi
!              1 + 0.5 * i * deltim * T^

!     in 3 Dimensions with 5- point precision (Numerov)
!     with boundary conditions zero at doth ends (reflecting)

!   . The "Laplacian"
!     operator is here dimensionless with diagonal element -2
!     and off-diagonals 1.


!     list parameters:


COMPLEX(DP), INTENT(IN OUT)                  :: psi(minx:maxx, miny:maxy,minz:maxz)
COMPLEX(DP), INTENT(IN OUT)                  :: neupsi(minx:maxx, miny:maxy,minz:maxz)
REAL(DP), INTENT(IN OUT)                     :: deltim

!                      ! old Function

!                      ! new Function

! Timestep

!     variables:

INTEGER :: ix, iy, iz
! loop indices

!-----------------------------------------------------------------------

DO iz = minz, maxz
  DO iy = miny, maxy
    CALL kinprop_1d5(psi(minx,iy,iz),2*maxx+1,1,dx,deltim)
  END DO
END DO            ! x- direction

DO iz = minz, maxz
  DO ix = minx, maxx
    CALL kinprop_1d5(psi(ix,miny,iz),2*maxy+1,2*maxx+1, dy,deltim)
  END DO
END DO            ! y- direction

DO iy = miny, maxy
  DO ix = minx, maxx
    CALL kinprop_1d5(psi(ix,iy,minz),2*maxz+1,  &
        (2*maxx+1)*(2*maxy+1),dz,deltim)
  END DO
END DO            ! z- direction

RETURN
END SUBROUTINE d3mixpropag

!-----inv5p_ini-------------------------------------------------------------

SUBROUTINE inv5p_ini(deltim)

!     initialize Array of fixed values for kinprop__1d5.
!     Only for same delta in all directions, yet

!USE params
!      include 'pot3D.inc'


REAL(DP), INTENT(IN)                         :: deltim
COMPLEX(DP) :: invnum(2*(maxx+maxy+maxz)+3)
!inverse diagonal Elements
COMPLEX(DP) :: diag,offd
!matrix-elements for forward-matrix
COMMON /invnum5c/ invnum,diag,offd


REAL(DP):: fac1,fac2   !diagonal elements
COMPLEX(DP) :: fac
INTEGER :: n,dim   !loop



!---------------------------------------------------------------------------

WRITE(*,*) 'Inititalisierung der Koeffizienten fuer Numerov'

IF(dx /= dy .OR. dy /= dz .OR. dz /= dy)  &
    STOP 'only same gridstep in every direction'

fac1=(-72.0D0*deltim*deltim + 10.0D0*dx**4)/ (36.0D0*deltim*deltim + dx**4)
fac2=72.0D0*deltim*deltax*dx / (36.0D0*deltim*deltim + dx**4)
fac = CMPLX(fac1,fac2,DP)

invnum(1) = one/fac

DO n = 2, 2*(maxx+maxy+maxz)+3
  invnum(n) = one/(fac - invnum(n-1))
END DO


fac1 = (72.0D0*deltim*deltim + 10.0D0*dx**4)/ (36.0D0*deltim*deltim + dx**4)
fac2 = 48.0D0*deltim*deltax*dx / (36.0D0*deltim*deltim + dx**4)
diag = CMPLX(fac1,fac2,DP)

fac1=(-36.0D0*deltim*deltim + dx**4)/ (36.0D0*deltim*deltim + dx**4)
fac2=12.0D0*deltim*deltax*dx / (deltax**4 + 36.0D0*deltim*deltim)
offd = CMPLX(fac1,fac2,DP)



RETURN
END SUBROUTINE inv5p_ini


!-----kinprop_1d5--------------------------------------------------------------

SUBROUTINE kinprop_1d5 (psi, ndiml, inc, deltax, deltim)


!             1 - 0.5 * i * deltim * T
!     Solve = -------------------------  * Rhs
!             1 + 0.5 * i * deltim * T^


!     T^ is kinetic energy operator gives in 3 point precision
!     (finite differenc)
!     with boundary conditions ZERO at both ends (reflecting)

!USE params
!      include 'work3D.inc'


COMPLEX(DP), INTENT(IN OUT)                  :: psi(ndiml)
INTEGER, INTENT(IN)                      :: ndiml
INTEGER, INTENT(IN)                      :: inc
REAL(DP), INTENT(IN OUT)                     :: deltax
REAL(DP), INTENT(IN OUT)                     :: deltim
COMPLEX(DP) :: invnum(2*(maxx+maxy+maxz)+3)
COMPLEX(DP) :: diag,offd
COMMON /invnum5c/ invnum,diag,offd


!INTEGER :: ! array size
!COMPLEX(DP) :: ! wave function to be propagated
COMPLEX(DP) :: solve
COMPLEX(DP) :: psip,psim,psi0
!REAL(DP):: ! mesh size

!INTEGER ::               ! Work array size
INTEGER, PARAMETER :: ndimx = 2 * (maxx + maxy + maxz) + 3

COMPLEX(DP) :: reff (ndimx)      ! effective r.h.s.
!c     $                  offdiag,           ! constant off-diagonal el.
!c     $                  offdiag2           ! squared off-diagonal el.

INTEGER :: i, n, ninc            ! loop index
REAL(DP):: fac1,fac2,fac3,fac4

!----------------------------------------------------------------------------

!     direct step (1-i dt/2 H) on psi




ninc=1
psi0 = psi(ninc)
ninc=ninc+inc
psip = psi(ninc)
reff(1) = diag*psi0 + offd*psip

DO i = 2,ndiml-1
  psim = psi0
  psi0 = psip
  ninc = ninc+inc
  psip = psi(ninc)
  reff(i) = diag*psi0 + offd*(psip+psim)
END DO

psim = psi0
psi0 = psip
reff(ndiml) = diag*psi0 +  offd*psim


!     forward loop (eliminate lower off- diagonal elements)

DO n = 2, ndiml
  reff(n) = reff(n) - reff(n-1)*invnum(n-1)
END DO

!     backward loop: eliminate upper off-diagonal and solve

ninc  = (ndiml-1)*inc + 1
psi(ninc) = reff(ndiml)*invnum(ndiml)
solve = psi(ninc)
DO n = ndiml-1, 1, -1
  ninc = ninc - inc
  solve = (reff(n) - solve)*invnum(n)
  psi(ninc) = solve
END DO

RETURN
END SUBROUTINE kinprop_1d5
#endif
END MODULE kinetic



