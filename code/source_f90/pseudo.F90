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

!     the "periodic table"
!       (negative numbers stand for the local electron cloud of an element)

REAL(DP) :: ch(-99:99)=0D0,amu(-99:99)=0D0,cc1(-99:99)=0D0,cc2(-99:99)=0D0
REAL(DP) :: crloc(-99:99)=0D0,crs(-99:99)=0D0,chs(-99:99)=0D0,chg1(-99:99)=0D0,chg2(-99:99)=0D0
REAL(DP) :: sgm1(-99:99)=0D0,sgm2(-99:99)=0D0
REAL(DP) :: r0g(-99:99)=0D0,r1g(-99:99)=0D0,r2g(-99:99)=0D0,h0_11g(-99:99)=0D0,h0_22g(-99:99)=0D0
REAL(DP) :: h0_33g(-99:99)=0D0,h1_11g(-99:99)=0D0,h1_22g(-99:99)=0D0
REAL(DP) :: h2_11g(-99:99)=0D0,radiong(-99:99)=0D0
REAL(DP) :: h0_12g(-99:99)=-1D20               ! default signal "automatic"
!INTEGER :: nrow(-99:99)
!fix! INTEGER :: np(0:ng)

INTEGER,PARAMETER :: knl=18000       ! storage for PsP projectors
!fix! LOGICAL :: tblock(0:ng)
!fix! REAL(DP) :: p0_1(knl,0:ng),p0_2(knl,0:ng),p1_1(knl,0:ng),p1_1x(knl,0:ng)
!fix! REAL(DP) :: p1_1y(knl,0:ng),p1_1z(knl,0:ng) 
!fix! REAL(DP) :: p0_3(knl,0:ng),p1_2(knl,0:ng),p1_2x(knl,0:ng) 
!fix! REAL(DP) :: p1_2y(knl,0:ng),p1_2z(knl,0:ng) 
!fix! REAL(DP) :: p2_1(knl,0:ng),p2_xy(knl,0:ng),p2_xz(knl,0:ng) 
!fix! REAL(DP) :: p2_yz(knl,0:ng),p2_xy2(knl,0:ng),p2_z2(knl,0:ng)

LOGICAL :: tnonlocany                          ! flag to invoke non-local part
LOGICAL,ALLOCATABLE :: tnonloc(:)              ! flag to invoke non-local part
INTEGER,ALLOCATABLE :: ifin(:),icount(:,:),np(:)
REAL(DP),ALLOCATABLE :: p0_1(:,:),p0_2(:,:),p1_1(:,:),p1_1x(:,:)
REAL(DP),ALLOCATABLE :: p1_1y(:,:),p1_1z(:,:) 
REAL(DP),ALLOCATABLE :: p0_3(:,:),p1_2(:,:),p1_2x(:,:) 
REAL(DP),ALLOCATABLE :: p1_2y(:,:),p1_2z(:,:) 
REAL(DP),ALLOCATABLE :: p2_1(:,:),p2_xy(:,:),p2_xz(:,:) 
REAL(DP),ALLOCATABLE :: p2_yz(:,:),p2_xy2(:,:),p2_z2(:,:)
