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

PROGRAM master_MPI

INCLUDE 'mpif.h'
INTEGER :: is(mpi_status_size)
INTEGER,PARAMETER :: knode=2

WRITE(*,*) ' before mpi_init'
CALL mpi_init(icode)

WRITE(*,*) ' before mpi_comm_size'
CALL  mpi_comm_size(mpi_comm_world,nprocs,icode)

WRITE(*,*) ' before mpi_comm_rank'
CALL  mpi_comm_rank(mpi_comm_world,myn,icode)
WRITE(*,*) ' nprocs,myn=',nprocs,myn

IF(nprocs<knode) STOP " not enough nodes"

IF(myn==0) THEN
  WRITE(*,*) ' submitted myn=',myn
  call system('./testsub.bat 001')
ELSE IF(myn==1) THEN
  WRITE(*,*) ' submitted myn=',myn
  call system('./testsub.bat 010')
ELSE IF(myn==2) THEN
  WRITE(*,*) ' submitted myn=',myn
  call system('./testsub.bat 011')
ELSE IF(myn==3) THEN
  WRITE(*,*) ' submitted myn=',myn
  call system('./testsub.bat 100')
ELSE IF(myn==4) THEN
  WRITE(*,*) ' submitted myn=',myn
  call system('./testsub.bat 101')
ELSE IF(myn==5) THEN
  WRITE(*,*) ' submitted myn=',myn
  call system('./testsub.bat 110')
END IF
  

CALL mpi_barrier (mpi_comm_world, mpi_ierror)
CALL mpi_finalize(icode)

END PROGRAM master_MPI
