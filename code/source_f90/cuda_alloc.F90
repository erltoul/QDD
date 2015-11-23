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

MODULE cuda_alloc
 INTERFACE
! cudaMallocHost
 integer (C_INT) function cudaMallocHost(buffer, size)  bind(C,name="cudaMallocHost")
  use iso_c_binding
  implicit none
  type (C_PTR)  :: buffer
  integer (C_LONG), value :: size
 end function cudaMallocHost
! cudaFreeHost
 integer (C_INT) function cudaFreeHost(buffer)  bind(C,name="cudaFreeHost")
  use iso_c_binding
  implicit none
  type (C_PTR), value :: buffer
 end function cudaFreeHost
! cudaMalloc
 integer (C_INT) function cudaMalloc(buffer, size)  bind(C,name="cudaMalloc")
  use iso_c_binding
  implicit none
  type (C_PTR)  :: buffer
  integer (C_LONG), value :: size
 end function cudaMalloc
! cudaFree
 integer (C_INT) function cudaFree(buffer)  bind(C,name="cudaFree")
  use iso_c_binding
  implicit none
  type (C_PTR), value :: buffer
 end function cudaFree
 END INTERFACE
END MODULE cuda_alloc
