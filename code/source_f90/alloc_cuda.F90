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
! cudaFreeHost
 integer (C_INT) function cudaFree(buffer)  bind(C,name="cudaFree")
  use iso_c_binding
  implicit none
  type (C_PTR), value :: buffer
 end function cudaFree
 END INTERFACE
END MODULE cuda_alloc
