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
 END INTERFACE
END MODULE cuda_alloc
