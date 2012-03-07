MODULE CUDAFFT
USE, intrinsic :: iso_c_binding
implicit none

INTERFACE
  type(C_PTR) function cuda_plan_1d(plan,N,batch) bind(C,name='cuda_plan_1d')
    import
    INTEGER*8 :: plan
    INTEGER(C_INT), value :: N,batch
  end function cuda_plan_1d

  function cuda_plan_3d(plan,n1,n2,n3) bind(C,name='cuda_plan_3d')
    import
    INTEGER*8 :: plan
    INTEGER(C_INT), value :: n1,n2
  end function cuda_plan_3d

  type(C_PTR) function kill_plan(plan) bind(C,name='kill_plan')
    import
    INTEGER*8 :: plan
  end function kill_plan

  type(C_PTR) function run_fft_for(plan,a,b,Np) bind(C,name='run_fft_for')
    import
    INTEGER*8 :: plan
    INTEGER(C_INT), value :: Np
    COMPLEX(C_DOUBLE_COMPLEX),intent(in) :: a
    COMPLEX(C_DOUBLE_COMPLEX),intent(out) :: b
  end function run_fft_for

  type(C_PTR) function run_fft_back(plan,a,b,Np) bind(C,name='run_fft_back')
    import
    INTEGER*8 :: plan
    INTEGER(C_INT), value :: Np
    COMPLEX(C_DOUBLE_COMPLEX),intent(in) :: a
    COMPLEX(C_DOUBLE_COMPLEX),intent(out) :: b
  end function run_fft_back

  function run_fft_for3d(plan,a,b,Npx,Npy,Npz) bind(C,name='run_fft_for3d')
    import
    INTEGER*8 :: plan
    INTEGER(C_INT), value :: Npx,Npy,Npz
    COMPLEX(C_DOUBLE_COMPLEX),intent(in) :: a
    COMPLEX(C_DOUBLE_COMPLEX),intent(out) :: b
  end function run_fft_for3d

  function run_fft_back3d(plan,a,b,Npx,Npy,Npz) bind(C,name='run_fft_back3d')
    import
    INTEGER*8 :: plan
    INTEGER(C_INT), value :: Npx,Npy,Npz
    COMPLEX(C_DOUBLE_COMPLEX),intent(in) :: a
    COMPLEX(C_DOUBLE_COMPLEX),intent(out) :: b
  end function run_fft_back3d
END INTERFACE

END MODULE CUDAFFT
