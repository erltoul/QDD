# Choose compilation parameters

# CF90: Fortran compiler
# Available options:
#   * GFORT (serial version only)
#   * IFORT (serial & parallel & debug version)
#   * XLF_R (serial & parallel version)
#   * MPIF90 (parallel version only)
CF90 = IFORT

# MPI_PARALLEL: Variable for parallel processing (only for some compilers)
# Available options:
#   * YES
#   * SIM: simulate parallelization (simultaneous mono jobs)
#   * NO
MPI_PARALLEL = NO

# OMP_THREADS: Invoke OpenMP
# Available options:
#   * DYN: wave function parallelization (threads)
#   * YES: use threads for FFT
#   * NO
OMP_THREADS = NO

# TYPE_FFT: FFT solver
# Available options:
#   * NETLIB
#   * FFTW
#   * MKL
#   * cuFFT
TYPE_FFT = cuFFT

# DEBUG: enable debugging
# Available options:
#   * YES
#   * NO
DEBUG  = NO

# LINK_STATIC: select static linkage of the binary
# Available options:
#   * YES
#   * NO
LINK_STATIC = NO

# MKL_THREADS: Enable MKL threading (only used for TYPE_FFT = cuFFT
# Available options:
#   * YES
#   * NO
MKL_THREADS = NO

# MKL_WRAPPERS: Path to the MKL FFTW wrappers (only used for TYPE_FFT = cuFFT
# Default value:
MKL_WRAPPERS =

# MKLPATH: Path to the MKL libraries (only used for TYPE_FFT = cuFFT
# Default value:
MKLPATH =
# Path for hyperion:
# MKLPATH = /divers/intel/Compiler/12.1/mkl/lib/intel64

# MKLINCLUDE: Include path for the MKL (only used for TYPE_FFT = cuFFT
# Default value:
MKLINCLUDE =

# For large boxe sizes, try option
# '-i-dynamic -mcmodel=medium' or '-shared-intel -mcmodel=medium'


