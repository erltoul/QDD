#include 'define.h'

#if(gridfft)
#include 'fft.F90'
#endif 

#if(findiff|numerov)
INCLUDE "findiff.F90"
#endif

