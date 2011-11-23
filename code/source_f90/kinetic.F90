#include"define.h"

#if(gridfft)
INCLUDE "fft.F90"
#endif

#if(findiff|numerov)
INCLUDE "findiff.F90"
#endif

