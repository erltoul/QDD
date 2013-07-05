#include "define.h"
 
#if(gridfft)
#if(coufou)

#include "falr.F90"

#endif !coufou

#if(coudoub)

#include "coulex.F90"

#endif !coudoub
#endif !gridfft

#if(findiff|numerov)
!INCLUDE "gridcoul.F90"
INCLUDE "findiff-sinft.F90"
#endif
