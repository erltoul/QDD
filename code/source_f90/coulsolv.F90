#include"define.h"
 
#if(gridfft)
#if(coufou)
INCLUDE 'falr.F90'
#endif
 
#if(coudoub)
INCLUDE 'coulex.F90'
#endif
#endif

#if(findiff|numerov)
!INCLUDE 'gridcoul.F90'
INCLUDE "findiff-sinft.F90"
#endif
