!     the "periodic table"
!       (negative numbers stand for the local electron cloud of an element)

REAL(DP) :: ch(-92:92),amu(-92:92),cc1(-92:92),cc2(-92:92)
REAL(DP) :: crloc(-92:92),crs(-92:92),chs(-92:92),chg1(-92:92),chg2(-92:92)
REAL(DP) :: sgm1(-92:92),sgm2(-92:92)
REAL(DP) :: r0g(-92:92),r1g(-92:92),r2g(-92:92),h0_11g(-92:92),h0_22g(-92:92)
REAL(DP) :: h0_33g(-92:92),h1_11g(-92:92),h1_22g(-92:92)
REAL(DP) :: h2_11g(-92:92),radiong(-92:92)
REAL(DP) :: h0_12g(-92:92)=-1D20               ! default signal "automatic"
!INTEGER :: nrow(-92:92)
!fix! INTEGER :: np(0:ng)

INTEGER,PARAMETER :: knl=4000       ! storage for PsP projectors
!fix! LOGICAL :: tblock(0:ng)
!fix! REAL(DP) :: p0_1(knl,0:ng),p0_2(knl,0:ng),p1_1(knl,0:ng),p1_1x(knl,0:ng)
!fix! REAL(DP) :: p1_1y(knl,0:ng),p1_1z(knl,0:ng) 
!fix! REAL(DP) :: p0_3(knl,0:ng),p1_2(knl,0:ng),p1_2x(knl,0:ng) 
!fix! REAL(DP) :: p1_2y(knl,0:ng),p1_2z(knl,0:ng) 
!fix! REAL(DP) :: p2_1(knl,0:ng),p2_xy(knl,0:ng),p2_xz(knl,0:ng) 
!fix! REAL(DP) :: p2_yz(knl,0:ng),p2_xy2(knl,0:ng),p2_z2(knl,0:ng)

LOGICAL,SAVE :: tnonlocany                     ! flag to invoke non-local part
LOGICAL,ALLOCATABLE :: tnonloc(:)              ! flag to invoke non-local part
INTEGER,ALLOCATABLE :: ifin(:),icount(:,:),np(:)
REAL(DP),ALLOCATABLE :: p0_1(:,:),p0_2(:,:),p1_1(:,:),p1_1x(:,:)
REAL(DP),ALLOCATABLE :: p1_1y(:,:),p1_1z(:,:) 
REAL(DP),ALLOCATABLE :: p0_3(:,:),p1_2(:,:),p1_2x(:,:) 
REAL(DP),ALLOCATABLE :: p1_2y(:,:),p1_2z(:,:) 
REAL(DP),ALLOCATABLE :: p2_1(:,:),p2_xy(:,:),p2_xz(:,:) 
REAL(DP),ALLOCATABLE :: p2_yz(:,:),p2_xy2(:,:),p2_z2(:,:)
