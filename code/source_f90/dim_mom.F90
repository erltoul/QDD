PARAMETER (pi=3.1415926)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2010-04-16  Time: 12:34:36
 
DO i=1,999
  READ(5,*,ERR=99,END=99) n,xx,yy,zz
  rms   = (xx+yy+zz)
  alpha = (zz+zz-xx-yy)*SQRT(5.0/(16.0*pi))/(rms)
  beta  = (xx-yy)*SQRT(15.0/(8.0*pi))/(rms)
  beta2 = SQRT(alpha*alpha+2.0*beta*beta)
  gamma = ATAN(1.4142136*beta/alpha)*180.0/pi
  IF(alpha < 0.0) gamma=180.0+gamma
  gamma = ABS(gamma)
  gamma = MOD(gamma,120.0)
  IF(gamma > 60.0 .AND. gamma <= 120.0) gamma = 120.0-gamma
  rms   = SQRT(rms)
  WRITE(6,'(1x,i5,f8.2,f8.3,f8.2)')  n,rms,beta2,gamma
END DO
99   CONTINUE
STOP
END SUBROUTINE pois

