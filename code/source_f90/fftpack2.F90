SUBROUTINE cosqb (n,x,wsave)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2010-04-16  Time: 12:36:26

INTEGER, INTENT(IN)                      :: n
REAL(8), INTENT(OUT)                        :: x(10)
REAL(8), INTENT(IN OUT)                     :: wsave(1)

DATA tsqrt2 /2.82842712474619/

IF (n-2 < 0) THEN
  GO TO   101
ELSE IF (n-2 == 0) THEN
  GO TO   102
ELSE
  GO TO   103
END IF
101 x(1) = 4.*x(1)
RETURN
102 x1 = 4.*(x(1)+x(2))
x(2) = tsqrt2*(x(1)-x(2))
x(1) = x1
RETURN
103 CALL cosqb1 (n,x,wsave,wsave(n+1))
RETURN
END SUBROUTINE cosqb

SUBROUTINE cosqb1 (n,x,w,xh)


INTEGER, INTENT(IN)                      :: n
REAL(8), INTENT(IN OUT)                     :: x(1)
REAL(8), INTENT(IN)                         :: w(1)
REAL(8), INTENT(OUT)                        :: xh(1)

ns2 = (n+1)/2
np2 = n+2
DO  i=3,n,2
  xim1 = x(i-1)+x(i)
  x(i) = x(i)-x(i-1)
  x(i-1) = xim1
END DO
x(1) = x(1)+x(1)
modn = MOD(n,2)
IF (modn == 0) x(n) = x(n)+x(n)
CALL rfftb (n,x,xh)
DO  k=2,ns2
  kc = np2-k
  xh(k) = w(k-1)*x(kc)+w(kc-1)*x(k)
  xh(kc) = w(k-1)*x(k)-w(kc-1)*x(kc)
END DO
IF (modn == 0) x(ns2+1) = w(ns2)*(x(ns2+1)+x(ns2+1))
DO  k=2,ns2
  kc = np2-k
  x(k) = xh(k)+xh(kc)
  x(kc) = xh(k)-xh(kc)
END DO
x(1) = x(1)+x(1)
RETURN
END SUBROUTINE cosqb1

SUBROUTINE radb2 (ido,l1,cc,ch,wa1)


INTEGER, INTENT(IN)                      :: ido
INTEGER, INTENT(IN)                      :: l1
REAL(8), INTENT(IN)                         :: cc(ido,2,l1)
REAL(8), INTENT(OUT)                        :: ch(ido,l1,2)
REAL(8), INTENT(IN)                         :: wa1(1)

DO  k=1,l1
  ch(1,k,1) = cc(1,1,k)+cc(ido,2,k)
  ch(1,k,2) = cc(1,1,k)-cc(ido,2,k)
END DO
IF (ido-2 < 0) THEN
  GO TO   107
ELSE IF (ido-2 == 0) THEN
  GO TO   105
END IF
102 idp2 = ido+2
DO  k=1,l1
  DO  i=3,ido,2
    ic = idp2-i
    ch(i-1,k,1) = cc(i-1,1,k)+cc(ic-1,2,k)
    tr2 = cc(i-1,1,k)-cc(ic-1,2,k)
    ch(i,k,1) = cc(i,1,k)-cc(ic,2,k)
    ti2 = cc(i,1,k)+cc(ic,2,k)
    ch(i-1,k,2) = wa1(i-2)*tr2-wa1(i-1)*ti2
    ch(i,k,2) = wa1(i-2)*ti2+wa1(i-1)*tr2
  END DO
END DO
IF (MOD(ido,2) == 1) RETURN
105 DO  k=1,l1
  ch(ido,k,1) = cc(ido,1,k)+cc(ido,1,k)
  ch(ido,k,2) = -(cc(1,2,k)+cc(1,2,k))
END DO
107 RETURN
END SUBROUTINE radb2

SUBROUTINE radb3 (ido,l1,cc,ch,wa1,wa2)

INTEGER, INTENT(IN)                      :: ido
INTEGER, INTENT(IN)                      :: l1
REAL(8), INTENT(IN)                         :: cc(ido,3,l1)
REAL(8), INTENT(OUT)                        :: ch(ido,l1,3)
REAL(8), INTENT(IN)                         :: wa1(1)
REAL(8), INTENT(IN)                         :: wa2(1)

DATA taur,taui /-.5,.866025403784439/

DO  k=1,l1
  tr2 = cc(ido,2,k)+cc(ido,2,k)
  cr2 = cc(1,1,k)+taur*tr2
  ch(1,k,1) = cc(1,1,k)+tr2
  ci3 = taui*(cc(1,3,k)+cc(1,3,k))
  ch(1,k,2) = cr2-ci3
  ch(1,k,3) = cr2+ci3
END DO
IF (ido == 1) RETURN
idp2 = ido+2
DO  k=1,l1
  DO  i=3,ido,2
    ic = idp2-i
    tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
    cr2 = cc(i-1,1,k)+taur*tr2
    ch(i-1,k,1) = cc(i-1,1,k)+tr2
    ti2 = cc(i,3,k)-cc(ic,2,k)
    ci2 = cc(i,1,k)+taur*ti2
    ch(i,k,1) = cc(i,1,k)+ti2
    cr3 = taui*(cc(i-1,3,k)-cc(ic-1,2,k))
    ci3 = taui*(cc(i,3,k)+cc(ic,2,k))
    dr2 = cr2-ci3
    dr3 = cr2+ci3
    di2 = ci2+cr3
    di3 = ci2-cr3
    ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
    ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
    ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
    ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
  END DO
END DO
RETURN
END SUBROUTINE radb3

SUBROUTINE radb4 (ido,l1,cc,ch,wa1,wa2,wa3)

INTEGER, INTENT(IN)                      :: ido
INTEGER, INTENT(IN)                      :: l1
REAL(8), INTENT(IN)                         :: cc(ido,4,l1)
REAL(8), INTENT(OUT)                        :: ch(ido,l1,4)
REAL(8), INTENT(IN)                         :: wa1(1)
REAL(8), INTENT(IN)                         :: wa2(1)
REAL(8), INTENT(IN)                         :: wa3(1)

DATA sqrt2 /1.414213562373095/

DO  k=1,l1
  tr1 = cc(1,1,k)-cc(ido,4,k)
  tr2 = cc(1,1,k)+cc(ido,4,k)
  tr3 = cc(ido,2,k)+cc(ido,2,k)
  tr4 = cc(1,3,k)+cc(1,3,k)
  ch(1,k,1) = tr2+tr3
  ch(1,k,2) = tr1-tr4
  ch(1,k,3) = tr2-tr3
  ch(1,k,4) = tr1+tr4
END DO
IF (ido-2 < 0) THEN
  GO TO   107
ELSE IF (ido-2 == 0) THEN
  GO TO   105
END IF
102 idp2 = ido+2
DO  k=1,l1
  DO  i=3,ido,2
    ic = idp2-i
    ti1 = cc(i,1,k)+cc(ic,4,k)
    ti2 = cc(i,1,k)-cc(ic,4,k)
    ti3 = cc(i,3,k)-cc(ic,2,k)
    tr4 = cc(i,3,k)+cc(ic,2,k)
    tr1 = cc(i-1,1,k)-cc(ic-1,4,k)
    tr2 = cc(i-1,1,k)+cc(ic-1,4,k)
    ti4 = cc(i-1,3,k)-cc(ic-1,2,k)
    tr3 = cc(i-1,3,k)+cc(ic-1,2,k)
    ch(i-1,k,1) = tr2+tr3
    cr3 = tr2-tr3
    ch(i,k,1) = ti2+ti3
    ci3 = ti2-ti3
    cr2 = tr1-tr4
    cr4 = tr1+tr4
    ci2 = ti1+ti4
    ci4 = ti1-ti4
    ch(i-1,k,2) = wa1(i-2)*cr2-wa1(i-1)*ci2
    ch(i,k,2) = wa1(i-2)*ci2+wa1(i-1)*cr2
    ch(i-1,k,3) = wa2(i-2)*cr3-wa2(i-1)*ci3
    ch(i,k,3) = wa2(i-2)*ci3+wa2(i-1)*cr3
    ch(i-1,k,4) = wa3(i-2)*cr4-wa3(i-1)*ci4
    ch(i,k,4) = wa3(i-2)*ci4+wa3(i-1)*cr4
  END DO
END DO
IF (MOD(ido,2) == 1) RETURN
105 CONTINUE
DO  k=1,l1
  ti1 = cc(1,2,k)+cc(1,4,k)
  ti2 = cc(1,4,k)-cc(1,2,k)
  tr1 = cc(ido,1,k)-cc(ido,3,k)
  tr2 = cc(ido,1,k)+cc(ido,3,k)
  ch(ido,k,1) = tr2+tr2
  ch(ido,k,2) = sqrt2*(tr1-ti1)
  ch(ido,k,3) = ti2+ti2
  ch(ido,k,4) = -sqrt2*(tr1+ti1)
END DO
107 RETURN
END SUBROUTINE radb4

SUBROUTINE radb5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)

INTEGER, INTENT(IN)                      :: ido
INTEGER, INTENT(IN)                      :: l1
REAL(8), INTENT(IN)                         :: cc(ido,5,l1)
REAL(8), INTENT(OUT)                        :: ch(ido,l1,5)
REAL(8), INTENT(IN)                         :: wa1(1)
REAL(8), INTENT(IN)                         :: wa2(1)
REAL(8), INTENT(IN)                         :: wa3(1)
REAL(8), INTENT(IN)                         :: wa4(1)

DATA tr11,ti11,tr12,ti12 /.309016994374947,.951056516295154,  &
    -.809016994374947,.587785252292473/

DO  k=1,l1
  ti5 = cc(1,3,k)+cc(1,3,k)
  ti4 = cc(1,5,k)+cc(1,5,k)
  tr2 = cc(ido,2,k)+cc(ido,2,k)
  tr3 = cc(ido,4,k)+cc(ido,4,k)
  ch(1,k,1) = cc(1,1,k)+tr2+tr3
  cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
  cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
  ci5 = ti11*ti5+ti12*ti4
  ci4 = ti12*ti5-ti11*ti4
  ch(1,k,2) = cr2-ci5
  ch(1,k,3) = cr3-ci4
  ch(1,k,4) = cr3+ci4
  ch(1,k,5) = cr2+ci5
END DO
IF (ido == 1) RETURN
idp2 = ido+2
DO  k=1,l1
  DO  i=3,ido,2
    ic = idp2-i
    ti5 = cc(i,3,k)+cc(ic,2,k)
    ti2 = cc(i,3,k)-cc(ic,2,k)
    ti4 = cc(i,5,k)+cc(ic,4,k)
    ti3 = cc(i,5,k)-cc(ic,4,k)
    tr5 = cc(i-1,3,k)-cc(ic-1,2,k)
    tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
    tr4 = cc(i-1,5,k)-cc(ic-1,4,k)
    tr3 = cc(i-1,5,k)+cc(ic-1,4,k)
    ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
    ch(i,k,1) = cc(i,1,k)+ti2+ti3
    cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
    ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
    cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
    ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
    cr5 = ti11*tr5+ti12*tr4
    ci5 = ti11*ti5+ti12*ti4
    cr4 = ti12*tr5-ti11*tr4
    ci4 = ti12*ti5-ti11*ti4
    dr3 = cr3-ci4
    dr4 = cr3+ci4
    di3 = ci3+cr4
    di4 = ci3-cr4
    dr5 = cr2+ci5
    dr2 = cr2-ci5
    di5 = ci2-cr5
    di2 = ci2+cr5
    ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
    ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
    ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
    ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
    ch(i-1,k,4) = wa3(i-2)*dr4-wa3(i-1)*di4
    ch(i,k,4) = wa3(i-2)*di4+wa3(i-1)*dr4
    ch(i-1,k,5) = wa4(i-2)*dr5-wa4(i-1)*di5
    ch(i,k,5) = wa4(i-2)*di5+wa4(i-1)*dr5
  END DO
END DO
RETURN
END SUBROUTINE radb5

SUBROUTINE radbg (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)

INTEGER, INTENT(IN)                      :: ido
INTEGER, INTENT(IN)                      :: ip
INTEGER, INTENT(IN)                      :: l1
INTEGER, INTENT(IN)                      :: idl1
REAL(8), INTENT(IN)                         :: cc(ido,ip,l1)
REAL(8), INTENT(IN OUT)                     :: c1(ido,l1,ip)
REAL(8), INTENT(OUT)                        :: c2(idl1,ip)
REAL(8), INTENT(OUT)                        :: ch(ido,l1,ip)
REAL(8), INTENT(IN OUT)                     :: ch2(idl1,ip)
REAL(8), INTENT(IN)                         :: wa(1)

DATA tpi/6.28318530717959/

arg = tpi/FLOAT(ip)
dcp = COS(arg)
dsp = SIN(arg)
idp2 = ido+2
nbd = (ido-1)/2
ipp2 = ip+2
ipph = (ip+1)/2
IF (ido < l1) GO TO 103
DO  k=1,l1
  DO  i=1,ido
    ch(i,k,1) = cc(i,1,k)
  END DO
END DO
GO TO 106
103 DO  i=1,ido
  DO  k=1,l1
    ch(i,k,1) = cc(i,1,k)
  END DO
END DO
106 DO  j=2,ipph
  jc = ipp2-j
  j2 = j+j
  DO  k=1,l1
    ch(1,k,j) = cc(ido,j2-2,k)+cc(ido,j2-2,k)
    ch(1,k,jc) = cc(1,j2-1,k)+cc(1,j2-1,k)
  END DO
END DO
IF (ido == 1) GO TO 116
IF (nbd < l1) GO TO 112
DO  j=2,ipph
  jc = ipp2-j
  DO  k=1,l1
    DO  i=3,ido,2
      ic = idp2-i
      ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
      ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
      ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
      ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
    END DO
  END DO
END DO
GO TO 116
112 DO  j=2,ipph
  jc = ipp2-j
  DO  i=3,ido,2
    ic = idp2-i
    DO  k=1,l1
      ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
      ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
      ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
      ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
    END DO
  END DO
END DO
116 ar1 = 1.
ai1 = 0.
DO  l=2,ipph
  lc = ipp2-l
  ar1h = dcp*ar1-dsp*ai1
  ai1 = dcp*ai1+dsp*ar1
  ar1 = ar1h
  DO  ik=1,idl1
    c2(ik,l) = ch2(ik,1)+ar1*ch2(ik,2)
    c2(ik,lc) = ai1*ch2(ik,ip)
  END DO
  dc2 = ar1
  ds2 = ai1
  ar2 = ar1
  ai2 = ai1
  DO  j=3,ipph
    jc = ipp2-j
    ar2h = dc2*ar2-ds2*ai2
    ai2 = dc2*ai2+ds2*ar2
    ar2 = ar2h
    DO  ik=1,idl1
      c2(ik,l) = c2(ik,l)+ar2*ch2(ik,j)
      c2(ik,lc) = c2(ik,lc)+ai2*ch2(ik,jc)
    END DO
  END DO
END DO
DO  j=2,ipph
  DO  ik=1,idl1
    ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
  END DO
END DO
DO  j=2,ipph
  jc = ipp2-j
  DO  k=1,l1
    ch(1,k,j) = c1(1,k,j)-c1(1,k,jc)
    ch(1,k,jc) = c1(1,k,j)+c1(1,k,jc)
  END DO
END DO
IF (ido == 1) GO TO 132
IF (nbd < l1) GO TO 128
DO  j=2,ipph
  jc = ipp2-j
  DO  k=1,l1
    DO  i=3,ido,2
      ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
      ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
      ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
      ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
    END DO
  END DO
END DO
GO TO 132
128 DO  j=2,ipph
  jc = ipp2-j
  DO  i=3,ido,2
    DO  k=1,l1
      ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
      ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
      ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
      ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
    END DO
  END DO
END DO
132 CONTINUE
IF (ido == 1) RETURN
DO  ik=1,idl1
  c2(ik,1) = ch2(ik,1)
END DO
DO  j=2,ip
  DO  k=1,l1
    c1(1,k,j) = ch(1,k,j)
  END DO
END DO
IF (nbd > l1) GO TO 139
is = -ido
DO  j=2,ip
  is = is+ido
  idij = is
  DO  i=3,ido,2
    idij = idij+2
    DO  k=1,l1
      c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
      c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
    END DO
  END DO
END DO
GO TO 143
139 is = -ido
DO  j=2,ip
  is = is+ido
  DO  k=1,l1
    idij = is
    DO  i=3,ido,2
      idij = idij+2
      c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
      c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
    END DO
  END DO
END DO
143 RETURN
END SUBROUTINE radbg

SUBROUTINE rfftb (n,r,wsave)


INTEGER, INTENT(IN)                      :: n
REAL(8), INTENT(IN OUT)                     :: r(1)
REAL(8), INTENT(IN OUT)                     :: wsave(1)

IF (n == 1) RETURN
CALL rfftb1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
RETURN
END SUBROUTINE rfftb

SUBROUTINE rfftb1 (n,c,ch,wa,ifac)


INTEGER, INTENT(IN)                      :: n
REAL(8), INTENT(OUT)                        :: c(1)
REAL(8), INTENT(IN)                         :: ch(1)
REAL(8), INTENT(IN OUT)                     :: wa(1)
INTEGER, INTENT(IN)                      :: ifac(10)

nf = ifac(2)
na = 0
l1 = 1
iw = 1
DO  k1=1,nf
  ip = ifac(k1+2)
  l2 = ip*l1
  ido = n/l2
  idl1 = ido*l1
  IF (ip /= 4) GO TO 103
  ix2 = iw+ido
  ix3 = ix2+ido
  IF (na /= 0) GO TO 101
  CALL radb4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
  GO TO 102
  101    CALL radb4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
  GO TO 115
  103    IF (ip /= 2) GO TO 106
  IF (na /= 0) GO TO 104
  CALL radb2 (ido,l1,c,ch,wa(iw))
  GO TO 105
  104    CALL radb2 (ido,l1,ch,c,wa(iw))
  105    na = 1-na
  GO TO 115
  106    IF (ip /= 3) GO TO 109
  ix2 = iw+ido
  IF (na /= 0) GO TO 107
  CALL radb3 (ido,l1,c,ch,wa(iw),wa(ix2))
  GO TO 108
  107    CALL radb3 (ido,l1,ch,c,wa(iw),wa(ix2))
  108    na = 1-na
  GO TO 115
  109    IF (ip /= 5) GO TO 112
  ix2 = iw+ido
  ix3 = ix2+ido
  ix4 = ix3+ido
  IF (na /= 0) GO TO 110
  CALL radb5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  GO TO 111
  110    CALL radb5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
  GO TO 115
  112    IF (na /= 0) GO TO 113
  CALL radbg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
  GO TO 114
  113    CALL radbg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  114    IF (ido == 1) na = 1-na
  115    l1 = l2
  iw = iw+(ip-1)*ido
END DO
IF (na == 0) RETURN
DO  i=1,n
  c(i) = ch(i)
END DO
RETURN
END SUBROUTINE rfftb1










SUBROUTINE cosqf (n,x,wsave)

INTEGER, INTENT(IN)                      :: n
REAL(8), INTENT(IN OUT)                     :: x(10)
REAL(8), INTENT(IN OUT)                     :: wsave(1)

DATA sqrt2 /1.4142135623731/

IF (n-2 < 0) THEN
  GO TO   102
ELSE IF (n-2 == 0) THEN
  GO TO   101
ELSE
  GO TO   103
END IF
101 tsqx = sqrt2*x(2)
x(2) = x(1)-tsqx
x(1) = x(1)+tsqx
102 RETURN
103 CALL cosqf1 (n,x,wsave,wsave(n+1))
RETURN
END SUBROUTINE cosqf

SUBROUTINE cosqf1 (n,x,w,xh)


INTEGER, INTENT(IN)                      :: n
REAL(8), INTENT(IN OUT)                     :: x(1)
REAL(8), INTENT(IN)                         :: w(1)
REAL(8), INTENT(OUT)                        :: xh(1)

ns2 = (n+1)/2
np2 = n+2
DO  k=2,ns2
  kc = np2-k
  xh(k) = x(k)+x(kc)
  xh(kc) = x(k)-x(kc)
END DO
modn = MOD(n,2)
IF (modn == 0) xh(ns2+1) = x(ns2+1)+x(ns2+1)
DO  k=2,ns2
  kc = np2-k
  x(k) = w(k-1)*xh(kc)+w(kc-1)*xh(k)
  x(kc) = w(k-1)*xh(k)-w(kc-1)*xh(kc)
END DO
IF (modn == 0) x(ns2+1) = w(ns2)*xh(ns2+1)
CALL rfftf (n,x,xh)
DO  i=3,n,2
  xim1 = x(i-1)-x(i)
  x(i) = x(i-1)+x(i)
  x(i-1) = xim1
END DO
RETURN
END SUBROUTINE cosqf1

SUBROUTINE radf2 (ido,l1,cc,ch,wa1)


INTEGER, INTENT(IN)                      :: ido
INTEGER, INTENT(IN)                      :: l1
REAL(8), INTENT(IN)                         :: cc(ido,l1,2)
REAL(8), INTENT(OUT)                        :: ch(ido,2,l1)
REAL(8), INTENT(IN)                         :: wa1(1)

DO  k=1,l1
  ch(1,1,k) = cc(1,k,1)+cc(1,k,2)
  ch(ido,2,k) = cc(1,k,1)-cc(1,k,2)
END DO
IF (ido-2 < 0) THEN
  GO TO   107
ELSE IF (ido-2 == 0) THEN
  GO TO   105
END IF
102 idp2 = ido+2
DO  k=1,l1
  DO  i=3,ido,2
    ic = idp2-i
    tr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
    ti2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
    ch(i,1,k) = cc(i,k,1)+ti2
    ch(ic,2,k) = ti2-cc(i,k,1)
    ch(i-1,1,k) = cc(i-1,k,1)+tr2
    ch(ic-1,2,k) = cc(i-1,k,1)-tr2
  END DO
END DO
IF (MOD(ido,2) == 1) RETURN
105 DO  k=1,l1
  ch(1,2,k) = -cc(ido,k,2)
  ch(ido,1,k) = cc(ido,k,1)
END DO
107 RETURN
END SUBROUTINE radf2

SUBROUTINE radf3 (ido,l1,cc,ch,wa1,wa2)

INTEGER, INTENT(IN)                      :: ido
INTEGER, INTENT(IN)                      :: l1
REAL(8), INTENT(IN)                         :: cc(ido,l1,3)
REAL(8), INTENT(OUT)                        :: ch(ido,3,l1)
REAL(8), INTENT(IN)                         :: wa1(1)
REAL(8), INTENT(IN)                         :: wa2(1)

DATA taur,taui /-.5,.866025403784439/

DO  k=1,l1
  cr2 = cc(1,k,2)+cc(1,k,3)
  ch(1,1,k) = cc(1,k,1)+cr2
  ch(1,3,k) = taui*(cc(1,k,3)-cc(1,k,2))
  ch(ido,2,k) = cc(1,k,1)+taur*cr2
END DO
IF (ido == 1) RETURN
idp2 = ido+2
DO  k=1,l1
  DO  i=3,ido,2
    ic = idp2-i
    dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
    di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
    dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
    di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
    cr2 = dr2+dr3
    ci2 = di2+di3
    ch(i-1,1,k) = cc(i-1,k,1)+cr2
    ch(i,1,k) = cc(i,k,1)+ci2
    tr2 = cc(i-1,k,1)+taur*cr2
    ti2 = cc(i,k,1)+taur*ci2
    tr3 = taui*(di2-di3)
    ti3 = taui*(dr3-dr2)
    ch(i-1,3,k) = tr2+tr3
    ch(ic-1,2,k) = tr2-tr3
    ch(i,3,k) = ti2+ti3
    ch(ic,2,k) = ti3-ti2
  END DO
END DO
RETURN
END SUBROUTINE radf3

SUBROUTINE radf4 (ido,l1,cc,ch,wa1,wa2,wa3)

INTEGER, INTENT(IN)                      :: ido
INTEGER, INTENT(IN)                      :: l1
REAL(8), INTENT(IN)                         :: cc(ido,l1,4)
REAL(8), INTENT(OUT)                        :: ch(ido,4,l1)
REAL(8), INTENT(IN)                         :: wa1(1)
REAL(8), INTENT(IN)                         :: wa2(1)
REAL(8), INTENT(IN)                         :: wa3(1)

DATA hsqt2 /.7071067811865475/

DO  k=1,l1
  tr1 = cc(1,k,2)+cc(1,k,4)
  tr2 = cc(1,k,1)+cc(1,k,3)
  ch(1,1,k) = tr1+tr2
  ch(ido,4,k) = tr2-tr1
  ch(ido,2,k) = cc(1,k,1)-cc(1,k,3)
  ch(1,3,k) = cc(1,k,4)-cc(1,k,2)
END DO
IF (ido-2 < 0) THEN
  GO TO   107
ELSE IF (ido-2 == 0) THEN
  GO TO   105
END IF
102 idp2 = ido+2
DO  k=1,l1
  DO  i=3,ido,2
    ic = idp2-i
    cr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
    ci2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
    cr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
    ci3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
    cr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
    ci4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
    tr1 = cr2+cr4
    tr4 = cr4-cr2
    ti1 = ci2+ci4
    ti4 = ci2-ci4
    ti2 = cc(i,k,1)+ci3
    ti3 = cc(i,k,1)-ci3
    tr2 = cc(i-1,k,1)+cr3
    tr3 = cc(i-1,k,1)-cr3
    ch(i-1,1,k) = tr1+tr2
    ch(ic-1,4,k) = tr2-tr1
    ch(i,1,k) = ti1+ti2
    ch(ic,4,k) = ti1-ti2
    ch(i-1,3,k) = ti4+tr3
    ch(ic-1,2,k) = tr3-ti4
    ch(i,3,k) = tr4+ti3
    ch(ic,2,k) = tr4-ti3
  END DO
END DO
IF (MOD(ido,2) == 1) RETURN
105 CONTINUE
DO  k=1,l1
  ti1 = -hsqt2*(cc(ido,k,2)+cc(ido,k,4))
  tr1 = hsqt2*(cc(ido,k,2)-cc(ido,k,4))
  ch(ido,1,k) = tr1+cc(ido,k,1)
  ch(ido,3,k) = cc(ido,k,1)-tr1
  ch(1,2,k) = ti1-cc(ido,k,3)
  ch(1,4,k) = ti1+cc(ido,k,3)
END DO
107 RETURN
END SUBROUTINE radf4

SUBROUTINE radf5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)

INTEGER, INTENT(IN)                      :: ido
INTEGER, INTENT(IN)                      :: l1
REAL(8), INTENT(IN)                         :: cc(ido,l1,5)
REAL(8), INTENT(OUT)                        :: ch(ido,5,l1)
REAL(8), INTENT(IN)                         :: wa1(1)
REAL(8), INTENT(IN)                         :: wa2(1)
REAL(8), INTENT(IN)                         :: wa3(1)
REAL(8), INTENT(IN)                         :: wa4(1)

DATA tr11,ti11,tr12,ti12 /.309016994374947,.951056516295154,  &
    -.809016994374947,.587785252292473/

DO  k=1,l1
  cr2 = cc(1,k,5)+cc(1,k,2)
  ci5 = cc(1,k,5)-cc(1,k,2)
  cr3 = cc(1,k,4)+cc(1,k,3)
  ci4 = cc(1,k,4)-cc(1,k,3)
  ch(1,1,k) = cc(1,k,1)+cr2+cr3
  ch(ido,2,k) = cc(1,k,1)+tr11*cr2+tr12*cr3
  ch(1,3,k) = ti11*ci5+ti12*ci4
  ch(ido,4,k) = cc(1,k,1)+tr12*cr2+tr11*cr3
  ch(1,5,k) = ti12*ci5-ti11*ci4
END DO
IF (ido == 1) RETURN
idp2 = ido+2
DO  k=1,l1
  DO  i=3,ido,2
    ic = idp2-i
    dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
    di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
    dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
    di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
    dr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
    di4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
    dr5 = wa4(i-2)*cc(i-1,k,5)+wa4(i-1)*cc(i,k,5)
    di5 = wa4(i-2)*cc(i,k,5)-wa4(i-1)*cc(i-1,k,5)
    cr2 = dr2+dr5
    ci5 = dr5-dr2
    cr5 = di2-di5
    ci2 = di2+di5
    cr3 = dr3+dr4
    ci4 = dr4-dr3
    cr4 = di3-di4
    ci3 = di3+di4
    ch(i-1,1,k) = cc(i-1,k,1)+cr2+cr3
    ch(i,1,k) = cc(i,k,1)+ci2+ci3
    tr2 = cc(i-1,k,1)+tr11*cr2+tr12*cr3
    ti2 = cc(i,k,1)+tr11*ci2+tr12*ci3
    tr3 = cc(i-1,k,1)+tr12*cr2+tr11*cr3
    ti3 = cc(i,k,1)+tr12*ci2+tr11*ci3
    tr5 = ti11*cr5+ti12*cr4
    ti5 = ti11*ci5+ti12*ci4
    tr4 = ti12*cr5-ti11*cr4
    ti4 = ti12*ci5-ti11*ci4
    ch(i-1,3,k) = tr2+tr5
    ch(ic-1,2,k) = tr2-tr5
    ch(i,3,k) = ti2+ti5
    ch(ic,2,k) = ti5-ti2
    ch(i-1,5,k) = tr3+tr4
    ch(ic-1,4,k) = tr3-tr4
    ch(i,5,k) = ti3+ti4
    ch(ic,4,k) = ti4-ti3
  END DO
END DO
RETURN
END SUBROUTINE radf5

SUBROUTINE radfg (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)

INTEGER, INTENT(IN)                      :: ido
INTEGER, INTENT(IN)                      :: ip
INTEGER, INTENT(IN)                      :: l1
INTEGER, INTENT(IN)                      :: idl1
REAL(8), INTENT(OUT)                        :: cc(ido,ip,l1)
REAL(8), INTENT(IN OUT)                     :: c1(ido,l1,ip)
REAL(8), INTENT(IN OUT)                     :: c2(idl1,ip)
REAL(8), INTENT(OUT)                        :: ch(ido,l1,ip)
REAL(8), INTENT(OUT)                        :: ch2(idl1,ip)
REAL(8), INTENT(IN)                         :: wa(1)

DATA tpi/6.28318530717959/

arg = tpi/FLOAT(ip)
dcp = COS(arg)
dsp = SIN(arg)
ipph = (ip+1)/2
ipp2 = ip+2
idp2 = ido+2
nbd = (ido-1)/2
IF (ido == 1) GO TO 119
DO  ik=1,idl1
  ch2(ik,1) = c2(ik,1)
END DO
DO  j=2,ip
  DO  k=1,l1
    ch(1,k,j) = c1(1,k,j)
  END DO
END DO
IF (nbd > l1) GO TO 107
is = -ido
DO  j=2,ip
  is = is+ido
  idij = is
  DO  i=3,ido,2
    idij = idij+2
    DO  k=1,l1
      ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
      ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
    END DO
  END DO
END DO
GO TO 111
107 is = -ido
DO  j=2,ip
  is = is+ido
  DO  k=1,l1
    idij = is
    DO  i=3,ido,2
      idij = idij+2
      ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
      ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
    END DO
  END DO
END DO
111 IF (nbd < l1) GO TO 115
DO  j=2,ipph
  jc = ipp2-j
  DO  k=1,l1
    DO  i=3,ido,2
      c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
      c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
      c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
      c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
    END DO
  END DO
END DO
GO TO 121
115 DO  j=2,ipph
  jc = ipp2-j
  DO  i=3,ido,2
    DO  k=1,l1
      c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
      c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
      c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
      c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
    END DO
  END DO
END DO
GO TO 121
119 DO  ik=1,idl1
  c2(ik,1) = ch2(ik,1)
END DO
121 DO  j=2,ipph
  jc = ipp2-j
  DO  k=1,l1
    c1(1,k,j) = ch(1,k,j)+ch(1,k,jc)
    c1(1,k,jc) = ch(1,k,jc)-ch(1,k,j)
  END DO
END DO

ar1 = 1.
ai1 = 0.
DO  l=2,ipph
  lc = ipp2-l
  ar1h = dcp*ar1-dsp*ai1
  ai1 = dcp*ai1+dsp*ar1
  ar1 = ar1h
  DO  ik=1,idl1
    ch2(ik,l) = c2(ik,1)+ar1*c2(ik,2)
    ch2(ik,lc) = ai1*c2(ik,ip)
  END DO
  dc2 = ar1
  ds2 = ai1
  ar2 = ar1
  ai2 = ai1
  DO  j=3,ipph
    jc = ipp2-j
    ar2h = dc2*ar2-ds2*ai2
    ai2 = dc2*ai2+ds2*ar2
    ar2 = ar2h
    DO  ik=1,idl1
      ch2(ik,l) = ch2(ik,l)+ar2*c2(ik,j)
      ch2(ik,lc) = ch2(ik,lc)+ai2*c2(ik,jc)
    END DO
  END DO
END DO
DO  j=2,ipph
  DO  ik=1,idl1
    ch2(ik,1) = ch2(ik,1)+c2(ik,j)
  END DO
END DO

IF (ido < l1) GO TO 132
DO  k=1,l1
  DO  i=1,ido
    cc(i,1,k) = ch(i,k,1)
  END DO
END DO
GO TO 135
132 DO  i=1,ido
  DO  k=1,l1
    cc(i,1,k) = ch(i,k,1)
  END DO
END DO
135 DO  j=2,ipph
  jc = ipp2-j
  j2 = j+j
  DO  k=1,l1
    cc(ido,j2-2,k) = ch(1,k,j)
    cc(1,j2-1,k) = ch(1,k,jc)
  END DO
END DO
IF (ido == 1) RETURN
IF (nbd < l1) GO TO 141
DO  j=2,ipph
  jc = ipp2-j
  j2 = j+j
  DO  k=1,l1
    DO  i=3,ido,2
      ic = idp2-i
      cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
      cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
      cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
      cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
    END DO
  END DO
END DO
RETURN
141 DO  j=2,ipph
  jc = ipp2-j
  j2 = j+j
  DO  i=3,ido,2
    ic = idp2-i
    DO  k=1,l1
      cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
      cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
      cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
      cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
    END DO
  END DO
END DO
RETURN
END SUBROUTINE radfg

SUBROUTINE rfftf (n,r,wsave)


INTEGER, INTENT(IN)                      :: n
REAL(8), INTENT(IN OUT)                     :: r(1)
REAL(8), INTENT(IN OUT)                     :: wsave(1)

IF (n == 1) RETURN
CALL rfftf1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
RETURN
END SUBROUTINE rfftf

SUBROUTINE rfftf1 (n,c,ch,wa,ifac)


INTEGER, INTENT(IN)                      :: n
REAL(8), INTENT(OUT)                        :: c(1)
REAL(8), INTENT(IN)                         :: ch(1)
REAL(8), INTENT(IN OUT)                     :: wa(1)
INTEGER, INTENT(IN)                      :: ifac(10)

nf = ifac(2)
na = 1
l2 = n
iw = n
DO  k1=1,nf
  kh = nf-k1
  ip = ifac(kh+3)
  l1 = l2/ip
  ido = n/l2
  idl1 = ido*l1
  iw = iw-(ip-1)*ido
  na = 1-na
  IF (ip /= 4) GO TO 102
  ix2 = iw+ido
  ix3 = ix2+ido
  IF (na /= 0) GO TO 101
  CALL radf4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
  GO TO 110
  101    CALL radf4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
  GO TO 110
  102    IF (ip /= 2) GO TO 104
  IF (na /= 0) GO TO 103
  CALL radf2 (ido,l1,c,ch,wa(iw))
  GO TO 110
  103    CALL radf2 (ido,l1,ch,c,wa(iw))
  GO TO 110
  104    IF (ip /= 3) GO TO 106
  ix2 = iw+ido
  IF (na /= 0) GO TO 105
  CALL radf3 (ido,l1,c,ch,wa(iw),wa(ix2))
  GO TO 110
  105    CALL radf3 (ido,l1,ch,c,wa(iw),wa(ix2))
  GO TO 110
  106    IF (ip /= 5) GO TO 108
  ix2 = iw+ido
  ix3 = ix2+ido
  ix4 = ix3+ido
  IF (na /= 0) GO TO 107
  CALL radf5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  GO TO 110
  107    CALL radf5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  GO TO 110
  108    IF (ido == 1) na = 1-na
  IF (na /= 0) GO TO 109
  CALL radfg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
  na = 1
  GO TO 110
  109    CALL radfg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  na = 0
  110    l2 = l1
END DO
IF (na == 1) RETURN
DO  i=1,n
  c(i) = ch(i)
END DO
RETURN
END SUBROUTINE rfftf1

SUBROUTINE cosqi (n,wsave)

INTEGER, INTENT(IN)                      :: n
REAL(8), INTENT(OUT)                        :: wsave(1)

DATA pih /1.57079632679491/

dt = pih/FLOAT(n)
fk = 0.
DO  k=1,n
  fk = fk+1.
  wsave(k) = COS(fk*dt)
END DO
CALL rffti (n,wsave(n+1))
RETURN
END SUBROUTINE cosqi

SUBROUTINE cost (n,x,wsave)


INTEGER, INTENT(IN)                      :: n
REAL(8), INTENT(IN OUT)                     :: x(10)
REAL(8), INTENT(IN)                         :: wsave(1)

nm1 = n-1
np1 = n+1
ns2 = n/2
IF (n-2 < 0) THEN
  GO TO   106
ELSE IF (n-2 == 0) THEN
  GO TO   101
ELSE
  GO TO   102
END IF
101 x1h = x(1)+x(2)
x(2) = x(1)-x(2)
x(1) = x1h
RETURN
102 IF (n > 3) GO TO 103
x1p3 = x(1)+x(3)
tx2 = x(2)+x(2)
x(2) = x(1)-x(3)
x(1) = x1p3+tx2
x(3) = x1p3-tx2
RETURN
103 c1 = x(1)-x(n)
x(1) = x(1)+x(n)
DO  k=2,ns2
  kc = np1-k
  t1 = x(k)+x(kc)
  t2 = x(k)-x(kc)
  c1 = c1+wsave(kc)*t2
  t2 = wsave(k)*t2
  x(k) = t1-t2
  x(kc) = t1+t2
END DO
modn = MOD(n,2)
IF (modn /= 0) x(ns2+1) = x(ns2+1)+x(ns2+1)
CALL rfftf (nm1,x,wsave(n+1))
xim2 = x(2)
x(2) = c1
DO  i=4,n,2
  xi = x(i)
  x(i) = x(i-2)-x(i-1)
  x(i-1) = xim2
  xim2 = xi
END DO
IF (modn /= 0) x(n) = xim2
106 RETURN
END SUBROUTINE cost






SUBROUTINE costi (n,wsave)

INTEGER, INTENT(IN)                      :: n
REAL(8), INTENT(OUT)                        :: wsave(1)

DATA pi /3.14159265358979/

IF (n <= 3) RETURN
nm1 = n-1
np1 = n+1
ns2 = n/2
dt = pi/FLOAT(nm1)
fk = 0.
DO  k=2,ns2
  kc = np1-k
  fk = fk+1.
  wsave(k) = 2.*SIN(fk*dt)
  wsave(kc) = 2.*COS(fk*dt)
END DO
CALL rffti (nm1,wsave(n+1))
RETURN
END SUBROUTINE costi

SUBROUTINE rffti (n,wsave)


INTEGER, INTENT(IN)                      :: n
REAL(8), INTENT(IN OUT)                     :: wsave(1)

IF (n == 1) RETURN
CALL rffti1 (n,wsave(n+1),wsave(2*n+1))
RETURN
END SUBROUTINE rffti

SUBROUTINE sinqb (n,x,wsave)


INTEGER, INTENT(IN)                      :: n
REAL(8), INTENT(OUT)                        :: x(1)
REAL(8), INTENT(IN OUT)                     :: wsave(1)

IF (n > 1) GO TO 101
x(1) = 4.*x(1)
RETURN
101 ns2 = n/2
DO  k=2,n,2
  x(k) = -x(k)
END DO
CALL cosqb (n,x,wsave)
DO  k=1,ns2
  kc = n-k
  xhold = x(k)
  x(k) = x(kc+1)
  x(kc+1) = xhold
END DO
RETURN
END SUBROUTINE sinqb

SUBROUTINE sinqf (n,x,wsave)


INTEGER, INTENT(IN)                      :: n
REAL(8), INTENT(IN OUT)                     :: x(1)
REAL(8), INTENT(IN OUT)                     :: wsave(1)

IF (n == 1) RETURN
ns2 = n/2
DO  k=1,ns2
  kc = n-k
  xhold = x(k)
  x(k) = x(kc+1)
  x(kc+1) = xhold
END DO
CALL cosqf (n,x,wsave)
DO  k=2,n,2
  x(k) = -x(k)
END DO
RETURN
END SUBROUTINE sinqf

SUBROUTINE sinqi (n,wsave)


INTEGER, INTENT(IN OUT)                  :: n
REAL(8), INTENT(IN OUT)                     :: wsave(1)

CALL cosqi (n,wsave)
RETURN
END SUBROUTINE sinqi

SUBROUTINE sint (n,x,wsave)


INTEGER, INTENT(IN)                      :: n
REAL(8), INTENT(IN OUT)                     :: x(1)
REAL(8), INTENT(IN OUT)                     :: wsave(1)

np1 = n+1
iw1 = n/2+1
iw2 = iw1+np1
iw3 = iw2+np1
CALL sint1(n,x,wsave,wsave(iw1),wsave(iw2),wsave(iw3))
RETURN
END SUBROUTINE sint

SUBROUTINE sint1(n,war,was,xh,x,ifac)

INTEGER, INTENT(IN)                      :: n
REAL(8), INTENT(IN OUT)                     :: war(1)
REAL(8), INTENT(IN)                         :: was(1)
REAL(8), INTENT(OUT)                        :: xh(10)
REAL(8), INTENT(IN OUT)                     :: x(1)
INTEGER, INTENT(IN OUT)                  :: ifac(1)

DATA sqrt3 /1.73205080756888/

DO  i=1,n
  xh(i) = war(i)
  war(i) = x(i)
END DO
IF (n-2 < 0) THEN
  GO TO   101
ELSE IF (n-2 == 0) THEN
  GO TO   102
ELSE
  GO TO   103
END IF
101 xh(1) = xh(1)+xh(1)
GO TO 106
102 xhold = sqrt3*(xh(1)+xh(2))
xh(2) = sqrt3*(xh(1)-xh(2))
xh(1) = xhold
GO TO 106
103 np1 = n+1
ns2 = n/2
x(1) = 0.
DO  k=1,ns2
  kc = np1-k
  t1 = xh(k)-xh(kc)
  t2 = was(k)*(xh(k)+xh(kc))
  x(k+1) = t1+t2
  x(kc+1) = t2-t1
END DO
modn = MOD(n,2)
IF (modn /= 0) x(ns2+2) = 4.*xh(ns2+1)
CALL rfftf1 (np1,x,xh,war,ifac)
xh(1) = .5*x(1)
DO  i=3,n,2
  xh(i-1) = -x(i)
  xh(i) = xh(i-2)+x(i-1)
END DO
IF (modn /= 0) GO TO 106
xh(n) = -x(n+1)
106 DO  i=1,n
  x(i) = war(i)
  war(i) = xh(i)
END DO
RETURN
END SUBROUTINE sint1

SUBROUTINE sinti (n,wsave)

INTEGER, INTENT(IN)                      :: n
REAL(8), INTENT(OUT)                        :: wsave(1)

DATA pi /3.14159265358979/

IF (n <= 1) RETURN
ns2 = n/2
np1 = n+1
dt = pi/FLOAT(np1)
DO  k=1,ns2
  wsave(k) = 2.*SIN(k*dt)
END DO
CALL rffti (np1,wsave(ns2+1))
RETURN
END SUBROUTINE sinti

SUBROUTINE rffti1 (n,wa,ifac)

INTEGER, INTENT(IN)                      :: n
REAL(8), INTENT(OUT)                        :: wa(1)
INTEGER, INTENT(OUT)                     :: ifac(10)
INTEGER :: ntryh(4)
DATA ntryh(1),ntryh(2),ntryh(3),ntryh(4)/4,2,3,5/

nl = n
nf = 0
j = 0
101 j = j+1
IF (j-4 > 0) THEN
  GO TO   103
END IF
102 ntry = ntryh(j)
GO TO 104
103 ntry = ntry+2
104 nq = nl/ntry
nr = nl-ntry*nq
IF (nr == 0) THEN
  GO TO   105
ELSE
  GO TO   101
END IF
105 nf = nf+1
ifac(nf+2) = ntry
nl = nq
IF (ntry /= 2) GO TO 107
IF (nf == 1) GO TO 107
DO  i=2,nf
  ib = nf-i+2
  ifac(ib+2) = ifac(ib+1)
END DO
ifac(3) = 2
107 IF (nl /= 1) GO TO 104
ifac(1) = n
ifac(2) = nf
tpi = 6.28318530717959
argh = tpi/FLOAT(n)
is = 0
nfm1 = nf-1
l1 = 1
IF (nfm1 == 0) RETURN
DO  k1=1,nfm1
  ip = ifac(k1+2)
  ld = 0
  l2 = l1*ip
  ido = n/l2
  ipm = ip-1
  DO  j=1,ipm
    ld = ld+l1
    i = is
    argld = FLOAT(ld)*argh
    fi = 0.
    DO  ii=3,ido,2
      i = i+2
      fi = fi+1.
      arg = fi*argld
      wa(i-1) = COS(arg)
      wa(i) = SIN(arg)
    END DO
    is = is+ido
  END DO
  l1 = l2
END DO
RETURN
END SUBROUTINE rffti1






