!This file is a part of PW-TELEMAN project.
!PW-TELEMAN is a Time-Dependent Electronic Dynamics in Molecules And Nanosystems library.
!Copyright (C) 2011-2015  Paul-Gerhard Reinhard, Eric Suraud, Florent Calvayrac,
!Phuong Mai Dinh, David Brusson, Philipp Wopperer, José María Escartín Esteban.
!
!PW-Teleman is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!PW-Teleman is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with PW-Teleman.  If not, see <http://www.gnu.org/licenses/>.


! Code converted using TO_F90 by Alan Miller
! Date: 2010-04-16  Time: 12:36:26

!References
!
!              P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!                 Computations (G. Rodrigue, ed.), Academic Press,
!                 1982, pp. 51-83.

SUBROUTINE cosqb (n,x,wsave)
!  Subroutine COSQB computes the fast Fourier transform of quarter
!  wave data. That is, COSQB computes a sequence from its
!  representation in terms of a cosine series with odd wave numbers.
!  The transform is defined below at output parameter X.
!
!  COSQB is the unnormalized inverse of COSQF since a call of COSQB
!  followed by a call of COSQF will multiply the input sequence X
!  by 4*N.
!
!  The array WSAVE which is used by subroutine COSQB must be
!  initialized by calling subroutine COSQI(N,WSAVE).
!
!  Input Parameters
!
!  N       the length of the array X to be transformed.  The method
!          is most efficient when N is a product of small primes.
!
!  X       an array which contains the sequence to be transformed
!
!  WSAVE   a work array which must be dimensioned at least 3*N+15
!          in the program that calls COSQB.  The WSAVE array must be
!          initialized by calling subroutine COSQI(N,WSAVE), and a
!          different WSAVE array must be used for each different
!          value of N.  This initialization does not have to be
!          repeated so long as N remains unchanged.  Thus subsequent
!          transforms can be obtained faster than the first.
!
!  Output Parameters
!
!  X       For I=1,...,N
!
!               X(I)= the sum from K=1 to K=N of
!
!                  2*X(K)*COS((2*K-1)*(I-1)*PI/(2*N))
!
!               A call of COSQB followed by a call of
!               COSQF will multiply the sequence X by 4*N.
!               Therefore COSQF is the unnormalized inverse
!               of COSQB.
!
!  WSAVE   contains initialization calculations which must not
!          be destroyed between calls of COSQB or COSQF.
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: n
REAL(DP), INTENT(OUT)                        :: x(10)
REAL(DP), INTENT(IN OUT)                     :: wsave(1)

DATA tsqrt2 /2.82842712474619D0/

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
!  Subroutine COSQB1 computes the fast Fourier transform of quarter
!  wave data. That is, COSQB1 computes a sequence from its
!  representation in terms of a cosine series with odd wave numbers.
!  The transform is defined below at output parameter X.
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: n
REAL(DP), INTENT(IN OUT)                     :: x(1)
REAL(DP), INTENT(IN)                         :: w(1)
REAL(DP), INTENT(OUT)                        :: xh(1)

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
!  Calculate the fast Fourier transform of subvectors of
!  length two.
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: ido
INTEGER, INTENT(IN)                      :: l1
REAL(DP), INTENT(IN)                         :: cc(ido,2,l1)
REAL(DP), INTENT(OUT)                        :: ch(ido,l1,2)
REAL(DP), INTENT(IN)                         :: wa1(1)

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
!  Calculate the fast Fourier transform of subvectors of
!  length three.
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: ido
INTEGER, INTENT(IN)                      :: l1
REAL(DP), INTENT(IN)                         :: cc(ido,3,l1)
REAL(DP), INTENT(OUT)                        :: ch(ido,l1,3)
REAL(DP), INTENT(IN)                         :: wa1(1)
REAL(DP), INTENT(IN)                         :: wa2(1)

DATA taur,taui /-.5D0,.866025403784439D0/

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
!  Calculate the fast Fourier transform of subvectors of
!  length four.
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: ido
INTEGER, INTENT(IN)                      :: l1
REAL(DP), INTENT(IN)                         :: cc(ido,4,l1)
REAL(DP), INTENT(OUT)                        :: ch(ido,l1,4)
REAL(DP), INTENT(IN)                         :: wa1(1)
REAL(DP), INTENT(IN)                         :: wa2(1)
REAL(DP), INTENT(IN)                         :: wa3(1)

DATA sqrt2 /1.414213562373095D0/

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
!  Calculate the fast Fourier transform of subvectors of
!  length five.
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: ido
INTEGER, INTENT(IN)                      :: l1
REAL(DP), INTENT(IN)                         :: cc(ido,5,l1)
REAL(DP), INTENT(OUT)                        :: ch(ido,l1,5)
REAL(DP), INTENT(IN)                         :: wa1(1)
REAL(DP), INTENT(IN)                         :: wa2(1)
REAL(DP), INTENT(IN)                         :: wa3(1)
REAL(DP), INTENT(IN)                         :: wa4(1)

DATA tr11,ti11,tr12,ti12 /.309016994374947D0,.951056516295154D0,  &
    -.809016994374947D0,.587785252292473D0/

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
!  Calculate the fast Fourier transform of subvectors of
!  arbitrary length.
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: ido
INTEGER, INTENT(IN)                      :: ip
INTEGER, INTENT(IN)                      :: l1
INTEGER, INTENT(IN)                      :: idl1
REAL(DP), INTENT(IN)                         :: cc(ido,ip,l1)
REAL(DP), INTENT(IN OUT)                     :: c1(ido,l1,ip)
REAL(DP), INTENT(OUT)                        :: c2(idl1,ip)
REAL(DP), INTENT(OUT)                        :: ch(ido,l1,ip)
REAL(DP), INTENT(IN OUT)                     :: ch2(idl1,ip)
REAL(DP), INTENT(IN)                         :: wa(1)

DATA tpi/6.28318530717959D0/

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
ai1 = 0D0
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
!Calls for rfftb1
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: n
REAL(DP), INTENT(IN OUT)                     :: r(1)
REAL(DP), INTENT(IN OUT)                     :: wsave(1)

IF (n == 1) RETURN
CALL rfftb1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
RETURN
END SUBROUTINE rfftb

SUBROUTINE rfftb1 (n,c,ch,wa,ifac)
! Called by rfftb
!   Compute the backward fast Fourier transform of a real
!   coefficient array
!   Subroutine RFFTB computes the real periodic sequence from its
!   Fourier coefficients (Fourier synthesis).  The transform is defined
!   below at output parameter R.
!
!   Input Arguments
!
!   N       the length of the array R to be transformed.  The method
!           is most efficient when N is a product of small primes.
!           N may change so long as different work arrays are provided.
!
!   R       a real array of length N which contains the sequence
!           to be transformed.
!
!   WSAVE   a work array which must be dimensioned at least 2*N+15
!           in the program that calls RFFTB.  The WSAVE array must be
!           initialized by calling subroutine RFFTI, and a different
!           WSAVE array must be used for each different value of N.
!           This initialization does not have to be repeated so long as
!           remains unchanged.  Thus subsequent transforms can be
!           obtained faster than the first.  Moreover, the same WSAVE
!           array can be used by RFFTF and RFFTB as long as N remains
!           unchanged.
!
!   Output Argument
!
!   R       For N even and for I = 1,...,N
!
!                R(I) = R(1)+(-1)**(I-1)*R(N)
!
!                     plus the sum from K=2 to K=N/2 of
!
!                      2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
!
!                     -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
!
!           For N odd and for I = 1,...,N
!
!                R(I) = R(1) plus the sum from K=2 to K=(N+1)/2 of
!
!                     2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
!
!                    -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
!
!   Note:  This transform is unnormalized since a call of RFFTF
!          followed by a call of RFFTB will multiply the input
!          sequence by N.
!
!   WSAVE  contains results which must not be destroyed between
!          calls of RFFTB or RFFTF.
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: n
REAL(DP), INTENT(OUT)                        :: c(1)
REAL(DP), INTENT(IN)                         :: ch(1)
REAL(DP), INTENT(IN OUT)                     :: wa(1)
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
!  Subroutine COSQF computes the fast Fourier transform of quarter
!  wave data. That is, COSQF computes the coefficients in a cosine
!  series representation with only odd wave numbers.  The transform
!  is defined below at Output Parameter X
!
!  COSQF is the unnormalized inverse of COSQB since a call of COSQF
!  followed by a call of COSQB will multiply the input sequence X
!  by 4*N.
!
!  The array WSAVE which is used by subroutine COSQF must be
!  initialized by calling subroutine COSQI(N,WSAVE).
!
!  Input Parameters
!
!  N       the length of the array X to be transformed.  The method
!          is most efficient when N is a product of small primes.
!
!  X       an array which contains the sequence to be transformed
!
!  WSAVE   a work array which must be dimensioned at least 3*N+15
!          in the program that calls COSQF.  The WSAVE array must be
!          initialized by calling subroutine COSQI(N,WSAVE), and a
!          different WSAVE array must be used for each different
!          value of N.  This initialization does not have to be
!          repeated so long as N remains unchanged.  Thus subsequent
!          transforms can be obtained faster than the first.
!
!  Output Parameters
!
!  X       For I=1,...,N
!
!               X(I) = X(1) plus the sum from K=2 to K=N of
!
!                  2*X(K)*COS((2*I-1)*(K-1)*PI/(2*N))
!
!               A call of COSQF followed by a call of
!               COSQB will multiply the sequence X by 4*N.
!               Therefore COSQB is the unnormalized inverse
!               of COSQF.
!
!  WSAVE   contains initialization calculations which must not
!          be destroyed between calls of COSQF or COSQB.
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: n
REAL(DP), INTENT(IN OUT)                     :: x(10)
REAL(DP), INTENT(IN OUT)                     :: wsave(1)

DATA sqrt2 /1.4142135623731D0/

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
!  Subroutine COSQF1 computes the fast Fourier transform of quarter
!  wave data. That is, COSQF1 computes the coefficients in a cosine
!  series representation with only odd wave numbers.  The transform
!  is defined below at Output Parameter X
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: n
REAL(DP), INTENT(IN OUT)                     :: x(1)
REAL(DP), INTENT(IN)                         :: w(1)
REAL(DP), INTENT(OUT)                        :: xh(1)

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
!  Calculate the fast Fourier transform of subvectors of
!  length two.
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: ido
INTEGER, INTENT(IN)                      :: l1
REAL(DP), INTENT(IN)                         :: cc(ido,l1,2)
REAL(DP), INTENT(OUT)                        :: ch(ido,2,l1)
REAL(DP), INTENT(IN)                         :: wa1(1)

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
!  Calculate the fast Fourier transform of subvectors of
!  length three.
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: ido
INTEGER, INTENT(IN)                      :: l1
REAL(DP), INTENT(IN)                         :: cc(ido,l1,3)
REAL(DP), INTENT(OUT)                        :: ch(ido,3,l1)
REAL(DP), INTENT(IN)                         :: wa1(1)
REAL(DP), INTENT(IN)                         :: wa2(1)

DATA taur,taui /-.5D0,.866025403784439D0/

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
!  Calculate the fast Fourier transform of subvectors of
!  length four.
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: ido
INTEGER, INTENT(IN)                      :: l1
REAL(DP), INTENT(IN)                         :: cc(ido,l1,4)
REAL(DP), INTENT(OUT)                        :: ch(ido,4,l1)
REAL(DP), INTENT(IN)                         :: wa1(1)
REAL(DP), INTENT(IN)                         :: wa2(1)
REAL(DP), INTENT(IN)                         :: wa3(1)

DATA hsqt2 /.7071067811865475D0/

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
!  Calculate the fast Fourier transform of subvectors of
!  length five.
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: ido
INTEGER, INTENT(IN)                      :: l1
REAL(DP), INTENT(IN)                         :: cc(ido,l1,5)
REAL(DP), INTENT(OUT)                        :: ch(ido,5,l1)
REAL(DP), INTENT(IN)                         :: wa1(1)
REAL(DP), INTENT(IN)                         :: wa2(1)
REAL(DP), INTENT(IN)                         :: wa3(1)
REAL(DP), INTENT(IN)                         :: wa4(1)

DATA tr11,ti11,tr12,ti12 /.309016994374947D0,.951056516295154D0,  &
    -.809016994374947D0,.587785252292473D0/

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
!  Calculate the fast Fourier transform of subvectors of
!  arbitrary length.
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: ido
INTEGER, INTENT(IN)                      :: ip
INTEGER, INTENT(IN)                      :: l1
INTEGER, INTENT(IN)                      :: idl1
REAL(DP), INTENT(OUT)                        :: cc(ido,ip,l1)
REAL(DP), INTENT(IN OUT)                     :: c1(ido,l1,ip)
REAL(DP), INTENT(IN OUT)                     :: c2(idl1,ip)
REAL(DP), INTENT(OUT)                        :: ch(ido,l1,ip)
REAL(DP), INTENT(OUT)                        :: ch2(idl1,ip)
REAL(DP), INTENT(IN)                         :: wa(1)

DATA tpi/6.28318530717959D0/

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
! Calls for rfftf1 if n is not 1
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: n
REAL(DP), INTENT(IN OUT)                     :: r(1)
REAL(DP), INTENT(IN OUT)                     :: wsave(1)

IF (n == 1) RETURN
CALL rfftf1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
RETURN
END SUBROUTINE rfftf

SUBROUTINE rfftf1 (n,c,ch,wa,ifac)
!   Subroutine RFFTF1 computes the Fourier coefficients of a real
!   periodic sequence (Fourier analysis).  The transform is defined
!   below at output parameter C.
!
!   The arrays WA and IFAC which are used by subroutine RFFTB1 must be
!   initialized by calling subroutine RFFTI1.
!
!   Input Arguments
!
!   N       the length of the array R to be transformed.  The method
!           is most efficient when N is a product of small primes.
!           N may change so long as different work arrays are provided.
!
!   C       a real array of length N which contains the sequence
!           to be transformed.
!
!   CH      a real work array of length at least N.
!
!   WA      a real work array which must be dimensioned at least N.
!
!   IFAC    an integer work array which must be dimensioned at least 15.
!
!           The WA and IFAC arrays must be initialized by calling
!           subroutine RFFTI1, and different WA and IFAC arrays must be
!           used for each different value of N.  This initialization
!           does not have to be repeated so long as N remains unchanged.
!           Thus subsequent transforms can be obtained faster than the
!           first.  The same WA and IFAC arrays can be used by RFFTF1
!           and RFFTB1.
!
!   Output Argument
!
!   C       C(1) = the sum from I=1 to I=N of R(I)
!
!           If N is even set L = N/2; if N is odd set L = (N+1)/2
!
!             then for K = 2,...,L
!
!                C(2*K-2) = the sum from I = 1 to I = N of
!
!                     C(I)*COS((K-1)*(I-1)*2*PI/N)
!
!                C(2*K-1) = the sum from I = 1 to I = N of
!
!                    -C(I)*SIN((K-1)*(I-1)*2*PI/N)
!
!           If N is even
!
!                C(N) = the sum from I = 1 to I = N of
!
!                     (-1)**(I-1)*C(I)
!
!   Notes:  This transform is unnormalized since a call of RFFTF1
!           followed by a call of RFFTB1 will multiply the input
!           sequence by N.
!
!           WA and IFAC contain initialization calculations which must
!           not be destroyed between calls of subroutine RFFTF1 or
!           RFFTB1.
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: n
REAL(DP), INTENT(OUT)                        :: c(1)
REAL(DP), INTENT(IN)                         :: ch(1)
REAL(DP), INTENT(IN OUT)                     :: wa(1)
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
!  Subroutine COSQI initializes the work array WSAVE which is used in
!  both COSQF1 and COSQB1.  The prime factorization of N together with
!  a tabulation of the trigonometric functions are computed and
!  stored in WSAVE.
!
!  Input Parameter
!
!  N       the length of the array to be transformed.  The method
!          is most efficient when N is a product of small primes.
!
!  Output Parameter
!
!  WSAVE   a work array which must be dimensioned at least 3*N+15.
!          The same work array can be used for both COSQF1 and COSQB1
!          as long as N remains unchanged.  Different WSAVE arrays
!          are required for different values of N.  The contents of
!          WSAVE must not be changed between calls of COSQF1 or COSQB1.
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: n
REAL(DP), INTENT(OUT)                        :: wsave(1)

DATA pih /1.57079632679491D0/

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
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

!  Subroutine COST computes the discrete Fourier cosine transform
!  of an even sequence X(I).  The transform is defined below at output
!  parameter X.
!
!  COST is the unnormalized inverse of itself since a call of COST
!  followed by another call of COST will multiply the input sequence
!  X by 2*(N-1).  The transform is defined below at output parameter X.
!
!  The array WSAVE which is used by subroutine COST must be
!  initialized by calling subroutine COSTI(N,WSAVE).
!
!  Input Parameters
!
!  N       the length of the sequence X.  N must be greater than 1.
!          The method is most efficient when N-1 is a product of
!          small primes.
!
!  X       an array which contains the sequence to be transformed
!
!  WSAVE   a work array which must be dimensioned at least 3*N+15
!          in the program that calls COST.  The WSAVE array must be
!          initialized by calling subroutine COSTI(N,WSAVE), and a
!          different WSAVE array must be used for each different
!          value of N.  This initialization does not have to be
!          repeated so long as N remains unchanged.  Thus subsequent
!          transforms can be obtained faster than the first.
!
!  Output Parameters
!
!  X       For I=1,...,N
!
!             X(I) = X(1)+(-1)**(I-1)*X(N)
!
!               + the sum from K=2 to K=N-1
!
!                 2*X(K)*COS((K-1)*(I-1)*PI/(N-1))
!
!               A call of COST followed by another call of
!               COST will multiply the sequence X by 2*(N-1).
!               Hence COST is the unnormalized inverse
!               of itself.
!
!  WSAVE   contains initialization calculations which must not be
!          destroyed between calls of COST.
!

INTEGER, INTENT(IN)                      :: n
REAL(DP), INTENT(IN OUT)                     :: x(10)
REAL(DP), INTENT(IN)                         :: wsave(1)

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
!  Subroutine COSTI initializes the array WSAVE which is used in
!  subroutine COST.  The prime factorization of N together with
!  a tabulation of the trigonometric functions are computed and
!  stored in WSAVE.
!
!  Input Parameter
!
!  N       the length of the sequence to be transformed.  The method
!          is most efficient when N-1 is a product of small primes.
!
!  Output Parameter
!
!  WSAVE   a work array which must be dimensioned at least 3*N+15.
!          Different WSAVE arrays are required for different values
!          of N.  The contents of WSAVE must not be changed between
!          calls of COST.
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: n
REAL(DP), INTENT(OUT)                        :: wsave(1)

DATA pi /3.14159265358979D0/

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
! Calls rffti1 if n is not 1 
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: n
REAL(DP), INTENT(IN OUT)                     :: wsave(1)

IF (n == 1) RETURN
CALL rffti1 (n,wsave(n+1),wsave(2*n+1))
RETURN
END SUBROUTINE rffti

SUBROUTINE sinqb (n,x,wsave)
!  Subroutine SINQB computes the fast Fourier transform of quarter
!  wave data.  That is, SINQB computes a sequence from its
!  representation in terms of a sine series with odd wave numbers.
!  the transform is defined below at output parameter X.
!
!  SINQF is the unnormalized inverse of SINQB since a call of SINQB
!  followed by a call of SINQF will multiply the input sequence X
!  by 4*N.
!
!  The array WSAVE which is used by subroutine SINQB must be
!  initialized by calling subroutine SINQI(N,WSAVE).
!
!  Input Parameters
!
!  N       the length of the array X to be transformed.  The method
!          is most efficient when N is a product of small primes.
!
!  X       an array which contains the sequence to be transformed
!
!  WSAVE   a work array which must be dimensioned at least 3*N+15
!          in the program that calls SINQB.  The WSAVE array must be
!          initialized by calling subroutine SINQI(N,WSAVE), and a
!          different WSAVE array must be used for each different
!          value of N.  This initialization does not have to be
!          repeated so long as N remains unchanged.  Thus subsequent
!          transforms can be obtained faster than the first.
!
!  Output Parameters
!
!  X       For I=1,...,N
!
!               X(I)= the sum from K=1 to K=N of
!
!                 4*X(K)*SIN((2*K-1)*I*PI/(2*N))
!
!               a call of SINQB followed by a call of
!               SINQF will multiply the sequence X by 4*N.
!               Therefore SINQF is the unnormalized inverse
!               of SINQB.
!
!  WSAVE   contains initialization calculations which must not
!          be destroyed between calls of SINQB or SINQF.
!
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: n
REAL(DP), INTENT(OUT)                        :: x(1)
REAL(DP), INTENT(IN OUT)                     :: wsave(1)

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
!  Subroutine SINQF computes the fast Fourier transform of quarter
!  wave data.  That is, SINQF computes the coefficients in a sine
!  series representation with only odd wave numbers.  The transform
!  is defined below at output parameter X.
!
!  SINQB is the unnormalized inverse of SINQF since a call of SINQF
!  followed by a call of SINQB will multiply the input sequence X
!  by 4*N.
!
!  The array WSAVE which is used by subroutine SINQF must be
!  initialized by calling subroutine SINQI(N,WSAVE).
!
!  Input Parameters
!
!  N       the length of the array X to be transformed.  The method
!          is most efficient when N is a product of small primes.
!
!  X       an array which contains the sequence to be transformed
!
!  WSAVE   a work array which must be dimensioned at least 3*N+15
!          in the program that calls SINQF.  The WSAVE array must be
!          initialized by calling subroutine SINQI(N,WSAVE), and a
!          different WSAVE array must be used for each different
!          value of N.  This initialization does not have to be
!          repeated so long as N remains unchanged.  Thus subsequent
!          transforms can be obtained faster than the first.
!
!  Output Parameters
!
!  X       For I=1,...,N
!
!               X(I) = (-1)**(I-1)*X(N)
!
!                  + the sum from K=1 to K=N-1 of
!
!                  2*X(K)*SIN((2*I-1)*K*PI/(2*N))
!
!               A call of SINQF followed by a call of
!               SINQB will multiply the sequence X by 4*N.
!               Therefore SINQB is the unnormalized inverse
!               of SINQF.
!
!  WSAVE   contains initialization calculations which must not
!          be destroyed between calls of SINQF or SINQB.
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: n
REAL(DP), INTENT(IN OUT)                     :: x(1)
REAL(DP), INTENT(IN OUT)                     :: wsave(1)

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
! Initialize a work array for SINQF and SINQB  
! Same as COSQI
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN OUT)                  :: n
REAL(DP), INTENT(IN OUT)                     :: wsave(1)

CALL cosqi (n,wsave)
RETURN
END SUBROUTINE sinqi

SUBROUTINE sint (n,x,wsave)
!  Subroutine SINT computes the discrete Fourier sine transform
!  of an odd sequence X(I).  The transform is defined below at
!  output parameter X.
!
!  SINT is the unnormalized inverse of itself since a call of SINT
!  followed by another call of SINT will multiply the input sequence
!  X by 2*(N+1).
!
!  The array WSAVE which is used by subroutine SINT must be
!  initialized by calling subroutine SINTI(N,WSAVE).
!
!  Input Parameters
!
!  N       the length of the sequence to be transformed.  The method
!          is most efficient when N+1 is the product of small primes.
!
!  X       an array which contains the sequence to be transformed
!
!
!  WSAVE   a work array with dimension at least INT(3.5*N+16)
!          in the program that calls SINT.  The WSAVE array must be
!          initialized by calling subroutine SINTI(N,WSAVE), and a
!          different WSAVE array must be used for each different
!          value of N.  This initialization does not have to be
!          repeated so long as N remains unchanged.  Thus subsequent
!          transforms can be obtained faster than the first.
!
!  Output Parameters
!
!  X       For I=1,...,N
!
!               X(I)= the sum from K=1 to K=N
!
!                    2*X(K)*SIN(K*I*PI/(N+1))
!
!               A call of SINT followed by another call of
!               SINT will multiply the sequence X by 2*(N+1).
!               Hence SINT is the unnormalized inverse
!               of itself.
!
!  WSAVE   contains initialization calculations which must not be
!          destroyed between calls of SINT.
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: n
REAL(DP), INTENT(IN OUT)                     :: x(1)
REAL(DP), INTENT(IN OUT)                     :: wsave(1)

np1 = n+1
iw1 = n/2+1
iw2 = iw1+np1
iw3 = iw2+np1
CALL sint1(n,x,wsave,wsave(iw1),wsave(iw2),wsave(iw3))
RETURN
END SUBROUTINE sint

SUBROUTINE sint1(n,war,was,xh,x,ifac)
! See sint
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: n
REAL(DP), INTENT(IN OUT)                     :: war(1)
REAL(DP), INTENT(IN)                         :: was(1)
REAL(DP), INTENT(OUT)                        :: xh(10)
REAL(DP), INTENT(IN OUT)                     :: x(1)
INTEGER, INTENT(IN OUT)                  :: ifac(1)

DATA sqrt3 /1.73205080756888D0/

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
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: n
REAL(DP), INTENT(OUT)                        :: wsave(1)

DATA pi /3.14159265358979D0/

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
!   Subroutine RFFTI1 initializes the work arrays WA and IFAC which are
!   used in both RFFTF1 and RFFTB1.  The prime factorization of N and a
!   tabulation of the trigonometric functions are computed and stored in
!   IFAC and WA, respectively.
!
!   Input Argument
!
!   N       the length of the sequence to be transformed.
!
!   Output Arguments
!
!   WA      a real work array which must be dimensioned at least N.
!
!   IFAC    an integer work array which must be dimensioned at least 15.
!
!   The same work arrays can be used for both RFFTF1 and RFFTB1 as long
!   as N remains unchanged.  Different WA and IFAC arrays are required
!   for different values of N.  The contents of WA and IFAC must not be
!   changed between calls of RFFTF1 or RFFTB1.
USE params, ONLY:DP
IMPLICIT REAL(DP) (A-H,O-Z)

INTEGER, INTENT(IN)                      :: n
REAL(DP), INTENT(OUT)                        :: wa(1)
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
tpi = 6.28318530717959D0
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






