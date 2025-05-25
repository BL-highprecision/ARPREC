subroutine f_main
! program mathinit

!   This program pre-computes data to be read by the mathtool program.  By
!   default, the parameter nquadopt is set below to 1, so that only data for
!   the tanh-sinh quadrature routine is generated.  Setting nquadopt to 2
!   generates data for two quadrature routines (tanh-sinh and error function). 
!   Setting nquadopt to 3 generates data for all three quadrature routines 
!   (tanh-sinh, error function and Gaussian).  Typical run times on a 2004-era
!   computer (assuming ndigmx1 = 1000, ndigmx2 = 2000 and nquadx = 10):
!     Tanh-sinh        2 minutes
!     Error function   10 minutes
!     Gaussian         100 minutes
!     Total            112 minutes

!   See the comments at the start of globdata.f regarding this program.
!   If the message "table space parameter is too small" appears when this 
!   program is run, the value of nquadz in globdata.f must be increased -- 
!   typically change the "18" in "nquadz = 18 ** nquadx", the formula for 
!   nquadz, to some higher value, such as 25.

!   David H Bailey    2004-07-01

use mpmodule
use globdata
integer i, nq3, nquadopt
parameter (nquadopt = 1)
real*8 second, tm0, tm1, tm2
type (mp_real) catalan, eulergamma, t1, t2
external catalan, eulergamma, second

call mpinit (ndigmx2)
ndebug = 2
write (6, *) 'mathinit: start'
tm0 = second ()

!   Generate several constants.

tm1 = second ()
open (11, file = 'const.dat', form = 'unformatted')
t1 = catalan ()
t2 = eulergamma ()
write (6, *) 'const complete'
write (11) t1, t2
rewind (11)
close (11)
write (6, *) 'file const.dat written'
tm2 = second ()
write (6, *) 'cpu time =', tm2 - tm1

call mpsetprec (ndigmx1)

if (nquadopt >= 1) then

!   Generate data for tanh-sinh quadrature.

  tm1 = second ()
  open (11, file = 'quadts.dat', form = 'unformatted')
  call initqts (nquadx, nquadz, nq3, quadwk, quadxk)
  write (6, *) 'initqts complete'
  write (11) nq3
  write (11) (quadwk(i), i = -1, nq3)
  write (11) (quadxk(i), i = -1, nq3)
  rewind (11)
  close (11)
  write (6, *) 'file quadts.dat written'
  tm2 = second ()
  write (6, *) 'cpu time =', tm2 - tm1
endif

if (nquadopt >= 2) then

!   Generate data for erf quadrature.

  tm1 = second ()
  open (11, file = 'quaderf.dat', form = 'unformatted')
  call initqerf (nquadx, nquadz, nq3, quadwk, quadxk)
  write (6, *) 'initqerf complete'
  write (11) nq3
  write (11) (quadwk(i), i = -1, nq3)
  write (11) (quadxk(i), i = -1, nq3)
  rewind (11)
  close (11)
  write (6, *) 'file quaderf.dat written'
  tm2 = second ()
  write (6, *) 'cpu time =', tm2 - tm1
endif

if (nquadopt >= 3) then

!   Generate data for Gaussian quadrature.

  tm1 = second ()
  open (11, file = 'quadgs.dat', form = 'unformatted')
  call initqgs (nquadx, nquadz, nq3, quadwk, quadxk)
  write (6, *) 'initqt complete'
  write (11) nq3
  write (11) (quadwk(i), i = -1, nq3)
  write (11) (quadxk(i), i = -1, nq3)
  rewind (11)
  close (11)
  write (6, *) 'file quadgs.dat written'
  tm2 = second ()
  write (6, *) 'cpu time =', tm2 - tm1
endif

write (6, *) 'total cpu time =', tm2 - tm0
stop
end

subroutine initqgs (nq1, nq2, nq3, wk, xk)

!   This subroutine initializes the quadrature arays xk and wk for Gaussian
!   quadrature.  It employs a Newton iteration scheme with a dynamic precision
!   level.  The argument nq2, which is the space allocated for wk and xk in
!   the calling program, should be at least 8 * 2^nq1 + 100, although a higher
!   value may be required, depending on precision level.  Monitor the space
!   figure given in the message below during initialization to be certain.
!   Modified for the Toolkit
!   David H Bailey    2004-05-27

use mpmodule
use globdata
implicit none
integer i, ierror, ik0, is, j, j1, k, n, nq1, nq2, nq3, nwp, nws
double precision pi
parameter (pi = 3.141592653589793238d0)
type (mp_real) eps, r, t1, t2, t3, t4, t5, wk(-1:nq2), xk(-1:nq2)
parameter (ik0 = 100)

if (ndebug >= 1) then
  write (6, 1)
1 format ('initqgs: Gaussian quadrature initialization')
endif

call mpgetprecwords (nws)
eps = 2.d0 ** (-96)
wk(-1) = dble (nq1)
wk(0) = 0.d0
xk(0) = 0.d0
wk(1) = dble (nq1)
xk(1) = dble (ik0)
i = ik0

do j = 2, ik0
  wk(j) = 0.d0
  xk(j) = 0.d0
enddo

do k = 1, nq1
  if (ndebug >= 2) write (6, *) k, i, nq2
  n = 3 * 2 ** (k + 1)

  do j = 1, n / 2

!   Compute a double precision estimate of the root.

    is = 0
    nwp = 3
    call mpsetprecwords (nwp)
    r = cos ((pi * (j - 0.25d0)) / (n + 0.5d0))

!   Compute the j-th root of the n-degree Legendre polynomial using Newton's
!   iteration.

100 continue

    t1 = 1.d0
    t2 = 0.d0

    do j1 = 1, n
      t3 = t2
      t2 = t1
      t1 = ((2 * j1 - 1) * r * t2 - (j1 - 1) * t3) / j1
    enddo

    t4 = n * (r * t1 - t2) / (r ** 2  - 1.d0)
    t5 = r
    r = r - t1 / t4

!   Once convergence is achieved at nwp = 3, then start doubling (almost) the
!   working precision level at each iteration until full precision is reached.

    if (nwp == 3) then
      if (abs (r - t5) > eps) goto 100
      nwp = min (2 * nwp - 1, nws)
      call mpsetprecwords (nwp)
      goto 100
    elseif (nwp < nws) then
      nwp = min (2 * nwp - 1, nws)
      call mpsetprecwords (nwp)
      goto 100
    endif

    i = i + 1
    if (i > nq2) goto 110
    xk(i) = r
    t4 = n * (r * t1 - t2) / (r ** 2  - 1.d0)
    wk(i) = 2.d0 / ((1.d0 - r ** 2) * t4 ** 2)
    call mpgetpar ('mpier', ierror)
    if (ierror > 0) goto 120
  enddo

  xk(k+1) = i
enddo

xk(-1) = dble (i)
nq3 = i
if (ndebug >= 2) then
  write (6, 2) i
2 format ('initqerf: Table spaced used =',i8)
endif
goto 130

110 continue

nq3 = 0
write (6, 3) nq2
3 format ('initqgs: Table space parameter is too small; value =',i8)
nerror = 92
goto 130

120 continue

nerror = ierror + 100
write (6, 4) nerror
4 format ('initqgs: Error in quadrature initialization; code =',i5)

130 continue

call mpsetprecwords (nws)
return
end

subroutine initqerf (nq1, nq2, nq3, wk, xk)

!   This subroutine initializes the quadrature arays xk and wk using the
!   function x(t) = erf(t) = 1 - erfc(t).  The argument nq2 is the space
!   allocated for wk and xk in the calling program.  By default it is set to 
!   12 * 2^nq1.  Increase nq2 if directed by a message produced below.
!   Upon completion, wk(-1) = nq1, and xk(-1) = n, the maximum space parameter
!   for these arrays.  In other words, the arrays occupy (wk(i), i = -1 to n)
!   and (xk(i), i = -1 to n), where n = xk(-1).  The array x_k contain 1 minus
!   the abscissas; the wk array are the weights at these abscissas.
!   Modified for the toolkit.

!   David H Bailey    2004-05-27

use mpmodule
use globdata
implicit none
integer i, ierror, iprint, j, k, k1, nq1, nq2, nq3
logical log1
real*8 h
parameter (iprint = 1000)
type (mp_real) eps2, p2, spi, t1, t2, t3, t4, t5, wk(-1:nq2), xk(-1:nq2)

if (ndebug >= 1) then
  write (6, 1)
1 format ('initqerf: Error function quadrature initialization')
endif

eps2 = mpreal (10.d0) ** (-ndigmx2)
p2 = 0.5d0 * mppic
spi = 2.d0 / sqrt (mppic)
h = 0.5d0 ** (nq1 - 2)
wk(-1) = dble (nq1)

do k = 0, nq2
  if (ndebug >= 2 .and. mod (k, iprint) == 0) write (6, *) k, nq2
  t1 = dble (k) * h
  xk(k) = erfc (t1)
  wk(k) = spi * exp (- t1 ** 2)
  call mpgetpar ('mpier', ierror)
  if (ierror > 0) goto 120
  if (wk(k) < eps2) goto 100
enddo

nq3 = 0
write (6, 2) nq2
2 format ('initqerf: Table space parameter is too small; value =',i8)
nerror = 91
goto 130

100 continue

nq3 = k
xk(-1) = dble (k)
if (ndebug >= 2) then
  write (6, 3) k
3 format ('initqerf: Table spaced used =',i8)
endif
goto 130

120 continue

nerror = ierror + 100
write (6, 4) nerror
4 format ('initqerf: Error in quadrature initialization; code =',i5)

130  continue

return
end

subroutine initqts (nq1, nq2, nq3, wk, xk)

!   This subroutine initializes the quadrature arays xk and wk using the
!   function x(t) = tanh (pi/2*sinh(t)).  The argument nq2 is the space
!   allocated for wk and xk in the calling program.  By default it is set to 
!   12 * 2^nq1.  Increase nq2 if directed by a message produced below.
!   Upon completion, wk(-1) = nq1, and xk(-1) = n, the maximum space parameter
!   for these arrays.  In other words, the arrays occupy (wk(i), i = -1 to n)
!   and (xk(i), i = -1 to n), where n = xk(-1).   The array x_k contain 1 minus
!   the abscissas; the wk array are the weights at these abscissas.
!   Modified for the toolkit.
!   David H Bailey    2004-05-27

use mpmodule
use globdata
implicit none
integer i, ierror, iprint, j, k, k1, nq1, nq2, nq3
logical log1
real*8 h
parameter (iprint = 1000)
type (mp_real) eps2, p2, spi, t1, t2, t3, t4, u1, u2, wk(-1:nq2), xk(-1:nq2)

if (ndebug >= 1) then
  write (6, 1)
1 format ('initqts: Tanh-sinh quadrature initialization')
endif

eps2 = mpreal (10.d0) ** (-ndigmx2)
p2 = 0.5d0 * mppic
h = 0.5d0 ** nq1
wk(-1) = dble (nq1)

do k = 0, nq2
  if (ndebug >= 2 .and. mod (k, iprint) == 0) write (6, *) k, nq2
  t1 = dble (k) * h

!   xk(k) = 1 - tanh (u1) = 1 /(e^u1 * cosh (u1))
!   wk(k) = u2 / cosh (u1)^2
!   where u1 = pi/2 * cosh (t1), u2 = pi/2 * sinh (t1)

  t2 = exp (t1)
  u1 = 0.5d0 * p2 * (t2 + 1.d0 / t2)
  u2 = 0.5d0 * p2 * (t2 - 1.d0 / t2)
  t3 = exp (u2)
  t4 = 0.5d0 * (t3 + 1.d0 / t3)
  xk(k) = 1.d0 / (t3 * t4)
  wk(k) = u1 / t4 ** 2
  call mpgetpar ('mpier', ierror)
  if (ierror > 0) goto 120
  if (wk(k) < eps2) goto 100
enddo

nq3 = 0
write (6, 2) nq2
2 format ('initqts: Table space parameter is too small; value =',i8)
nerror = 91
goto 130

100 continue

nq3 = k
xk(-1) = dble (k)
if (ndebug >= 2) then
  write (6, 3) k
3 format ('initqts: Table spaced used =',i8)
endif
goto 130

120 continue

nerror = ierror + 100
write (6, 4) nerror
4 format ('initqts: Error in quadrature initialization; code =',i5)

130  continue

return
end

function catalan ()
use mpmodule
use globdata
implicit none
integer k
real*8 dk
type (mp_real) catalan, c1, c2, c4, c8, r16, t1, t2, t3
type (mp_real) x1, x2, x3, x4, x5, x6, eps

c1 = 1.d0
c2 = 2.d0
c4 = 4.d0
c8 = 8.d0
r16 = 1.d0 / 16.d0
t1 = 0.d0
t2 = 1.d0
eps = mpreal (10.d0) ** (-ndigmx1-1)

do k = 0, 10000000
  dk = k
  t3 = t2 * (c8 / (8.d0 * dk + 1.d0) ** 2 + c8 / (8.d0 * dk + 2.d0) ** 2 &
       + c4 / (8.d0 * dk + 3.d0) ** 2 - c2 / (8.d0 * dk + 5.d0) ** 2 &
       - c2 / (8.d0 * dk + 6.d0) ** 2 - c1 / (8.d0 * dk + 7.d0) ** 2)
  t1 = t1 + t3
  t2 = r16 * t2
  if (t3 < eps) goto 100
enddo

write (6, *) 'catalan: error - contact author'

100 continue

catalan = 1.d0 / 8.d0 * mppic * mpl02 + 1.d0 / 16.d0 * t1
return
end

function eulergamma ()
use mpmodule
use globdata
implicit none
integer k, m, n
real*8 d1, dplog10q
type (mp_real) eulergamma, c1, t1, t2, t3, t4, t5
external dplog10q

n = log (ndigmx1 * log (10.d0)) / log (2.d0) + 1.d0
d1 = 2.d0 ** n
c1 = 1.d0
t1 = 0.d0
t3 = 1.d0 / d1
t4 = 1.d0
t5 = 1.d0

do m = 0, 1000000000
  t3 = d1 * t3
  t4 = (m + 1) * t4
  t5 = t5 + c1 / (m + 1)
  t2 = t3 / t4 * t5
  t1 = t1 + t2
  if (t2 < 1.d0) goto 100
enddo

write (6, *) 'eulergamma: error - contact author'

100 continue

t5 = d1
eulergamma = d1 / exp (t5) * t1 - n * mpl02 - 1.d0
return
end

function dplog10q (a)

!   For input MP value a, this routine returns a DP approximation to log10 (a).

use mpmodule
implicit none
integer ia
double precision da, dplog10q, t1
type (mp_real) a

call mpmdc (a%mpr, da, ia)
if (da .eq. 0.d0) then
  dplog10q = -9999.d0
else
  dplog10q = log10 (abs (da)) + ia * log10 (2.d0)
endif

100 continue
return
end

