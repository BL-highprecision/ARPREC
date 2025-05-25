function quadgs (kv, str1, x1, x2, tmp1, nq1, nq2, wk, xk)

!   This routine computes the integral of the function fun on the interval
!   [0, 1], with up to nq1 iterations, with a target tolerance of 10^nepsilon1.
!   wk and xk are precomputed tables of weights and abscissas, each of size
!   nq2 and of type mp_real.  The function fun is not evaluated at x1 or x2.
!   Modified for the Toolkit.

!   David H Bailey    2004-07-27
!   David H Bailey    2006-11-03   Modifications made for toolkit usage

use mpmodule
use globdata
implicit none
integer i, ierror, ik0, k, j, n, nds, nq1, nq2, nerr, nwp
double precision d1, d2, d3, d4, dplog10q, dpw
type (mp_real) ax, bx, c10, eps1, epsilon1, epsilon2, err, fun, &
  quadgs, tsum, s1, s2, s3, t1, t2, t3, tw1, tw2, twmx, &
  wk(-1:nq2), xk(-1:nq2), x1, x2, xx1, xx2
integer kv, ntmp1
character*2048 str1
character*1 stout(72)
type (mp_real) tmp1(ntmpx)
external dplog10q
parameter (ik0 = 100, dpw = 14.44943979d0)

! nds = mpoud
! mpoud = 56
ax = 0.5d0 * (x2 - x1)
bx = 0.5d0 * (x2 + x1)
epsilon1 = mpreal (10.d0) ** nepsilon1
epsilon2 = mpreal (10.d0) ** nepsilon2
s1 = 0.d0
s2 = 0.d0
c10 = 10.d0

if (wk(-1) < dble (nq1)) then
  write (6, 1) nq1
1 format ('quadgs: quadrature arrays have not been initialized; nq1 =',i6)
  nerror = 70
  goto 130
endif

nwp = 8
call mpsetprecwords (nwp)

do k = 1, nq1
  n = 3 * 2 ** (k + 1)
  s3 = s2
  s2 = s1

100 continue

  twmx = 0.d0
  tsum = 0.d0
  i = dble (xk(k))

  do j = 1, n / 2
    i = i + 1
    xx1 = - ax * xk(i) + bx
    xx2 = ax * xk(i) + bx

    if (xx1 > x1) then
      var(kv) = xx1
      call parse (str1, ntmp1, tmp1(1))  
      t1 = tmp1(1)
      call mpgetpar ('mpier', ierror)
      if (ierror > 0 .or. nerror > 0) goto 120
      tw1 = t1 * wk(i)
    else
      t1 = 0.d0
      tw1 = 0.d0
    endif

    if (xx2 < x2 .and. j + k > 2) then
      var(kv) = xx2
      call parse (str1, ntmp1, tmp1(1)) 
      t2 = tmp1(1)  
      call mpgetpar ('mpier', ierror)
      if (ierror > 0 .or. nerror > 0) goto 120
      tw2 = t2 * wk(i)
    else
      t2 = 0.d0
      tw2 = 0.d0
    endif

    tsum = tsum + tw1 + tw2
    twmx = max (twmx, abs (tw1), abs (tw2))
  enddo

!   Compute s1 = current integral approximation and err = error estimate.
!   Tsum is the sum of all tw1 and tw2 from the loop above.
!   Twmx is the largest absolute value of tw1 and tw2 from the loop above.

  s1 =  ax * tsum
  eps1 = twmx * mpreal (10.d0) ** nint (22 - dpw * nwp)
  d1 = dplog10q (abs (s1 - s2))
  d2 = dplog10q (abs (s1 - s3))
  d3 = dplog10q (eps1) - 1.d0

  if (k <= 2) then
    err = 1.d0
  elseif (d1 .eq. -9999.d0) then
    err = 0.d0
  else
    nerr = nint (min (0.d0, max (d1 ** 2 / d2, 2.d0 * d1, d3, d4)))
    err = c10 ** nerr
  endif

!   Output current integral approximation and error estimate, to 56 dp.

  if (ndebug >= 2) then
    write (6, 2) k, nq1, nint (dplog10q (abs (err)))
2   format ('quadgs: Iteration',i3,' of',i3,'; est error = 10^',i5, &
      '; approx value =')
!    call mpwrite (6, s1)
    call mpeform (s1, 72, 60, stout)
    write (6, '(72a1)') stout
  endif
  if (k >= 3 .and. err < eps1) goto 130

  if (-2 * nerr > nint (dpw * nwp)) then
    nwp = min (2 * nwp, nwords1)
    call mpsetprecwords (nwp)
  endif
enddo

write (6, 3) nint (dplog10q (abs (err))), nquadl
3 format ('quadgs: Estimated error = 10^',i5/&
  'Increase Quadlevel for greater accuracy. Current Quadlevel =',i4)
if (err > 1.d-20) then
  write (6, 4)
4 format ('quadgs: Poor results may be due to singularities at endpoints.'/&
  'If so, try the erf or tanh-sinh quadrature routines (Quadtype = 2 or 3).')
endif
goto 130

120 continue

nerror = ierror + 100
write (6, 6) nerror
6 format ('quadgs: Error in quadrature calculation; code =',i5)
s1 = 0.d0

130 continue

quadgs = s1
! mpoud = nds
call mpsetprecwords (nwords1)
return
end

function quaderf (kv, str1, x1, x2, tmp1, nq1, nq2, wk, xk)

!   This routine computes the integral of the function in fun on the interval
!   [x1, x2], with up to nq1 iterations, with a target tolerance of 10^nepsilon1.
!   xk and wk are precomputed tables of abscissas and weights.  The function
!   fun is not evaluated at x = x1 or x2.  The array x_k contains 1 minus
!   the abscissas; the wk array contains the weights at these abscissas.
!   Modified for the toolkit.

!   David H. Bailey     2004-07-27

use mpmodule
use globdata
implicit none
integer i, ierror, ip(0:100), iz1, iz2, izx, j, k, k1, k2, n, nds, nq1, nq2, &
  nqq1
parameter (izx = 5)
logical log1, log2
real*8 d1, d2, d3, d4, dplog10q, h
type (mp_real) ax, bx, c10, eps1, eps2, epsilon1, epsilon2, err, fun, &
  quaderf, tsum, s1, s2, s3, t1, t2, t3, tw1, tw2, twi1, twi2, twmx, &
  wk(-1:nq2), xk(-1:nq2), x1, x2, xki, xt1, xx1, xx2
integer kv, ntmp1
character*1 stout(72)
character*2048 str1
type (mp_real) tmp1(ntmpx)
external dplog10q

! nds = mpoud
! mpoud = 56
call mpsetprecwords (nwords2)
ax = 0.5d0 * (x2 - x1)
bx = 0.5d0 * (x2 + x1)
call mpsetprecwords (nwords1)
epsilon1 = mpreal (10.d0) ** nepsilon1
epsilon2 = mpreal (10.d0) ** nepsilon2
tsum = 0.d0
s1 = 0.d0
s2 = 0.d0
h = 4.d0
c10 = 10.d0

if (wk(-1) < dble (nq1)) then
  write (6, 1) nq1
1 format ('quaderf: quadrature arrays have not been initialized; nq1 =',i6)
  nerror = 70
  goto 140
endif
nqq1 = dble (wk(-1))
n = dble (xk(-1))

do k = 0, nqq1
  ip(k) = 2 ** k
enddo

do k = 1, nq1
  h = 0.5d0 * h
  s3 = s2
  s2 = s1
  k1 = ip(nqq1-k)
  k2 = ip(nqq1-k+1)
  iz1 = 0
  iz2 = 0
  twmx = 0.d0

!   Evaluate function at level k in x, avoiding unnecessary computation.

  do i = 0, n, k1
    if (mod (i, k2) /= 0 .or. k == 1) then
      call mpsetprecwords (nwords2)
      xki = xk(i)
      xt1 = 1.d0 - xki
      xx1 = - ax * xt1 + bx
      xx2 = ax * xt1 + bx
      log1 = xx1 > x1
      log2 = xx2 < x2
      call mpsetprecwords (nwords1)

      if (log1 .and. iz1 < izx) then
        var(kv) = xx1
        call parse (str1, ntmp1, tmp1(1))  
        t1 = tmp1(1)  
        call mpgetpar ('mpier', ierror)
        if (ierror > 0 .or. nerror > 0) goto 130
        tw1 = t1 * wk(i)
        twi1 = abs (tw1)
        if (twi1 < epsilon1) then
          iz1 = iz1 + 1
        else
          iz1 = 0
        endif
      else
        t1 = 0.d0
        tw1 = 0.d0
      endif

      if (i > 0 .and. log2 .and. iz2 < izx) then
        var(kv) = xx2
        call parse (str1, ntmp1, tmp1(1)) 
        t2 = tmp1(1)  
        call mpgetpar ('mpier', ierror)
        if (ierror > 0 .or. nerror > 0) goto 130
        tw2 = t2 * wk(i)
        twi2 = abs (tw2)
        if (abs (tw2) < epsilon1) then
          iz2 = iz2 + 1
        else
          iz2 = 0
        endif
      else
        t2 = 0.d0
        tw2 = 0.d0
      endif

      tsum = tsum + tw1 + tw2
      twmx = max (twmx, abs (tw1), abs (tw2))
    endif
  enddo

!   Compute s1 = current integral approximation and err = error estimate.
!   Tsum is the sum of all tw1 and tw2 from the loop above.
!   Twmx is the largest absolute value of tw1 and tw2 from the loop above.
!   Twi1 and twi2 are the final nonzero values of abs(tw1) and abs(tw2).

  s1 =  ax * h * tsum
  eps1 = twmx * epsilon1
  eps2 = max (twi1, twi2)
  d1 = dplog10q (abs (s1 - s2))
  d2 = dplog10q (abs (s1 - s3))
  d3 = dplog10q (eps1) - 1.d0
  d4 = dplog10q (eps2) - 1.d0

  if (k <= 2) then
    err = 1.d0
  elseif (d1 .eq. -9999.d0) then
    err = 0.d0
  else
    err = c10 ** nint (min (0.d0, max (d1 ** 2 / d2, 2.d0 * d1, d3, d4)))
  endif

!   Output current integral approximation and error estimate, to 56 dp.

  if (ndebug >= 2) then
    write (6, 2) k, nq1, nint (dplog10q (abs (err)))
2   format ('quaderf: Iteration',i3,' of',i3,'; est error = 10^',i5, &
      '; approx value =')
!    call mpwrite (6, s1)
    call mpeform (s1, 72, 60, stout)
    write (6, '(72a1)') stout
  endif

  if (k >= 3 .and. err < eps1) goto 140
  if (k >= 3 .and. err < eps2) goto 120
enddo

write (6, 3) nint (dplog10q (abs (err))), nquadl
3 format ('quaderf: Estimated error = 10^',i5/&
  'Increase Quadlevel for greater accuracy. Current Quadlevel =',i4)
goto 140

120 continue

write (6, 4) nint (dplog10q (abs (err))), ndigits2
4 format ('quaderf: Estimated error = 10^',i5/&
  'Increase secondary prec (Ndigits2) for greater accuracy. Current value =',i4)
goto 140

130 continue

if (ierror > 0) nerror = ierror + 100
write (6, 5) nerror
5 format ('quaderf: Error in quadrature calculation; code =',i5)
s1 = 0.d0

140 continue

quaderf = s1
! mpoud = nds
return
end

function quadts (kv, str1, x1, x2, tmp1, nq1, nq2, wk, xk)

!   This routine computes the integral of the function in fun on the interval
!   [x1, x2], with up to nq1 iterations, with a target tolerance of 10^nepsilon1.
!   wk and xk are precomputed tables of abscissas and weights.  The function
!   fun is not evaluated at x = x1 or x2.  The array x_k contains 1 minus
!   the abscissas; the wk array contains the weights at these abscissas.
!   Modified for the Toolkit.

!   David H. Bailey     2004-07-27

use mpmodule
use globdata
implicit none
integer i, ierror, ip(0:100), iz1, iz2, izx, j, k, k1, k2, n, nds, nq1, nq2, &
  nqq1
parameter (izx = 5)
logical log1, log2
real*8 d1, d2, d3, d4, dplog10q, h
type (mp_real) ax, bx, c10, eps1, eps2, epsilon1, epsilon2, err, fun, &
  quadts, tsum, s1, s2, s3, t1, t2, t3, tw1, tw2, twi1, twi2, twmx, &
  wk(-1:nq2), xk(-1:nq2), x1, x2, xki, xt1, xx1, xx2
integer kv, ntmp1
character*1 stout(72)
character*2048 str1
type (mp_real) tmp1(ntmpx)
external dplog10q

! nds = mpoud
! mpoud = 56
call mpsetprecwords (nwords2)
ax = 0.5d0 * (x2 - x1)
bx = 0.5d0 * (x2 + x1)
call mpsetprecwords (nwords1)
epsilon1 = mpreal (10.d0) ** nepsilon1
epsilon2 = mpreal (10.d0) ** nepsilon2
tsum = 0.d0
s1 = 0.d0
s2 = 0.d0
h = 1.d0
c10 = 10.d0

if (wk(-1) < dble (nq1)) then
  write (6, 1) nq1
1 format ('quadts: quadrature arrays have not been initialized; nq1 =',i6)
  nerror = 70
  goto 140
endif
nqq1 = dble (wk(-1))
n = dble (xk(-1))

do k = 0, nqq1
  ip(k) = 2 ** k
enddo

do k = 1, nq1
  h = 0.5d0 * h
  s3 = s2
  s2 = s1
  k1 = ip(nqq1-k)
  k2 = ip(nqq1-k+1)
  iz1 = 0
  iz2 = 0
  twmx = 0.d0

!   Evaluate function at level k in x, avoiding unnecessary computation.

  do i = 0, n, k1
    if (mod (i, k2) /= 0 .or. k == 1) then
      call mpsetprecwords (nwords2)
      xki = xk(i)
      xt1 = 1.d0 - xki
      xx1 = - ax * xt1 + bx
      xx2 = ax * xt1 + bx
      log1 = xx1 > x1
      log2 = xx2 < x2
      call mpsetprecwords (nwords1)

      if (log1 .and. iz1 < izx) then
        var(kv) = xx1
        call parse (str1, ntmp1, tmp1(1))  
        t1 = tmp1(1)  
        call mpgetpar ('mpier', ierror)
        if (ierror > 0 .or. nerror > 0) goto 130
        tw1 = t1 * wk(i)
        twi1 = abs (tw1)
        if (twi1 < epsilon1) then
          iz1 = iz1 + 1
        else
          iz1 = 0
        endif
      else
        t1 = 0.d0
        tw1 = 0.d0
      endif

      if (i > 0 .and. log2 .and. iz2 < izx) then
        var(kv) = xx2
        call parse (str1, ntmp1, tmp1(1)) 
        t2 = tmp1(1)  
        call mpgetpar ('mpier', ierror)
        if (ierror > 0 .or. nerror > 0) goto 130
        tw2 = t2 * wk(i)
        twi2 = abs (tw2)
        if (abs (tw2) < epsilon1) then
          iz2 = iz2 + 1
        else
          iz2 = 0
        endif
      else
        t2 = 0.d0
        tw2 = 0.d0
      endif

      tsum = tsum + tw1 + tw2
      twmx = max (twmx, abs (tw1), abs (tw2))
    endif
  enddo

!   Compute s1 = current integral approximation and err = error estimate.
!   Tsum is the sum of all tw1 and tw2 from the loop above.
!   Twmx is the largest absolute value of tw1 and tw2 from the loop above.
!   Twi1 and twi2 are the final nonzero values of abs(tw1) and abs(tw2).

  s1 =  ax * h * tsum
  eps1 = twmx * epsilon1
  eps2 = max (twi1, twi2)
  d1 = dplog10q (abs (s1 - s2))
  d2 = dplog10q (abs (s1 - s3))
  d3 = dplog10q (eps1) - 1.d0
  d4 = dplog10q (eps2) - 1.d0

  if (k <= 2) then
    err = 1.d0
  elseif (d1 .eq. -9999.d0) then
    err = 0.d0
  else
    err = c10 ** nint (min (0.d0, max (d1 ** 2 / d2, 2.d0 * d1, d3, d4)))
  endif

!   Output current integral approximation and error estimate, to 56 dp.

  if (ndebug >= 2) then
    write (6, 2) k, nq1, nint (dplog10q (abs (err)))
2   format ('quadts: Iteration',i3,' of',i3,'; est error = 10^',i5, &
      '; approx value =')
!    call mpwrite (6, s1)
    call mpeform (s1, 72, 60, stout)
    write (6, '(72a1)') stout
  endif
  if (k >= 3 .and. err < eps1) goto 140
  if (k >= 3 .and. err < eps2) goto 120
enddo

write (6, 3) nint (dplog10q (abs (err))), nquadl
3 format ('quadts: Estimated error = 10^',i5/&
  'Increase Quadlevel for greater accuracy. Current Quadlevel =',i4)
goto 140

120 continue

write (6, 4) nint (dplog10q (abs (err))), ndigits2
4 format ('quadts: Estimated error = 10^',i5/&
  'Increase secondary prec (Ndigits2) for greater accuracy. Current value =',i4)
goto 140

130 continue

if (ierror > 0) nerror = ierror + 100
write (6, 5) nerror
5 format ('quadts: Error in quadrature calculation; code =',i5)
s1 = 0.d0

140 continue

quadts = s1
! mpoud = nds
return
end
