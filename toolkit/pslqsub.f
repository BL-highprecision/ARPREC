!   Combined pslqm1/pslqm2 subroutine file.
!   Modified for the toolkit.
!   David H Bailey 2004-06-11

!------------------------------

!   The following code performs the one-level, multi-pair PSLQ algorithm.
!   David H. Bailey     2004-06-10

subroutine pslqm1 (idb, gam, n, nsq, rb, x, is0, s1, s2, iq, r)

!   This routine allocates the scratch space arrays is0, s1 and s2.  Handling
!   scratch space in this manner is not really necessary in Fortran-90, but 
!   this design facilitates straightforward translation to Fortran-77.

!   Arguments are as follows:
!     idb = Debug flag (0-3).
!     gam = PSLQ gamma parameter, normally set to sqrt(4/3).
!     n   = Dimension of input vector.
!     nsq = Size of tables used in itermpw.
!     rb  = Log10 of max size (Euclidean norm) of acceptable relation (type DP).
!     x   = Input real vector (type mp_real).
!     is0 = Scratch array (type integer -- see size below).
!     s1  = Scratch array (type double precision -- see size below).
!     s2  = Scratch array (type mp_real -- see size below).
!     iq  = Output flag: 0 (unsuccessful) or 1 (successful).
!     r   = Output integer coefficient vector (type mp_real).

use mpmodule
implicit none
integer idb, iq, n, nsq
double precision gam, rb
integer is0(3*n)
double precision s1(n)
type (mp_real) s2(3*n*n+nsq*n+2*n), x(n), r(n)

call pslqm1x (idb, gam, n, nsq, rb, x, is0(1), is0(n+1), is0(2*n+1), &
  s1(1), s2(1), s2(n*n+1), s2(2*n*n+1), s2(3*n*n+1), &
  s2(3*n*n+nsq*n+1), s2(3*n*n+nsq*n+n+1), iq, r)
return
end

subroutine pslqm1x (idb, gam, n, nsq, rb, x, ip, ir, is, dx, b, h, t, &
  syq, s, y, iq, r)

!   The following parameters are set in this routine:
!     ipi = Iteration print interval when idb >= 2.
!     ipm = Iteration check interval for MP iterations.
!     itm = Maximum iteration count.  Run is aborted when this is exceeded.

use mpmodule
implicit none
integer i, idb, imq, ipi, ipm, iq, it, itm, izm, j, j1, n, n1, n2, n3, n4, nsq
parameter (ipi = 25, ipm = 100, itm = 100000)
double precision d1, d2, d3, d4, dplog10, gam, rb, second, tm0, tm1, times(2)
integer ip(n), ir(n), is(n)
double precision dx(n)
type (mp_real) b(n,n), h(n,n), s(n), syq(n,nsq), t(n,n), r(n), x(n), &
  y(n), bound, rn, t1, t2, t3, t4
external bound, dplog10, second

!   Initialize.

if (idb .ge. 2) write (6, 1) n
1 format ('PSLQM1 integer relation detection: n =',i5)
iq = 0
it = 0
imq = 0
rn = 0.d0

do i = 1, 2
  times(i) = 0.d0
enddo

if (idb .ge. 2) write (6, 2) it
2 format ('Iteration',i7,3x,'MP initialization')
tm0 = second ()
call initmp (idb, n, nsq, ip, dx, b, h, s, syq, x, y)
tm1 = second ()
times(1) = tm1 - tm0

!   Perform MP iterations.

if (idb .ge. 2) write (6, 3) it
3 format ('Iteration',i7,3x,'Start MP iterations')

100 continue

it = it + 1
if (idb .eq. 3 .or. idb .ge. 2 .and. mod (it, ipi) .eq. 0) write (6, 4) it
4 format ('Iteration',i7)
tm0 = second ()
call itermp (idb, gam, it, n, nsq, ip, ir, is, dx, b, h, s, syq, t, y, &
  imq, izm)
tm1 = second ()
times(2) = times(2) + (tm1 - tm0)

if (mod (it, ipm) .eq. 0) then

!   Find min and max absolute value of y vector.

  t1 = 1.d300
  t2 = 0.d0

  do i = 1, n
    t1 = min (t1, abs (y(i)))
    t2 = max (t2, abs (y(i)))
  enddo

  if (idb .ge. 2) then
    call decmd (t1, d1, n1)
    call decmd (t2, d2, n2)
    write (6, 5) it, d1, n1, d2, n2
5   format ('Iteration',i7,3x,'Min, max of y =',f11.6,'D',i5, &
      f11.6,'D',i5)
  endif

!   Compute norm bound.

  t1 = bound (n, h)
  rn = max (rn, t1)
  if (idb .ge. 2) then
    call decmd (t1, d1, n1)
    call decmd (rn, d2, n2)
    write (6, 6) it, d1, n1, d2, n2
6   format ('Iteration',i7,3x,'Norm bound =',f11.6,'D',i5,4x,'Max. bound =', &
      f11.6,'D',i5)
  endif
  if (it .gt. itm) then
    if (idb .ge. 1) write (6, 7) itm
7   format ('Iteration limit exceeded',i7)
    goto 120
  endif
  if (dplog10 (rn) .gt. rb) then
    if (idb .ge. 1) write (6, 8) rb
8   format ('Norm bound limit exceeded.',1p,d15.6)
    goto 120
  endif
endif

if (izm .eq. 0) then
  goto 100
elseif (izm .eq. 1) then
  goto 110
elseif (izm .eq. 2) then
  goto 120
endif

110 continue

!   A relation has been detected.

tm0 = second ()
t1 = 1.d300
t2 = 0.d0
t3 = 0.d0

!   Select the relation corresponding to the smallest y entry and compute norm.

do j = 1, n
  if (abs (y(j)) < t1) then
    j1 = j
    t1 = abs (y(j))
  endif
  t2 = max (t2, abs (y(j)))
enddo

do i = 1, n
  r(i) = b(j1,i)
  t3 = t3 + r(i) ** 2
enddo

t3 = sqrt (t3)
call decmd (t3, d3, n3)
d3 = d3 * 10.d0 ** n3

!   Output the final norm bound and other info.

if (idb .ge. 1) then
  t4 = bound (n, h)
  call decmd (t1, d1, n1)
  call decmd (t2, d2, n2)
  call decmd (t3, d3, n3)
  call decmd (t4, d4, n4)
  write (6, 9) it, d1, n1, d2, n2, d4, n4
9 format ('Iteration',i7,3x,'Relation detected'/ &
  'Min, max of y =',0p,f11.6,'D',i5,f11.6,'D',i5/'Max. bound =',f11.6,'D',i5)
  write (6, 10) j1, d3, n3, d1, n1
10 format ('Index of relation =',i4,3x,'Norm =',f11.6,'D',i5,3x, &
  'Residual =',f11.6,'D',i5)
endif

if (dplog10 (t3) .le. rb) then
  iq = 1
else
  if (idb .ge. 2) write (6, 11)
11 format ('Relation is too large.')
endif

120 continue

if (idb .ge. 2) write (6, 12) times
12 format ('CPU times:'/(5f12.2))

return
end

!------------------------------

!   First-level subroutines.

subroutine initmp (idb, n, nsq, ix, dx, b, h, s, syq, x, y)

!   This initializes MP arrays at the beginning.

use mpmodule
implicit none
integer i, idb, j, n, nsq
double precision d1
integer ix(n)
double precision dx(n)
type (mp_real) b(n,n), h(n,n), s(n), syq(n,nsq), x(n), y(n), t1

if (idb .ge. 3) then
  write (6, 1)
1 format ('initmp: Input x vector:')
  call matoutmd (1, n, ix, dx, x)
endif

!   Set b to the identity matrix.

do j = 1, n
  do i = 1, n
    b(i,j) = 0.d0
  enddo

  b(j,j) = 1.d0
enddo

t1 = 0.d0

!   Compute the x vector, the square root of the partial sum of squares of x,
!   and the y vector, which is the normalized x vector.

do i = n, 1, -1
  t1 = t1 + x(i) ** 2
  s(i) = sqrt (t1)
enddo

t1 = 1.d0 / s(1)

do i = 1, n
  y(i) = t1 * x(i)
  s(i) = t1 * s(i)
enddo

!   Compute the initial h matrix.

!$omp parallel do private (i, j, t1) copyin (mpnw)
do j = 1, n - 1
  do i = 1, j - 1
    h(i,j) = 0.d0
  enddo

  h(j,j) = s(j+1) / s(j)
  t1 = y(j) / (s(j) * s(j+1))

  do i = j + 1, n
    h(i,j) = - y(i) * t1
  enddo
enddo
!$omp end parallel do

!   Zero the syq array.

do j = 1, nsq
  do i = 1, n
    syq(i,j) = 0.d0
  enddo
enddo

if (idb .ge. 3) then
  write (6, 2)
2 format ('initmp: Initial y vector:')
  call matoutmd (1, n, ix, dx, y)
  write (6, 3)
3 format ('initmp: Initial h matrix:')
  call matoutmd (n, n - 1, ix, dx, h)
endif

return
end

subroutine itermp (idb, gam, it, n, nsq, ip, ir, is, dx, b, h, q, syq, t, &
  y, imq, izm)

!   This performs one iteration of the PSLQM algorithm using MP arithmetic.

use mpmodule
use globdata
implicit none
integer i, idb, ii, ij, im, im1, imq, it, izm, j, j1, j2, k, mpr, mq, n, &
  n1, nsq, ntl
parameter (ntl = 72)
double precision d1, d2, gam
integer ip(n), ir(n), is(n)
double precision dx(n)
type (mp_real) b(n,n), h(n,n), q(n), syq(n,nsq), t(n,n), y(n)
type (mp_real) eps, t1, t2, t3, t4, teps

eps = mpreal (10.d0) ** (nepsilon1+10)
teps = 2.d0 ** ntl * eps
izm = 0
mpr = nint (0.4d0 * n)

!   Compute q vector = {gam^i * |h(i,i)|}, then sort in ascending order.

do i = 1, n - 1
  q(i) = gam ** i * abs (h(i,i))
enddo

call qsortmp (n - 1, q, ip)

!   Select up to mpr disjoint pairs of indices (m,m+1), where m is an index
!   from the list of the largest q(i).

do i = 1, n
  is(i) = 0
enddo

if (imq .eq. 0) then
  mq = mpr
else
  mq = 1
  imq = 0
endif
ii = n

do i = 1, mq
100 continue
  ii = ii - 1
  if (ii .eq. 0) then
    mq = i - 1
    goto 110
  endif
  j1 = ip(ii)
  j2 = j1 + 1
  if (is(j1) .ne. 0 .or. is(j2) .ne. 0) goto 100
  ir(i) = j1
  is(j1) = 1
  is(j2) = 1
enddo

110 continue

!   Exchange the pairs of entries of y, and rows of b and h.

!$omp parallel do private (i, j, im, im1, t1) copyin (mpnw)
do j = 1, mq
  im = ir(j)
  im1 = im + 1
  t1 = y(im)
  y(im) = y(im1)
  y(im1) = t1

  do i = 1, n
    t1 = b(im,i)
    b(im,i) = b(im1,i)
    b(im1,i) = t1
  enddo

  do i = 1, n - 1
    t1 = h(im,i)
    h(im,i) = h(im1,i)
    h(im1,i) = t1
  enddo
enddo
!$omp end parallel do

!   Eliminate the "corners" produced by the above permutation in h.

!$omp parallel do private (i, j, im, im1, t1, t2, t3, t4) copyin (mpnw)
do j = 1, mq
  im = ir(j)
  im1 = im + 1
  if (im .le. n - 2) then
    t1 = h(im,im)
    t2 = h(im,im1)
    t3 = sqrt (t1 ** 2 + t2 ** 2)
    t1 = t1 / t3
    t2 = t2 / t3

    do i = im, n
      t3 = h(i,im)
      t4 = h(i,im1)
      h(i,im) = t1 * t3 + t2 * t4
      h(i,im1) = - t2 * t3 + t1 * t4
    enddo
  endif
enddo
!$omp end parallel do

!   Perform reduction on h, using the diagonal scheme.  Multipliers are
!   saved in the t array.

do i = 2, n
!$omp parallel do private (ij, j, k) copyin (mpnw)
  do j = 1, n - i + 1
    ij = i + j - 1

    do k = j + 1, ij - 1
      h(ij,j) = h(ij,j) - t(ij,k) * h(k,j)
    enddo

    t(ij,j) = anint (h(ij,j) / h(j,j))
    h(ij,j) = h(ij,j) - t(ij,j) * h(j,j)
  enddo
!$omp end parallel do
enddo

!   Update y, using the t array.  Find min absolute value of y.

t1 = abs (y(n))
j1 = n

do j = 1, n - 1
  do i = j + 1, n
    y(j) = y(j) + t(i,j) * y(i)
  enddo

  if (abs (y(j)) < t1) then
    j1 = j
    t1 = abs (y(j))
  endif
enddo

!   Update b, using the t array.

!$omp parallel do private (i, j, k) copyin (mpnw)
do k = 1, n
  do j = 1, n - 1
    do i = j + 1, n
      b(j,k) = b(j,k) + t(i,j) * b(i,k)
    enddo
  enddo
enddo
!$omp end parallel do

!  Find the largest entry of b in the same row as the smallest y.

t2 = 0.d0

do i = 1, n
  t2 = max (t2, abs (b(j1,i)))
enddo

if (t1 .le. t2 * teps) then
  if (idb .ge. 2) then
    call decmd (t1, d1, n1) 
    write (6, 1) it, d1, n1
1   format ('Iteration',i7,3x,'itermp: Small value in y =',f11.6,'D',i5)
  endif
  if (t1 .le. t2 * eps) then
    izm = 1
  else
    if (idb .ge. 1) write (6, 2) it
2   format ('Iteration',i7,3x,'itermp: Precision exhausted.'/ &
      'Try increasing epsilon to a less negative value.')
    izm = 2
  endif
endif

!   Compare the y vector with those of recent iterations.  If a duplicate is
!   found, then the next iteration must be performed with mq = 1.

do j = 1, nsq
  t1 = 0.d0

  do i = 1, n
    t1 = max (t1, abs (y(i) - syq(i,j)))
  enddo

  if (t1 .le. t2 * teps) then
    if (idb .ge. 2) write (6, 3) it, j
 3  format ('Iteration',i7,3x,'itermp: Duplicate found, j =',i6)
    imq = 1
    goto 120
  endif
enddo

!   Place the vector y in the table syq.

120   continue
k = 1 + mod (it, nsq)

do i = 1, n
  syq(i,k) = y(i)
enddo

if (idb .ge. 3) then
  write (6, 4)
4 format ('itermp: Updated y:')
!  call matoutmd (1, n, ip, dx, y)
  call matoutmp (1, n, y)
  write (6, 5)
5 format ('itermp: Updated b matrix:')
  call matoutmd (n, n, ip, dx, b)
  write (6, 6)
6 format ('itermp: Updated h matrix:')
  call matoutmd (n, n - 1, ip, dx, h)
endif

return
end

!------------------------------

!   Second- and third-level subroutines.

function bound (n, h)

!   This computes the norm bound using MP arithmetic.

use mpmodule
implicit none
integer i, n
type (mp_real) bound,  h(n,n), t1, t2

t1 = 0.d0

do i = 1, n - 1
  t1 = max (t1, abs (h(i,i)))
enddo

bound = 1.d0 / t1

return
end

function dplog10 (a)

!   For input MP value a, this routine returns a DP approximation to log10 (a).

use mpmodule
implicit none
integer ia
double precision da, dplog10, t1
type (mp_real) a

call mpmdc (a%mpr, da, ia)
if (da .le. 0.d0) then
  dplog10 = -1d300
else
  dplog10 = log10 (abs (da)) + ia * log10 (2.d0)
endif

100 continue
return
end

subroutine decmd (a, b, ib)

!   For input MP value a, this routine returns DP b and integer ib such that 
!   a = b * 10^ib, with 1 <= abs (b) < 10 for nonzero a.

use mpmodule
implicit none
integer ia, ib
double precision da, b, t1, xlt
parameter (xlt = 0.3010299956639812d0)
type (mp_real) a

call mpmdc (a%mpr, da, ia)
if (da .ne. 0.d0) then
  t1 = xlt * ia + log10 (abs (da))
  ib = t1
  if (t1 .lt. 0.d0) ib = ib - 1
  b = sign (10.d0 ** (t1 - ib), da)
else
  b = 0.d0
  ib = 0
endif

return
end

subroutine matoutmd (n1, n2, ix, dx, a)

!   This outputs the MP matrix a as a DP matrix.

use mpmodule
implicit none
integer i, j, n1, n2
integer ix(n2)
double precision dx(n2)
type (mp_real) a(n1,n2)

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)

  do j = 1, n2
    call decmd (a(i,j), dx(j), ix(j))
  enddo

  write (6, 2) (dx(j), ix(j), j = 1, n2)
2 format (4(f13.8,'D',i5))
enddo

return
end

subroutine matoutmp (n1, n2, a)

!   This outputs the MP matrix a.  It may be used in place of calls to matoutmd
!   in the code above if greater accuracy is desired in debug output.

use mpmodule
implicit none
integer i, j, n1, n2
type (mp_real) a(n1,n2)

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)

  do j = 1, n2
    call mpwrite (6, a(i,j))
  enddo
enddo

return
end

subroutine qsortmp (n, a, ip)

!   This routine sorts the entries of the N-long MP vector A into ascending
!   order using the quicksort algorithm.  The permutation vector that would
!   sort the vector is returned in IP.

use mpmodule
implicit none
integer i, iq, it, j, jq, jz, k, l, n
integer ip(n), ik(50), jk(50)
type (mp_real) a(n), s0, s1, s2

do i = 1, n
  ip(i) = i
enddo

if (n .eq. 1) return

k = 1
ik(1) = 1
jk(1) = n

130 i = ik(k)
j = jk(k)
iq = i
jq = j
it = (i + j + 1) / 2
l = ip(j)
ip(j) = ip(it)
ip(it) = l
s0 = a(ip(j))
j = j - 1

140 continue

do l = i, j
  if (s0 .lt. a(ip(l))) goto 160
enddo

i = j
goto 190

160 i = l

do l = j, i, -1
  if (s0 .gt. a(ip(l))) goto 180
enddo

j = i
goto 190

180 j = l
if (i .ge. j) goto 190
l = ip(i)
ip(i) = ip(j)
ip(j) = l
goto 140

190 continue
if (s0 .ge. a(ip(i))) goto 200
l = ip(jq)
ip(jq) = ip(i)
ip(i) = l

200 k = k - 1
jz = 0
if (j .eq. iq) goto 210
k = k + 1
jk(k) = j
jz = 1

210 i = i + 1
if (i .eq. jq) goto 220
k = k + 1
ik(k) = i
jk(k) = jq
if (jz .eq. 0) goto 220
if (j - iq .ge. jq - i) goto 220
ik(k-1) = i
jk(k-1) = jq
ik(k) = iq
jk(k) = j

220 if (k .gt. 0) goto 130

return
end

!------------------------------

!   The following code performs the two-level, multi-pair PSLQ algorithm.
!   David H. Bailey     2004-06-10

subroutine pslqm2 (idb, gam, n, nsq, rb, x, is0, s1, s2, iq, r)

!   This routine allocates the scratch space arrays is0, s1 and s2.  Handling
!   scratch space in this manner is not really necessary in Fortran-90, but 
!   this design facilitates straightforward translation to Fortran-77.

!   Arguments are as follows:
!     idb = Debug flag (0-3).
!     gam = PSLQ gamma parameter, normally set to sqrt(4/3).
!     n   = Dimension of input vector.
!     nsq = Size of tables used in iterdp and itermp.
!     rb  = Log10 of max size (Euclidean norm) of acceptable relation (type DP).
!     x   = Input real vector (type mp_real).
!     is0 = Scratch array (type integer -- see size below).
!     s1  = Scratch array (type double precision -- see size below).
!     s2  = Scratch array (type mp_real -- see size below).
!     iq  = Output flag: 0 (unsuccessful) or 1 (successful).
!     r   = Output integer coefficient vector (type mp_real).

use mpmodule
implicit none
integer idb, iq, n, nsq
double precision gam, rb
integer is0(3*n)
double precision s1(7*n*n+nsq*n+3*n)
type (mp_real) s2(3*n*n+nsq*n+2*n), x(n), r(n)

call pslqm2x (idb, gam, n, nsq, rb, x, is0(1), is0(n+1), is0(2*n+1), &
  s1(1), s1(n*n+1), s1(2*n*n+1), s1(3*n*n+1), s1(4*n*n+1), s1(5*n*n+1), &
  s1(6*n*n+1), s1(7*n*n+1), s1(7*n*n+nsq*n+1), s1(7*n*n+nsq*n+n+1), &
  s1(7*n*n+nsq*n+2*n+1), s2(1), s2(n*n+1), s2(2*n*n+1), s2(3*n*n+1), &
  s2(3*n*n+nsq*n+1), s2(3*n*n+nsq*n+n+1), iq, r)
return
end

subroutine pslqm2x (idb, gam, n, nsq, rb, x, ip, ir, is, da, db, dh, dsa, &
  dsb, dsh, dt, dsyq, dq, dy, dsy, b, h, t, syq, s, y, iq, r)

!   The following parameters are set in this routine:
!     ipi = Iteration print interval when idb >= 2.
!     ipm = Iteration check interval for MP iterations; also save interval.
!     itm = Maximum iteration count.  Run is aborted when this is exceeded.

use mpmodule
implicit none
integer i, idb, imq, ipi, ipm, iq, it, itm, its, izd, izm, j, j1, n, n1, n2, &
  n3, nsq
parameter (ipi = 100, ipm = 10, itm = 1000000)
double precision bounddp, d1, d2, d3, d4, gam, rb, rn, second, tm0, tm1, &
  times(4)
integer ip(n), ir(n), is(n)
double precision da(n,n), db(n,n), dh(n,n), dq(n), dsa(n,n), dsb(n,n), &
  dsh(n,n), dsyq(n,nsq), dt(n,n), dy(n), dsy(n)
type (mp_real) b(n,n), h(n,n), s(n), syq(n,nsq), t(n,n), r(n), x(n), &
  y(n), t1, t2, t3
external bounddp, second

!   Initialize.

if (idb .ge. 2) write (6, 1) n
1 format ('PSLQM2 integer relation detection: n =',i5)
iq = 0
it = 0
imq = 0
rn = 0.d0

do i = 1, 4
  times(i) = 0.d0
enddo

if (idb .ge. 2) write (6, 2) it
2 format ('Iteration',i7,3x,'MP initialization')
tm0 = second ()
call initmp (idb, n, nsq, ip, dt, b, h, s, syq, x, y)
tm1 = second ()
times(1) = tm1 - tm0

100 continue

!   Initialize DP arrays from MP arrays.

if (idb .ge. 3) write (6, 3) it
3 format ('Iteration',i7,3x,'Start DP iterations')
call initdp (idb, it, n, nsq, da, db, dh, dy, dsyq, h, y, izd)
call savedp (n, da, db, dh, dy, dsa, dsb, dsh, dsy)
its = it

if (izd .eq. 0) then

!   Perform an LQ decomposition on dh, prior to DP iterations.

  call lqdp (n, n - 1, dh)

!   Perform DP iterations.

110 continue

  it = it + 1
  if (idb .ge. 3 .or. idb .ge. 2 .and. mod (it, ipi) .eq. 0) write (6, 4) it
4 format ('Iteration',i7)
  tm0 = second ()
  call iterdp (idb, gam, it, n, nsq, ip, ir, is, da, db, dh, dq, dsyq, &
    dt, dy, imq, izd)
  tm1 = second ()
  times(2) = times(2) + (tm1 - tm0)

  if (izd .eq. 0) then
    if (mod (it - its, ipm) .eq. 0) then
      call savedp (n, da, db, dh, dy, dsa, dsb, dsh, dsy)
      its = it
    endif
    goto 110
  else

!   Update MP arrays from the DP arrays.

    if (izd .eq. 2) then

!   DP iteration was aborted -- revert to previous data.

      it = its
      call savedp (n, dsa, dsb, dsh, dsy, da, db, dh, dy)
    endif

    if (idb .ge. 3) write (6, 5) it
5   format ('Iteration',i7,3x,'MP update')
    tm0 = second ()
    call updtmp (idb, it, n, ip, dt, da, db, b, h, t, y, izm)
    tm1 = second ()
    times(3) = times(3) + (tm1 - tm0)

!   Compute norm bound.

    do j = 1, n - 1
      do i = 1, n
        dh(i,j) = h(i,j)
      enddo
    enddo
    
    call lqdp (n, n - 1, dh)
    d3 = bounddp (n, dh)
    if (d3 == -1.d0) goto 150
    rn = max (rn, d3)
    if (idb .ge. 2) then
      write (6, 6) it, d3, rn
6     format ('Iteration',i7,3x,'Norm bound =',1p,d15.6,4x,'Max. bound =', &
        1p,d15.6)
    endif
    if (it .gt. itm) then
      if (idb .ge. 1) write (6, 7) itm
7     format ('Iteration limit exceeded',i7)
      goto 150
    endif
    if (log10 (rn) .gt. rb) then
      if (idb .ge. 1) write (6, 8) rb
8     format ('Norm bound limit exceeded.',1p,d15.6)
      goto 150
    endif
    if (izm .eq. 0) then
      if (izd .eq. 2) then
        goto 120
      else
        goto 100
      endif
    elseif (izm .eq. 1) then
      goto 140
    elseif (izm .eq. 2) then
      goto 150
    endif
  endif
endif

120 continue

!   Perform an LQ decomposition on h, prior to MP iterations.

if (idb .ge. 2) write (6, 9) it
9 format ('Iteration',i7,3x,'Start MP iterations')
tm0 = second ()
call lqmp (n, n - 1, h)
tm1 = second ()
times(1) = times(1) + (tm1 - tm0)

!   Perform MP iterations.

130 continue

it = it + 1
if (idb .ge. 2) write (6, 10) it
10 format ('Iteration',i7)
tm0 = second ()
call itermp (idb, gam, it, n, nsq, ip, ir, is, dt, b, h, s, syq, t, y, &
  imq, izm)
tm1 = second ()
times(4) = times(4) + (tm1 - tm0)

if (izm .eq. 0) then
  if (mod (it - its, ipm) .eq. 0) then

!   Check to see if DP iterations can be resumed.

    call check (idb, it, n, h, y, izd)
    if (izd .eq. 0) then
      if (idb .ge. 2) write (6, 11) it
11    format ('Iteration',i7,3x,'Return to DP iterations')
      goto 100
    endif
  endif
  goto 130
elseif (izm .eq. 1) then
  goto 140
elseif (izm .eq. 2) then
  goto 150
endif

!   A relation has been detected.

140 continue

tm0 = second ()
t1 = 1.d300
t2 = 0.d0
t3 = 0.d0

!   Select the relation corresponding to the smallest y entry and compute norm.

do j = 1, n
  if (abs (y(j)) < t1) then
    j1 = j
    t1 = abs (y(j))
  endif
  t2 = max (t2, abs (y(j)))
enddo

do i = 1, n
  r(i) = b(j1,i)
  t3 = t3 + r(i) ** 2
enddo

t3 = sqrt (t3)
call decmd (t3, d3, n3)
d3 = d3 * 10.d0 ** n3

!   Output the final norm bound and other info.

if (idb .ge. 1) then
  do j = 1, n - 1
    do i = 1, n
      dh(i,j) = h(i,j)
    enddo
  enddo

  call lqdp (n, n - 1, dh)
  d4 = bounddp (n, dh)
  rn = max (rn, d4)
  call decmd (t1, d1, n1)
  call decmd (t2, d2, n2)
  write (6, 12) it, d1, n1, d2, n2, rn
12 format ('Iteration',i7,3x,'Relation detected'/ &
  'Min, max of y =',0p,f11.6,'D',i5,f11.6,'D',i5/'Max. bound =',1p,d15.6)
  write (6, 13) j1, d3, d1, n1
13 format ('Index of relation =',i4,3x,'Norm =',1p,d15.6,3x, &
  'Residual =',0p,f11.6,'D',i5)
endif

if (log10 (d3) .le. rb) then
  iq = 1
else
  if (idb .ge. 2) write (6, 14)
14 format ('Relation is too large.')
endif

150 continue

if (idb .ge. 2) write (6, 15) times
15 format ('CPU times:'/(4f12.2))

return
end

!------------------------------

!   First-level subroutines.

subroutine check (idb, it, n, h, y, izd)

!   This checks h and y to see if DP iterations can be resumed.

use mpmodule
use globdata
implicit none
integer i, idb, it, izd, n, n1
double precision d1, deps
parameter (deps = 1d-10)
type (mp_real) h(n,n), y(n), t1, t2, t3, t4

call mpsetprecwords (5)
izd = 0
t1 = 1.d300
t2 = 0.d0

!   Find the min and max absolute value in the y vector.

do i = 1, n
  t1 = min (t1, abs (y(i)))
  t2 = max (t2, abs (y(i)))
enddo

t3 = t1 / t2
t4 = deps
if (t3 .lt. t4) then
  if (idb .ge. 2) then
    call decmd (t3, d1, n1)
    write (6, 1) it, d1, n1
1   format ('Iteration',i7,3x,'check: Small min/max ratio in y =',f11.6, &
    'D',i5)
  endif
  izd = 1
endif

call mpsetprecwords (nwords1)
return
end

subroutine initdp (idb, it, n, nsq, da, db, dh, dy, dsyq, h, y, izd)

!   This re-initializes the DP arrays from the MP arrays.

use mpmodule
use globdata
implicit none
integer i, idb, it, izd, j, n, n1, nsq
double precision da(n,n), db(n,n), dh(n,n), dy(n), dsyq(n,nsq), d1, deps
type (mp_real) h(n,n), y(n), t1, t2, t3, t4
parameter (deps = 1d-10)

call mpsetprecwords (5)
izd = 0
t1 = 1.d300
t2 = 0.d0

!   Find the min and max absolute value in the y vector.

do i = 1, n
  t1 = min (t1, abs (y(i)))
  t2 = max (t2, abs (y(i)))
enddo

!   If the dynamic range of the y vector is too great, set the izd flag and
!   abort this initialization.

t3 = t1 / t2
t4 = deps
if (t3 .lt. t4) then
  if (idb .ge. 3) then
    call decmd (t3, d1, n1)
    write (6, 1) it, d1, n1
1   format ('Iteration',i7,3x,'initdp: Small min/max ratio in y =',f11.6, &
    'D',i5)
  endif
  izd = 1
  goto 100
endif

!   Set dy to be the scaled y vector.

t1 = 1.d0 / t2

do i = 1, n
  dy(i) = t1 * y(i)
enddo

!   Find the maximum absolute value of the h matrix diagonals.

t2 = 0.d0

do j = 1, n - 1
  t2 = max (t2, abs (h(j,j)))
enddo

!   Set dh to be the scaled h matrix.

t1 = 1.d0 / t2

do j = 1, n - 1
  do i = 1, n
    dh(i,j) = t1 * h(i,j)
  enddo
enddo

!   Set da and db to the identity.

do j = 1, n
  do i = 1, n
    da(i,j) = 0.d0
    db(i,j) = 0.d0
  enddo

  da(j,j) = 1.d0
  db(j,j) = 1.d0
enddo

!   Zero the dsyq array.

do j = 1, nsq
  do i = 1, n
    dsyq(i,j) = 0.d0
  enddo
enddo

if (idb .ge. 3) then
  write (6, 2)
2 format ('initdp: Scaled dy vector:')
  call matoutdp (1, n, dy)
  write (6, 3)
3 format ('initdp: Scaled dh matrix:')
  call matoutdp (n, n - 1, dh)
endif

100 continue

call mpsetprecwords (nwords1)
return
end

subroutine iterdp (idb, gam, it, n, nsq, ip, ir, is, da, db, dh, dq, dsyq, &
  dt, dy, imq, izd)

!   This performs one iteration of the PSLQ algorithm using DP arithmetic.

use mpmodule
implicit none
integer i, idb, ii, ij, im, im1, imq, it, izd, j, j1, j2, k, mpr, mq, n, &
  n1, nsq
double precision deps, gam, t1, t2, t3, t4, tmx1, tmx2
integer ip(n), ir(n), is(n)
double precision da(n,n), db(n,n), dh(n,n), dq(n), dsyq(n,nsq), &
  dt(n,n), dy(n)
parameter (tmx1 = 1.d13, tmx2 = 2.d0**52, deps = 1.d-14)

izd = 0
mpr = nint (0.4d0 * n)

!   Compute dq vector = {gam^i * |dh(i,i)|}, then sort in ascending order.

do i = 1, n - 1
  dq(i) = gam ** i * abs (dh(i,i))
enddo

call qsortdp (n - 1, dq, ip)

!   Select up to mpr disjoint pairs of indices (m,m+1), where m is an index
!   from the list of the largest dq(i).

do i = 1, n
  is(i) = 0
enddo

if (imq .eq. 0) then
  mq = mpr
else
  mq = 1
  imq = 0
endif
ii = n

do i = 1, mq
100 continue
  ii = ii - 1
  if (ii .eq. 0) then
    mq = i - 1
    goto 110
  endif
  j1 = ip(ii)
  j2 = j1 + 1
  if (is(j1) .ne. 0 .or. is(j2) .ne. 0) goto 100
  ir(i) = j1
  is(j1) = 1
  is(j2) = 1
enddo

110 continue

!   Exchange the pairs of entries of dy, and rows of da, db and dh.

do j = 1, mq
  im = ir(j)
  im1 = im + 1
  t1 = dy(im)
  dy(im) = dy(im1)
  dy(im1) = t1

  do i = 1, n
    t1 = da(im,i)
    da(im,i) = da(im1,i)
    da(im1,i) = t1
    t1 = db(im,i)
    db(im,i) = db(im1,i)
    db(im1,i) = t1
  enddo

  do i = 1, n - 1
    t1 = dh(im,i)
    dh(im,i) = dh(im1,i)
    dh(im1,i) = t1
  enddo
enddo

!   Eliminate the "corners" produced by the above permutation in dh.

do j = 1, mq
  im = ir(j)
  im1 = im + 1
  if (im .le. n - 2) then
    t1 = dh(im,im)
    t2 = dh(im,im1)
    t3 = sqrt (t1 ** 2 + t2 ** 2)
    t1 = t1 / t3
    t2 = t2 / t3

    do i = im, n
      t3 = dh(i,im)
      t4 = dh(i,im1)
      dh(i,im) = t1 * t3 + t2 * t4
      dh(i,im1) = - t2 * t3 + t1 * t4
    enddo
  endif
enddo

!   Perform reduction on dh, using the diagonal scheme.  Multipliers are
!   saved in the dt array.

do i = 2, n
  do j = 1, n - i + 1
    ij = i + j - 1

    do k = j + 1, ij - 1
      dh(ij,j) = dh(ij,j) - dt(ij,k) * dh(k,j)
    enddo

    dt(ij,j) = anint (dh(ij,j) / dh(j,j))
    dh(ij,j) = dh(ij,j) - dt(ij,j) * dh(j,j)
  enddo
enddo

!   Update dy, using the dt array.  Find min absolute value of dy.

t1 = abs (dy(n))

do j = 1, n - 1
  do i = j + 1, n
    dy(j) = dy(j) + dt(i,j) * dy(i)
  enddo

  t1 = min (t1, abs (dy(j)))
enddo

!   Update da and db, using the dt array.  Find the max absolute value of
!   da and db entries as they are calculated (not merely at the end).

t2 = 0.d0

!$omp parallel do private (i, j, k, t3)
do k = 1, n
  t3 = 0.d0

  do j = 1, n - 1
    do i = j + 1, n
      da(i,k) = da(i,k) - dt(i,j) * da(j,k)
      db(j,k) = db(j,k) + dt(i,j) * db(i,k)
      t3 = max (t3, abs (da(i,k)), abs (db(j,k)))
    enddo
  enddo

  dq(k) = t3
enddo
!$omp end parallel do

do k = 1, n
  t2 = max (t2, dq(k))
enddo

if (t1 .le. deps) then
  if (idb .ge. 3) write (6, 1) it, t1
1 format ('Iteration',i7,3x,'iterdp: Small value in dy =',1pd15.6)
  izd = 1
endif

if (t2 .gt. tmx1 .and. t2 .le. tmx2) then
  if (idb .ge. 3) write (6, 2) it, t2
2 format ('Iteration',i7,3x,'iterdp: Large value in da or db =',1pd15.6)
  izd = 1
elseif (t2 .gt. tmx2) then
  if (idb .ge. 2) write (6, 3) it, t2
3 format ('Iteration',i7,3x,'iterdp: Very large value in da or db =',1pd15.6)
  izd = 2
  return
endif

!   Compare the dy vector with those of recent iterations.  If a duplicate is
!   found, then the next iteration must be performed with mq = 1.

do j = 1, nsq
  t1 = 0.d0

  do i = 1, n
    t1 = max (t1, abs (dy(i) - dsyq(i,j)))
  enddo

  if (t1 .le. deps) then
    if (idb .ge. 2) write (6, 4) it, j
4   format ('Iteration',i7,3x,'iterdp: Duplicate found, j =',i6)
    imq = 1
    goto 120
  endif
enddo

!   Place the vector dy in the table dsyq.

120   continue
k = 1 + mod (it, nsq)

do i = 1, n
  dsyq(i,k) = dy(i)
enddo

if (idb .ge. 3) then
  write (6, 5)
5 format ('iterdp: Updated dy:')
  call matoutdp (1, n, dy)
  write (6, 6)
6 format ('iterdp: Updated da matrix:')
  call matoutdp (n, n, da)
  write (6, 7)
7 format ('iterdp: Updated db matrix:')
  call matoutdp (n, n, db)
  write (6, 8)
8 format ('iterdp: Updated dh matrix:')
  call matoutdp (n, n - 1, dh)
endif

return
end

subroutine lqdp (n, m, dh)

!   This performs an LQ decomposition on the DP matrix dh.  It is a simplified
!   and transposed adaptation of the subroutine dqrdc from Linpack.

implicit none
integer i, j, jp, l, lup, m, ml, n
double precision dh(n,m), nrmxl, one, t, zero

zero = 0.d0
one = 1.d0
lup = min (m,n)

!   Perform the householder reduction of dh.

do l = 1, lup
  if (l .eq. m) go to 280

!   Compute the householder transformation for column l.

  ml = m - l
  t = zero

  do i = 0, ml
    t = t + dh(l,l+i) ** 2
  enddo

  nrmxl = sqrt (t)
  if (nrmxl .eq. zero) go to 270
  if (dh(l,l) .ne. zero) nrmxl = sign (nrmxl, dh(l,l))
  t = one / nrmxl

  do i = 0, ml
    dh(l,l+i) = t * dh(l,l+i)
  enddo

  dh(l,l) = one + dh(l,l)

!   Apply the transformation to the remaining columns, updating the norms.

  do j = l + 1, n
    t = zero

    do i = 0, ml
      t = t + dh(l,l+i) * dh(j,l+i)
    enddo

    t = - t / dh(l,l)

    do i = 0, ml
      dh(j,l+i) = dh(j,l+i) + t * dh(l,l+i)
    enddo
  enddo

!   Save the transformation.

  dh(l,l) = - nrmxl
270 continue
280 continue
enddo

!   Zero dh above the diagonal.

do j = 1, m
  do i = 1, j - 1
    dh(i,j) = 0.d0
  enddo
enddo

return
end

subroutine lqmp (n, m, h)

!   This performs an LQ decomposition on the MP matrix h.  It is a simplified
!   and transposed adaptation of the subroutine dqrdc from Linpack.

use mpmodule
implicit none
integer i, j, jp, l, lup, m, ml, n
type (mp_real) h(n,m), nrmxl, one, t, zero

zero = 0.d0
one = 1.d0
lup = min (m,n)

!   Perform the householder reduction of h.

do l = 1, lup
  if (l .eq. m) go to 280

!   Compute the householder transformation for column l.

  ml = m - l
  t = zero

  do i = 0, ml
    t = t + h(l,l+i) ** 2
  enddo

  nrmxl = sqrt (t)
  if (nrmxl .eq. zero) go to 270
  if (h(l,l) .ne. zero) nrmxl = sign (nrmxl, h(l,l))
  t = one / nrmxl

  do i = 0, ml
    h(l,l+i) = t * h(l,l+i)
  enddo

  h(l,l) = one + h(l,l)

!   Apply the transformation to the remaining columns, updating the norms.

!$omp parallel do private (i, j, t) copyin (mp5)
  do j = l + 1, n
    t = zero

    do i = 0, ml
      t = t + h(l,l+i) * h(j,l+i)
    enddo

    t = - t / h(l,l)

    do i = 0, ml
      h(j,l+i) = h(j,l+i) + t * h(l,l+i)
    enddo
  enddo
!$omp end parallel do

!   Save the transformation.

  h(l,l) = - nrmxl
270 continue
280 continue
enddo

!   Zero h above the diagonal.

do j = 1, m
  do i = 1, j - 1
    h(i,j) = 0.d0
  enddo
enddo

return
end

subroutine savedp (n, da, db, dh, dy, dsa, dsb, dsh, dsy)

!   This saves the arrays dy, da, db, dh in case dp iterations must be aborted.
!   A call to the same routine, with (da,db,dh,dy) and (dsa,dsb,dsh,dsy)
!   exchanged, serves to restore these arrays.

implicit none
integer i, j, n
double precision da(n,n), db(n,n), dh(n,n), dy(n), dsa(n,n), dsb(n,n), &
  dsh(n,n), dsy(n)

do i = 1, n
  dsy(i) = dy(i)
enddo

do j = 1, n
  do i = 1, n
    dsa(i,j) = da(i,j)
    dsb(i,j) = db(i,j)
  enddo
enddo

do j = 1, n - 1
  do i = 1, n
    dsh(i,j) = dh(i,j)
  enddo
enddo

return
end

subroutine updtmp (idb, it, n, ix, dx, da, db, b, h, t, y, izm)

!   This updates the MP arrays from the DP arrays.

use mpmodule
use globdata
implicit none 
integer i, i1, idb, it, izm, j, n, n1, n2, ntl
parameter (ntl = 72)
integer ix(n)
double precision da(n,n), db(n,n), dx(n), d1, d2
type (mp_real) b(n,n), h(n,n), t(n), y(n), t1, t2, eps, teps

if (idb .ge. 3) write (6, 1) it
1 format ('Iteration',i7,3x,'updtmp: MP update from DP arrays.')
eps = mpreal (10.d0) ** (nepsilon1+10)
teps = 2.d0 ** ntl * eps
izm = 0
t1 = 1.d300
t2 = 0.d0

!   Update y with db.

call mxmdm (n, 1, db, y, t)

do i = 1, n
  if (abs (y(i)) < t1) then
    i1 = i
    t1 = abs (y(i))
  endif
  t2 = max (t2, abs (y(i)))
enddo

if (idb .ge. 2) then
  call decmd (t1, d1, n1)
  call decmd (t2, d2, n2)
  write (6, 2) it, d1, n1, d2, n2
2  format ('Iteration',i7,3x,'updtmp: Min, max of y =',f11.6,'D',i5, &
     f11.6,'D',i5)
endif

!   Update b with db.

call mxmdm (n, n, db, b, t)

!   Update h with da.  There is no need to perform a LQ decomposition on h.

call mxmdm (n, n - 1, da, h, t)

!   Find the largest entry of b in the same row as the smallest y.

t2 = 0.d0

do i = 1, n
  t2 = max (t2, abs (b(i1,i)))
enddo

if (t1 .le. t2 * teps) then
  if (idb .ge. 2) then
    call decmd (t1, d1, n1)
    write (6, 3) it, d1, n1
3   format ('Iteration',i7,3x,'updtmp: Small value in y =',f11.6,'D',i5)
  endif
  if (t1 .le. t2 * eps) then
    izm = 1
  else
    if (idb .ge. 1) write (6, 4) it
4   format ('Iteration',i7,3x,'updtmp: Precision exhausted.'/ &
      'Try increasing epsilon to a less negative value.')
    izm = 2
  endif
endif

if (idb .ge. 3) then
  write (6, 5)
5 format ('updtmp: Updated y:')
  call matoutmd (1, n, ix, dx, y)
  write (6, 6)
6 format ('updtmp: Updated b matrix:')
  call matoutmd (n, n, ix, dx, b)
  write (6, 7)
7 format ('updtmp: Updated h matrix:')
  call matoutmd (n, n - 1, ix, dx, h)
endif

return
end

!------------------------------

!   Second- and third-level subroutines.

function bounddp (n, dh)

!   This computes the norm bound using DP arithmetic.

implicit none
integer i, n
double precision bounddp, t1, t2
double precision dh(n,n)

t1 = 0.d0

do i = 1, n - 1
  t1 = max (t1, abs (dh(i,i)))
enddo

if (t1 < 1.d-300) then
  write (6, 1)
1 format ('bound: h matrix too small -- use 1-level or 3-level pslqm program.')
  bounddp = -1.d0
else
  bounddp = 1.d0 / t1
endif

return
end

subroutine matoutdp (n1, n2, a)

!   This outputs the DP matrix a.

implicit none
integer i, j, n1, n2
double precision a(n1,n2)

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)
  write (6, 2) (a(i,j), j = 1, n2)
2 format (1p5d15.5)
enddo

return
end

subroutine mxmdm (n1, n2, a, b, c)

!  This multiplies the DP square matrix a by the MP matrix b, and the result
!  is placed in b.  c is a MP scratch array with n1 cells.  n1, n2 are
!  the matrix dimensions as indicated below.

use mpmodule
implicit none
integer i, j, n1, n2
double precision a(n1,n1)
type (mp_real) b(n1,n2), c(n1), t1, t2

!$omp parallel do private (i, j, c, t1) copyin (mp5)
do j = 1, n2
  do i = 1, n1
    call mpdotd (n1, 1, b(1,j), n1, a(i,1), t1)
    c(i) = t1
  enddo

  do i = 1, n1
    b(i,j) = c(i)
  enddo
enddo
!$omp end parallel do

return
end

subroutine qsortdp (n, a, ip)

!   This routine sorts the entries of the N-long DP vector A into ascending
!   order using the quicksort algorithm.  The permutation vector that would
!   sort the vector is returned in IP.

implicit none
integer i, iq, it, j, jq, jz, k, l, n
integer ip(n), ik(50), jk(50)
double precision a(n), s0, s1, s2

do i = 1, n
  ip(i) = i
enddo

if (n .eq. 1) return

k = 1
ik(1) = 1
jk(1) = n

130 i = ik(k)
j = jk(k)
iq = i
jq = j
it = (i + j + 1) / 2
l = ip(j)
ip(j) = ip(it)
ip(it) = l
s0 = a(ip(j))
j = j - 1

140 continue

do l = i, j
  if (s0 .lt. a(ip(l))) goto 160
enddo

i = j
goto 190

160 i = l

do l = j, i, -1
  if (s0 .gt. a(ip(l))) goto 180
enddo

j = i
goto 190

180 j = l
if (i .ge. j) goto 190
l = ip(i)
ip(i) = ip(j)
ip(j) = l
goto 140

190 continue
if (s0 .ge. a(ip(i))) goto 200
l = ip(jq)
ip(jq) = ip(i)
ip(i) = l

200 k = k - 1
jz = 0
if (j .eq. iq) goto 210
k = k + 1
jk(k) = j
jz = 1

210 i = i + 1
if (i .eq. jq) goto 220
k = k + 1
ik(k) = i
jk(k) = jq
if (jz .eq. 0) goto 220
if (j - iq .ge. jq - i) goto 220
ik(k-1) = i
jk(k-1) = jq
ik(k) = iq
jk(k) = j

220 if (k .gt. 0) goto 130

return
end

