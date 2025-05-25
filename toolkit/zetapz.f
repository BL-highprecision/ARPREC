!---------------------------------------------

!   This code evaluates the zetaz and zetap functions.
!   David H Bailey    2004-07-14

function zetaz (n, izeta)

!   This evaluate zeta sums via Formula 44 in the paper "Special Values of
!   Multidimensional Polylogarithms", by J. M. Borwein et al, CECM 98-106.

use mpmodule
use globdata
implicit none
integer i, ii, j, k, ka, kb, kt1, kt2, m, n, nx
parameter (nx = 100)
integer ia(nx,0:nx+1), ial(0:nx+1), ib(nx,0:nx+1), ibl(0:nx+1), &
  ir(nx), is(nx), it1(nx), it2(nx), izeta(n)
type (mp_real) zetap, t1, t2, t3, zetaz
external zetap

!   Form the ir and is arrays.

k = 0

do m = 1, n
  k = k + 1
  if (k > n) goto 110
  if (izeta(k) <= 1) then
    nerror = 301
    goto 400
  endif
  ir(m) = 0
  is(m) = izeta(k) - 2

100 continue

  k = k + 1
  if (k > n) goto 110
  if (izeta(k) == 1) then
    ir(m) = ir(m) + 1
    goto 100
  else
    k = k - 1
  endif
enddo

nerror = 302
goto 400

110 continue

! write (6, *) 'm =', m
! write (6, *) 'ir array:'
! write (6, '(3i6)') (i, ir(i), i = 1, m)
! write (6, *) 'is array:'
! write (6, '(3i6)') (i, is(i), i = 1, m)

!   Form the ia and ib arrays.

t1 = 0.d0
ial(0) = 0
ial(m+1) = 0
ibl(0) = 0
ibl(m+1) = 0

do j = 1, m
  ka = 0

  do i = j, m
    ka = ka + 1
    if (ka > nx) then
      nerror = 303
      goto 400
    endif
    ia(ka,j) = is(i) + 2

    do ii = 1, ir(i)
      ka = ka + 1
      if (ka > nx) then
        nerror = 304
        goto 400
      endif
      ia(ka,j) = 1
    enddo
  enddo

  ial(j) = ka
  kb = 0

  do i = j, 1, -1
    kb = kb + 1
    if (kb > nx) then
      nerror = 305
      goto 400
    endif
    ib(kb,j) = ir(i) + 2

    do ii = 1, is(i)
      kb = kb + 1
      if (kb > nx) then
        nerror = 306
        goto 400
      endif
      ib(kb,j) = 1
    enddo
  enddo

  ibl(j) = kb
enddo

! write (6, *) 'ia array:'
! write (6, '(3i6)') ((i, j, ia(i,j), i = 1, ial(j)), j = 1, m)
! write (6, *) 'ib array:'
! write (6, '(3i6)') ((i, j, ib(i,j), i = 1, ibl(j)), j = 1, m)

!   Compute the sum.

t3 = 0.d0

do j = 1, m
  do i = 0, is(j) + 1

!    write (6, *) 'i, j =', i, j

    kt1 = 1
    it1(1) = is(j) + 2 - i

    do ii = 1, ir(j)
      kt1 = kt1 + 1
      if (kt1 > nx) then
        nerror = 307
        goto 400
      endif
      it1(kt1) = 1
    enddo

    do ii = 1, ial(j+1)
      kt1 = kt1 + 1
      if (kt1 > nx) then
        nerror = 308
        goto 400
      endif
      it1(kt1) = ia(ii,j+1)
    enddo

    kt2 = 0

    do ii = 1, i
      kt2 = kt2 + 1
      if (kt2 > nx) then
        nerror = 309
        goto 400
      endif
      it2(kt2) = 1
    enddo

    do ii = 1, ibl(j-1)
      kt2 = kt2 + 1
      if (kt2 > nx) then
        nerror = 310
        goto 400
      endif
      it2(kt2) = ib(ii,j-1)
    enddo

    t1 = zetap (kt1, it1)
    if (nerror > 0) goto 500
    t2 = zetap (kt2, it2)
    if (nerror > 0) goto 500

!    write (6, *) 'it1 =', (it1(k), k = 1, kt1)
!    write (6, *) 't1 ='
!    call mpwrite (6, t1)
!    write (6, *) 'it2 =', (it2(k), k = 1, kt2)
!    write (6, *) 't2 ='
!    call mpwrite (6, t2)

    t3 = t3 + t1 * t2
  enddo

  do i = 1, ir(j)
    kt1 = 0

    do ii = 1, i
      kt1 = kt1 + 1
      if (kt1 > nx) then
        nerror = 311
        goto 400
      endif
      it1(kt1) = 1
    enddo

    do ii = 1, ial(j+1)
      kt1 = kt1 + 1
      if (kt1 > nx) then
        nerror = 312
        goto 400
      endif
      it1(kt1) = ia(ii,j+1)
    enddo

    kt2 = 1
    it2(kt2) = ir(j) + 2 - i

    do ii = 1, is(j)
      kt2 = kt2 + 1
      if (kt2 > nx) then
        nerror = 313
        goto 400
      endif
      it2(kt2) = 1
    enddo

    do ii = 1, ibl(j-1)
      kt2 = kt2 + 1
      if (kt2 > nx) then
        nerror = 314
        goto 400
      endif
      it2(kt2) = ib(ii,j-1)
    enddo

    t1 = zetap (kt1, it1)
    if (nerror > 0) goto 500
    t2 = zetap (kt2, it2)
    if (nerror > 0) goto 500

!    write (6, *) 'it1 =', (it1(k), k = 1, kt1)
!    write (6, *) 't1 ='
!    call mpwrite (6, t1)
!    write (6, *) 'it2 =', (it2(k), k = 1, kt2)
!    write (6, *) 't2 ='
!    call mpwrite (6, t2)

    t3 = t3 + t1 * t2
  enddo
enddo

kt1 = 0

do ii = 1, ibl(m)
  kt1 = kt1 + 1
  if (kt1 > nx) then
    nerror = 315
    goto 400
  endif
  it1(kt1) = ib(ii,m)
enddo

t1 = zetap (kt1, it1)
if (nerror > 0) goto 500

! write (6, *) 'final evaluation'
! write (6, *) 'it1 =', (it1(k), k = 1, kt1)
! write (6, *) 't1 ='
! call mpwrite (6, t1)

zetaz = t3 + t1

! write (6, *) 'final result ='
! call mpwrite (6, zetaz)

goto 500

400 continue

write (6, 1) 
1 format ('zetaz: program error.')

500 continue

return
end

function zetap (ka, ia)
use mpmodule
use globdata
implicit none
integer i, ka, ia(0:ka-1), j, jmax, ierror
parameter (jmax = 1000000)
type (mp_real) eps, st(0:ka-1), tj, zetap, t1, t2, t3, t4, t5, t6, t7

if (ka == 0) then
  t1 = 1.d0
  goto 100
endif
eps = 0.5d0 ** 24 * mpreal (10.d0) ** nepsilon1

do i = 0, ka - 2
  st(i) = 0.d0
enddo

st(ka-1) = 1.d0
t1 = 0.d0
t2 = 1.d0

do j = 1, jmax
  tj = j
  t2 = 0.5d0 * t2
  t3 = t2 * st(0) / tj ** ia(0)
  t1 = t1 + t3

  do i = 1, ka - 1
    st(i-1) = st(i-1) + st(i) / tj ** ia(i)
    call mpgetpar ('mpier', ierror)
    if (ierror > 0) then
      nerror = 316
      goto 100
    endif
  enddo

  if (abs (t3) < eps .and. j > 20 * ka) goto 100
enddo

100 continue

zetap = t1
return
end
