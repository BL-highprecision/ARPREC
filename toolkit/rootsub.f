subroutine rroot (n, r0, a, nr, r)

!   Finds a single arbitrary-precision real root, beginning at r0.
!   Modified for the toolkit.
!   David H Bailey   2004-05-27

use mpmodule
use globdata
implicit none
integer i, j, k, n, nit, nr, nwp
parameter (nit = 100)
type (mp_real) a(0:n), ad(0:n-1), eps1
type (mp_real) t1, t2, t3, t4, r, r0

!   Compute derivative.

do i = 0, n - 1
  ad(i) = (i + 1) * a(i+1)
enddo

nwp = 3
call mpsetprecwords (nwp)
t1 = r0
t2 = 0.d0
if (r0 /= t2) then
  eps1 = 2.d0 ** (-96) * abs (t1)
else
  eps1 = 2.d0 ** (-96)
endif

!   Perform up to nit double-precision Newton iterations.

do j = 1, nit

!   Evaluate polynomial at t1.

  t2 = 0.d0

  do i = n, 0, -1
    t2 = t1 * t2 + a(i)
  enddo

!   Evaluate polynomial derivative at t1.

  t3 = 0.d0

  do i = n - 1, 0, -1
    t3 = t1 * t3 + ad(i)
  enddo

  t4 = t1 - t2 / t3
  if (nwp == 3) then
    if (abs (t1 - t4) < eps1) then
      nwp = min (2 * nwp - 1, nwords1)
      call mpsetprecwords (nwp)
    endif
  elseif (nwp < nwords1) then
    nwp = min (2 * nwp - 1, nwords1)
    call mpsetprecwords (nwp)
  else
    t1 = t4
    goto 100
  endif
  t1 = t4
enddo

write (6, 1)
1 format ('rroot: failed to find arbitrary precision real root.')
nr = 0
r = 0.d0
goto 110

100 continue

nr = 1
r = t1

110 continue

call mpsetprecwords (nwords1)
return
end

subroutine croot (n, r0, a, nr, r)

!   Finds a single arbitrary-precision complex root, beginning at r0.

use mpmodule
use globdata
implicit none
integer i, j, k, n, nit, nr, nwp
parameter (nit = 100)
type (mp_real) a(0:n), ad(0:n-1), eps1
type(mp_complex) t1, t2, t3, t4, r, r0

!   Compute derivative.

do i = 0, n - 1
  ad(i) = (i + 1) * a(i+1)
enddo

nwp = 3
call mpsetprecwords (nwp)
t1 = r0
t2 = (0.d0, 0.d0)
if (r0 /= t2) then
  eps1 = 2.d0 ** (-96) * abs (t1)
else
  eps1 = 2.d0 ** (-96)
endif

!   Perform up to nit double-precision Newton iterations.

do j = 1, nit

!   Evaluate polynomial at t1.

  t2 = (0.d0, 0.d0)

  do i = n, 0, -1
    t2 = t1 * t2 + a(i)
  enddo

!   Evaluate polynomial derivative at t1.

  t3 = (0.d0, 0.d0)

  do i = n - 1, 0, -1
    t3 = t1 * t3 + ad(i)
  enddo

  t4 = t1 - t2 / t3
  if (nwp == 3) then
    if (abs (t1 - t4) < eps1) then
      nwp = min (2 * nwp - 1, nwords1)
      call mpsetprecwords (nwp)
    endif
  elseif (nwp < nwords1) then
    nwp = min (2 * nwp - 1, nwords1)
    call mpsetprecwords (nwp)
  else
    t1 = t4
    goto 100
  endif
  t1 = t4
enddo

write (6, 1)
1 format ('croot: failed to find multiple precision complex root.')
nr = 0
r = (0.d0, 0.d0)
goto 110

100 continue

nr = 1
r = t1

110 continue

call mpsetprecwords (nwords1)
return
end


