      subroutine f_main

! A simple test of the fortran wrappers

      use mpmodule

      type (mp_real) x, y, q2, q3, q4
      real*8 dx, dy
      real rx
      type (mp_integer) jx, jy
      integer i, j
      type (mp_complex) z1, z2, z3
      complex (kind(1.d0)) x1, x2
!
!      double precision, intrinsic :: log10

      print*, "before mpinit"
      call mpinit(100)
      print*, "after mpinit"
      call mpgetpar('mpndb', i)
      print*, 'MPNDB = ', i
      call mpsetpar('mpndb', 3)
      call mpgetpar('mpndb', i)
      print*, 'MPNDB = ', i

!      x = qdpi()
      dx = 4.0
      rx = 4.3
      x = mpreal(-4.6)
      x = aint(x)
      print*, 'x = '
      call mpwrite(6, x)
      dy = 3.0
      y = dy

      q2 = sinh(y)
      dx = sinh(dy)
      q4 = q2 - dx
      q4 = (y - dy * y)  + x
      q4 = sinh(y)  - sinh(dy)  ! cause seg fault
      q4 = log10(y) - log10(dy) ! cause seg fault, also with +, *
      q3 = log10(y)
      dx = log10(dy)
      q3 = mod(x,y)     ! dies at the place with "new"

      i = 2
      jx = i
      j = 30
      jy = j

      x1 = cmplx(1.0, 2.0, kind(1.d0))
      x2 = cmplx(1.0, 2.0, kind(1.d0))
      z1 = x1
      z2 = x2
      z3 = sin(z1) - sin(x1)   ! cause seg fault
      z3 = log(z2) - log(x2)   ! cause seg fault
      z3 = sqrt(z2)
      z3 = x - z1

      print*, 'x1 = ', x1
      print*, 'z1 ='
      call mpwrite(6, z1)
      print*, 'z3 ='
      call mpwrite(6, z3)
      print*, 'jx = '
      call mpwrite(6, jx)
      print*, 'q3 = '
      call mpwrite(6, q3)
      print*, 'q4 = '
      call mpwrite(6, q4)

      print*, 'mpl02 ='
      call mpwrite(6, mpl02)
      print*, 'mpl10 ='
      call mpwrite(6, mpl10)
      print*, 'mppic ='
      call mpwrite(6, mppic)
      print*, 'mpeps ='
      call mpwrite(6, mpeps)
      print*, 'mplrg ='
      call mpwrite(6, mplrg)
      print*, 'mpsml ='
      call mpwrite(6, mpsml)

      end
