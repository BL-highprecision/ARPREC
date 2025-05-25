!  mp_modx.f

!  This Fortran-90 code provides a basic translation facility for an "mp_realx"
!  datatype, which is an extra-high precision (typically double-long) variant of 
!  the mp_real datatype.  The four arithmetic operations are defined between 
!  two mp_realx variables, between a mp_realx variable and a double precision 
!  variable, and between a mp_realx variable and a mp_real variable.  
!  Comparison operations between two mp_realx variables, and a few basic 
!  initrinsics are also defined here.  This satisfies the needs of the F-90
!  tquaderf and tquadts codes only -- it is not a complete package.

!  Note that these routines do NOT automatically change the working precision
!  level to extra-high precision, because of the overhead of saving, setting and
!  restoring precision with each individual call.  Instead, this is done in the
!  user program at the beginning and end of a section of extra-high precision
!  computation, by using calls to mpsetprec or mpsetprecwords.  See the programs
!  tquaderf.f and tquadts.f for some examples of this usage.

!  If order to use these routines, subroutine mpinitx (see below) must be
!  called immediately after the call to mpinit in the user's main program.

!   David H Bailey    2010-08-27

!  This work was supported by the Director, Office of Science, Division
!  of Mathematical, Information, and Computational Sciences of the
!  U.S. Department of Energy under contract number DE-AC03-76SF00098.

module mpmodulex
use mpmodule
implicit none

!   mpiplx is the maximum precision level, in digits.

integer mpiplx
parameter (mpiplx = 4000)

!   NOTE:  This code should not be changed below this point.

integer mpwdsx, mpworkx5
parameter (mpwdsx = (mpiplx-1) / digits_per_word + 2.d0)

type mp_realx
  sequence
  real*8 mpr(mpwdsx+5)
end type

type (mp_realx), public:: mpl02x, mpl10x, mppicx, mpepsx

! mp_realx operator extension interface blocks.
  
  interface operator (+)
      module procedure mpadd_xx
      module procedure mpadd_xq
      module procedure mpadd_qx
      module procedure mpadd_xd
      module procedure mpadd_dx
  end interface

  interface operator (-)
      module procedure mpsub_xx
      module procedure mpsub_xq
      module procedure mpsub_qx
      module procedure mpsub_xd
      module procedure mpsub_dx
! negation
      module procedure mpneg_x
  end interface

  interface operator (*)
      module procedure mpmul_xx
      module procedure mpmul_xq
      module procedure mpmul_qx
      module procedure mpmul_xd
      module procedure mpmul_dx
  end interface

  interface operator (/)
      module procedure mpdiv_xx
      module procedure mpdiv_xq
      module procedure mpdiv_qx
      module procedure mpdiv_xd
      module procedure mpdiv_dx
  end interface

  interface operator (**)
      module procedure mpexp_xi
  end interface

  interface assignment (=)
      module procedure mpeq_xx
      module procedure mpeq_xq
      module procedure mpeq_qx
      module procedure mpeq_xd
      module procedure mpeq_dx
  end interface

  interface operator (.eq.)
      module procedure mpeqt_xx
  end interface

  interface operator (.ne.)
      module procedure mpnet_xx
  end interface

  interface operator (.le.)
      module procedure mplet_xx
  end interface

  interface operator (.ge.)
      module procedure mpget_xx
  end interface

  interface operator (.lt.)
      module procedure mpltt_xx
  end interface

  interface operator (.gt.)
      module procedure mpgtt_xx
  end interface

  interface abs
    module procedure mp_absx
  end interface

  interface atan
    module procedure mp_atanx
  end interface

  interface dble
    module procedure mp_xtod
  end interface

  interface exp
    module procedure mp_expx
  end interface

  interface max
    module procedure mp_maxx
  end interface

  interface min
    module procedure mp_minx
  end interface

  interface mpread
    module procedure mp_inpx
  end interface

  interface mpwrite
    module procedure mp_outx
  end interface

  interface mpreal
    module procedure mp_xtoq
  end interface

  interface mprealx
    module procedure mp_dtox
    module procedure mp_qtox
  end interface

  interface sqrt
    module procedure mp_sqrtx
  end interface
contains

  subroutine mpinitx (n_mpiplx)
!  MPINIT must be called at the start of execution in the user's main program,
!  immediately after the call to mpinit.  It sets the extra-high precision level
!  and epsilon level, and computes extra-high precision versions of the 
!  constants Pi, Log(2) and Log(10).

!  The arguments are as follows:
!  n_mpiplx: integer, optional.  Extra-high working precision level, in digits.
!     Must not exceed mpiplx, which is set near the start of this file.
!     If n_mpiplx is not present, then precision level is set to mpiplx.

      implicit none
      integer, intent(in), optional :: n_mpiplx
      integer iconst_temp

      if (present (n_mpiplx)) then
         if (n_mpiplx .gt. mpiplx) then
            write(mpldb, *) '*** MPINITX: new precision level is too high'
            stop
         endif
         iconst_temp = 1
         call f_mpinit (n_mpiplx, iconst_temp, mpwdsx, mpworkx5, &
           mpepsx%mpr, mpl02x%mpr, mpl10x%mpr, mppicx%mpr)
         new_mpipl = n_mpiplx
         new_mpwork = mpworkx5 - 5
      else
         iconst_temp = 1
         call f_mpinit (mpiplx, iconst_temp, mpwdsx, mpworkx5, &
           mpepsx%mpr, mpl02x%mpr, mpl10x%mpr, mppicx%mpr)
         new_mpipl = mpiplx
         new_mpwork = mpworkx5 - 5
      endif
  end subroutine

! Additions
  type (mp_realx) function mpadd_xx (a, b)
      type (mp_realx), intent(in) :: a, b
      mpadd_xx%mpr(1) = mpworkx5
      call f_mpadd (a%mpr, b%mpr, mpadd_xx%mpr)
  end function

  type (mp_realx) function mpadd_xq (a, b)
      type (mp_realx), intent(in) :: a
      type (mp_real), intent(in) :: b
      mpadd_xq%mpr(1) = mpworkx5
      call f_mpadd (a%mpr, b%mpr, mpadd_xq%mpr)
  end function

  type (mp_realx) function mpadd_qx (a, b)
      type (mp_real), intent(in) :: a
      type (mp_realx), intent(in) :: b
      mpadd_qx%mpr(1) = mpworkx5
      call f_mpadd (a%mpr, b%mpr, mpadd_qx%mpr)
  end function

  type (mp_realx) function mpadd_xd (a, b)
      type (mp_realx), intent(in) :: a
      real*8, intent(in) :: b
      mpadd_xd%mpr(1) = mpworkx5
      call f_mpadd_d (a%mpr, b, mpadd_xd%mpr)
  end function

  type (mp_realx) function mpadd_dx (a, b)
      real*8, intent(in) :: a
      type (mp_realx), intent(in) :: b
      mpadd_dx%mpr(1) = mpworkx5
      call f_mpadd_d (b%mpr, a, mpadd_dx%mpr)
  end function

! Subtractions
  type (mp_realx) function mpsub_xx (a, b)
      type (mp_realx), intent(in) :: a, b
      mpsub_xx%mpr(1) = mpworkx5
      call f_mpsub (a%mpr, b%mpr, mpsub_xx%mpr)
  end function

  type (mp_realx) function mpsub_xq (a, b)
      type (mp_realx), intent(in) :: a
      type (mp_real), intent(in) :: b
      mpsub_xq%mpr(1) = mpworkx5
      call f_mpsub (a%mpr, b%mpr, mpsub_xq%mpr)
  end function

  type (mp_realx) function mpsub_qx (a, b)
      type (mp_real), intent(in) :: a
      type (mp_realx), intent(in) :: b
      mpsub_qx%mpr(1) = mpworkx5
      call f_mpsub (a%mpr, b%mpr, mpsub_qx%mpr)
  end function

  type (mp_realx) function mpsub_xd (a, b)
      type (mp_realx), intent(in) :: a
      real*8, intent(in) :: b
      mpsub_xd%mpr(1) = mpworkx5
      call f_mpsub_d (a%mpr, b, mpsub_xd%mpr)
  end function

  type (mp_realx) function mpsub_dx (a, b)
      real*8, intent(in) :: a
      type (mp_realx), intent(in) :: b
      mpsub_dx%mpr(1) = mpworkx5
      call f_mpsub_dq (a, b%mpr, mpsub_dx%mpr)
  end function

! Unary Minus
  type (mp_realx) function mpneg_x (a)
    type (mp_realx), intent(in) :: a
    mpneg_x%mpr(1) = mpworkx5
    call f_mpneg_q (a%mpr, mpneg_x%mpr);
  end function

! Multiplications
  type (mp_realx) function mpmul_xx (a, b)
      type (mp_realx), intent(in) :: a, b
      mpmul_xx%mpr(1) = mpworkx5
      call f_mpmul (a%mpr, b%mpr, mpmul_xx%mpr)
  end function

  type (mp_realx) function mpmul_xq (a, b)
      type (mp_realx), intent(in) :: a
      type (mp_real), intent(in) :: b
      mpmul_xq%mpr(1) = mpworkx5
      call f_mpmul (a%mpr, b%mpr, mpmul_xq%mpr)
  end function

  type (mp_realx) function mpmul_qx (a, b)
      type (mp_real), intent(in) :: a
      type (mp_realx), intent(in) :: b
      mpmul_qx%mpr(1) = mpworkx5
      call f_mpmul (a%mpr, b%mpr, mpmul_qx%mpr)
  end function

  type (mp_realx) function mpmul_xd (a, b)
      type (mp_realx), intent(in) :: a
      real*8, intent(in) :: b
      mpmul_xd%mpr(1) = mpworkx5
      call f_mpmul_qd (a%mpr, b, mpmul_xd%mpr)
  end function

  type (mp_realx) function mpmul_dx (a, b)
      real*8, intent(in) :: a
      type (mp_realx), intent(in) :: b
      mpmul_dx%mpr(1) = mpworkx5
      call f_mpmul_qd (b%mpr, a, mpmul_dx%mpr)
  end function

! Divisions
  type (mp_realx) function mpdiv_xx (a, b)
      type (mp_realx), intent(in) :: a, b
      mpdiv_xx%mpr(1) = mpworkx5
      call f_mpdiv (a%mpr, b%mpr, mpdiv_xx%mpr)
  end function

  type (mp_realx) function mpdiv_xq (a, b)
      type (mp_realx), intent(in) :: a
      type (mp_real), intent(in) :: b
      mpdiv_xq%mpr(1) = mpworkx5
      call f_mpdiv (a%mpr, b%mpr, mpdiv_xq%mpr)
  end function

  type (mp_realx) function mpdiv_qx (a, b)
      type (mp_real), intent(in) :: a
      type (mp_realx), intent(in) :: b
      mpdiv_qx%mpr(1) = mpworkx5
      call f_mpdiv (a%mpr, b%mpr, mpdiv_qx%mpr)
  end function

  type (mp_realx) function mpdiv_xd (a, b)
      type (mp_realx), intent(in) :: a
      real*8, intent(in) :: b
      mpdiv_xd%mpr(1) = mpworkx5
      call f_mpdiv_qd (a%mpr, b, mpdiv_xd%mpr)
  end function

  type (mp_realx) function mpdiv_dx (a, b)
      real*8, intent(in) :: a
      type (mp_realx), intent(in) :: b
      mpdiv_dx%mpr(1) = mpworkx5
      call f_mpdiv_dq (a, b%mpr, mpdiv_dx%mpr)
  end function

! Powers
  type (mp_realx) function mpexp_xi (a, ib)
      type (mp_realx), intent(in) :: a
      integer, intent(in) :: ib 
      mpexp_xi%mpr(1) = mpworkx5
      call f_mppwr_qi (a%mpr, ib, mpexp_xi%mpr)
  end function

! Assignments
  subroutine mpeq_xx (a, b)
      type (mp_realx), intent(inout) :: a
      type (mp_realx), intent(in) :: b
      a%mpr(1) = mpworkx5
      call f_mpeq (b%mpr, a%mpr)
  end subroutine

  subroutine mpeq_xq (a, b)
      type (mp_realx), intent(inout) :: a
      type (mp_real), intent(in) :: b
      a%mpr(1) = mpworkx5
      call f_mpeq (b%mpr, a%mpr)
  end subroutine

  subroutine mpeq_qx (a, b)
      type (mp_real), intent(inout) :: a
      type (mp_realx), intent(in) :: b
      a%mpr(1) = mpwork5
      call f_mpeq (b%mpr, a%mpr)
  end subroutine

  subroutine mpeq_xd (a, b)
      type (mp_realx), intent(inout) :: a
      real*8, intent(in) :: b
      a%mpr(1) = mpworkx5
      call f_mpeq_d (b, a%mpr)
  end subroutine

  subroutine mpeq_dx (a, b)
      real*8, intent(out) :: a
      type (mp_realx), intent(in) :: b
      double precision db
      integer ib
      call f_mpmdc (b%mpr, db, ib)
      a = db * 2.d0 ** ib
  end subroutine

! Equality
  logical function mpeqt_xx (a, b)
      type (mp_realx), intent(in) :: a, b
      integer ic
      call f_mpcpr (a%mpr, b%mpr, ic)
      if (ic .eq. 1) then
         mpeqt_xx = .true.
      else
         mpeqt_xx = .false.
      endif
      return
  end function

! Inequality
  logical function mpnet_xx (a, b)
      type (mp_realx), intent(in) :: a, b
      integer ic
      call f_mpcpr (a%mpr, b%mpr, ic)
      if (ic .ne. 1) then
         mpnet_xx = .true.
      else
         mpnet_xx = .false.
      endif
      return
  end function
      
! Less-Than-Or-Equal-To
  logical function mplet_xx (a, b)
      type (mp_realx), intent(in) :: a, b
      integer ic
      call f_mplet (a%mpr, b%mpr, ic)
      if (ic .eq. 1) then
         mplet_xx = .true.
      else
         mplet_xx = .false.
      endif
      return
  end function
      
! Greater-Than-Or-Equal-To
  logical function mpget_xx (a, b)
      type (mp_realx), intent(in) :: a, b
      integer ic
      call f_mpget (a%mpr, b%mpr, ic)
      if (ic .eq. 1) then
         mpget_xx = .true.
      else
         mpget_xx = .false.
      endif
      return
  end function

! Less-Than
  logical function mpltt_xx (a, b)
      type (mp_realx), intent(in) :: a, b
      integer ic
      call f_mpltt (a%mpr, b%mpr, ic)
      if (ic .eq. 1) then
         mpltt_xx = .true.
      else
         mpltt_xx = .false.
      endif
      return
  end function
      
! Greater-Than
  logical function mpgtt_xx (a, b)
      type (mp_realx), intent(in) :: a, b
      integer ic
      call f_mpgtt (a%mpr, b%mpr, ic)
      if (ic .eq. 1) then
         mpgtt_xx = .true.
      else
         mpgtt_xx = .false.
      endif
      return
  end function
      
! Absolute value
  type (mp_realx) function mp_absx (a)
    type (mp_realx), intent(in) :: a
    mp_absx%mpr(1) = mpworkx5
    call f_mpabs (a%mpr, mp_absx%mpr)
  end function

  type (mp_realx) function mp_atanx(a)
    type (mp_realx), intent(in) :: a
    mp_atanx%mpr(1) = mpworkx5
    call f_mpatan(a%mpr, mp_atanx%mpr)
  end function

  type (mp_realx) function mp_expx(a)
    type (mp_realx), intent(in) :: a
    mp_expx%mpr(1) = mpworkx5
    call f_mpexp(a%mpr, mp_expx%mpr)
  end function

  double precision function mp_xtod (a)
      type (mp_realx), intent (in):: a
      call f_mpdble (a%mpr, mp_xtod)
  end function

  type (mp_realx) function mp_maxx (a, b)
      type (mp_realx), intent(in) :: a, b
      integer ic
      mp_maxx%mpr(1) = mpworkx5
      call f_mpget (a%mpr, b%mpr, ic)
      if (ic .eq. 1) then
         call f_mpeq (a%mpr, mp_maxx%mpr)
      else
         call f_mpeq (b%mpr, mp_maxx%mpr)
      endif
  end function

  type (mp_realx) function mp_minx (a, b)
      type (mp_realx), intent(in) :: a, b
      integer ic
      mp_minx%mpr(1) = mpworkx5
      call f_mplet (a%mpr, b%mpr, ic)
      if (ic .eq. 1) then
         call f_mpeq (a%mpr, mp_minx%mpr)
      else
         call f_mpeq (b%mpr, mp_minx%mpr)
      endif
  end function

! SQRT, etc.
  type (mp_realx) function mp_sqrtx (a)
    type (mp_realx), intent(in) :: a
    mp_sqrtx%mpr(1) = mpworkx5
    call f_mpsqrt (a%mpr, mp_sqrtx%mpr)
  end function

!  Conversion
  type (mp_real) function mp_xtoq (a)
      type (mp_realx), intent(in) :: a
      mp_xtoq%mpr(1) = mpwork5
      call f_mpeq (a%mpr, mp_xtoq%mpr)
  end function

  type (mp_realx) function mp_qtox (a)
      type (mp_real), intent(in) :: a
      mp_qtox%mpr(1) = mpworkx5
      call f_mpeq (a%mpr, mp_qtox%mpr)
  end function

  type (mp_realx) function mp_dtox (da)
      real*8, intent(in) :: da
      mp_dtox%mpr(1) = mpworkx5
      call f_mpdmc (da, mp_dtox%mpr)
  end function

! Input

  subroutine mp_inpx (iu, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: iu
    type (mp_realx), intent (out) :: q1, q2, q3, q4, q5, q6, q7, q8, q9
    optional :: q2, q3, q4, q5, q6, q7, q8, q9

    call mpinpx (iu, q1)
    if (present (q2)) call mpinpx (iu, q2)
    if (present (q3)) call mpinpx (iu, q3)
    if (present (q4)) call mpinpx (iu, q4)
    if (present (q5)) call mpinpx (iu, q5)
    if (present (q6)) call mpinpx (iu, q6)
    if (present (q7)) call mpinpx (iu, q7)
    if (present (q8)) call mpinpx (iu, q8)
    if (present (q9)) call mpinpx (iu, q9)
    return
  end subroutine

  subroutine mpinpx (iu, qa)
    integer, intent(in) :: iu
    type (mp_realx), intent(out) :: qa
    character*1 az(new_mpipl+100)
    integer i, l, l1, nn
    character*200 line

    l = 0
    nn = new_mpipl + 100
    qa%mpr(1) = mpworkx5

100  continue
    read(iu, '(a)', end = 200) line

    do i = 200, 1, -1
      if (line(i:i) /= ' ') goto 110
    enddo

    i = 0

110  continue
    l1 = i

    do i = 1, l1
      if (line(i:i) == ',') goto 150
      if (l + 1 <= nn) then
        l = l + 1
        az(l) = line(i:i)
      endif
    enddo

    goto 100

150  continue

    call mpinpc (az, l, qa%mpr)
    goto 300

200  continue

     write (mpldb, 1)
1    format (&
     'mpinpx: end of file or no comma terminating multiprecion input.')
     stop

300  continue
    return
  end subroutine

! Output
  subroutine mp_outx (iu, q1, q2, q3, q4, q5, q6, q7, q8, q9)
      integer, intent(in) :: iu
      type (mp_realx), intent(in) :: q1, q2, q3, q4, q5, q6, q7, q8, q9
      optional :: q2, q3, q4, q5, q6, q7, q8, q9
      call mpoutx (iu, q1)
      if (present (q2)) call mpoutx (iu, q2)
      if (present (q3)) call mpoutx (iu, q3)
      if (present (q4)) call mpoutx (iu, q4)
      if (present (q5)) call mpoutx (iu, q5)
      if (present (q6)) call mpoutx (iu, q6)
      if (present (q7)) call mpoutx (iu, q7)
      if (present (q8)) call mpoutx (iu, q8)
      if (present (q9)) call mpoutx (iu, q9)
      return
  end subroutine

  subroutine mpoutx (iu, q)
      integer, intent(in) :: iu
      type (mp_realx), intent(in) :: q
      character*1 az(new_mpipl+100)
      integer i, l
      call mpoutc (q%mpr, az, l)
      az(l+1) = ','
      write(iu, '(78A1)') (az(i), i = 1, l+1)
  end subroutine

end module mpmodulex
