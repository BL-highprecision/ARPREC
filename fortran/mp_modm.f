!  mp_modm.f

!  This Fortran-90 code provides a basic translation facility for an "mp_realm"
!  datatype, which is a moderate precision (typically 125 digit) variant of 
!  the mp_real datatype.  The four arithmetic operations are defined between 
!  two mp_realm variables, between a mp_realm variable and a double precision 
!  variable, and between a mp_realm variable and a mp_real variable.  
!  Comparison operations between two mp_realm variables, and a few basic 
!  initrinsics are also defined here.  This satisfies the needs of the F-90
!  PSLQ3 and PSLQM3 codes only -- it is not a complete package.  The purpose
!  of these routines is to save both computation time and memory (since this
!  datatype uses only a few words).

!  Note that these routines do NOT automatically change the working precision
!  level to moderate precision, because of the overhead of saving, setting and
!  restoring precision with each individual call.  Instead, this is done in the
!  user program at the beginning and end of a section of moderate precision
!  computation, by using calls to mpsetprec or mpsetprecwords.  See the programs
!  tpslq3.f and tpslqm3.f for some examples of this usage.  Indeed, the working
!  precision should NOT be lowered for many of the routines below (example:
!  mpmul_mq, which performs moderate x full precision multiplication), since
!  such routines are normally used to return full precision results.

!   David H Bailey    2004-04-28

!  This work was supported by the Director, Office of Science, Division
!  of Mathematical, Information, and Computational Sciences of the
!  U.S. Department of Energy under contract number DE-AC03-76SF00098.

module mpmodulem
use mpmodule
implicit none

!   mpiplm is the maximum precision level, in digits.

integer mpiplm
parameter (mpiplm = 125)

!   NOTE:  This code should not be changed below this point.

integer mpwdsm, mpm5
parameter (mpwdsm = (mpiplm-1) / digits_per_word + 2.d0, mpm5 = mpwdsm + 5)
character*1, private :: az(mpiplm+100)

type mp_realm
  sequence
  real*8 mpr(mpm5)
end type

! mp_realm operator extension interface blocks.
  
  interface operator (+)
      module procedure mpadd_mm
      module procedure mpadd_mq
      module procedure mpadd_qm
      module procedure mpadd_md
      module procedure mpadd_dm
  end interface

  interface operator (-)
      module procedure mpsub_mm
      module procedure mpsub_mq
      module procedure mpsub_qm
      module procedure mpsub_md
      module procedure mpsub_dm
! negation
      module procedure mpneg_m
  end interface

  interface operator (*)
      module procedure mpmul_mm
      module procedure mpmul_mq
      module procedure mpmul_qm
      module procedure mpmul_md
      module procedure mpmul_dm
  end interface

  interface operator (/)
      module procedure mpdiv_mm
      module procedure mpdiv_mq
      module procedure mpdiv_qm
      module procedure mpdiv_md
      module procedure mpdiv_dm
  end interface

  interface operator (**)
      module procedure mpexp_mi
  end interface

  interface assignment (=)
      module procedure mpeq_mm
      module procedure mpeq_mq
      module procedure mpeq_qm
      module procedure mpeq_md
      module procedure mpeq_dm
  end interface

  interface operator (.eq.)
      module procedure mpeqt_mm
  end interface

  interface operator (.ne.)
      module procedure mpnet_mm
  end interface

  interface operator (.le.)
      module procedure mplet_mm
  end interface

  interface operator (.ge.)
      module procedure mpget_mm
  end interface

  interface operator (.lt.)
      module procedure mpltt_mm
  end interface

  interface operator (.gt.)
      module procedure mpgtt_mm
  end interface

  interface abs
    module procedure mp_absm
  end interface

  interface aint
     module procedure mp_aintm
  end interface

  interface anint
    module procedure mp_anintm
  end interface

  interface dble
    module procedure mp_mtod
  end interface

  interface max
    module procedure mp_maxm
    module procedure mp_maxm3
  end interface

  interface min
    module procedure mp_minm
    module procedure mp_minm3
  end interface

  interface mpread
    module procedure mp_inpm
  end interface

  interface mpwrite
    module procedure mp_outm
  end interface

  interface mpreal
    module procedure mp_mtoq
  end interface

  interface mprealm
    module procedure mp_dtom
    module procedure mp_qtom
  end interface

  interface sign
    module procedure mp_signm
  end interface

  interface sqrt
    module procedure mp_sqrtm
  end interface
contains

  subroutine mpdotdm (n, isa, a, isb, db, c)
!   This routine computes the dot product of the MPM vector A with the DP
!   vector DB, returning the MPM result in C.  This routine is used in the
!   author's customized PSLQ routine, resulting in substantial speedup.
!   The length of both the A and DB vectors is N, and ISA and ISB are the 
!   skip distances between successive elements of A and DB, measured in 
!   MPM words and DP words, respectively.  The DP values in DB must be
!   whole numbers, so for example they cannot be larger than 2^53.

      integer n, isa, isb
      double precision db(isb*n)
      type(mp_realm) a(isa*n), c
      c%mpr(1) = mpm5
      call f_mpdotd (n, isa * (mpwdsm + 5), a(1)%mpr, isb, db, c%mpr)
  end subroutine
      
! Additions
  type (mp_realm) function mpadd_mm (a, b)
      type (mp_realm), intent(in) :: a, b
      mpadd_mm%mpr(1) = mpm5
      call f_mpadd (a%mpr, b%mpr, mpadd_mm%mpr)
  end function

  type (mp_real) function mpadd_mq (a, b)
      type (mp_realm), intent(in) :: a
      type (mp_real), intent(in) :: b
      mpadd_mq%mpr(1) = mpwork5
      call f_mpadd (a%mpr, b%mpr, mpadd_mq%mpr)
  end function

  type (mp_real) function mpadd_qm (a, b)
      type (mp_real), intent(in) :: a
      type (mp_realm), intent(in) :: b
      mpadd_qm%mpr(1) = mpwork5
      call f_mpadd (a%mpr, b%mpr, mpadd_qm%mpr)
  end function

  type (mp_realm) function mpadd_md (a, b)
      type (mp_realm), intent(in) :: a
      real*8, intent(in) :: b
      mpadd_md%mpr(1) = mpm5
      call f_mpadd_d (a%mpr, b, mpadd_md%mpr)
  end function

  type (mp_realm) function mpadd_dm (a, b)
      real*8, intent(in) :: a
      type (mp_realm), intent(in) :: b
      mpadd_dm%mpr(1) = mpm5
      call f_mpadd_d (b%mpr, a, mpadd_dm%mpr)
  end function

! Subtractions
  type (mp_realm) function mpsub_mm (a, b)
      type (mp_realm), intent(in) :: a, b
      mpsub_mm%mpr(1) = mpm5
      call f_mpsub (a%mpr, b%mpr, mpsub_mm%mpr)
  end function

  type (mp_real) function mpsub_mq (a, b)
      type (mp_realm), intent(in) :: a
      type (mp_real), intent(in) :: b
      mpsub_mq%mpr(1) = mpwork5
      call f_mpsub (a%mpr, b%mpr, mpsub_mq%mpr)
  end function

  type (mp_real) function mpsub_qm (a, b)
      type (mp_real), intent(in) :: a
      type (mp_realm), intent(in) :: b
      mpsub_qm%mpr(1) = mpwork5
      call f_mpsub (a%mpr, b%mpr, mpsub_qm%mpr)
  end function

  type (mp_realm) function mpsub_md (a, b)
      type (mp_realm), intent(in) :: a
      real*8, intent(in) :: b
      mpsub_md%mpr(1) = mpm5
      call f_mpsub_d (a%mpr, b, mpsub_md%mpr)
  end function

  type (mp_realm) function mpsub_dm (a, b)
      real*8, intent(in) :: a
      type (mp_realm), intent(in) :: b
      mpsub_dm%mpr(1) = mpm5
      call f_mpsub_dq (a, b%mpr, mpsub_dm%mpr)
  end function

! Unary Minus
  type (mp_realm) function mpneg_m (a)
    type (mp_realm), intent(in) :: a
    mpneg_m%mpr(1) = mpm5
    call f_mpneg_q (a%mpr, mpneg_m%mpr);
  end function

! Multiplications
  type (mp_realm) function mpmul_mm (a, b)
      type (mp_realm), intent(in) :: a, b
      mpmul_mm%mpr(1) = mpm5
      call f_mpmul (a%mpr, b%mpr, mpmul_mm%mpr)
  end function

  type (mp_real) function mpmul_mq (a, b)
      type (mp_realm), intent(in) :: a
      type (mp_real), intent(in) :: b
      mpmul_mq%mpr(1) = mpwork5
      call f_mpmul (a%mpr, b%mpr, mpmul_mq%mpr)
  end function

  type (mp_real) function mpmul_qm (a, b)
      type (mp_real), intent(in) :: a
      type (mp_realm), intent(in) :: b
      mpmul_qm%mpr(1) = mpwork5
      call f_mpmul (a%mpr, b%mpr, mpmul_qm%mpr)
  end function

  type (mp_realm) function mpmul_md (a, b)
      type (mp_realm), intent(in) :: a
      real*8, intent(in) :: b
      mpmul_md%mpr(1) = mpm5
      call f_mpmul_qd (a%mpr, b, mpmul_md%mpr)
  end function

  type (mp_realm) function mpmul_dm (a, b)
      real*8, intent(in) :: a
      type (mp_realm), intent(in) :: b
      mpmul_dm%mpr(1) = mpm5
      call f_mpmul_qd (b%mpr, a, mpmul_dm%mpr)
  end function

! Divisions
  type (mp_realm) function mpdiv_mm (a, b)
      type (mp_realm), intent(in) :: a, b
      mpdiv_mm%mpr(1) = mpm5
      call f_mpdiv (a%mpr, b%mpr, mpdiv_mm%mpr)
  end function

  type (mp_real) function mpdiv_mq (a, b)
      type (mp_realm), intent(in) :: a
      type (mp_real), intent(in) :: b
      mpdiv_mq%mpr(1) = mpwork5
      call f_mpdiv (a%mpr, b%mpr, mpdiv_mq%mpr)
  end function

  type (mp_real) function mpdiv_qm (a, b)
      type (mp_real), intent(in) :: a
      type (mp_realm), intent(in) :: b
      mpdiv_qm%mpr(1) = mpwork5
      call f_mpdiv (a%mpr, b%mpr, mpdiv_qm%mpr)
  end function

  type (mp_realm) function mpdiv_md (a, b)
      type (mp_realm), intent(in) :: a
      real*8, intent(in) :: b
      mpdiv_md%mpr(1) = mpm5
      call f_mpdiv_qd (a%mpr, b, mpdiv_md%mpr)
  end function

  type (mp_realm) function mpdiv_dm (a, b)
      real*8, intent(in) :: a
      type (mp_realm), intent(in) :: b
      mpdiv_dm%mpr(1) = mpm5
      call f_mpdiv_dq (a, b%mpr, mpdiv_dm%mpr)
  end function

! Powers
  type (mp_realm) function mpexp_mi (a, ib)
      type (mp_realm), intent(in) :: a
      integer, intent(in) :: ib 
      mpexp_mi%mpr(1) = mpm5
      call f_mppwr_qi (a%mpr, ib, mpexp_mi%mpr)
  end function

! Assignments
  subroutine mpeq_mm (a, b)
      type (mp_realm), intent(inout) :: a
      type (mp_realm), intent(in) :: b
      a%mpr(1) = mpm5
      call f_mpeq (b%mpr, a%mpr)
  end subroutine

  subroutine mpeq_mq (a, b)
      type (mp_realm), intent(inout) :: a
      type (mp_real), intent(in) :: b
      a%mpr(1) = mpm5
      call f_mpeq (b%mpr, a%mpr)
  end subroutine

  subroutine mpeq_qm (a, b)
      type (mp_real), intent(inout) :: a
      type (mp_realm), intent(in) :: b
      a%mpr(1) = mpwork5
      call f_mpeq (b%mpr, a%mpr)
  end subroutine

  subroutine mpeq_md (a, b)
      type (mp_realm), intent(inout) :: a
      real*8, intent(in) :: b
      a%mpr(1) = mpm5
      call f_mpeq_d (b, a%mpr)
  end subroutine

  subroutine mpeq_dm (a, b)
      real*8, intent(out) :: a
      type (mp_realm), intent(in) :: b
      double precision db
      integer ib
      call f_mpmdc (b%mpr, db, ib)
      a = db * 2.d0 ** ib
  end subroutine

! Equality
  logical function mpeqt_mm (a, b)
      type (mp_realm), intent(in) :: a, b
      integer ic
      call f_mpcpr (a%mpr, b%mpr, ic)
      if (ic .eq. 1) then
         mpeqt_mm = .true.
      else
         mpeqt_mm = .false.
      endif
      return
  end function

! Inequality
  logical function mpnet_mm (a, b)
      type (mp_realm), intent(in) :: a, b
      integer ic
      call f_mpcpr (a%mpr, b%mpr, ic)
      if (ic .ne. 1) then
         mpnet_mm = .true.
      else
         mpnet_mm = .false.
      endif
      return
  end function
      
! Less-Than-Or-Equal-To
  logical function mplet_mm (a, b)
      type (mp_realm), intent(in) :: a, b
      integer ic
      call f_mplet (a%mpr, b%mpr, ic)
      if (ic .eq. 1) then
         mplet_mm = .true.
      else
         mplet_mm = .false.
      endif
      return
  end function
      
! Greater-Than-Or-Equal-To
  logical function mpget_mm (a, b)
      type (mp_realm), intent(in) :: a, b
      integer ic
      call f_mpget (a%mpr, b%mpr, ic)
      if (ic .eq. 1) then
         mpget_mm = .true.
      else
         mpget_mm = .false.
      endif
      return
  end function

! Less-Than
  logical function mpltt_mm (a, b)
      type (mp_realm), intent(in) :: a, b
      integer ic
      call f_mpltt (a%mpr, b%mpr, ic)
      if (ic .eq. 1) then
         mpltt_mm = .true.
      else
         mpltt_mm = .false.
      endif
      return
  end function
      
! Greater-Than
  logical function mpgtt_mm (a, b)
      type (mp_realm), intent(in) :: a, b
      integer ic
      call f_mpgtt (a%mpr, b%mpr, ic)
      if (ic .eq. 1) then
         mpgtt_mm = .true.
      else
         mpgtt_mm = .false.
      endif
      return
  end function
      
! Absolute value
  type (mp_realm) function mp_absm (a)
    type (mp_realm), intent(in) :: a
    mp_absm%mpr(1) = mpm5
    call f_mpabs (a%mpr, mp_absm%mpr)
  end function

  type (mp_realm) function mp_aintm (a)
      type (mp_realm), intent(in) :: a
      mp_aintm%mpr(1) = mpm5
      call f_mpaint (a%mpr, mp_aintm%mpr)
  end function

  type (mp_realm) function mp_anintm (a)
      type (mp_realm), intent(in) :: a
      mp_anintm%mpr(1) = mpm5
      call f_mpnint (a%mpr, mp_anintm%mpr)
  end function

  double precision function mp_mtod (a)
      type (mp_realm), intent (in):: a
      call f_mpdble (a%mpr, mp_mtod)
  end function

  type (mp_realm) function mp_maxm (a, b)
      type (mp_realm), intent(in) :: a, b
      integer ic
      mp_maxm%mpr(1) = mpm5
      call f_mpget (a%mpr, b%mpr, ic)
      if (ic .eq. 1) then
         call f_mpeq (a%mpr, mp_maxm%mpr)
      else
         call f_mpeq (b%mpr, mp_maxm%mpr)
      endif
  end function

  type (mp_realm) function mp_maxm3 (a, b, c)
      type (mp_realm), intent(in) :: a, b, c
      integer ic
      mp_maxm3%mpr(1) = mpm5
      call f_mpget (a%mpr, b%mpr, ic)
      if (ic .eq. 1) then
         call f_mpeq (a%mpr, mp_maxm3%mpr)
      else
         call f_mpeq (b%mpr, mp_maxm3%mpr)
      endif
      call f_mpget (c%mpr, mp_maxm3%mpr, ic)
      if (ic .eq. 1) call f_mpeq (c%mpr, mp_maxm3%mpr)
  end function

  type (mp_realm) function mp_minm (a, b)
      type (mp_realm), intent(in) :: a, b
      integer ic
      mp_minm%mpr(1) = mpm5
      call f_mplet (a%mpr, b%mpr, ic)
      if (ic .eq. 1) then
         call f_mpeq (a%mpr, mp_minm%mpr)
      else
         call f_mpeq (b%mpr, mp_minm%mpr)
      endif
  end function

  type (mp_realm) function mp_minm3 (a, b, c)
      type (mp_realm), intent(in) :: a, b, c
      integer ic
      mp_minm3%mpr(1) = mpm5
      call f_mplet (a%mpr, b%mpr, ic)
      if (ic .eq. 1) then
         call f_mpeq (a%mpr, mp_minm3%mpr)
      else
         call f_mpeq (b%mpr, mp_minm3%mpr)
      endif
      call f_mplet (c%mpr, mp_minm3%mpr, ic)
      if (ic .eq. 1) call f_mpeq (c%mpr, mp_minm3%mpr)
  end function

  type (mp_realm) function mp_signm (a, b)
      type (mp_realm), intent (in) :: a, b
      intrinsic dsign
      mp_signm%mpr(1) = mpm5
      call f_mpeq (a%mpr, mp_signm%mpr)
      mp_signm%mpr(2) = dsign (mp_signm%mpr(2), b%mpr(2))
  end function

! SQRT, etc.
  type (mp_realm) function mp_sqrtm (a)
    type (mp_realm), intent(in) :: a
    mp_sqrtm%mpr(1) = mpm5
    call f_mpsqrt (a%mpr, mp_sqrtm%mpr)
  end function

!  Conversion
  type (mp_real) function mp_mtoq (a)
      type (mp_realm), intent(in) :: a
      mp_mtoq%mpr(1) = mpwork5
      call f_mpeq (a%mpr, mp_mtoq%mpr)
  end function

  type (mp_realm) function mp_qtom (a)
      type (mp_real), intent(in) :: a
      mp_qtom%mpr(1) = mpm5
      call f_mpeq (a%mpr, mp_qtom%mpr)
  end function

  type (mp_realm) function mp_dtom (da)
      real*8, intent(in) :: da
      mp_dtom%mpr(1) = mpm5
      call f_mpdmc (da, mp_dtom%mpr)
  end function

! Input

  subroutine mp_inpm (iu, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: iu
    type (mp_realm), intent (out) :: q1, q2, q3, q4, q5, q6, q7, q8, q9
    optional :: q2, q3, q4, q5, q6, q7, q8, q9

    call mpinpm (iu, q1)
    if (present (q2)) call mpinpm (iu, q2)
    if (present (q3)) call mpinpm (iu, q3)
    if (present (q4)) call mpinpm (iu, q4)
    if (present (q5)) call mpinpm (iu, q5)
    if (present (q6)) call mpinpm (iu, q6)
    if (present (q7)) call mpinpm (iu, q7)
    if (present (q8)) call mpinpm (iu, q8)
    if (present (q9)) call mpinpm (iu, q9)
    return
  end subroutine

  subroutine mpinpm (iu, a)
    integer, intent(in) :: iu
    type (mp_realm), intent(out) :: a
    integer i, l, l1
    character*200 line

    l = 0
    a%mpr(1) = mpm5

100 continue
    read(iu, '(200a1)', end = 200) line

    do i = 200, 1, -1
      if (line(i:i) /= ' ') goto 110
    enddo

    i = 0

110 continue
    l1 = i

    do i = 1, l1
      az(l+i) = line(i:i)
    enddo

    if (line(l1:l1) /= ',') then
      l = l + l1
      goto 100
    endif
    call mpinpc (az, l, a%mpr)
    goto 300

200 continue
    write (6, 1)
1   format ('mpinpm: no comma terminating multiprecion input.')

300 continue
    return
  end subroutine

! Output
  subroutine mp_outm (iu, q1, q2, q3, q4, q5, q6, q7, q8, q9)
      integer, intent(in) :: iu
      type (mp_realm), intent(in) :: q1, q2, q3, q4, q5, q6, q7, q8, q9
      optional :: q2, q3, q4, q5, q6, q7, q8, q9
      call mpoutm (iu, q1)
      if (present (q2)) call mpoutm (iu, q2)
      if (present (q3)) call mpoutm (iu, q3)
      if (present (q4)) call mpoutm (iu, q4)
      if (present (q5)) call mpoutm (iu, q5)
      if (present (q6)) call mpoutm (iu, q6)
      if (present (q7)) call mpoutm (iu, q7)
      if (present (q8)) call mpoutm (iu, q8)
      if (present (q9)) call mpoutm (iu, q9)
      return
  end subroutine

  subroutine mpoutm (iu, q)
      integer, intent(in) :: iu
      type (mp_realm), intent(in) :: q
      integer i, l
      call mpoutc (q%mpr, az, l)
      az(l+1) = ','
      write(iu, '(78A1)') (az(i), i = 1, l+1)
  end subroutine

end module mpmodulem
