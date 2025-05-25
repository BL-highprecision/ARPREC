!  mp_mod.f

!  ARPREC Fortran-90 translation modules.
!  Version date 2010-08-25
!  Copyright 2002-2010
  
!  Authors: Sherry Li (LBNL) and David H Bailey (LBNL)
!  Contact email:  dhbailey@lbl.gov

!  This work was supported by the Director, Office of Science, Division
!  of Mathematical, Information, and Computational Sciences of the
!  U.S. Department of Energy under contract number DE-AC03-76SF00098.

module mpdefmod

!   mpipl is the maximum precision level, in digits.

integer mpipl
parameter (mpipl = 2000)

!   NOTE:  This code should not be changed below this point in ordinary usage.

integer kdb, mpldb, mpnbt, mpoud, mpwds, mpwds51, mpwork5, new_mpipl, &
  new_mpwork
real*8 digits_per_word
parameter (mpldb = 6, mpnbt = 48)
parameter (digits_per_word = 14.44943979187109d0)     ! mpnbt*log10(2.d0)
parameter (mpwds = (mpipl-1) / digits_per_word + 2.d0)
parameter (kdb = kind (0.d0), mpwds51 = mpwds + 5 + 1)

! Arbitrary precision datatypes

type mp_integer
  sequence
  real*8 mpi(mpwds+5)
end type
type mp_real
  sequence
  real*8 mpr(mpwds+5)
end type
type mp_complex
  sequence
  real*8 mpc(2*(mpwds+5))
end type

type (mp_real), public:: mpl02, mpl10, mppic, mpeps, mplrg, mpsml

contains
      
  subroutine mpinit (n_mpipl, iconst, filename)
!  MPINIT must be called at the start of execution in the user's main program.
!  It sets the numeric precision level, the MP epsilon level, the output
!  precision level, and computes the constants Pi, Log(2) and Log(10).

!  The arguments are as follows:
!  n_mpipl: integer, optional.  Working precision level, in digits.
!     Must not exceed mpipl, which is set a few lines above.
!     If n_mpipl is not present, then precision level is set to mpipl.
!  iconst: integer, optional.  Flag for computing constants log(2), log(10), pi.
!     0: Do not compute constants, and do not read or write to the file.
!     1: Compute the constants if filename is not present; if filename is
!          present, then read constants from file specified by filename.
!     2: Compute the constants; if filename is present, then write the
!          constants to file specified by filename.
!     If iconst is not present, then compute constants, as if iconst = 1.
!  filename: character, optional.  Name of file containing log(2), log(10), pi.
!     utilized only if iconst = 1 or iconst = 2 (see above).
!     If filename is not present, then no attempt is made to read or write file.

      implicit none
      integer, intent(in), optional :: n_mpipl, iconst
      character*(*), intent(in), optional :: filename
      integer iconst_temp

      if (present (n_mpipl)) then
         if (n_mpipl .gt. mpipl) then
            write(mpldb, *) &
              '*** MPINIT: new precision level is too high'
            stop
         endif
         if (present (iconst)) then
            if (iconst == 0) then
               call f_mpinit (n_mpipl, iconst, mpwds, mpwork5, &
                  mpeps%mpr, mpl02%mpr, mpl10%mpr, mppic%mpr)
               new_mpipl = n_mpipl
               new_mpwork = mpwork5 - 5
            elseif (iconst == 1 .and. present (filename)) then
               if (filename == ' ') then
                  write (mpldb, *) '*** MPINIT: filename is blank'
                  stop
               endif
               iconst_temp = 0
               call f_mpinit (n_mpipl, iconst_temp, mpwds, mpwork5, &
                 mpeps%mpr, mpl02%mpr, mpl10%mpr, mppic%mpr)
               new_mpipl = n_mpipl
               new_mpwork = mpwork5 - 5
               open (51, file = filename, form = 'unformatted')
               rewind (51)
               read (51) mpl02
               read (51) mpl10
               read (51) mppic
               close (51)
            elseif (iconst == 2 .and. present (filename)) then
               if (filename == ' ') then
                  write (mpldb, *) '*** MPINIT: filename is blank'
                  stop
               endif
               call f_mpinit (n_mpipl, iconst, mpwds, mpwork5, &
                 mpeps%mpr, mpl02%mpr, mpl10%mpr, mppic%mpr)
               new_mpipl = n_mpipl
               new_mpwork = mpwork5 - 5
               open (51, file = filename, form = 'unformatted')
               rewind (51)
               write (51) mpl02
               write (51) mpl10
               write (51) mppic
               close (51)
            else
               call f_mpinit (n_mpipl, iconst, mpwds, mpwork5, &
                 mpeps%mpr, mpl02%mpr, mpl10%mpr, mppic%mpr)
               new_mpipl = n_mpipl
               new_mpwork = mpwork5 - 5
            endif
         else
           iconst_temp = 1
           call f_mpinit (n_mpipl, iconst_temp, mpwds, mpwork5, &
             mpeps%mpr, mpl02%mpr, mpl10%mpr, mppic%mpr)
           new_mpipl = n_mpipl
           new_mpwork = mpwork5 - 5
         endif
      else
         iconst_temp = 1
         call f_mpinit (mpipl, iconst_temp, mpwds, mpwork5, &
           mpeps%mpr, mpl02%mpr, mpl10%mpr, mppic%mpr)
         new_mpipl = mpipl
         new_mpwork = mpwork5 - 5
      endif

      mpoud = 56
      mplrg%mpr(1) = mpwork5
      mplrg%mpr(2) = 1.
      mplrg%mpr(3) = 2**27
      mplrg%mpr(4) = 1.
      mplrg%mpr(5) = 0.
      mpsml%mpr(1) = mpwork5
      mpsml%mpr(2) = 1.
      mpsml%mpr(3) = -2**27
      mpsml%mpr(4) = 1.
      mpsml%mpr(5) = 0.
  end subroutine

  subroutine mpsetprec (num_digits)
    integer num_digits

    if (num_digits > new_mpipl) then
      write (mpldb, *) 'mpsetprec: invalid argument; precision set to ', &
        new_mpipl, ' digits'
      call f_mpsetprec (new_mpipl)
    else
      call f_mpsetprec (num_digits)
    endif
  end subroutine

  subroutine mpgetprec (num_digits)
    integer num_digits
    call f_mpgetprec (num_digits)
  end subroutine

  subroutine mpsetprecwords (num_words)
    integer num_words
    if (num_words > new_mpwork) then
      write (mpldb, *) &
        'mpsetprecwords: invalid argument; precision set to ', &
        new_mpwork, ' words'
      call f_mpsetprecwords (new_mpwork)
    else
      call f_mpsetprecwords (num_words)
    endif
  end subroutine

  subroutine mpgetprecwords (num_words)
    integer num_words
    call f_mpgetprecwords (num_words)
  end subroutine

  subroutine mpsetoutputprec (num_digits)
    integer num_digits
    if (num_digits > new_mpipl) then
      write (mpldb, *) &
        'mpsetoutputprec: invalid argument; output precision set to ', &
          new_mpipl, ' digits'
      mpoud = new_mpipl
    else
      mpoud = num_digits
    endif
  end subroutine

  subroutine mpgetoutputprec (num_digits)
    integer num_digits
    num_digits = mpoud
  end subroutine

  subroutine mpgetpar (s, n, k)
!  MPGETPAR retrieves some ARPREC C++ integer parameters.
      character*(*), intent(in) :: s
      integer, intent(out) :: n
      integer, intent(in), optional :: k
      
      if (s == 'mpnw') then 
        call f_mpgetpar (1, n, 0)
      elseif (s == 'mpidb') then
        call f_mpgetpar (2, n, 0)
      elseif (s == 'mpndb') then
        call f_mpgetpar (3, n, 0)
      elseif (s == 'mpmcr') then
        call f_mpgetpar (4, n, 0)
      elseif (s == 'mpird') then
        call f_mpgetpar (5, n, 0)
      elseif (s == 'mpier') then
        call f_mpgetpar (6, n, 0)
      elseif (s == 'mpker') then
        call f_mpgetpar (7, n, k)
      else
        write (mpldb, 1) s
1       format ('mpgetpar: invalid parameter name: ',a)
        n = 0
      endif
  end subroutine

  subroutine mpsetpar (s, n, k)
!  MPSETPAR sets some ARPREC C++ integer parameters.
      character*(*), intent(in) :: s
      integer, intent(in) :: n
      integer, intent(in), optional :: k
      
      if (s == 'mpnw') then 
         call f_mpsetpar(1, n, 0)
      elseif (s == 'mpidb') then
         call f_mpsetpar(2, n, 0)
      elseif (s == 'mpndb') then
         call f_mpsetpar(3, n, 0)
      elseif (s == 'mpmcr') then
         call f_mpsetpar(4, n, 0)
      elseif (s == 'mpird') then
         call f_mpsetpar(5, n, 0)
      elseif (s == 'mpier') then
         call f_mpsetpar(6, n, 0)
      elseif (s == 'mpker') then
        call f_mpsetpar (7, n, k)
      else
        write (mpldb, 1) s
1       format ('mpsetpar: invalid parameter name: ',a)
      endif
  end subroutine

  subroutine mpinpc (a, n, b)

!   Converts the CHARACTER*1 array A of length N into the MP number B.  The
!   string A must be in the format '10^s a x tb.c' where a, b and c are digit
!   strings; s and t are '-', '+' or blank; x is either 'x' or '*'.  Blanks may
!   be embedded anywhere.  The exponent string a is limited to nine digits and
!   80 total characters, including blanks.  The exponent portion (i.e. the
!   portion up to and including x) and the period may optionally be omitted.
!   Debug output starts with MPIDB = 7.

!   Max SP space for B: MPNW + 4 cells.

!   The following example shows how this routine may be used to input a MP
!   number:

!   CHARACTER*1 CX(800)
!   READ (1, '(80A1)') (CX(I), I = 1, ND)
!   CALL MPINPC (CX, ND, B)

  integer i, ib, id, ier, ip, it, i1, i2, k0, k1, k2, l1, n, nb, nn, n5
  double precision bi
  character*1 a(n), ai
  character*10 dig
  character*80 ca
  parameter (dig = '0123456789')
  real*8 b(0:new_mpwork+5), f(0:8), s(0:3*(new_mpwork+5))
  integer mpidb, mpier, mpndb

  call mpgetpar ('mpier', mpier)
  if (mpier .ne. 0) goto 220
  call mpgetpar ('mpidb', mpidb)
  if (mpidb .ge. 7) then
    call mpgetpar ('mpndb', mpndb)
    no = min (n, int (7.225 * mpndb) + 20)
    write (mpldb, 1) (a(i), i = 1, no)
1   format ('MPINPC I'/(78a1))
  endif

  n5 = new_mpwork + 5
  k0 = 0
  k1 = k0 + n5
  k2 = k1 + n5
  s(k0) = n5
  s(k1) = n5
  s(k2) = n5
  b(0) = n5
  b(1) = 0.d0
  b(2) = 0.d0
  i1 = 1
  nn = 0

!   Find the carat, period, plus or minus sign, whichever comes first.

  do i = 1, n
    ai = a(i)
    if (ai .eq. '^') goto 110
    if (ai .eq. '.' .or. ai .eq. '+' .or. ai .eq. '-') goto 160
  enddo

  goto 160

!   Make sure number preceding the carat is 10.

110 continue

  i2 = i - 1
  if (i2 .gt. 80) then
    ier = 1
    goto 210
  endif
  ca = ' '

  do i = 1, i2
    ai = a(i)
    if (ai .eq. ' ') then
      goto 120
    elseif (index (dig, ai) .eq. 0) then
      ier = 2
      goto 210
    endif
    ca(i:i) = ai
120 continue
  enddo

  nn = mpdigin (ca, 80)
  if (nn .ne. 10) then
    ier = 3
    goto 210
  endif
  i1 = i2 + 2

!   Find the x or *.

  do i = i1, n
    ai = a(i)
    if (ai .eq. 'x' .or. ai .eq. '*') goto 140
  enddo

  ier = 4
  goto 210

!   Convert the exponent.

140 continue

i2 = i - 1
  l1 = i2 - i1 + 1
  if (l1 .gt. 80) then
    ier = 5
    goto 210
  endif
  ca = ' '
  id = 0
  is = 1

  do i = 1, l1
    ai = a(i+i1-1)
    if (ai .eq. ' ' .or. ai .eq. '+') then
      goto 150
    elseif (ai .eq. '-' .and. id .eq. 0) then
      id = 1
      is = -1
      ca(i:i) = ' '
    else
      if (index (dig, ai) .eq. 0) then
        ier = 6
        goto 210
      endif
      id = 1
      ca(i:i) = ai
    endif
150 continue
  enddo

  nn = is * mpdigin (ca, 80)
  i1 = i2 + 2

!   Find the next nonblank character.

160 continue

  do i = i1, n
    if (a(i) .ne. ' ') goto 180
  enddo

  ier = 7
  goto 210

!   Check if the nonblank character is a plus or minus sign.

180 continue

  i1 = i
  if (a(i1) .eq. '+') then
    i1 = i1 + 1
    is = 1
  elseif (a(i1) .eq. '-') then
    i1 = i1 + 1
    is = -1
  else
    is = 1
  endif
  nb = 0
  ib = 0
  id = 0
  ip = 0
  s(k2+1) = 0.
  s(k2+2) = 0.
  f(0) = 8.d0
  f(1) = 1.d0
  f(2) = 0.d0
  it = 0

190 continue

  ip = 0
  ca(1:12) = '000000000000'

!   Scan for digits, looking for the period also.  On the first pass we just
!   count, so that on the second pass it will come out right.

  do i = i1, n
    ai = a(i)
    if (ai .eq. ' ') then
    elseif (ai .eq. '.') then
      if (ip .ne. 0) then
        ier = 8
        goto 210
      endif
      ip = id
    elseif (index (dig, ai) .eq. 0) then
      ier = 9
      goto 210
    else
      ib = ib + 1
      id = id + 1
      ca(ib:ib) = ai
    endif
    if (ib .eq. 12 .or. i .eq. n .and. ib .ne. 0) then
      if (it .ne. 0) then
        nb = nb + 1
        bi = mpdigin (ca(1:12), 12)
        call f_mpmul_qd (s(k2), 1.d12, s(k0))
        if (bi .ne. 0) then
          f(1) = 1.d0
          f(3) = bi
        else
          f(1) = 0.d0
        endif
        call f_mpadd (s(k0), f, s(k2))
        ca(1:12) = '000000000000'
      endif
      if (i .ne. n) ib = 0
    endif
  enddo

  if (it .eq. 0) then
    ib = 12 - ib
    if (ib .eq. 12) ib = 0
    it = 1
    goto 190
  endif

  if (is .eq. -1) s(k2+1) = - s(k2+1)
  if (ip .eq. 0) ip = id
  nn = nn + ip - id
  f(1) = 1.d0
  f(3) = 10.d0
  
  call f_mppwr_qi (f, nn, s(k0))
  call f_mpmul (s(k2), s(k0), s(k1))
  call f_mpeq (s(k1), b)

  if (mpidb .ge. 7) then
    no = min (int (abs (b(1))), mpndb) + 2
    write (mpldb, 2) (b(i), i = 1, no)
2   format ('MPINPC O'/(6f12.0))
  endif
  goto 220

210 continue

  call mpgetpar ('mpker', i, 41)
  if (i > 0) then
    write (mpldb, 3) ier
3   format ('*** MPINPC: Syntax error in literal string; ier =',i4)
    call mpsetpar ('mpier', 41)
    if (i == 2) stop
  endif

220 continue
  return
  end subroutine

  subroutine mpoutc (a, b, n)

!   Converts the MP number A into character form in the CHARACTER*1 array B.
!   N (an output parameter) is the length of the output.  In other words, B is
!   contained in B(1), ..., B(N).  The format is analogous to the Fortran
!   exponential format (E format), except that the exponent is placed first.
!   Debug output starts with MPIDB = 7.

!   Max CHARACTER*1 space for B: 7.225 * MPNW + 30 cells.

!   This routine is called by MPOUT, but it may be directly called by the user
!   if desired for custom output.  Example:

!   CHARACTER*1 CX(800)
!   CALL MPOUTC (A, CX, ND)
!   WRITE (1, '(20A1/(72A1))') (CX(I), I = 1, ND)

  double precision aa, al2, t1
  character*1 b(*)
  character*16 ca
  parameter (al2 = 0.301029995663981195d0, con = 0.8304820235d0)
  real*8 an, a(0:*), f(0:8), s(0:2*(new_mpwork+5))
  real*8 mpbdx, mprdx
  parameter (mpbdx = 2.d0 ** mpnbt, mprdx = 0.5d0 ** mpnbt)
  integer mpidb, mpier, mpndb

  call mpgetpar ('mpier', mpier)
  if (mpier .ne. 0) goto 200
  call mpgetpar ('mpidb', mpidb)

  if (mpidb .ge. 7) then
    call mpgetpar ('mpndb', mpndb)
    no = min (int (abs (a(1))), mpndb) + 2
    write (mpldb, 1) (a(i), i = 1, no)
1   format ('MPOUTC I'/(6f12.0))
  endif

  ia = sign (1.d0, a(1))
  na = min (int (abs (a(1))), new_mpwork)
  n5 = new_mpwork + 5
  k0 = 0
  k1 = k0 + n5
  s(0) = n5
  s(n5) = n5
  f(0) = 8.d0
  f(1) = 1.d0
  f(2) = 0.d0
  f(3) = 10.d0
  f(4) = 0.d0
  f(5) = 0.d0

!   Determine exact power of ten for exponent.

  if (na .ne. 0) then
    aa = a(3)
    if (na .ge. 2) aa = aa + mprdx * a(4)
    t1 = al2 * mpnbt * a(2) + log10 (aa)
    if (t1 .ge. 0.d0) then
      nx = t1
    else
      nx = t1 - 1.d0
    endif
    call f_mppwr_qi (f, nx, s(k0))
    call f_mpdiv (a, s(k0), s(k1))

!   If we didn't quite get it exactly right, multiply or divide by 10 to fix.

100 continue

    if (s(k1+2) .lt. 0.d0) then
      nx = nx - 1
      call f_mpmul_qd (s(k1), 10.d0, s(k0))
      call f_mpeq (s(k0), s(k1))
      goto 100
    elseif (s(k1+3) .ge. 10.d0) then
      nx = nx + 1
      call f_mpdiv_qd (s(k1), 10.d0, s(k0))
      call f_mpeq (s(k0), s(k1))
      goto 100
    endif
    s(k1+1) = abs (s(k1+1))
  else
    nx = 0
  endif

!   Place exponent first instead of at the very end as in Fortran.

  b(1) = '1'
  b(2) = '0'
  b(3) = ' '
  b(4) = '^'
  ca = mpdigout (dble (nx), 10)

  do i = 1, 10
    b(i+4) = ca(i:i)
  enddo

  b(15) = ' '
  b(16) = 'x'
  b(17) = ' '

!   Insert sign and first digit.

  if (ia .eq. -1) then
    b(18) = '-'
  else
    b(18) = ' '
  endif
  if (na .ne. 0) then
    an = s(k1+3)
  else
    an = 0
  endif

  ca = mpdigout (an, 1)

  b(19) = ca(1:1)
  b(20) = '.'
  ix = 20
  if (na .eq. 0) goto 190
  f(3) = an
  call f_mpsub (s(k1), f, s(k0))
  if (s(k0+1) .eq. 0.d0) goto 190
  call f_mpmul_qd (s(k0), 1.d12, s(k1))

  nl = max (new_mpwork * digits_per_word / 12.d0 - 1.d0, 1.d0)
  nl = min (nl, mpoud / 12 + 1)

!   Insert the digits of the remaining words.

  do j = 1, nl
    if (s(k1+2) .eq. 0.d0) then
      an = s(k1+3)
      f(1) = 1.d0
      f(3) = an
    else
      f(1) = 0.d0
      an = 0.d0
    endif
    ca = mpdigout (an, 12)

    do i = 1, 12
      if (ca(i:i) == ' ') ca(i:i) = '0'
      b(i+ix) = ca(i:i)
    enddo

    ix = ix + 12
    call f_mpsub (s(k1), f, s(k0))
    call f_mpmul_qd (s(k0), 1.d12, s(k1))
    if (s(k1+1) .eq. 0.d0) goto 140
  enddo

!   Check if trailing zeroes should be trimmed.

  j = nl + 1

140 continue

  l = ix
  if (b(l) .eq. '0' .and. b(l-1) .eq. '0' .or. (j .gt. nl .and. &
    b(l-2) .eq. '0' .and. b(l-3) .eq. '0')) then
    b(l) = ' '
    b(l-1) = ' '

    do i = l - 2, 21, -1
      if (b(i) .ne. '0') then
        ix = i
        goto 190
      endif
      b(i) = ' '
    enddo

    ix = 20

!   Check if trailing nines should be rounded up.

  elseif (j .gt. nl .and. b(l-2) .eq. '9' .and. b(l-3) .eq. '9') then
    b(l) = ' '
    b(l-1) = ' '

    do i = l - 2, 21, -1
      if (b(i) .ne. '9') goto 180
      b(i) = ' '
    enddo

!   We have rounded away all digits to the right of the decimal point, and the
!   digit to the left of the digit is a 9.  Set the digit to 1 and increase
!   the exponent by one.

    ix = 20
    if (b(19) .eq. '9') then
      b(19) = '1'
      ca = mpdigout (dble (nx + 1), 10)

      do i = 1, 10
        b(i+4) = ca(i:i)
      enddo
    else
      ca = b(19)
      an = mpdigin (ca, 1)
      ca = mpdigout (an + 1.d0, 1)
      b(19) = ca(1:1)
    endif
    goto 190

180 continue

    ca = b(i)
    an = mpdigin (ca, 1)
    ca = mpdigout (an + 1.d0, 1)
    b(i) = ca(1:1)
    ix = i
  endif

190 continue

  n = min (ix, mpoud + 20)

  if (mpidb .ge. 7) then
    no = min (n, 12 * mpndb + 20)
    write (mpldb, 2) (b(i), i = 1, no)
2   format ('MPOUTC O'/(78a1))
  endif

200 continue
  return
  end subroutine

  real*8 function mpdigin (ca, n)
    implicit none
    real*8 d1
    character*(*), ca
    character*16 digits
    integer i, is, k, n
    parameter (digits = '0123456789-')

    is = 1
    d1 = 0.d0

    do i = 1, n
      k = index (digits, ca(i:i)) - 1
      if (k < 0) then
        write (mpldb, *) 'mpdigin: non-digit in character string'
      elseif (k <= 9) then
        d1 = 10.d0 * d1 + k
      elseif (k == 10) then
        is = -1
      endif
    enddo

    mpdigin = is * d1
  end function

  character*16 function mpdigout (a, n)
    implicit none
    real*8 a, d1, d2
    character*16 ca, digits
    parameter (digits = '0123456789')
    integer i, is, k, n

    ca = ' '
    is = sign (1.d0, a)
    d1 = abs (a)

    do i = n, 1, -1
      d2 = aint (d1 / 10.d0)
      k = 1.d0 + (d1 - 10.d0 * d2)
      d1 = d2
      ca(i:i) = digits(k:k)
      if (d1 == 0.d0) goto 100
    enddo

    i = 0

100 continue

    if (is < 0 .and. i > 1) then
      ca(i-1:i-1) = '-'
    elseif (i == 0 .or. is < 0 .and. i == 1) then
      ca = '****************'
    endif

    mpdigout = ca
    return
  end function

  subroutine mpeform (a, n1, n2, b)
    type (mp_real) a
    integer n1, n2
    character*1 b(n1)
    call mpeformx (a%mpr, n1, n2, b)
    return
  end subroutine
    
  subroutine mpeformx (a, n1, n2, b)

!   This routine converts the MP number A to E format, i.e. E N1.N2.
!   B is the output array (type CHARACTER*1) of size N1.

      integer i, j, k, lex, n1, n2
      real*8 a(mpwds+5)
      character*1 b(n1), c(new_mpipl+100)

      if (n1 > mpoud) then
        write (mpldb, '("*** mpeformx: mpoud must exceed n1")')
        goto 110
      endif
      call mpoutc (a, c, n)

!   Find length of exponent field.

      do i = 5, 14
        if (c(i) /= ' ') goto 100
      enddo

100   continue

      lex = 15 - i
      k = n1 - lex - n2 - 4

!   Check for overflow of field length.

      if (k < 0) then
         do j = 1, n1
            b(j) = '*'
         enddo

         goto 110
      endif

!   Copy characters to appropriate positions.

      do j = 1, k
        b(j) = ' '
      enddo

      do j = 1, min (n2 + 3, n - 17)
        b(j+k) = c(j+17)
      enddo

      do j = n - 16, n2 + 3
        b(j+k) = '0'
      enddo

      b(k+n2+4) = 'e'

      do j = 1, lex
        b(j+k+n2+4) = c(i+j-1)
      enddo

110   continue

  return
  end subroutine

  subroutine mpfform (a, n1, n2, b)
    type (mp_real) a
    integer n1, n2
    character*1 b(n1)
    call mpfformx (a%mpr, n1, n2, b)
    return
  end subroutine
    
  subroutine mpfformx (a, n1, n2, b)

!   This routine converts the MP number A to F format, i.e. F N1.N2.
!   B is the output array (type CHARACTER*1) of size N1.
!
      real*8 a(mpwds+5), a1(mpwds+5), a2(mpwds+5)
      character*1 b(n1), c(new_mpipl+100)
      character*16 chr16

      if (n1 > mpoud) then
        write (mpldb, '("*** mpfformx: mpoud must exceed n1")')
        goto 200
      endif

!   Add a small "fuzz" in the direction of larger magnitude.

      a1(1) = mpwork5
      a2(1) = mpwork5
      a1(2) = 1.d0
      a1(3) = 0.d0
      a1(4) = 10.d0
      a1(5) = 0.d0
      a1(6) = 0.d0
      a1(7) = 0.d0
      i1 = - n2 - 1
      call f_mppwr_qi(a1, i1, a2)
      if (a(2) >= 0.d0) then
        call f_mpadd (a, a2, a1)
      else
        call f_mpsub (a, a2, a1)
      endif

!   Convert 

      call mpoutc (a1, c, n)
      chr16 = ' '

      do i = 1, 10
        chr16(i:i) = c(i+4)
      enddo

      ix = mpdigin (chr16, 16)

      if (a(2) .ge. 0.) then
         ls = 0
      else
         ls = 1
      endif
      if (ix .ge. 0 .and. a(2) .ne. 0.) then
         lz = 0
      else
         lz = 1
      endif
      mx = max (ix, 0)

!   Check for overflow of field length.

      if (ls + lz + mx + n2 + 2 .gt. n1) then
         do i = 1, n1
            b(i) = '*'
         enddo

         goto 200
      endif

!   Check if a zero should be output.

      if (a(2) .eq. 0 .or. -ix .gt. n2) then
         do i = 1, n1 - n2 - 2
            b(i) = ' '
         enddo

         b(n1-n2-1) = '0'
         b(n1-n2) = '.'

         do i = 1, n2
            b(i+n1-n2) = '0'
         enddo

         goto 200
      endif

!   Process other cases.
      
      do i = 1, n1 - n2 - mx - 2
         b(i) = ' '
      enddo

      if (a(2) .lt. 0.) b(n1-n2-mx-2) = '-'
      if (ix .ge. 0) then
         b(n1-n2-ix-1) = c(19)
         kx = min (n - 20, ix)

         do i = 1, kx
            b(i+n1-n2-ix-1) = c(i+20)
         enddo
         
         do i = kx + 1, ix
            b(i+n1-n2-ix-1) = '0'
         enddo

         b(n1-n2) = '.'
         kx = max (min (n - ix - 20, n2), 0)

         do i = 1, kx
            b(i+n1-n2) = c(i+ix+20)
         enddo

         do i = kx + 1, n2
            b(i+n1-n2) = '0'
         enddo
      else
         nx = - ix
         b(n1-n2-1) = '0'
         b(n1-n2) = '.'

         do i = 1, nx - 1
            b(i+n1-n2) = '0'
         enddo

         b(n1-n2+nx) = c(19)
         kx = min (n - 20, n2 - nx)

         do i = 1, kx
            b(i+n1-n2+nx) = c(i+20)
         enddo

         do i = kx + 1, n2 - nx
            b(i+n1-n2+nx) = '0'
         enddo
      endif

200   continue

      return
  end subroutine

  subroutine mpdexc (a, l, b)

!   This routine converts the character*1 string A, which
!   represents a multiprecision number in Fortran style, i.e.
!   '1234567890' or '1.23456789D-21', into standard MP binary format.
!   This routine is not intended to be called directly by the user.

    implicit none
    integer i, i1, l, l1, l2
    character*1 a(l), c(new_mpipl+100)
    real*8 b(*)

    do i = 1, l
      if (a(i) .eq. 'D' .or. a(i) .eq. 'E' .or. a(i) .eq. 'd' &
        .or. a(i) .eq. 'e') goto 100
    enddo

    call mpinpc (a, l, b)
    goto 110

100 i1 = i
    l1 = i - 1
    l2 = l - i
    c(1) = '1'
    c(2) = '0'
    c(3) = '^'

    do i = 1, l2
      c(i+3) = a(i+i1)
    enddo

    c(l2+4) = 'x'

    do i = 1, l1
      c(i+l2+4) = a(i)
    enddo

    call mpinpc (c, l1 + l2 + 4, b)
110 return
  end subroutine

  subroutine mpmdc(a, b, n)
!   This converts the MP number A to the DPE form (B, N), accurate to between
!   14 and 17 digits, depending on system.  B will be between 1 and MPBDX.
!   Debug output starts with MPIDB = 9.
      double precision a(mpwork5), b
      integer n
      call f_mpmdc(a, b, n)
  end subroutine

  subroutine mpdotd (n, isa, a, isb, db, c)
!   This routine computes the dot product of the MP vector A with the DP
!   vector DB, returning the MP result in C.  This routine is used in the
!   author's customized PSLQ routine, resulting in substantial speedup.
!   The length of both the A and DB vectors is N, and ISA and ISB are the 
!   skip distances between successive elements of A and DB, measured in 
!   MP words and DP words, respectively.  The DP values in DB must be
!   whole numbers, so for example they cannot be larger than 2^53.

      integer n, isa, isb
      double precision db(isb*n)
      type (mp_real) a(isa*n), c
      c%mpr(1) = mpwork5
      call f_mpdotd (n, isa * (mpwds + 5), a(1)%mpr, isb, db, c%mpr)
  end subroutine
      
end module mpdefmod


module mpintmod
!
!  This Fortran-90 module defines operator extensions involving the
!  MP_INTEGER datatype.  For operations involving two MP data types,
!  those whose first argument is MP_INTEGER are included here.
!  Others are handled in other modules.
!
!  The subroutines and functions defined in this module are private
!  and not intended to be called directly by the user.
!
  use mpdefmod
  implicit none
  private mpadd, mpsub, mpmul, mpdiv, mpeq, mpeq_str

! MP Integer operator extension interface blocks.

  interface operator (+)
      module procedure mpadd
      module procedure mpadd_ji
      module procedure mpadd_ij
      module procedure mpadd_jd
      module procedure mpadd_dj
  end interface

   interface operator (-)
      module procedure mpsub
      module procedure mpsub_ji
      module procedure mpsub_ij
      module procedure mpsub_jd
      module procedure mpsub_dj
! negation
      module procedure mpneg_j
   end interface

  interface operator (*)
      module procedure mpmul
      module procedure mpmul_ji
      module procedure mpmul_ij
      module procedure mpmul_jd
      module procedure mpmul_dj
  end interface

  interface operator (/)
      module procedure mpdiv
      module procedure mpdiv_ji
      module procedure mpdiv_ij
      module procedure mpdiv_jd
      module procedure mpdiv_dj
  end interface

  interface assignment (=)
      module procedure mpeq_str
      module procedure mpeq
      module procedure mpeq_ji
      module procedure mpeq_ij
  end interface

  interface operator (**)
      module procedure mpexp_jj
      module procedure mpexp_ji
!      module procedure mpexp_ij
  end interface

  interface operator (.eq.)
      module procedure mpeqt_jj
      module procedure mpeqt_jq
!      module procedure mpeqt_ji
  end interface

  interface operator (.ne.)
      module procedure mpnet_jj
      module procedure mpnet_jq
      module procedure mpnet_qj
  end interface

  interface operator (.le.)
      module procedure mplet_jj
      module procedure mplet_jq
      module procedure mplet_qj
  end interface

  interface operator (.ge.)
      module procedure mpget_jj
      module procedure mpget_jq
      module procedure mpget_qj
  end interface

  interface operator (.lt.)
      module procedure mpltt_jj
      module procedure mpltt_jq
      module procedure mpltt_qj
  end interface
  
  interface operator (.gt.)
      module procedure mpgtt_jj
      module procedure mpgtt_jq
      module procedure mpgtt_qj
  end interface

contains

! Additions
  type (mp_integer) function mpadd(ja, jb)
      type (mp_integer), intent(in) :: ja, jb
      mpadd%mpi(1) = mpwork5
      call f_mpadd(ja%mpi, jb%mpi, mpadd%mpi)
      call f_ovcheck (mpadd%mpi)
  end function

  type (mp_integer) function mpadd_ji(ja, ib)
      type (mp_integer), intent(in) :: ja
      integer, intent(in) :: ib
      mpadd_ji%mpi(1) = mpwork5
      call f_mpadd_ji(ja%mpi, ib, mpadd_ji%mpi)
      call f_ovcheck (mpadd_ji%mpi)
  end function

  type (mp_integer) function mpadd_ij(ia, jb)
      integer, intent(in) :: ia
      type (mp_integer), intent(in) :: jb
      mpadd_ij%mpi(1) = mpwork5
      call f_mpadd_ji(jb%mpi, ia, mpadd_ij%mpi)
      call f_ovcheck (mpadd_ij%mpi)
  end function

  type (mp_real) function mpadd_jd(ja, db)
      type (mp_integer), intent(in) :: ja
      real*8, intent(in) :: db
      mpadd_jd%mpr(1) = mpwork5
      call f_mpadd_jd(ja%mpi, db, mpadd_jd%mpr)
  end function

  type (mp_real) function mpadd_dj(da, jb)
      real*8, intent(in) :: da
      type (mp_integer), intent(in) :: jb
      mpadd_dj%mpr(1) = mpwork5
      call f_mpadd_jd(jb%mpi, da, mpadd_dj%mpr)
  end function

! Subtractions
  type (mp_integer) function mpsub(ja, jb)
      type (mp_integer), intent(in) :: ja, jb
      mpsub%mpi(1) = mpwork5
      call f_mpsub(ja%mpi, jb%mpi, mpsub%mpi)
      call f_ovcheck (mpsub%mpi)
  end function

  type (mp_integer) function mpsub_ji(ja, ib)
      type (mp_integer), intent(in) :: ja
      integer, intent(in) :: ib
      mpsub_ji%mpi(1) = mpwork5
      call f_mpsub_ji(ja%mpi, ib, mpsub_ji%mpi)
      call f_ovcheck (mpsub_ji%mpi)
  end function

  type (mp_integer) function mpsub_ij(ia, jb)
      integer, intent(in) :: ia
      type (mp_integer), intent(in) :: jb
      mpsub_ij%mpi(1) = mpwork5
      call f_mpsub_ij(ia, jb%mpi, mpsub_ij%mpi)
      call f_ovcheck (mpsub_ij%mpi)
  end function

  type (mp_real) function mpsub_jd(ja, db)
      type (mp_integer), intent(in) :: ja
      real*8, intent(in) :: db
      mpsub_jd%mpr(1) = mpwork5
      call f_mpsub_jd(ja%mpi, db, mpsub_jd%mpr)
  end function

  type (mp_real) function mpsub_dj(da, jb)
      real*8, intent(in) :: da
      type (mp_integer), intent(in) :: jb
      mpsub_dj%mpr(1) = mpwork5
      call f_mpsub_dj(da, jb%mpi, mpsub_dj%mpr)
  end function

! Unary Minus
  type (mp_integer) function mpneg_j(ja)
    type (mp_integer), intent(in) :: ja
    mpneg_j%mpi(1) = mpwork5
    call f_mpneg_q(ja%mpi, mpneg_j%mpi);
      call f_ovcheck (mpneg_j%mpi)
  end function

! Multiplications
  type (mp_integer) function mpmul(ja, jb)
      type (mp_integer), intent(in) :: ja, jb
      mpmul%mpi(1) = mpwork5
      call f_mpmul(ja%mpi, jb%mpi, mpmul%mpi)
      call f_ovcheck (mpmul%mpi)
  end function

  type (mp_integer) function mpmul_ji(ja, ib)
      type (mp_integer), intent(in) :: ja
      integer, intent(in) :: ib
      mpmul_ji%mpi(1) = mpwork5
      call f_mpmul_ji(ja%mpi, ib, mpmul_ji%mpi)
      call f_ovcheck (mpmul_ji%mpi)
  end function

  type (mp_integer) function mpmul_ij(ia, jb)
      integer, intent(in) :: ia
      type (mp_integer), intent(in) :: jb
      mpmul_ij%mpi(1) = mpwork5
      call f_mpmul_ji(jb%mpi, ia, mpmul_ij%mpi)
      call f_ovcheck (mpmul_ij%mpi)
  end function

  type (mp_real) function mpmul_jd(ja, db)
      type (mp_integer), intent(in) :: ja
      real*8, intent(in) :: db
      mpmul_jd%mpr(1) = mpwork5
      call f_mpmul_qd(ja%mpi, db, mpmul_jd%mpr)
  end function

  type (mp_real) function mpmul_dj(da, jb)
      real*8, intent(in) :: da
      type (mp_integer), intent(in) :: jb
      mpmul_dj%mpr(1) = mpwork5
      call f_mpmul_qd(jb%mpi, da, mpmul_dj%mpr)
  end function

! Divisions
  type (mp_integer) function mpdiv(a, b)
      type (mp_integer), intent(in) :: a, b
      mpdiv%mpi(1) = mpwork5
      call f_mpdiv_jj(a%mpi, b%mpi, mpdiv%mpi)
      call f_ovcheck (mpdiv%mpi)
  end function

  type (mp_integer) function mpdiv_ji(ja, ib)
      type (mp_integer), intent(in) :: ja
      integer, intent(in) :: ib
      mpdiv_ji%mpi(1) = mpwork5
      call f_mpdiv_ji(ja%mpi, ib, mpdiv_ji%mpi)
      call f_ovcheck (mpdiv_ji%mpi)
  end function

  type (mp_integer) function mpdiv_ij(ia, jb)
      integer, intent(in) :: ia
      type (mp_integer), intent(in) :: jb
      mpdiv_ij%mpi(1) = mpwork5
      call f_mpdiv_ij(ia, jb%mpi, mpdiv_ij%mpi)
      call f_ovcheck (mpdiv_ij%mpi)
  end function

  type (mp_real) function mpdiv_jd(ja, db)
      type (mp_integer), intent(in) :: ja
      real*8, intent(in) :: db
      mpdiv_jd%mpr(1) = mpwork5
      call f_mpdiv_qd(ja%mpi, db, mpdiv_jd%mpr)
  end function

  type (mp_real) function mpdiv_dj(da, jb)
      real*8, intent(in) :: da
      type (mp_integer), intent(in) :: jb
      mpdiv_dj%mpr(1) = mpwork5
      call f_mpdiv_dq(da, jb%mpi, mpdiv_dj%mpr)
  end function

! Assignments
  subroutine mpeq_str (ja, ab)
    type (mp_integer), intent(out):: ja
    character*(*), intent (in):: ab
    type (mp_integer):: jtmp1, jtmp2
    character*1 az(len(ab))
    integer i, l
    ja%mpi(1) = mpwork5
    l = len (ab)
    do i = 1, l
      az(i) = ab(i:i)
    enddo
    call mpdexc (az, l, jtmp1%mpi)
    call f_mpinfr (jtmp1%mpi, ja%mpi, jtmp2%mpi)
    return
  end subroutine

  subroutine mpeq(ja, jb)
      type (mp_integer), intent(inout) :: ja
      type (mp_integer), intent(in) :: jb
      ja%mpi(1) = mpwork5
      call f_mpeq(jb%mpi, ja%mpi)
      call f_ovcheck (ja%mpi)      
  end subroutine

  subroutine mpeq_ji(ja, ib)
      type (mp_integer), intent(inout) :: ja
      integer, intent(in) :: ib
      ja%mpi(1) = mpwork5
      call f_mpeq_ji(ib, ja%mpi)
      call f_ovcheck (ja%mpi)
  end subroutine

  subroutine mpeq_ij(ia, jb)
      integer, intent(out) :: ia
      type (mp_integer), intent(in) :: jb
      double precision db
      integer ib
      call f_mpmdc(jb%mpi, db, ib)
      ia = db * 2.d0 ** ib
  end subroutine

! Powers
  type (mp_integer) function mpexp_jj(ja, jb)
      type (mp_integer), intent(in) :: ja, jb
      mpexp_jj%mpi(1) = mpwork5
      call f_mppwr_jj(ja%mpi, jb%mpi, mpexp_jj%mpi)
      call f_ovcheck (mpexp_jj%mpi)
  end function

  type (mp_integer) function mpexp_ji(ja, ib)
      type (mp_integer), intent(in) :: ja
      integer, intent(in) :: ib 
      mpexp_ji%mpi(1) = mpwork5
      call f_mppwr_ji(ja%mpi, ib, mpexp_ji%mpi)
      call f_ovcheck (mpexp_ji%mpi)
  end function

! Equalities
  logical function mpeqt_jj(ja, jb)
      type (mp_integer), intent(in) :: ja, jb
      integer ic
      call f_mpcpr(ja%mpi, jb%mpi, ic)
      if (ic .eq. 1) then
         mpeqt_jj = .true.
      else
         mpeqt_jj = .false.
      endif
      return
  end function
      
  logical function mpeqt_jq(ja, b)
      type (mp_integer), intent(in) :: ja
      type (mp_real), intent(in) :: b
      integer ic
      call f_mpcpr(ja%mpi, b%mpr, ic)
      if (ic .eq. 1) then
         mpeqt_jq = .true.
      else
         mpeqt_jq = .false.
      endif
      return
  end function

! Inequalities  
  logical function mpnet_jj(ja, jb)
      type (mp_integer), intent(in) :: ja, jb
      integer ic
      call f_mpcpr(ja%mpi, jb%mpi, ic)
      if (ic .ne. 1) then
         mpnet_jj = .true.
      else
         mpnet_jj = .false.
      endif
      return
  end function
      
  logical function mpnet_jq(ja, b)
      type (mp_integer), intent(in) :: ja
      type (mp_real), intent(in) :: b
      integer ic
      call f_mpcpr(ja%mpi, b%mpr, ic)
      if (ic .ne. 1) then
         mpnet_jq = .true.
      else
         mpnet_jq = .false.
      endif
      return
  end function

  logical function mpnet_qj(a, jb)
      type (mp_real), intent(in) :: a
      type (mp_integer), intent(in) :: jb
      integer ic
      call f_mpcpr(a%mpr, jb%mpi, ic)
      if (ic .ne. 1) then
         mpnet_qj = .true.
      else
         mpnet_qj = .false.
      endif
      return
  end function

! Less-Than-Or-Equal-To
  logical function mplet_jj(ja, jb)
      type (mp_integer), intent(in) :: ja, jb
      integer ic
      call f_mplet(ja%mpi, jb%mpi, ic)
      if (ic .eq. 1) then
         mplet_jj = .true.
      else
         mplet_jj = .false.
      endif
      return
  end function
      
  logical function mplet_jq(ja, b)
      type (mp_integer), intent(in) :: ja
      type (mp_real), intent(in) :: b
      integer ic
      call f_mplet(ja%mpi, b%mpr, ic)
      if (ic .eq. 1) then
         mplet_jq = .true.
      else
         mplet_jq = .false.
      endif
      return
  end function

  logical function mplet_qj(a, jb)
      type (mp_real), intent(in) :: a
      type (mp_integer), intent(in) :: jb
      integer ic
      call f_mplet(a%mpr, jb%mpi, ic)
      if (ic .eq. 1) then
         mplet_qj = .true.
      else
         mplet_qj = .false.
      endif
      return
  end function

! Greater-Than-Or-Equal-To
  logical function mpget_jj(ja, jb)
      type (mp_integer), intent(in) :: ja, jb
      integer ic
      call f_mpget(ja%mpi, jb%mpi, ic)
      if (ic .eq. 1) then
         mpget_jj = .true.
      else
         mpget_jj = .false.
      endif
      return
  end function
      
  logical function mpget_jq(ja, b)
      type (mp_integer), intent(in) :: ja
      type (mp_real), intent(in) :: b
      integer ic
      call f_mpget(ja%mpi, b%mpr, ic)
      if (ic .eq. 1) then
         mpget_jq = .true.
      else
         mpget_jq = .false.
      endif
      return
  end function

  logical function mpget_qj(a, jb)
      type (mp_real), intent(in) :: a
      type (mp_integer), intent(in) :: jb
      integer ic
      call f_mpget(a%mpr, jb%mpi, ic)
      if (ic .eq. 1) then
         mpget_qj = .true.
      else
         mpget_qj = .false.
      endif
      return
  end function

! Less-Than
  logical function mpltt_jj(ja, jb)
      type (mp_integer), intent(in) :: ja, jb
      integer ic
      call f_mpltt(ja%mpi, jb%mpi, ic)
      if (ic .eq. 1) then
         mpltt_jj = .true.
      else
         mpltt_jj = .false.
      endif
      return
  end function
      
  logical function mpltt_jq(ja, b)
      type (mp_integer), intent(in) :: ja
      type (mp_real), intent(in) :: b
      integer ic
      call f_mpltt(ja%mpi, b%mpr, ic)
      if (ic .eq. 1) then
         mpltt_jq = .true.
      else
         mpltt_jq = .false.
      endif
      return
  end function

  logical function mpltt_qj(a, jb)
      type (mp_real), intent(in) :: a
      type (mp_integer), intent(in) :: jb
      integer ic
      call f_mpltt(a%mpr, jb%mpi, ic)
      if (ic .eq. 1) then
         mpltt_qj = .true.
      else
         mpltt_qj = .false.
      endif
      return
  end function

! Greater-Than
  logical function mpgtt_jj(ja, jb)
      type (mp_integer), intent(in) :: ja, jb
      integer ic
      call f_mpgtt(ja%mpi, jb%mpi, ic)
      if (ic .eq. 1) then
         mpgtt_jj = .true.
      else
         mpgtt_jj = .false.
      endif
      return
  end function
      
  logical function mpgtt_jq(ja, b)
      type (mp_integer), intent(in) :: ja
      type (mp_real), intent(in) :: b
      integer ic
      call f_mpgtt(ja%mpi, b%mpr, ic)
      if (ic .eq. 1) then
         mpgtt_jq = .true.
      else
         mpgtt_jq = .false.
      endif
      return
  end function

  logical function mpgtt_qj(a, jb)
      type (mp_real), intent(in) :: a
      type (mp_integer), intent(in) :: jb
      integer ic
      call f_mpgtt(a%mpr, jb%mpi, ic)
      if (ic .eq. 1) then
         mpgtt_qj = .true.
      else
         mpgtt_qj = .false.
      endif
      return
  end function

end module mpintmod


module mprealmod
!
!  This Fortran-90 module defines operator extensions involving the
!  MP_REAL datatype.  For operations involving two MP data types,
!  those whose first argument is MP_REAL are included here.
!  Others are handled in other modules.
!
!  The subroutines and functions defined in this module are private
!  and not intended to be called directly by the user.
!
  use mpdefmod
  implicit none
  integer, parameter :: out_str_len = 72
  private mpadd, mpsub, mpmul, mpdiv, mpeq, mpeq_qd, mpeq_str

! MP Real operator extension interface blocks.
  
  interface operator (+)
      module procedure mpadd
      module procedure mpadd_qi
      module procedure mpadd_iq
      module procedure mpadd_qd
      module procedure mpadd_dq
!      module procedure mpadd_qz
  end interface

  interface operator (-)
      module procedure mpsub
      module procedure mpsub_qi
      module procedure mpsub_iq
      module procedure mpsub_qd
      module procedure mpsub_dq
! negation
      module procedure mpneg_q
  end interface

  interface operator (*)
      module procedure mpmul
      module procedure mpmul_qi
      module procedure mpmul_iq
      module procedure mpmul_qd
      module procedure mpmul_dq
  end interface

  interface operator (/)
      module procedure mpdiv
      module procedure mpdiv_qi
      module procedure mpdiv_iq
      module procedure mpdiv_qd
      module procedure mpdiv_dq
  end interface

  interface assignment (=)
      module procedure mpeq_str
      module procedure mpeq
      module procedure mpeq_qi
      module procedure mpeq_iq
      module procedure mpeq_qd
      module procedure mpeq_dq
  end interface

  interface operator (**)
      module procedure mpexp
      module procedure mpexp_dq
      module procedure mpexp_qd
      module procedure mpexp_qi
  end interface

  interface operator (.eq.)
      module procedure mpeqt_qq
      module procedure mpeqt_qi
      module procedure mpeqt_iq
      module procedure mpeqt_qd
      module procedure mpeqt_dq
  end interface

  interface operator (.ne.)
      module procedure mpnet_qq
      module procedure mpnet_qi
      module procedure mpnet_iq
      module procedure mpnet_qd
      module procedure mpnet_dq
  end interface

  interface operator (.le.)
      module procedure mplet_qq
      module procedure mplet_qi
      module procedure mplet_iq
      module procedure mplet_qd
      module procedure mplet_dq
  end interface

  interface operator (.ge.)
      module procedure mpget_qq
      module procedure mpget_qi
      module procedure mpget_iq
      module procedure mpget_qd
      module procedure mpget_dq
  end interface

  interface operator (.lt.)
      module procedure mpltt_qq
      module procedure mpltt_qi
      module procedure mpltt_iq
      module procedure mpltt_qd
      module procedure mpltt_dq
  end interface

  interface operator (.gt.)
      module procedure mpgtt_qq
      module procedure mpgtt_qi
      module procedure mpgtt_iq
      module procedure mpgtt_qd
      module procedure mpgtt_dq
  end interface

contains

! Additions
  type (mp_real) function mpadd(a, b)
      type (mp_real), intent(in) :: a, b
      mpadd%mpr(1) = mpwork5
      call f_mpadd(a%mpr, b%mpr, mpadd%mpr)
  end function

  type (mp_real) function mpadd_qi(a, ib)
      type (mp_real), intent(in) :: a
      integer, intent(in) :: ib
      double precision db
      mpadd_qi%mpr(1) = mpwork5
      db = ib
      call f_mpadd_d(a%mpr, db, mpadd_qi%mpr)
  end function

  type (mp_real) function mpadd_iq(ia, b)
      integer, intent(in) :: ia
      type (mp_real), intent(in) :: b
      double precision da
      mpadd_iq%mpr(1) = mpwork5
      da = ia
      call f_mpadd_d(b%mpr, da, mpadd_iq%mpr)
  end function

  type (mp_real) function mpadd_qd(a, b)
      type (mp_real), intent(in) :: a
      real*8, intent(in) :: b
      mpadd_qd%mpr(1) = mpwork5
      call f_mpadd_d(a%mpr, b, mpadd_qd%mpr)
  end function

  type (mp_real) function mpadd_dq(a, b)
      real*8, intent(in) :: a
      type (mp_real), intent(in) :: b
      mpadd_dq%mpr(1) = mpwork5
      call f_mpadd_d(b%mpr, a, mpadd_dq%mpr)
  end function

! Subtractions
  type (mp_real) function mpsub(a, b)
      type (mp_real), intent(in) :: a, b
      mpsub%mpr(1) = mpwork5
      call f_mpsub(a%mpr, b%mpr, mpsub%mpr)
  end function

  type (mp_real) function mpsub_qi(a, ib)
      type (mp_real), intent(in) :: a
      integer, intent(in) :: ib
      real*8 db
      mpsub_qi%mpr(1) = mpwork5
      db = ib
      call f_mpsub_d(a%mpr, db, mpsub_qi%mpr)
  end function

  type (mp_real) function mpsub_iq(ia, b)
      integer, intent(in) :: ia
      type (mp_real), intent(in) :: b
      real*8 da
      mpsub_iq%mpr(1) = mpwork5
      da = ia
      call f_mpsub_dq(da, b%mpr, mpsub_iq%mpr)
  end function

  type (mp_real) function mpsub_qd(a, b)
      type (mp_real), intent(in) :: a
      real* 8, intent(in) :: b
      mpsub_qd%mpr(1) = mpwork5
      call f_mpsub_d(a%mpr, b, mpsub_qd%mpr)
  end function

  type (mp_real) function mpsub_dq(a, b)
      real* 8, intent(in) :: a
      type (mp_real), intent(in) :: b
      mpsub_dq%mpr(1) = mpwork5
      call f_mpsub_dq(a, b%mpr, mpsub_dq%mpr)
  end function

! Unary Minus
  type (mp_real) function mpneg_q(a)
    type (mp_real), intent(in) :: a
    mpneg_q%mpr(1) = mpwork5
    call f_mpneg_q(a%mpr, mpneg_q%mpr);
  end function

! Multiplications
  type (mp_real) function mpmul(a, b)
      type (mp_real), intent(in) :: a, b
      mpmul%mpr(1) = mpwork5
      call f_mpmul(a%mpr, b%mpr, mpmul%mpr)
  end function

  type (mp_real) function mpmul_qi(a, ib)
      type (mp_real), intent(in) :: a
      integer, intent(in) :: ib
      mpmul_qi%mpr(1) = mpwork5
      call f_mpmul_qi(a%mpr, ib, mpmul_qi%mpr)
  end function

  type (mp_real) function mpmul_iq(ia, b)
      integer, intent(in) :: ia
      type (mp_real), intent(in) :: b
      mpmul_iq%mpr(1) = mpwork5
      call f_mpmul_qi(b%mpr, ia, mpmul_iq%mpr)
  end function

  type (mp_real) function mpmul_qd(a, b)
      type (mp_real), intent(in) :: a
      real*8, intent(in) :: b
      mpmul_qd%mpr(1) = mpwork5
      call f_mpmul_qd(a%mpr, b, mpmul_qd%mpr)
  end function

  type (mp_real) function mpmul_dq(a, b)
      real*8, intent(in) :: a
      type (mp_real), intent(in) :: b
      mpmul_dq%mpr(1) = mpwork5
      call f_mpmul_qd(b%mpr, a, mpmul_dq%mpr)
  end function

! Divisions
  type (mp_real) function mpdiv(a, b)
      type (mp_real), intent(in) :: a, b
      mpdiv%mpr(1) = mpwork5
      call f_mpdiv(a%mpr, b%mpr, mpdiv%mpr)
  end function

  type (mp_real) function mpdiv_qi(a, ib)
      type (mp_real), intent(in) :: a
      integer, intent(in) :: ib
      mpdiv_qi%mpr(1) = mpwork5
      call f_mpdiv_qi(a%mpr, ib, mpdiv_qi%mpr)
  end function

  type (mp_real) function mpdiv_iq(ia, b)
      integer, intent(in) :: ia
      type (mp_real), intent(in) :: b
      mpdiv_iq%mpr(1) = mpwork5
      call f_mpdiv_iq(ia, b%mpr, mpdiv_iq%mpr)
  end function

  type (mp_real) function mpdiv_qd(a, b)
      type (mp_real), intent(in) :: a
      real*8, intent(in) :: b
      mpdiv_qd%mpr(1) = mpwork5
      call f_mpdiv_qd(a%mpr, b, mpdiv_qd%mpr)
  end function

  type (mp_real) function mpdiv_dq(a, b)
      real*8, intent(in) :: a
      type (mp_real), intent(in) :: b
      mpdiv_dq%mpr(1) = mpwork5
      call f_mpdiv_dq(a, b%mpr, mpdiv_dq%mpr)
  end function

! Assignments
  subroutine mpeq_str (qa, ab)
    type (mp_real), intent(out):: qa
    character*(*), intent (in):: ab
    character*1 az(len(ab))
    integer i, l
    qa%mpr(1) = mpwork5
    l = len (ab)
    do i = 1, l
      az(i) = ab(i:i)
    enddo
    call mpdexc (az, l, qa%mpr)
    return
  end subroutine

  subroutine mpeq(a, b)
      type (mp_real), intent(inout) :: a
      type (mp_real), intent(in) :: b
      a%mpr(1) = mpwork5
      call f_mpeq(b%mpr, a%mpr)
  end subroutine

  subroutine mpeq_qi(a, ib)
      type (mp_real), intent(inout) :: a
      integer, intent(in) :: ib
      real*8 db
      a%mpr(1) = mpwork5
      db = ib
      call f_mpeq_d(db, a%mpr)
  end subroutine

  subroutine mpeq_iq(ia, b)
      integer, intent(out) :: ia
      type (mp_real), intent(in) :: b
      double precision db
      integer ib
      call f_mpmdc(b%mpr, db, ib)
      ia = db * 2.d0 ** ib
  end subroutine

  subroutine mpeq_qd(a, b)
      type (mp_real), intent(inout) :: a
      real*8, intent(in) :: b
      a%mpr(1) = mpwork5
      call f_mpeq_d(b, a%mpr)
  end subroutine

  subroutine mpeq_dq(a, b)
      double precision, intent(out) :: a
      type (mp_real), intent(in) :: b
      double precision db
      integer ib
      call f_mpmdc(b%mpr, db, ib)
      a = db * 2.d0 ** ib
  end subroutine

! Equalities
  logical function mpeqt_qq(a, b)
      type (mp_real), intent(in) :: a, b
      integer ic
      call f_mpcpr(a%mpr, b%mpr, ic)
      if (ic .eq. 1) then
         mpeqt_qq = .true.
      else
         mpeqt_qq = .false.
      endif
      return
  end function
      
  logical function mpeqt_qi(a, ib)
      type (mp_real), intent(in) :: a
      integer, intent(in) :: ib
      integer ic
      call f_mpcpr_i(a%mpr, ib, ic)
      if (ic .eq. 1) then
         mpeqt_qi = .true.
      else
         mpeqt_qi = .false.
      endif
      return
  end function
      
  logical function mpeqt_iq(ia, b)
      integer, intent(in) :: ia
      type (mp_real), intent(in) :: b
      integer ic
      call f_mpcpr_i(b%mpr, ia, ic)
      if (ic .eq. 1) then
         mpeqt_iq = .true.
      else
         mpeqt_iq = .false.
      endif
      return
  end function
      
  logical function mpeqt_qd(a, b)
      type (mp_real), intent(in) :: a
      real*8, intent(in) :: b
      integer ic
      call f_mpcpr_d(a%mpr, b, ic)
      if (ic .eq. 1) then
         mpeqt_qd = .true.
      else
         mpeqt_qd = .false.
      endif
      return
  end function

  logical function mpeqt_dq(a, b)
      real*8, intent(in) :: a
      type (mp_real), intent(in) :: b
      integer ic
      call f_mpcpr_d(b%mpr, a, ic)
      if (ic .eq. 1) then
         mpeqt_dq = .true.
      else
         mpeqt_dq = .false.
      endif
      return
  end function

! Inequalities
  logical function mpnet_qq(a, b)
      type (mp_real), intent(in) :: a, b
      integer ic
      call f_mpcpr(a%mpr, b%mpr, ic)
      if (ic .ne. 1) then
         mpnet_qq = .true.
      else
         mpnet_qq = .false.
      endif
      return
  end function
      
  logical function mpnet_qi(a, ib)
      type (mp_real), intent(in) :: a
      integer, intent(in) :: ib
      integer ic
      call f_mpcpr_i(a%mpr, ib, ic)
      if (ic .ne. 1) then
         mpnet_qi = .true.
      else
         mpnet_qi = .false.
      endif
      return
  end function
      
  logical function mpnet_iq(ia, b)
      integer, intent(in) :: ia
      type (mp_real), intent(in) :: b
      integer ic
      call f_mpcpr_i(b%mpr, ia, ic)
      if (ic .ne. 1) then
         mpnet_iq = .true.
      else
         mpnet_iq = .false.
      endif
      return
  end function
      
  logical function mpnet_qd(a, b)
      type (mp_real), intent(in) :: a
      real*8, intent(in) :: b
      integer ic
      call f_mpcpr_d(a%mpr, b, ic)
      if (ic .ne. 1) then
         mpnet_qd = .true.
      else
         mpnet_qd = .false.
      endif
      return
  end function

  logical function mpnet_dq(a, b)
      real*8, intent(in) :: a
      type (mp_real), intent(in) :: b
      integer ic
      call f_mpcpr_d(b%mpr, a, ic)
      if (ic .ne. 1) then
         mpnet_dq = .true.
      else
         mpnet_dq = .false.
      endif
      return
  end function

! Less-Than-Or-Equal-To
  logical function mplet_qq(a, b)
      type (mp_real), intent(in) :: a, b
      integer ic
      call f_mplet(a%mpr, b%mpr, ic)
      if (ic .eq. 1) then
         mplet_qq = .true.
      else
         mplet_qq = .false.
      endif
      return
  end function
      
  logical function mplet_qi(a, ib)
      type (mp_real), intent(in) :: a
      integer, intent(in) :: ib
      integer ic
      call f_mplet_i(a%mpr, ib, ic)
      if (ic .eq. 1) then
         mplet_qi = .true.
      else
         mplet_qi = .false.
      endif
      return
  end function
      
  logical function mplet_iq(ia, b)
      integer, intent(in) :: ia
      type (mp_real), intent(in) :: b
      integer ic
      call f_mpget_i(b%mpr, ia, ic)
      if (ic .eq. 1) then
         mplet_iq = .true.
      else
         mplet_iq = .false.
      endif
      return
  end function
      
  logical function mplet_qd(a, b)
      type (mp_real), intent(in) :: a
      real*8, intent(in) :: b
      integer ic
      call f_mplet_d(a%mpr, b, ic)
      if (ic .eq. 1) then
         mplet_qd = .true.
      else
         mplet_qd = .false.
      endif
      return
  end function

  logical function mplet_dq(a, b)
      real*8, intent(in) :: a
      type (mp_real), intent(in) :: b
      integer ic
      call f_mpget_d(b%mpr, a, ic)
      if (ic .eq. 1) then
         mplet_dq = .true.
      else
         mplet_dq = .false.
      endif
      return
  end function

! Greater-Than-Or-Equal-To
  logical function mpget_qq(a, b)
      type (mp_real), intent(in) :: a, b
      integer ic
      call f_mpget(a%mpr, b%mpr, ic)
      if (ic .eq. 1) then
         mpget_qq = .true.
      else
         mpget_qq = .false.
      endif
      return
  end function

  logical function mpget_qi(a, ib)
      type (mp_real), intent(in) :: a
      integer, intent(in) :: ib
      integer ic
      call f_mpget_i(a%mpr, ib, ic)
      if (ic .eq. 1) then
         mpget_qi = .true.
      else
         mpget_qi = .false.
      endif
      return
  end function
      
  logical function mpget_iq(ia, b)
      integer, intent(in) :: ia
      type (mp_real), intent(in) :: b
      integer ic
      call f_mplet_i(b%mpr, ia, ic)
      if (ic .eq. 1) then
         mpget_iq = .true.
      else
         mpget_iq = .false.
      endif
      return
  end function
      
  logical function mpget_qd(a, b)
      type (mp_real), intent(in) :: a
      real*8, intent(in) :: b
      integer ic
      call f_mpget_d(a%mpr, b, ic)
      if (ic .eq. 1) then
         mpget_qd = .true.
      else
         mpget_qd = .false.
      endif
      return
  end function

  logical function mpget_dq(a, b)
      real*8, intent(in) :: a
      type (mp_real), intent(in) :: b
      integer ic
      call f_mplet_d(b%mpr, a, ic)
      if (ic .eq. 1) then
         mpget_dq = .true.
      else
         mpget_dq = .false.
      endif
      return
  end function

! Less-Than
  logical function mpltt_qq(a, b)
      type (mp_real), intent(in) :: a, b
      integer ic
      call f_mpltt(a%mpr, b%mpr, ic)
      if (ic .eq. 1) then
         mpltt_qq = .true.
      else
         mpltt_qq = .false.
      endif
      return
  end function
      
  logical function mpltt_qi(a, ib)
      type (mp_real), intent(in) :: a
      integer, intent(in) :: ib
      integer ic
      call f_mpltt_i(a%mpr, ib, ic)
      if (ic .eq. 1) then
         mpltt_qi = .true.
      else
         mpltt_qi = .false.
      endif
      return
  end function
      
  logical function mpltt_iq(ia, b)
      integer, intent(in) :: ia
      type (mp_real), intent(in) :: b
      integer ic
      call f_mpgtt_i(b%mpr, ia, ic)
      if (ic .eq. 1) then
         mpltt_iq = .true.
      else
         mpltt_iq = .false.
      endif
      return
  end function
      
  logical function mpltt_qd(a, b)
      type (mp_real), intent(in) :: a
      real*8, intent(in) :: b
      integer ic
      call f_mpltt_d(a%mpr, b, ic)
      if (ic .eq. 1) then
         mpltt_qd = .true.
      else
         mpltt_qd = .false.
      endif
      return
  end function

  logical function mpltt_dq(a, b)
      real*8, intent(in) :: a
      type (mp_real), intent(in) :: b
      integer ic
      call f_mpgtt_d(b%mpr, a, ic)
      if (ic .eq. 1) then
         mpltt_dq = .true.
      else
         mpltt_dq = .false.
      endif
      return
  end function

! Greater-Than
  logical function mpgtt_qq(a, b)
      type (mp_real), intent(in) :: a, b
      integer ic
      call f_mpgtt(a%mpr, b%mpr, ic)
      if (ic .eq. 1) then
         mpgtt_qq = .true.
      else
         mpgtt_qq = .false.
      endif
      return
  end function
      
  logical function mpgtt_qi(a, ib)
      type (mp_real), intent(in) :: a
      integer, intent(in) :: ib
      integer ic
      call f_mpgtt_i(a%mpr, ib, ic)
      if (ic .eq. 1) then
         mpgtt_qi = .true.
      else
         mpgtt_qi = .false.
      endif
      return
  end function
      
  logical function mpgtt_iq(ia, b)
      integer, intent(in) :: ia
      type (mp_real), intent(in) :: b
      integer ic
      call f_mpltt_i(b%mpr, ia, ic)
      if (ic .eq. 1) then
         mpgtt_iq = .true.
      else
         mpgtt_iq = .false.
      endif
      return
  end function
      
  logical function mpgtt_qd(a, b)
      type (mp_real), intent(in) :: a
      real*8, intent(in) :: b
      integer ic
      call f_mpgtt_d(a%mpr, b, ic)
      if (ic .eq. 1) then
         mpgtt_qd = .true.
      else
         mpgtt_qd = .false.
      endif
      return
  end function

  logical function mpgtt_dq(a, b)
      real*8, intent(in) :: a
      type (mp_real), intent(in) :: b
      integer ic
      call f_mpltt_d(b%mpr, a, ic)
      if (ic .eq. 1) then
         mpgtt_dq = .true.
      else
         mpgtt_dq = .false.
      endif
      return
  end function

! Powers
  type (mp_real) function mpexp(a, b)
      type (mp_real), intent(in) :: a, b
      mpexp%mpr(1) = mpwork5
      call f_mppwr(a%mpr, b%mpr, mpexp%mpr)
  end function

  type (mp_real) function mpexp_dq(a, b)
      real*8, intent(in) :: a
      type (mp_real), intent(in) :: b
      type (mp_real) t1
      t1%mpr(1) = mpwork5
      call f_mpdmc (a, t1%mpr)
      mpexp_dq%mpr(1) = mpwork5
      call f_mppwr(t1%mpr, b%mpr, mpexp_dq%mpr)
  end function

  type (mp_real) function mpexp_qd(a, b)
      type (mp_real), intent(in) :: a
      real*8, intent(in) :: b 
      mpexp_qd%mpr(1) = mpwork5
      call f_mppwr_d(a%mpr, b, mpexp_qd%mpr)
  end function

  type (mp_real) function mpexp_qi(a, ib)
      type (mp_real), intent(in) :: a
      integer, intent(in) :: ib 
      mpexp_qi%mpr(1) = mpwork5
      call f_mppwr_qi(a%mpr, ib, mpexp_qi%mpr)
  end function

end module mprealmod


module mpcomplexmod

!  This Fortran-90 module defines operator extensions involving the
!  MP_COMPLEX datatype.  For operations involving two MP data types,
!  those whose first argument is MP_COMPLEX are included here.
!  Others are handled in other modules.
!
!  The subroutines and functions defined in this module are private
!  and not intended to be called directly by the user.

  use mpdefmod
  implicit none

!  MP Complex operator extension interface blocks.

  interface operator (+)
      module procedure mpadd_zj
      module procedure mpadd_jz
      module procedure mpadd_zq
      module procedure mpadd_qz
      module procedure mpadd_zx
      module procedure mpadd_xz
      module procedure mpadd_zz
  end interface

  interface operator (-)
      module procedure mpsub_zj
      module procedure mpsub_jz
      module procedure mpsub_zq
      module procedure mpsub_qz
      module procedure mpsub_zx
      module procedure mpsub_xz
      module procedure mpsub_zz
! negation
      module procedure mpneg_z
  end interface

  interface operator (*)
      module procedure mpmul_zj
      module procedure mpmul_jz
      module procedure mpmul_zq
      module procedure mpmul_qz
      module procedure mpmul_zz
      module procedure mpmul_zd
      module procedure mpmul_dz
  end interface

  interface operator (/)
      module procedure mpdiv_iz
      module procedure mpdiv_zj
      module procedure mpdiv_zq
      module procedure mpdiv_qz
      module procedure mpdiv_zz
      module procedure mpdiv_zd
  end interface

  interface assignment (=)
!      module procedure mpeq_zj
      module procedure mpeq_zq
      module procedure mpeq_zz
      module procedure mpeq_zx
      module procedure mpeq_xz
  end interface

  interface operator (**)
      module procedure mpexp_zi
      module procedure mpexp_zq
!     module procedure mpexp_zz
  end interface

  interface operator (.eq.)
!      module procedure mpeq_tzj
!      module procedure mpeq_tzq
      module procedure mpeqt_zz
  end interface

  interface operator (.ne.)
!      module procedure mpnet_zj
!      module procedure mpnet_zq
      module procedure mpnet_zz
  end interface

contains 

! Additions
  type (mp_complex) function mpadd_zj(za, jb)
      type (mp_complex), intent(in) :: za
      type (mp_integer), intent(in) :: jb
      mpadd_zj%mpc(1) = mpwork5
      mpadd_zj%mpc(mpwds51) = mpwork5
      call f_mpadd_zq(za%mpc, jb%mpi, mpadd_zj%mpc)
  end function

  type (mp_complex) function mpadd_jz(ja, zb)
      type (mp_integer), intent(in) :: ja
      type (mp_complex), intent(in) :: zb
      mpadd_jz%mpc(1) = mpwork5
      mpadd_jz%mpc(mpwds51) = mpwork5
      call f_mpadd_zq(zb%mpc, ja%mpi, mpadd_jz%mpc)
  end function

  type (mp_complex) function mpadd_zq(za, qb)
      type (mp_complex), intent(in) :: za
      type (mp_real), intent(in) :: qb
      mpadd_zq%mpc(1) = mpwork5
      mpadd_zq%mpc(mpwds51) = mpwork5
      call f_mpadd_zq(za%mpc, qb%mpr, mpadd_zq%mpc)
  end function

  type (mp_complex) function mpadd_qz(qa, zb)
      type (mp_real), intent(in) :: qa
      type (mp_complex), intent(in) :: zb
      mpadd_qz%mpc(1) = mpwork5
      mpadd_qz%mpc(mpwds51) = mpwork5
      call f_mpadd_zq(zb%mpc, qa%mpr, mpadd_qz%mpc)
  end function

  type (mp_complex) function mpadd_zx(za, xb)
      type (mp_complex), intent(in) :: za
      complex (kdb), intent(in) :: xb
      double precision br, bi
      mpadd_zx%mpc(1) = mpwork5
      mpadd_zx%mpc(mpwds51) = mpwork5
      br = xb
      bi = aimag(xb)
      call f_mpadd_zx(za%mpc, br, bi, mpadd_zx%mpc)
  end function

  type (mp_complex) function mpadd_xz(xa, zb)
      complex (kdb), intent(in) :: xa
      type (mp_complex), intent(in) :: zb
      double precision ar, ai
      mpadd_xz%mpc(1) = mpwork5
      mpadd_xz%mpc(mpwds51) = mpwork5
      ar = xa
      ai = aimag(xa)
      call f_mpadd_zx(zb%mpc, ar, ai, mpadd_xz%mpc)
  end function

  type (mp_complex) function mpadd_zz(za, zb)
      type (mp_complex), intent(in) :: za, zb
      mpadd_zz%mpc(1) = mpwork5
      mpadd_zz%mpc(mpwds51) = mpwork5
      call f_mpadd_zz(za%mpc, zb%mpc, mpadd_zz%mpc)
  end function

! Subtractions
  type (mp_complex) function mpsub_zj(za, jb)
      type (mp_complex), intent(in) :: za
      type (mp_integer), intent(in) :: jb
      mpsub_zj%mpc(1) = mpwork5
      mpsub_zj%mpc(mpwds51) = mpwork5
      call f_mpsub_zq(za%mpc, jb%mpi, mpsub_zj%mpc)
  end function

  type (mp_complex) function mpsub_jz(ja, zb)
      type (mp_integer), intent(in) :: ja
      type (mp_complex), intent(in) :: zb
      mpsub_jz%mpc(1) = mpwork5
      mpsub_jz%mpc(mpwds51) = mpwork5
      call f_mpsub_qz(ja%mpi, zb%mpc, mpsub_jz%mpc)
  end function

  type (mp_complex) function mpsub_zq(za, qb)
      type (mp_complex), intent(in) :: za
      type (mp_real), intent(in) :: qb
      mpsub_zq%mpc(1) = mpwork5
      mpsub_zq%mpc(mpwds51) = mpwork5
      call f_mpsub_zq(za%mpc, qb%mpr, mpsub_zq%mpc)
  end function

  type (mp_complex) function mpsub_qz(qa, zb)
      type (mp_real), intent(in) :: qa
      type (mp_complex), intent(in) :: zb
      mpsub_qz%mpc(1) = mpwork5
      mpsub_qz%mpc(mpwds51) = mpwork5
      call f_mpsub_qz(qa%mpr, zb%mpc, mpsub_qz%mpc)
  end function

  type (mp_complex) function mpsub_zx(za, xb)
      type (mp_complex), intent(in) :: za
      complex (kdb), intent(in) :: xb
      double precision br, bi
      mpsub_zx%mpc(1) = mpwork5
      mpsub_zx%mpc(mpwds51) = mpwork5
      br = xb
      bi = aimag(xb)
      call f_mpsub_zx(za%mpc, br, bi, mpsub_zx%mpc)
  end function

  type (mp_complex) function mpsub_xz(xa, zb)
      complex (kdb), intent(in) :: xa
      type (mp_complex), intent(in) :: zb
      double precision ar, ai
      mpsub_xz%mpc(1) = mpwork5
      mpsub_xz%mpc(mpwds51) = mpwork5
      ar = xa
      ai = aimag(xa)
      call f_mpsub_xz(ar, ai, zb%mpc, mpsub_xz%mpc)
  end function

  type (mp_complex) function mpsub_zz(za, zb)
      type (mp_complex), intent(in) :: za, zb
      mpsub_zz%mpc(1) = mpwork5
      mpsub_zz%mpc(mpwds51) = mpwork5
      call f_mpsub_zz(za%mpc, zb%mpc, mpsub_zz%mpc)
  end function

!  negation routine.
  type (mp_complex) function mpneg_z (za)
    type (mp_complex), intent(in) :: za
    mpneg_z%mpc(1) = mpwork5
    mpneg_z%mpc(mpwds51) = mpwork5
    call f_mpneg_z(za%mpc, mpneg_z%mpc);
    return
  end function

! Multiplications
  type (mp_complex) function mpmul_zj(za, jb)
      type (mp_complex), intent(in) :: za
      type (mp_integer), intent(in) :: jb
      mpmul_zj%mpc(1) = mpwork5
      mpmul_zj%mpc(mpwds51) = mpwork5
      call f_mpmul_zq(za%mpc, jb%mpi, mpmul_zj%mpc)
  end function

  type (mp_complex) function mpmul_jz(ja, zb)
      type (mp_integer), intent(in) :: ja
      type (mp_complex), intent(in) :: zb
      mpmul_jz%mpc(1) = mpwork5
      mpmul_jz%mpc(mpwds51) = mpwork5
      call f_mpmul_zq(zb%mpc, ja%mpi, mpmul_jz%mpc)
  end function

  type (mp_complex) function mpmul_zq(za, qb)
      type (mp_complex), intent(in) :: za
      type (mp_real), intent(in) :: qb
      mpmul_zq%mpc(1) = mpwork5
      mpmul_zq%mpc(mpwds51) = mpwork5
      call f_mpmul_zq(za%mpc, qb%mpr, mpmul_zq%mpc)
  end function

  type (mp_complex) function mpmul_qz(qa, zb)
      type (mp_real), intent(in) :: qa
      type (mp_complex), intent(in) :: zb
      mpmul_qz%mpc(1) = mpwork5
      mpmul_qz%mpc(mpwds51) = mpwork5
      call f_mpmul_zq(zb%mpc, qa%mpr, mpmul_qz%mpc)
  end function

  type (mp_complex) function mpmul_zz(za, zb)
      type (mp_complex), intent(in) :: za, zb
      mpmul_zz%mpc(1) = mpwork5
      mpmul_zz%mpc(mpwds51) = mpwork5
      call f_mpmul_zz(za%mpc, zb%mpc, mpmul_zz%mpc)
  end function

  type (mp_complex) function mpmul_zd(za, db)
      type (mp_complex), intent(in) :: za
      double precision, intent(in) :: db
      mpmul_zd%mpc(1) = mpwork5
      mpmul_zd%mpc(mpwds51) = mpwork5
      call f_mpmul_zd(za%mpc, db, mpmul_zd%mpc)
  end function

  type (mp_complex) function mpmul_dz(da, zb)
      double precision, intent(in) :: da
      type (mp_complex), intent(in) :: zb
      mpmul_dz%mpc(1) = mpwork5
      mpmul_dz%mpc(mpwds51) = mpwork5
      call f_mpmul_zd(zb%mpc, da, mpmul_dz%mpc)
  end function

! Divisions
  type (mp_complex) function mpdiv_iz(ia, zb)
      integer, intent(in) :: ia
      type (mp_complex), intent(in) :: zb
      double precision t
      mpdiv_iz%mpc(1) = mpwork5
      mpdiv_iz%mpc(mpwds51) = mpwork5
      t = ia
      call f_mpdiv_dz(t, zb%mpc, mpdiv_iz%mpc)
  end function

  type (mp_complex) function mpdiv_zj(za, jb)
      type (mp_complex), intent(in) :: za
      type (mp_integer), intent(in) :: jb
      mpdiv_zj%mpc(1) = mpwork5
      mpdiv_zj%mpc(mpwds51) = mpwork5
      call f_mpdiv_zq(za%mpc, jb%mpi, mpdiv_zj%mpc)
  end function

  type (mp_complex) function mpdiv_zq(za, qb)
      type (mp_complex), intent(in) :: za
      type (mp_real), intent(in) :: qb
      mpdiv_zq%mpc(1) = mpwork5
      mpdiv_zq%mpc(mpwds51) = mpwork5
      call f_mpdiv_zq(za%mpc, qb%mpr, mpdiv_zq%mpc)
  end function

  type (mp_complex) function mpdiv_qz(qa, zb)
      type (mp_real), intent(in) :: qa
      type (mp_complex), intent(in) :: zb
      mpdiv_qz%mpc(1) = mpwork5
      mpdiv_qz%mpc(mpwds51) = mpwork5
      call f_mpdiv_qz(qa%mpr, zb%mpc, mpdiv_qz%mpc)
  end function

  type (mp_complex) function mpdiv_zz(za, zb)
      type (mp_complex), intent(in) :: za, zb
      mpdiv_zz%mpc(1) = mpwork5
      mpdiv_zz%mpc(mpwds51) = mpwork5
      call f_mpdiv_zz(za%mpc, zb%mpc, mpdiv_zz%mpc)
  end function

  type (mp_complex) function mpdiv_zd(za, db)
      type (mp_complex), intent(in) :: za
      double precision, intent(in) :: db
      mpdiv_zd%mpc(1) = mpwork5
      mpdiv_zd%mpc(mpwds51) = mpwork5
      call f_mpdiv_zd(za%mpc, db, mpdiv_zd%mpc)
  end function

! Assignments
  subroutine mpeq_zq(za, qb)
      type (mp_complex), intent(inout) :: za
      type (mp_real), intent(in) :: qb
      za%mpc(1) = mpwork5
      za%mpc(mpwds51) = mpwork5
      call f_mpeq_zq(qb%mpr, za%mpc)
  end subroutine

  subroutine mpeq_zz(za, zb)
      type (mp_complex), intent(inout) :: za
      type (mp_complex), intent(in) :: zb
      za%mpc(1) = mpwork5
      za%mpc(mpwds51) = mpwork5
      call f_mpeq_zz(zb%mpc, za%mpc)
  end subroutine

  subroutine mpeq_zx(za, xb)
      type (mp_complex), intent(inout) :: za
      complex (kdb), intent(in) :: xb
      real*8 dbr, dbi
      za%mpc(1) = mpwork5
      za%mpc(mpwds51) = mpwork5
      dbr = xb
      dbi = aimag(xb)
      call f_mpeq_zx(dbr, dbi, za%mpc)
  end subroutine

  subroutine mpeq_xz(xa, zb)
      complex (kdb), intent(inout) :: xa
      type (mp_complex), intent(in) :: zb
      real*8 dbr, dbi
      integer ibr, ibi
      call f_mpmdc (zb%mpc(1), dbr, ibr)
      call f_mpmdc (zb%mpc(mpwds51), dbi, ibi)
      xa = cmplx (dbr * 2.d0 ** ibr, dbi * 2.d0 ** ibi, kdb)
  end subroutine

! Powers
  type (mp_complex) function mpexp_zi(za, ib)
      type (mp_complex), intent(in) :: za
      integer, intent(in) :: ib
      mpexp_zi%mpc(1) = mpwork5
      mpexp_zi%mpc(mpwds51) = mpwork5
      call f_mppwr_zi(za%mpc, ib, mpexp_zi%mpc)
  end function

  type (mp_complex) function mpexp_zq(za, qb)
      type (mp_complex), intent(in) :: za
      type (mp_real), intent(in) :: qb
      mpexp_zq%mpc(1) = mpwork5
      mpexp_zq%mpc(mpwds51) = mpwork5
      call f_mppwr_zq(za%mpc, qb%mpr, mpexp_zq%mpc)
  end function

! type (mp_complex) function mpexp_zz(za, zb)
!     type (mp_complex), intent(in) :: za
!     type (mp_complex), intent(in) :: zb
!     mpexp_zz%mpc(1) = mpwork5
!     mpexp_zz%mpc(mpwds51) = mpwork5
!     call f_mppwr_zz(za%mpc, zb%mpc, mpexp_zz%mpc)
! end function

! Equalities
  logical function mpeqt_zz(za, zb)
      type (mp_complex), intent(in) :: za, zb
      integer ic
      call f_mpcpr_z(za%mpc, zb%mpc, ic)
      if (ic .eq. 1) then
         mpeqt_zz = .true.
      else
         mpeqt_zz = .false.
      endif
      return
  end function

! Inequalities
  logical function mpnet_zz(za, zb)
      type (mp_complex), intent(in) :: za, zb
      integer ic
      call f_mpcpr_z(za%mpc, zb%mpc, ic)
      if (ic .ne. 1) then
         mpnet_zz = .true.
      else
         mpnet_zz = .false.
      endif
      return
  end function

end module mpcomplexmod


module mpgenmod

!  This Fortran-90 module defines generic functions involving all
!  MP datatypes.

!  The subroutines and functions defined in this module are private
!  and not intended to be called directly by the user.  The generic
!  names (i.e. interface block names) are publicly accessible, though.
!
  use mpdefmod
!  implicit none

! MP generic interface blocks

  interface abs
      module procedure mp_absj
      module procedure mp_absq
      module procedure mp_absz
  end interface

  interface acos
      module procedure mp_acos
  end interface

  interface aimag
      module procedure mp_imag
  end interface

  interface aint
      module procedure mp_aint
  end interface

  interface anint
      module procedure mp_anint
  end interface

  interface arg
      module procedure mp_arg
  end interface

  interface asin
      module procedure mp_asin
  end interface

  interface atan
      module procedure mp_atan
  end interface

  interface atan2
      module procedure mp_atan2
  end interface

  interface bessel
      module procedure mp_bessel
  end interface

  interface besselexp
      module procedure mp_besselexp
  end interface

  interface cmplx
    module procedure mp_cmplx
  end interface

  interface conjg
      module procedure mp_conjg
  end interface

  interface cos
      module procedure mp_cos
      module procedure mp_cosz
  end interface

  interface cosh
      module procedure mp_cosh
  end interface

  interface dble
      module procedure mp_jtod
      module procedure mp_qtod
      module procedure mp_ztod
  end interface

  interface erf
      module procedure mp_erf
  end interface

  interface erfc
      module procedure mp_erfc
  end interface

  interface exp
      module procedure mp_exp
      module procedure mp_expz
  end interface

  interface gamma
      module procedure mp_gamma
  end interface

  interface int
      module procedure mp_int
  end interface

  interface log
      module procedure mp_log
      module procedure mp_logz
  end interface

  interface log10
      module procedure mp_log10
  end interface

  interface max
      module procedure mp_maxj
      module procedure mp_maxq
      module procedure mp_maxj3
      module procedure mp_maxq3
  end interface

  interface min
      module procedure mp_minj
      module procedure mp_minq
      module procedure mp_minj3
      module procedure mp_minq3
  end interface

  interface mod
      module procedure mp_modj
      module procedure mp_modq
  end interface

  interface mpcsshf
      module procedure mp_cssh
  end interface

  interface mpcssnf
      module procedure mp_cssn
  end interface

  interface mpnrtf
      module procedure mp_nrt
  end interface

! Random Number Generator
  interface mpranf
      module procedure mp_rand
  end interface

  interface mpread
      module procedure mp_inpj
      module procedure mp_inpq
      module procedure mp_inpz
  end interface

! Conversion
  interface mpcmpl
    module procedure mp_qqtoz
  end interface

  interface mpint
      module procedure mp_qtoj
!      module procedure mp_ztoj
      module procedure mp_itoj
      module procedure mp_rtoj
      module procedure mp_dtoj
!      module procedure mp_xtoj
      module procedure mp_atoj
  end interface

  interface mpreal
      module procedure mp_jtoq
      module procedure mp_ztoq
      module procedure mp_itoq
      module procedure mp_rtoq
      module procedure mp_dtoq
!      module procedure mp_xtoq
      module procedure mp_atoq
  end interface

  interface mpwrite
      module procedure mp_outj
      module procedure mp_outq
      module procedure mp_outz
  end interface

  interface sign
      module procedure mp_signj
      module procedure mp_signq
  end interface

  interface sin
      module procedure mp_sin
      module procedure mp_sinz
  end interface

  interface sinh
      module procedure mp_sinh
  end interface

  interface sqrt
      module procedure mp_sqrt
      module procedure mp_sqrtz
  end interface

  interface tan
      module procedure mp_tan
  end interface

  interface tanh
      module procedure mp_tanh
  end interface

contains

! Absolute value
  type (mp_integer) function mp_absj(ja)
    type (mp_integer), intent(in) :: ja
    mp_absj%mpi(1) = mpwork5
    call f_mpabs(ja%mpi, mp_absj%mpi)
    call f_ovcheck (mp_absj%mpi)
  end function

  type (mp_real) function mp_absq(a)
    type (mp_real), intent(in) :: a
    mp_absq%mpr(1) = mpwork5
    call f_mpabs(a%mpr, mp_absq%mpr)
  end function

  type (mp_real) function mp_absz(za)
    type (mp_complex), intent(in) :: za
    mp_absz%mpr(1) = mpwork5
    call f_mpabs_z(za%mpc, mp_absz%mpr)
  end function

  type (mp_real) function mp_acos(a)
    type (mp_real), intent(in) :: a
    mp_acos%mpr(1) = mpwork5
    call f_mpacos(a%mpr, mp_acos%mpr)
  end function

  type (mp_real) function mp_imag(za)
    type (mp_complex), intent(in) :: za
    mp_imag%mpr(1) = mpwork5
    call f_mpeq(za%mpc(mpwds51), mp_imag%mpr)
  end function

  type (mp_real) function mp_aint(qa)
      type (mp_real), intent(in) :: qa
      mp_aint%mpr(1) = mpwork5
      call f_mpaint(qa%mpr, mp_aint%mpr)
  end function

  type (mp_real) function mp_anint(qa)
      type (mp_real), intent(in) :: qa
      mp_anint%mpr(1) = mpwork5
      call f_mpnint(qa%mpr, mp_anint%mpr)
  end function

  type (mp_real) function mp_arg(za)
    type (mp_complex), intent(in) :: za
    mp_arg%mpr(1) = mpwork5
    call f_mparg(za%mpc, mp_arg%mpr)
  end function

  type (mp_real) function mp_asin(a)
    type (mp_real), intent(in) :: a
    mp_asin%mpr(1) = mpwork5
    call f_mpasin(a%mpr, mp_asin%mpr)
  end function

  type (mp_real) function mp_atan(a)
    type (mp_real), intent(in) :: a
    mp_atan%mpr(1) = mpwork5
    call f_mpatan(a%mpr, mp_atan%mpr)
  end function

  type (mp_real) function mp_atan2(a, b)
    type (mp_real), intent(in) :: a, b
    mp_atan2%mpr(1) = mpwork5
    call f_mpatan2(a%mpr, b%mpr, mp_atan2%mpr)
  end function

  type (mp_real) function mp_bessel(qa)
    type (mp_real), intent(in) :: qa
    mp_bessel%mpr(1) = mpwork5
    call f_mpbessel(qa%mpr, mp_bessel%mpr)
  end function  

  type (mp_real) function mp_besselexp(qa)
    type (mp_real), intent(in) :: qa
    mp_besselexp%mpr(1) = mpwork5
    call f_mpbesselexp(qa%mpr, mp_besselexp%mpr)
  end function  

  type (mp_complex) function mp_cmplx (qa, qb)
    type (mp_real), intent(in):: qa, qb
    mp_cmplx%mpc(1) = mpwork5
    mp_cmplx%mpc(mpwds51) = mpwork5
    call f_mpeq(qa%mpr, mp_cmplx%mpc)
    call f_mpeq(qb%mpr, mp_cmplx%mpc(mpwds51))
    return
  end function

  type (mp_complex) function mp_conjg(za)
    type (mp_complex), intent(in):: za
    mp_conjg%mpc(1) = mpwork5
    mp_conjg%mpc(mpwds51) = mpwork5
    call f_mpeq(za%mpc, mp_conjg%mpc)
    call f_mpneg_q(za%mpc(mpwds51), mp_conjg%mpc(mpwds51))
    return
  end function

  type (mp_real) function mp_cos(a)
    type (mp_real), intent(in) :: a
    mp_cos%mpr(1) = mpwork5
    call f_mpcos(a%mpr, mp_cos%mpr)
  end function

  type (mp_complex) function mp_cosz(za)
    type (mp_complex), intent(in) :: za
    mp_cosz%mpc(1) = mpwork5
    mp_cosz%mpc(mpwds51) = mpwork5
    call f_mpcos_z(za%mpc, mp_cosz%mpc)
  end function

  type (mp_real) function mp_cosh(a)
    type (mp_real), intent(in) :: a
    mp_cosh%mpr(1) = mpwork5
    call f_mpcosh(a%mpr, mp_cosh%mpr)
  end function

  double precision function mp_jtod(ja)
      type (mp_integer), intent (in):: ja
      call f_mpdble(ja%mpi, mp_jtod)
  end function

  double precision function mp_qtod(qa)
      type (mp_real), intent (in):: qa
      call f_mpdble(qa%mpr, mp_qtod)
  end function

  double precision function mp_ztod(za)
      type (mp_complex), intent (in):: za
      call f_mpdble(za%mpc, mp_ztod)
  end function

  type (mp_real) function mp_exp(qa)
    type (mp_real), intent(in) :: qa
    mp_exp%mpr(1) = mpwork5
    call f_mpexp(qa%mpr, mp_exp%mpr)
  end function

  type (mp_real) function mp_erf(qa)
    type (mp_real), intent(in) :: qa
    mp_erf%mpr(1) = mpwork5
    call f_mperf(qa%mpr, mp_erf%mpr)
  end function

  type (mp_real) function mp_erfc(qa)
    type (mp_real), intent(in) :: qa
    mp_erfc%mpr(1) = mpwork5
    call f_mperfc(qa%mpr, mp_erfc%mpr)
  end function

  type (mp_complex) function mp_expz(za)
    type (mp_complex), intent(in) :: za
    mp_expz%mpc(1) = mpwork5
    mp_expz%mpc(mpwds51) = mpwork5
    call f_mpexp_z(za%mpc, mp_expz%mpc)
  end function

  type (mp_real) function mp_gamma(qa)
    type (mp_real), intent(in) :: qa
    mp_gamma%mpr(1) = mpwork5
    call f_mpgamma(qa%mpr, mp_gamma%mpr)
  end function  

  integer function mp_int(ja)
      type (mp_integer), intent(in) :: ja
      double precision da
      integer ia
      call f_mpmdc(ja%mpi, da, ia)
      mp_int = da * 2.d0 ** ia
  end function

  type (mp_real) function mp_log(qa)
    type (mp_real), intent(in) :: qa
    mp_log%mpr(1) = mpwork5
    call f_mplog(qa%mpr, mp_log%mpr)
  end function

  type (mp_complex) function mp_logz(za)
    type (mp_complex), intent(in) :: za
    mp_logz%mpc(1) = mpwork5
    mp_logz%mpc(mpwds51) = mpwork5
    call f_mplog_z(za%mpc, mp_logz%mpc)
  end function

  type (mp_real) function mp_log10(qa)
    type (mp_real), intent(in) :: qa
    mp_log10%mpr(1) = mpwork5
    call f_mplog10(qa%mpr, mp_log10%mpr)
  end function

  type (mp_integer) function mp_maxj(ja, jb)
      type (mp_integer), intent(in) :: ja, jb
      integer ic
      mp_maxj%mpi(1) = mpwork5
      call f_mpget(ja%mpi, jb%mpi, ic)
      if (ic .eq. 1) then
         call f_mpeq(ja%mpi, mp_maxj%mpi)
      else
         call f_mpeq(jb%mpi, mp_maxj%mpi)
      endif
      call f_ovcheck (mp_maxj%mpi)
  end function

  type (mp_real) function mp_maxq(qa, qb)
      type (mp_real), intent(in) :: qa, qb
      integer ic
      mp_maxq%mpr(1) = mpwork5
      call f_mpget(qa%mpr, qb%mpr, ic)
      if (ic .eq. 1) then
         call f_mpeq(qa%mpr, mp_maxq%mpr)
      else
         call f_mpeq(qb%mpr, mp_maxq%mpr)
      endif
  end function

  type (mp_integer) function mp_maxj3(ja, jb, jc)
      type (mp_integer), intent(in) :: ja, jb, jc
      integer ic
      mp_maxj3%mpi(1) = mpwork5
      call f_mpget(ja%mpi, jb%mpi, ic)
      if (ic .eq. 1) then
         call f_mpeq(ja%mpi, mp_maxj3%mpi)
      else
         call f_mpeq(jb%mpi, mp_maxj3%mpi)
      endif
      call f_mpget(jc%mpi, mp_maxj3%mpi, ic)
      if (ic .eq. 1) call f_mpeq (jc%mpi, mp_maxj3%mpi)
      call f_ovcheck (mp_maxj3%mpi)
  end function

  type (mp_real) function mp_maxq3(qa, qb, qc)
      type (mp_real), intent(in) :: qa, qb, qc
      integer ic
      mp_maxq3%mpr(1) = mpwork5
      call f_mpget(qa%mpr, qb%mpr, ic)
      if (ic .eq. 1) then
         call f_mpeq(qa%mpr, mp_maxq3%mpr)
      else
         call f_mpeq(qb%mpr, mp_maxq3%mpr)
      endif
      call f_mpget(qc%mpr, mp_maxq3%mpr, ic)
      if (ic .eq. 1) call f_mpeq (qc%mpr, mp_maxq3%mpr)
  end function

  type (mp_integer) function mp_minj(ja, jb)
      type (mp_integer), intent(in) :: ja, jb
      integer ic
      mp_minj%mpi(1) = mpwork5
      call f_mplet(ja%mpi, jb%mpi, ic)
      if (ic .eq. 1) then
         call f_mpeq(ja%mpi, mp_minj%mpi)
      else
         call f_mpeq(jb%mpi, mp_minj%mpi)
      endif
      call f_ovcheck (mp_minj%mpi)
  end function

  type (mp_real) function mp_minq(qa, qb)
      type (mp_real), intent(in) :: qa, qb
      integer ic
      mp_minq%mpr(1) = mpwork5
      call f_mplet(qa%mpr, qb%mpr, ic)
      if (ic .eq. 1) then
         call f_mpeq(qa%mpr, mp_minq%mpr)
      else
         call f_mpeq(qb%mpr, mp_minq%mpr)
      endif
  end function

  type (mp_integer) function mp_minj3(ja, jb, jc)
      type (mp_integer), intent(in) :: ja, jb, jc
      integer ic
      mp_minj3%mpi(1) = mpwork5
      call f_mplet(ja%mpi, jb%mpi, ic)
      if (ic .eq. 1) then
         call f_mpeq(ja%mpi, mp_minj3%mpi)
      else
         call f_mpeq(jb%mpi, mp_minj3%mpi)
      endif
      call f_mplet(jc%mpi, mp_minj3%mpi, ic)
      if (ic .eq. 1) call f_mpeq (jc%mpi, mp_minj3%mpi)
      call f_ovcheck (mp_minj3%mpi)
  end function

  type (mp_real) function mp_minq3(qa, qb, qc)
      type (mp_real), intent(in) :: qa, qb, qc
      integer ic
      mp_minq3%mpr(1) = mpwork5
      call f_mplet(qa%mpr, qb%mpr, ic)
      if (ic .eq. 1) then
         call f_mpeq(qa%mpr, mp_minq3%mpr)
      else
         call f_mpeq(qb%mpr, mp_minq3%mpr)
      endif
      call f_mplet(qc%mpr, mp_minq3%mpr, ic)
      if (ic .eq. 1) call f_mpeq (qc%mpr, mp_minq3%mpr)
  end function

  type (mp_integer) function mp_modj (ja, jb)
      type (mp_integer), intent (in):: ja, jb
      mp_modj%mpi(1) = mpwork5
      call f_mpmod(ja%mpi, jb%mpi, mp_modj%mpi)
      call f_ovcheck (mp_modj%mpi)
  end function

  type (mp_real) function mp_modq (qa, qb)
      type (mp_real), intent (in):: qa, qb
      mp_modq%mpr(1) = mpwork5
      call f_mpmod(qa%mpr, qb%mpr, mp_modq%mpr)
  end function

  type (mp_integer) function mp_signj(ja, jb)
      type (mp_integer), intent (in) :: ja, jb
      mp_signj%mpi(1) = mpwork5
      call f_mpeq(ja%mpi, mp_signj%mpi)
      mp_signj%mpi(2) = sign(mp_signj%mpi(2), jb%mpi(2))
      call f_ovcheck (mp_signj%mpi)
  end function
  
  type (mp_real) function mp_signq(qa, qb)
      type (mp_real), intent (in) :: qa, qb
      mp_signq%mpr(1) = mpwork5
      call f_mpeq(qa%mpr, mp_signq%mpr)
      mp_signq%mpr(2) = sign(mp_signq%mpr(2), qb%mpr(2))
  end function

  type (mp_real) function mp_sin(a)
    type (mp_real), intent(in) :: a
    mp_sin%mpr(1) = mpwork5
    call f_mpsin(a%mpr, mp_sin%mpr)
  end function

  type (mp_complex) function mp_sinz(za)
    type (mp_complex), intent(in) :: za
    mp_sinz%mpc(1) = mpwork5
    mp_sinz%mpc(mpwds51) = mpwork5
    call f_mpsin_z(za%mpc, mp_sinz%mpc)
  end function

  type (mp_real) function mp_sinh(a)
    type (mp_real), intent(in) :: a
    mp_sinh%mpr(1) = mpwork5
    call f_mpsinh(a%mpr, mp_sinh%mpr)
  end function

! SQRT, etc.
  type (mp_real) function mp_sqrt(a)
    type (mp_real), intent(in) :: a
    mp_sqrt%mpr(1) = mpwork5
    call f_mpsqrt(a%mpr, mp_sqrt%mpr)
  end function

  type (mp_complex) function mp_sqrtz(za)
    type (mp_complex), intent(in) :: za
    mp_sqrtz%mpc(1) = mpwork5
    mp_sqrtz%mpc(mpwds51) = mpwork5
    call f_mpsqrt_z(za%mpc, mp_sqrtz%mpc)
  end function

  type (mp_real) function mp_tan(a)
    type (mp_real), intent(in) :: a
    mp_tan%mpr(1) = mpwork5
    call f_mptan(a%mpr, mp_tan%mpr)
  end function

  type (mp_real) function mp_tanh(a)
    type (mp_real), intent(in) :: a
    mp_tanh%mpr(1) = mpwork5
    call f_mptanh(a%mpr, mp_tanh%mpr)
  end function

  subroutine mp_cssh (qa, qb, qc)
      type (mp_real), intent (in):: qa
      type (mp_real), intent (out):: qb, qc
      qb%mpr(1) = mpwork5
      qc%mpr(1) = mpwork5
      call f_mpcsshf(qa%mpr, qb%mpr, qc%mpr)
  end subroutine

  subroutine mp_cssn (qa, qb, qc)
      type (mp_real), intent (in):: qa
      type (mp_real), intent (out):: qb, qc
      qb%mpr(1) = mpwork5
      qc%mpr(1) = mpwork5
      call f_mpcssnf(qa%mpr, qb%mpr, qc%mpr)
  end subroutine

  type (mp_real) function mp_nrt(qa, ib)
      type (mp_real), intent (in):: qa
      integer, intent(in) :: ib
      mp_nrt%mpr(1) = mpwork5
      call f_mpnrt(qa%mpr, ib, mp_nrt%mpr)
  end function

  type (mp_real) function mp_rand ()
      mp_rand%mpr(1) = mpwork5
      call f_mprand (mp_rand%mpr)
  end function

! Conversion
  type (mp_complex) function mp_qqtoz (qa, qb)
    type (mp_real), intent(in):: qa, qb
    mp_qqtoz%mpc(1) = mpwork5
    mp_qqtoz%mpc(mpwds51) = mpwork5
    call f_mpeq (qa%mpr(1), mp_qqtoz%mpc(1))
    call f_mpeq (qb%mpr(1), mp_qqtoz%mpc(mpwds51))
    return
  end function

  type (mp_integer) function mp_qtoj (qa)
      type (mp_real), intent(in) :: qa
      mp_qtoj%mpi(1) = mpwork5
      call f_mpeq(qa%mpr, mp_qtoj%mpi)
      call f_ovcheck (mp_qtoj%mpi)
  end function

  type (mp_real) function mp_jtoq (ja)
      type (mp_integer), intent(in) :: ja
      mp_jtoq%mpr(1) = mpwork5
      call f_mpeq(ja%mpi, mp_jtoq%mpr)
  end function

  type (mp_real) function mp_ztoq (za)
      type (mp_complex), intent(in) :: za
      mp_ztoq%mpr(1) = mpwork5
      call f_mpeq(za%mpc, mp_ztoq%mpr)
  end function

  type (mp_integer) function mp_itoj (ia)
      integer, intent(in) :: ia
      double precision da
      da = ia
      mp_itoj%mpi(1) = mpwork5
      call f_mpdmc(da, mp_itoj%mpi)
      call f_ovcheck (mp_itoj%mpi)
  end function

  type (mp_integer) function mp_rtoj(ra)
      real, intent(in) :: ra
      double precision da
      da = ra
      mp_rtoj%mpi(1) = mpwork5
      call f_mpdmc(da, mp_rtoj%mpi)
      call f_ovcheck (mp_rtoj%mpi)
  end function

  type (mp_integer) function mp_dtoj (da)
      real*8, intent(in) :: da
      mp_dtoj%mpi(1) = mpwork5
      call f_mpdmc(da, mp_dtoj%mpi)
      call f_ovcheck (mp_dtoj%mpi)
  end function

  type (mp_integer) function mp_atoj (ab)
    character*(*), intent (in):: ab
    type (mp_integer) jtmp1, jtmp2
    character*1 az(len(ab))
    integer i, l
    mp_atoj%mpi(1) = mpwork5
    l = len (ab)
    do i = 1, l
      az(i) = ab(i:i)
    enddo
    call mpdexc (az, l, jtmp1%mpi)
    call f_mpinfr (jtmp1%mpi, mp_atoj%mpi, jtmp2%mpi)
    return
  end function

  type (mp_real) function mp_itoq (ia)
      integer, intent(in) :: ia
      double precision da
      da = ia
      mp_itoq%mpr(1) = mpwork5
      call f_mpdmc(da, mp_itoq%mpr)
  end function

  type (mp_real) function mp_rtoq(ra)
      real, intent(in) :: ra
      double precision da
      da = ra
      mp_rtoq%mpr(1) = mpwork5
      call f_mpdmc(da, mp_rtoq%mpr)
  end function

  type (mp_real) function mp_dtoq (da)
      real*8, intent(in) :: da
      mp_dtoq%mpr(1) = mpwork5
      call f_mpdmc(da, mp_dtoq%mpr)
  end function

  type (mp_real) function mp_atoq (ab)
    character*(*), intent (in):: ab
    character*1 az(len(ab))
    integer i, l
    mp_atoq%mpr(1) = mpwork5
    l = len (ab)
    do i = 1, l
      az(i) = ab(i:i)
    enddo
    call mpdexc (az, l, mp_atoq%mpr)
    return
  end function

! Input
  subroutine mp_inpj (iu, j1, j2, j3, j4, j5, j6, j7, j8, j9)
    integer, intent(in) :: iu
    type (mp_integer), intent (out) :: j1, j2, j3, j4, j5, j6, j7, j8, j9
    optional :: j2, j3, j4, j5, j6, j7, j8, j9

    call mpinpj (iu, j1)
    if (present (j2)) call mpinpj (iu, j2)
    if (present (j3)) call mpinpj (iu, j3)
    if (present (j4)) call mpinpj (iu, j4)
    if (present (j5)) call mpinpj (iu, j5)
    if (present (j6)) call mpinpj (iu, j6)
    if (present (j7)) call mpinpj (iu, j7)
    if (present (j8)) call mpinpj (iu, j8)
    if (present (j9)) call mpinpj (iu, j9)
    return
  end subroutine

  subroutine mpinpj (iu, ja)
    integer, intent(in) :: iu
    type (mp_integer), intent(out) :: ja
    type (mp_integer) jtmp1
    type (mp_real) qa

    call mpinpq (iu, qa)
    call f_mpinfr (qa%mpr, ja%mpi, jtmp1%mpi)
    return
  end subroutine

  subroutine mp_inpq (iu, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: iu
    type (mp_real), intent (out) :: q1, q2, q3, q4, q5, q6, q7, q8, q9
    optional :: q2, q3, q4, q5, q6, q7, q8, q9

    call mpinpq (iu, q1)
    if (present (q2)) call mpinpq (iu, q2)
    if (present (q3)) call mpinpq (iu, q3)
    if (present (q4)) call mpinpq (iu, q4)
    if (present (q5)) call mpinpq (iu, q5)
    if (present (q6)) call mpinpq (iu, q6)
    if (present (q7)) call mpinpq (iu, q7)
    if (present (q8)) call mpinpq (iu, q8)
    if (present (q9)) call mpinpq (iu, q9)
    return
  end subroutine

  subroutine mpinpq (iu, qa)
    integer, intent(in) :: iu
    type (mp_real), intent(out) :: qa
    character*1 az(new_mpipl+100)
    integer i, l, l1, nn
    character*200 line

    l = 0
    nn = new_mpipl + 100
    qa%mpr(1) = mpwork5

100 continue
    read(iu, '(a)', end = 200) line

    do i = 200, 1, -1
      if (line(i:i) /= ' ') goto 110
    enddo

    i = 0

110 continue
    l1 = i

    do i = 1, l1
      if (line(i:i) == ',') goto 150
      if (l + 1 <= nn) then
        l = l + 1
        az(l) = line(i:i)
      endif
    enddo

    goto 100

150 continue

    call mpinpc (az, l, qa%mpr)
    goto 300

200 continue

    call mpgetpar ('mpker', i, 72)
    if (i > 0) then
      write (mpldb, 1)
1     format (&
       'mpinpq: end of file or no comma terminating multiprecion input.')
      call mpsetpar ('mpier', 72)
      if (i == 2) stop
    endif

300 continue
    return
  end subroutine

  subroutine mp_inpz (iu, z1, z2, z3, z4, z5, z6, z7, z8, z9)
    integer, intent(in) :: iu
    type (mp_complex), intent (out) :: z1, z2, z3, z4, z5, z6, z7, z8, z9
    optional :: z2, z3, z4, z5, z6, z7, z8, z9

    call mpinpz (iu, z1)
    if (present (z2)) call mpinpz (iu, z2)
    if (present (z3)) call mpinpz (iu, z3)
    if (present (z4)) call mpinpz (iu, z4)
    if (present (z5)) call mpinpz (iu, z5)
    if (present (z6)) call mpinpz (iu, z6)
    if (present (z7)) call mpinpz (iu, z7)
    if (present (z8)) call mpinpz (iu, z8)
    if (present (z9)) call mpinpz (iu, z9)
    return
  end subroutine

  subroutine mpinpz (iu, z)
    integer, intent(in) :: iu
    type (mp_complex), intent(out):: z
    type (mp_real) q1, q2

    call mpinpq (iu, q1)
    call mpinpq (iu, q2)
    z%mpc(1) = mpwork5
    z%mpc(mpwds51) = mpwork5
    call f_mpeq (q1%mpr(1), z%mpc(1))
    call f_mpeq (q2%mpr(1), z%mpc(mpwds51))
    return
  end subroutine

! Output
  subroutine mp_outj (iu, j1, j2, j3, j4, j5, j6, j7, j8, j9)
      integer, intent(in) :: iu
      type (mp_integer), intent (in) :: j1, j2, j3, j4, j5, j6, j7, j8, j9
      optional :: j2, j3, j4, j5, j6, j7, j8, j9
      call mpoutj (iu, j1)
      if (present (j2)) call mpoutj (iu, j2)
      if (present (j3)) call mpoutj (iu, j3)
      if (present (j4)) call mpoutj (iu, j4)
      if (present (j5)) call mpoutj (iu, j5)
      if (present (j6)) call mpoutj (iu, j6)
      if (present (j7)) call mpoutj (iu, j7)
      if (present (j8)) call mpoutj (iu, j8)
      if (present (j9)) call mpoutj (iu, j9)
      return
  end subroutine

  subroutine mpoutj(iu, ja)
      integer, intent(in) :: iu
      type (mp_integer), intent(in) :: ja
      character*1 az(new_mpipl+100)
      integer l
      call mpoutc (ja%mpi, az, l)
      az(l+1) = ','
      write(iu, '(78A1)') (az(i), i = 1, l+1)
  end subroutine

  subroutine mp_outq (iu, q1, q2, q3, q4, q5, q6, q7, q8, q9)
      integer, intent(in) :: iu
      type (mp_real), intent(in) :: q1, q2, q3, q4, q5, q6, q7, q8, q9
      optional :: q2, q3, q4, q5, q6, q7, q8, q9
      call mpoutq (iu, q1)
      if (present (q2)) call mpoutq (iu, q2)
      if (present (q3)) call mpoutq (iu, q3)
      if (present (q4)) call mpoutq (iu, q4)
      if (present (q5)) call mpoutq (iu, q5)
      if (present (q6)) call mpoutq (iu, q6)
      if (present (q7)) call mpoutq (iu, q7)
      if (present (q8)) call mpoutq (iu, q8)
      if (present (q9)) call mpoutq (iu, q9)
      return
  end subroutine

  subroutine mpoutq(iu, q)
      integer, intent(in) :: iu
      type (mp_real), intent(in) :: q
      character*1 az(new_mpipl+100)
      integer l
      call mpoutc (q%mpr, az, l)
      az(l+1) = ','
      write(iu, '(78A1)') (az(i), i = 1, l+1)
  end subroutine

  subroutine mp_outz (iu, z1, z2, z3, z4, z5, z6, z7, z8, z9)
      integer, intent(in) :: iu
      type (mp_complex), intent(in) :: z1, z2, z3, z4, z5, z6, z7, z8, z9
      optional :: z2, z3, z4, z5, z6, z7, z8, z9
      call mpoutz (iu, z1)
      if (present (z2)) call mpoutz (iu, z2)
      if (present (z3)) call mpoutz (iu, z3)
      if (present (z4)) call mpoutz (iu, z4)
      if (present (z5)) call mpoutz (iu, z5)
      if (present (z6)) call mpoutz (iu, z6)
      if (present (z7)) call mpoutz (iu, z7)
      if (present (z8)) call mpoutz (iu, z8)
      if (present (z9)) call mpoutz (iu, z9)
      return
  end subroutine

  subroutine mpoutz(iu, z)
      integer, intent(in) :: iu
      type (mp_complex), intent(in) :: z
      character*1 az(new_mpipl+100)
      integer l
      call mpoutc (z%mpc(1), az, l)
      az(l+1) = ','
      write(iu, '(78A1)') (az(i), i = 1, l+1)
      call mpoutc (z%mpc(mpwds51), az, l)
      az(l+1) = ','
      write(iu, '(78A1)') (az(i), i = 1, l+1)
  end subroutine

end module mpgenmod


module mpmodule
use mpintmod
use mprealmod
use mpcomplexmod
use mpgenmod
end module mpmodule
