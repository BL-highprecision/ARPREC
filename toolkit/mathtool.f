subroutine f_main
! program mathtool

!   The Mathematician's Toolkit
!   Version date:  2008-08-13
!   Copyright 2002, 2003, 2004, 2005, 2006, 2007, 2008

!   Author: David H Bailey (LBNL)
!   Contact email:  dhbailey@lbl.gov

!   This work was supported by the Director, Office of Science, Division
!   of Mathematical, Information, and Computational Sciences of the
!   U.S. Department of Energy under contract number DE-AC02-05CH11231.

use mpmodule
use globdata
implicit none
integer i, i1, i2, ires, j, k, lfname, lnam1, lnam2, lnblk, lstr, lstx, &
  lst1, lst2, ndig, nres, nunit
parameter (lstx = 2048)
character*1 st1(lstx), stout(lstx+100)
character*16 lcase, strx1, strx2
character*40 fname
character*160 str1, str2
character*2048 string
type (mp_real) res(ntmpx)
double precision tm0, tm1, second
external lnblk, lcase

write (6, 1)
1 format ('Welcome to the Experimental Mathematician''s Toolkit'/&
  'Initializing...')
call mpinit (ndigmx2)
call toolinit
ires = 0
nunit = 5
lfname = 0
fname = ' '

100 continue

i1 = lnblk (quadn(nquadt))
write (6, 2) ndebug, ndigits1, ndigits2, nepsilon1, nepsilon2, npslqb, npslql, &
  nquadl, nquadt, quadn(nquadt)(1:i1)
2 format (/'Current settings:'/&
  'Debug level =',i4/&
  'Primary precision level =',i6,' digits'/&
  'Secondary precision level =',i6,' digits'/&
  'Primary epsilon = 10^',i6/&
  'Secondary epsilon = 10^',i6/&
  'PSLQ bound =',i4/&
  'PSLQ level =',i4/&
  'Quadrature level =',i4/&
  'Quadrature type =',i4,'  (',a,')')
write (6, 3)
3 format (/'Sample expressions (case insensitive):'/ &
  'e + pi + log2 + log10 + catalan  Adds these pre-defined constants.'/&
  'result[1] + result[2]            Adds result #1 to result #2.'/&
  'alpha = arctan[3] - 3*log[2]     Defines or sets user variable alpha.'/&
  'fun1[x,y] = 2*sqrt[x]*erf[y]     Defines user function fun1.'/&
  'clear[nam1]                      Clears definition of variable or function.'/&
  'integrate[1/(1+x^2), {x, 0, 1}]  Integrates 1/(1+x^2) from x=0 to 1.'/&
  'sum[1/2^k, {k, 0, infinity}]     Sums 1/2^k from k=0 to infinity.'/&
  'binomial[20,10]*factorial[10]    Evaluates binomial coeff and factorial.'/&
  'zeta[3] + zetaz[1,1,2]           Evaluates zeta and multi-zeta functions.'/&
  'table[x^k, {k, 1, 4}]            Forms the list [x^1, x^2, x^3, x^4].'/&
  'pslq[table[x^k, {k, 0, n}]]      Finds coeffs of degree-n poly for x.'/&
  'polyroot[1,-1,-1,{0.618}]        Finds real root of 1-x-x^2=0 near 0.618.'/&
  'polyroot[1,2,3,{-0.33, 0.47}]    Complex root of 1+2x+3x^2 near -0.33+0.47i.'/&
  'digits = 200                     Sets working precision to 200 digits.'/&
  '                                 This does NOT change display precision.'/&
  'epsilon = -190                   Sets epsilon level to 10^(-190).'/&
  'eformat[190,180]                 Display using E format with 180 digits.'/&
  'fformat[60,50]                   Display using F format with 50 digits.'/&
  'input file.dat                   Inputs commands from file file.dat.'/&
  'output file.dat                  Outputs user vars and funs to file.dat.'/&
  'help polyroot                    Displays a brief explanation of polyroot.'/&
  'functions                        Displays a list of all defined functions.'/&
  'variables                        Displays a list of all defined variables.'/&
  'prompt                           Displays this message.'/&
  'exit                             Exits this program.'/&
  'Expressions may be continued on next line by typing \ at end of line.')

!   Input expression string.

110 continue

nerror = 0
lstr = 0
string = ' '
write (6, '()')

120 continue

read (nunit, '(a160)', end = 300) str1
lst1 = lnblk (str1)

if (lst1 == 0 .and. lstr == 0) then

!   Null string.

  goto 110
elseif (lstr + lst1 > lstx) then

!   Input string exceeds maximum length.

  write (6, 4) lstx
4 format ('Combined input string too long; max length = ',i5)
  goto 110
elseif (str1(lst1:lst1) == '\') then

!   Input string is continued on next line.

  string(lstr+1:lstr+lst1-1) = str1(1:lst1-1)
  lstr = lstr + lst1 - 1
  goto 120
elseif (str1 == 'prompt') then

!   Prompt.

  goto 100
elseif (str1 == 'exit') then

!   Exit

  goto 400
elseif (str1(1:5) == 'input') then

!  Input definitions from file.

  do i = 6, lst1
    if (str1(i:i) /= ' ') goto 130
  enddo

130 continue

  lfname = lst1 - i + 1
  fname = str1(i:lst1)
  if (lfname <= 0 .or. lfname > 40 .or. fname == ' ') then
    write (6, 5)
5   format ('Input file name missing or too long; max 40 chars.')
    goto 110
  endif
  nunit = 10
  open (nunit, file = fname, status = 'old', form = 'formatted', err = 140)
  rewind (nunit)
  goto 110

140 continue

  write (6, 6) fname(1:lfname)
6 format ('Input file does not exist: ',a)
  nunit = 5
  goto 110
elseif (str1(1:6) == 'output') then

!   Output variables and function definitions to file.

  do i = 7, lst1  
    if (str1(i:i) /= ' ') goto 150
  enddo
  
150 continue 
  
  lfname = lst1 - i + 1 
  fname = str1(i:lst1) 
  if (lfname <= 0 .or. lfname > 40 .or. fname == ' ') then 
    write (6, 7)  
7   format ('Output file name missing or too long; max 40 chars.') 
    goto 110 
  endif 
  open (10, file = fname, form = 'formatted')
  rewind 10

!   Output variable names and numerical values to file.

  do i = nvara + nvarb + nvarc + 1, nvar
    lnam1 = lnblk (varn(i))
    lst2 = lnam1
    str2(1:lst2) = varn(i)(1:lst2)
    str2(lst2+1:lst2+4) = ' = \'
    lst2 = lst2 + 4
    write (10, '(a)') str2(1:lst2)
    call mpeform (var(i), ndigits1+20, ndigits1, stout)
    lnam2 = ndigits1 + 20

    do k = 0, lnam2 - 1, 80
      lst2 = min (lnam2 - k, 80)
      str2 = ' '

      do j = 1, lst2
        str2(j:j) = stout(j+k)
      enddo

      if (k + lst2 < lnam2 - 1) then
        str2(lst2+1:lst2+1) = '\'
        lst2 = lst2 + 1
      endif
      write (10, '(a)') str2(1:lst2)
    enddo
  enddo

!   Output function definitions to file.

  do i = nfuna + 1, nfun
    lnam1 = lnblk (funn(i))
    lnam2 = lnblk (fund(i))
    str2(1:lnam1) = funn(i)(1:lnam1)
    str2(lnam1+1:lnam1+1) = '['
    lst2 = lnam1 + 1

    do j = 1, narg(i)
      write (strx1, '("arg",i1,",")') j
      str2(lst2+1:lst2+5) = strx1(1:5)
      lst2 = lst2 + 5
    enddo

    str2(lst2:lst2+4) = '] = \'
    lst2 = lst2 + 4
    write (10, '(a)') str2(1:lst2)
    lst2 = 0
    str2 = ' '

    do k = 0, lnam2 - 1, 80
      lst2 = min (lnam2 - k, 80)
      str2(1:lst2) = fund(i)(k+1:k+lst2)
      if (k + lst2 < lnam2 - 1) then
        str2(lst2+1:lst2+1) = '\'
        lst2 = lst2 + 1
      endif
      write (10, '(a)') str2(1:lst2)
    enddo
  enddo

  rewind 10
  close (10)
  write (6, 8) fname(1:lfname)
8 format ('User variables and function definitions output to file: ',a)
  goto 110 
elseif (str1(1:5) == 'clear') then

!   Clear a function or variable name.

  i1 = index (str1, ']')
  if (str1(6:6) /= '[' .or. i1 == 0) then
    write (6, 9)
9   format ('Syntax eror in clear command.')
    goto 110
  endif
  strx1 = str1(7:i1-1)

  do j = nvara + nvarb + nvarc + 1, nvar
    if (varn(j) == strx1) then
      do i = j, nvar - 1
        varn(i) = varn(i+1)
        var(i) = var(i+1)
      enddo

      nvar = nvar - 1
      write (6, 10) strx1(1:16)
10    format ('Variable name removed: ',a)
      goto 110
    endif
  enddo

  do j = nfuna + 1, nfun
    if (funn(j) == strx1) then
      do i = j, nfun - 1
        narg(i) = narg(i+1)
        funn(i) = funn(i+1)
        fund(i) = fund(i+1)
      enddo

      nfun = nfun - 1
      write (6, 11) strx1(1:16)
11    format ('Function name removed: ',a)
      goto 110
    endif
  enddo

  write (6, 12) strx1
12 format ('No user variable or function by this name: ',a)
  goto 110
elseif (str1(1:7) == 'eformat') then

!   Switch to E Format.

  i1 = index (str1, ',')
  i2 = index (str1, ']')
  if (str1(8:8) /= '[' .or. i1 == 0 .or. i2 == 0 .or. i2 < i1 .or. lst1 > i2) &
    then
    write (6, 13)
13  format ('Syntax error in eformat command.')
    goto 110
  endif
  strx1 = str1(9:i1-1)
  strx2 = str1(i1+1:i2-1)
  read (strx1, '(i16)') i1
  read (strx2, '(i16)') i2
  if (i1 < 0 .or. i1 > lstx .or. i2 > ndigits1 .or. i2 < 0 .or. i2 > i1 - 5) then
    write (6, 14) i1, i2
14  format ('Improper arguments in eformat command:',2i7)
    goto 110
  endif
  nef = 1
  nef1 = i1
  nef2 = i2
  goto 110
elseif (str1(1:7) == 'fformat') then

!   Switch to F Format.

  i1 = index (str1, ',')
  i2 = index (str1, ']')
  if (str1(8:8) /= '[' .or. i1 == 0 .or. i2 == 0 .or. i2 < i1 .or. lst1 > i2) &
    then
    write (6, 15)
15  format ('Syntax error in fformat command.')
    goto 110
  endif
  strx1 = str1(9:i1-1)
  strx2 = str1(i1+1:i2-1)
  read (strx1, '(i16)') i1
  read (strx2, '(i16)') i2
  if (i1 < 0 .or. i1 > lstx .or. i2 > ndigits1 .or. i2 < 0 .or. i2 > i1 - 5) then
    write (6, 16) i1, i2
16  format ('Improper arguments in fformat command:',2i7)
    goto 110
  endif
  nef = 2
  nef1 = i1
  nef2 = i2
  goto 110
elseif (str1(1:4) == 'help') then

  do i = 5, lst1
    if (str1(i:i) /= ' ') goto 160
  enddo

160 continue

  i1 = lst1 - i + 1
  strx1 = str1(i:lst1)
  if (i1 <= 0 .or. strx1 == ' ') strx1 = '[blank]'
  i2 = 0

  do j = 1, nfuna
    if (lcase (strx1) == lcase (funn(j))) then
      i2 = j

      do i = 1, 4
        if (funhelp(i,j) /= ' ') write (6, '(a)') funhelp(i,j)
      enddo
    endif
  enddo

  if (i2 == 0) write (6, 17) strx1, (funn(i), i = 1, nfuna)
17 format ('This is not a valid predefined function name: ',a/ &
     'Help information available for:'/ (a,' ',a,' ',a,' ',a))
   goto 110
elseif (str1(1:9) == 'variables') then
  write (6, '(a," ",a," ",a," ",a)') (varn(i), i = 1, nvar)
  goto 110
elseif (str1(1:9) == 'functions') then
  write (6, '(a," ",a," ",a," ",a)') (funn(i), i = 1, nfun)
  goto 110
endif

!   Form full input string.

if (lst1 > 0) then
  string(lstr+1:lstr+lst1) = str1(1:lst1)
  lstr = lstr + lst1
endif

!   Check if string is function name without brackets.

if (lstr <= 16) then
  do j = nfuna + 1, nfun
    if (string(1:lstr) == funn(j)) then
      i1 = lnblk (fund(j))

      do i = 1, i1, 76
        if (i + 76 < i1) then
          write (6, '(a,"\")') fund(j)(i:i+75)
        else
          write (6, '(a)') fund(j)(i:i1)
        endif
      enddo

      goto 110
    endif

  enddo
endif

!   Parse input string.

call mpsetpar ('mpier', 0)
tm0 = second ()
call parse (string, nres, res)
tm1 = second ()

if (nerror > 0) then

!   Error code resulted during execution.

  write (6, 18) nerror
18 format ('Error; code =',i3)
else

!   Output results, and post results to % and to Result array.

  if (nres > 0) then
    if (ires + nres > mxres) ires = 0
    if (nres == 1) then
      write (6, 19) ires + 1
19    format ('Result[',i3,'] =')
    elseif (nres > 1) then
      write (6, 20) ires + 1, ires + nres
20    format ('Result[',i3,'] through Result[',i3,'] =')
    endif

    do i = 1, nres
      ires = ires + 1
      result(ires) = res(i)
      if (nef == 1) then
        call mpeform (res(i), nef1, nef2, stout)
      elseif (nef == 2) then
        call mpfform (res(i), nef1, nef2, stout)
      elseif (nef == 3) then
        if (abs (res(i)) > 0.0001d0 .and. abs (res(i)) < 10000.d0) then
          call mpfform (res(i), nef1, nef2, stout)
        else
          call mpeform (res(i), nef1, nef2, stout)
        endif
      endif
      lnam2 = nef1

      do k = 0, lnam2 - 1, 80
        lst2 = min (lnam2 - k, 80)
        str2 = ' '

        do j = 1, lst2
          str2(j:j) = stout(j+k)
        enddo

        if (k + lst2 < lnam2 - 1) then
          str2(lst2+1:lst2+1) = '\'
          lst2 = lst2 + 1
        endif
        write (6, '(a)') str2(1:lst2)
      enddo
    enddo
  endif
  if (ndebug > 0) then
    write (6, 21) tm1 - tm0
21  format ('CPU time =',f12.4)
  endif
endif

goto 110

!   End of file for input data from file.

300 continue

write (6, 22) fname(1:lfname)
22 format ('Input commands read from file: ',a)
rewind (10)
close (10)
nunit = 5
lfname = 0
fname = ' '
goto 110

400 continue

stop
end

subroutine toolinit

!   This subroutine initializes data for the Toolkit.

use mpmodule
use globdata
implicit none
integer i, i1, i2, ndeb0, nq0, nq3
type (mp_real) catalan, eulergamma, t1, t2, t3
parameter (ndeb0 = 2, nq0 = 3)

!   Set basic parameters.

i1 = 1

do i = 1, 72
  call mpsetpar ('mpker', i1, i)
enddo

nfun = nfuna
nvar = nvara + nvarb + nvarc
ndebug = ndeb0
ndigits1 = min (ndigi, ndigmx1)
ndigits2 = min (2 * ndigits1, ndigmx2)
nepsilon1 = -ndigits1
nepsilon2 = -ndigits2
mpoud = ndigits1 + 20
nquadl = int (log (dble (ndigits1)) / log (2.d0))
npslqb = 100
npslql = 1
nquadt = nq0
nef = 3
nef1 = 78
nef2 = 68

!   Read constants and quadrature arrays from files.

open (11, file = 'const.dat', form = 'unformatted')
rewind (11)
read (11) t1, t2
rewind (11)
close (11)

if (nquadt == 1) then
  open (11, file = 'quadgs.dat', form = 'unformatted')
elseif (nquadt == 2) then
  open (11, file = 'quaderf.dat', form = 'unformatted')
elseif (nquadt == 3) then
  open (11, file = 'quadts.dat', form = 'unformatted')
endif
rewind (11)
read (11) nq3
read (11) (quadwk(i), i = -1, nq3)
read (11) (quadxk(i), i = -1, nq3)
rewind (11)
close (11)

!  Initialize var array with basic constants.

t3 = 1.d0
var(1) = exp (t3)
var(2) = mpl02
var(3) = mpl10
var(4) = mppic
var(5) = t1
var(6) = t2
var(7) = mplrg

do i = 1, nvarb
  var(nvara+i) = 0.d0
enddo

do i = nvara + nvarb + 1, nvara + nvarb + nvarc
  if (varn(i) == 'Debug') then
    var(i) = dble (ndebug)
  elseif (varn(i) == 'Digits') then
    var(i) = dble (ndigits1)
  elseif (varn(i) == 'Digits2') then
    var(i) = dble (ndigits2)
  elseif (varn(i) == 'Epsilon') then
    var(i) = dble (nepsilon1)
  elseif (varn(i) == 'Epsilon2') then
    var(i) = dble (nepsilon2)
  elseif (varn(i) == 'Pslqbound') then
    var(i) = dble (npslqb)
  elseif (varn(i) == 'Pslqlevel') then
    var(i) = dble (npslql)
  elseif (varn(i) == 'Quadlevel') then
    var(i) = dble (nquadl)
  elseif (varn(i) == 'Quadtype') then
    var(i) = dble (nquadt)
  else
    write (6, 1) 
1   format ('toolinit: impossible case -- please contact author.')
    stop
  endif
enddo

do i = 1, nvarz
  var(nvara+nvarb+nvarc+i) = 0.d0
enddo

!   Set precision to ndigits1.  This must be done last, so as to insure
!   that the constants in var (see above) are stored to full precision.

call mpsetprec (ndigits2)
call mpgetprecwords (nwords2)
call mpsetprec (ndigits1)
call mpgetprecwords (nwords1)

return
end

recursive subroutine parse (string, nres, res)

!   This subroutine parses and processes the input character string.

use mpmodule
use globdata
implicit none
integer i, inam, inum, ioper, iop(8), iprec(0:7), iasgn, i1, i2, i3, j, &
  k, ks, k1, k2, lnam, lnam1, lnam2, lnamx, lnblk, lnum, lnum1, lnum2, &
  lnumx, lstr, nop, nres, ntmp1, ntmp2
character*27 alphal, alphau
character*11 digits
character*13 delim
parameter (lnamx = 16, lnumx = 2048, &
  alphal = 'abcdefghijklmnopqrstuvwxyz_', &
  alphau = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ_', &
  digits = '0123456789.', delim = ' +-*/^:=,()[]')
character*2048 num, num1, num2
character*2048 string, str1, str2, str3
character*16 lcase, nam, nam1, nam2, argn(nvarb)
character*1 char
type (mp_real) res(ntmpx), tmp1(ntmpx), tmp2(ntmpx), t1, t2, t3
type (mp_real) t300
external lnblk, lcase
data iprec /3, 1, 1, 2, 2, 4, 0, 0/

t300 = 1.d300
lstr = lnblk (string)

! write (6, *) 'parse: enter; string ='
! write (6, *) string(1:lstr)

ioper = 0
iasgn = 0
k = 0
lnam = 0
lnum = 0
nam = ' '
num = ' '
nop = 0
nres = 0
ntmp1 = 0
ntmp2 = 0

100 continue

k = k + 1
if (k > lstr) goto 200
char = string(k:k)

!   Check for alphabetic.

! write (6, *) 'parse: k, lstr, char =', k, lstr, '#'//char//'#'
! write (6, *) 'nop, iop array =', nop, (iop(i), i = 1, nop), &
!   (delim(iop(i)+1:iop(i)+1), i = 1, nop)
! write (6, *) 'nres, res array =', nres
! write (6, '(1p,d20.10)') (dble(min(max(res(i),-t300),t300)), i=1,nres)

i1 = index (alphau, char)
i2 = index (alphal, char)
i3 = max (i1, i2)
if (i3 > 0 .and. (char /= 'e' .or. lnum == 0)) then
  if (ioper /= 0) then
    write (6, 1) 
1   format ('parse: operator expected, variable name found.'/ &
      'If multiplication is intended, use * between operands.')
    nerror = 1
    goto 400
  endif
  if (lnum > 0) then
    write (6, 2) 
2   format ('parse: alphabetic character in digit string.')
    nerror = 2
    goto 400
  endif
  lnam = lnam + 1
  if (lnam > lnamx) then
    write (6, 3) lnamx
3   format ('parse: name too long; max chars =',i4)
    nerror = 3
    goto 400
  endif
  nam(lnam:lnam) = char
  goto 100
endif

!   Check for numeric.

i1 = index (digits, char)

if (i1 > 0 .or. i3 > 0) then
  if (ioper /= 0) then
    write (6, 4) 
4   format ('parse: operator expected, numeric constant found.'/ &
      'If multiplication is intended, use * between operands.')
    nerror = 4
    goto 400
  endif
  if (lnam > 0) then
    if (i1 == 11) then
      write (6, 31) nam(1:lnam)
31    format ('parse: period after name =',a)
      goto 400
    endif
    lnam = lnam + 1
    if (lnam > lnamx) then
      write (6, 3) lnamx
      nerror = 5
      goto 400
    endif
    nam(lnam:lnam) = char
  else
    lnum = lnum + 1
    if (lnum > lnumx) then
      write (6, 5) lnumx
5     format ('parse: digit string too long; max chars =',i6)
      nerror = 6
      goto 400
    endif
    num(lnum:lnum) = char

!   If this is an 'e' (in an e-format number), check if next char is '-'.

    if (char == 'e' .and. k < lstr .and. lnum < lnumx) then
      if (string(k+1:k+1) == '-') then
        k = k + 1
        lnum = lnum + 1
        num(lnum:lnum) = string(k:k)
      endif
    endif
  endif
  goto 100
endif

200 continue

!   Check if this is end of name.

if (lnam > 0) then

!   Look for brackets containing dummy argument list.

  call findbrackets (k, string, k1, k2)
  if (nerror > 0) goto 400
  if (k1 == 0) then

!   Check table of defined variables.

    do j = 1, nvar
      if (lcase (nam) == lcase (varn(j))) then

!   Add variable's value to res array.

        nres = nres + 1
        if (nres > ntmpx) then
          write (6, 6) ntmpx
6         format ('parse: vector too long; max =',i4)
          nerror = 7
          goto 400
        endif
        res(nres) = var(j)
        if (nres == 1) iasgn = j
        lnam = 0
        nam = ' '
        ioper = 1
        goto 300
      endif
    enddo

!   New name is not found in variable name table.  Check if it appears in 
!   function name table.

    do i = 1, nfun
      if (lcase (nam) == lcase (funn(i))) then
        write (6, 7) nam
7       format ('parse: function name is not followed by brackets = ',a)
        nerror = 8
        goto 400
      endif
    enddo

!   If new name appears after first variable or numeric string, it is undefined.

    if (nres > 0) then
      write (6, 8) nam
8     format ('parse: variable name not found: ',a/ &
      'Currently defined variable names =')
      write (6, '(4a18)') (varn(i), i = 1, nvar)
      nerror = 9
      goto 400
    endif

!   Add new name to variable table.

    nvar = nvar + 1
    if (nvar > nvarx) then
      write (6, 9) nvarx
9     format ('parse: too many variable names; max =',i4)
      nerror = 10
      goto 400
    endif
    var(nvar) = 0.d0
    varn(nvar) = nam
    nres = 1
    res(1) = 0.d0
    iasgn = nvar
    write (6, 10) nam
10  format ('parse: new variable name = ',a)
    lnam = 0
    nam = ' '
    ioper = 1

!   Make sure that new variable name is followed by = or :=

    do i = k, lstr - 1
      if (string(i:i) /= ' ') then
        if (string(i:i) /= '=' .and. string(i:i+1) /= ':=') then
          goto 210
        else
          goto 300
        endif
      endif
    enddo

210  continue

    write (6, 11) 
11  format ('parse: new variable name is not followed by = or :=')
    nerror = 11
    goto 400
  else

!   Pair of brackets found -- check table of defined functions for name.

    do j = 1, nfun
      if (lcase (nam) == lcase (funn(j))) then

!   Evaluate function arguments, then evaluate function.

        str1 = string(k1+1:k2-1)
        call evalfun (j, str1, tmp1, ntmp2, tmp2)
        if (nerror > 0) goto 400

!   Append function value(s) to res array.

        if (nres + ntmp2 > ntmpx) then
          write (6, 6) ntmpx
          nerror = 12
          goto 400
        endif

        do i = 1, ntmp2
          nres = nres + 1
          res(nres) = tmp2(i)
        enddo

        k = k2 + 1
        ioper = ntmp2
        if (k > lstr) goto 300
        char = string(k:k)
        lnam = 0
        nam = ' '
        goto 300
      endif
    enddo

!   New name is not found in function name table.  Check if it appears in 
!   variable name table.

    do i = 1, nvar
      if (lcase (nam) == lcase (varn(i))) then
        write (6, 13) nam
13      format ('parse: variable name is followed by brackets = ',a)
        nerror = 13
        goto 400
      endif
    enddo

!   If new name appears after first variable or numeric string, it is undefined.

    if (nres > 0) then
      write (6, 14) nam
14    format ('parse: function name not found: ',a/ &
      'Currently defined function names =')
      write (6, '(4a18)') (funn(i), i = 1, nfun)
      nerror = 14
      goto 400
    endif

!   Add new name to function table.

    nfun = nfun + 1
    if (nfun > nfunx) then
      write (6, 15) nfunx
15    format ('parse: too many function names; max =',i4)
      nerror = 15
      goto 400
    endif
    funn(nfun) = nam
    write (6, 16) nam
16  format ('parse: new function name = ',a)
    k = k1 + 1

!   Make sure the new function is followed by a valid dummy argument list.

    do j = 1, nvarb
      do i = k, k2
        if (string(i:i) /= ' ') goto 220
      enddo

220   continue

      k = i
      if (k == k2) then
        write (6, 17) nvarb
17      format ('parse: syntax error in dummy argument list; max args =',i2)
        nerror = 16
        goto 400
      endif

      i1 = index (string(k:k2), ',')
      if (i1 == 0) i1 = k2 - k + 1
      i2 = k - 1 + i1
      nam1 = string (k:i2-1)
      i3 = lnblk (nam1)

      do i = 1, i3
        if (index (delim, nam1(i:i)) > 0) then
          write (6, 17) nvarb
          nerror = 17
          goto 400
        endif
      enddo

      argn(j) = nam1
      k = i2 + 1
      if (k >= k2) goto 240
    enddo

    write (6, 17) nvarb
    nerror = 18
    goto 400

240 continue

    narg(nfun) = j
    k = k + 1
    
!   Make sure that new function argument list is followed by = or :=

    do i = k, lstr - 1
      if (string(i:i) /= ' ') then
        if (string(i:i) == '=') then
          k = i + 1
          goto 260
        elseif (string(i:i+1) == ':=') then
          k = i + 2
          goto 260
        else
          goto 250
        endif
      endif
    enddo

250  continue

    write (6, 19)
19  format ('parse: new function name is not followed by = or :=')
    nerror = 19
    goto 400

260 continue

!   Replace dummy argument names in function definition with arg1, arg2, etc.


    do i = k, lstr
      if (string(i:i) /= ' ') goto 270
    enddo

270 continue

    k = i
    str1 = string(k:lstr)

    do i = 1, j
      write (nam1, '("arg",i1)') i
      call replace (argn(i), nam1, str1, str2)
      str1 = str2
    enddo

!   Insert remainder of input string into function definition array.

    fund(nfun) = str1
    goto 400
  endif
elseif (lnum > 0) then

!   End of digit string.  Append value to res array.

  nres = nres + 1
  if (nres > ntmpx) then
    write (6, 6) ntmpx
    nerror = 20
    goto 400
  endif
  res(nres) = num(1:lnum)
  lnum = 0
  num = ' '
  ioper = 1
  goto 300
endif

300 continue

! write (6, *) 'parse: after 300:'
! write (6, *) 'nop, iop array =', nop, (iop(i), i = 1, nop), &
!   (delim(iop(i)+1:iop(i)+1), i = 1, nop)
! write (6, '(1p,d20.10)') (dble(min(max(res(i),-t300),t300)), i=1,nres)

!   Check for end-of-line -- if so, clear pending arithmetic operation(s).

if (k > lstr) then

!   write (6, *) 'parse: end of line'

  if (ioper == 0) then
    write (6, 20)
20  format ('parse: premature end of line.')
    nerror = 21
    goto 400
  elseif (ioper > 1 .and. nop > 0) then
    write (6, 102)
102 format ('parse: illegal operator with array result.')
    nerror = 22
    goto 400
  endif

  do i = nop, 1, -1
    i1 = iop(i)
    call oper (i1, iasgn, nres, res)
    if (nerror > 0) goto 400
  enddo

  goto 400
endif

!   Check for delimiters.

i1 = index (delim, char)

!   Handle other delimiters.

if (i1 == 0) then
  write (6, 21) char
21 format ('parse: illegal character =',a)
  nerror = 23
  goto 400
elseif (i1 == 1) then

!   Delimiter is a blank.

  goto 100
elseif (ioper == 0 .and. (i1 == 2 .or. i1 >= 4 .and. i1 <= 9)) then

!   Consecutive operators.

  write (6, 22)
22 format ('parse: illegal operator syntax.')
  nerror = 24
  goto 400
elseif (i1 >= 2 .and. i1 <= 8) then
  if (ioper > 1) then
    write (6, 22)
    nerror = 25
    goto 400
  endif

! write (6, *) 'parse: start of arith op block:'
! write (6, *) 'nop, iop array =', nop, (iop(i), i = 1, nop), &
!   (delim(iop(i)+1:iop(i)+1), i = 1, nop)
! write (6, '(1p,d20.10)') (dble(min(max(res(i),-t300),t300)), i=1,nres)

!   Delimiter is arithmetic operator or equal.

  if (i1 == 7) then
    if (k == lstr .or. string(k:k+1) /= ':=') then
      write (6, 23)
23    format ('parse: colon is not followed by equal sign.')
      nerror = 26
      goto 400
    endif
  endif
  if (ioper == 0 .and. i1 == 3) then

!   Unary minus case.

    i2 = 0
  else
    i2 = i1 - 1
  endif

!   Compare precedence with pending operations, and clear as needed.

  do i = nop, 1, -1
    if (iprec(iop(i)) < iprec(i2)) goto 350
    i3 = iop(i)
    call oper (i3, iasgn, nres, res)
    if (nerror > 0) goto 400
  enddo

  i = 0

350 continue

  nop = i + 1
  iop(nop) = i2
  if (i2 == 6) k = k + 1
  ioper = 0

! write (6, *) 'parse: end of arith op block:'
! write (6, *) 'nop, iop array =', nop, (iop(i), i = 1, nop), &
!   (delim(iop(i)+1:iop(i)+1), i = 1, nop)
! write (6, '(1p,d20.10)') (dble(min(max(res(i),-t300),t300)), i=1,nres)

  goto 100
elseif (i1 == 9) then

!   Delimiter is a comma -- clear pending arithmetic operation(s).

  do i = nop, 1, -1
    i1 = iop(i)
    call oper (i1, iasgn, nres, res)
    if (nerror > 0) goto 400
  enddo

  nop = 0
  ioper = 0
  goto 100
elseif (i1 == 10) then

!   Delimiter is left parentheses.

  if (ioper > 0) then
    write (6, 24)
24  format ('parse: illegal parenthesis syntax.')
    nerror = 27
    goto 400
  endif

!   find matching right parenthesis.

  call findpars (k, string, k1, k2)
  if (nerror > 0) goto 400

!   Evaluate expression inside parentheses.

  str1 = string(k1+1:k2-1)

!   write (6, *) 'parse: before recursive call to parse'

  call parse (str1, ntmp1, tmp1)

!   write (6, *) 'parse: after recursive call to parse'

  if (nerror > 0) goto 400
  if (ntmp1 > 1) then
    write (6, 25)
25  format ('parse: expression inside parentheses has multiple values.')
    nerror = 28
    goto 400
  endif

!   Append value to res array.

  nres = nres + 1
  if (nres > ntmpx) then
    write (6, 6) ntmpx
    nerror = 29
    goto 400
  endif
  res(nres) = tmp1(1)
  k = k2
  ioper = 1
  goto 100
elseif (i1 == 11) then

!   Delimiter is right parenthesis.

  write (6, 26)
26 format ('parse: mismatched right parenthesis.')
  nerror = 30
  goto 400
elseif (i1 == 12) then

!   Delimiter is left bracket.

  write (6, 27)
27 format ('parse: illegal left bracket.')
  nerror = 31
  goto 400
elseif (i1 == 13) then

!   Delimiter is right bracket.

  write (6, 28)
28 format ('parse: illegal right bracket.')
  nerror = 32
  goto 400
else
  write (6, 29) char
29 format ('parse: illegal delimiter =',a)
  nerror = 33
  goto 400
endif

400 continue

! write (6, *) 'parse: exit; nerror, nres =', nerror, nres
! write (6, '(1p,d20.10)') (dble(min(max(res(i),-t300),t300)), i=1,nres)

return
end

subroutine oper (ix, iasgn, nres, res)

!   This evaluates arithmetic operations.

use mpmodule
use globdata
implicit none
integer i, i1, i2, ix, iasgn, k, k1, lnblk, lnm, ierror, nq3, nres
character*16 delim, lcase, nam
parameter (delim = ' +-*/^:=,()[]')
type (mp_real) res(ntmpx), t1, t2
! type (mp_real) t300
external lnblk, lcase

! t300 = 1.d300

!   Handle the five basic arithmetic operations.

! write (6, *) 'oper input: ix, iasgn, symbol =', ix, iasgn, delim(ix+1:ix+1)
! write (6, *) 'nres, res array =', nres
! write (6, '(1p,d20.10)') (dble(min(max(res(i),-t300),t300)), i=1,nres)

if (ix == 0) then
  if (nres <= 0) goto 100
  res(nres) = - res(nres)
  goto 200
elseif (ix == 1) then
  if (nres <= 1) goto 100
  t1 = res(nres-1) + res(nres)
  nres = nres - 1
  res(nres) = t1
  goto 200
elseif (ix == 2) then
  if (nres <= 1) goto 100
  t1 = res(nres-1) - res(nres)
  nres = nres - 1
  res(nres) = t1
  goto 200
elseif (ix == 3) then
  if (nres <= 1) goto 100
  t1 = res(nres-1) * res(nres)
  nres = nres - 1
  res(nres) = t1
  goto 200
elseif (ix == 4) then
  if (nres <= 1) goto 100
  t1 = res(nres-1) / res(nres)
  nres = nres - 1
  res(nres) = t1
  goto 200
elseif (ix == 5) then
  if (nres <= 1) goto 100
  t1 = aint (res(nres))
  if (t1 == res(nres) .and. abs (t1) < 1.d9) then
    k = t1
    t1 = res(nres-1) ** k
  else
    t1 = res(nres-1) ** res(nres)
  endif
  nres = nres - 1
  res(nres) = t1
  goto 200
elseif (ix == 6 .or. ix == 7) then

!   Assignment.

  if (nres <= 1 .or. nres > 2) then
    goto 100
  elseif (iasgn <= 0) then
    write (6, 1) iasgn
1   format ('oper: illegal assignment; iasgn =',i3)
    nerror = 41
  elseif (iasgn <= nvara + nvarb) then
    write (6, 2)  varn(iasgn)
2   format ('oper: illegal assignment to reserved name = ',a)
    nerror = 42
  elseif (iasgn <= nvara + nvarb + nvarc) then

!   Variable is one of the special names: Debug, Digits, Epsilon, etc.

    nam = lcase (varn(iasgn))
    lnm = lnblk (nam)
    t1 = res(2)
    if (abs (t1) > 1.d6) then
      write (6, 3) varn(iasgn)(1:lnm)
3     format ('oper: improper value assigned to ',a)
      nerror = 43
      goto 200
    endif
    var(iasgn) = res(2)
    nres = 0
    k1 = dble (var(iasgn))

    select case (nam)
    case ('debug')
      if (k1 < 0 .or. k1 > 10) then
        write (6, 4) 0, 10
4       format ('oper: improper value assigned to Debug; min, max =',2i8)
        nerror = 44
        goto 200
      endif
      ndebug = k1
      var(iasgn) = dble (ndebug)
    case ('digits')
      if (k1 < 30 .or. k1 > ndigmx1) then
        write (6, 5) 30, ndigmx1
5       format ('oper: improper value assigned to Digits; min, max =',2i8)
        nerror = 45
        goto 200
      elseif (k1 > ndigits1 .and. nvar > nvara + nvarb + nvarc) then
        write (6, 6)
6       format ('oper: user variables should be recalculated before use.')
      endif
      ndigits1 = k1
      var(iasgn) = dble (ndigits1)
      call mpsetprec (ndigits1)
      mpoud = ndigits1 + 20
      ndigits2 = min (2 * ndigits1, ndigmx2)
      nepsilon1 = - ndigits1
      nepsilon2 = - ndigits2
      nquadl = int (log (dble (ndigits1)) / log (2.d0))
      call mpsetprec (ndigits2)
      call mpgetprecwords (nwords2)
      call mpsetprec (ndigits1)
      call mpgetprecwords (nwords1)

      do i = nvara + nvarb + 1, nvara + nvarb + nvarc
        if (varn(i) == 'Digits2') then
          var(i) = dble (ndigits2)
        elseif (varn(i) == 'Epsilon') then
          var(i) = dble (nepsilon1)
        elseif (varn(i) == 'Epsilon2') then
          var(i) = dble (nepsilon2)
        elseif (varn(i) == 'Quadlevel') then
          var(i) = dble (nquadl)
        endif
      enddo

      write (6, 7) ndigits1, ndigits2, nepsilon1, nepsilon2, nquadl
7     format ('oper: Digits changed to ',i5/'Digits2 changed to ',i5/ &
        'Epsilon changed to ',i5/'Epsilon2 changed to ',i5/ &
        'Quadlevel changed to ',i2)
    case ('digits2')
      if (k1 < 30 .or. k1 > ndigmx2) then
        write (6, 25) 30, ndigmx2
25       format ('oper: improper value assigned to Digits2; min, max =',2i8)
        nerror = 53
        goto 200
      endif
      ndigits2 = k1
      var(iasgn) = dble (ndigits2)
      nepsilon2 = - ndigits2
      call mpsetprec (ndigits2)
      call mpgetprecwords (nwords2)
      call mpsetprec (ndigits1)

      do i = nvara + nvarb + 1, nvara + nvarb + nvarc
        if (varn(i) == 'Epsilon2') var(i) = dble (nepsilon2)
      enddo

      write (6, 27) ndigits2, nepsilon2
27     format ('oper: digits2 changed to ',i5/'epsilon2 changed to ',i5)
    case ('epsilon')
      if (k1 > -10) then
        write (6, 8) -10
8       format ('oper: improper value assigned to Epsilon; max =',i8)
        nerror = 46
        goto 200
      elseif (k1 < - ndigits1) then
        write (6, 9) ndigits1
9       format ('oper: warning - Epsilon exceeds -Digits; Digits =',i5)
      endif
      nepsilon1 = k1
      var(iasgn) = dble (nepsilon1)
    case ('epsilon2')
      if (k1 > -10) then
        write (6, 8) k1
28       format ('oper: improper value assigned to Epsilon2; max =',i8)
        nerror = 54
        goto 200
      elseif (k1 < - ndigits2) then
        write (6, 9) ndigits2
29       format ('oper: warning - Epsilon2 exceeds -Digits2; Digits2 =',i5)
      endif
      nepsilon2 = k1
      var(iasgn) = dble (nepsilon2)
    case ('pslqbound')
      if (k1 < 1 .or. k1 > 1000) then
        write (6, 10) 1, 1000
10      format ('oper: improper value assigned to Pslqbound; min, max =',2i8)
        nerror = 47
        goto 200
      endif
      npslqb = k1
      var(iasgn) = dble (npslqb)
    case ('pslqlevel')
      if (k1 < 1 .or. k1 > 3) then
        write (6, 11) 1, 3
11      format ('oper: improper value assigned to Pslqlevel; min, max =',2i8)
        nerror = 48
        goto 200
      endif
      npslql = k1
      var(iasgn) = dble (npslql)
    case ('quadlevel')
      if (k1 < 2 .or. k1 > nquadx) then
        write (6, 12) 2, nquadx
12      format ('oper: improper value assigned to QuadLevel; min, max =',2i8)
        nerror = 49
        goto 200
      endif
      nquadl = k1
      var(iasgn) = dble (nquadl)
    case ('quadtype')
      if (k1 < 1 .or. k1 > 3) then
        write (6, 13) 1, 3
13      format ('oper: improper value assigned to Quadtype; min, max =',2i8)
        nerror = 50
        goto 200
      endif
      nquadt = k1
      var(iasgn) = dble (nquadt)
      if (ndebug > 0) write (6, 14) quadn(k1)
14    format ('New quadrature type: ',a)
      if (k1 == 1) then
        open (11, file = 'quadgs.dat', form = 'unformatted')
      elseif (k1 == 2) then
        open (11, file = 'quaderf.dat', form = 'unformatted')
      elseif (k1 == 3) then
        open (11, file = 'quadts.dat', form = 'unformatted')
      else
        write (6, 15)
15      format ('oper: impossible case; please contact author.')
        stop
      endif
      rewind (11)
      read (11, end = 110) nq3
      read (11, end = 110) (quadwk(i), i = -1, nq3)
      read (11, end = 110) (quadxk(i), i = -1, nq3)
    case default
      goto 100
    end select
    goto 200
  elseif (iasgn <= nvar) then
    var(iasgn) = res(2)
    nres = 0
    goto 200
  else
    write (6, 1) iasgn
    nerror = 51
  endif
  goto 200
else
  goto 100
endif

100 continue

write (6, 16) ix, nres
16 format ('oper: stack error or assignment error; ix, nres =',2i4,';'/&
  'please contact author.')
stop

110 continue

write (6, 17)
17 format (&
  'End-of-file encountered when reading quadrature data. Most likely the'/ &
  'quadrature type selected has not been initialized.  Select another type.')
nerror = 52

200  continue

call mpgetpar ('mpier', ierror)
if (ierror > 0) nerror = ierror + 1000

! write (6, *) 'oper exit:'
! write (6, *) 'nres, res array =', nres
! write (6, '(1p,d20.10)') (dble(min(max(res(i),-t300),t300)), i=1,nres)

210 continue

return
end

recursive subroutine evalfun (ix, string, tmp, nvalue, value)

!   This evaluates function references.

use mpmodule
use globdata
implicit none
integer i, idb, iq, ix, izeta(100), k, k1, k2, lnblk, ierror, ntmp, nvalue
character*16 lcase, nam
character*2048 string, stx1
type (mp_real) binomial, factorial, qinteg, qsum, tmp(ntmpx), t1, t2, t3, &
  value(ntmpx), zetap, zetaz
external binomial, factorial, lcase, lnblk, qinteg, qsum, zetap, zetaz

! write (6, *) 'evalfun: enter; ix, name =', ix, funn(ix)
! write (6, *) 'string =', string(1:64)

nvalue = 1
nam = lcase (funn(ix))

if (nam /= 'polyroot' .and. nam /= 'integrate' .and. nam /= 'sum' &
  .and. nam /= 'table') then
  call parse (string, ntmp, tmp)
  if (nerror > 0) goto 300
  if (ntmp /= narg(ix) .and. narg(ix) /= 1000) then
    write (6, 1) nam, narg(ix)
1   format ('parse: wrong number of arguments for function ',a,i2)
    nerror = 61
    goto 300
  endif
endif

select case (nam)
case ('abs')
  value(1) = abs (tmp(1))
case ('arccos')
  value(1) = acos (tmp(1))
case ('arcsin')
  value(1) = asin (tmp(1))
case ('arctan')
  value(1) = atan (tmp(1))
case ('arctan2')
  value(1) = atan2 (tmp(1), tmp(2))
case ('bessel')
  value(1) = bessel (tmp(1))
case ('besselexp')
  value(1) = besselexp (tmp(1))
case ('binomial')
  value(1) = binomial (tmp(1), tmp(2))
case ('cos')
  value(1) = cos (tmp(1))
case ('erf')
  value(1) = erf (tmp(1))
case ('exp')
  value(1) = exp (tmp(1))
case ('factorial')
  value(1) = factorial (tmp(1))
case ('gamma')
  value(1) = gamma (tmp(1))
case ('integrate')
  value(1) = qinteg (string, tmp(1))
case ('log')
  value(1) = log (tmp(1))
case ('max')
  value(1) = max (tmp(1), tmp(2))
case ('min')
  value(1) = min (tmp(1), tmp(2))
case ('polyroot')
  call polyroot (string, tmp(1), nvalue, value)
case ('pslq')
  call qpslq (ntmp, tmp(1), string, value)
  nvalue = ntmp
case ('result')
  k = dble (tmp(1))
  if (k <= 0 .or. k > mxres) then
    write (6, 3) dble (tmp(1)), mxres
3   format ('evalfun: illegal argument to Result =',1pd15.6/ &
      'Must be integer between zero and ',i4)
    nerror = 62
    goto 300
  endif
  value(1) = result(k)
case ('sin')
  value(1) = sin (tmp(1))
case ('sqrt')
  value(1) = sqrt (tmp(1))
case ('sum')
  value(1) = qsum (string, tmp(1))
case ('table')
  call table (string, tmp(1), nvalue, value)
case ('tan')
  value(1) = tan (tmp(1))
case ('zeta')
  do i = 1, ntmp
    t1 = tmp(i)
    if (i > 1 .or. t1 > 100.d0 .or. t1 /= aint (t1)) then
      write (6, 4) 1, 100
4     format ('evalfun: improper argument to Zeta; integer, min, max =',2i8)
      call mpwrite (6, t1)
      nerror = 63
      goto 300
    endif
    izeta(i) = t1
  enddo

  nvalue = 1
  value(1) = zetaz (ntmp, izeta)
  if (nerror > 0) goto 300
case ('zetap')
  do i = 1, ntmp
    t1 = tmp(i)
    if (t1 > 100.d0 .or. t1 /= aint (t1)) then
      write (6, 5) 100
5     format ('evalfun: improper argument to Zetap; integer, max =',i8)
      call mpwrite (6, t1)
      nerror = 64
      goto 300
    endif
    izeta(i) = t1
  enddo

  nvalue = 1
  value(1) = zetap (ntmp, izeta)
  if (nerror > 0) goto 300
case ('zetaz')
  do i = 1, ntmp
    t1 = tmp(i)
    if (t1 > 100.d0 .or. t1 /= aint (t1)) then
      write (6, 6) 100
6     format ('evalfun: improper argument to Zetaz; integer, max =',i8)
      call mpwrite (6, t1)
      nerror = 65
      goto 300
    endif
    izeta(i) = t1
  enddo

  nvalue = 1
  value(1) = zetaz (ntmp, izeta)
  if (nerror > 0) goto 300
case default
  if (ix <= nfuna) then
    write (6, 7) 
7   format ('evalfun: impossible case; contact author.')
    stop
  endif

!   Evaluate one of the user-defined functions.

  stx1 = fund(ix)

  do i = 1, ntmp
    var(i+nvara) = tmp(i)
  enddo

  call parse (stx1, nvalue, value)
  if (nerror > 0) goto 300
end select

call mpgetpar ('mpier', ierror)
if (ierror > 0) nerror = ierror + 1000

300 continue

! write (6, *) 'evalfun: exit; nvalue, value array =', nvalue
! write (6, '(1p,4d19.12)') (dble (value(i)), i = 1, nvalue)

return
end
  
function binomial (a1, a2)

!   This evaluates the binomial function.  Computation is accelerated by
!   accumulating products in double precision before updating the multiprecision
!   value.

use mpmodule
use globdata
type (mp_real) a1, a2, binomial, t1, t2, t3
real*8 d1, d2, t53
parameter (t53 = 2.d0 ** 53)
integer i, k1, k2, n1, n2

if (a1 == 0.d0 .and. a2 == 0.d0) then
  binomial = 1.d0
  goto 100
elseif (a1 < 1.d0 .or. a1 > 1.d5 .or. a1 /= aint (a1) .or. a2 < 1.d0 &
  .or. a2 > 1.d5 .or. a2 /= aint (a2) .or. a2 > a1) then
  call decmd (a1, d1, n1)
  call decmd (a2, d2, n2)
  write (6, 1) d1, n1, d2, n2
1 format ('binomial: improper arguments; values =', &
    f10.6,' x 10^',i5,3x,f10.6,'x 10^',i5)
  nerror = 71
  goto 100
endif

k1 = a1
k2 = min (a2, a1 - a2)
t1 = 1.d0
d1 = 1.d0

do k = 1, k2
  d2 = (k1 + 1 - k) * d1
  if (d2 > t53) then
    t1 = d1 * t1
    d2 = k1 + 1 - k
  endif
  d1 = d2
enddo

t1 = d1 * t1
d1 = 1.d0

do k = 1, k2
  d2 = k * d1
  if (d2 > t53) then
    t1 = t1 / d1
    d2 = k
  endif
  d1 = d2
enddo

binomial = t1 / d1

100 continue

! write (6, *) 'binomial: exit'

return
end

function factorial (a1)

!   This evaluates the factorial function.  Computation is accelerated by
!   accumulating products in double precision before updating the multiprecision
!   value.

use mpmodule
use globdata
type (mp_real) a1, factorial, t1, t2, t3
real*8 d1, d2, t53
parameter (t53 = 2.d0 ** 53)
integer i, k1, n1

! write (6, *) 'factorial: enter'

if (a1 == 0.d0) then
  factorial = 1.d0
  goto 100
elseif (a1 < 1.d0 .or. a1 > 1.d5 .or. a1 /= aint (a1)) then
  call decmd (a1, d1, n1)
  write (6, 1) d1, n1
1 format ('factorial: improper argument; value =', &
    f10.6,' x 10^',i5,3x,f10.6,'x 10^',i5)
  nerror = 81
  goto 100
endif

k1 = anint (a1)
t1 = 1.d0
d1 = 1.d0

do k = 1, k1
  d2 = k * d1
  if (d2 > t53) then
    t1 = d1 * t1
    d2 = k
  endif
  d1 = d2
enddo

factorial = d1 * t1

100 continue

! write (6, *) 'factorial: exit'

return
end

recursive subroutine polyroot (string, tmp1, nvalue, value)

!   This handles calls to polyroot, which finds real or complex roots of
!   polynomials.

use mpmodule
use globdata
implicit none
integer i, j, k, k1, k2, lnblk, lstr, lx1, lx2, ierror, ntmp1, ntmp2, nr, nvalue
character*2048 string, stx1, stx2
type (mp_real) tmp1(ntmpx), tmp2(ntmpx), t1, value(ntmpx)
external lcase, lnblk

! write (6, *) 'polyroot: enter, string ='
! write (6, '(a)') string(1:76)

lstr = lnblk (string)

!  Look for starting value braces.

k = index (string, '{')
if (k == 0) goto 300
call findbraces (k, string, k1, k2)
if (nerror > 0) goto 300
if (k1 == 0) goto 200

!  Find comma before braces.

do i = k - 1, 1, -1
  if (string(i:i) /= ' ') then
    if (string(i:i) == ',') then
      goto 100
    else
      goto 200
    endif
  endif
enddo

goto 200

100 continue

stx1 = string(1:i-1)
lx1 = lnblk (stx1)
call parse (stx1, ntmp1, tmp1)
if (nerror > 0) goto 300

! write (6, *) 'polyroot: coefficients ='
! write (6, '(1p,4d19.11)') (dble (tmp1(i)), i = 1, ntmp1)

!   Find starting value(s).

stx2 = string(k1+1:k2-1)
lx2 = lnblk (stx2)
call parse (stx2, ntmp2, tmp2)
if (nerror > 0) goto 300
if (ntmp2 == 0 .or. ntmp2 > 2) goto 200

! write (6, *) 'polyroot: starting value(s) ='
! write (6, '(1p,2d25.15)') (dble (tmp2(i)), i = 1, ntmp2)

!   Find root, real or complex according to ntmp2 = 1 or 2.

if (ntmp2 == 1) then
  call rroot (ntmp1 - 1, tmp2(1), tmp1, nr, value)
  if (nerror > 0) goto 250
  nvalue = 1
else
  call croot (ntmp1 - 1, tmp2(1), tmp1, nr, value)
  if (nerror > 0) goto 250
  nvalue = 2
endif

goto 300

200 continue

write (6, 5)
5 format ('polyroot: syntax error in arguments to polyroot.'/ &
 'Examples: "polyroot[1, -1, -1, {0.618}]" (real root of 1-x-x^2=0 near 0.618)'/&
 '"polyroot[1, 1, 1, {-0.5, 0.866}]" (complex root of 1+x+x^2 near .5+.866i)')
nerror = 91
goto 300

250 continue

write (6, 6)
6 format ('polyroot: root not found')
nerror = 92

300 continue

! write (6, *) 'polyroot: exit'

return
end

subroutine qpslq (n, x, str, r)

!   This handles calls to pslq.

use mpmodule
use globdata
implicit none
integer i, i1, i2, j, k, ks, k1, k2, k3, lnblk, lstr, n, nsq, iq, lnamx(ntmpx), &
  nstring1
parameter (nstring1 = 1024)
character*1 string1(nstring1)
character*60 nam1, nam2, namx(ntmpx)
character*2048 str, str1
type (mp_real) t300, x(n), r(n)
real*8 d1, d2, dplog10q, gam, rb
parameter (nsq = 8, gam = 1.1547005438d0)
integer is0(3*n)
real*8 s1(7*n*n+nsq*n+3*n)
type (mp_real) s2(3*n*n+nsq*n+2*n)
external dplog10q, lnblk

t300 = 1.d300
rb = npslqb
if (npslql == 1) then
  call pslqm1 (ndebug, gam, n, nsq, rb, x, is0, s1, s2, iq, r)
elseif (npslql == 2) then
  call pslqm2 (ndebug, gam, n, nsq, rb, x, is0, s1, s2, iq, r)
elseif (npslql == 3) then
  write (6, *) 'qpslq: 3-level PSLQ not available yet.'
  iq = 0
endif
if (iq == 0) then
  write (6, 1) 
1 format ('evalfun: integer relation not found.')
  nerror = 101
elseif (ndebug > 0) then
  write (6, 2)
2 format ('Relation:  0 =')
  k = 0
  ks = 0
  lstr = lnblk (str)

!   Parse the string to identify each term.

  do j = 1, n

100 continue

    k = k + 1
    if (k > lstr) then
      goto 110
    elseif (str(k:k) == ',') then
      goto 110
    elseif (str(k:k) == '(') then
      call findpars (k, str, k1, k2)
      k = k2
      goto 100
    elseif (str(k:k) == '[') then
      call findbrackets (k, str, k1, k2)
      k = k2
      goto 100
    elseif (str(k:k) == '{') then
      call findbraces (k, str, k1, k2)
      k = k2
      goto 100
    else
      goto 100
    endif

110 continue
    if (k - 1 < ks + 1) goto 120
    nam1 = str(ks+1:k-1)
    k1 = k - ks - 1
    ks = k
    nam2 = ' '
    k2 = 0

    do i = 1, k1
      if (nam1(i:i) /= ' ') then
        k2 = k2 + 1
        nam2(k2:k2) = nam1(i:i)
      endif
    enddo

    lnamx(j) = k2
    namx(j) = nam2
  enddo

  if (k <= lstr) goto 120
  goto 200

120 continue

!   Parsing failed -- construct a series of simple names, e.g. pslqnnn.

  do j = 1, n
    lnamx(j) = 7
    write (nam1, '(i3.3)') j
    namx(j) = 'pslq'//nam1(1:3)
  enddo

200 continue

!   Find the length, in digits, of the largest r(j).

  d1 = 0.d0

  do j = 1, n
    if (r(j) == 0.d0) then
      d2 = 1.d0
    else
      d2 = 1.d0 + abs (dplog10q (r(j)))
    endif
    d1 = max (d1, d2)
  enddo

  k1 = min (int (d1 + 3.d0), nstring1)

!   Convert each r(j) to decimal (Fortran F format), then append * and name.
!   Output resulting string in 80-character blocks (almost always it will be
!   less than 80 characters total).

  do j = 1, n
    call mpfform (r(j), k1, 0, string1)
    k2 = lnamx(j)
    str1 = '+'

    do i = 1, k1
      str1(i+1:i+1) = string1(i)
    enddo

    str1(k1+2:k1+3) = '* '
    str1(k1+4:k1+k2+3) = namx(j)(1:k2)
    
    do i1 = 1, k1 + k2 + 3, 80
      i2 = min (i1 + 79, k1 + k2 + 3)
      write (6, '(80a1)') (str1(i:i), i = i1, i2)
    enddo
  enddo
endif
return
end

recursive function qinteg (string, tmp1)

!   This handles calls to integrate (numerical quadrature).

use mpmodule
use globdata
implicit none
integer i, i1, i2, i3, j, k, kvar, k1, k2, lnam1, lnblk, lstr, lr1, lr2, &
  lx1, lx2, lx3, ntmp1
character*2048 string, stx1, stx2, stx3
character*16 lcase, nam1, str1, str2
character*27 alphal, alphau
character*10 digits
parameter (digits = '0123456789', &
  alphal = 'abcdefghijklmnopqrstuvwxyz_', &
  alphau = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ_')
type (mp_real) qinteg, quad, tmp1(ntmpx), t1, t2, t3, x1, x2
type (mp_real) t300
external lcase, lnblk, quad

t300 = 1.d300

! write (6, *) 'qinteg: enter, string ='
! write (6, '(a)') string(1:76)

qinteg = 0.d0
k = 0
lstr = lnblk (string)

!   Look for function definition.

100 continue

k = k + 1
if (k > lstr) then
  goto 200
elseif (string(k:k) == '[') then
  call findbrackets (k, string, k1, k2)
  if (nerror > 0) goto 300
  k = k2
  goto 100
elseif (string(k:k) == '(') then
  call findpars (k, string, k1, k2)
  if (nerror > 0) goto 300
  k = k2
  goto 100
elseif (string(k:k) == ',') then
  goto 110
else
  goto 100
endif

110 continue

stx1 = string(1:k-1)
lx1 = lnblk (stx1)

! write (6, *) 'qinteg: stx1 ='
! write (6, '(a)') stx1(1:lx1)

!  Look for variable of integration and limits of integration.

k = k + 1
call findbraces (k, string, k1, k2)
if (nerror > 0) goto 300
if (k1 == 0) goto 200
k = k1 + 1

do i = k, lstr
  if (string(i:i) /= ' ') goto 120
enddo

120 continue

k = i
k1 = k + index (string(k+1:k2-1), ',')
if (k1 <= k) goto 200
lnam1 = lnblk (string(k:k1-1))
if (lnam1 > 16) goto 200
nam1 = string(k:k1-1)

do i = 1, lnam1
  i1 = index (alphal, nam1(i:i))
  i2 = index (alphau, nam1(i:i))
  i3 = index (digits, nam1(i:i))
  if (i == 1 .and. i1 + i2 == 0 .or. i1 + i2 + i3 == 0) goto 200
enddo

do i = 1, nvar
  if (lcase (nam1) == lcase (varn(i))) goto 130
enddo

if (nvar + 1 > nvarx) then
  write (6, 1) nvarz
1 format ('qinteg: too many user variables; max =',i5)
  nerror = 111
  goto 300
endif
nvar = nvar + 1
i = nvar
varn(i) = nam1

130 continue

kvar = i
stx2 = string(k1+1:k2-1)
lx2 = lnblk (stx2)
 
! write (6, *) 'qinteg: before parse; stx2 =', stx2(1:lx2)

call mpsetprecwords (nwords2)
call parse (stx2, ntmp1, tmp1)
if (nerror > 0) goto 300
if (ntmp1 /= 2) goto 200
x1 = tmp1(1)
x2 = tmp1(2)
call mpsetprec (nwords1)

! write (6, *) 'x1, x2 =', dble(min(max(x1,-t300),t300)), &
!   dble(min(max(x2,-t300),t300))

!  Perform quadrature.

if (x2 <= x1) then
  qinteg = 0.d0
elseif (x1 > - mplrg .and. x2 < mplrg) then

!  Neither limit is infinite.

  t1 = quad (kvar, stx1, x1, x2, tmp1)
  if (nerror > 0) goto 300
  qinteg = t1
else

!   Handle cases when one or both of the limits is infinite.

  str1 = varn(kvar)
  lr1 = lnblk (str1)
  lr2 = lnblk (varn(kvar))
  str2 = '(1/' // varn(kvar)(1:lr2) // ')'
  lr2 = lnblk (str2)
  call replace (str1, str2, stx1, stx2)
  lx2 = lnblk (stx2)
  stx3 = '(' // stx2(1:lx2) // ')/' // str1(1:lr1) // '^2'
  lx3 = lnblk (stx3)

  if (x1 > 0.d0 .and. x2 == mplrg) then
    x2 = 1.d0 / x1
    x1 = 0.d0
    t1 = quad (kvar, stx3, x1, x2, tmp1)
    if (nerror > 0) goto 300
    qinteg = t1
  elseif (x1 <= 0.d0 .and. x1 > -mplrg .and. x2 == mplrg) then
    x2 = 1.d0
    t1 = quad (kvar, stx1, x1, x2, tmp1)
    if (nerror > 0) goto 300
    x1 = 0.d0
    x2 = 1.d0
    t2 = quad (kvar, stx3, x1, x2, tmp1)
    if (nerror > 0) goto 300
    qinteg = t1 + t2
  elseif (x1 == -mplrg .and. x2 < 0.d0) then
    x1 = 1.d0 / x2
    x2 = 0.d0
    t1 = quad (kvar, stx3, x1, x2, tmp1)
    if (nerror > 0) goto 300
    qinteg = t1
  elseif (x1 == -mplrg .and. x2 >= 0.d0 .and. x2 < mplrg) then
    x1 = -1.d0
    t1 = quad (kvar, stx1, x1, x2, tmp1)
    if (nerror > 0) goto 300
    x2 = -1.d0
    x1 = 0.d0
    t2 = quad (kvar, stx3, x1, x2, tmp1)
    if (nerror > 0) goto 300
    qinteg = t1 + t2
  elseif (x1 == -mplrg .and. x2 == mplrg) then
    x1 = -1.d0
    x2 = 1.d0
    t1 = quad (kvar, stx1, x1, x2, tmp1)
    if (nerror > 0) goto 300
    x1 = -1.d0
    x2 = 0.d0
    t2 = quad (kvar, stx3, x1, x2, tmp1)
    if (nerror > 0) goto 300
    x1 = 0.d0
    x2 = 1.d0
    t3 = quad (kvar, stx3, x1, x2, tmp1)
    if (nerror > 0) goto 300
    qinteg = t1 + t2 + t3
  endif
endif
goto 300
  
200 continue

write (6, 2)
2 format ('qinteg: syntax error in arguments to integrate.'/ &
    'Example of proper call: "integrate[cos[x]*exp[-x], {x, 0, 1}]"')
nerror = 112

300 continue

return
end

function quad (kvar, stx1, x1, x2, tmp1)
use mpmodule
use globdata
implicit none
integer kvar
character*2048 stx1
type (mp_real) quad, quadgs, quaderf, quadts, tmp1(ntmpx), x1, x2
external quadgs, quaderf, quadts

if (nquadt == 1) then
  quad = quadgs (kvar, stx1, x1, x2, tmp1, nquadl, nquadz, quadwk, quadxk)
elseif (nquadt == 2) then
  quad = quaderf (kvar, stx1, x1, x2, tmp1, nquadl, nquadz, quadwk, quadxk)
elseif (nquadt == 3) then
  quad = quadts (kvar, stx1, x1, x2, tmp1, nquadl, nquadz, quadwk, quadxk)
else
  write (6, *) 'quad: impossible case; contact author'
  stop
endif
return
end

recursive function qsum (string, tmp1)

!   This handles calls to sum (numerical summation).

use mpmodule
use globdata
implicit none
integer i, i1, i2, i3, j, k, kvar, k1, k2, lnam1, lnblk, lstr, lx1, lx2, &
  ierror, ntmp1
character*2048 string, stx1, stx2
character*16 lcase, nam1
character*27 alphal, alphau
character*10 digits
parameter (digits = '0123456789', &
  alphal = 'abcdefghijklmnopqrstuvwxyz_', &
  alphau = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ_')
real*8 d1, d2
type (mp_real) eps, qsum, tmp1(ntmpx), t1, t2, t3
external lcase, lnblk

! write (6, *) 'qsum: enter, string ='
! write (6, '(a)') string(1:76)

qsum = 0.d0
k = 0
lstr = lnblk (string)

!   Look for function definition.

100 continue

k = k + 1
if (k > lstr) then
  goto 200
elseif (string(k:k) == '[') then
  call findbrackets (k, string, k1, k2)
  if (nerror > 0) goto 300
  k = k2
  goto 100
elseif (string(k:k) == '(') then
  call findpars (k, string, k1, k2)
  if (nerror > 0) goto 300
  k = k2
  goto 100
elseif (string(k:k) == ',') then
  goto 110
else
  goto 100
endif

110 continue

stx1 = string(1:k-1)
lx1 = lnblk (stx1)

! write (6, *) 'qsum: stx1 ='
! write (6, '(a)') stx1(1:lx1)

!  Look for index of summation and limits of summation.

k = k + 1
call findbraces (k, string, k1, k2)
if (nerror > 0) goto 300
if (k1 == 0) goto 200
k = k1 + 1

do i = k, lstr
  if (string(i:i) /= ' ') goto 120
enddo

goto 200

120 continue

k = i
k1 = k + index (string(k+1:k2-1), ',')
if (k1 <= k) goto 200
lnam1 = lnblk (string(k:k1-1))
if (lnam1 > 16) goto 200
nam1 = string(k:k1-1)

do i = 1, lnam1
  i1 = index (alphal, nam1(i:i))
  i2 = index (alphau, nam1(i:i))
  i3 = index (digits, nam1(i:i))
  if (i == 1 .and. i1 + i2 == 0 .or. i1 + i2 + i3 == 0) goto 200
enddo

do i = 1, nvar
  if (lcase (nam1) == lcase (varn(i))) goto 130
enddo

if (nvar + 1 > nvarx) then
  write (6, 1) nvarz
1 format ('qsum: too many user variables; max =',i5)
  nerror = 121
  goto 300
endif
nvar = nvar + 1
i = nvar
varn(i) = nam1

130 continue

kvar = i
stx2 = string(k1+1:k2-1)
lx2 = lnblk (stx2)

! write (6, *) 'qsum: before parse; stx2 =', stx2(1:lx2)

call parse (stx2, ntmp1, tmp1)

! write (6, *) 'qsum: after parse'

if (nerror > 0) goto 300
if (ntmp1 /= 2) goto 200

d1 = min (max (tmp1(1), mpreal (-1.d300)), mpreal (1.d300))
d2 = min (max (tmp1(2), mpreal (-1.d300)), mpreal (1.d300))

if (abs (d1) >= 1.d9 .or. abs (d2) >= 1.d9 .and. d2 < 1.d300 .or. &
  d1 /= anint (d1) .or. d2 /= anint (d2) .and. d2 < 1.d300) then
  write (6, 2) d1, d2
2 format ('qsum: improper limits of sum; values =',1p,2d25.15)
  nerror = 122
  goto 300
endif

k1 = anint (d1)
k2 = anint (min (d2, 1.d9))
t1 = 0.d0
j = 0
eps = mpreal (10.d0) ** nepsilon1

! write (6, *) 'k1, k2 =', k1, k2

!   Perform summation.

do k = k1, k2

!  write (6, *) 'qsum: before call to parse; k =', k

  var(kvar) = k
  call parse (stx1, ntmp1, tmp1(1))
  call mpgetpar ('mpier', ierror)
  if (ierror > 0 .or. nerror > 0) goto 300
  if (ntmp1 > 1) then
    write (6, 3)
3   format ('qsum: function returns more than one value per call.')
    nerror = 123
    goto 300
  endif
  t1 = t1 + tmp1(1)
  
!  write (6, *) 'qsum: after call to parse; tmp1 =', dble (tmp1(1))

!   For infinite sums, check for convergence (ten consecutive terms < eps).

  if (k2 == 1000000000) then
    if (abs (tmp1(1)) < eps) then
      j = j + 1
      if (j == 10) goto 140
    else
      j = 0
    endif
  endif
enddo

140 continue

qsum = t1
goto 300
  
200 continue

write (6, 5)
5 format ('qsum: syntax error in arguments to sum.'/ &
    'Example of proper call: "sum[n/(1+2^n), {n, 0, Infinity}]"')
nerror = 124

300 continue

! write (6, *) 'qsum: exit'

return
end

recursive subroutine table (string, tmp1, nvalue, value)

!   This handles calls to table, which forms lists.

use mpmodule
use globdata
implicit none
integer i, i1, i2, i3, j, k, kvar, k1, k2, lnam1, lnblk, lstr, lx1, lx2, &
  ierror, ntmp1, nvalue
character*2048 string, stx1, stx2
character*16 lcase, nam1
character*27 alphal, alphau
character*10 digits
parameter (digits = '0123456789', &
  alphal = 'abcdefghijklmnopqrstuvwxyz_', &
  alphau = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ_')
real*8 d1, d2
type (mp_real) tmp1(ntmpx), t1, value(ntmpx)
external lcase, lnblk

! write (6, *) 'table: enter, string ='
! write (6, '(a)') string(1:76)

k = 0
lstr = lnblk (string)

!   Look for function definition.

100 continue

k = k + 1
if (k > lstr) then
  goto 200
elseif (string(k:k) == '[') then
  call findbrackets (k, string, k1, k2)
  if (nerror > 0) goto 300
  k = k2
  goto 100
elseif (string(k:k) == '(') then
  call findpars (k, string, k1, k2)
  if (nerror > 0) goto 300
  k = k2
  goto 100
elseif (string(k:k) == ',') then
  goto 110
else
  goto 100
endif

110 continue

stx1 = string(1:k-1)
lx1 = lnblk (stx1)

! write (6, *) 'table: stx1 ='
! write (6, '(a)') stx1(1:lx1)

!  Look for index of table and limits of table.

k = k + 1
call findbraces (k, string, k1, k2)
if (nerror > 0) goto 300
if (k1 == 0) goto 200
k = k1 + 1

do i = k, lstr
  if (string(i:i) /= ' ') goto 120
enddo

goto 200

120 continue

k = i
k1 = k + index (string(k+1:k2-1), ',')
if (k1 <= k) goto 200
lnam1 = lnblk (string(k:k1-1))
if (lnam1 > 16) goto 200
nam1 = string(k:k1-1)

do i = 1, lnam1
  i1 = index (alphal, nam1(i:i))
  i2 = index (alphau, nam1(i:i))
  i3 = index (digits, nam1(i:i))
  if (i == 1 .and. i1 + i2 == 0 .or. i1 + i2 + i3 == 0) goto 200
enddo

do i = 1, nvar
  if (lcase (nam1) == lcase (varn(i))) goto 130
enddo

if (nvar + 1 > nvarx) then
  write (6, 1) nvarz
1 format ('table: too many user variables; max =',i5)
  nerror = 131
  goto 300
endif
nvar = nvar + 1
i = nvar
varn(i) = nam1

130 continue

kvar = i
stx2 = string(k1+1:k2-1)
lx2 = lnblk (stx2)

! write (6, *) 'table: before parse; stx2 =', stx2(1:lx2)

call parse (stx2, ntmp1, tmp1)

! write (6, *) 'table: after parse'

if (nerror > 0) goto 300
if (ntmp1 /= 2) goto 200

d1 = min (max (tmp1(1), mpreal (-1.d300)), mpreal (1.d300))
d2 = min (max (tmp1(2), mpreal (-1.d300)), mpreal (1.d300))

if (abs (d1) > 1.d9 .or. abs (d2) > 1.d9) then
  write (6, 2) d1, d2
2  format ('table: index limit is too large; values =',1p,2d25.15)
  nerror = 132
  goto 300
endif

k1 = anint (d1)
k2 = anint (d2)
if (k2 - k1 + 1 > ntmpx) then
  write (6, 3) ntmpx
3 format ('table: too many elements in list; max =',i4)
  nerror = 133
  goto 300
endif
i = 0

!   Form list.

do k = k1, k2

!    write (6, *) 'table: before call to parse; k =', k

  var(kvar) = k
  call parse (stx1, ntmp1, tmp1(1))
  call mpgetpar ('mpier', ierror)
  if (ierror > 0 .or. nerror > 0) goto 300
  if (ntmp1 > 1) then
    write (6, 4)
4   format ('table: function returns more than one value per call.')
    nerror = 134
    goto 300
  endif
  
!    write (6, *) 'table: after call to parse; t1 =', dble (t1)

  i = i + 1
  value(i) = tmp1(1)
enddo

nvalue = i

goto 300
  
200 continue

write (6, 5)
5 format ('table: syntax error in arguments to table.'/ &
    'Example of proper call: "table[alpha^n, {n, 0, 10}]"')
nerror = 135

300 continue

! write (6, *) 'table: exit'

return
end

function lnblk (string)

!   This finds the index of the last non-blank character in string.

integer i, lnblk
character*(*) string

do i = len (string), 1, -1
  if (string(i:i) .ne. ' ') goto 110
enddo

i = 0

110  continue

lnblk = i
return
end

function lcase (string)

!   This converts string to lower case.

integer i, k
character*16 string, lcase
character*27 alphal, alphau
parameter (alphal = 'abcdefghijklmnopqrstuvwxyz%', &
  alphau = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ%')

lcase = string

do i = 1, 16
  k = index (alphau, string(i:i))
  if (k > 0) lcase(i:i) = alphal(k:k)
enddo

return
end

subroutine findbraces (k, string, k1, k2)

!   This routine finds a set of matching braces beginning at position k
!   in string.  The positions are returned in k1 and k2.

use globdata
implicit none
integer i1, i2, j, k, k1, k2, k3, lnblk, lstr
character*2048 string
external lnblk

i1 = 0
i2 = 0
k1 = 0
k2 = 0
lstr = lnblk (string)

do j = k, lstr
  if (string(j:j) == '{') then
    i1 = i1 + 1
    i2 = 1
    if (i1 == 1) k1 = j
  elseif (string(j:j) == '}') then
    i1 = i1 - 1
    i2 = 1
    if (i1 == 0) then
      k2 = j
      goto 200
    elseif (i1 < 0) then
      write (6, 1) 
1     format ('findbraces: mismatched right brace.')
      nerror = 141
      goto 200
    endif
  elseif (i2 == 0 .and. string(j:j) /= ' ') then
    goto 200
  endif
enddo

if (i1 > 0) then
  write (6, 2) 
2 format ('findbraces: mismatched left brace.')
  nerror = 142
endif

200 continue

return
end

subroutine findbrackets (k, string, k1, k2)

!   This routine finds a set of matching brackets beginning at position k
!   in string.  The positions are returned in k1 and k2.

use globdata
implicit none
integer i1, i2, j, k, k1, k2, k3, lnblk, lstr
character*2048 string
external lnblk

i1 = 0
i2 = 0
k1 = 0
k2 = 0
lstr = lnblk (string)

do j = k, lstr
  if (string(j:j) == '[') then
    i1 = i1 + 1
    i2 = 1
    if (i1 == 1) k1 = j
  elseif (string(j:j) == ']') then
    i1 = i1 - 1
    i2 = 1
    if (i1 == 0) then
      k2 = j
      goto 200
    elseif (i1 < 0) then
      write (6, 1) 
1     format ('findbrackets: mismatched right bracket.')
      nerror = 151
      goto 200
    endif
  elseif (i2 == 0 .and. string(j:j) /= ' ') then
    goto 200
  endif
enddo

if (i1 > 0) then
  write (6, 2) 
2 format ('findbrackets: mismatched left bracket.')
  nerror = 152
endif

200 continue

return
end

subroutine findpars (k, string, k1, k2)

!   This routine finds a set of matching parentheses beginning at position k
!   in string.  The positions are returned in k1 and k2.

use globdata
implicit none
integer i1, i2, j, k, k1, k2, k3, lnblk, lstr
character*2048 string
external lnblk

i1 = 0
i2 = 0
k1 = 0
k2 = 0
lstr = lnblk (string)

do j = k, lstr
  if (string(j:j) == '(') then
    i1 = i1 + 1
    i2 = 1
    if (i1 == 1) k1 = j
  elseif (string(j:j) == ')') then
    i1 = i1 - 1
    i2 = 1
    if (i1 == 0) then
      k2 = j
      goto 200
    elseif (i1 < 0) then
      write (6, 1) 
1     format ('findpars: mismatched right parenthesis.')
      nerror = 161
      goto 200
    endif
  elseif (i2 == 0 .and. string(j:j) /= ' ') then
    goto 200
  endif
enddo

if (i1 > 0) then
  write (6, 2) 
2 format ('findpars: mismatched left parenthesis.')
  nerror = 162
endif

200 continue

return
end

subroutine replace (str1, str2, stx1, stx2)

!   This routine replaces string str1 with str2 in stx1, returning stx2,
!   provided str1 is delimited by one of the characters in delim.

implicit none
integer i, j, k, koff, lnblk, lr1, lr2, lx1, lx2
character*16 str1, str2
character*2048 stx1, stx2
character*27 alphal, alphau
character*11 digits
character*13 delim
parameter ( &
  alphal = 'abcdefghijklmnopqrstuvwxyz%', &
  alphau = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ%', &
  digits = '0123456789.', delim = ' +-*/^:=,()[]')
external lnblk

lr1 = lnblk (str1)
lr2 = lnblk (str2)
lx1 = lnblk (stx1)
lx2 = lx1
stx2 = stx1
k = 0
koff = 0

100 continue

if (k > lx1) goto 200
i = index (stx1(k:lx1), str1(1:lr1))
if (i == 0) goto 200
j = i + k - 1

if (j == 1) then
  if (j + lr1 - 1 == lx1 .or. index (delim, stx1(j+lr1:j+lr1)) > 0) then
    goto 110
  endif
elseif (index (delim, stx1(j-1:j-1)) > 0 .and. (j + lr1 - 1 == lx1 .or. &
  index (delim, stx1(j+lr1:j+lr1)) > 0)) then
  goto 110
endif
k = j + lr1
goto 100

110 continue

stx2 = stx2(1:j-1+koff) // str2(1:lr2) // stx2(j+lr1+koff:lx2)
koff = koff + (lr2 - lr1)
lx2 = lx2 + (lr2 - lr1)
k = j + lr1
goto 100

200 continue

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
