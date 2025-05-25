module globdata

!   This module defines global parameters and variables for the mathinit and
!   mathtool programs.

!   David H Bailey    2004-06-14

!     Parameter   Default
!     Name        Value  Explanation

!     mxres       999    Size of result table, ie when one types "result[34]".
!     nfuna        27    Number of predefined functions.
!     nfunx       100    Maximum number of functions (including predefined).
!     nfunz        75    nfunx - nfuna.
!     ndigi       100    Initial precision level, in digits, not exceed ndigmx1.
!     ndigmx1    1000    Maximum primary precision level, in digits (see note).
!     ndigmx2    2000    Maximum secondary precision level, in digits (see note).
!     nquadx       10    Maximum number of quadrature "levels".
!     nquadz    18432    Size of quadrature table, in MP values (see note).
!     ntmpx       100    Maximum number of temporary variables.
!     nvara         7    Number of predefined constants (e, pi, etc).
!     nvarb         9    Number of predefined arguments (arg1, arg2, etc).
!     nvarc         7    Number of predefined variables (digits, epsilon, etc).
!     nvarx       200    Maximum number of variables (including predefined).
!     nvarz        76    nvarx - nvara - nvarb - nvarc.

!     Variable
!     Name        Explanation

!     ndebug      Debug printout level: 0 = none; 3 = max.
!     ndigits1    Primary working precision, in digits.
!     ndigits2    Secondary precision, in digits.
!     nef         Output format: 1 = E format; 2 = F format.
!     nef1        Width of output format.
!     nef2        Digits after period in output format.
!     nerror      Error number (error has occurred if nonzero).
!     nepsilon1   Primary epsilon variable (epsilon value = 10^(10-nepsilon)).
!     nepsilon2   Secondary epsilon variable (epsilon2 value = 10^(10-nepsilon2).
!     nfun        Number of defined functions (including predefined).
!     npslqb      PSLQ bound (iterations are terminated when bound > 10^npslqb).
!     npslql      PSLQ level: 1, 2 or 3.  PSLQ level 3 is not yet implemented.
!     nquadl      Quadrature level (at most nquadl levels are used).
!     nquadt      Quadrature type: 1 = Gaussian; 2 = erf; 3 = tanh-sinh.
!     nvar        Number of defined variables (including predefined).
!     nwords1     Number of words of primary precision.
!     nwords2     Number of words of secondary precision.

!     Array
!     Name    Type     Explanation

!     narg    integer  Number of arguments for defined functions.
!     fund    char     Function definitions.
!     funn    char     Function names.
!     varn    char     Variable names.
!     quadn   char     Names of the three quadrature types.
!     quadwk  mp_real  Quadrature weights.
!     quadxk  mp_real  Quadrature abscissas.
!     result  mp_real  Array of results, ie when one types "result[34]".
!     var     mp_real  Variable values.

!   If any change is made to this file, then at the least this file must be
!   re-compiled, then all the mathtool files must be re-compiled and re-linked.
!   If either ndigmx1, ndigmx2 or nquadx is changed, then in addition the 
!   mathinit.f program must be re-compiled and re-run.  Additionally, if
!   ndigmx1 or ndigmx2 is changed, make sure that neither exceeds mpipl in
!   mp_mod.f.  If either exceeds, then mpipl must be increased to at least this
!   level and mp_mod.f must be re-compiled; then all mathtool files must be
!   re-compiled, and the mathinit.f program must be re-compiled and re-run.
!   Finally, if ndigmx1 exceeds 2000, then all instances of "2048" in this file,
!   mathtool.f and quadsub.f must be changed to at least ndigmx1 + 40, and then
!   all programs must be re-compiled and re-linked.

use mpmodule
implicit none
integer i, j
private i, j
integer mxres, ndebug, ndigi, ndigmx1, ndigmx2, ndigits1, ndigits2, nef, &
  nef1, nef2, nepsilon1, nepsilon2, nerror, nfun, nfuna, nfunx, nfunz, &
  npslqb, npslql, nquadl, nquadt, nquadx, nquadz, ntmpx, nvar, nvara, nvarb, &
  nvarc, nvarx, nvarz, nwords1, nwords2
parameter (mxres = 999, nfuna = 28, nfunx = 100, nfunz = nfunx - nfuna, &
  ndigi = 100, ndigmx1 = 1000, ndigmx2 = 2000, nquadx = 10, &
  nquadz = 18 * 2 ** nquadx + 200, ntmpx = 100, nvara = 7, nvarb = 9, &
  nvarc = 9, nvarx = 200, nvarz = nvarx - nvara - nvarb - nvarc)
integer narg(nfunx)
character*2048 fund(nfunx)
character*76 funhelp(4,nfuna)
character*16 funn(nfunx)
character*16 varn(nvarx)
character*16 quadn(3)
type (mp_real) quadwk(-1:nquadz), quadxk(-1:nquadz), result(mxres), var(nvarx)
data narg /1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 4, 1, 2, 2, 1000, 1000, &
  1, 1, 1, 4, 4, 1, 1, 1000, 1000, nfunz*0/
data fund /100 * ' '/
data funn /'Abs', 'Arccos', 'Arcsin', 'Arctan', 'Arctan2', 'Bessel', &
  'Besselexp', 'Binomial', 'Cos', 'Erf', 'Exp', 'Factorial', 'Gamma', &
  'Integrate', 'Log', 'Max', 'Min', 'Polyroot', 'Pslq', 'Result', 'Sin', &
  'Sqrt', 'Sum', 'Table', 'Tan', 'Zeta', 'Zetap', 'Zetaz', nfunz*' '/
data varn /'E', 'Log2', 'Log10', 'Pi', 'Catalan', 'Eulergamma', 'Infinity', &
  'Arg1', 'Arg2', 'Arg3', 'Arg4', 'Arg5', 'Arg6', 'Arg7', 'Arg8', 'Arg9', &
  'Debug', 'Digits', 'Digits2', 'Epsilon', 'Epsilon2', 'Pslqbound', &
  'Pslqlevel', 'Quadlevel', 'Quadtype', nvarz*' '/
data quadn /'Gaussian', 'Error function', 'Tanh-Sinh'/
data ((funhelp(i,j), i = 1, 4), j = 1, 14) / &
  'Abs[x] computes the absolute value of the x.', 3*' ', &
  'Arccos[x] computes arccos of x, with result in [0,pi].', 3*' ',&
  'Arcsin[x] computes arcsin of x, with result in [-pi/2,pi/2].', 3*' ', &
  'Arctan[x] computes arctan of x, with result in [-pi/2,pi/2].', 3*' ',&
  'Arctan2[y,x] computes arctangent of y/x, ie placing result in (-pi,pi]', &
  'according to the coordinates (x,y).  x or y may be 0 but not both.', 2*' ', &
  'Bessel[x] computes BesselI[0,x].', 3*' ', &
  'Besselexp[x] computes BesselI[0,x] / Exp[x].', 3*' ',&
  'Binomial[m,n] computes the binomial coefficient of integers (m,n).', 3*' ', &
  'Cos[x] computes the cosine of x.', 3*' ', &
  'Erf[x] computes the error function of x.', 3*' ', &
  'Exp[x] computes the exponential function of x.', 3*' ', &
  'Factorial[n] computes the factorial of the integer n.', 3*' ', &
  'Gamma[x] computes the Gamma function.', 3*' ', &
  'Integrate[fun[x], {x, a, b}] computes the integral of function fun[x] (which', &
  'may be an expression involving x), from x=a to x=b. Either a or b or both', &
  'may be "infinity" or "-infinity".  fun[x] may have a singularity at a or b,',&
  'but not between. quadlevel (4-10) controls level; quadtype (1-3) is type.'/
data ((funhelp(i,j), i = 1, 4), j = 15, 20) / &
  'Log[x] computes the natural logarithm of x.', 3*' ', &
  'Max[x,y] returns the maximum of x and y.', 3*' ', &
  'Min[x,y] returns the minimum of x and y.', 3*' ', &
  'Polyroot[1,2,3,4,{-6}] finds the real root of polynomial 1+2x+3x^2+4x^3', &
  'near -6. Polyroot[1, 2, 3, {-0.33, 0.47}] finds the complex root of', &
  '1+2x+3x^2 near -0.33+0.47i.  Starting value should be close to root,', &
  'otherwise no root may be found.  Some experimentation may be necessary.', &
  'Pslq[x1,x2,...,xn] finds integers ai such that a1*x1+a2*x2+...+an*xn = 0,', &
  'to within tolerance 10^epsilon.  The xi may be expressions; n <= 50. The', &
  'search is terminated when the relation bound exceeds 10^pslqbound. A two-', &
  'level PSLQ (faster for large n) may be selected by typing "pslqlevel=2".', &
  'Result[n] is used to recall the n-th result from the Toolkit.', 3*' '/
data ((funhelp(i,j), i = 1, 4), j = 21, nfuna) / &
  'Sin[x] computes the sine of x.', 3*' ', &
  'Sqrt[x] computes the square root of x.', 3*' ', &
  'Sum[fun[n], {n,a,b}] computes the sum of the function fun[n] (which may be', &
  'an expression involving n), from n=a to n=b (integer a and b). The', &
  'parameter b may be "infinity".  The function fun should rapidly go to 0;', &
  'otherwise "infinite" summations may take an unreasonably long run time.', &
  'Table[fun[n], {n,a,b}] generates the list (fun[a], fun[a+1], ..., fun[b]),', &
  'where a and b are integers.  No more than 50 entries may generated.', 2*' ', &
  'Tan[x] computes the tangent of x.', 3*' ', &
  'Zeta[n] computes the Riemann zeta function of positive integer n.', 3*' ', &
  'Zetap[n1, n2, ..., nk] computes the multi-zeta function zetap for positive', &
  'integers ni. See http://www.cecm.sfu.ca/projects/EZFace for definition.', &
  2*' ', &
  'Zetaz[n1, n2, ..., nk] computes the multi-zeta function zetaz for positive', &
  'integers ni. See http://www.cecm.sfu.ca/projects/EZFace for definition.', &
  2*' '/
end
