#  configure.ac
#
#  This work was supported by the Director, Office of Science, Division
#  of Mathematical, Information, and Computational Sciences of the
#  U.S. Department of Energy under contract number DE-AC03-76SF00098.
# 
#  Copyright (c) 2001-2002
#  
#  Creation Date:  May 1, 2002
#  Modified:       
#  Version:        1.0
#
#
# Process this file with autoconf to produce a configure script.
define([ARPREC_PATCHLEVEL], 1)
AC_PREREQ(2.60)
AC_INIT(arprec, 2.2.ARPREC_PATCHLEVEL, dhbailey@lbl.gov)
AC_CONFIG_SRCDIR([src/init.cpp])

PACKAGE="arprec"
MAJOR_VERSION="2"
MINOR_VERSION="8"
PATCH_VERSION="ARPREC_PATCHLEVEL"
VERSION=$MAJOR_VERSION.$MINOR_VERSION.$PATCH_VERSION
AC_DEFINE_UNQUOTED([MAJOR_VERSION],$MAJOR_VERSION,[arprec major version number])
AC_DEFINE_UNQUOTED([MINOR_VERSION],$MINOR_VERSION,[arprec minor version number])
AC_DEFINE_UNQUOTED([PATCH_VERSION],$PATCH_VERSION,[arprec patch number (sub minor version)])
echo "configuring $PACKAGE $VERSION..."

AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE([nostdinc])
AM_CONFIG_HEADER([config.h include/arprec/arprec_config.h])

AC_CANONICAL_HOST

# Checks for command-line options
# --enable-inline option.
AC_ARG_ENABLE(enable_inline, AC_HELP_STRING(--enable-inline, [inline commonly used functions.  [[default=yes]]]))
AC_MSG_CHECKING([if commonly used function is to be inlined])
if test "x$enable_inline" != "xno"; then
  enable_inline="yes"
  AC_DEFINE([ARPREC_INLINE], [1], [Define to 1 to inline commonly used functions.])
fi
AC_MSG_RESULT($enable_inline)

# --enable-debug
AC_ARG_ENABLE(enable_debug, AC_HELP_STRING(--enable-debug, [enable debugging code. [[default=no]]]))
AC_MSG_CHECKING([if debugging code is to be enabled])
if test "x$enable_debug" = "xyes"; then
  AC_DEFINE([ARPREC_DEBUG], [1], [Define to 1 to enable debugging code.])
else
  enable_debug="no"
fi
AC_MSG_RESULT($enable_debug)

# --enable-always-fix-fpu
AC_ARG_ENABLE(enable_always_fix_fpu, AC_HELP_STRING(--enable-always-fix-fpu, [Use FPU fix whenever C interface is called.  [[default=no]]]))
if test "x$enable_always_fix_fpu" = "xyes"; then
  AC_DEFINE([ALWAYS_FIX_FPU], [1], [Define to 1 to use FPU fix whenever C interface is called.])
else
  enable_always_fix_fpu="no"
fi

# Checks for programs.
AC_LANG(C++)
if test "x$CXXFLAGS" = "x"; then
  if test "x$enable_debug" = "xyes"; then
    CXXFLAGS="-O"
  else
    CXXFLAGS="-O2"
  fi
fi



# Set up compiler search list.  DUe to possible case insensitive filesystems, 
# Cygwin and Apple systems sometimes gets confused between cc and CC.
# Try to use vendor compilers first, as they are often the fastest.
case $host in 
  *cygwin* | *mingw* ) cxx_list="" ;;
  *apple* ) cxx_list="xlC_r xlC icpc pathCC c++ g++ gpp" ;;
  * ) cxx_list="xlC_r xlC icpc pathCC pgCC cxx aCC CC cc++ FCC KCC RCC c++ g++ gpp" ;;
esac

AC_PROG_CXX($cxx_list)
if test "x$CC" = "x"; then
  CC="$CXX"
fi
AC_PROG_CC

# Set CXX_PROG_NAME to the program name only.
CXX_PROG_NAME=`basename $CXX`

case "$CXX_PROG_NAME" in
  icpc|icc )
    CONFIG_CXXFLAGS="-mp"
    CXXFLAGS="$CXXFLAGS -mp -wd1572"
  ;;
  xlC* )
    CONFIG_CXXFLAGS="-qstrict"
    CXXFLAGS="$CXXFLAGS -qstrict"
  ;;
  CC )
    case "$host" in
      sparc-* )
        CXXFLAGS="$CXXFLAGS -dalign"
      ;;
    esac
  ;;
  * )
    if test "x$GXX" = "xyes"; then
      CXXFLAGS="$CXXFLAGS -Wall"
    fi
  ;;
esac
if test "x$enable_debug" = "xyes"; then
  CXXFLAGS="$CXXFLAGS -g"
fi

# --enable-qd
AC_ARG_ENABLE(enable-qd, AC_HELP_STRING(--enable-qd, [Use quad-double (QD)
library. [[default=no]]]))
if test "x$enable_qd" = "xyes"; then
  AC_CHECK_QD([2.3], [CXXFLAGS="$CXXFLAGS $QD_CXXFLAGS"])
  AC_DEFINE([ARPREC_QD], [1], [Define to 1 to use the QD library.])
else
  enable_qd="no"
fi

# --enable-fma
AC_ARG_ENABLE(enable_fma, AC_HELP_STRING(--enable-fma, [use fused multiply-add/subtract (auto,gnu,c99,ibm,compiler).  Use this option with care.  [[default=auto]]]))
if test "x$enable_fma" = "x"; then
  enable_fma="auto"
fi
case "$enable_fma" in
  yes)
    AX_CXX_FMA
    if test "x$ax_cxx_fma" = "x"; then
      AC_MSG_ERROR([Cannot find working fused multiply add.])
    fi
  ;;
  auto)
    case $host in
      powerpc*-*-*)
        AX_CXX_FMA(["ibm gnu"])
      ;;
    esac
  ;;
  no) 
    ax_cxx_fma =
  ;;
  *) 
    AX_CXX_FMA([$enable_fma])
    if test "x$ax_cxx_fma" = "x"; then
      AC_MSG_ERROR([Cannot find working fused multiply add.])
    fi
  ;;
esac
if test "x$ax_cxx_fma" = "x"; then
  ax_cxx_fma="none"
  ax_cxx_fms="none"
else
  AC_DEFINE_UNQUOTED([ARPREC_FMA(x,y,z)], [$ax_cxx_fma], 
    [If fused multiply-add is available, define correct macro for using it.])
  AC_DEFINE_UNQUOTED([ARPREC_FMS(x,y,z)], [$ax_cxx_fms], 
    [If fused multiply-subtract is available, define correct macro for using it.])
fi

# Check for F90 compiler
if test "x$FCFLAGS" = "x"; then
  if test "x$enable_debug" = "xyes"; then
    FCFLAGS="-O"
  else
    FCFLAGS="-O2"
  fi
fi
AC_ARG_ENABLE(enable_fortran, AC_HELP_STRING(--enable-fortran, [build Fortran interfaces [[default=auto]]]))
if test "x$enable_fortran" != "xno"; then
  AC_LANG_PUSH(Fortran)
  AC_PROG_FC([xlf95 ifort pathf95 f95 gfortran g95 pgf95 lf95 fort ifc efc pathf90 xlf90 pgf90 epcf90 xlf f90])
  if test "x$FC" != "x"; then
    enable_fortran="yes"
    AC_DEFINE([HAVE_FORTRAN], [1], [Define to 1 if Fortran interface is to be compiled.])

    AC_FC_LIBRARY_LDFLAGS
    AC_FC_FREEFORM
    AX_CXX_FC_LIB

    AC_FC_WRAPPERS
    AC_FC_MAIN
    
    fcname=`basename $FC`

    case $host in
      i?86-*-* | k?*-*-* | athlon-*-* | x86_64-*-*)
        case $fcname in
          ifc|ifort) 
            FCFLAGS="$FCFLAGS -mp" 
            CONFIG_FCFLAGS="-mp"
          ;;
        esac
      ;;
      *-*-aix* | powerpc-apple-*)
        case $fcname in
          xlf*) 
            FCFLAGS="$FCFLAGS -qstrict"
            CONFIG_FCFLAGS="-qstrict"
          ;;
        esac
      ;;
    esac

    AX_F90_MODULE_FLAG
    AX_FC_BACKSLASH_FLAG
    AC_SUBST(ax_f90_modflag)

    AX_F90_MODULE_STYLE
    module_case=`echo "$ax_f90_module_style" | sed 's/,.*$//'`
    module_ext=`echo "$ax_f90_module_style" | sed 's/^.*, //'`
    AC_SUBST(module_ext)

    AX_FC_ETIME
    AC_SUBST(ETIME, $ax_fc_etime)
  else
    enable_fortran="no"
  fi
  if test "x$enable_debug" = "xyes"; then
    FCFLAGS="$FCFLAGS -g"
  fi

  AC_LANG_POP(Fortran)
fi
AM_CONDITIONAL(HAVE_FORTRAN, [test "x$enable_fortran" = "xyes"])
AM_CONDITIONAL(UPCASE_MODULE, [test "$module_case" = "uppercase"])

# Check ranlib
AC_PROG_RANLIB

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([fp.h sys/time.h])
if test "x$enable_debug" = "xyes"; then
  AC_CHECK_HEADER([valgrind/memcheck.h], 
    [AC_DEFINE([ARPREC_VALGRIND], [1], [Define to 1 to use valgrind debugging code.  Requires ARPREC_DEBUG to be set.])])
fi
AC_CXX_TEMPLATES

# Checks for typedefs, structures, and compiler characteristics.

# Check for x86 FPU fix
AC_MSG_CHECKING([whether x86 FPU control code is needed])
x86_fix="no"
case "$host" in
  i?86-*-* | k?*-*-* | athlon-*-* | x86_64-*-*)
    AC_DEFINE(X86, [1], [Whether to use x86 fpu fix.])
    x86_fix="yes"
  ;;
esac
AC_MSG_RESULT($x86_fix)

# Check for fpu_control.h
if test "x$x86_fix" = "xyes"; then
  AC_CHECK_HEADERS(fpu_control.h)
fi

# Checks for library functions.
# drand48 function
AC_CHECK_FUNCS([drand48 copysign])

AC_DEFINE([ARPREC_API], [], [Any special symbols needed for exporting APIs.])

# Output
AC_SUBST(CONFIGURE_ARGS, $ac_configure_args)
AC_SUBST(CONFIG_CXXFLAGS, $CONFIG_CXXFLAGS)
AC_SUBST(CONFIG_FCFLAGS, $CONFIG_FCFLAGS)
AC_SUBST(F77, $FC)
AC_SUBST(FFLAGS, $FCFLAGS)

AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile include/Makefile
                 fortran/Makefile toolkit/Makefile arprec-config 
                 fortran/second.f toolkit/second.f])

AC_OUTPUT

# Inform user of important options
echo
echo "Using "
echo "   C++ Compiler   = $CXX"
echo "   C++ Flags      = $CXXFLAGS $DEFS"
if test "x$enable_fortran" = "xyes"; then
  echo "   F90 Compiler   = $FC"
  echo "   F90 Flags      = $FCFLAGS"
  echo "   F90 Libs       = $FCLIBS"
fi
echo "   enable fortran = $enable_fortran"
echo "   enable inline  = $enable_inline"
echo "   enable x86 fix = $x86_fix"
echo "   enable debug   = $enable_debug"
echo "   fma            = $ax_cxx_fma"
echo "   fms            = $ax_cxx_fms"
if test "x$ac_cv_cxx_templates" != "xyes"; then
  echo
  echo  "Warning: the C++ compiler does not seem to support templates."
  echo  "         Some C++ example programs may not compile correctly."
fi

