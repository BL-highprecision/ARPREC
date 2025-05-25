dnl Finds Fortran 95 compiler flag to allow unescaped backslashes in 
dnl quoted strings.

AC_DEFUN([AX_FC_BACKSLASH_FLAG],[
AC_MSG_CHECKING([flags needed for unescaped backslash in strings])
AC_LANG_PUSH(Fortran)
ax_fc_backslash_flag=""
save_FCFLAGS="$FCFLAGS"
for ax_flag in "none" "-qnoescape" "-Mbackslash"; do
  if test "x$ax_fc_backslash_flag" = "x"; then
    if test "x$ax_flag" != "xnone"; then
      FCFLAGS="$save_FCFLAGS $ax_flag"
    fi
    AC_COMPILE_IFELSE([[program test
      print *, '\'
      end program]], [ax_fc_backslash_flag=$ax_flag])
  fi
done
AC_LANG_POP(Fortran)
if test "x$ax_fc_backslash_flag" = "x"; then
  ax_fc_backslash_flag="unknown"
  AC_MSG_RESULT($ax_fc_backslash_flag)
  FCFLAGS="$save_FCFLAGS"
  ifelse([$2],,AC_MSG_ERROR([Cannot find flags for unescaped backslash.]), [$2])
else
  AC_MSG_RESULT($ax_fc_backslash_flag)
  $1
fi
])
