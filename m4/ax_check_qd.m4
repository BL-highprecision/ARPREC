dnl AC_CHECK_QD(version, action-if, action-not)
dnl Copyright © 2007 Yozo Hida <yozo@cs.berkeley.edu>
dnl
dnl File template taken from ac_check_curl.m4 in autoconf macro archive:
dnl Copyright © 2005 Akos Maroy <darkeye@tyrell.hu>
dnl
dnl Copying and distribution of this file, with or without modification, 
dnl are permitted in any medium without royalty provided the copyright 
dnl notice and this notice are preserved.
dnl
AC_DEFUN([AC_CHECK_QD], [
  succeeded=no

  if test -z "$QD_CONFIG"; then
    AC_PATH_PROG(QD_CONFIG, qd-config, no)
  fi

  if test "$QD_CONFIG" = "no" ; then
    echo "*** The qd-config script could not be found. Make sure it is"
    echo "*** in your path, and that the QD library is properly installed."
  else
    dnl drop the patchlevel (may not be a number).
    QD_VERSION=`$QD_CONFIG --version | cut -d. -f1,2`
    AC_MSG_CHECKING(for qd library version >= $1)
        VERSION_CHECK=`expr $QD_VERSION \>\= $1`
        if test "$VERSION_CHECK" = "1" ; then
            AC_MSG_RESULT(yes)
            succeeded=yes

            AC_MSG_CHECKING(for qd library C++ flags)
            QD_CXXFLAGS=`$QD_CONFIG --cxxflags`
            AC_MSG_RESULT($QD_CXXFLAGS)
        else
            QD_CXXFLAGS=""
            ## If we have a custom action on failure, don't print errors, but
            ## do set a variable so people can do so.
            ifelse([$3], ,echo "can't find qd library version >= $1",)
        fi
  fi

  if test $succeeded = yes; then
     ifelse([$2], , :, [$2])
  else
     ifelse([$3], , AC_MSG_ERROR([Library requirements (qd) not met.]), [$3])
  fi
])
