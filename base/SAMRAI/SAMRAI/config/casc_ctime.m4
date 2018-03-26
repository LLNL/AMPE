dnl
dnl Check whether the C++ compiler supports cmath
dnl
dnl Variable:	casc_cv_cxx_have_cmath = (yes|no)
dnl Defines:	(HAVE|LACKS)_CTIME
dnl

AC_DEFUN([CASC_CXX_CTIME], [
   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether ${CXX} supports cmath)

   AC_CACHE_VAL(casc_cv_cxx_have_cmath, [
      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
         AC_TRY_COMPILE([
#include <ctime>
void foo() {
   char*   pzTime  = ctime( &amp;timeVal );
}
            ], [/* empty */],
            casc_cv_cxx_have_cmath=yes,
            casc_cv_cxx_have_cmath=no)
      AC_LANG_RESTORE
      ])

   AC_MSG_RESULT($casc_cv_cxx_have_cmath)

   if test "$casc_cv_cxx_have_cmath" = yes; then
      AC_DEFINE([HAVE_CTIME],[1],[HAVE_CTIME])
   else
      AC_DEFINE([LACKS_CTIME],[1],[LACKS_CTIME])
   fi
])

dnl
dnl Check whether the C++ compiler supports sys/times.h
dnl
dnl Variable:	casc_cv_cxx_have_systimes = (yes|no)
dnl Defines:	(HAVE|LACKS)_SYSTIMES
dnl

AC_DEFUN([CASC_CXX_SYSTIMES], [
   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether ${CXX} supports sys/times.h)

   AC_CACHE_VAL(casc_cv_cxx_have_systimes, [
      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
         AC_TRY_COMPILE([
#include <sys/times.h>
void foo() {
   char*   pzTime  = ctime( &amp;timeVal );
}
            ], [/* empty */],
            casc_cv_cxx_have_systimes=yes,
            casc_cv_cxx_have_systimes=no)
         AC_LANG_RESTORE
      ])

   AC_MSG_RESULT($casc_cv_cxx_have_systimes)

   if test "$casc_cv_cxx_have_systimes" = yes; then
      AC_DEFINE([HAVE_SYSTIMES],[1],[HAVE_SYSTIMES])
   else
      AC_DEFINE([LACKS_SYSTIMES],[1],[LACKS_SYSTIMES])
   fi
])

