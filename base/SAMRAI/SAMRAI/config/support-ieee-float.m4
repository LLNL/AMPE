dnl $Id$

AC_DEFUN([CASC_C_IEEE_FLOAT],[
dnl Check on certain declarations in the float.h file:
dnl FLT_SNAN DBL_SNAN
dnl
dnl ac_define ..._IS_BROKEN for symbols that are not defined.
dnl
# Begin macro CASC_IEEE_FLOAT

AC_LANG_C

AC_EGREP_CPP([^nan is broken],
[#include <float.h>
#ifndef NAN
nan is broken
#endif],
AC_DEFINE([NAN_IS_BROKEN],[1],[Define if NAN is not in float.h])
CASC_AC_LOG(["NAN is broken (not in float.h)"]),
CASC_AC_LOG(["NAN is ok (in float.h)"])
)

AC_EGREP_CPP([^flt snan is broken],
[#include <float.h>
#ifndef FLT_SNAN
flt snan is broken
#endif],
AC_DEFINE([FLT_SNAN_IS_BROKEN],[1],[Define if FLT_SNAN is not in float.h])
CASC_AC_LOG(["FLT_NAN is broken (not in float.h)"]),
CASC_AC_LOG(["FLT_NAN is ok (in float.h)"])
)

AC_EGREP_CPP([^dbl snan is broken],
[#include <float.h>
#ifndef DBL_SNAN
dbl snan is broken
#endif],
AC_DEFINE([DBL_SNAN_IS_BROKEN],[1],[Define if DBL_SNAN is not in float.h])
CASC_AC_LOG(["DBL_NAN is broken (not in float.h)"]),
CASC_AC_LOG(["DBL_NAN is ok (in float.h)"])
)

AC_EGREP_CPP([^flt snan is broken],
[#include <float.h>
#ifndef FLT_MAX
flt snan is broken
#endif],
AC_DEFINE([FLT_MAX_IS_BROKEN],[1],[Define if FLT_MAX is not in float.h])
CASC_AC_LOG(["FLT_MAX is broken (not in float.h)"]),
CASC_AC_LOG(["FLT_MAX is ok (in float.h)"])
)

AC_EGREP_CPP([^dbl snan is broken],
[#include <float.h>
#ifndef DBL_MAX
dbl snan is broken
#endif],
AC_DEFINE([DBL_MAX_IS_BROKEN],[1],[Define if DBL_MAX is not in float.h])
CASC_AC_LOG(["DBL_MAX is broken (not in float.h)"]),
CASC_AC_LOG(["DBL_MAX is ok (in float.h)"])
)

# End macro CASC_IEEE_FLOAT
])dnl
