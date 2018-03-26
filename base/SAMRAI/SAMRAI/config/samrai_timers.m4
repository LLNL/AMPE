dnl
dnl  File:           $HeadURL$
dnl  Package:        SAMRAI
dnl  Copyright:      (c) 1997-2012 Lawrence Livermore National Security, LLC
dnl  Date:           $Date$
dnl  Revision:       $LastChangedRevision$
dnl  Description:    Macro to control whether timers are compile into or
dnl                  out of SAMRAI
dnl               
dnl                  Variable: samrai_enable_timers
dnl                  DEFINES:  ENABLE_SAMRAI_TIMERS
dnl

AC_DEFUN_ONCE([SAMRAI_TIMERS],[

# Begin SAMRAI_TIMERS
# Defines ENABLE_SAMRAI_TIMERS if --enable-timers is specified.  This is
# turned on by default 

AC_MSG_CHECKING([if SAMRAI Timers are enabled])
AC_ARG_ENABLE([timers],
[AS_HELP_STRING([--disable-timers],
   [Disable SAMRAI Timers.])],
   [
      if test "x$enableval" = "xno"; then
         samrai_enable_timers="$enableval"
      fi
   ], [
      samrai_enable_timers="yes"
   ])

dnl By default, the timers are enabled.  Explicitly disabling them
dnl means that we will not compile the timers into the code.  All
dnl timer calls will essentially be no-ops.
dnl
dnl However, if timers are enabled, we'll define this symbol and
dnl timers will be built into the code. 

AS_IF([test "x$samrai_enable_timers" = "xyes"], [
   AC_DEFINE([ENABLE_SAMRAI_TIMERS],[1],[ENABLE_SAMRAI_TIMERS])
])
AC_MSG_RESULT([$samrai_enable_timers])
]
)
