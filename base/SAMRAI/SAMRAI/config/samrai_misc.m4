dnl
dnl  File:           $HeadURL$
dnl  Package:        SAMRAI
dnl  Copyright:      (c) 1997-2010 Lawrence Livermore National Security, LLC
dnl  Date:           $Date$
dnl  Revision:       $LastChangedRevision$
dnl  Description:    Misc SAMRAI macros
dnl

AC_DEFUN_ONCE([SAMRAI_MISC],[

AC_ARG_ENABLE([box_counting],
[AS_HELP_STRING([--enable-box_counting],
   [Turns on Box and telemetry.])],
   [
      if test "x$enableval" = "xyes"; then
         AC_DEFINE(BOX_TELEMETRY,1,Enable Box counting)
      elif test "x$enableval" = "x"; then
         AC_DEFINE(BOX_TELEMETRY,1,Enable Box counting)
      fi
   ],)
]
)
