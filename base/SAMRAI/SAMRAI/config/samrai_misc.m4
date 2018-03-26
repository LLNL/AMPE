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
   [Turns on Box and MappedBox telemetry.])],
   [
      if test "x$enableval" = "xyes"; then
         CPPFLAGS_EXTRA="-DBOX_TELEMETRY $CPPFLAGS_EXTRA"
      elif test "x$enableval" = "x"; then
         CPPFLAGS_EXTRA="-DBOX_TELEMETRY $CPPFLAGS_EXTRA"
      fi
   ],)
]
)
