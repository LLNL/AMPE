dnl Define macros for supporting z compression library.
dnl $Id$


AC_DEFUN([CASC_SUPPORT_SPOOLES],[
dnl Support spooles library by setting the variables
dnl spooles_PREFIX, spooles_INCLUDES, and spooles_LIBS.
dnl Arg1: non-empty if you want the default to be on.
dnl
# Begin macro CASC_SUPPORT_SPOOLES

ifelse($1,,unset spooles_PREFIX,spooles_PREFIX=)

AC_ARG_WITH(spooles,
ifelse($1,,
[  --with-spooles[=PATH]
			Use spooles library and optionally specify where
			it is installed.],
[  --without-spooles	Do not use the spooles library.]),
if test "${with_spooles+set}" = yes ; then
  spooles_PREFIX=
elif test "${with_spooles}" = no; then
  unset spooles_PREFIX;
else
  spooles_PREFIX="${with_spooles}"
fi
)

if test "${spooles_PREFIX+set}" = set; then
  # Set spooles_LIBS and spooles_INCLUDE if they are not already set.
  # Note that we expect library archives and headers to be
  # directly under spooles_PREFIX rather than subdirectories
  # lib and include of spooles_PREFIX.
  if test ! "${spooles_LIBS+set}" = set; then
    if test ! "${spooles_PREFIX}" = ''; then
      spooles_LIBS="-L${spooles_PREFIX}"
    fi
    spooles_LIBS="${spooles_LIBS} -lspoolesMPI -lspooles"
  fi
  if test ! "${spooles_INCLUDES+set}" = set ; then
    if test ! "${spooles_PREFIX}" = ''; then
      spooles_INCLUDES="-I${spooles_PREFIX}"
    fi
  fi
fi

CASC_AC_LOG_VAR(spooles_PREFIX spooles_LIBS spooles_INCLUDES)

# End macro CASC_SUPPORT_SPOOLES
])
