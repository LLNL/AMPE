dnl Define macros for supporting z compression library.
dnl $Id: support-libz.m4,v 1.1 2006/01/14 00:29:08 jafh Exp $

AC_DEFUN([BTNG_VAR_SET_LIBZ],[
dnl Provides support for the z compression library.
dnl
dnl Arguments are:
dnl 1. Name of variable to set to path where z compression library is installed.
dnl    Nothing is done if this variable is unset.
dnl    If you only want to look in default locations, set it to blank.
dnl 2. Name of the INCLUDES variable similar to the automake INCLUDES variable.
dnl    This variable is modified ONLY if it is NOT set.
dnl 3. Name of the LIBS variable similar to the automake LIBS variable.
dnl    This variable is modified ONLY if it is NOT set.
dnl
dnl If arg1 is defined, assume that the user wants z compression
dnl support.  Do so by assigning arg2 and arg3 if they are not defined.
dnl
# End macro BTNG_VAR_SET_LIBZ
if test "${$1+set}" = set ; then
  # Modify the output INCLUDES variable, if it is not set.
  if test ! "${$2+set}" = set ; then
    test -n "${$1}" && $2="-I${$1}/include"
  fi
  # Modify the output LIBS variable, if it is not set.
  if test ! "${$3+set}" = set ; then
    # Save LIBS for later recovery.
    btng_save_LIBS="$LIBS";
    # Extra libraries, if any, required by this check.
    btng_extra_libs="$libz_LIBS -lm"
    # If path is given, add path to extra flag for library search.
    test -n "${$1}" && btng_extra_libs="-L${$1}/lib $btng_extra_libs"
    # Look for library.
    AC_SEARCH_LIBS([zlibVersion],z,[
      BTNG_AC_LOG_VAR(LIBS,After finding libz flag)
      # Action if found ...
      # Extract modifications to LIB into library-specific LIBS variable.
      $3=`echo " $LIBS" | sed "s! $btng_save_LIBS!!"`;
      test -n "${$1}" && $3="-L${$1}/lib ${$3}"
      BTNG_AC_LOG_VAR($3, Found libz library flag)
      ],[
      # Action if NOT found ...
      BTNG_AC_LOG_VAR($3, Did not find libz library flag)
      ],[$btng_extra_libs])
    LIBS="$btng_save_LIBS";	# Restore global-use variable.
    unset btng_extra_libs
    unset btng_save_LIBS
  else
    BTNG_AC_LOG(Not looking for libz because $3 is already set)
  fi
fi
# End macro BTNG_VAR_SET_LIBZ
])dnl



AC_DEFUN([BTNG_SUPPORT_LIBZ],[
dnl Support z compression library by setting the variables
dnl libz_PREFIX, libz_INCLUDES, and libz_LIBS.
dnl Arg1: non-empty if you want the default to be on.
dnl
# Begin macro BTNG_SUPPORT_LIBZ

BTNG_ARG_WITH_ENV_WRAPPER(libz, libz_PREFIX,
ifelse($1,,
[  --with-libz[=PATH]
			Use z compression and optionally specify where
			it is installed.],
[  --without-libz	Do not use the z compression library.]),
ifelse($1,,unset libz_PREFIX; test "${with_libz+set}" = set && libz_PREFIX=,libz_PREFIX=))

BTNG_ARG_WITH_PREFIX(libz-includes,libz_INCLUDES,
[  --with-libz-includes=STRING
			Specify the INCLUDES flags for z compression.
			If not specified, and --with-libz=PATH is,
			this defaults to "-IPATH/include".])dnl

BTNG_ARG_WITH_PREFIX(libz-libs,libz_LIBS,
[  --with-libz-libs=STRING
			Specify LIBS flags for z compression.
			If not specified, and --with-libz=PATH is,
			this defaults to "-LPATH/lib -llapack -lblas".])dnl

BTNG_VAR_SET_LIBZ(libz_PREFIX,libz_INCLUDES,libz_LIBS)
# End macro BTNG_SUPPORT_LIBZ
])
