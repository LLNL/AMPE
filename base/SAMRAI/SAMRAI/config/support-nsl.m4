dnl $Id$

AC_DEFUN([CASC_VAR_SET_NSL],[
dnl Provides support for the nsl library.
dnl
dnl Arguments are:
dnl 1. Name of variable to set to path where nsl are installed.
dnl    Nothing is done if this variable is unset.
dnl    If you only want to look in default locations, set it to blank.
dnl 2. Name of the INCLUDES variable similar to the automake INCLUDES variable.
dnl    This variable is modified ONLY if it is NOT set and the path
dnl    is non-blank.
dnl 3. Name of the LIBS variable similar to the automake LIBS variable.
dnl    This variable is modified ONLY if it is NOT set.
dnl    If the library cannot be found, this remains unset.
dnl
dnl If arg1 is defined, assume that the user wants nsl support.
dnl Do so by assigning arg2 and arg3 if they are not defined.
dnl
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
    AC_SEARCH_LIBS([getnetname],nsl,[
      CASC_AC_LOG_VAR(LIBS,After finding nsl flag)
      # Action if found ...
      # Extract modifications to LIB into library-specific LIBS variable.
      $3=`echo " $LIBS" | sed "s! $btng_save_LIBS!!"`;
      test -n "${$1}" && $3="-L${$1}/lib ${$3}"
      CASC_AC_LOG_VAR($3, Found nsl library flag)
      ],[
      # Action if NOT found ...
      CASC_AC_LOG_VAR($3, Did not find nsl library flag)
      ],[$btng_extra_libs])
    LIBS="$btng_save_LIBS";	# Restore global-use variable.
    unset btng_extra_libs
    unset btng_save_LIBS
  else
    CASC_AC_LOG(Not looking for nsl because $3 is already set)
  fi
fi
])dnl



AC_DEFUN([CASC_SUPPORT_NSL],[
dnl Support nsl library by setting the variables
dnl nsl_PREFIX, nsl_INCLUDES, and nsl_LIBS.
dnl Arg1: non-empty if you want the default to be on.
dnl
# Begin macro CASC_SUPPORT_NSL

CASC_ARG_WITH_ENV_WRAPPER(nsl, nsl_PREFIX,
ifelse($1,,
[  --with-nsl[=PATH]
			Use nsl and optionally specify where
			it is installed.],
[  --without-nsl	Do not use the nsl library.]),
ifelse($1,,unset nsl_PREFIX; test "${with_nsl+set}" = set && nsl_PREFIX=,nsl_PREFIX=))
CASC_AC_LOG_VAR(nsl_PREFIX nsl_INCLUDES nsl_LIBS, before looking)

CASC_ARG_WITH_PREFIX(nsl-includes,nsl_INCLUDES,
[  --with-nsl-includes=STRING
			Specify the INCLUDES flags for nsl.
			If not specified, and --with-nsl=PATH is,
			this defaults to "-IPATH/include".])dnl

CASC_ARG_WITH_PREFIX(nsl-libs,nsl_LIBS,
[  --with-nsl-libs=STRING
			Specify LIBS flags for nsl.
			If not specified, and --with-nsl=PATH is,
			this defaults to "-LPATH/lib -lnsl".])dnl

CASC_VAR_SET_NSL(nsl_PREFIX,nsl_INCLUDES,nsl_LIBS)
# End macro CASC_SUPPORT_NSL
])
