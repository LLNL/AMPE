dnl $Id$

AC_DEFUN([CASC_VAR_SET_BLAS],[
dnl Provides support for the blas library.
dnl
dnl Arguments are:
dnl 1. Name of variable to set to path where blas are installed.
dnl    Nothing is done if this variable is unset.
dnl    If you only want to look in default locations, set it to blank.
dnl 2. Name of the INCLUDES variable similar to the automake INCLUDES variable.
dnl    This variable is modified ONLY if it is NOT set.
dnl 3. Name of the LIBS variable similar to the automake LIBS variable.
dnl    This variable is modified ONLY if it is NOT set.
dnl
dnl If arg1 is defined, assume that the user wants blas
dnl support.  Do so by assigning arg2 and arg3 if they are not defined.
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
    AC_SEARCH_LIBS([daxpy_],blas,[
      CASC_AC_LOG_VAR(LIBS,After finding blas flag)
      # Action if found ...
      # Extract modifications to LIB into library-specific LIBS variable.
      $3=`echo " $LIBS" | sed "s! $btng_save_LIBS!!"`;
      test -n "${$1}" && $3="-L${$1}/lib ${$3}"
      CASC_AC_LOG_VAR($3, Found blas library flag)
      ],[
      # Action if NOT found ...
      CASC_AC_LOG_VAR($3, Did not find blas library flag)
      AC_MSG_WARN(
[I could not systematically find the name of
the blas library so I am using -lblas instead.])
      $3="-lblas"
      test -n "${$1}" &&	\
	$3="-L${$1}/lib ${$3}"	# Add path flag to output variable.
      ],[$btng_extra_libs])
    LIBS="$btng_save_LIBS";	# Restore global-use variable.
    unset btng_extra_libs
    unset btng_save_LIBS
  else
    CASC_AC_LOG(Not looking for blas because $3 is already set)
  fi
fi
])dnl




AC_DEFUN([CASC_SUPPORT_BLAS],[
dnl Support blas library by setting the variables
dnl blas_PREFIX, blas_INCLUDES, and blas_LIBS.
dnl Arg1: non-empty if you want the default to be on.
dnl
# Begin macro CASC_SUPPORT_BLAS

CASC_ARG_WITH_ENV_WRAPPER(blas, blas_PREFIX,
ifelse($1,,
[  --with-blas[=PATH]
			Use blas and optionally specify where
			they are installed.],
[  --without-blas	Do not use the blas library.]),
if test "${with_blas+set}" = set; then
  blas_PREFIX=
else
ifelse($1,,unset blas_PREFIX,blas_PREFIX=)
fi
)

CASC_ARG_WITH_PREFIX(blas-includes,blas_INCLUDES,
[  --with-blas-includes=STRING
			Specify the INCLUDES flags for blas.
			If not specified, and --with-blas=PATH is,
			this defaults to "-IPATH/include".])dnl

CASC_ARG_WITH_PREFIX(blas-libs,blas_LIBS,
[  --with-blas-libs=STRING
			Specify LIBS flags for blas.
			If not specified, and --with-blas=PATH is,
			this defaults to "-LPATH/lib -lblas".])dnl

CASC_VAR_SET_BLAS(blas_PREFIX,blas_INCLUDES,blas_LIBS)

CASC_AC_LOG_VAR(blas_PREFIX blas_INCLUDES blas_LIBS)
# End macro CASC_SUPPORT_BLAS
])
