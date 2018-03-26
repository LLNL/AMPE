dnl $Id: support-blaslapack.m4,v 1.1 2006/01/14 00:29:08 jafh Exp $

AC_DEFUN([BTNG_VAR_SET_BLAS],[
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
      BTNG_AC_LOG_VAR(LIBS,After finding blas flag)
      # Action if found ...
      # Extract modifications to LIB into library-specific LIBS variable.
      $3=`echo " $LIBS" | sed "s! $btng_save_LIBS!!"`;
      test -n "${$1}" && $3="-L${$1}/lib ${$3}"
      BTNG_AC_LOG_VAR($3, Found blas library flag)
      ],[
      # Action if NOT found ...
      BTNG_AC_LOG_VAR($3, Did not find blas library flag)
      ],[$btng_extra_libs])
    LIBS="$btng_save_LIBS";	# Restore global-use variable.
    unset btng_extra_libs
    unset btng_save_LIBS
  else
    BTNG_AC_LOG(Not looking for blas because $3 is already set)
  fi
fi
])dnl

AC_DEFUN([BTNG_VAR_SET_LAPACK],[
dnl Provides support for the lapack library.
dnl
dnl Arguments are:
dnl 1. Name of variable to set to path where lapack are installed.
dnl    Nothig is done if this variable is unset.
dnl 2. Name of the INCLUDES variable similar to the automake INCLUDES variable.
dnl    This variable is modified ONLY if it is NOT set.
dnl 3. Name of the LIBS variable similar to the automake LIBS variable.
dnl    This variable is modified ONLY if it is NOT set.
dnl
dnl If arg1 is defined, assume that the user wants lapack
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
    AC_SEARCH_LIBS([xerbla_],lapack,[
      BTNG_AC_LOG_VAR(LIBS,After finding lapack flag)
      # Action if found ...
      # Extract modifications to LIB into library-specific LIBS variable.
      $3=`echo " $LIBS" | sed "s! $btng_save_LIBS!!"`;
      test -n "${$1}" && $3="-L${$1}/lib ${$3}"
      BTNG_AC_LOG_VAR($3, Found lapack library flag)
      ],[
      # Action if NOT found ...
      BTNG_AC_LOG_VAR($3, Did not find lapack library flag)
      ],[$btng_extra_libs])
    LIBS="$btng_save_LIBS";	# Restore global-use variable.
    unset btng_extra_libs
    unset btng_save_LIBS
  else
    BTNG_AC_LOG(Not looking for lapack because $3 is already set)
  fi
fi
])dnl




AC_DEFUN([BTNG_SUPPORT_BLASLAPACK],[
dnl Support blas and lapack libraries by setting the variables
dnl blaslapack_PREFIX, blaslapack_INCLUDES, and blaslapack_LIBS.
dnl Arg1: non-empty if you want the default to be on.
dnl
# Begin macro BTNG_SUPPORT_BLASLAPACK

BTNG_ARG_WITH_ENV_WRAPPER(blaslapack, blaslapack_PREFIX,
ifelse($1,,
[  --with-blaslapack[=PATH]
			Use blas and lapack and optionally specify where
			they are installed.],
[  --without-blaslapack	Do not use the blas and lapack libraries.]),
if test "${with_blaslapack+set}" = set; then
  blaslapack_PREFIX=
else
ifelse($1,,unset blaslapack_PREFIX,blaslapack_PREFIX=)
fi
)

BTNG_ARG_WITH_PREFIX(blaslapack-includes,blaslapack_INCLUDES,
[  --with-blaslapack-includes=STRING
			Specify the INCLUDES flags for blas and lapack.
			If not specified, and --with-blaslapack=PATH is,
			this defaults to "-IPATH/include".])dnl

BTNG_ARG_WITH_PREFIX(blaslapack-libs,blaslapack_LIBS,
[  --with-blaslapack-libs=STRING
			Specify LIBS flags for blas and lapack.
			If not specified, and --with-blaslapack=PATH is,
			this defaults to "-LPATH/lib -llapack -lblas".])dnl

# Look for blas and lapack flags separately.
# Then combine them.
if test "${blaslapack_PREFIX+set}" = set; then
  blas_PREFIX="$blaslapack_PREFIX"
  lapack_PREFIX="$blaslapack_PREFIX"
fi
BTNG_VAR_SET_BLAS(blas_PREFIX,blas_INCLUDES,blas_LIBS)
BTNG_VAR_SET_LAPACK(lapack_PREFIX,lapack_INCLUDES,lapack_LIBS)

# Combine include flags from blas and lapack without repeating -I flags.
if test -n "$blas_INCLUDES" || test -n "$lapack_INCLUDES" ; then
unset blaslapack_INCLUDES
for i in $lapack_INCLUDES $blas_INCLUDES ; do
  case $i in
  -I*) 
    echo "${blaslapack_INCLUDES}" | egrep "($i\$|$i[[ 	]] *)" > /dev/null	\
      && continue	# already contains the flag we want to add.
    blaslapack_INCLUDES="$blaslapack_INCLUDES $i"
    ;;
  esac
done
fi

# Combine lib flags from blas and lapack without repeating -L flags.
if test -n "$blas_LIBS" || test -n "$lapack_LIBS" ; then
unset blaslapack_LIBS
for i in $lapack_LIBS $blas_LIBS ; do
  case $i in
  -L*) 
    echo "${blaslapack_LIBS}" | egrep "($i\$|$i[[ 	]] *)" > /dev/null	\
      && continue	# already contains the flag we want to add.
    blaslapack_LIBS="$blaslapack_LIBS $i"
    ;;
  esac
done
for i in $lapack_LIBS $blas_LIBS ; do
  case $i in
  -L*) : ;;
  *) blaslapack_LIBS="$blaslapack_LIBS $i";
  esac
done
fi

BTNG_AC_LOG_VAR(blaslapack_PREFIX blaslapack_INCLUDES blaslapack_LIBS)
# End macro BTNG_SUPPORT_BLASLAPACK
])
