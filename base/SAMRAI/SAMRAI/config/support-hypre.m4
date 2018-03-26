dnl $Id$

dnl Define macros for supporting HYPRE.


AC_DEFUN([CASC_SUPPORT_HYPRE],[
dnl Support hypre libraries by setting the variables
dnl hypre_PREFIX, hypre_INCLUDES, and hypre_LIBS.
dnl Arg1: empty if you want the default to be off.
dnl
# Begin macro CASC_SUPPORT_HYPRE
CASC_ARG_WITH_ENV_WRAPPER(hypre, hypre_PREFIX,
ifelse($1,,
[  --with-hypre[=PATH]	Use HYPRE and optionally specify where it is installed.],
[  --without-hypre	Do not use the HYPRe library.]),
ifelse($1,,if test "$with_hypre" = '' ; then unset hypre_PREFIX; else hypre_PREFIX=; fi, hypre_PREFIX=)
)
CASC_VAR_SET_HYPRE(hypre_PREFIX,hypre_INCLUDES,hypre_LIBS)
CASC_AC_LOG_VAR(hypre_PREFIX hypre_INCLUDES hypre_LIBS)
if test "${hypre_PREFIX+set}" = set; then
  btng_save_cppflags=$CPPFLAGS

  # Add hypre include flags to cpp so we can examine its header file.
  CPPFLAGS="$hypre_INCLUDES $CPPFLAGS"
  CASC_AC_LOG_VAR(hypre_INCLUDES CPPFLAGS)

  # Check if HYPRE header is ok.
  AC_CHECK_HEADER(HYPRE_config.h,:,
    [AC_MSG_ERROR(Problems checking HYPRE_config.h)])

  # Check if HYPRE was compiled with parallelism.
  AC_MSG_CHECKING(if hypre is serial or parallel)
  AC_EGREP_CPP([^HYPRE_SEQUENTIAL_IS_DEFINED$], [
#include <HYPRE_config.h>
#ifdef HYPRE_SEQUENTIAL
HYPRE_SEQUENTIAL_IS_DEFINED
#endif
    ],
    hypre_PARALLELISM=serial,
    hypre_PARALLELISM=parallel)
  AC_MSG_RESULT($hypre_PARALLELISM)

  # Reset cpp after checking hypre header file.
  CPPFLAGS=$btng_save_cppflags
  unset btng_save_cppflags

  CASC_AC_LOG_VAR(CPPFLAGS)
  CASC_AC_LOG_VAR(hypre_config_file hypre_PARALLELISM)
fi
# End macro CASC_SUPPORT_HYPRE
])


AC_DEFUN([CASC_VAR_SET_HYPRE],[
dnl Provides support for the blas and lapack libraries.
dnl
dnl Arguments are:
dnl 1. Name of variable to set to path where hypre is installed.
dnl    Nothig is done if this variable is unset.
dnl 2. Name of the INCLUDES variable similar to the automake INCLUDES variable.
dnl    This variable is modified ONLY if it is NOT set.
dnl 3. Name of the LIBS variable similar to the automake LIBS variable.
dnl    This variable is modified ONLY if it is NOT set.
dnl
dnl If arg1 is defined, assume that the user wants blas and lapack
dnl support.  Do so by assigning arg2 and arg3 if they are not defined.
dnl
# Begin macro CASC_VAR_SET_HYPRE
if test "${$1+set}" = set ; then
  if test ! "${$2+set}" = set ; then
    test -n "${$1}" && $2="-I${$1}/include"
  fi
  if test ! "${$3+set}" = set ; then
    $3='-lHYPRE'
    if test -n "${$1}" ; then
      for i in ${$3} ; do
	tmp_name=`echo $i | sed 's/^-l//'`
        if test ! -f "${$1}/lib/lib${tmp_name}.a" && \
          test ! -f "${$1}/lib/lib${tmp_name}.so"; then
          AC_MSG_WARN(Library file for ${tmp_name} is missing from ${$1}/lib.)
        fi
      done
      $3="-L${$1}/lib ${$3}"
    fi
  fi
fi
# End macro CASC_VAR_SET_HYPRE
])dnl
