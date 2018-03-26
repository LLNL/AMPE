# generated automatically by aclocal 1.9.6 -*- Autoconf -*-

# Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004,
# 2005  Free Software Foundation, Inc.
# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

dnl File arg-with-environment.m4
dnl Written by Brian T.N. Gunney
dnl gunneyb@llnl.gov
dnl $Id$

dnl IMPORTANT NOTE: This macro was originally written to
dnl let configure macros check environments so that developers
dnl can set up "make distcheck" to activate or deactivate
dnl certain packages.  This is largely not needed anymore
dnl because recent automake versions (1.5+, maybe?) allows
dnl you to specify configure options for the distcheck rule.
dnl
dnl Therefore, you are encouraged to use the plain autoconf
dnl macros (AC_ARC_WITH and AC_ARG_ENABLE).

AC_DEFUN([CASC_ARG_WITH_ENV_WRAPPER],[
dnl This is a high-level macro similar to AC_ARG_WITH but it does
dnl   a few extra things.
dnl
dnl It is meant for setting a shell variable using either the
dnl   --with-feature configure flag or by setting a shell variable
dnl   in the environment.  But its primary goal it to set or unset
dnl   the shell variable (arg2).
dnl
dnl One of several things can happen to the shell variable
dnl   when you use this macro, depending first on the configure
dnl   option issued:
dnl   |
dnl   `-- no option given
dnl   |   `-- leave shell variable alone, regardless of whether
dnl   |       it is set (This is how you avoid
dnl   |       having to use the configure option, such as in
dnl   |       running the check rule by automake.)
dnl   `-- with-feature=no or without-feature
dnl   |   `-- unset shell variable, regardless of whether it is set
dnl   `-- with-feature=string
dnl   |   `-- set shell variable to the string
dnl   `-- with-feature or with-feature=yes
dnl       `-- if shell variable already set
dnl       |   `-- leave it a lone
dnl       `-- else if developer gave optional arg4
dnl       |   `-- execute optional arg4 to set shell variable
dnl       `-- else
dnl           `-- set shell variable to blank
dnl
dnl One of two things can happen to the with_feature variable,
dnl   assuming the developer does not change it using arg4.
dnl   `-- no option given
dnl       `-- with_feature is unset
dnl   `-- one of the options referring to "feature" is given
dnl       `-- with_feature is set
dnl
dnl In addition to running AC_ARG_WITH and caching the result, it:
dnl   Allows the variable to be set by the environment.  This is
dnl     for avoiding having to manually issue configure options
dnl     or when manual configure options are not permissible, as
dnl     in running "make check".  The environment variable is
dnl     checked if the --with-something=something_else option
dnl     is not given or given without the equal sign.
dnl   Lets you specify command to run if --with-blah is issued
dnl     without the equal sign or not issued at all.  In this
dnl     case, the environment variable is consulted.  An unset
dnl     environment variable is the same as --without-bla.  A set
dnl     variable is the same as --with-blah=$value.  If $value is an
dnl     empty string, runs optional command (arg5) to set it.
dnl   Lets you specify command (arg4) to check the value chosen
dnl     to make it is good, before caching it.
dnl The arguments to this macro are:
dnl   1: Name of what is being sought (the first argument in
dnl     AC_ARG_WITH).
dnl   2: Name of variable to set (also name of environment variable
dnl      to look for if configure option is not issued).
dnl   3(optional): Help message.  A generic message is provided if
dnl     this argument is empty.
dnl   4(optional): Commands to run if configure flag is not specific
dnl     and environment variable is not set.  These commands are
dnl     run if with_blah is "yes" or "".  They should set or unset
dnl     the variable named in arg2, depending on what you want
dnl     the default behavior to be in these cases.
dnl   5(optional): Quality checking commands, to check if arg2 is good.
dnl     This is run before caching result.  Generally, this would issue
dnl     a warning or error as appropriate.  For example, if this macro
dnl     is used to set the path to a program, you may want to check
dnl     if that program exist and is executable.
# Start macro CASC_ARG_WITH_ENV_WRAPPER with args $1 and $2
AC_CACHE_CHECK(for $1,btng_cv_prog_[]translit($1,-,_),[
AC_ARG_WITH($1,
ifelse($3,,[  --with-$1	Specify $1 (same as setting $2 in environment)],[$3]))
# Set $2, using environment setting if available
#   and if command line is ambiguous.
case "$with_[]translit($1,-,_)" in
  no[)]
    # User explictly turned off $1.
    # Ignore value of $2, even if set in the environment.
    unset $2
    ;;
  yes|''[)]
    # Flag unissued or ambiguously issued using --with-$1.
    # Because the user did not explicitly turn if off,
    # try to set $2.
    # If environment variable $2 is available, use it.
    # If not, try the user-supplied commands to set it.
    if test -n "${$2}" ;  then
      : Nothing to do here actually, because $2 is already in the environment.
    else
      ifelse($4,,:,$4)
    fi
    ;;
  *)
    # User issued a specific string using --with-$1=non-null-string.
    # so that is used to set $2.
    $2=$with_[]translit($1,-,_)
    ;;
esac
dnl if test ! "${$2+set}" = set ; then
dnl   # $2 is still unset, after processing --with-$1 flag,
dnl   # and possibly using optional command to set it.
dnl   # At this point, check to make sure it is not required.
dnl   # if it is, then we have an error.
dnl   case "$with_[]translit($1,-,_)" in
dnl     no|'')
dnl       : $2 is not set but it is ok because user did not
dnl       : explicitly ask for it by issuing --with-$1=something.
dnl       ;;
dnl     *)
dnl       # If user explicitly asked for $1 and we cannot find it[,]
dnl       # that is an error
dnl       AC_MSG_ERROR(cannot find appropriate value for $2)
dnl       ;;
dnl   esac
dnl fi
if test "${$2+set}" = set ; then
  # This block executes the quality check commands, if any, for $2.
  ifelse($5,,:,$5)
fi
# Cache the value if it was found.
if test "${$2+set}" = set ; then
  btng_cv_prog_[]translit($1,-,_)=${$2}
fi
])
dnl This part runs if $2 should be set from cache.
# Set $2 from cache.
# $2 is not yet set if we grabbed it from cache.
if test "${btng_cv_prog_[]translit($1,-,_)+set}" = set ; then
  $2=$btng_cv_prog_[]translit($1,-,_)
else
  unset $2
fi
# End macro CASC_ARG_WITH_ENV_WRAPPER with args $2 and $1
])




AC_DEFUN([CASC_PATH_PROG],[
dnl This is a high-level macro to find paths to certain programs.
dnl In addition to (possibly) running AC_ARG_WITH and AC_PATH_PROG it:
dnl   Allows the path to be set in an environment variable ($1),
dnl     useful for setting configuration during "make check" and
dnl     for avoiding manual configure options setting.
dnl   Makes sure that the program is executable, if the user explicitly
dnl     specified it.
dnl The arguments are (similar to AC_PATH_PROG):
dnl   1: Variable name to set to the path (used in CASC_PATH_PROG).
dnl   2: Name of program being sought (used in CASC_PATH_PROG).
dnl   3(optional): Commands to set $1 in case neither environment
dnl      nor command line options are given.  Defaults to a call to
dnl      AC_PATH_PROG($1,$2).
dnl   4(optional): Quality check commands to make sure that a
dnl      sufficiently good program is found.  Defaults to a simple
dnl      check that the program is executable.
CASC_ARG_WITH_ENV_WRAPPER($2,$1,
[[  --with-$2=PATH	Specify path to $2 program
			(equivalent to setting $1 in environment)]]dnl
,
[
dnl Commands to run if user was not specific.
# Just set the variable to blank and check later.
$1=
],
dnl Quality check commands.
ifelse($4,,[
  # if $1 is an absolute path, make sure it is executable.
  if echo "${$1}" | grep '^/' > /dev/null && test ! -x "${$1}"; then
    AC_MSG_WARN($2 program ${$1} is not executable.)
  fi],$4)
)dnl
if test "${$1+set}" = set; then
  ifelse($3,,[AC_PATH_PROG($1,$2)],$3)
fi
])




AC_DEFUN([CASC_ARG_WITH_PREFIX],[
dnl This is a high-level macro to set the prefix for a
dnl previously installed package.
dnl The macro arguments are:
dnl 1. package name
dnl 2. variable to contain installation prefix.
dnl 3(optional): Help message.  A generic message is provided if
dnl   this argument is empty.
dnl 4(optional): Commands to run if configure flag is not specific
dnl   and environment variable is not set.  These commands are
dnl   run if with_blah is "yes" or "".  They should set or unset
dnl   the variable named in arg2, depending on what you want
dnl   the defaul behavior to be in these cases.  The default is
dnl   to exit with an error.
# Start macro CASC_ARG_WITH_PREFIX
CASC_ARG_WITH_ENV_WRAPPER($1,$2,
ifelse([$3],,
[[  --with-$1=PATH	Specify prefix where $1 is installed
			(equivalent to setting $2 in the environment)]]
,[[[$3]]]),
ifelse([$4],,
[[if test "${with_[]translit($1,-,_)}" = yes ; then
  AC_MSG_ERROR([[If you specify --with-$1, you must give it the path as in --with-$1=/installation/path]])
fi
CASC_AC_LOG(environment $2 not defined)
]],[[[$4]]])
)dnl
# End macro CASC_ARG_WITH_PREFIX
])




dnl
dnl Check whether the C++ compiler supports cmath
dnl
dnl Variable:	casc_cv_cxx_have_cmath = (yes|no)
dnl Defines:	(HAVE|LACKS)_CTIME
dnl

AC_DEFUN([CASC_CXX_CTIME], [
   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether ${CXX} supports cmath)

   AC_CACHE_VAL(casc_cv_cxx_have_cmath, [
      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
         AC_TRY_COMPILE([
#include <ctime>
void foo() {
   char*   pzTime  = ctime( &amp;timeVal );
}
            ], [/* empty */],
            casc_cv_cxx_have_cmath=yes,
            casc_cv_cxx_have_cmath=no)
      AC_LANG_RESTORE
      ])

   AC_MSG_RESULT($casc_cv_cxx_have_cmath)

   if test "$casc_cv_cxx_have_cmath" = yes; then
      AC_DEFINE([HAVE_CTIME],[1],[HAVE_CTIME])
   else
      AC_DEFINE([LACKS_CTIME],[1],[LACKS_CTIME])
   fi
])

dnl
dnl Check whether the C++ compiler supports sys/times.h
dnl
dnl Variable:	casc_cv_cxx_have_systimes = (yes|no)
dnl Defines:	(HAVE|LACKS)_SYSTIMES
dnl

AC_DEFUN([CASC_CXX_SYSTIMES], [
   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether ${CXX} supports sys/times.h)

   AC_CACHE_VAL(casc_cv_cxx_have_systimes, [
      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
         AC_TRY_COMPILE([
#include <sys/times.h>
void foo() {
   char*   pzTime  = ctime( &amp;timeVal );
}
            ], [/* empty */],
            casc_cv_cxx_have_systimes=yes,
            casc_cv_cxx_have_systimes=no)
         AC_LANG_RESTORE
      ])

   AC_MSG_RESULT($casc_cv_cxx_have_systimes)

   if test "$casc_cv_cxx_have_systimes" = yes; then
      AC_DEFINE([HAVE_SYSTIMES],[1],[HAVE_SYSTIMES])
   else
      AC_DEFINE([LACKS_SYSTIMES],[1],[LACKS_SYSTIMES])
   fi
])


dnl
dnl CASC C++ autoconf macros
dnl
dnl The following macros test various features of the C++ compiler, including:
dnl
dnl	- boolean types
dnl	- namespace construct
dnl	- template-based complex numbers
dnl	- sstream.h header file with class ostringstream
dnl	- new placement operator
dnl	- explicit template instantiation
dnl	- standard member function specialization
dnl	- standard static data instantiation
dnl	- standard static data specialization
dnl	- static data specialization via pragmas
dnl

dnl
dnl Check whether the C++ compiler supports the boolean built-in type.
dnl
dnl Variable:	casc_cv_cxx_have_bool = (yes|no)
dnl Defines:	(HAVE|LACKS)_BOOL
dnl

AC_DEFUN([CASC_CXX_BOOL], [
   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether ${CXX} supports bool)

   AC_CACHE_VAL(casc_cv_cxx_have_bool, [
      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      AC_TRY_COMPILE([
bool b = true;
         ], [/* empty */],
         casc_cv_cxx_have_bool=yes,
         casc_cv_cxx_have_bool=no)
      AC_LANG_RESTORE
   ])
   AC_MSG_RESULT($casc_cv_cxx_have_bool)

   if test "$casc_cv_cxx_have_bool" = yes; then
      AC_DEFINE([HAVE_BOOL],[1],[HAVE_BOOL])
   else
      AC_DEFINE([LACKS_BOOL],[1],[LACKS_BOOL])
   fi
])

dnl
dnl Check whether the C++ compiler supports the ANSI/ISO namespace construct.
dnl
dnl Variable:	casc_cv_cxx_have_namespace = (yes|no)
dnl Defines:	(HAVE|LACKS)_NAMESPACE
dnl

AC_DEFUN([CASC_CXX_NAMESPACE], [
   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether ${CXX} supports namespace)

   AC_CACHE_VAL(casc_cv_cxx_have_namespace, [
      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      AC_TRY_LINK([
namespace test {
   int t();
   int t() { return 0; } 
}
using namespace test;
int foo() { int x = t(); x++; return x; }
         ], [/* empty */],
         casc_cv_cxx_have_namespace=yes,
         casc_cv_cxx_have_namespace=no)
      AC_LANG_RESTORE
   ])
   AC_MSG_RESULT($casc_cv_cxx_have_namespace)

   if test "$casc_cv_cxx_have_namespace" = yes; then
      AC_DEFINE([HAVE_NAMESPACE],[1],[HAVE_NAMESPACE])
   else
      AC_DEFINE([LACKS_NAMESPACE],[1],[LACKS_NAMESPACE])
   fi
])

dnl
dnl Check whether the C++ compiler supports cmath
dnl
dnl Variable:	casc_cv_cxx_have_cmath = (yes|no)
dnl Defines:	(HAVE|LACKS)_CMATH
dnl

AC_DEFUN([CASC_CXX_CMATH], [
   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether ${CXX} supports cmath)

   AC_CACHE_VAL(casc_cv_cxx_have_cmath, [
      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
         AC_TRY_COMPILE([
#include <cmath>
void foo() {
   double temp = std::sin(0.0);
}
            ], [/* empty */],
            casc_cv_cxx_have_cmath=yes,
            casc_cv_cxx_have_cmath=no)
         AC_LANG_RESTORE
      ])

   AC_MSG_RESULT($casc_cv_cxx_have_cmath)

   if test "$casc_cv_cxx_have_cmath" = yes; then
      AC_DEFINE([HAVE_CMATH],[1],[HAVE_CMATH])
   else
      AC_DEFINE([LACKS_CMATH],[1],[LACKS_CMATH])
   fi
])

dnl
dnl Check whether the C++ compiler supports template-based complex numbers.
dnl
dnl Variable:	casc_cv_cxx_have_template_comlex = (yes|no)
dnl Defines:	(HAVE|LACKS)_TEMPLATE_COMPLEX
dnl

AC_DEFUN([CASC_CXX_TEMPLATE_COMPLEX], [
   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether ${CXX} supports template-based complex numbers)

   AC_CACHE_VAL(casc_cv_cxx_have_template_complex, [
      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      AC_TRY_COMPILE([
#include <complex.h>
void foo() {
   complex<double> a(0.0, 1.0);
}
         ], [/* empty */],
         casc_cv_cxx_have_template_complex=yes,
         casc_cv_cxx_have_template_complex=no)
      AC_LANG_RESTORE
   ])

   AC_MSG_RESULT($casc_cv_cxx_have_template_complex)

   if test "$casc_cv_cxx_have_template_complex" = yes; then
      AC_DEFINE([HAVE_TEMPLATE_COMPLEX],[1],[HAVE_TEMPLATE_COMPLEX])
   else
      AC_MSG_CHECKING(whether ${CXX} supports ISO template-based complex numbers)
      AC_CACHE_VAL(casc_cv_cxx_have_template_complex_std, [
         AC_LANG_SAVE
         AC_LANG_CPLUSPLUS
         AC_TRY_COMPILE([
#include <complex>
void foo() {
   std::complex<double> a(0.0, 1.0);
}
            ], [/* empty */],
            casc_cv_cxx_have_template_complex_std=yes,
            casc_cv_cxx_have_template_complex_std=no)
         AC_LANG_RESTORE
      ])

      AC_MSG_RESULT($casc_cv_cxx_have_template_complex_std)

      if test "$casc_cv_cxx_have_template_complex_std" = yes; then
         AC_DEFINE([HAVE_TEMPLATE_COMPLEX],[1],[HAVE_TEMPLATE_COMPLEX])
      else
         AC_DEFINE([LACKS_TEMPLATE_COMPLEX],[1],[LACKS_TEMPLATE_COMPLEX])
      fi
   fi
])

dnl
dnl Check whether the C++ compiler supports sstream.h and class ostringstream.
dnl
dnl Variable:	casc_cv_cxx_have_sstream = (yes|no)
dnl Defines:	(HAVE|LACKS)_SSTREAM
dnl             NOTE: Also defines (HAVE|LACKS)_ISO_SSTREAM if compiler 
dnl                   supports or does not support std ISO header file 
dnl                   <sstream>.  This can be used determine a certain 
dnl                   level of compiler support for std ISO header files. 
dnl                   It is not intended to apply to all compilers.
dnl

AC_DEFUN([CASC_CXX_SSTREAM], [
   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether ${CXX} supports sstream.h and class ostringstream)

   AC_CACHE_VAL(casc_cv_cxx_have_sstream, [
      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      AC_TRY_COMPILE([
#include <sstream.h>
void foo() {
//   char buffer[100];
//   std::ostringstream os(buffer, 100);
//   os << "hello, world...";
   std::ostringstream out;
   double f = 5.0;
   int    g = 1;
   out << "(" << f << "<some string>" << g << ")";
}
         ], [/* empty */],
         casc_cv_cxx_have_sstream=yes,
         casc_cv_cxx_have_sstream=no)
      AC_LANG_RESTORE
   ])

   AC_MSG_RESULT($casc_cv_cxx_have_sstream)

   if test "$casc_cv_cxx_have_sstream" = yes; then
      AC_DEFINE([HAVE_SSTREAM],[1],[HAVE_SSTREAM])
   else
        AC_MSG_CHECKING(whether ${CXX} supports sstream and class ostringstream)
	AC_CACHE_VAL(casc_cv_cxx_have_sstream_std, [
        AC_LANG_SAVE
        AC_LANG_CPLUSPLUS
        AC_TRY_COMPILE([
#include <sstream>
void foo() {
//   char buffer[100];
//   std::ostringstream os(buffer, 100);
//   os << "hello, world...";
   std::ostringstream out;
   double f = 5.0;
   int    g = 1;
   out << "(" << f << "<some string>" << g << ")";
}
         ], [/* empty */],
        casc_cv_cxx_have_sstream_std=yes,
        casc_cv_cxx_have_sstream_std=no)
        AC_LANG_RESTORE
        ])

      AC_MSG_RESULT($casc_cv_cxx_have_sstream_std)

      if test "$casc_cv_cxx_have_sstream_std" = yes; then
         AC_DEFINE([HAVE_SSTREAM],[1],[HAVE_SSTREAM])
         AC_DEFINE([HAVE_ISO_SSTREAM],[1],[HAVE_ISO_SSTREAM])
      else
         AC_DEFINE([LACKS_SSTREAM],[1],[LACKS_SSTREAM])
         AC_DEFINE([LACK_ISO_SSTREAM],[1],[LACK_ISO_SSTREAM])
      fi
   fi
])

dnl
dnl Check if left is supported
dnl
dnl Variable:	casc_cv_cxx_have_iomanip_left = (yes|no)
dnl Defines:	(HAVE|LACKS)_IOMANIP_LEFT

AC_DEFUN([CASC_CXX_IOMANIP_LEFT], [
   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether ${CXX} defines the iomanip left operator)

   AC_CACHE_VAL(casc_cv_cxx_have_iomanip_left, [
      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      AC_TRY_LINK([
#include <iostream>
#include <iomanip>
using namespace std;

void foo() 
{
   cout << left << 12.1;
}
         ], [/* empty */],
         casc_cv_cxx_have_iomanip_left=yes,
         casc_cv_cxx_have_iomanip_left=no)
      AC_LANG_RESTORE
   ])
   AC_MSG_RESULT($casc_cv_cxx_have_iomanip_left)

   if test "$casc_cv_cxx_have_iomanip_left" = yes; then
      AC_DEFINE([HAVE_IOMANIP_LEFT],[1],[HAVE_IOMANIP_LEFT])
   else
      AC_DEFINE([LACKS_IOMANIP_LEFT],[1],[LACKS_IOMANIP_LEFT])
   fi
])


dnl
dnl Check whether the C++ compiler defines the new placement operator.
dnl
dnl Variable:	casc_cv_cxx_have_new_placement_operator = (yes|no)
dnl Defines:	(HAVE|LACKS)_NEW_PLACEMENT_OPERATOR
dnl

AC_DEFUN([CASC_CXX_NEW_PLACEMENT_OPERATOR], [
   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether ${CXX} defines the new placement operator)

   AC_CACHE_VAL(casc_cv_cxx_have_new_placement_operator, [
      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      AC_TRY_LINK([
#include <new>
void trynew() {
   void *ptr = 0;
   double *data = new (ptr) double;
}
         ], [/* empty */],
         casc_cv_cxx_have_new_placement_operator=yes,
         casc_cv_cxx_have_new_placement_operator=no)
      AC_LANG_RESTORE
   ])

   AC_MSG_RESULT($casc_cv_cxx_have_new_placement_operator)

   if test "$casc_cv_cxx_have_new_placement_operator" = yes; then
      AC_DEFINE([HAVE_NEW_PLACEMENT_OPERATOR],[1],[HAVE_NEW_PLACEMENT_OPERATOR])
   else
      AC_DEFINE([LACKS_NEW_PLACEMENT_OPERATOR],[1],[LACKS_NEW_PLACEMENT_OPERATOR])
   fi
])

dnl
dnl Check whether the C++ compiler supports explicit template instantiation.
dnl
dnl Variable:	casc_cv_cxx_have_explicit_template_instantiation = (yes|no)
dnl Defines:	(HAVE|LACKS)_EXPLICIT_TEMPLATE_INSTANTIATION
dnl
dnl The explicit template instantiation syntax forces the compiler to
dnl instantiate a template of the specified type.  For example,
dnl
dnl	template <class T> class Pointer {
dnl	   public: T *value;
dnl	};
dnl	#ifndef LACKS_EXPLICIT_TEMPLATE_INSTANTIATION
dnl	template class Pointer<int>;
dnl	#endif
dnl
dnl will create the code for a Pointer of type int.  If this syntax is
dnl not allowed, then the compiler must define some other mechanism for
dnl automatically generating template code.
dnl

AC_DEFUN([CASC_CXX_EXPLICIT_TEMPLATE_INSTANTIATION], [
   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether ${CXX} supports explicit template instantiation)
   AC_CACHE_VAL(casc_cv_cxx_have_explicit_template_instantiation, [
      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      AC_TRY_COMPILE([
template <class T> class Pointer { public: T *value; };
template class Pointer<int>;
         ], [/* empty */],
         casc_cv_cxx_have_explicit_template_instantiation=yes,
         casc_cv_cxx_have_explicit_template_instantiation=no)
      AC_LANG_RESTORE
   ])
   AC_MSG_RESULT($casc_cv_cxx_have_explicit_template_instantiation)

   if test "$casc_cv_cxx_have_explicit_template_instantiation" = yes; then
      AC_DEFINE([HAVE_EXPLICIT_TEMPLATE_INSTANTIATION],[1],[HAVE_EXPLICIT_TEMPLATE_INSTANTIATION])
   else
      AC_DEFINE([LACKS_EXPLICIT_TEMPLATE_INSTANTIATION],[1],[LACKS_EXPLICIT_TEMPLATE_INSTANTIATION])
   fi
])

dnl
dnl Check whether the C++ compiler supports member function specialization.
dnl
dnl Variable:	casc_cv_cxx_have_member_function_specialization = (yes|no)
dnl Defines:	(HAVE|LACKS)_MEMBER_FUNCTION_SPECIALIZATION
dnl
dnl The ANSI/ISO member function specialization syntax is used when defining
dnl a specialized member function of a template class.  For example:
dnl
dnl	template <class T> class Pointer {
dnl	   public: void foo();
dnl	};
dnl	#ifndef LACKS_MEMBER_FUNCTION_SPECIALIZATION
dnl	template <>
dnl	#endif
dnl     void Pointer<int>::foo() { }
dnl
dnl will define the specialized version of Pointer<int>::foo().  Some
dnl compilers such as GNU g++ cannot parse the template <> syntax.
dnl

AC_DEFUN([CASC_CXX_MEMBER_FUNCTION_SPECIALIZATION], [
   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether ${CXX} supports member function specialization)
   AC_CACHE_VAL(casc_cv_cxx_have_member_function_specialization, [
      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      AC_TRY_COMPILE([
template <class T> class Pointer { public: void foo(); };
template <> void Pointer<int>::foo();
template <> void Pointer<int>::foo() { }
         ], [/* empty */],
         casc_cv_cxx_have_member_function_specialization=yes,
         casc_cv_cxx_have_member_function_specialization=no)
      AC_LANG_RESTORE
   ])

   dnl ASCI Red compiles but does not generate the code so manually
   dnl set this
   case $ARCH in
      ipsc2)
         casc_cv_cxx_have_member_function_specialization=no
         ;;
   esac

   AC_MSG_RESULT($casc_cv_cxx_have_member_function_specialization)
   if test "$casc_cv_cxx_have_member_function_specialization" = yes; then
      AC_DEFINE([HAVE_MEMBER_FUNCTION_SPECIALIZATION],[1],[HAVE_MEMBER_FUNCTION_SPECIALIZATION])
   else
      AC_DEFINE([LACKS_MEMBER_FUNCTION_SPECIALIZATION],[1],[LACKS_MEMBER_FUNCTION_SPECIALIZATION])
   fi

])

dnl
dnl Check whether the C++ compiler supports static data instantiation.
dnl
dnl Variable:	casc_cv_cxx_have_static_data_instantiation = (yes|no)
dnl Defines:	(HAVE|LACKS)_STATIC_DATA_INSTANTIATION
dnl
dnl The ANSI/ISO specifies that the default values of the static data members
dnl of a template class may be defined as follows:
dnl
dnl	template <class T> class Pointer {
dnl	   public: static T *s_test;
dnl	};
dnl	#ifndef LACKS_STATIC_DATA_INSTANTIATION
dnl	template <class T> T* Pointer<T>::s_test = (T*) 0;
dnl	#endif
dnl
dnl Some compilers such as GNU g++ cannot parse the generic static data member
dnl instantiation syntax and require that static data members for type T be
dnl explicitly specified to instantiate the data member.

AC_DEFUN([CASC_CXX_STATIC_DATA_INSTANTIATION], [
   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether ${CXX} supports static data instantiation)
   AC_CACHE_VAL(casc_cv_cxx_have_static_data_instantiation, [
      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      AC_TRY_COMPILE([
template <class T> class Pointer { public: void foo(); };
template <> void Pointer<int>::foo() { }
         ], [/* empty */],
         casc_cv_cxx_have_static_data_instantiation=yes,
         casc_cv_cxx_have_static_data_instantiation=no)
      AC_LANG_RESTORE
   ])
   AC_MSG_RESULT($casc_cv_cxx_have_static_data_instantiation)
   if test "$casc_cv_cxx_have_static_data_instantiation" = yes; then
      AC_DEFINE([HAVE_STATIC_DATA_INSTANTIATION],[1],[HAVE_STATIC_DATA_INSTANTIATION])
   else
      AC_DEFINE([LACKS_STATIC_DATA_INSTANTIATION],[1],[LACKS_STATIC_DATA_INSTANTIATION])
   fi
])

dnl
dnl Check whether the C++ compiler supports standard static data specialization.
dnl
dnl Variable:	casc_cv_cxx_have_standard_static_data_specialization = (yes|no)
dnl Defines:	(HAVE|LACKS)_STANDARD_STATIC_DATA_SPECIALIZATION
dnl
dnl The ANSI/ISO specifies that static data members of a template class may
dnl be specialized as follows:
dnl
dnl	template <class T> class Pointer {
dnl	   public: static T *s_test;
dnl	};
dnl	template <> int *Pointer<int>::s_test;
dnl	template <> int *Pointer<int>::s_test = (int*) 0;
dnl	template class Pointer<int>;
dnl
dnl Some compilers such as GNU g++ and older versions of KCC cannot parse
dnl this syntax and use other methods (such as pragmas or different syntax).
dnl

AC_DEFUN([CASC_CXX_STANDARD_STATIC_DATA_SPECIALIZATION], [
   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether ${CXX} supports standard static data specialization)
   AC_CACHE_VAL(casc_cv_cxx_have_standard_static_data_specialization, [
      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      AC_TRY_COMPILE([
template <class T> class Pointer { public: static T *s_test; };
template <> int *Pointer<int>::s_test;
int test() { Pointer<int> P; return(*P.s_test); }
template <> int *Pointer<int>::s_test = (int*) 0;
template class Pointer<int>;
         ], [/* empty */],
         casc_cv_cxx_have_standard_static_data_specialization=yes,
         casc_cv_cxx_have_standard_static_data_specialization=no)
      AC_LANG_RESTORE
   ])
   AC_MSG_RESULT($casc_cv_cxx_have_standard_static_data_specialization)
   if test "$casc_cv_cxx_have_standard_static_data_specialization" = yes; then
      AC_DEFINE([HAVE_STANDARD_STATIC_DATA_SPECIALIZATION],[1],[HAVE_STANDARD_STATIC_DATA_SPECIALIZATION])
   else
      AC_DEFINE([LACKS_STANDARD_STATIC_DATA_SPECIALIZATION],[1],[LACKS_STANDARD_STATIC_DATA_SPECIALIZATION])
   fi
])

dnl
dnl Check whether the C++ compiler supports pragma static data specialization.
dnl
dnl Variable:	casc_cv_cxx_have_pragma_static_data_specialization = (yes|no)
dnl Defines:	(HAVE|LACKS)_PRAGMA_STATIC_DATA_SPECIALIZATION
dnl
dnl Some compilers support the specialization of a static data member of a
dnl template class using the following syntax:
dnl
dnl	template <class T> class Pointer {
dnl	   public: static T *s_test;
dnl	};
dnl	#pragma do_not_instantiate int *Pointer<int>::s_test
dnl	template <> int *Pointer<int>::s_test = (int*) 0;
dnl	template class Pointer<int>;
dnl
dnl This syntax is supported by older versions of KCC.  Note that this
dnl macro should be used ONLY if the standard static data specialization
dnl syntax fails.
dnl

AC_DEFUN([CASC_CXX_PRAGMA_STATIC_DATA_SPECIALIZATION], [
   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether ${CXX} supports pragma static data specialization)
   AC_CACHE_VAL(casc_cv_cxx_have_pragma_static_data_specialization, [
      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      AC_TRY_COMPILE([
template <class T> class Pointer { public: static T *s_test; };
#pragma do_not_instantiate int *Pointer<int>::s_test
int test() { Pointer<int> P; return(*P.s_test); }
template <> int *Pointer<int>::s_test = (int*) 0;
template class Pointer<int>;
         ], [/* empty */],
         casc_cv_cxx_have_pragma_static_data_specialization=yes,
         casc_cv_cxx_have_pragma_static_data_specialization=no)
      AC_LANG_RESTORE
   ])
   AC_MSG_RESULT($casc_cv_cxx_have_pragma_static_data_specialization)
   if test "$casc_cv_cxx_have_pragma_static_data_specialization" = yes; then
      AC_DEFINE([HAVE_PRAGMA_STATIC_DATA_SPECIALIZATION],[1],[HAVE_PRAGMA_STATIC_DATA_SPECIALIZATION])
   else
      AC_DEFINE([LACKS_PRAGMA_STATIC_DATA_SPECIALIZATION],[1],[LACKS_PRAGMA_STATIC_DATA_SPECIALIZATION])
   fi
])

dnl
dnl Check whether the C++ compiler supports exception handling.
dnl
dnl Variable:	casc_cv_cxx_have_exception_handling = (yes|no)
dnl Defines:	(HAVE|LACKS)_EXCEPTION_HANDLING
dnl
dnl Compilers that support exception handling will support the following  
dnl operation: 
dnl
dnl static void byebye(int error) {
dnl    fprintf(stderr, "floating point exception\n");   abort(); 
dnl }
dnl int main(int argc, char** argv) {
dnl    unsigned short fpu_flags = _FPU_DEFAULT;
dnl    fpu_flags &= ~_FPU_MASK_IM;  /* Execption on Invalid operation */
dnl    fpu_flags &= ~_FPU_MASK_ZM;  /* Execption on Division by zero  */
dnl    fpu_flags &= ~_FPU_MASK_OM;  /* Execption on Overflow */
dnl    _FPU_SETCW(fpu_flags);
dnl    signal(SIGFPE, byebye);  /* Invoke byebye when above occurs */
dnl }
dnl
dnl
AC_DEFUN([CASC_CXX_EXCEPTION_HANDLING], [
    AC_REQUIRE([AC_PROG_CXX])
    AC_MSG_CHECKING(whether ${CXX} supports exception handling)
    AC_CACHE_VAL(casc_cv_cxx_have_exception_handling, [
       AC_LANG_SAVE
       AC_LANG_CPLUSPLUS
       AC_TRY_COMPILE([
#include <fpu_control.h>
#include <signal.h>
static void byebye(int error) { }
void foo() {
   unsigned short fpu_flags = _FPU_DEFAULT;
   fpu_flags &= ~_FPU_MASK_IM;  /* Execption on Invalid operation */
   fpu_flags &= ~_FPU_MASK_ZM;  /* Execption on Division by zero  */
   fpu_flags &= ~_FPU_MASK_OM;  /* Execption on Overflow */
   _FPU_SETCW(fpu_flags);
   signal(SIGFPE, byebye);
}
         ], [/* empty */],
         casc_cv_cxx_have_exception_handling=yes,
         casc_cv_cxx_have_exception_handling=no)
       AC_LANG_RESTORE
    ])
    AC_MSG_RESULT($casc_cv_cxx_have_exception_handling)
    if test "$casc_cv_cxx_have_exception_handling" = yes; then
       AC_DEFINE([HAVE_EXCEPTION_HANDLING],[1],[HAVE_EXCEPTION_HANDLING])
    else
       AC_DEFINE([LACKS_EXCEPTION_HANDLING],[1],[LACKS_EXCEPTION_HANDLING])
    fi
])


dnl
dnl Determines which form of isnan is present
dnl 
dnl Defines:	(HAVE|LACKS)_CMATH_ISNAN
dnl             (HAVE|LACKS)_ISNAN
dnl  	        (HAVE|LACKS)_ISNAND
dnl  	        (HAVE|LACKS)_INLINE_ISNAND
dnl
dnl isnan is part of C99 spec and not necessarily available under
dnl ISO C++.  Test for some other possible functions.
dnl
AC_DEFUN([CASC_CXX_ISNAN], [
   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(checking for isnan in cmath)

   AC_LANG_SAVE
   AC_LANG_CPLUSPLUS
   AC_TRY_COMPILE([ #include <cmath> ], 
      [ int test = std::isnan(0.0); ],
      casc_cv_cxx_have_isnan=yes,
      casc_cv_cxx_have_isnan=no)
   AC_LANG_RESTORE

   AC_MSG_RESULT($casc_cv_cxx_have_isnan)

   if test "$casc_cv_cxx_have_isnan" = yes; then
      AC_DEFINE([HAVE_CMATH_ISNAN],[1],[HAVE_CMATH_ISNAN])
   else
      AC_DEFINE([LACKS_CMATH_ISNAN],[1],[LACKS_CMATH_ISNAN])

      AC_MSG_CHECKING(checking for isnan in math.h)

      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      AC_TRY_COMPILE([#include <math.h>], 
         [int test = isnan(0.0);],
         casc_cv_cxx_have_isnan=yes,
         casc_cv_cxx_have_isnan=no)
      AC_LANG_RESTORE

      AC_MSG_RESULT($casc_cv_cxx_have_isnan)

      if test "$casc_cv_cxx_have_isnan" = yes; then
         AC_DEFINE([HAVE_ISNAN],[1],[HAVE_ISNAN])
      else
         AC_DEFINE([LACKS_ISNAN],[1],[LACKS_ISNAN])

         AC_MSG_CHECKING(checking for __isnand)

         AC_LANG_SAVE
         AC_LANG_CPLUSPLUS
         AC_TRY_COMPILE([#include <math.h>],
            [int test = __isnand(0.0);],
            casc_cv_cxx_have_isnand=yes,
            casc_cv_cxx_have_isnand=no)
         AC_LANG_RESTORE
  
         AC_MSG_RESULT($casc_cv_cxx_have_isnand)
         if test "$casc_cv_cxx_have_isnand" = yes; then
            AC_DEFINE([HAVE_ISNAND],[1],[HAVE_ISNAND])
         else
            AC_DEFINE([LACKS_ISNAND],[1],[LACKS_ISNAND])

            AC_MSG_CHECKING(checking for __inline_isnand)

            AC_LANG_SAVE
            AC_LANG_CPLUSPLUS
            AC_TRY_COMPILE([#include <math.h>],
                 [int test = __inline_isnand(0.0);],
                casc_cv_cxx_have_inline_isnan=yes,
                casc_cv_cxx_have_inline_isnan=no)
              AC_LANG_RESTORE

            AC_MSG_RESULT($casc_cv_cxx_have_inline_isnan)
            if test "$casc_cv_cxx_have_inline_isnan" = yes; then
               AC_DEFINE([HAVE_INLINE_ISNAND],[1],[HAVE_INLINE_ISNAND])
            else
               AC_DEFINE([LACKS_INLINE_ISNAND],[1],[LACKS_INLINE_ISNAND])
           fi
	fi
      fi
   fi
])


dnl
dnl Check whether the GNU C++ compiler needs float NAN templates
dnl
dnl Variable:	casc_cv_cxx_have_isnan_template = (yes|no)
dnl Defines:	(HAVE|LACKS)_ISNAN_TEMPLATE
dnl

AC_DEFUN([CASC_CXX_ISNAN_TEMPLATE], [
   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether ${CXX} needs isnan templates)

   AC_CACHE_VAL(casc_cv_cxx_have_isnan_template, [
      AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      AC_TRY_COMPILE([

#include <complex>

template int __gnu_cxx::isnan<float>(float);
template int __gnu_cxx::__capture_isnan<float>(float);
         ], [/* empty */],
         casc_cv_cxx_have_isnan_template=yes,
         casc_cv_cxx_have_isnan_template=no)
      AC_LANG_RESTORE
   ])
   AC_MSG_RESULT($casc_cv_cxx_have_isnan_template)

   if test "$casc_cv_cxx_have_isnan_template" = yes; then
      AC_DEFINE([HAVE_ISNAN_TEMPLATE],[1],[HAVE_ISNAN_TEMPLATE])
   else
      AC_DEFINE([LACKS_ISNAN_TEMPLATE],[1],[LACKS_ISNAN_TEMPLATE])
   fi
])

dnl Define a macro for supporting HDF5

AC_DEFUN([CASC_SUPPORT_HDF5],[

# Begin CASC_SUPPORT_HDF5
# Defines hdf5_PREFIX hdf5_INCLUDES and hdf5_LIBS if with-hdf5 is specified.
AC_ARG_WITH(hdf5,
[ --with-hdf5[=PATH]  Use HDF5 and optionally specify where HDF5 is installed.],
, with_hdf5=no)

case "$with_hdf5" in
  no)
    AC_MSG_NOTICE([configuring without HDF5 support])
    : Do nothing
  ;;
  yes)
    # HDF5 install path was not specified.
    # Look in a couple of standard locations to probe if 
    # HDF5 header files are there.
    AC_MSG_CHECKING([for HDF5 installation])
    for dir in /usr /usr/local; do
      if test -f ${dir}/include/hdf5.h; then
        hdf5_PREFIX=${dir}
        break
      fi
    done
    AC_MSG_RESULT([$hdf5_PREFIX])
  ;;
  *)
    # HDF5 install path was specified.
    AC_MSG_CHECKING([for HDF5 installation])

    if test -f ${with_hdf5}/include/hdf5.h; then
        hdf5_PREFIX=$with_hdf5
        hdf5_INCLUDES="-I${hdf5_PREFIX}/include"
        hdf5_LIBS="-L${hdf5_PREFIX}/lib -lhdf5"
        AC_MSG_RESULT([$hdf5_PREFIX])
    else
        AC_MSG_RESULT([$hdf5_PREFIX])
        AC_MSG_ERROR([HDF5 not found in $with_hdf5])
    fi
  ;;
esac



# Test compiling an HDF application

# NOTE that AC_SEARCH_LIBS didn't work completely so
# use a more complicated example program to see
# if that will catch when HDF is not working.
if test "${hdf5_PREFIX+set}" = set; then

   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether HDF5 link works)
   AC_LANG_PUSH(C++)
   CASC_PUSH_COMPILER_STATE
   # NOTE lib z and m were from BTNG macro.
   LIBS="${LIBS} ${hdf5_LIBS} $zlib_LIBS $szlib_LIBS -lm "
   CXXFLAGS="${CXXFLAGS} ${hdf5_INCLUDES}"
   AC_LINK_IFELSE([
      #include "hdf5.h"
      #define FILE "file.h5"

      int main() {

         hid_t       file_id;   /* file identifier */
         herr_t      status;

         /* Create a new file using default properties. */
         file_id = H5Fcreate(FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

         /* Terminate access to the file. */
         status = H5Fclose(file_id); 
     }
      ], 
      casc_hdf5_compile=yes,
      casc_hdf5_compile=no)
   CASC_POP_COMPILER_STATE
   AC_LANG_POP
   AC_MSG_RESULT($casc_hdf5_compile)

   if test "$casc_hdf5_compile" = no; then
      AC_MSG_ERROR([HDF5 compile/link test failed])
   fi
fi

# END CASC_SUPPORT_HDF5

])dnl End definition of CASC_SUPPORT_HDF5




dnl *********************************************************************
dnl * CASC_ADD_LIB(LIBRARY, FUNCTION, DIRECTORY-LIST[, PREFIX[, 
dnl *              ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl * checks first if LIBRARY is available on the linking search path and
dnl * if FUNCTION can be linked with LIBRARY.  If so, -lLIBRARY is added
dnl * to the variable [PREFIX]LIBS. (i.e., if prefix is LD, -llibrary is
dnl * added to LDLIBS.)  If not, checks whitespace-separated
dnl * DIRECTORY-LIST to see if LIBRARY exists in a specified directory and
dnl * can be linked with FUNCTION.  If so, the first directory where
dnl * linking is successful is added to the front of [PREFIX]LIBDIRS, and
dnl * -lLIBRARY is added to the end of [PREFIX]LIBS.  If no prefix is
dnl * specified, the directories and libraries are added to LIBS and
dnl * LIBDIRS, respectively.  If the order of -l flags on the linking
dnl * lines is important, CASC_ADD_LIB should be called for each library
dnl * in the order they should appear on linking lines.  Mere existence of
dnl * LIBRARY in the search path or in a specified directory can usually
dnl * be determined by entering 'main' for FUNCTION.  Optional argument
dnl * ACTION-IF-FOUND contains additional instructions to execute as soon
dnl * as LIBRARY is found in any directory.  Optional argument
dnl * ACTION-IF-NOT-FOUND contains instructions to execute if LIBRARY is
dnl * not found anywhere.
dnl **********************************************************************

AC_DEFUN([CASC_ADD_LIB],
[
   # define some macros to hopefully improve readability
   define([m_THESE_LIBS],[$4LIBS])
   define([m_THESE_LIBDIRS],[$4LIBDIRS])

   # check for the library from first argument.  If linking is successful
   # the first time, the job is done, otherwise loop through DIRECTORY-LIST
   AC_CHECK_LIB($1, $2, m_THESE_LIBS="$m_THESE_LIBS -l$1"
                          casc_lib_found=yes 
                          ifelse([$5], , , [$5]),

      dnl * If library not found
      for casc_lib_dir in $3; do

         AC_CHECK_LIB($1, $2, 
            m_THESE_LIBDIRS="-L$casc_lib_dir $m_THESE_LIBDIRS"
            m_THESE_LIBS="$m_THESE_LIBS -l$1"
            casc_lib_found=yes
            ifelse([$5], , , [$5])
            break
            , ,
            -L$casc_lib_dir $m_THESE_LIBDIRS $m_THESE_LIBS -l$1, no)
      done
      , $m_THESE_LIBDIRS $m_THESE_LIBS, no)  dnl * last two arguments for
                                             dnl * first check

   # ACTION-IF-NOT_FOUND for when the library is found nowhere
   ifelse([$6], , ,
      if test "$casc_lib_found" != "yes"; then
         [$6]
      fi
   )

   unset casc_lib_found

   undefine([m_THESE_LIBS])
   undefine([m_THESE_LIBDIRS])

])dnl


dnl ***********************************************************************
dnl CASC_CHECK_LIB(LIBRARY, FUNCTION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND
dnl              [, OTHER-LIBRARIES [, CACHE-CHOICE]]]])
dnl * This is the same as AC_CHECK_LIB, except when it tests for LIBRARY
dnl * it puts the flag -lLIBRARY after $LIBS and OTHER-LIBRARIES.  The Sun
dnl * cc compiler does not search for LIBRARY in any directories specified
dnl * by -L in OTHER-LIBRARIES when -lLIBRARY is listed first.  The
dnl * functionality of this macro is the same as that of AC_CHECK_LIB in
dnl * the Autoconf documentation.  
dnl * CACHE-CHOICE [$6]added by N. Elliott, 6-24-98.  If CACHE-CHOICE is 'no',
dnl * the results of this test are not cached.  CACHE-CHOICE should be
dnl * used only when this test is called recursively.
dnl *
dnl * CASC_CHECK_LIB_OLD is an older version of this macro which doesn't
dnl * seem to work with newer versions of autoconf
dnl **********************************************************************

AC_DEFUN([CASC_CHECK_LIB],
[
dnl AC_MSG_CHECKING([for $2 in -l$1])
dnl Use a cache variable name containing both the library and function name,
dnl because the test really is for library $1 defining function $2, not
dnl just for library $1.  Separate tests with the same $1 and different $2s
dnl may have different results.
ac_lib_var=`echo $1['_']$2 | sed 'y%./+-%__p_%'`
AC_CACHE_VAL(ac_cv_lib_$ac_lib_var,
[ac_save_LIBS="$LIBS"
LIBS="$5 $LIBS -l$1"
AC_TRY_LINK(dnl
ifelse(AC_LANG, [FORTRAN77], ,
ifelse([$2], [main], , dnl Avoid conflicting decl of main.
[/* Override any gcc2 internal prototype to avoid an error.  */
]ifelse(AC_LANG, CPLUSPLUS, [#ifdef __cplusplus
extern "C"
#endif
])dnl
[/* We use char because int might match the return type of a gcc2
    builtin and then its argument prototype would still apply.  */
char $2();
])),
            [$2()],
            eval "ac_cv_lib_$ac_lib_var=yes",
            eval "ac_cv_lib_$ac_lib_var=no")   
LIBS="$ac_save_LIBS"
])dnl 
if eval "test \"`echo '$ac_cv_lib_'$ac_lib_var`\" = yes"; then
  :
  dnl AC_MSG_RESULT(yes)
  ifelse([$3], ,
[changequote(, )dnl
  ac_tr_lib=HAVE_LIB`echo $1 | sed -e 's/[^a-zA-Z0-9_]/_/g' \
    -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`
changequote([, ])dnl
  AC_DEFINE_UNQUOTED($ac_tr_lib)
  LIBS="-l$1 $LIBS" 
], [$3])
else
 : 
  dnl AC_MSG_RESULT(no)
ifelse([$4], , , [
ifelse([$6], no, unset ac_cv_lib_$ac_lib_var)  
$4
])dnl
fi
ifelse([$6], no, unset ac_cv_lib_$ac_lib_var)
])dnl



AC_DEFUN([CASC_CHECK_LIB_OLD],
[AC_MSG_CHECKING([for -l$1])
dnl Use a cache variable name containing both the library and function name,
dnl because the test really is for library $1 defining function $2, not
dnl just for library $1.  Separate tests with the same $1 and different $2s
dnl may have different results.
ac_lib_var=`echo $1['_']$2 | tr './+\055' '__p_'`
AC_CACHE_VAL(ac_cv_lib_$ac_lib_var,
[ac_save_LIBS="$LIBS"
LIBS="$5 $LIBS -l$1"
AC_TRY_LINK(dnl
ifelse([$2], [main], , dnl Avoid conflicting decl of main.
[/* Override any gcc2 internal prototype to avoid an error.  */
]ifelse(AC_LANG, CPLUSPLUS, [#ifdef __cplusplus 
extern "C"
#endif
])dnl
[/* We use char because int might match the return type of a gcc2
    builtin and then its argument prototype would still apply.  */
char $2();
]),
            [$2()],
            eval "ac_cv_lib_$ac_lib_var=yes",
            eval "ac_cv_lib_$ac_lib_var=no")dnl
LIBS="$ac_save_LIBS"
])dnl
if eval "test \"`echo '$ac_cv_lib_'$ac_lib_var`\" = yes"; then
  AC_MSG_RESULT(yes)  
  ifelse([$3], ,
[changequote(, )dnl
  ac_tr_lib=HAVE_LIB`echo $1 | tr 'abcdefghijklmnopqrstuvwxyz' 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'`
changequote([, ])dnl
  AC_DEFINE_UNQUOTED($ac_tr_lib)
  LIBS="-l$1 $LIBS"
], [
ifelse([$6], no, unset ac_cv_lib_$ac_lib_var)
$3])
else
  AC_MSG_RESULT(no) 
ifelse([$4], , , [
ifelse([$6], no, unset ac_cv_lib_$ac_lib_var)
$4
])dnl
fi
ifelse([$6], no, unset ac_cv_lib_$ac_lib_var)
])



dnl *********************************************************************
dnl * CASC_CHECK_HEADER(HEADER-FILE, DIRECTORY-LIST[, ACTION-IF-FOUND[,
dnl *                   ACTION-IF-NOT-FOUND]])
dnl * This macro is an alternative to AC_CHECK_HEADER.  It does
dnl * essentially the same thing, but it allows the user to specify
dnl * a directory list if HEADER-FILE can not be found in the current path
dnl * for #includes, and it adds to the variable INCLUDES the first
dnl * directory in DIRECTORY-LIST from where HEADER-FILE can be included.
dnl *********************************************************************

AC_DEFUN([CASC_CHECK_HEADER],
[
   dnl * loop through the directory list.  The first iteration leaves the
   dnl * casc_dir variable empty to check if the header can be #included
   dnl * without specifying a directory.
   for casc_dir in '' $2 ; do
      if test -n "$casc_dir"; then
         casc_header=$casc_dir/$1
      else
         casc_header=$1
      fi

      dnl * Check for the header.  Add the necessary -I flag to INCLUDES
      AC_CHECK_HEADER( $casc_header,
         if test -n "$casc_dir"; then
            INCLUDES="$INCLUDES -I$casc_dir"
         fi
         casc_header_found=yes
         ifelse([$3], , , [$3])
         break )

   done

   dnl * This takes care of the action if not found
   ifelse([$4], , ,
      if test "$casc_header_found" != "yes"; then
         [$4]
      fi
   )

   unset casc_header_found
])dnl


dnl **********************************************************************
dnl * CASC_CREATE_PACKAGE_OPTION(PACKAGE-NAME[, DIR-LIST[, FILE]])
dnl * This is a general macro that creates a configure command-line option
dnl * called `--with-PACKAGE-NAME-dir' which will allow the user to
dnl * specify the location of the installation of an outside software
dnl * package, such as PETSc or ISIS++.  After a check to make sure the
dnl * given directory is valid (see below for discussion of validity), the
dnl * directory's path is stored in the shell variable PACKAGE-NAME_DIR.
dnl * For example, to allow the user to specify the location of PETSc,
dnl * place `CASC_CREATE_PACKAGE_OPTION(PETSC)' in configure.in.  Then the
dnl * user, if configuring on the CASC Sun cluster, would type `configure
dnl * --with-PETSC-dir=/home/casc/petsc', and the directory's path would
dnl * be stored in PETSC_DIR.  With this macro, the user is also permitted
dnl * to set the variable PACKAGE-NAME_DIR in the environment before
dnl * running configure, but any choice made on the command line would
dnl * override any preset values.  
dnl *
dnl * This macro takes an optional second argument, DIR-LIST, which is a
dnl * whitespace-separated list of directories where the developer thinks
dnl * PACKAGE-NAME might be installed.  If DIR-LIST is given, and the user
dnl * does not use the `--with' option to give the location of
dnl * PACKAGE-NAME (or if the directory given by the user does not exist),
dnl * then configure will assign to PACKAGE-NAME_DIR the path of the first
dnl * directory in DIR-LIST that is valid.
dnl *
dnl * Validity:  The optional third argument to this macro is FILE, which
dnl * should be either the name of a file in the top directory of the
dnl * package in question or the relative path of a file in a subdirectory
dnl * of the package.  If the argument FILE is given, then configure will
dnl * consider a user specified directory or a directory from DIR-LIST 
dnl * valid only if FILE exists in the directory.  If this argument is not
dnl * given, then configure will consider a directory valid simply if it
dnl * is indeed a directory.  FILE should be a file with a unique name
dnl * that can be expected to exist in the same location in any 
dnl * installation of the package in question.  If you know of no such
dnl * file, do not include a third argument when invoking this macro.
dnl * 
dnl * This macro also gives the user the command-line option
dnl * `--without-PACKAGE-NAME-dir', which, when invoked, will leave the
dnl * variable PACKAGE-NAME_DIR empty.  This option should be invoked when
dnl * the user wants to exclude a package from the configuration.
dnl * 
dnl * NOTE:  Since PACKAGE-NAME is used as part of both a command-line
dnl * option and a variable name, it MUST consist of only alphanumeric
dnl * characters.  PACKAGE-NAME is only a label, so it need not conform to
dnl * any existing directory or file name.  I would recommend that it be
dnl * all caps, as it becomes part of the name of a variable that is
dnl * substituted into the Makefile.
dnl **********************************************************************

AC_DEFUN([CASC_CREATE_PACKAGE_OPTION],
[
   AC_MSG_CHECKING([for $1 directory])

   dnl * $1 stands for the PACKAGE-NAME.  If [$1]_DIR has been set in the
   dnl * environment, give its value to casc_env_[$1]_dir, and clear
   dnl * [$1]_DIR.  The environmental value will ultimately be reassigned
   dnl * to [$1]_DIR if it is valid and no command-line options are able
   dnl * to change [$1]_DIR to a valid directory.  The environmental value
   dnl * will also be used even if it is invalid, if the command-line
   dnl * options and the DIRECTORY-LIST are both unable to generate a
   dnl * valid value.
   casc_result=
   casc_env_[$1]_dir=$[$1]_DIR
   [$1]_DIR=

   AC_ARG_WITH($1-dir, 
[  --with-$1-dir=DIR    $1 is installed in directory DIR
  --without-$1-dir     do not look for $1],

               if test "$withval" = "no"; then
                  casc_result="configuring without [$1]"
                  [$1]_DIR=
               fi
               , )

   dnl * If "--without-$1-dir" was given, then [$1]_DIR is left blank.
   dnl * Otherwise there is the following procedure to try to give
   dnl * [$1]_DIR a valid value:
   dnl *
   dnl * if "--with-$1-dir" was given
   dnl *    if the argument to "--with-$1-dir" is valid
   dnl *       assign the argument to [$1]_DIR
   dnl *    endif
   dnl * endif
   dnl *
   dnl * if a value for [$1]_DIR has not yet been found
   dnl *    if [$1]_DIR from the environment exists and is valid
   dnl *       assign the environmental value to [$1]_DIR
   dnl *    endif
   dnl * endif
   dnl *
   dnl * if [$1]_DIR still has no value
   dnl *    if the macro was given a DIRECTORY-LIST argument
   dnl *       for each directory in the list
   dnl *          if the directory is valid
   dnl *             assign the directory to [$1]_DIR
   dnl *             break loop
   dnl *          else
   dnl *             continue loop
   dnl *          endif
   dnl *       end loop
   dnl *       if [$1]_DIR still doesn't have a value
   dnl *          casc_result="none"
   dnl *       else
   dnl *          casc_result=$[$1]_DIR
   dnl *       endif
   dnl *    else
   dnl *       casc_result="none"
   dnl *    endif
   dnl * endif

   if test "$with_[$1]_dir" != "no"; then

      if test -d "$with_[$1]_dir"; then

         ifelse([$3], , ,
            if test -f $with_[$1]_dir/[$3]; then)

               casc_result="$with_[$1]_dir"
               [$1]_DIR="$casc_result"

         ifelse([$3], , ,
            fi)
      fi

      if test -z "$casc_result"; then

         if test -d "$casc_env_[$1]_dir"; then

            ifelse([$3], , ,
               if test -f $casc_env_[$1]_dir/[$3]; then)

                  casc_result="$casc_env_[$1]_dir"
                  [$1]_DIR="$casc_result"

            ifelse([$3], , ,
               fi)
         fi
      fi



      if test -z "$casc_result"; then
         [$1]_DIR=
   
         ifelse([$2], ,
            casc_result="none" ,

            for casc_dir in $2; do

               if test -d "$casc_dir"; then

                  ifelse([$3], , ,
                     if test -f $casc_dir/[$3]; then)

                        $1_DIR=$casc_dir

                        break

                  ifelse([$3], , ,
                     fi)
               fi
            done

            if test -z "$[$1]_DIR"; then
               casc_result="none"

            else
               casc_result="$[$1]_DIR"
            fi
         )
      fi
   fi

   dnl * $casc_result either is a valid value for [$1]_DIR or "none".
   dnl * if none, then assign the original environmental value of
   dnl * [$1]_DIR, whatever it may be, to casc_result and [$1]_DIR.  If
   dnl * there was no environmental value, then $casc_result remains
   dnl * "none" and [$1]_DIR is left empty.

   if test "$casc_result" = "none"; then

      if test -n "$casc_env_[$1]_dir"; then

         casc_result="$casc_env_[$1]_dir"
         [$1]_DIR="$casc_result"
      fi
   fi

   AC_MSG_RESULT($casc_result)
   AC_SUBST([$1]_DIR)
])


dnl smr_ARG_WITHLIB from FVWM by S. Robbins 
dnl Allow argument for optional libraries; wraps AC_ARG_WITH, to
dnl provide a "--with-foo-lib" option in the configure script, where foo
dnl is presumed to be a library name.  The argument given by the user
dnl (i.e. "bar" in ./configure --with-foo-lib=bar) may be one of four 
dnl things:
dnl     * boolean (no, yes or blank): whether to use library or not
dnl     * file: assumed to be the name of the library
dnl     * directory: assumed to *contain* the library
dnl     * a quoted, space-separated list of linker flags needed to link
dnl       with this library.  (To be used if this library requires
dnl       linker flags other than the normal `-L' and `-l' flags.)
dnl 
dnl The argument is sanity-checked.  If all is well, two variables are
dnl set: "with_foo" (value is yes, no, or maybe), and "foo_LIBFLAGS" (value
dnl is either blank, a file, -lfoo, '-L/some/dir -lfoo', or whatever 
dnl linker flags the user gives). The idea is: the first tells you whether
dnl the library is to be used or not (or the user didn't specify one way
dnl or the other) and the second to put on the command line for linking
dnl with the library.
dnl
dnl Usage:
dnl smr_ARG_WITHLIB(name, libname, description)
dnl 
dnl name                name for --with argument ("foo" for libfoo)
dnl libname             (optional) actual name of library,
dnl                     if different from name
dnl description         (optional) used to construct help string
dnl 
dnl Changes:  Changed some identifier names.
dnl           --with-foo-library is now --with-foo-lib
dnl           foo_LIBS is now foo_LIBFLAGS
dnl           Fourth posibility for argument to --with-foo-lib added
dnl           Documentation above changed to reflect these changes
dnl           Noah Elliott, October 1998


AC_DEFUN([CASC_SMR_ARG_WITHLIB],
[
   smr_ARG_WITHLIB([$1],[$2],[$3])
])dnl

AC_DEFUN([smr_ARG_WITHLIB], [

ifelse($2, , smr_lib=[$1], smr_lib=[$2]) 
    
AC_ARG_WITH([$1]-lib,
ifelse($3, ,
[  --with-$1-lib[=PATH]       use $1 library], 
[  --with-$1-lib[=PATH]       use $1 library ($3)]),
[
    if test "$withval" = yes; then
        with_[$1]=yes
        [$1]_LIBFLAGS="-l${smr_lib}"
    elif test "$withval" = no; then
        with_[$1]=no
        [$1]_LIBFLAGS=
    else
        with_[$1]=yes
        if test -f "$withval"; then
            [$1]_LIBFLAGS=$withval
        elif test -d "$withval"; then
            [$1]_LIBFLAGS="-L$withval -l${smr_lib}"
        else
            case $withval in
            -*)
               [$1]_LIBFLAGS="$withval"
            ;;
            *)
               AC_MSG_ERROR(
                  [argument must be boolean, file, directory, or compiler flags]
                           )
            ;;
            esac
        fi
    fi
], [
    with_[$1]=maybe
    [$1]_LIBFLAGS="-l${smr_lib}"
])])

    
dnl smr_ARG_WITHINCLUDES from FVWM by S. Robbins
dnl Check if the include files for a library are accessible, and
dnl define the variable "name_INCLUDE" with the proper "-I" flag for
dnl the compiler.  The user has a chance to specify the includes
dnl location, using "--with-foo-include".
dnl 
dnl This should be used *after* smr_ARG_WITHLIB *and* AC_CHECK_LIB are
dnl successful.
dnl 
dnl Usage:
dnl smr_ARG_WITHINCLUDES(name, header, extra-flags)
dnl 
dnl name                library name, MUST same as used with smr_ARG_WITHLIB
dnl header              a header file required for using the lib
dnl extra-flags         (optional) flags required when compiling the
dnl                     header, typically more includes; for ex. X_CFLAGS
dnl
dnl Changes:  Changed some identifier names.
dnl           --with-foo-includes is now --with-foo-include
dnl           name_CFLAGS is now name_INCLUDE
dnl           Documentation above changed to reflect these changes
dnl           Noah Elliott, October 1998

AC_DEFUN([CASC_SMR_ARG_WITHINCLUDES],
[
   smr_ARG_WITHINCLUDES([$1], [$2], [$3])
])dnl

AC_DEFUN([smr_ARG_WITHINCLUDES], [

AC_ARG_WITH([$1]-include,
[  --with-$1-include=DIR  set directory for $1 headers],
[
    if test -d "$withval"; then
        [$1]_INCLUDE="-I${withval}"
    else
        AC_MSG_ERROR(argument must be a directory)
    fi])

dnl This bit of logic comes from autoconf's AC_PROG_CC macro.  We need
dnl to put the given include directory into CPPFLAGS temporarily, but
dnl then restore CPPFLAGS to its old value.
dnl 
smr_test_CPPFLAGS="${CPPFLAGS+set}"
smr_save_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS ${[$1]_CFLAGS}"

    ifelse($3, , , CPPFLAGS="$CPPFLAGS [$3]")
    AC_CHECK_HEADERS($2)
   
if test "$smr_test_CPPFLAGS" = set; then
    CPPFLAGS=$smr_save_CPPFLAGS
else
    unset CPPFLAGS
fi
])
    
        
dnl smr_CHECK_LIB from FVWM by S. Robbins
dnl Probe for an optional library.  This macro creates both
dnl --with-foo-lib and --with-foo-include options for the configure
dnl script.  If --with-foo-lib is *not* specified, the default is to
dnl probe for the library, and use it if found.
dnl
dnl Usage:
dnl smr_CHECK_LIB(name, libname, desc, func, header, x-libs, x-flags)
dnl 
dnl name        name for --with options
dnl libname     (optional) real name of library, if different from
dnl             above
dnl desc        (optional) short descr. of library, for help string
dnl func        function of library, to probe for
dnl header      (optional) header required for using library
dnl x-libs      (optional) extra libraries, if needed to link with lib
dnl x-flags     (optional) extra flags, if needed to include header files
dnl
dnl Changes:  identifier names and documentation modified to reflect
dnl           changes to smr_ARG_WITHLIB and smr_ARG_WITHINCLUDES
dnl           Noah Elliott, October 1998

AC_DEFUN([CASC_SMR_CHECK_LIB],
[
   smr_CHECK_LIB([$1], [$2], [$3], [$4], [$5], [$6], [$7])
])dnl

AC_DEFUN([smr_CHECK_LIB],
[   
ifelse($2, , smr_lib=[$1], smr_lib=[$2])
ifelse($5, , , smr_header=[$5])
smr_ARG_WITHLIB($1,$2,$3)
if test "$with_$1" != no; then
    AC_CHECK_LIB($smr_lib, $4,
        smr_havelib=yes, smr_havelib=no,
        ifelse($6, , ${$1_LIBFLAGS}, [${$1_LIBFLAGS} $6]))
    if test "$smr_havelib" = yes -a "$smr_header" != ""; then
        smr_ARG_WITHINCLUDES($1, $smr_header, $7)
        smr_safe=`echo "$smr_header" | sed 'y%./+-%__p_%'`
        if eval "test \"`echo '$ac_cv_header_'$smr_safe`\" != yes"; then
            smr_havelib=no
        fi
    fi
    if test "$smr_havelib" = yes; then
        AC_MSG_RESULT(Using $1 library)
    else
        $1_LIBFLAGS=
        $1_INCLUDE=
        if test "$with_$1" = maybe; then
            AC_MSG_RESULT(Not using $1 library)
        else
            AC_MSG_WARN(Requested $1 library not found!)
        fi
    fi
fi])


dnl **********************************************************************
dnl * CASC_CONFIG_OUTPUT_LIST(DIR-LIST[, OUTPUT-FILE])
dnl *
dnl * The intent of this macro is to make configure handle the possibility
dnl * that a portion of the directory tree of a project may not be
dnl * present.  This will modify the argument list of AC_OUTPUT to contain
dnl * only output file names for which corresponding input files exist.
dnl * If you are not concerned about the possible absence of the necessary
dnl * input (.in) files, it is better to not use this macro and to
dnl * explicitly list all of the output files in a call to AC_OUTPUT.
dnl * Also, If you wish to create a file Foo from a file with a name
dnl * other than Foo.in, this macro will not work, and you must use
dnl * AC_OUTPUT.
dnl *
dnl * This macro checks for the existence of the file OUTPUT-FILE.in in
dnl * each directory specified in the whitespace-separated DIR-LIST.  
dnl * (Directories should be specified by relative path from the directory 
dnl * containing configure.in.) If OUTPUT-FILE is not specified, the
dnl * default is 'Makefile'.  For each directory that contains 
dnl * OUTPUT-FILE.in, the relative path of OUTPUT-FILE is added to the 
dnl * shell variable OUTPUT-FILE_list.  When AC_OUTPUT is called,
dnl * '$OUTPUT-FILE_list' should be included in the argument list.  So if
dnl * you have a directory tree and each subdirectory contains a 
dnl * Makefile.in, DIR-LIST should be a list of every subdirectory and
dnl * OUTPUT-FILE can be omitted, because 'Makefile' is the default.  When
dnl * configure runs, it will check for the existence of a Makefile.in in
dnl * each directory in DIR-LIST, and if so, the relative path of each
dnl * intended Makefile will be added to the variable Makefile_list.
dnl *
dnl * This macro can be called multiple times, if there are files other
dnl * than Makefile.in with a .in suffix other that are intended to be 
dnl * processed by configure. 
dnl *
dnl * Example
dnl *     If directories dir1 and dir2 both contain a file named Foo.in, 
dnl *     and you wish to use configure to create a file named Foo in each
dnl *     directory, then call 
dnl *     CASC_CONFIG_OUTPUT_LIST(dir1 dir2, Foo)
dnl *     If you also called this macro for Makefile as described above,
dnl *     you should call
dnl *     AC_OUTPUT($Makefile_list $Foo_list)
dnl *     at the end of configure.in .
dnl *********************************************************************


AC_DEFUN([CASC_CONFIG_OUTPUT_LIST],
[
   dnl * m_OUTPUT_LIST is a macro to store the name of the variable
   dnl * which will contain the list of output files
   define([m_OUTPUT_LIST], ifelse([$2], , Makefile_list, [$2_list]))

   if test -z "$srcdir"; then
      srcdir=.
   fi

   dnl * use "Makefile" if second argument not given
   if test -n "$2"; then
      casc_output_file=$2
   else   
      casc_output_file=Makefile
   fi   
      
   dnl * Add a file to the output list if its ".in" file exists.
   for casc_dir in $1; do
      if test -f $srcdir/$casc_dir/$casc_output_file.in; then
         m_OUTPUT_LIST="$m_OUTPUT_LIST $casc_dir/$casc_output_file"
      fi
   done
])dnl


dnl **********************************************************************
dnl * CASC_GUESS_ARCH
dnl * Guesses a one-word name for the current architecture, unless ARCH
dnl * has been preset.  This is an alternative to the built-in macro
dnl * AC_CANONICAL_HOST, which gives a three-word name.  Uses the utility
dnl * 'tarch', which is a Bourne shell script that should be in the same  
dnl * directory as the configure script.  If tarch is not present or if it
dnl * fails, ARCH is set to the value, if any, of shell variable HOSTTYPE,
dnl * otherwise ARCH is set to "unknown".
dnl **********************************************************************

AC_DEFUN([CASC_GUESS_ARCH],
[
   AC_MSG_CHECKING(the architecture)

   dnl * $ARCH could already be set in the environment or earlier in configure
   dnl * Use the preset value if it exists, otherwise go throug the procedure
   if test -z "$ARCH"; then

      dnl * configure searches for the tool "tarch".  It should be in the
      dnl * same directory as configure.in, but a couple of other places
      dnl * will be checked.  casc_tarch stores a relative path for "tarch".
      casc_tarch_dir=
      for casc_dir in $srcdir $srcdir/.. $srcdir/../.. $srcdir/config; do
         if test -f $casc_dir/tarch; then
            casc_tarch_dir=$casc_dir
            casc_tarch=$casc_tarch_dir/tarch
            break
         fi
      done

      dnl * if tarch was not found or doesn't work, try using env variable
      dnl * $HOSTTYPE
      if test -z "$casc_tarch_dir"; then
         AC_MSG_WARN(cannot find tarch, using \$HOSTTYPE as the architecture)
         ARCH=$HOSTTYPE
      else
         ARCH="`$casc_tarch`"

         if test -z "$ARCH" || test "$ARCH" = "unknown"; then
            ARCH=$HOSTTYPE
         fi
      fi

      dnl * if $ARCH is still empty, give it the value "unknown".
      if test -z "$ARCH"; then
         ARCH=unknown
         AC_MSG_WARN(architecture is unknown)
      else
         AC_MSG_RESULT($ARCH)
      fi    
   else
      AC_MSG_RESULT($ARCH)
   fi

   AC_SUBST(ARCH)

])dnl


dnl **********************************************************************
dnl * CASC_SET_SUFFIX_RULES is not like the other macros in aclocal.m4
dnl * because it does not run any kind of test on the system on which it
dnl * is running.  All it does is create several variables which contain
dnl * the text of some simple implicit suffix rules that can be
dnl * substituted into Makefile.in.  The suffix rules that come from the
dnl * macro all deal with compiling a source file into an object file.  If
dnl * this macro is called in configure.in, then if `@CRULE@' is placed in
dnl * Makefile.in, the following will appear in the generated Makefile:
dnl *
dnl * .c.o:
dnl *         @echo "Making (c) " $@ 
dnl *         @${CC} -o $@ -c ${CFLAGS} $<	
dnl *
dnl * The following is a list of the variables created by this macro and
dnl * the corresponding suffixes of the files that each implicit rule 
dnl * deals with.
dnl *
dnl * CRULE       --   .c
dnl * CXXRULE     --   .cxx
dnl * CPPRULE     --   .cpp
dnl * CCRULE      --   .cc
dnl * CAPCRULE    --   .C
dnl * F77RULE     --   .f
dnl *
dnl * There are four suffix rules for C++ files because of the different
dnl * suffixes that can be used for C++.  Only use the one which
dnl * corresponds to the suffix you use for your C++ files.
dnl *
dnl * The rules created by this macro require you to use the following
dnl * conventions for Makefile variables:
dnl *
dnl * CC        = C compiler
dnl * CXX       = C++ compiler
dnl * F77       = Fortran 77 compiler
dnl * CFLAGS    = C compiler flags
dnl * CXXFLAGS  = C++ compiler flags
dnl * FFLAGS    = Fortran 77 compiler flags
dnl **********************************************************************

AC_DEFUN([CASC_SET_SUFFIX_RULES],
[
   dnl * Things weren't working whenever "$@" showed up in the script, so
   dnl * I made the symbol $at_sign to signify '@'
   at_sign=@

   dnl * All of the backslashes are used to handle the $'s and the
   dnl * newlines which get passed through echo and sed.

   CRULE=`echo ".c.o:\\\\
\t@echo \"Making (c) \" \\$$at_sign \\\\
\t@\\${CC} -o \\$$at_sign -c \\${CFLAGS} \$<"`

   AC_SUBST(CRULE)

   CXXRULE=`echo ".cxx.o:\\\\
\t@echo \"Making (c++) \" \\$$at_sign \\\\
\t@\\${CXX} -o \\$$at_sign -c \\${CXXFLAGS} \$<"`

   AC_SUBST(CXXRULE)

   CPPRULE=`echo ".cpp.o:\\\\
\t@echo \"Making (c++) \" \\$$at_sign \\\\
\t@\\${CXX} -o \\$$at_sign -c \\${CXXFLAGS} \$<"`

   AC_SUBST(CPPRULE)

   CCRULE=`echo ".cc.o:\\\\
\t@echo \"Making (c++) \" \\$$at_sign \\\\
\t@\\${CXX} -o \\$$at_sign -c \\${CXXFLAGS} \$<"`

   AC_SUBST(CCRULE)

   CAPCRULE=`echo ".C.o:\\\\
\t@echo \"Making (c++) \" \\$$at_sign \\\\
\t@\\${CXX} -o \\$$at_sign -c \\${CXXFLAGS} \$<"`

   AC_SUBST(CAPCRULE)

   F77RULE=`echo ".f.o:\\\\
\t@echo \"Making (f) \" \\$$at_sign \\\\
\t@\\${F77} -o \\$$at_sign -c \\${FFLAGS} \$<"`

   AC_SUBST(F77RULE)

])

dnl Macro to save compiler state flags for invoking dnl compiler tests
dnl NOTE that this is NOT currently a stack so can dnl only be called
dnl in push/pop order.  push push pop pop dnl will fail
AC_DEFUN([CASC_PUSH_COMPILER_STATE],[
   casc_save_LIBS=$LIBS
   casc_save_CXXFLAGS=$CXXFLAGS
])

dnl Macro to restore compiler state flags for invoking
dnl compiler tests
AC_DEFUN([CASC_POP_COMPILER_STATE],[
   LIBS=$casc_save_LIBS
   unset casc_save_LIBS
   CXXFLAGS=$casc_save_CXXFLAGS
   unset casc_save_CXXFLAGS
])

dnl ********************************************************************
dnl * CASC_PROG_MPICC searches the PATH for an available MPI C compiler
dnl * wraparound.  It assigns the name to MPICC.
dnl ********************************************************************

AC_DEFUN([CASC_PROG_MPICC],
[
   AC_CHECK_PROGS(MPICC, mpcc mpicc tmcc hcc)
   test -z "$MPICC" && AC_MSG_ERROR([no acceptable mpicc found in \$PATH])
])dnl


dnl ********************************************************************
dnl * CASC_PROG_MPICXX searches the PATH for an available MPI C++
dnl * compiler wraparound.  It assigns the name to MPICXX.
dnl ********************************************************************

AC_DEFUN([CASC_PROG_MPICXX],
[
   AC_CHECK_PROGS(MPICXX, mpKCC mpCC mpig++ mpiCC hcp)
   test -z "$MPICXX" && AC_MSG_ERROR([no acceptable mpic++ found in \$PATH])
])dnl


dnl **********************************************************************
dnl * CASC_PROG_MPIF77 searches the PATH for an available MPI Fortran 77
dnl * compiler wraparound.  It assigns the name to MPIF77.
dnl **********************************************************************

AC_DEFUN([CASC_PROG_MPIF77],
[
   AC_CHECK_PROGS(MPIF77, mpf77 mpxlf mpif77 mpixlf tmf77 hf77)
   test -z "$MPIF77" && AC_MSG_ERROR([no acceptable mpif77 found in \$PATH])
])dnl


dnl ***********************************************************************
dnl * CASC_CHECK_MPIF77_PP checks whether the preprocessor needs to
dnl * be called before calling the compiler for Fortran files with
dnl * preprocessor directives and MPI function calls.  If the preprocessor
dnl * is necessary, MPIF77NEEDSPP is set to "yes", otherwise it is set to
dnl * "no"
dnl ***********************************************************************

AC_DEFUN([CASC_CHECK_MPIF77_PP],
[
   AC_REQUIRE([CASC_PROG_MPIF77])

   rm -f testppmp.o

   AC_MSG_CHECKING(whether $FPP needs to be called before $MPIF77)

   # This follows the same procedur as CASC_CHECK_F77_PP, except it tests
   # $MPIF77 using a test program that includes MPI functions.

   cat > testppmp.F <<EOF
#define FOO 3
	program testppmp
	include 'mpif.h'
	integer rank,size,mpierr,sum
	call MPI_INIT(mpierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,size,mpierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD,rank,mpierr)
#ifdef FORTRAN_NO_UNDERSCORE
        sum = rank + size
#else
        sum = rank + rank
#endif
        call MPI_FINALIZE(mpierr)
        end
EOF

   $MPIF77 -DBAR -c testppmp.F
   if test -f testppmp.o; then
      MPIF77NEEDSPP=no
   else
      MPIF77NEEDSPP=yes
   fi

   echo $MPIF77NEEDSPP
   rm -f testppmp.o testppmp.F
   AC_SUBST(MPIF77NEEDSPP)
])dnl


dnl *********************************************************************
dnl * CASC_SET_MPI sets up the needed MPI library and directory flags.
dnl * The location of the file mpi.h is put into the variable MPIINCLUDE
dnl * as a -I flag.  The -l flags that specify the needed libraries and
dnl * the -L flags that specify the paths of those libraries are placed in
dnl * the variables MPILIBS and MPILIBDIRS, respectively.  To set the MPI
dnl * libraries and directories manually, use the --with-mpi-include,
dnl * --with-mpi-libs, and --with-mpi-lib-dirs command-line options when
dnl * invoking configure.  Only one directory should be specified with
dnl * --with-mpi-include, while any number of directories can be specified
dnl * by --with-mpi-lib-dirs.  Any number of libraries can be specified
dnl * with --with-mpi-libs, and the libraries must be referred to by their
dnl * base names, so libmpi.a is just mpi.  It is adviseable to use all
dnl * three --with flags whenever one is used, because it is likely that
dnl * when one is chosen it will mess up the automatic choices for the
dnl * other two.  If the architecture is unknown, or if the needed MPI
dnl * settings for the current architecture are not known, then the naive
dnl * settings of MPILIBS="-lmpi" and MPILIBDIRS="-L/usr/local/mpi/lib"
dnl * are tested, and if they exist they are used, otherwise the MPILIB*
dnl * variables are left blank.  In the case of rs6000, the variable
dnl * MPIFLAGS is also set.
dnl **********************************************************************

AC_DEFUN([CASC_SET_MPI],
        [

   dnl * If called from within CASC_FIND_MPI, then the configure-line
   dnl * options will already exist.  This ifdef creates them otherwise.
   ifdef([AC_PROVIDE_CASC_FIND_MPI], ,
      [AC_ARG_WITH(mpi-include, [  --with-mpi-include=DIR  mpi.h is in DIR],
                  casc_mpi_include_dir=$withval)

      AC_ARG_WITH(mpi-libs,
[  --with-mpi-libs=LIBS    LIBS is space-separated list of library names
                          needed for MPI, e.g. \"nsl socket mpi\"],
                  casc_mpi_libs=$withval)

      AC_ARG_WITH(mpi-lib-dirs,
[  --with-mpi-lib-dirs=DIRS
                          DIRS is space-separated list of directories
                          containing the libraries specified by
                          \`--with-mpi-libs', e.g \"/usr/lib /usr/local/mpi/lib\"],
                  casc_mpi_lib_dirs=$withval)]
   )

   if test -z "$casc_mpi_libs"; then
      AC_REQUIRE([CASC_GUESS_ARCH])

      dnl * Set everything to known values
      case $ARCH in

         sun4 | solaris)
            case $F77 in
               *g77)
                   if test -z "$casc_mpi_include_dir"; then
                      casc_mpi_include_dir=/usr/local/mpi/lam/h
                   fi

                   if test -z "$casc_mpi_lib_dirs"; then
                      casc_mpi_lib_dirs="/usr/local/mpi/lam/lib"
                   fi

                   casc_mpi_libs="socket mpi trillium args tstdio t";;

               *)

                  if test -z "$casc_mpi_include_dir"; then
                     MPIINCLUDE="-I/usr/local/mpi/mpich/include \
                                 -I/usr/local/mpi/mpich/lib/solaris/ch_p4"
                  fi

                  if test -z "$casc_mpi_lib_dirs"; then
                     casc_mpi_lib_dirs="/usr/local/mpi/mpich/lib/solaris/ch_p4 \
                                       /usr/lib"
                  fi

               casc_mpi_libs="nsl socket mpi";;
               esac

            if test -z "$MPIINCLUDE"; then
               AC_CHECK_HEADER($casc_mpi_include_dir/mpi.h,
                               MPIINCLUDE="-I$casc_mpi_include_dir")
            fi
         ;;

         alpha)
            if test -z "$casc_mpi_include_dir"; then
               casc_mpi_include_dir=/usr/local/mpi/include
            fi
            AC_CHECK_HEADER($casc_mpi_include_dir/mpi.h,
                               MPIINCLUDE="-I$casc_mpi_include_dir")

            if test -z "$casc_mpi_lib_dirs"; then
               casc_mpi_lib_dirs="/usr/local/mpi/lib/alpha/ch_shmem \
                                  /usr/local/lib"
            fi

            casc_mpi_libs="mpich gs";;

         rs6000)

dnl            if test -z "$casc_mpi_include_dir"; then
dnl               casc_mpi_include_dir=/usr/lpp/ppe.poe/include
dnl            fi
dnl            AC_CHECK_HEADER($casc_mpi_include_dir/mpi.h,
dnl                               MPIINCLUDE="-I$casc_mpi_include_dir")

dnl            if test -z "$casc_mpi_lib_dirs"; then
dnl               casc_mpi_lib_dirs=/usr/lpp/ppe.poe/lib
dnl            fi

            casc_mpi_libs=mpi

            MPIFLAGS="-binitfini:poe_remote_main";;

         IRIX64 | iris4d)
            if test -z "$casc_mpi_include_dir"; then
               casc_mpi_include_dir=/usr/local/mpi/include
            fi
            AC_CHECK_HEADER($casc_mpi_include_dir/mpi.h,
                               MPIINCLUDE="-I$casc_mpi_include_dir")

            if test -z "$casc_mpi_lib_dirs"; then
               casc_mpi_lib_dirs=/usr/local/mpi/lib/IRIX64/ch_p4
            fi

            casc_mpi_libs=mpi;;

         *)
AC_MSG_WARN([trying naive MPI settings - can use --with flags to change])
            if test -z "$casc_mpi_include_dir"; then
               casc_mpi_include_dir=/usr/local/mpi/include
            fi
            AC_CHECK_HEADER($casc_mpi_include_dir/mpi.h,
                               MPIINCLUDE="-I$casc_mpi_include_dir")

            if test -z "$casc_mpi_lib_dirs"; then
               casc_mpi_lib_dirs=/usr/local/mpi/lib
            fi
            casc_mpi_libs=mpi ;;
      esac

      for casc_lib in $casc_mpi_libs; do
         CASC_ADD_LIB($casc_lib, main, $casc_mpi_lib_dirs, MPI)
      done

   else
      if test -n "$casc_mpi_include_dir"; then
         MPIINCLUDE="-I$casc_mpi_include_dir"
      else
         MPIINCLUDE=
      fi

      if test -n "$casc_mpi_lib_dirs"; then
         for casc_lib_dir in $casc_mpi_lib_dirs; do
            MPILIBDIRS="-L$casc_lib_dir $MPILIBDIRS"
         done
      else
         MPILIBDIRS=
      fi

      for casc_lib in $casc_mpi_libs; do
         MPILIBS="$MPILIBS -l$casc_lib"
      done
   fi
])dnl


dnl ********************************************************************
dnl * CASC_FIND_MPI will determine the libraries, directories, and other
dnl * flags needed to compile and link programs with MPI function calls.
dnl * This macro runs tests on the script found by the CASC_PROG_MPICC
dnl * macro.  If there is no such mpicc-type script in the PATH and
dnl * MPICC is not set manually, then this macro will not work.
dnl *
dnl * One may question why these settings would need to be determined if
dnl * there already is mpicc available, and that is a valid question.  I
dnl * can think of a couple of reasons one may want to use these settings
dnl * rather than using mpicc directly.  First, these settings allow you
dnl * to choose the C compiler you wish to use rather than using whatever
dnl * compiler is written into mpicc.  Also, the settings determined by
dnl * this macro should also work with C++ and Fortran compilers, so you
dnl * won't need to have mpiCC and mpif77 alongside mpicc.  This is
dnl * especially helpful on systems that don't have mpiCC.  The advantage
dnl * of this macro over CASC_SET_MPI is that this one doesn't require
dnl * a test of the machine type and thus will hopefully work on unknown
dnl * architectures.  The main disadvantage is that it relies on mpicc.
dnl *
dnl * --with-mpi-include, --with-mpi-libs, and --with-mpi-lib-dirs can be
dnl * used to manually override the automatic test, just as with
dnl * CASC_SET_MPI.  If any one of these three options are used, the
dnl * automatic test will not be run, so it is best to call all three
dnl * whenever one is called.  In addition, the option --with-mpi-flags is
dnl * available here to set any other flags that may be needed, but it
dnl * does not override the automatic test.  Flags set by --with-mpi-flags
dnl * will be added to the variable MPIFLAGS.  This way, if the macro, for
dnl * whatever reason, leaves off a necessary flag, the flag can be added
dnl * to MPIFLAGS without eliminating anything else.  The other variables
dnl * set are MPIINCLUDE, MPILIBS, and MPILIBDIRS, just as in
dnl * CASC_SET_MPI.  This macro also incorporates CASC_SET_MPI as a backup
dnl * plan, where if there is no mpicc, it will use the settings
dnl * determined by architecture name in CASC_SET_MPI
dnl ********************************************************************

AC_DEFUN([CASC_FIND_MPI],
[

   casc_find_mpi_cache_used=yes

   AC_CACHE_VAL(casc_cv_mpi_include, casc_find_mpi_cache_used=no)
   AC_CACHE_VAL(casc_cv_mpi_libs, casc_find_mpi_cache_used=no)
   AC_CACHE_VAL(casc_cv_mpi_lib_dirs, casc_find_mpi_cache_used=no)
   AC_CACHE_VAL(casc_cv_mpi_flags, casc_find_mpi_cache_used=no)

   if test "$casc_find_mpi_cache_used" = "yes"; then
      AC_MSG_CHECKING(for location of mpi.h)
      MPIINCLUDE=$casc_cv_mpi_include
      AC_MSG_RESULT("\(cached\) $MPIINCLUDE")

      AC_MSG_CHECKING(for MPI library directories)
      MPILIBDIRS=$casc_cv_mpi_lib_dirs
      AC_MSG_RESULT("\(cached\) $MPILIBDIRS")

      AC_MSG_CHECKING(for MPI libraries)
      MPILIBS=$casc_cv_mpi_libs
      AC_MSG_RESULT("\(cached\) $MPILIBS")

      AC_MSG_CHECKING(for other MPI-related flags)
      MPIFLAGS=$casc_cv_mpi_flags
      AC_MSG_RESULT("\(cached\) $MPIFLAGS")
   else


      dnl * Set up user options.  If user uses any of the fist three options,
      dnl * then automatic tests are not run.

      casc_user_chose_mpi=no
      AC_ARG_WITH(mpi-include, [  --with-mpi-include=DIR  mpi.h is in DIR],
                  for mpi_dir in $withval; do
                     MPIINCLUDE="$MPIINCLUDE -I$withval"
                  done; casc_user_chose_mpi=yes)

      AC_ARG_WITH(mpi-libs,
[  --with-mpi-libs=LIBS    LIBS is space-separated list of library names
                          needed for MPI, e.g. \"nsl socket mpi\"],
                  for mpi_lib in $withval; do
                     MPILIBS="$MPILIBS -l$mpi_lib"
                  done; casc_user_chose_mpi=yes)


      AC_ARG_WITH(mpi-lib-dirs,
[  --with-mpi-lib-dirs=DIRS
                          DIRS is space-separated list of directories
                          containing the libraries specified by
                          \`--with-mpi-libs', e.g \"/usr/lib /usr/local/mpi/lib\"],
                  for mpi_lib_dir in $withval; do
                     MPILIBDIRS="-L$mpi_lib_dir $MPILIBDIRS"
                  done; casc_user_chose_mpi=yes)

      dnl * --with-mpi-flags only adds to automatic selections,
      dnl * does not override

      AC_ARG_WITH(mpi-flags,
[  --with-mpi-flags=FLAGS  FLAGS is space-separated list of whatever flags other
                          than -l and -L are needed to link with mpi libraries],
                          MPIFLAGS=$withval)


      if test "$casc_user_chose_mpi" = "no"; then

      dnl * Find an MPICC.  If there is none, call CASC_SET_MPI to choose MPI
      dnl * settings based on architecture name.  If CASC_SET_MPI fails,
      dnl * print warning message.  Manual MPI settings must be used.

         AC_ARG_WITH(MPICC,
[  --with-MPICC=ARG        ARG is mpicc or similar MPI C compiling tool],
            MPICC=$withval,
            [AC_CHECK_PROGS(MPICC, mpcc mpicc tmcc hcc)])

         if test -z "$MPICC"; then
            AC_MSG_WARN([no acceptable mpicc found in \$PATH])
            CASC_SET_MPI
            if test -z "$MPILIBS"; then
             AC_MSG_WARN([MPI not found - must set manually using --with flags])
            fi

         dnl * When $MPICC is there, run the automatic test
         dnl * here begins the hairy stuff

         else

            dnl changequote(, )dnl

            dnl * Create a minimal MPI program.  It will be compiled using
            dnl * $MPICC with verbose output.
            cat > mpconftest.c << EOF
#include "mpi.h"

main(int argc, char **argv)
{
   int rank, size;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Finalize();
   return 0;
}
EOF

            casc_mplibs=
            casc_mplibdirs=
            casc_flags=
            casc_lmpi_exists=no

            dnl * These are various ways to produce verbose output from $MPICC
            dnl * All of their outputs are stuffed into variable
            dnl * $casc_mpoutput

            for casc_command in "$MPICC -show"\
                                "$MPICC -v"\
                                "$MPICC -#"\
                                "$MPICC"; do

               casc_this_output=`$casc_command mpconftest.c -o mpconftest 2>&1`

               dnl * If $MPICC uses xlc, then commas must be removed from output
               xlc_p=`echo $casc_this_output | grep xlcentry`
               if test -n "$xlc_p"; then
                  casc_this_output=`echo $casc_this_output | sed 's/,/ /g'`
               fi

               dnl * Turn on flag once -lmpi is found in output
               lmpi_p=`echo $casc_this_output | grep "\-lmpi"`
               if test -n "$lmpi_p"; then
                  casc_lmpi_exists=yes
               fi

               casc_mpoutput="$casc_mpoutput $casc_this_output"
               casc_this_output=

            done

            rm -rf mpconftest*

            dnl * little test to identify $CC as IBM's xlc
            echo "main() {}" > cc_conftest.c
            cc_output=`${CC-cc} -v -o cc_conftest cc_conftest.c 2>&1`
            xlc_p=`echo $cc_output | grep xlcentry`
            if test -n "$xlc_p"; then
               casc_compiler_is_xlc=yes
            fi
            rm -rf cc_conftest*

            dnl * $MPICC might not produce '-lmpi', but we still need it.
            dnl * Add -lmpi to $casc_mplibs if it was never found
            if test "$casc_lmpi_exists" = "no"; then
               casc_mplibs="-lmpi"
            else
               casc_mplibs=
            fi

            casc_want_arg=

            dnl * Loop through every word in output to find possible flags.
            dnl * If the word is the absolute path of a library, it is added
            dnl * to $casc_flags.  Any "-llib", "-L/dir", "-R/dir" and
            dnl * "-I/dir" is kept.  If '-l', '-L', '-R', '-I', '-u', or '-Y'
            dnl * appears alone, then the next word is checked.  If the next
            dnl * word is another flag beginning with '-', then the first
            dnl * word is discarded.  If the next word is anything else, then
            dnl * the two words are coupled in the $casc_arg variable.
            dnl * "-binitfini:poe_remote_main" is a flag needed especially
            dnl * for IBM MPI, and it is always kept if it is found.
            dnl * Any other word is discarded.  Also, after a word is found
            dnl * and kept once, it is discarded if it appears again

            for casc_arg in $casc_mpoutput; do

               casc_old_want_arg=$casc_want_arg
               casc_want_arg=

               if test -n "$casc_old_want_arg"; then
                  case "$casc_arg" in
                  [-*)]
                     casc_old_want_arg=
                  ;;
                  esac
               fi

               case "$casc_old_want_arg" in
               ['')]
                  case $casc_arg in
                  [/*.a)]
                     exists=false
                     for f in $casc_flags; do
                        if test x$casc_arg = x$f; then
                           exists=true
                        fi
                     done
                     if $exists; then
                        casc_arg=
                     else
                        casc_flags="$casc_flags $casc_arg"
                     fi
                  ;;
                  [-binitfini:poe_remote_main)]
                     exists=false
                     for f in $casc_flags; do
                        if test x$casc_arg = x$f; then
                           exists=true
                        fi
                     done
                     if $exists; then
                        casc_arg=
                     else
                        casc_flags="$casc_flags $casc_arg"
                     fi
                  ;;
                  [-lang*)]
                     casc_arg=
                  ;;
                  [-[lLR])]
                     casc_want_arg=$casc_arg
                     casc_arg=
                  ;;
                  [-[lLR]*)]
                     exists=false
                     for f in $casc_flags; do
                        if test x$casc_arg = x$f; then
                           exists=true
                        fi
                     done
                     if $exists; then
                        casc_arg=
                     else
                       casc_flags="$casc_flags $casc_arg"
                     fi
                  ;;
                 [-u)]
                     casc_want_arg=$casc_arg
                     casc_arg=
                  ;;
                  [-Y)]
                     casc_want_arg=$casc_arg
                     casc_arg=
                  ;;
                  [-I)]
                     casc_want_arg=$casc_arg
                     casc_arg=
                  ;;
                  [-I*)]
                     exists=false
                     for f in $casc_flags; do
                        if test x$casc_arg = x$f; then
                           exists=true
                        fi
                     done
                     if $exists; then
                        casc_arg=
                     else
                        casc_flags="$casc_flags $casc_arg"
                     fi
                  ;;
                  [*)]
                     casc_arg=
                  ;;
                  esac

               ;;
               [-[lLRI])]
                  casc_arg="casc_old_want_arg $casc_arg"
               ;;
               [-u)]
                  casc_arg="-u $casc_arg"
               ;;
               [-Y)]
                  casc_arg=`echo $casc_arg | sed -e 's%^P,%%'`
                  SAVE_IFS=$IFS
                  IFS=:
                  casc_list=
                  for casc_elt in $casc_arg; do
                     casc_list="$casc_list -L$casc_elt"
                  done
                  IFS=$SAVE_IFS
                  casc_arg="$casc_list"
               ;;
               esac

               dnl * Still inside the big for loop, we separate each flag
               dnl * into includes, libdirs, libs, flags
               if test -n "$casc_arg"; then
                  case $casc_arg in
                  [-I*)]

                     dnl * if the directory given in this flag contains mpi.h
                     dnl * then the flag is assigned to $MPIINCLUDE
                     if test -z "$MPIINCLUDE"; then
                        casc_cppflags="$casc_cppflags $casc_arg"
                        casc_include_dir=`echo "$casc_arg" | sed 's/-I//g'`

                        SAVE_CPPFLAGS="$CPPFLAGS"
                        CPPFLAGS="$casc_cppflags"
                        dnl changequote([, ])dnl

                        unset ac_cv_header_mpi_h
                        AC_CHECK_HEADER(mpi.h,
                                        MPIINCLUDE="$casc_cppflags")

                        dnl changequote(, )dnl
                        CPPFLAGS="$SAVE_CPPFLAGS"

                     else
                        casc_arg=
                     fi
                  ;;
                  [-[LR]*)]

                     dnl * These are the lib directory flags
                     casc_mplibdirs="$casc_mplibdirs $casc_arg"
                  ;;
                  [-l* | /*)]

                     dnl * These are the libraries
                     casc_mplibs="$casc_mplibs $casc_arg"
                  ;;
                  [-binitfini:poe_remote_main)]
                     if test "$casc_compiler_is_xlc" = "yes"; then
                        casc_mpflags="$casc_mpflags $casc_arg"
                     fi
                  ;;
                  [*)]
                     dnl * any other flag that has been kept goes here
                     casc_mpflags="$casc_mpflags $casc_arg"
                  ;;
                  esac

                  dnl * Upcoming test needs $LIBS to contain the flags
                  dnl * we've found
                  LIBS_SAVE=$LIBS
                  LIBS="$MPIINCLUDE $casc_mpflags $casc_mplibdirs $casc_mplibs"

                  if test -n "`echo $LIBS | grep '\-R/'`"; then
                     LIBS=`echo $LIBS | sed 's/-R\//-R \//'`
                  fi

                  dnl changequote([, ])dnl


                  dnl * Test to see if flags found up to this point are
                  dnl * sufficient to compile and link test program.  If not,
                  dnl * the loop keeps going to the next word
                  AC_LANG_PUSH(C)
                  AC_TRY_LINK(
dnl                     ifelse(AC_LANG, [C++],

dnl [#ifdef __cplusplus
dnl extern "C"
dnl #endif
dnl ])dnl
[#include "mpi.h"
], [int rank, size;
   int argc;
   char **argv;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Finalize();
],
                     casc_result=yes)
                  AC_LANG_POP(C)
                  LIBS=$LIBS_SAVE

                  if test "$casc_result" = yes; then
                     casc_result=
                     break
                  fi
               fi
            done

            dnl * After loop is done, set variables to be substituted
            MPILIBS=$casc_mplibs
            MPILIBDIRS=$casc_mplibdirs
            MPIFLAGS="$MPIFLAGS $casc_mpflags"

            dnl * IBM MPI uses /usr/lpp/ppe.poe/libc.a instead of /lib/libc.a
            dnl * so we need to make sure that -L/lib is not part of the
            dnl * linking line when we use IBM MPI.  This only appears in
            dnl * configure when CASC_FIND_MPI is called first.
	    dnl            ifdef([AC_PROVIDE_CASC_FIND_F77LIBS],
            dnl               if test -n "`echo $F77LIBFLAGS | grep '\-L/lib '`"; then
            dnl                  if test -n "`echo $F77LIBFLAGS | grep xlf`"; then
            dnl                     F77LIBFLAGS=`echo $F77LIBFLAGS | sed 's/-L\/lib //g'`
            dnl                  fi
            dnl               fi
            dnl            )

            if test -n "`echo $MPILIBS | grep pmpich`" &&
               test -z "`echo $MPILIBS | grep pthread`"; then
                  LIBS_SAVE=$LIBS
                  LIBS="$MPIINCLUDE $MPIFLAGS $MPILIBDIRS $MPILIBS -lpthread"
                  AC_LANG_PUSH(C)
                  AC_TRY_LINK(
dnl                     ifelse(AC_LANG, [C++],

dnl [#ifdef __cplusplus
dnl extern "C"
dnl #endif
dnl ])dnl
[#include "mpi.h"
], [int rank, size;
   int argc;
   char **argv;
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Finalize();
],
                     MPILIBS="$MPILIBS -lpthread")
                  AC_LANG_POP(C)
                  LIBS=$LIBS_SAVE
            fi

            AC_MSG_CHECKING(for MPI include directories)
            AC_MSG_RESULT($MPIINCLUDE)
            AC_MSG_CHECKING(for MPI library directories)
            AC_MSG_RESULT($MPILIBDIRS)
            AC_MSG_CHECKING(for MPI libraries)
            AC_MSG_RESULT($MPILIBS)
            AC_MSG_CHECKING(for other MPI-related flags)
            AC_MSG_RESULT($MPIFLAGS)
         fi
      fi

      AC_CACHE_VAL(casc_cv_mpi_include, casc_cv_mpi_include=$MPIINCLUDE)
      AC_CACHE_VAL(casc_cv_mpi_lib_dirs, casc_cv_mpi_lib_dirs=$MPILIBDIRS)
      AC_CACHE_VAL(casc_cv_mpi_libs, casc_cv_mpi_libs=$MPILIBS)
      AC_CACHE_VAL(casc_cv_mpi_flags, casc_cv_mpi_flags=$MPIFLAGS)
   fi

   AC_SUBST(MPIINCLUDE)
   AC_SUBST(MPILIBDIRS)
   AC_SUBST(MPILIBS)
   AC_SUBST(MPIFLAGS)

])dnl

dnl ********************************************************************
dnl * CASC_FIND_MPI_ALPHA is a special case of CASC_FIND_MPI for the
dnl * compass cluster.  The original CASC_FIND_MPI looks for existence
dnl * of mpCC and mpiCC.  If the former is found it uses native (proprietary)
dnl * mpi and if the latter is found, it uses mpich.  The DECs are a
dnl * special case because mpCC does not exist and mpiCC does, but we want
dnl * to use the native version by default.  Therefore, the original macro
dnl * did not work for this case so I added this one to deal with it.
dnl * AMW 9/00
dnl ********************************************************************

AC_DEFUN([CASC_FIND_MPI_ALPHA],
[

   casc_find_mpi_cache_used=yes

   AC_CACHE_VAL(casc_cv_mpi_include, casc_find_mpi_cache_used=no)
   AC_CACHE_VAL(casc_cv_mpi_libs, casc_find_mpi_cache_used=no)
   AC_CACHE_VAL(casc_cv_mpi_lib_dirs, casc_find_mpi_cache_used=no)
   AC_CACHE_VAL(casc_cv_mpi_flags, casc_find_mpi_cache_used=no)

   if test "$casc_find_mpi_cache_used" = "yes"; then
      AC_MSG_CHECKING(for location of mpi.h)
      MPIINCLUDE=$casc_cv_mpi_include
      AC_MSG_RESULT("\(cached\) $MPIINCLUDE")

      AC_MSG_CHECKING(for MPI library directories)
      MPILIBDIRS=$casc_cv_mpi_lib_dirs
      AC_MSG_RESULT("\(cached\) $MPILIBDIRS")

      AC_MSG_CHECKING(for MPI libraries)
      MPILIBS=$casc_cv_mpi_libs
      AC_MSG_RESULT("\(cached\) $MPILIBS")

      AC_MSG_CHECKING(for other MPI-related flags)
      MPIFLAGS=$casc_cv_mpi_flags
      AC_MSG_RESULT("\(cached\) $MPIFLAGS")
   else


      dnl * Set up user options.  If user uses any of the fist three options,
      dnl * then automatic tests are not run.

      casc_user_chose_mpi=no
      AC_ARG_WITH(mpi-include, [  --with-mpi-include=DIR  mpi.h is in DIR],
                  for mpi_dir in $withval; do
                     MPIINCLUDE="$MPIINCLUDE -I$withval"
                  done; casc_user_chose_mpi=yes)

      AC_ARG_WITH(mpi-libs,
[  --with-mpi-libs=LIBS    LIBS is space-separated list of library names
                          needed for MPI, e.g. \"nsl socket mpi\"],
                  for mpi_lib in $withval; do
                     MPILIBS="$MPILIBS -l$mpi_lib"
                  done; casc_user_chose_mpi=yes)


      AC_ARG_WITH(mpi-lib-dirs,
[  --with-mpi-lib-dirs=DIRS
                          DIRS is space-separated list of directories
                          containing the libraries specified by
                          \`--with-mpi-libs', e.g \"/usr/lib /usr/local/mpi/lib\"],
                  for mpi_lib_dir in $withval; do
                     MPILIBDIRS="-L$mpi_lib_dir $MPILIBDIRS"
                  done; casc_user_chose_mpi=yes)

      dnl * --with-mpi-flags only adds to automatic selections,
      dnl * does not override

      AC_ARG_WITH(mpi-flags,
[  --with-mpi-flags=FLAGS  FLAGS is space-separated list of whatever flags other
                          than -l and -L are needed to link with mpi libraries],
                          MPIFLAGS=$withval)


      if test "$casc_user_chose_mpi" = "no"; then

         dnl * Set defaults for Compass cluster here.  This is the point where
         dnl * we call CASC_SET_MPI in CASC_FIND_MPI macro.

         casc_mpi_include_dir=
         casc_mpi_lib_dirs=
         casc_mpi_libs="mpi rt rpc gs pthread"

         for casc_incl_dir in $casc_mpi_include_dir; do
            MPIINCLUDE="-I$casc_incl_dir $MPIINCLUDE"
         done
         for casc_lib_dir in $casc_mpi_lib_dirs; do
            MPILIBDIRS="-L$casc_lib_dir $MPILIBDIRS"
         done
         for casc_lib in $casc_mpi_libs; do
            MPILIBS="$MPILIBS -l$casc_lib"
         done
      fi


      AC_MSG_CHECKING(for MPI include directories)
      AC_MSG_RESULT($MPIINCLUDE)
      AC_MSG_CHECKING(for MPI library directories)
      AC_MSG_RESULT($MPILIBDIRS)
      AC_MSG_CHECKING(for MPI libraries)
      AC_MSG_RESULT($MPILIBS)
      AC_MSG_CHECKING(for other MPI-related flags)
      AC_MSG_RESULT($MPIFLAGS)

   fi

])dnl



dnl Define a macro for supporting SILO

AC_DEFUN([CASC_SUPPORT_SILO],[

# Begin CASC_SUPPORT_SILO
# Defines silo_PREFIX silo_INCLUDES and silo_LIBS if with-silo is specified.
AC_ARG_WITH(silo,
[ --with-silo[=PATH]  Use SILO and optionally specify where SILO is installed.],
, with_silo=no)

case "$with_silo" in
  no)
    AC_MSG_NOTICE([configuring without SILO support])
    : Do nothing
  ;;
  yes)
    # SILO install path was not specified.
    # Look in a couple of standard locations to probe if 
    # SILO header files are there.
    AC_MSG_CHECKING([for SILO installation])
    for dir in /usr /usr/local; do
      if test -f ${dir}/include/silo.h; then
        silo_PREFIX=${dir}
        break
      fi
    done
    AC_MSG_RESULT([$silo_PREFIX])
  ;;
  *)
    # SILO install path was specified.
    AC_MSG_CHECKING([for SILO installation])
    silo_PREFIX=$with_silo
    silo_INCLUDES="-I${silo_PREFIX}/include"
    if test -f ${silo_PREFIX}/include/silo.h; then
        AC_MSG_RESULT([$silo_PREFIX])
    else
        AC_MSG_RESULT([$silo_PREFIX])
        AC_MSG_ERROR([SILO not found in $with_silo])
    fi
  ;;
esac

# Determine which SILO library is built
if test "${silo_PREFIX+set}" = set; then
   AC_MSG_CHECKING([for SILO library])
   if test -f ${silo_PREFIX}/lib/libsilo.a; then
      silo_LIBS='-lsilo'
      AC_MSG_RESULT([using $silo_LIBS])
   elif test -f ${silo_PREFIX}/lib/libsiloh5.a; then
      silo_LIBS='-lsiloh5'
      AC_MSG_RESULT([using $silo_LIBS])
   else
      AC_MSG_RESULT([using $silo_LIBS])
      AC_MSG_ERROR([Could not fine silo library in $silo_PREFIX])
   fi

   silo_LIBS="-L${silo_PREFIX}/lib ${silo_LIBS}"
fi

# END CASC_SUPPORT_SILO

])dnl End definition of CASC_SUPPORT_SILO

dnl Define a macro for supporting VALGRIND

AC_DEFUN([CASC_SUPPORT_VALGRIND],[

# Begin CASC_SUPPORT_VALGRIND
# Defines valgrind_EXE
AC_ARG_WITH(valgrind,
[ --with-valgrind[=PATH]  Use VALGRIND and optionally specify where VALGRIND is installed.],
, with_valgrind=no)

case "$with_valgrind" in
  no)
    AC_MSG_NOTICE([configuring without VALGRIND support])
    : Do nothing
  ;;
  yes)
    # VALGRIND install path was not specified.
    # Look in a couple of standard locations to probe if 
    # VALGRIND header files are there.
    AC_MSG_CHECKING([for VALGRIND installation])
    for dir in /usr /usr/local; do
      if test -f ${dir}/bin/valgrind; then
        valgrind_PREFIX=${dir}
        break
      fi
    done
    AC_MSG_RESULT([$valgrind_PREFIX])
  ;;
  *)
    # VALGRIND install path was specified.
    AC_MSG_CHECKING([for VALGRIND installation])
    valgrind_PREFIX=$with_valgrind
    ;;
esac

if test "${valgrind_PREFIX+set}" = set 
then
   valgrind_EXE="${valgrind_PREFIX}/bin/valgrind"
   if test -f ${valgrind_PREFIX}/bin/valgrind; then
      AC_MSG_RESULT([$valgrind_PREFIX])
   else
      AC_MSG_RESULT([$valgrind_PREFIX])
      AC_MSG_ERROR([VALGRIND not found in $with_valgrind])
   fi
fi

# END CASC_SUPPORT_VALGRIND

])dnl End definition of CASC_SUPPORT_VALGRIND

dnl Define macros for supporting XLC 

AC_DEFUN([CASC_CXX_STD_FILL_N_RETURNS_VOID],[

# Begin CASC_CXX_STD_FILL_N_RETURNS_VOID
# Defines CASC_STD_FILL_N_RETURNS_VOID

# Check if std::fill_n returns a void, older XLC compilers do this.
   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether std:fill_n returns void)
   AC_LANG_PUSH(C++)
   CASC_PUSH_COMPILER_STATE
   AC_COMPILE_IFELSE([

      #include <vector>

      template void std::fill_n<unsigned int*, int, int>(unsigned int*, int, int const&);      

      ], 
      casc_std_fill_n_returns_void=yes,
      casc_std_fill_n_returns_void=no)
   CASC_POP_COMPILER_STATE
   AC_LANG_POP
   AC_MSG_RESULT($casc_std_fill_n_returns_void)

   if test "$casc_std_fill_n_returns_void" = yes; then
      AC_DEFINE([CASC_STD_FILL_N_RETURNS_VOID], 1, [Define if std::fill_n returns void])
   fi


# END CASC_CXX_STD_FILL_N_RETURNS_VOID

])dnl End definition of CASC_CXX_STD_FILL_N_RETURNS_VOID






# ===========================================================================
#
# SYNOPSIS
#
#   CHECK_SZLIB()
#
# DESCRIPTION
#
#   This macro searches for an installed szlib library. If nothing was
#   specified when calling configure, it searches first in /usr/local and
#   then in /usr. If the --with-szlib=DIR is specified, it will try to find
#   it in DIR/include/szlib.h and DIR/lib/libsz.a. If --without-szlib is
#   specified, the library is not searched at all.
#
#   If either the header file (szlib.h) or the library (libsz) is not found,
#   the configuration exits on error, asking for a valid szlib installation
#   directory or --without-szlib.
#
#   The macro defines the symbol HAVE_LIBSZ if the library is found. You
#   should use autoheader to include a definition for this symbol in a
#   config.h file. Sample usage in a C/C++ source is as follows:
#
#     #ifdef HAVE_LIBSZ
#     #include <szlib.h>
#     #endif /* HAVE_LIBSZ */
#

AC_DEFUN([CHECK_SZLIB],
#
# DEFINES :
#	        szlib_PREFIX
#		szlib_INCLUDES
#		szlib_LIBS
#
[AC_MSG_CHECKING(if szlib is wanted)
AC_ARG_WITH(szlib,
[  --with-szlib=DIR root directory path of szlib installation [DIR defaults to
                    /usr/local or /usr if not found in /usr/local]
  --without-szlib to disable szlib usage completely [the default]],
[if test "$withval" != no ; then
  AC_MSG_RESULT(yes)
  if test "$withval" == yes ;
  then
     SZLIB_HOME=/usr/local
  else
     SZLIB_HOME="$withval"
  fi
  if test ! -d "$SZLIB_HOME"
  then
    AC_MSG_WARN([Sorry, $SZLIB_HOME does not exist, checking usual places])
    SZLIB_HOME=/usr/local
    if test ! -f "${SZLIB_HOME}/include/szlib.h"
    then
       SZLIB_HOME=/usr
    fi
  fi
else
  AC_MSG_RESULT(no)
fi])

#
# Locate szlib, if wanted
#
if test -n "${SZLIB_HOME}"
then
        SZLIB_OLD_LDFLAGS=$LDFLAGS
        SZLIB_OLD_CPPFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS -L${SZLIB_HOME}/lib"
        CPPFLAGS="$CPPFLAGS -I${SZLIB_HOME}/include"
        AC_LANG_SAVE
        AC_LANG_C
        AC_CHECK_LIB(z, inflateEnd, [szlib_cv_libsz=yes], [szlib_cv_libsz=no])
        AC_CHECK_HEADER(szlib.h, [szlib_cv_szlib_h=yes], [szlib_cv_szlib_h=no])
        AC_LANG_RESTORE
        if test "$szlib_cv_libsz" = "yes" -a "$szlib_cv_szlib_h" = "yes"
        then
	        szlib_PREFIX="${SZLIB_HOME}"
		szlib_INCLUDES="-I${SZLIB_HOME}/include"
		szlib_LIBS="-L${SZLIB_HOME}/lib -lsz"
                #
                # If both library and header were found, use them
                #
                AC_CHECK_LIB(z, inflateEnd)
                AC_MSG_CHECKING(szlib in ${SZLIB_HOME})
                AC_MSG_RESULT(ok)
        else
                #
                # If either header or library was not found, revert and bomb
                #
                AC_MSG_CHECKING(szlib in ${SZLIB_HOME})
                LDFLAGS="$SZLIB_OLD_LDFLAGS"
                CPPFLAGS="$SZLIB_OLD_CPPFLAGS"
                AC_MSG_RESULT(failed)
                AC_MSG_ERROR(either specify a valid szlib installation with --with-szlib=DIR or disable szlib usage with --without-szlib)
        fi
fi

])

# ===========================================================================
#           http://www.nongnu.org/autoconf-archive/check_zlib.html
# ===========================================================================
#
# SYNOPSIS
#
#   CHECK_ZLIB()
#
# DESCRIPTION
#
#   This macro searches for an installed zlib library. If nothing was
#   specified when calling configure, it searches first in /usr/local and
#   then in /usr. If the --with-zlib=DIR is specified, it will try to find
#   it in DIR/include/zlib.h and DIR/lib/libz.a. If --without-zlib is
#   specified, the library is not searched at all.
#
#   If either the header file (zlib.h) or the library (libz) is not found,
#   the configuration exits on error, asking for a valid zlib installation
#   directory or --without-zlib.
#
#   The macro defines the symbol HAVE_LIBZ if the library is found. You
#   should use autoheader to include a definition for this symbol in a
#   config.h file. Sample usage in a C/C++ source is as follows:
#
#     #ifdef HAVE_LIBZ
#     #include <zlib.h>
#     #endif /* HAVE_LIBZ */
#
# LICENSE
#
#   Copyright (c) 2008 Loic Dachary <loic@senga.org>
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([CHECK_ZLIB],
#
# DEFINES :
#	        zlib_PREFIX
#		zlib_INCLUDES
#		zlib_LIBS
#
[AC_MSG_CHECKING(if zlib is wanted)
AC_ARG_WITH(zlib,
[  --with-zlib=DIR root directory path of zlib installation [DIR defaults to
                    /usr/local or /usr if not found in /usr/local]
  --without-zlib to disable zlib usage completely [the default]],
[if test "$withval" != no ; then
  AC_MSG_RESULT(yes)
  if test "$withval" == yes ;
  then
     ZLIB_HOME=/usr/local
  else
     ZLIB_HOME="$withval"
  fi
  if test ! -d "$ZLIB_HOME"
  then
    AC_MSG_WARN([Sorry, $ZLIB_HOME does not exist, checking usual places])
    ZLIB_HOME=/usr/local
    if test ! -f "${ZLIB_HOME}/include/zlib.h"
    then
       ZLIB_HOME=/usr
    fi
  fi
else
  AC_MSG_RESULT(no)
fi],
  [AC_MSG_RESULT(no)]
)


#
# Locate zlib, if wanted
#
if test -n "${ZLIB_HOME}"
then
        ZLIB_OLD_LDFLAGS=$LDFLAGS
        ZLIB_OLD_CPPFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS -L${ZLIB_HOME}/lib"
        CPPFLAGS="$CPPFLAGS -I${ZLIB_HOME}/include"
        AC_LANG_SAVE
        AC_LANG_C
        AC_CHECK_LIB(z, inflateEnd, [zlib_cv_libz=yes], [zlib_cv_libz=no])
        AC_CHECK_HEADER(zlib.h, [zlib_cv_zlib_h=yes], [zlib_cv_zlib_h=no])
        AC_LANG_RESTORE
        if test "$zlib_cv_libz" = "yes" -a "$zlib_cv_zlib_h" = "yes"
        then
	        zlib_PREFIX="${ZLIB_HOME}"
		zlib_INCLUDES="-I${ZLIB_HOME}/include"
		zlib_LIBS="-L${ZLIB_HOME}/lib -lz"
                #
                # If both library and header were found, use them
                #
                AC_CHECK_LIB(z, inflateEnd)
                AC_MSG_CHECKING(zlib in ${ZLIB_HOME})
                AC_MSG_RESULT(ok)
        else
                #
                # If either header or library was not found, revert and bomb
                #
                AC_MSG_CHECKING(zlib in ${ZLIB_HOME})
                LDFLAGS="$ZLIB_OLD_LDFLAGS"
                CPPFLAGS="$ZLIB_OLD_CPPFLAGS"
                AC_MSG_RESULT(failed)
                AC_MSG_ERROR(either specify a valid zlib installation with --with-zlib=DIR or disable zlib usage with --without-zlib)
        fi
fi

])

dnl $Id$

dnl Determines which compiler is being used.
dnl This check uses the compiler behavior when possible.
dnl For some compiler, we resort to a best guess,
dnl because we do not know a foolproof way to get the info.

dnl Much of the information used here came from the very
dnl helpful predef project (http://predef.sourceforge.net/).




dnl Simple wrappers to allow using CASC_INFO_CXX_ID_NAMES and
dnl CASC_INFO_CC_ID_NAMES without arguments.
dnl The names CC_ID and CC_VERSION are used for the C compiler id and version.
dnl The names CXX_ID and CXX_VERSION are used for the C++ compiler id and version.
AC_DEFUN([CASC_INFO_CXX_ID],[
  CASC_INFO_CXX_ID_NAMES(CXX_ID,CXX_VERSION)
])
AC_DEFUN([CASC_INFO_CC_ID],[
  CASC_INFO_CC_ID_NAMES(CC_ID,CC_VERSION)
])
AC_DEFUN([CASC_INFO_CC_CXX_ID],[
  AC_REQUIRE([CASC_INFO_CC_ID])
  AC_REQUIRE([CASC_INFO_CXX_ID])
])


dnl CASC_INFO_CXX_ID and CASC_INFO_C_ID determine which C or C++ compiler
dnl is being used.
# Set the variables CXX_ID or C_ID as follows:
# Gnu		-> gnu
# SUNWspro	-> sunpro
# Dec		-> dec
# KCC		-> kai
# Intel		-> intel
# SGI		-> sgi
# IBM xlc	-> xlc


AC_DEFUN([CASC_INFO_CXX_ID_NAMES],
dnl Arguments are:
dnl 1. Name of variable to set to the ID string.
dnl 2. Name of variable to set to the version number.
[
# Start macro CASC_INFO_CXX_ID_NAMES
  AC_REQUIRE([AC_PROG_CXXCPP])
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  CASC_AC_LOG(CXXP is $CXX)
  CASC_AC_LOG(CXXCPP is $CXXCPP)

  $1=unknown
  $2=unknown

dnl Do not change the following chain of if blocks into a case statement.
dnl We may eventually have a compiler that must be tested in a different
dnl method


  # Check if it is a Sun compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CXX is sunpro)
    AC_EGREP_CPP([^0x[0-9]+],__SUNPRO_CC,
      $1=sunpro
      # SUN compiler defines __SUNPRO_CC to the version number.
      echo __SUNPRO_CC > conftest.C
      $2=`${CXXCPP} conftest.C | sed -n 2p`
      rm -f conftest.C
    )
  fi


  # Check if it is a Intel compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CXX is intel)
    AC_EGREP_CPP(^yes,
#ifdef __INTEL_COMPILER
yes;
#endif
,
      $1=intel
      # Intel compiler defines __INTEL_COMPILER to the version number.
      echo __INTEL_COMPILER > conftest.C
      $2=`${CXXCPP} conftest.C | sed -n 2p`
      rm -f conftest.C
    )
  fi


  # Check if it is a GNU compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CXX is gnu)
    AC_EGREP_CPP(^yes,
#ifdef __GNUC__
yes;
#endif
,
    $1=gnu
    # GNU compilers output version number with option --version.
    # Alternatively, it also defines the macros __GNUC__,
    # GNUC_MINOR__ and __GNUC_PATCHLEVEL__
    [[$2=`$CXX --version | sed -e 's/[^0-9]\{0,\}\([^ ]\{1,\}\).\{0,\}/\1/' -e 1q`]]
    )
  fi


  # Check if it is a DEC compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CXX is dec)
    AC_EGREP_CPP(^1,__DECCXX,
      $1=dec
      # DEC compiler defines __DECCXX_VER to the version number.
      echo __DECCXX_VER > conftest.C
      $2=`${CXXCPP} conftest.C | sed -n 2p`
      rm -f conftest.C
    )
  fi


  # Check if it is a KAI compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CXX is kai)
    AC_EGREP_CPP(^1,__KCC,
      $1=kai
      # KCC compiler defines __KCC_VERSION to the version number.
      echo __KCC_VERSION > conftest.C
      $2=`${CXXCPP} conftest.C | sed -n 2p`
      rm -f conftest.C
    )
  fi


  # Check if it is a SGI compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CXX is sgi)
    AC_EGREP_CPP(^1,__sgi,
      $1=sgi
      # SGI compiler defines _COMPILER_VERSION to the version number.
      echo _COMPILER_VERSION > conftest.C
      $2=`${CXXCPP} conftest.C | sed /^\\#/d`
      rm -f conftest.C
    )
  fi


  # Check if it is a IBM compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CXX is xlc)
    AC_EGREP_CPP(^yes,
#ifdef __xlC__
yes;
#endif
,
    $1=xlc
    # IBM compiler defines __xlC__ to the version number.
    echo __xlC__ > conftest.C
    $2=`${CXXCPP} conftest.C | sed /^\\#/d`
    rm -f conftest.C
    )
  fi


  AC_LANG_RESTORE
  CASC_AC_LOG_VAR(CXX_ID CXX_VERSION)
# End macro CASC_INFO_CXX_ID_NAMES
])





AC_DEFUN([CASC_INFO_CC_ID_NAMES],
dnl Arguments are:
dnl 1. Name of variable to set to the ID string.
dnl 2. Name of variable to set to the version number.
[
# Start macro CASC_INFO_CC_ID_NAMES
  AC_REQUIRE([AC_PROG_CPP])
  AC_LANG_SAVE
  AC_LANG_C
  CASC_AC_LOG(CC is $CC)
  CASC_AC_LOG(CPP is $CPP)

  $1=unknown
  $2=unknown

dnl Do not change the following chain of if blocks into a case statement.
dnl We may eventually have a compiler that must be tested in a different
dnl method


  # Check if it is a Sun compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CC is sunpro)
    AC_EGREP_CPP([^ 0x[0-9]+],__SUNPRO_C,
      $1=sunpro
      # SUN compiler defines __SUNPRO_C to the version number.
      echo __SUNPRO_C > conftest.c
      $2=`${CPP} ${CPPFLAGS} conftest.c | sed -n -e 's/^ //' -e 2p`
      rm -f conftest.c
    )
  fi


  # Check if it is a Intel compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CC is intel)
    AC_EGREP_CPP(^yes,
#ifdef __INTEL_COMPILER
yes;
#endif
,
      $1=intel
      # Intel compiler defines __INTEL_COMPILER to the version number.
      echo __INTEL_COMPILER > conftest.C
      $2=`${CPP} conftest.C | sed -n 2p`
      rm -f conftest.C
    )
  fi


  # Check if it is a GNU compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CC is gnu)
    AC_EGREP_CPP(^yes,
#ifdef __GNUC__
yes;
#endif
,
    $1=gnu
    [[$2=`$CC --version | sed -e 's/[^0-9]\{0,\}\([^ ]\{1,\}\).\{0,\}/\1/' -e 1q`]]
    )
  fi


  # Check if it is a DEC compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CC is dec)
    AC_EGREP_CPP(^ 1,__DECC,
      $1=dec
      # DEC compiler defines __DECC_VER to the version number.
      echo __DECC_VER > conftest.c
      $2=`${CPP} ${CPPFLAGS} conftest.c | sed -n -e 's/^ //' -e 2p`
      rm -f conftest.c
    )
  fi


  # Check if it is a KAI compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CC is kai)
    AC_EGREP_CPP(^1,__KCC,
      $1=kai
      # KCC compiler defines __KCC_VERSION to the version number.
      echo __KCC_VERSION > conftest.c
      $2=`${CPP} ${CPPFLAGS} conftest.c | sed -n 2p`
      rm -f conftest.c
    )
  fi


  # Check if it is a SGI compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CC is sgi)
    AC_EGREP_CPP(^1,__sgi,
      $1=sgi
      # SGI compiler defines _COMPILER_VERSION to the version number.
      echo _COMPILER_VERSION > conftest.c
      $2=`${CPP} ${CPPFLAGS} conftest.c | sed /^\\#/d`
      rm -f conftest.c
    )
  fi


  # Check if it is a IBM compiler.
  if test $$1 = unknown; then
    CASC_AC_LOG(checking if $CC is xlc)
    if echo "$host_os" | grep "aix" >/dev/null ; then
      # The wretched IBM shell does not eval correctly,
      # so we have to help it with a pre-eval eval statement.
      ac_cpp=`eval "echo $ac_cpp"`
      save_ac_cpp=$ac_cpp
      CASC_AC_LOG(ac_cpp is temporarily set to $ac_cpp)
    else
      save_ac_cpp=
    fi
    CASC_AC_LOG(ac_cpp is $ac_cpp)
    AC_EGREP_CPP(^yes,
#ifdef __xlC__
yes;
#endif
,
    $1=xlc
    # IBM compiler defines __xlC__ to the version number.
    echo __xlC__ > conftest.C
    $2=`${CPP} conftest.C | sed /^\\#/d`
    rm -f conftest.C
    )
    test "$save_ac_cpp" && ac_cpp=$save_ac_cpp
    CASC_AC_LOG(ac_cpp is restored to $ac_cpp)
  fi


  AC_LANG_RESTORE
  CASC_AC_LOG_VAR(CC_ID CC_VERSION)
# End macro CASC_INFO_CC_ID_NAMES
])

dnl $Id$


AC_DEFUN([CASC_TYPE_BOOL],[

# Start macro CASC_TYPE_BOOL

AC_MSG_CHECKING(checking whether bool type is broken)

AC_CACHE_VAL(btng_cv_type_bool_broken, [

  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS

  AC_TRY_COMPILE(, bool b = true; ,
    # bool is not broken.
    btng_cv_type_bool_broken=no
    ,
    # bool is broken.
    btng_cv_type_bool_broken=yes
  )	dnl End AC_TRY_COMPILE call

  AC_LANG_RESTORE

])	dnl End AC_CACHE_VAL call

AC_MSG_RESULT($btng_cv_type_bool_broken)

if test "$btng_cv_type_bool_broken" = yes; then
  AC_DEFINE([BOOL_IS_BROKEN],1,Define if bool type is not properly supported)
fi


# End macro CASC_TYPE_BOOL

])	dnl End of COMPILE_BOOLEAN_MACRO definition.

AC_DEFUN([CASC_AC_LOG],[echo "configure:__oline__:" $1 >&AC_FD_CC])

AC_DEFUN([CASC_AC_LOG_VAR],[
dnl arg1 is list of variables to log.
dnl arg2 (optional) is a label.
dnl
dnl This macro makes code that write out at configure time
dnl label: x is '...'
dnl if x is set and
dnl label: x is unset
dnl otherwise.
define([btng_log_label],ifelse($2,,,[$2: ]))
btng_log_vars="$1"
for btng_log_vars_index in $btng_log_vars ; do
  eval "test \"\${${btng_log_vars_index}+set}\" = set"
  if test $? = 0; then
    btng_log_vars_value="'`eval echo \\${$btng_log_vars_index}`'";
  else
    btng_log_vars_value="unset";
  fi
  CASC_AC_LOG("btng_log_label$btng_log_vars_index is $btng_log_vars_value");
dnl
dnl This is a shorter version, but it does not work for some Bourne shells
dnl due to misinterpretation of the multiple backslashes
dnl CASC_AC_LOG("btng_log_label$btng_log_vars_index is `eval if test \\\"\$\{$btng_log_vars_index+set\}\\\"\; then echo \\\""'"\$\{$btng_log_vars_index\}"'"\\\"\; else echo 'unset'\; fi`")
done
undefine([btng_log_label])
])

dnl Define a macro for supporting generalized serial-parallel run.


AC_DEFUN([SAMRAI_SERPA],[
dnl Support a generalized way to run a program in serial or parallel mode.
dnl (serpa is serial + parallel).
dnl
dnl Arg1: name (may include relative path) to give to shell script created.
dnl   If omitted, the name will be serpa-run, in the top level directory.
dnl
dnl This macro brings out an automake bug.  See "IMPORTANT NOTE:" below
dnl for comments and how to circumvent.

# Begin macro $0


# The generalized script needs to know about the system
# so it can form the appropriate parallel run command line.
AC_REQUIRE([AC_CANONICAL_TARGET])


dnl The parallel run program (usually mpirun) is hardwired for now,
dnl but it should be more flexible.
dnl I plan to eventualy use the following variable,
dnl but it is not currently used.
AC_ARG_WITH(parallel-run-bin,[
  --with-parallel-run-bin=STRING
		Specify the parallel run binary (i.e., mpirun, dmpirun,
		poe, etc) to be used in the serial-parallel run script.
],[
case "$with_parallel_run_bin" in
  no)
    # Set PARALLEL_RUN_BIN to blank to remove it from the serpa-run script.
    PARALLEL_RUN_BIN=;;
  yes)
    # Unset the variable PARALLEL_RUN_BIN so it is guessed in the next step.
    unset PARALLEL_RUN_BIN;
    ;;
  *)
    # Set PARALLEL_RUN_BIN to user specification.
    PARALLEL_RUN_BIN="$with_parallel_run_bin"
    ;;
esac
CASC_AC_LOG_VAR(with_parallel_run_bin PARALLEL_RUN_BIN target_os, with-parallel-run-bin given)

],[
unset PARALLEL_RUN_BIN;
CASC_AC_LOG_VAR(with_parallel_run_bin PARALLEL_RUN_BIN target_os, with-parallel-run-bin NOT given)
])
# If PARALLEL_RUN_BIN is unset, guess it.
if test ! "${PARALLEL_RUN_BIN+set}" = set; then
  case "$target_os" in
    osf*) PARALLEL_RUN_BIN=dmpirun
    ;;
    *) PARALLEL_RUN_BIN=mpirun
    ;;
  esac
  CASC_AC_LOG_VAR(with_parallel_run_bin PARALLEL_RUN_BIN target_os, after setting PARALLEL_RUN_BIN)
fi
AC_SUBST(PARALLEL_RUN_BIN)


define(btng_serpa_run_fn,ifelse($1,,serpa-run,$1))
dnl Create serpa-run.in file (at autoconf time, before configure time)
dnl File serpa-run.in will be used at configure time to create serpa-run
dnl in the compile tree.
syscmd([ file_name=']btng_serpa_run_fn[' cat <<'__EOM__'> ${file_name}.in
#!/bin/sh

# This script runs an executable program serially or in parallel.

# The syntax for running programs in parallel differs from platform to platform,
# software to software.  This script automatically uses the correct syntax for
# the environments that it knows about.  For the environments it does not know,
# it has a default that may work.

# This script defines the set of tests to put the program through.
# This set is parametrized by these variables:
# nproc_list: space- or comma-separated list of number of processors to use.


# The name and directory of this script
# (Do not rely on existence of dirname and basename programs.)
script_name=`echo ${0} | sed -e 's:.*/::'`;
dir_name=`echo ${0} | sed -e 's:^\([^/]*\)$:./\1:' -e 's:/[^/]*$::'`;

# Define a way to gracefully die from within this script.
# The die function is similar to Perl's.
# Arguments are: <exit_value> <exit_message>
die () {
    echo "ERROR_MESSAGE_FROM ${script_name}:"
    if [ -n "${2}" ]; then echo "Error ${2}"; fi
    if [ -n "${1}" ]; then exit ${1}; fi
    exit 99
}


# When no argument is given, print the help message and exit.
if [ "${1}" = '-h' ] || [ "${1}" = '--help' ] || [ ${#} -lt 2 ] ; then
  cat <<-_EOF_
	Usage: ${script_name} <list of nproc> <program name>

	This is a generalized script to run a program in serial
	and/or parallel.

	The number of processors are given by <list of nproc>,
	which must be a comma- or space-delimited list of integers.
	Example: "${script_name} 1,2,5 parallel_program" runs
	parallel_program 3 times with 1, 2 and 5 processors in turn.

	The special case of number of processors = 0 means to run serially.
	This differs from one processor in that one processor means to run
	in parallel, but with one processor.

	This script exits with an error if any instance of
	running the program fails.

	If SERPA_REDIRECT_OUTPUT_TO is defined in the environment,
	standard output of the program is redirected there.
	If SERPA_REDIRECT_ERRORS_TO is defined in the environment,
	error output of the program is redirected there.
	Outputs directly from THIS script are NOT affected
	by these environment variables.

	If SERPA_MAX_FAILS is defined, ${script_name} will
	continue until that many failures before exiting.
	The default is 1 failure (${script_name} exits on the
	first failure).  The exit value will always be the
	failure count.

	When a serial run os made, if SERPA_PERFORMANCE_FILE
	is defined and is the name of a file, that file is searched
	for performance data.  The run is timed and the result
	is compared to the performance data found.
	_EOF_
exit
fi

# How to invoke a GNU compatible time command
if [ -z "${SERPA_GNUTIME}" ]; then
    SERPA_GNUTIME=/usr/bin/time
fi

# Percent difference for performance comparisons
if [ -z "${SERPA_DIFFERENCE}" ]; then
    SERPA_DIFFERENCE=0.1
fi


# Check the redirection environment variables.
# See the help message for how these are used.
if [ -n "${SERPA_REDIRECT_OUTPUT_TO}" ]; then
  output_redirection_string="1> ${SERPA_REDIRECT_OUTPUT_TO}"
else
  unset output_redirection_string
fi
if [ -n "${SERPA_REDIRECT_ERRORS_TO}" ]; then
  errors_redirection_string="2> ${SERPA_REDIRECT_ERRORS_TO}"
else
  unset errors_redirection_string
fi

# Save the arguments for later reference.
serpa_args="${@}"


# Get the number of processors from the first argument of this script.
nproc_list=${1}; shift;
# We allow comma-delimited lists, so we now remove those commas.
echo ${nproc_list} | grep '^[0-9 ,]\{1,\}$' > /dev/null	\
  || die 1 "Invalid number of processors list '${nproc_list}'"
nproc_list=`echo ${nproc_list} | sed 's/,/ /g'`

# Get the program name from the next argument of this script.
program=${1}; shift;
test -x ${program} || die 1 "No execute permission on '${program}'"


# Determine host name for use below.
if [ ${?}{HOST} ]; then
  HOST=`uname -n`
  export HOST
fi


# Use additional_env to specify additional environments that
# should be set before running.  Although you can, do not use
# this to set environments for programs.  It is meant to be
# environments for the parallel execution program such as
# mpirun.
additional_env=


# Variables required for using parallelrun.  These may not be needed,
# but if they are needed and unset, there will be an error.
# Some of these strings are determined at configure time.
parallelrun_prog="@PARALLEL_RUN_BIN@";	# Name of parallelrun program.


# Determine what platform we are on.
target_cpu=@target_cpu@
target_os=@target_os@
target_vendor=@target_vendor@


# Determine the machine file name.
serpa_machine_file='serpa.machines'


# Define a function to run a program in parallel.
# We have to do this because different environments require different syntaxes.
# The function arguments are:
#   program name
#   number of processor
#
# We define the function in a big if-else statement based on the
# target computer.  We assume that there is a one-to-one correspondence.
# If there is not, there may have to be a nested if-structure.
#
if echo "${HOST}" | egrep '^(blue|frost)[0-9]' > /dev/null ; then

  # For blue, a specific singleton platform.
  run_multiproc () {
    MP_RESD="YES"
    MP_HOSTFILE=""
    MP_EUILIB=us
    MP_EUIDEVICE=css0
    export MP_RESD MP_HOSTFILE MP_EUILIB MP_EUIDEVICE
  # IBM shell functions eat their parameters after the first assignment
  # from them so save the parameters first.
    case "${HOST}" in
      blue*) proc_per_node=4;;
      frost*) proc_per_node=16;;
    esac
    program=${1}; nproc=${2}; shift 2
    nodes=`expr 1 + \( ${nproc} - 1 \) / ${proc_per_node}`
    com="${additional_env} ${parallelrun_prog} ${program}"
    com="${com} -rmpool 0 -nodes ${nodes} -procs ${nproc} ${@}"
    # com="${program} ${@} -procs ${nproc}"
    echo ${com}
    eval ${com} ${output_redirection_string} ${errors_redirection_string}
    return ${?}
  }

elif echo "${HOST}" | grep '^tc2k' > /dev/null ; then

  # For tc2k, a specific singleton platform (contributed by Brian Miller).
  run_multiproc () {
    program=${1}; nproc=${2}; shift 2
    com="${additional_env} prun -n ${nproc} ${program} ${@}"
    com="${program} ${@} -p ${nproc}"
    echo ${com}
    eval ${com} ${output_redirection_string} ${errors_redirection_string}
    return ${?}
  }

elif echo "${HOST}" | egrep '^(mcr|pengra|thunder)[0-9]' > /dev/null ; then

  # For LC Linux clusters using srun.
  run_multiproc () {
    program=${1}; nproc=${2}; shift 2
    com="${additional_env} srun -n${nproc} -p pdebug ${program} ${@}"
    echo ${com}
    eval ${com} ${output_redirection_string} ${errors_redirection_string}
    return ${?}
  }

elif echo "${target_os}" | grep '^osf'	> /dev/null ; then

  # For Dec OSF.
  run_multiproc () {
    program=${1}; nproc=${2}; shift 2
    com="${additional_env} dmpirun -np ${nproc} ${program} ${@}"
    echo ${com}
    eval ${com} ${output_redirection_string} ${errors_redirection_string}
    return ${?}
  }

elif echo "${target_os}" | grep '^solaris' > /dev/null	\
  || echo "${target_os}" | grep '^irix' > /dev/null	\
  || echo "${target_os}" | grep '^linux' > /dev/null	\
  ; then

  # Most platforms fall into this case.
  run_multiproc () {
    program=${1}; nproc=${2}; shift 2
    com="${additional_env} ${parallelrun_prog}"
    com="${com} -machinefile ${dir_name}/${serpa_machine_file} -np ${nproc} ${program} ${@}"
    echo ${com}
    eval ${com} ${output_redirection_string} ${errors_redirection_string}
    return ${?}
  }

else

  # Simplest case.
  # This is generic.  It may not work, but it is our best guess
  # without knowledge of the system.
  run_multiproc () {
    program=${1}; nproc=${2}; shift 2
    com="${additional_env} ${program} -np ${nproc} ${@}"
    echo ${com}
    eval ${com} ${output_redirection_string} ${errors_redirection_string}
    return ${?}
  }

fi


# Initialize the failure count.
serpa_num_failures=0
test -z "${SERPA_MAX_FAILS}" && SERPA_MAX_FAILS=1

# Run the program.
for nproc in ${nproc_list}; do
  cat <<-_EOM_
	${script_name}================================================
	RUNNING: ${script_name} ${nproc} ${program} ${@}
	${script_name}::::::::::::::::::::::::::::::::::::::::::::::::
	_EOM_
  if test "${nproc}" = 0 ; then
    # Run serially.

    # If performance file exists then with performance testing
    # otherwise do a normal sequential run
    if test "${SERPA_PERFORMANCE_FILE+set}" = set &&
       test -f "${SERPA_PERFORMANCE_FILE}"; then

	# Run with performance check. 
	com="${additional_env} ${program} ${@}"
	echo "${com}"
	${SERPA_GNUTIME} -o $$.time -f "%e %t" ${com} \
	    ${output_redirection_string} ${errors_redirection_string}
	exit_value=${?}

	# Parse temporary file to get recorded time/memory usage
	read etime memory < $$.time
	rm $$.time
	
	# Report the run time/memory usage 
	echo "serpa-perf ${nproc} <${com}> ${etime} ${memory}"

	# Check if time is out of bounds
	search=`grep "serpa-perf ${nproc} <${com}>" ${SERPA_PERFORMANCE_FILE}`
	if test "$?" = "0";then
	    search=`echo $search | sed 's/.*>//'`
	    read compare_time compare_memory <<EOF	    
$search
EOF
        else
	    compare_time=999999
	    compare_memory=0
	fi

	difference=`bc -l <<EOF
a=${etime}
b=${compare_time}
t=${SERPA_DIFFERENCE}
d=a-b
if ( d < 0 ) {
	d = -d
}
p=d/b
print (p < t), "\n"
EOF`
	if test "${difference}" = "1"; then
	    echo "PERFORMANCE PASSED: target ${compare_time} (seconds) current ${etime} sec"
	else
	    echo "PERFORMANCE FAIL: target ${compare_time} (seconds) current ${etime}"
	fi

	# Currently don't do memory check since it is not reporting anything.

    else
	# Run without performance testing
	com="${additional_env} ${program} ${@}"
	echo "${com}"
	eval ${com} ${output_redirection_string} ${errors_redirection_string}
	exit_value=${?}
    fi
  else
    # Because parallel runs sometimes fails due to problems unrelated
    # to the program, we give it several tries before declaring failure.
    if test -z "${SERPA_PARALLEL_TRIES}"; then SERPA_PARALLEL_TRIES=1; fi
    c=1
    while test ${c} -le ${SERPA_PARALLEL_TRIES} ; do
      if test ${c} -gt 1 ; then
        echo "possibly failed.  TRY number ${c}"
      fi
      run_multiproc ${program} ${nproc} "${@}"
      exit_value=${?}
      if [ ${exit_value} = 0 ]; then
        break
      fi
      c=`expr ${c} + 1`
    done
  fi
  # Report pass or fail.
  # pf_string=FAILED; test "${exit_value}" = 0 && pf_string=PASSED
  pf_string=FAILED; test "${exit_value}" = 0 && pf_string=COMPLETED
  cat <<-_EOM_
	${script_name}::::::::::::::::::::::::::::::::::::::::::::::::
	${pf_string}: ${script_name} ${nproc} ${program} ${@}
	${script_name}================================================
	_EOM_
  # Count failures and possibly exit.
  if test ${exit_value} -ne 0 ; then
    serpa_num_failures=`expr ${serpa_num_failures} + 1`;
    if test ${serpa_num_failures} -eq ${SERPA_MAX_FAILS}; then
      die 1 "FAIL running ${program} with ${nproc} processors"
    fi
  fi
done

if test ${serpa_num_failures} -eq 0; then
  echo "${script_name} ${serpa_args} passed."
else
  echo "${script_name} ${serpa_args} had ${serpa_num_failures} failures."
fi
exit ${serpa_num_failures}

# End of script.
__EOM__
])dnl End of macro to create serpa-run.in


dnl Create the configured file from its .in version.
dnl AC_CONFIG_FILES( btng_serpa_run_fn, chmod +x btng_serpa_run_fn )
dnl IMPORTANT NOTE:
dnl Using the defined symbol btng_serpa_run_fn in AC_CONFIG_FILES
dnl causes automake to complain that btng_serpa_run_fn.in is not found.
dnl Since that symbol is a defined m4 macro, this is an automake deficiency.
dnl Until this is fixed, the above AC_CONFIG_FILES usage will not
dnl work!  You must call AC_CONFIG_FILES with the actual names of
dnl the file instead of the m4-defined symbol btng_serpa_run_fn.



AC_CONFIG_COMMANDS([create-serpa-machine-file],[
# Generate a machine file if needed.
if test ! -r "$serpa_machine_file" ; then
  # Set the serpa machine file name to be in the same directory as
  # the serpa-run script.
  btng_serpa_dir_name=`echo ']btng_serpa_run_fn[' | sed -e ['s:^\([^/]*\)$:./\1:'] -e ['s:/[^/]*$::']`;
  btng_serpa_machine_file="${btng_serpa_dir_name}/serpa.machines"
  hostname="`hostname`"
  cat <<-__EOF__>"$btng_serpa_machine_file"
	$hostname
	$hostname
	$hostname
	$hostname
	__EOF__
fi
],[
dnl Settings to make before running the above.
])


# Add the serpa machine file to the clean list.
btng_serpa_dir_name=`echo ']btng_serpa_run_fn[' | sed -e ['s:^\([^/]*\)$:./\1:'] -e ['s:/[^/]*$::']`;
btng_serpa_machine_file="${btng_serpa_dir_name}/serpa.machines"
DISTCLEANFILES="$DISTCLEANFILES ${btng_serpa_machine_file}"
CASC_AC_LOG_VAR(btng_serpa_machine_file DISTCLEANFILES, in support-serpa-run)


undefine([btng_serpa_run_fn])


# End macro $0
])dnl End CASC_SUPPORT_SERPA macro.

dnl Define a macro for supporting BOOST

AC_DEFUN([CASC_SUPPORT_BOOST],[

# Begin CASC_SUPPORT_BOOST
# Defines boost_PREFIX boost_INCLUDES and boost_LIBS.
AC_ARG_WITH(boost,
[ --with-boost[=PATH]  Use BOOST and specify where BOOST is installed.],
, with_boost=no)

case "$with_boost" in
  no)
    : Do nothing
  ;;
  yes)
    # BOOST install path was not specified.
    # Look in a couple of standard locations to probe if 
    # BOOST header files are there.
    AC_MSG_CHECKING([for BOOST installation])
    for dir in /usr /usr/local; do
      if test -f ${dir}/include/boost/tr1/unordered_map.hpp; then
        boost_PREFIX=${dir}
        break
      fi
    done
    AC_MSG_RESULT([$boost_PREFIX])
  ;;
  *)
    # BOOST install path was specified.
    AC_MSG_CHECKING([for BOOST installation])

    if test -f ${with_boost}/include/boost/tr1/unordered_map.hpp; then
        boost_PREFIX=$with_boost
        boost_INCLUDES="-I${boost_PREFIX}/include"
        dnl Do not have to set boost_LIBS because SAMRAI does not need it (yet).
        dnl boost_LIBS="-L${boost_PREFIX}/lib -lboost_regex"
        AC_MSG_RESULT([$boost_PREFIX])
    else
        AC_MSG_RESULT([$boost_PREFIX])
        AC_MSG_ERROR([BOOST not found in $with_boost])
    fi
  ;;
esac


# END CASC_SUPPORT_BOOST

])dnl End definition of CASC_SUPPORT_BOOST


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

dnl Define a macro for supporting Thread Building Blocks

AC_DEFUN([SAMRAI_SUPPORT_TBB],[

# Begin SAMRAI_SUPPORT_TBB
# Defines tbb_PREFIX tbb_INCLUDES and tbb_LIBS.
AC_ARG_WITH([tbb], 
   [AS_HELP_STRING([--with-tbb], 
      [Use Thread Building Blocks])],
   [],
   [with_tbb=no])

AC_ARG_WITH([tbb-include],
   AS_HELP_STRING([--with-tbb-include=DIR], [tbb.h is in DIR]), 

      for tbb_dir in $withval; do
         TBBINCLUDE="$TBBINCLUDE -I$withval"
      done;
      have_tbb=yes
)
   
AC_ARG_WITH( [tbb-lib-dirs],
   AS_HELP_STRING([--with-tbb-lib-dirs=DIRS],
                  [DIRS is a space-separated list of directories containing the libraries specified by \'--with-tbb-libs\']),

      for tbb_lib_dir in $withval; do
         TBBLIBDIRS="-L$tbb_lib_dir $TBBLIBDIRS"
      done; 
      have_tbb=yes
)

AC_ARG_WITH([tbb-libs],
   AS_HELP_STRING([--with-tbb-libs=LIBS], 
                  [LIBS is a space-separated list of library names needed for Thread Building Blocks, e.g., \"tbb tbbmalloc\"]),

      for tbb_lib in $withval; do
         TBBLIBS="$TBBLIBS -l$tbb_lib"
      done; 
      have_tbb=yes
)

dnl
dnl check if we're disabling tbb support
dnl
if test "$with_tbb" = "no" ; then
    AC_MSG_NOTICE([configuring without Thread Building Blocks support])
else
    AC_MSG_CHECKING([for Thread Building Blocks installation])

    if test -f ${with_tbb_include}/tbb.h; then
        AC_MSG_RESULT([$TBBINCLUDE])
    else
        AC_MSG_RESULT([$TBBINCLUDE])
        AC_MSG_ERROR([TBB not found in $with_tbb_include])
    fi

    if test -f ${with_tbb_lib_dirs}/libtbb.so; then
        AC_MSG_RESULT([$TBBLIBDIRS])
    else
        AC_MSG_RESULT([$TBBLIBDIRS])
        AC_MSG_ERROR([TBB libs not found in $with_tbb_lib_dirs])
    fi
fi

dnl Test compiling an TBB application
dnl
dnl NOTE that AC_SEARCH_LIBS didn't work completely so
dnl use a more complicated example program to see
dnl if that will catch when TBB is not working.  Using one of the TBB
dnl short tests for this.

if test "${TBBINCLUDE+set}" = set; then

   AC_REQUIRE([AC_PROG_CXX])
   AC_MSG_CHECKING(whether Thread Building Blocks link works)
   AC_LANG_PUSH(C++)
   CASC_PUSH_COMPILER_STATE
   LIBS="${LIBS} ${TBBLIBDIRS} ${TBBLIBS}"
   CXXFLAGS="${CXXFLAGS} ${TBBINCLUDE}"
   AC_LINK_IFELSE([
/*
    Copyright 2005-2011 Intel Corporation.  All Rights Reserved.

    This file is part of Threading Building Blocks.

    Threading Building Blocks is free software; you can redistribute it
    and/or modify it under the terms of the GNU General Public License
    version 2 as published by the Free Software Foundation.

    Threading Building Blocks is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Threading Building Blocks; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

    As a special exception, you may use this file as part of a free software
    library without restriction.  Specifically, if other files instantiate
    templates or use macros or inline functions from this file, or you compile
    this file and link it with other files to produce an executable, this
    file does not by itself cause the resulting executable to be covered by
    the GNU General Public License.  This exception does not however
    invalidate any other reasons why the executable file might be covered by
    the GNU General Public License.
*/

#include <iostream>
#include <string>
#include <algorithm>
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

using namespace tbb;
using namespace std;
static const size_t N = 23;

class SubStringFinder {
  const string str;
  size_t *max_array;
  size_t *pos_array;
public: 
  void operator() ( const blocked_range<size_t>& r ) const { 
    for ( size_t i = r.begin(); i != r.end(); ++i ) {
      size_t max_size = 0, max_pos = 0;
      for (size_t j = 0; j < str.size(); ++j)
      if (j != i) {
        size_t limit = str.size()-max(i,j);
        for (size_t k = 0; k < limit; ++k) {
          if (str[[i + k]] != str[[j + k]]) break;
          if (k > max_size) {
            max_size = k;
            max_pos = j;
          }
        }
      }
      max_array[[i]] = max_size;
      pos_array[[i]] = max_pos;
    }
  }
  SubStringFinder(string &s, size_t *m, size_t *p) :
    str(s), max_array(m), pos_array(p) { }
};

int main() {

  string str[[N]] = { string("a"), string("b") };
  for (size_t i = 2; i < N; ++i) str[[i]] = str[[i-1]]+str[[i-2]];
  string &to_scan = str[[N-1]]; 
  size_t num_elem = to_scan.size();

  size_t *max = new size_t[[num_elem]];
  size_t *pos = new size_t[[num_elem]];

  parallel_for(blocked_range<size_t>(0, num_elem ),
               SubStringFinder( to_scan, max, pos ) );

  for (size_t i = 0; i < num_elem; ++i)
    cout << " " << max[[i]] << "(" << pos[[i]] << ")" << endl;
  delete[[]] pos;
  delete[[]] max;
  return 0;
}
      ], 
      samrai_tbb_compile=yes,
      samrai_tbb_compile=no)
   CASC_POP_COMPILER_STATE
   AC_LANG_POP
   AC_MSG_RESULT($samrai_tbb_compile)

   if test "$samrai_tbb_compile" = no; then
      AC_MSG_ERROR([TBB compile/link test failed])
   fi
fi

# END SAMRAI_SUPPORT_TBB

])dnl End definition of SAMRAI_SUPPORT_TBB


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

AC_DEFUN([SAMRAI_SPLIT_LIBS_STRING],[
dnl
dnl Macro SAMRAI_SPLIT_LIBS_STRING
dnl Written by Brian Gunney.
dnl
dnl This macro takes an automake-style LIBS string (arg1) and
dnl splits it into the -L part (arg2, what SAMRAI usually calls
dnl LIB_PATH) and -l part (arg3, what SAMRAI calls LIB_NAME).
dnl The rest are also lumped into the LIB_NAME part, for lack
dnl of generality in the SAMRAI distinction.
dnl
dnl I think the SAMRAI format is limitting because not all parts
dnl if the LIBS strings can be recategorized into LIB_PATH and
dnl LIB_NAME parts.  But this macro allows us to use macros that
dnl conform to the standard automake format.  BTNG.
dnl
# Split $1 into the LIB_PATH part ($2) and the LIB_NAME part ($3)
if test -n "${$1}"; then
  for i in ${$1}; do
    case "$i" in
    -L*) $2="${$2} $i" ;;
    *) $3="${$3} $i" ;;
    esac
  done
fi
])

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

dnl Define macros for supporting HDF5.
dnl $Id$

AC_DEFUN([CASC_VAR_SET_DL],[
dnl Provides support for the dl (dynamic loading) library.
dnl
dnl Arguments are:
dnl 1. Name of variable to set to path where dl are installed.
dnl    Nothing is done if this variable is unset.
dnl    If you only want to look in default locations, set it to blank.
dnl 2. Name of the INCLUDES variable similar to the automake INCLUDES variable.
dnl    This variable is modified ONLY if it is NOT set and the path
dnl    is non-blank.
dnl 3. Name of the LIBS variable similar to the automake LIBS variable.
dnl    This variable is modified ONLY if it is NOT set.
dnl    If the library cannot be found, this remains unset.
dnl
dnl If arg1 is defined, assume that the user wants dl support.
dnl Do so by assigning arg2 and arg3 if they are not defined.
dnl
# Begin macro CASC_VAR_SET_DL
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
    AC_SEARCH_LIBS([dlopen],dl,[
      CASC_AC_LOG_VAR(LIBS,After finding dl flag)
      # Action if found ...
      # Extract modifications to LIB into library-specific LIBS variable.
      $3=`echo " $LIBS" | sed "s! $btng_save_LIBS!!"`;
      test -n "${$1}" && $3="-L${$1}/lib ${$3}"
      CASC_AC_LOG_VAR($3, Found dl library flag)
      ],[
      # Action if NOT found ...
      CASC_AC_LOG_VAR($3, Did not find dl library flag)
      ],[$btng_extra_libs])
    LIBS="$btng_save_LIBS";	# Restore global-use variable.
    unset btng_extra_libs
    unset btng_save_LIBS
  else
    CASC_AC_LOG(Not looking for dl because $3 is already set)
  fi
fi
# End macro CASC_VAR_SET_DL
])dnl



AC_DEFUN([CASC_SUPPORT_DL],[
dnl Support dl library by setting the variables
dnl dl_PREFIX, dl_INCLUDES, and dl_LIBS.
dnl Arg1: non-empty if you want the default to be on.
dnl
# Begin macro CASC_SUPPORT_DL

CASC_ARG_WITH_ENV_WRAPPER(dl, dl_PREFIX,
ifelse($1,,
[  --with-dl[=PATH]
			Use the dynamic loading library and optionally
			specify where it is installed.],
[  --without-dl		Do not use the dynamic loading library.]),
ifelse($1,,unset dl_PREFIX; test "${with_dl+set}" = set && dl_PREFIX=,dl_PREFIX=))

CASC_ARG_WITH_PREFIX(dl-includes,dl_INCLUDES,
[  --with-dl-includes=STRING
			Specify the INCLUDES flags for dl.
			If not specified, and --with-dl=PATH is,
			this defaults to "-IPATH/include".])dnl

CASC_ARG_WITH_PREFIX(dl-libs,dl_LIBS,
[  --with-dl-libs=STRING
			Specify LIBS flags for dl.
			If not specified, and --with-dl=PATH is,
			this defaults to "-LPATH/lib -ldl".])dnl

CASC_VAR_SET_DL(dl_PREFIX,dl_INCLUDES,dl_LIBS)
# End macro CASC_SUPPORT_DL
])

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

dnl $Id$

AC_DEFUN([CASC_VAR_SET_LAPACK],[
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
      CASC_AC_LOG_VAR(LIBS,After finding lapack flag)
      # Action if found ...
      # Extract modifications to LIB into library-specific LIBS variable.
      $3=`echo " $LIBS" | sed "s! $btng_save_LIBS!!"`;
      test -n "${$1}" && $3="-L${$1}/lib ${$3}"
      CASC_AC_LOG_VAR($3, Found lapack library flag)
      ],[
      # Action if NOT found ...
      CASC_AC_LOG_VAR($3, Did not find lapack library flag)
      AC_MSG_WARN(
[I could not systematically find the name of
the lapack library so I am using -llapack instead.])
      $3="-llapack"
      test -n "${$1}" &&	\
	$3="-L${$1}/lib ${$3}"	# Add path flag to output variable.
      ],[$btng_extra_libs])
    LIBS="$btng_save_LIBS";	# Restore global-use variable.
    unset btng_extra_libs
    unset btng_save_LIBS
  else
    CASC_AC_LOG(Not looking for lapack because $3 is already set)
  fi
fi
])dnl




AC_DEFUN([CASC_SUPPORT_LAPACK],[
dnl Support lapack library by setting the variables
dnl lapack_PREFIX, lapack_INCLUDES, and lapack_LIBS.
dnl Arg1: non-empty if you want the default to be on.
dnl
# Begin macro CASC_SUPPORT_LAPACK

CASC_ARG_WITH_ENV_WRAPPER(lapack, lapack_PREFIX,
ifelse($1,,
[  --with-lapack[=PATH]
			Use lapack and optionally specify where
			they are installed.],
[  --without-lapack	Do not use the lapack library.]),
if test "${with_lapack+set}" = set; then
  lapack_PREFIX=
else
ifelse($1,,unset lapack_PREFIX,lapack_PREFIX=)
fi
)

CASC_ARG_WITH_PREFIX(lapack-includes,lapack_INCLUDES,
[  --with-lapack-includes=STRING
			Specify the INCLUDES flags for lapack.
			If not specified, and --with-lapack=PATH is,
			this defaults to "-IPATH/include".])dnl

CASC_ARG_WITH_PREFIX(lapack-libs,lapack_LIBS,
[  --with-lapack-libs=STRING
			Specify LIBS flags for lapack.
			If not specified, and --with-lapack=PATH is,
			this defaults to "-LPATH/lib -llapack".])dnl

CASC_VAR_SET_LAPACK(lapack_PREFIX,lapack_INCLUDES,lapack_LIBS)

CASC_AC_LOG_VAR(lapack_PREFIX lapack_INCLUDES lapack_LIBS)
# End macro CASC_SUPPORT_LAPACK
])

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

dnl $Id$

AC_DEFUN([CASC_SUPPORT_PETSC],[
# Begin macro SUPPORT_PETSC
dnl Support PETSC by setting PETSC_DIR, PETSC_ARCH,
dnl petsc_INCLUDE and petsc_LIBS.
dnl Also set PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR and
dnl PETSC_VERSION_SUBMINOR to indicate PETSc version.
dnl
dnl Support --with-petsc-optimize to use optimized PETSC library.
dnl Support --with-petsc-mpiuni to use PETSC uniprocessor MPI library.
dnl
dnl Arg1: non-empty if you want the default to be on.
dnl
dnl This version supports PETSc-2.1.0 and later.

# Set PETSC_DIR to the PETSC root directory.
CASC_ARG_WITH_PREFIX(petsc,PETSC_DIR,
ifelse($1,,
[[  --with-petsc=PATH	Support PETSc, and specify PETSC top-level directory.
			Setting PETSC_DIR is equivalent to this.]],
[  --without-petsc	Do not support PETSc.])
,
# User was not specific about specifying PETSc.
ifelse($1,,
# PETSc should be off by default.
# So if user specified --with-petsc, it should include a path or it is an error.
if test "${with_petsc+set}" = set ; then
  # User specified --with-petsc ambiguously.
  AC_MSG_ERROR([You must specify a path with --with-petsc=...])
else
  # User did not specify --with-petsc, so turn it off.
  unset PETSC_DIR;
fi
,
# PETSc should be on by default.
# So require that user say where it is or turn it off.
AC_MSG_ERROR([You must specify either --with-petsc=... or --without-petsc.])
))

# Set version numbers (PETSC_VERSION_...) for use by configure.
# Users can directly access these from the PETSc header file.
if test "${PETSC_DIR+set}" = set; then
[
PETSC_VERSION_MAJOR=`sed -e '/^[ \t]\{0,\}#define PETSC_VERSION_MAJOR/!d' -e 's/.\{0,\}[ \t]\{1,\}//' $PETSC_DIR/include/petscversion.h`
PETSC_VERSION_MINOR=`sed -e '/^[ \t]\{0,\}#define PETSC_VERSION_MINOR/!d' -e 's/.\{0,\}[ \t]\{1,\}//' $PETSC_DIR/include/petscversion.h`
PETSC_VERSION_SUBMINOR=`sed -e '/^[ \t]\{0,\}#define PETSC_VERSION_SUBMINOR/!d' -e 's/.\{0,\}[ \t]\{1,\}//' $PETSC_DIR/include/petscversion.h`
]
fi
CASC_AC_LOG_VAR(PETSC_VERSION_MAJOR PETSC_VERSION_MINOR PETSC_VERSION_SUBMINOR)

# Set PETSC_ARCH.
CASC_ARG_WITH_PREFIX(petsc-arch,PETSC_ARCH,
[  --with-petsc-arch=PETSC_ARCH
			Specify the PETSC architecture.
			If omitted, the output of the petscarch script
			in the PETSc directory is used.])

# Set PETSC_OPTIMIZE.
CASC_ARG_WITH_ENV_WRAPPER(petsc-optimize,PETSC_OPTIMIZE,
[  --with-petsc-optimize
			Use the optimized PETSC libraries],
# By default, use the debug PETSC library.
PETSC_OPTIMIZE=g
CASC_AC_LOG_VAR(with_petsc_optimize)
test "$with_petsc_optimize" = yes && PETSC_OPTIMIZE=O
)

# Set PETSC_MPIUNI.
CASC_ARG_WITH_ENV_WRAPPER(petsc-mpiuni,PETSC_MPIUNI,
[  --with-petsc-mpiuni	Use the PETSC uniprocessor MPI library],
# By default, do not use the PETSC uniprocessor MPI library.
unset PETSC_MPIUNI
)

# Set PETSC_LIBFILES
unset PETSC_LIBFILES
AC_ARG_WITH(petsc-libfiles,
[  --with-petsc-libfiles
			Specify explit PETSc library files instead
			of -L and -l flags (may help some debuggers)],
test "$with_petsc_libfiles" = yes && PETSC_LIBFILES=yes)


if test "${PETSC_DIR+set}" = set; then
  # Set up PETSC only if PETSC_DIR is defined.

  if test ! -d "$PETSC_DIR"; then
    AC_MSG_WARN([PETSC directory ($PETSC_DIR) does not look right])
  fi
  export PETSC_DIR
  if test -z "$PETSC_ARCH"; then
    if test -f "$PETSC_DIR/bmake/petscconf"; then
	eval `grep PETSC_ARCH $PETSC_DIR/bmake/petscconf`
    elif test -x "$PETSC_DIR/bin/petscarch"; then
       PETSC_ARCH=`$PETSC_DIR/bin/petscarch`
    else
       AC_MSG_WARN([PETSC could not determine PETSC_ARCH])
    fi    
    export PETSC_ARCH
  fi
  CASC_AC_LOG_VAR(PETSC_ARCH)
  if test ! -d "$PETSC_DIR/bmake/$PETSC_ARCH"; then
    AC_MSG_WARN([PETSC architecture ($PETSC_ARCH) does not look right])
  fi
  if test ! "$PETSC_OPTIMIZE" = g && test ! "$PETSC_OPTIMIZE" = O; then
    AC_MSG_ERROR([PETSC optimize should be either g or O])
  fi

  petsc_INCLUDES="-I$PETSC_DIR/include -I$PETSC_DIR/bmake/$PETSC_ARCH"
  petsc_INCLUDES="$petsc_INCLUDES -I$PETSC_DIR/src/vec"
  # Currently, I'm not entirely sure why we have to explicitly specify
  # the src/vec directory in the include path.  But there is at least
  # one required file there that cannot be found in the include directory.

# SGS Support latter version of PETSc
# Try new structure and then old
  if test -d "${PETSC_DIR}/lib/${PETSC_ARCH}"; then
    petsc_LIBDIR="${PETSC_DIR}/lib/${PETSC_ARCH}"
  elif  test -d "${PETSC_DIR}/lib/lib${PETSC_OPTIMIZE}/${PETSC_ARCH}"; then
    petsc_LIBDIR="${PETSC_DIR}/lib/lib${PETSC_OPTIMIZE}/${PETSC_ARCH}"
  else 
    AC_MSG_WARN([PETSC lib directory does not look as expected])
  fi


  # Issue the -L flag if not specifying PETSc libraries by file names.
  test ! "${PETSC_LIBFILES+set}" = set && petsc_LIBS="-L${petsc_LIBDIR}"

  # Build up a list of PETSC library files.
# SGS
  petsc_libs_ls1=`cd ${petsc_LIBDIR} && echo lib*.*`

  if test -n "$petsc_libs_ls1"; then
    unset petsc_libs_ls
    for i in $petsc_libs_ls1; do
      j=`echo $i | sed -e 's/lib//' -e 's/\.a$//' -e 's/\.so$//'`
      if echo "$petsc_libs_ls" | grep -v " $j " > /dev/null; then # Note padding!
        petsc_libs_ls="$petsc_libs_ls $j ";	# Note space padding!
      fi
    done
  fi
  # Remove mpiuni from the list of PETSC libraries unless user asked for it.
  if test ! "${PETSC_MPIUNI}" = yes; then
    petsc_libs_ls=`echo "$petsc_libs_ls" | sed 's/ mpiuni //g'`
  fi
  # Move some low-level libraries to the end to ensure resolution
  # for linkers that only make one pass.
  for i in petscmat petscvec petsc; do
    petsc_libs_ls=`echo "$petsc_libs_ls" | sed 's/\(.*\)\( \{0,1\}'"$i"'\{0,1\} \)\(.*\)/\1 \3 \2/g'`
  done
  # Build up petsc_LIBS string using library names.
  CASC_AC_LOG_VAR(petsc_libs_ls1 petsc_libs_ls)
  if test -n "$petsc_libs_ls"; then
    if test "${PETSC_LIBFILES+set}" = set; then
      for i in $petsc_libs_ls; do
# SGS
        petsc_LIBS="$petsc_LIBS ${petsc_LIBDIR}/lib${i}.a"
      done
    else
      for i in $petsc_libs_ls; do
        petsc_LIBS="$petsc_LIBS -l$i"
      done
    fi
  fi

  CASC_AC_LOG_VAR(PETSC_DIR petsc_INCLUDES petsc_LIBS PETSC_OPTIMIZE PETSC_MPIUNI)

fi
# End macro SUPPORT_PETSC
])

dnl $Id$



AC_DEFUN([CASC_C_RESTRICT],[

# Start macro CASC_C_RESTRICT

AC_MSG_CHECKING(checking whether restrict is broken)

AC_CACHE_VAL(btng_cv_c_restrict_broken, [

  AC_LANG_PUSH([C++])
  AC_TRY_COMPILE([
struct array_test {
  double *ptr;
  int i0;
  array_test(double *p, int i);
  double &value(int i) const;
};
array_test::array_test(double *p, int i) : ptr(p), i0(i) {}
double &array_test::value(int i) __restrict__ const {
  return ptr[i-i0];
}
    ],[
double a[10];
array_test at(a,20);
at.value(5) = 5;
    ],
    # restrict is not broken.
    btng_cv_c_restrict_broken=no
    ,
    # restrict is broken.
    btng_cv_c_restrict_broken=yes
  )	dnl End AC_TRY_COMPILE call

  AC_LANG_POP([C++])

])	dnl End AC_CACHE_VAL call

AC_MSG_RESULT($btng_cv_c_restrict_broken)

if test "$btng_cv_c_restrict_broken" = yes; then
  AC_DEFINE([RESTRICT_IS_BROKEN],1,Define if restrict is not properly supported)
fi


# End macro CASC_C_RESTRICT

])	dnl End of CASC_C_RESTRICT definition.

dnl $Id$


AC_DEFUN([CASC_LIBS_ADD_RPATH],[
# Begin macro CASC_LIBS_ADD_RPATH
dnl Support RPATH by going in a LIBS string and, for each -L flag,
dnl add a flag immediately following it to set the RPATH, for
dnl paths that contain shared libraries.
dnl
dnl arg1 is a LIBS string.
dnl arg2 is the name of the variable to set to the new LIBS string.
dnl arg3 is non-empty to use id of the C++ compiler instead of the C compiler.


dnl Determine which compiler is being used, because
dnl the syntax of the RPATH flag depends on the compiler.
dnl Use the C++ compiler and assume the C compiler
dnl is from the same family.
AC_REQUIRE([CASC_INFO_CC_CXX_ID])


AC_ARG_ENABLE(rpath,
[  --enable-rpath=SYNTAX	When linking add syntax for rpath for every
			-L option that points to a directory with .so
			files in it.  If SYNTAX is omitted, an attempt
			is made to find out the correct rpath syntax for
			the compiler being used.]
,,enable_rpath=yes)

if test "$enable_rpath" = yes; then
  # Determine the proper rpath syntax.

  AC_LANG_SAVE

  ifelse([$3],,
  AC_LANG_C
  btng_rpath_compiler_id="$CC_ID",
  AC_LANG_CPLUSPLUS
  btng_rpath_compiler_id="$CXX_ID"
  )


  # Unset the rpath syntax variable so we can check on whether we
  # found a way to set it.
  unset btng_rpath_beginning;

  # Determine, based on the compiler, the syntax for specifying RPATH.
  # It should be of the form "$btng_rpath_beginning$the_path", where
  # btng_rpath_beginning is the compiler-dependent part.
  case "$btng_rpath_compiler_id" in
    gnu)
      # This compiler may use a variable rpath syntax because it may use
      # the native loader.
      CASC_LIBS_FIND_RPATH(btng_rpath_beginning,
	['---bogus-flag-meant-to-cause-error' '-Wl,-rpath ' '-Wl,-R' '-Wl,-R '])
    ;;
    intel)
      # This compiler may use a variable rpath syntax because it may use
      # the native loader.
      CASC_LIBS_FIND_RPATH(btng_rpath_beginning,
	['---bogus-flag-meant-to-cause-error' '-Wl,-rpath ' '-Wl,-R' '-Wl,-R '])
      if test "$btng_rpath_beginning" = "---bogus-flag-meant-to-cause-error"; then
        # Do not rely on the compiler return value to test for syntax
        # Guess the syntax assuming the native loader will be used.
        case "$host_os" in
          linux*) btng_rpath_beginning='-Wl,-rpath ' ;;
          sun*|solaris*) btng_rpath_beginning='-R' ;;
          osf*) btng_rpath_beginning='-rpath ' ;;
          *) btng_rpath_beginning='' ;;
        esac
        AC_MSG_WARN(
  [Your compiler ifelse($3,,$CC,$CXX) returns 0 even when it is
  given a bogus flag.  Therefore, I cannot find the proper syntax
  for the rpath for this compiler.  I have resorted to a guess that
  may not be correct: '$btng_rpath_beginning'.
  You can override this by using --enable-rpath=SYNTAX])
      fi
    ;;
    sunpro)
      # This compiler may use a variable rpath syntax.
      CASC_LIBS_FIND_RPATH(btng_rpath_beginning,['---bogus-flag-meant-to-cause-error' '-R' '-R '])
    ;;
    kai)
      # The KAI compilers use the system native loader.
      #
      # On some platforms (PC/Linux at least), this compiler seems
      # to return 0 even if it encounters error, thus it can return
      # the first guess for the rpath syntax, even if the guess is
      # wrong.  We try to catch this by making the first flag bogus.
      # If the compiler accepts this flag (by returning 0), we know
      # it is wrong and we resort to an alternative method for
      # getting the rpath syntax.
      CASC_LIBS_FIND_RPATH(btng_rpath_beginning,
	['---bogus-flag-meant-to-cause-error' '-R' '-R ' '-rpath ' '-Wl,-rpath ' '-Wl,-R' '-Wl,-R '])
      if test "$btng_rpath_beginning" = "---bogus-flag-meant-to-cause-error"; then
        # Do not rely on the compiler return value to test for syntax
        # Guess the syntax assuming the native loader will be used.
        case "$host_os" in
          linux*) btng_rpath_beginning='-Wl,-rpath ' ;;
          sun*|solaris*) btng_rpath_beginning='-R' ;;
          osf*) btng_rpath_beginning='-rpath ' ;;
          *) btng_rpath_beginning='' ;;
        esac
        AC_MSG_WARN(
  [Your compiler ifelse($3,,$CC,$CXX) returns 0 even when it is
  given a bogus flag.  Therefore, I cannot find the proper syntax
  for the rpath for this compiler.  I have resorted to a guess that
  may not be correct: '$btng_rpath_beginning'.
  You can override this by using --enable-rpath=SYNTAX])
      fi
    ;;
    *)
      CASC_LIBS_FIND_RPATH(btng_rpath_beginning)
    ;;
  esac
  CASC_AC_LOG_VAR(host_os CC_ID CXX_ID btng_rpath_compiler_id btng_rpath_beginning, forming rpaths)

  AC_LANG_RESTORE

  # It is valid to have btng_rpath_beginning be blank.
  # but if it is unset, we could not find a way to set it.
  if test ! "${btng_rpath_beginning+set}" = set; then
    AC_MSG_WARN(I cannot find a working syntax for setting relocatable paths)
  fi

elif test ! "${enable_rpath}" = no; then

  # User has provided the rpath syntax.
  btng_rpath_beginning=$enable_rpath

fi;	# End block determining the proper rpath syntax.


# Use the rpath syntax.
if test "${btng_rpath_beginning+set}" = set	\
  && test -n "${btng_rpath_beginning}" ; then
  # Add the RPATH flags only if we know the syntax for it,
  # and if it is needed as indicated by a non-empty btng_rpath_beginning.

  # Loop through the flags in $1, looking for the -L flag,
  # and append RPATH flag to each one found, if the the
  # path specified by the flag includes shared libraries.
  for i in ${$1}; do
    btng_new_$2="${btng_new_$2} ${i}"
    btng_tmp_addl_string=`echo $i | sed 's/^-L//'`
    test "$btng_tmp_addl_string" = "$i" && continue	# does not contain -L.
    test -d "$btng_tmp_addl_string" || continue;	# directory nonexistent.
    test "`echo $btng_tmp_addl_string/*.so`" = "$btng_tmp_addl_string/*.so" \
      && continue;	# does not contain shared libraries.
    echo "${btng_new_$2}"	\
      | grep ".*${btng_rpath_beginning}[[ 	]]*${btng_tmp_addl_string}"	\
      > /dev/null	\
      && continue	# already contains the flag we want to add.
    btng_new_$2="${btng_new_$2} ${btng_rpath_beginning}${btng_tmp_addl_string}"
  done
  $2="${btng_new_$2}"

fi

dnl Now, arg2 should be similar to arg1, but with the additional RPATH flags.

# End macro CASC_LIBS_ADD_RPATH
])

AC_DEFUN([CASC_LIBS_FIND_RPATH],[
# Begin macro CASC_LIBS_FIND_RPATH
dnl Find the correct rpath syntax from the list given in arg1.
dnl arg1: variable to set to the syntax string
dnl arg2: list of syntaxes to try;
dnl   if blank, a large number of syntaxes will be tried.
dnl
dnl arg1 is list of possible rpath syntaxes to try.
define(btng_possible_rpaths,dnl
[ifelse($2,,['-R ' '-R' '-rpath ' '-Wl,-rpath ' '-Wl,-R ' '-Wl,-R'],[[$2]])])
  btng_save_LIBS="$LIBS";
  for i in btng_possible_rpaths; do
    LIBS="${i}/usr/local"
    AC_TRY_LINK(,,$1="$i", unset $1)
    # Intel compiler does not fail on bad args but warning message is
    # created. If warning is found in the log then continue searching
    # for syntax as the current one is no good.  If warning is not
    # found use return status of link attempt to determine if
    # parameter was accepted by the compiler.
    SEARCH=`echo "ignoring unknown option '${LIBS}'" | sed -e "s/---/-f-/"`
    if ( grep "$SEARCH" config.log ) >/dev/null 2>&1
    then 
	:
    else 
        if test "${$1+set}" = set; then break; fi
    fi
  done
  LIBS="$btng_save_LIBS"
undefine([btng_possible_rpaths])
# End macro CASC_LIBS_FIND_RPATH
])

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

dnl $Id$

AC_DEFUN([CASC_FIND_CORRECT_HEADER_FILENAME],[
dnl There is no standard naming convention for STL header files.
dnl This macro helps to pick the right name out of a list.
dnl Arg1 is the variable to set to the found file name.
dnl Arg2 is the list of file names to search
dnl Arg3 are additional headers to include (for use by AC_TRY_COMPILE)
dnl Arg4 is the code body to test if the included file works.
# Start macro $0
  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  $1=
  AC_REQUIRE([CASC_TYPE_BOOL])
  CPPFLAGS_SAVE=$CPPFLAGS
  for file in $2; do
    AC_CHECK_HEADER($file, btng_header_found=1, unset btng_header_found)
    if test -n "$btng_header_found"; then
      AC_MSG_CHECKING(whether $file is the header sought)
      CASC_AC_LOG(found header file $file)
      CPPFLAGS="$CPPFLAGS_SAVE $CXX_OPTIONS"
      AC_TRY_COMPILE(
        [
/* macro $0 checking for $file */
#ifdef BOOL_IS_BROKEN
typedef int bool;
#define true 1
#define false 0
#endif
	$3
        #include <$file>
using namespace std;
],
        $4,
	AC_MSG_RESULT(yes)
        $1="$file",
	AC_MSG_RESULT(no)
      )
    fi
    if test -n "${$1}"; then break; fi
  done
  AC_LANG_RESTORE
  CPPFLAGS=$CPPFLAGS_SAVE
# End macro $0
])




AC_DEFUN([CASC_TREAT_VARIABLE_HEADER_FILENAME],[
dnl CASC_TREAT_VARIABLE_HEADER_FILENAME is a generic macro
dnl used by (and using) other macros in this file.
dnl It determines, from a given list, the correct name of
dnl a header file required to compile a test code body.
dnl It takes a list of possible of the header filenames.
dnl It reports whether each header file is the one sought
dnl until it finds the one that is.
dnl If none of the header filenames work:
dnl   It issues a warning.
dnl   It defines a ...IS_BROKEN C macro saying so.
dnl If it finds the first header filename that works:
dnl   It assigns a variable (..._HEADER_FILE) to the
dnl   correct filename and call AC_DEFINE for that variable.
dnl Arguments are:
dnl  1: a single name representing the header sought.
dnl  2: a list of possible header filenames.
dnl  3: other include lines (for use in AC_TRY_COMPILE).
dnl  4: code to test if the header file is the one being sought.
dnl
# Start macro $0
AC_CACHE_VAL(btng_cv_[]translit($1,[-],[_])[]_header_filename, [
  AC_ARG_WITH($1-header-file,
  [  --with-$1-header-file	Specify name of the $1 header file.],
  btng_cv_[]translit($1,[-],[_])[]_header_filename=$with_[]translit($1,[-],[_])[]_header_file,
  [CASC_FIND_CORRECT_HEADER_FILENAME(btng_cv_[]translit($1,[-],[_])[]_header_filename,$2,[$3],[[$4]])]
  )
])	dnl End AC_CACHE_VAL call
# We must be able to find the $1 header file or else.
translit($1,[-a-z],[_A-Z])[]_HEADER_FILE="$btng_cv_[]translit($1,[-],[_])[]_header_filename"
if test -z "$translit($1,[-a-z],[_A-Z])[]_HEADER_FILE"; then
  translit($1,[-],[_])[]_header_is_broken=1
  AC_MSG_WARN([cannot find a working $1 header file.
      Names tried: $2
      If you know the correct hame of this header file,
      use the option --with-[]$1[]-header-file=FILENAME
      with configure.])
  AC_DEFINE(translit($1,[-a-z],[_A-Z])[]_IS_BROKEN,[1],[The $1 header file is broken])
  CASC_AC_LOG(header file $1 is broken)
else
  unset translit($1,[-],[_])[]_header_is_broken
  AC_DEFINE_UNQUOTED(translit($1,[-a-z],[_A-Z])[]_HEADER_FILE,<$translit($1,[-a-z],[_A-Z])[]_HEADER_FILE>,
    [Header file for $1])
  CASC_AC_LOG(header file $1 is ok)
fi
# End macro $0
])	dnl end of CASC_TREAT_VARIABLE_HEADER_FILENAME definition.



dnl
dnl These are some STL headers with uncertain names.
dnl


AC_DEFUN([CASC_STL_STRING_HEADER_FILENAME],[
# Start macro $0
dnl dnl AC_MSG_CHECKING(name of the STL string header file)
CASC_TREAT_VARIABLE_HEADER_FILENAME([stl-string],
  [string strings string.h strings.h string.hxx strings.hxx],,
  [std::string s; s = "sample string";])
# End macro $0
])	dnl end of CASC_STL_STRING_HEADER_FILENAME definition.


AC_DEFUN([CASC_STL_SET_HEADER_FILENAME],[
# Start macro $0
dnl AC_MSG_CHECKING(name of the STL set header file)
CASC_TREAT_VARIABLE_HEADER_FILENAME([stl-set], [set set.h set.hxx],,
  [set<int> s; s.insert(1);])
# End macro $0
])	dnl end of CASC_STL_SET_HEADER_FILENAME definition.


AC_DEFUN([CASC_STL_STACK_HEADER_FILENAME],[
# Start macro $0
dnl AC_MSG_CHECKING(name of the STL stack header file)
AC_REQUIRE([CASC_INFO_CXX_ID])
AC_REQUIRE([CASC_STL_LIST_HEADER_FILENAME])
btng_stl_stack_test_body='[stack<int> s; s.push(1);]'
# The Sun compiler version 5.2 does not treat default template
# arguments correctly.  The STL standard states that for stack,
# only the first argument is required but this Sun compiler
# requires the second.
if test "$CXX_ID" = "sunpro" && echo "$CXX_VERSION" | grep '^0x420' > /dev/null || test "$CXX_ID" = "sunpro" && echo "$CXX_VERSION" | grep '^0x520' > /dev/null ; then
btng_stl_stack_test_body='[stack<int,list<int> > s; s.push(1);]'
fi
CASC_TREAT_VARIABLE_HEADER_FILENAME([stl-stack], [stack stack.h stack.hxx],,
  [$btng_stl_stack_test_body])
# End macro $0
])	dnl end of CASC_STL_STACK_HEADER_FILENAME definition.


AC_DEFUN([CASC_STL_VECTOR_HEADER_FILENAME],[
# Start macro $0
CASC_TREAT_VARIABLE_HEADER_FILENAME([stl-vector], [vector vector.h vector.hxx],,
[vector<int> v; v.insert(v.begin(),1);
vector<char> s; s.insert( s.end(), 10, '\0' );])
# End macro $0
])	dnl end of CASC_STL_VECTOR_HEADER_FILENAME definition.


AC_DEFUN([CASC_STL_LIST_HEADER_FILENAME],[
# Start macro $0
dnl AC_MSG_CHECKING(name of the STL list header file)
CASC_TREAT_VARIABLE_HEADER_FILENAME([stl-list], [list list.h list.hxx],,
  [list<int> v; v.insert(v.begin(),1);])
# End macro $0
])	dnl end of CASC_STL_LIST_HEADER_FILENAME definition.


AC_DEFUN([CASC_STL_MAP_HEADER_FILENAME],[
# Start macro $0
dnl AC_MSG_CHECKING(name of the STL map header file)
AC_REQUIRE([CASC_INFO_CXX_ID])
btng_stl_map_test_body='[map<int,int> v; v[0]=1;]'
# The Sun compiler version 4.2 does not treat default template
# arguments correctly.  The STL standard states that for map,
# only the first two arguments are required but the Sun compiler
# requires the third.
test "$CXX_ID" = "sunpro" && echo "$CXX_VERSION" | grep '^0x420' > /dev/null && \
btng_stl_map_test_body='[map<int,int,less<int> > v; v[0]=1;]'
CASC_TREAT_VARIABLE_HEADER_FILENAME([stl-map], [map map.h map.hxx],,
  [$btng_stl_map_test_body])
# End macro $0
])	dnl end of CASC_STL_MAP_HEADER_FILENAME definition.


AC_DEFUN([CASC_STL_ITERATOR_HEADER_FILENAME],[
# Start macro $0
dnl AC_MSG_CHECKING(name of the STL iterator header file)
CASC_TREAT_VARIABLE_HEADER_FILENAME([stl-iterator],
  [iterator iterator.h iterator.hxx],,
  [int a[10], size; size=distance(a,a+10);])
dnl  [ostream_iterator<int> v(cout," ");])
# End macro $0
])	dnl end of CASC_STL_ITERATOR_HEADER_FILENAME definition.


AC_DEFUN([CASC_STL_ALGO_HEADER_FILENAME],[
# Start macro $0
dnl AC_MSG_CHECKING(name of the STL algo header file)
CASC_TREAT_VARIABLE_HEADER_FILENAME([stl-algo],
  [algo algorithm algo.h algorithm.h algo.hxx algorithm.hxx] ,,
  [int n[10]; find(n,n+10,0);])
# End macro $0
])	dnl end of CASC_STL_ALGO_HEADER_FILENAME definition.


AC_DEFUN([CASC_STL_FUNCTION_HEADER_FILENAME],[
# Start macro $0
dnl AC_MSG_CHECKING(name of the STL numeric header file)
CASC_TREAT_VARIABLE_HEADER_FILENAME([stl-function],
  [function function.h function.hxx] ,,
  [int a=1, b=2, c; plus<int> adder; c=adder(a,b);])
# End macro $0
])	dnl end of CASC_STL_FUNCTION_HEADER_FILENAME definition.


AC_DEFUN([CASC_STL_NUMERIC_HEADER_FILENAME],[
# Start macro $0
dnl AC_MSG_CHECKING(name of the STL numeric header file)
CASC_TREAT_VARIABLE_HEADER_FILENAME([stl-numeric],
  [numeric numeric.h numeric.hxx] ,,
  [int n[10]; iota(n,n+10,0);])
# End macro $0
])	dnl end of CASC_STL_NUMERIC_HEADER_FILENAME definition.


AC_DEFUN([CASC_STL_SSTREAM_HEADER_FILENAME],[
# Start macro $0
dnl AC_MSG_CHECKING(name of the STL string stream header file)
btng_stl_sstream_test_body='/* New syntax */ istringstream ist("a string");'
dnl We think that the sun 4.2 compiler does not support the syntax,
dnl but we're not absolutely sure.
test "$CXX_ID" = "sunpro" && echo "$CXX_VERSION" | grep '^0x420' > /dev/null && \
btng_stl_sstream_test_body='/* Old syntax */ char i[[10]]; istrstream ist(i);'
test "$CXX_ID" = "gnu" && echo "$CXX_VERSION" | grep '^2.95.2' > /dev/null && \
btng_stl_sstream_test_body='/* Old syntax */ char i[[10]]; istrstream ist(i);'
CASC_TREAT_VARIABLE_HEADER_FILENAME([stl-sstream],
  [sstream stringstream strstream sstream.h stringstream.h strstream.h sstream.hxx stringstream.hxx strstream.hxx] ,,
  [$btng_stl_sstream_test_body] )
# End macro $0
])      dnl end of CASC_STL_SSTREAM_HEADER_FILENAME definition.


AC_DEFUN([CASC_STL_MULTIMAP_HEADER_FILENAME],[
# Start macro $0
dnl AC_MSG_CHECKING(name of the STL multimap header file)
AC_REQUIRE([CASC_INFO_CXX_ID])
btng_stl_multimap_test_body='[multimap<int,int > v; pair<const int,int> thePair(0,1); v.insert(thePair);]'
test "$CXX_ID" = "sunpro" && echo "$CXX_VERSION" | grep '^0x420' > /dev/null && \
btng_stl_multimap_test_body='[multimap<int,int,less<int> > v; pair<const int,int> thePair(0,1); v.insert(thePair);]'
CASC_TREAT_VARIABLE_HEADER_FILENAME([stl-multimap],
    [multimap mmap multimap.h mmap.h multimap.hxx mmap.hxx map map.h map.hxx],,
    [$btng_stl_multimap_test_body])
# End macro $0
])      dnl end of CASC_STL_MULTIMAP_HEADER_FILENAME definition.


AC_DEFUN([CASC_STL_PAIR_HEADER_FILENAME],[
# Start macro $0
dnl AC_MSG_CHECKING(name of the STL pair header file)
CASC_TREAT_VARIABLE_HEADER_FILENAME([stl-pair], [pair pair.h pair.hxx],,
  [pair<int,int> s(0,1);])
# End macro $0
])      dnl end of CASC_STL_PAIR_HEADER_FILENAME definition.




dnl
dnl These are some stream-related headers with uncertain names.
dnl


AC_DEFUN([CASC_IOSTREAM_HEADER_FILENAME],[
# Start macro $0
dnl AC_MSG_CHECKING(name of the iostream header file)
CASC_TREAT_VARIABLE_HEADER_FILENAME([iostream],
  [iostream iostream.h iostream.hxx],,
  [ostream &co=cout; // test ostream declaration
   istream &ci=cin; // test istream declaration
   cout<<"test"<<endl; // test extraction operator
   ])
# End macro $0
])	dnl end of CASC_IOSTREAM_HEADER_FILENAME definition.


AC_DEFUN([CASC_FSTREAM_HEADER_FILENAME],[
# Start macro $0
dnl AC_MSG_CHECKING(name of the fstream header file)
CASC_TREAT_VARIABLE_HEADER_FILENAME([fstream],
  [fstream fstream.h fstream.hxx],,
  [fstream iost("theStream",ios::app);])
# End macro $0
])	dnl end of CASC_FSTREAM_HEADER_FILENAME definition.


AC_DEFUN([CASC_IOMANIP_HEADER_FILENAME],[
# Start macro $0
dnl AC_MSG_CHECKING(name of the iomanip header file)
AC_REQUIRE([CASC_IOSTREAM_HEADER_FILENAME])
CASC_TREAT_VARIABLE_HEADER_FILENAME([iomanip],
  [iomanip iomanip.h iomanip.hxx],[#include IOSTREAM_HEADER_FILE],
  [cout<<setw(13)<<endl;])
# End macro $0
])	dnl end of CASC_IOMANIP_HEADER_FILENAME definition.


