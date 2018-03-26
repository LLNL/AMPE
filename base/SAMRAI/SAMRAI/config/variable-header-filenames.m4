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

