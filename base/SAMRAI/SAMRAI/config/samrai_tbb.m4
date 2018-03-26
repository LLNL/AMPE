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

