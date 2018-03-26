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





