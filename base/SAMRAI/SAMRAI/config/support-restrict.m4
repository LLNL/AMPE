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
