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
