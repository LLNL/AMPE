dnl Define a macro for supporting BOOST

AC_DEFUN([CASC_SUPPORT_BOOST],[

# Begin CASC_SUPPORT_BOOST
# Defines boost_PREFIX boost_INCLUDES and boost_LIBS.
AC_ARG_WITH(boost,
[  --with-boost[=PATH]  Use BOOST and specify where BOOST is installed.],
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

