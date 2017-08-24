dnl Define a macro for supporting VALGRIND

AC_DEFUN([CASC_SUPPORT_VALGRIND],[

# Begin CASC_SUPPORT_VALGRIND
# Defines valgrind_EXE
AC_ARG_WITH(valgrind,
[  --with-valgrind[=PATH]  Use VALGRIND and optionally specify where VALGRIND is installed.],
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
