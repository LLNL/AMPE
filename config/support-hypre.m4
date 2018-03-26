dnl Define a macro for supporting HYPRE

AC_DEFUN([BTNG_SUPPORT_HYPRE],[

# Begin BTNG_SUPPORT_HYPRE
# Defines hypre_PREFIX hypre_INCLUDES and hypre_LIBS.
AC_ARG_WITH(hypre,
[ --with-hypre[=PATH]  Use HYPRE and specify where HYPRE is installed.],
, with_hypre=no)

case "$with_hypre" in
  no)
    : Do nothing
  ;;
  yes)
    # HYPRE install path was not specified.
    # Look in a couple of standard locations to probe if 
    # HYPRE header files are there.
    AC_MSG_CHECKING([for HYPRE installation])
    for dir in /usr /usr/local; do
      if test -f ${dir}/include/HYPRE.h; then
        hypre_PREFIX=${dir}
        break
      fi
    done
    AC_MSG_RESULT([$hypre_PREFIX])
  ;;
  *)
    # HYPRE install path was specified.
    AC_MSG_CHECKING([for HYPRE installation])

    if test -f ${with_hypre}/include/HYPRE.h; then
        hypre_PREFIX=$with_hypre
        hypre_INCLUDES="-I${hypre_PREFIX}/include"
        dnl Do not have to set hypre_LIBS because SAMRAI does not need it (yet).
        dnl hypre_LIBS="-L${hypre_PREFIX}/lib -lhypre"
        AC_MSG_RESULT([$hypre_PREFIX])
    else
        AC_MSG_RESULT([$hypre_PREFIX])
        AC_MSG_ERROR([HYPRE not found in $with_hypre])
    fi
  ;;
esac


# END BTNG_SUPPORT_HYPRE

])dnl End definition of BTNG_SUPPORT_HYPRE

