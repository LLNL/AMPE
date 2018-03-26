dnl Define a macro for supporting NETCDF

AC_DEFUN([BTNG_SUPPORT_NETCDF],[

# Begin BTNG_SUPPORT_NETCDF
# Defines netcdf_PREFIX netcdf_INCLUDES and netcdf_LIBS.
AC_ARG_WITH(netcdf,
[ --with-netcdf[=PATH]  Use NETCDF and specify where NETCDF is installed.],
, with_netcdf=no)

case "$with_netcdf" in
  no)
    : Do nothing
  ;;
  yes)
    # NETCDF install path was not specified.
    # Look in a couple of standard locations to probe if 
    # NETCDF header files are there.
    AC_MSG_CHECKING([for NETCDF installation])
    for dir in /usr /usr/local; do
      if test -f ${dir}/include/netcdf.h; then
        netcdf_PREFIX=${dir}
        break
      fi
    done
    AC_MSG_RESULT([$netcdf_PREFIX])
  ;;
  *)
    # NETCDF install path was specified.
    AC_MSG_CHECKING([for NETCDF installation])

    if test -f ${with_netcdf}/include/netcdf.h; then
        netcdf_PREFIX=$with_netcdf
        netcdf_INCLUDES="-I${netcdf_PREFIX}/include"
        dnl Do not have to set netcdf_LIBS because SAMRAI does not need it (yet).
        dnl netcdf_LIBS="-L${netcdf_PREFIX}/lib -lnetcdf"
        AC_MSG_RESULT([$netcdf_PREFIX])
    else
        AC_MSG_RESULT([$netcdf_PREFIX])
        AC_MSG_ERROR([NETCDF not found in $with_netcdf])
    fi
  ;;
esac


# END BTNG_SUPPORT_NETCDF

])dnl End definition of BTNG_SUPPORT_NETCDF

