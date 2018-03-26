dnl $Id$
dnl
dnl This file was writen to kludge in fortran link flags
dnl when AC_F77_LIBRARY_LDFLAGS had some problems.  It is
dnl obsolete.  You should use AC_F77_LIBRARY_LDFLAGS instead.
dnl
dnl
dnl This file contains functions to crudely determine
dnl loader flags for linking C/C++ and Fortran.

AC_DEFUN([CASC_LIB_C_FORTRAN],
# Start macro $1
[
dnl Use this macro to set fortran_LIBS to the library flgas
dnl needed to link in fortran from C codes.
dnl
dnl The macro AC_F77_LIBRARY_LDFLAGS, which was installed
dnl with autoconf, did not work for me.  BTNG.

AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])

if test -z "$fortran_LIBS" ; then
  # Guessing the correct flags is just based on past experience.
  # But it seems to work well in this case.  We rely heavily on
  # the name of compiler, so compilers with strange names will
  # not successfily be guessed.
  AC_REQUIRE([CASC_INFO_CC_CXX_ID])
  CASC_AC_LOG(Guessing fortran_LIBS in $1)
  case "$CXX_ID" in
    gnu)	# The GNU compiler.
      CASC_AC_LOG(About to check for RedHat gnu-2.96 compiler)
      CASC_AC_LOG_VAR(CXX_ID CXX_VERSION)
      case "$host_os" in
	linux*) # Check which version to determine fortran libraries
	[ FOUND_G2C=`echo $FLIBS | sed -e 's/.*g2c.*//'` ]
        if  test -z "$FOUND_G2C" ; then
	    fortran_LIBS="-lg2c"
	else
       	    fortran_LIBS=""
	fi

        if test "$CXX_VERSION" = "2.96" ; then
	  # Fix problem with RedHat distribution of experimental
	  # gcc version 2.96 not knowing where it put the g2c library.
          [ extra_path=`which $CXX | sed 's:bin/[^/]*$:lib:'` ]
      	  CASC_AC_LOG_VAR(extra_path)
	  test -d "$extra_path/gcc-lib/i386-redhat-linux/2.96" && extra_path="$extra_path/gcc-lib/i386-redhat-linux/2.96"
      	  CASC_AC_LOG_VAR(extra_path)
	  if test -d "$extra_path" ; then
	    CASC_AC_LOG(Adding $extra_path to fortran_LIBS)
	    fortran_LIBS="-L$extra_path $fortran_LIBS"
	  fi
        fi
	;;
	solaris*)	fortran_LIBS='-lg2c' ;;
      esac
    ;;
    kai)	# The KCC compiler runs on many platforms.
      # This compiler usually uses the native loader so use the
      # flags from the native compiler/loader.
      case "$host_os" in
	osf*)		fortran_LIBS='-lfor' ;;
	linux*)		fortran_LIBS='-lg2c' ;;
	sun*|solaris*)	fortran_LIBS='-L/opt/SUNWspro/lib -lF77 -lM77 -lV77 -lsunmath' ;;
	sgi)		fortran_LIBS='-lftn' ;;
      esac
    ;;
    sgi)	# SGI compiler.
        fortran_LIBS='-lftn'
    ;;
    sunpro)		# Sunpro compiler.
	fortran_LIBS='-lF77 -lM77 -lV77 -lsunmath'
    ;;
    dec)	# The DEC compiler.
	fortran_LIBS='-lfor'
    ;;
  esac
fi

]
# End macro $1
)dnl	End definition of CASC_LIB_C_FORTRAN



AC_DEFUN([CASC_SUPPORT_FORTRAN_FROM_C],
# Begin macro $1
[
dnl Use this macro if you have to link with Fortran codes.
dnl The macro AC_F77_LIBRARY_LDFLAGS, which was installed
dnl with autoconf, did not work for me.  BTNG.
dnl
dnl Arg1: empty if you want the default to be off.
dnl
ifelse($1,,
[
CASC_ARG_WITH_ENV_WRAPPER(fortran-libs, fortran_LIBS,
[  --with-fortran-libs[=LIBS]
		Enable support for linking with FORTRAN code
		and optionally specify the FORTRAN library flags.]
)
],[
CASC_ARG_WITH_ENV_WRAPPER(fortran-libs, fortran_LIBS,
[  --without-fortran-libs
		Do not try to specify the FORTRAN library linking flags.
  --with-fortran-libs[=LIBS]
		Specify the FORTRAN library linking flags.],
[[CASC_LIB_C_FORTRAN]]
)
]
)
dnl CASC_ARG_WITH_ENV_WRAPPER(fortran-libs, fortran_LIBS,
dnl ifelse($1,,
dnl [  --with-fortran-libs[=LIBS]
dnl 		Enable support for linking with FORTRAN code
dnl 		and optionally specify the FORTRAN library flags.],
dnl [  --without-fortran-libs
dnl 		Do not try to specify the FORTRAN library linking flags.]),
dnl [[CASC_LIB_C_FORTRAN]]
dnl )
CASC_AC_LOG_VAR(fortran_LIBS)
]
# End macro $1
)


