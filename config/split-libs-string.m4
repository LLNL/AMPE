AC_DEFUN([SAMRAI_SPLIT_LIBS_STRING],[
dnl
dnl Macro SAMRAI_SPLIT_LIBS_STRING
dnl Written by Brian Gunney.
dnl
dnl This macro takes an automake-style LIBS string (arg1) and
dnl splits it into the -L part (arg2, what SAMRAI usually calls
dnl LIB_PATH) and -l part (arg3, what SAMRAI calls LIB_NAME).
dnl The rest are also lumped into the LIB_NAME part, for lack
dnl of generality in the SAMRAI distinction.
dnl
dnl I think the SAMRAI format is limitting because not all parts
dnl if the LIBS strings can be recategorized into LIB_PATH and
dnl LIB_NAME parts.  But this macro allows us to use macros that
dnl conform to the standard automake format.  BTNG.
dnl
# Split $1 into the LIB_PATH part ($2) and the LIB_NAME part ($3)
if test -n "${$1}"; then
  for i in ${$1}; do
    case "$i" in
    -L*) $2="${$2} $i" ;;
    *) $3="${$3} $i" ;;
    esac
  done
fi
])
