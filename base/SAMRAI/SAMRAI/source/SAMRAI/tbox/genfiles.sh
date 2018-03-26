#########################################################################
##
## This file is part of the SAMRAI distribution.  For full copyright 
## information, see COPYRIGHT and COPYING.LESSER. 
##
## Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
## Description:   simple shell script to generate flex and bison files 
##
#########################################################################

dir_name=`echo ${0} | sed -e 's:^\([^/]*\)$:./\1:' -e 's:/[^/]*$::'`;
cd $dir_name

#
# Use yacc since ASCI red does not support alloca() function used by bison
#

bison -d -p yy Grammar.y
perl grammer_fixup.pl Grammar.tab.c > Grammar.C
perl grammer_fixup.pl Grammar.tab.h > Grammar.h
rm Grammar.tab.c
rm Grammar.tab.h

#
# Scanner requires flex due to input reading method with MPI
#

flex -Pyy -otemp.$$ Scanner.l
perl scanner_fixup.pl temp.$$ > Scanner.C 
rm temp.$$

# Add some pragma's to ignore warnings from compilers for
# machine generated code.
cat >> temp.$$ <<-EOF 
#ifdef __GNUC__
#if __GNUC__ > 4 || \
              (__GNUC__ == 4 && (__GNUC_MINOR__ > 2 || \
                                 (__GNUC_MINOR__ == 2 && \
                                  __GNUC_PATCHLEVEL__ > 0)))
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wconversion"
#endif
#endif

#ifdef __INTEL_COMPILER
// Ignore Intel warnings about unreachable statements
#pragma warning (disable:177)
// Ignore Intel warnings about external declarations
#pragma warning (disable:1419)
// Ignore Intel warnings about type conversions
#pragma warning (disable:810)
// Ignore Intel remarks about non-pointer conversions
#pragma warning (disable:2259)
// Ignore Intel remarks about zero used for undefined preprocessor syms
#pragma warning (disable:193)
// Ignore Intel remarks about unreachable code
#pragma warning (disable:111)
#endif

#ifdef __xlC__
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif
EOF
cat temp.$$ Scanner.C > temp2.$$
mv temp2.$$ Scanner.C
cat temp.$$ Grammar.C > temp2.$$
mv temp2.$$ Grammar.C
rm temp.$$

