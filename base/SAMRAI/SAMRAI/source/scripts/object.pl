#########################################################################
##
## This file is part of the SAMRAI distribution.  For full copyright 
## information, see COPYRIGHT and COPYING.LESSER. 
##
## Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
## Description:   perl script to create Makefile.objects.tmp files from .C files 
##
#########################################################################

#
# Usage: object.pl <.C files>
#

#
# Set up global directory information and output depend file
#

$TARGET = shift(@ARGV);
@FILES   = @ARGV;
$DEPEND  = "Makefile.objects.tmp";

$TABLEN  = 8;
$LINLEN  = 72;

#
# Write an output object file that contains the .o files from the .C input
#

open(OUTFILE, ">>$DEPEND") || die "Cannot open output file $DEPEND...";

$EMPTY="";

if ($TARGET eq "libdefault") {
    print OUTFILE "##\n";
    print OUTFILE "## File:\tMakefile.objects\n";
    print OUTFILE "## Package:\tSAMRAI\n";
    print OUTFILE "## Copyright:\t(c) 1997-2012 Lawrence Livermore National Security, LLC\n";
    print OUTFILE "## Release:\t\$N${EMPTY}ame:  \$\n";
    print OUTFILE "## Revision:\t\$R${EMPTY}evision: \$\n";
    print OUTFILE "## Modified:\t\$D${EMPTY}ate: \$\n";
    print OUTFILE "## Description:\tmakefile objects\n";
    print OUTFILE "##\n\n";
}

print OUTFILE "\n\n";

@FILES=sort(@FILES);

# Reorder list of files by nested template level.  This is done so the
#"real" templates are compiled before things like Pointers, Lists,
#etc.

for(@FILES) {

    # Find number of "-" characters, this is the number of nested
    #templates.

    $string=$_;
    $lookfor="-";
    $template_depth=0;
    $pos = $[;
    while (( $pos = index($string,$lookfor,$pos)) >= $[) {
	$pos++;
	$template_depth++;
    }

    # Create 2D array ordered by template depth

    push @{$SORTED_FILES[$template_depth]}, $string;
}

# Rebuild list of files with new order

@FILES=();
for $template_depth (0 .. $#SORTED_FILES) {
    @FILES = (@FILES, @{$SORTED_FILES[$template_depth]});
}

printObjects(@FILES);
close(OUTFILE);

#
# Print out object data in a pleasing manner
#

sub printObjects {
   my @OBJS = @_;

   print OUTFILE "$TARGET:\t";
   $CURLEN = $TABLEN;
   for (@OBJS) {
      $OBJ = $_;
      $OBJ =~ s/^(.*)\.C/$1.o/;
      $OBJ =~ s/X.o/\$\{NDIM\}\.o/;
      if (length($OBJ)+$CURLEN >= $LINLEN) {
         $NTAB = ($LINLEN-$CURLEN)/$TABLEN;
         for ($i = 0; $i < $NTAB; $i++) {
            print OUTFILE "\t";
         }
         print OUTFILE "\\\n";
         print OUTFILE "\t";
         $CURLEN = $TABLEN;
      }
      if ($CURLEN == $TABLEN) {
         print OUTFILE "$OBJ";
         $CURLEN += length($OBJ);
      } else {
         print OUTFILE " $OBJ";
         $CURLEN += length($OBJ)+1;
      }
   }
   print OUTFILE "\n";
}
