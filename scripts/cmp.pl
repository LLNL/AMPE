#!/usr/local/bin/perl
##
## File:        cmp.pl
## Package:     SAMRAI scripts
## Copyright:   (c) 1997-2004 The Regents of the University of California
## Release:     $Name:  $
## Revision:    $Revision: 1.1.1.1 $
## Modified:    $Date: 2006/02/01 00:33:28 $
## Description: perl script to compare two files but ignore CVS comments
##

##
## Usage: cmp.pl <file1> <file2>
##

$ANAME = shift(@ARGV);
$BNAME = shift(@ARGV);

open(AFILE, "$ANAME") || die "Cannot open input file $ANAME...";
open(BFILE, "$BNAME") || die "Cannot open input file $BNAME...";

while (!eof(AFILE) && !eof(BFILE)) {
   $ALINE = <AFILE>;
   $BLINE = <BFILE>;
   $_ = $ALINE;

   if (!/^(\/\/|c|C|#|##| \*)[ ]*(Release:[\t ]*\$Name|Revision:[\t ]*\$Revision|Modified:[\t ]*\$Date):[^\$]*\$/o) {
      if ($ALINE ne $BLINE) {
         exit 1;
      }
   }
}

exit 0 if (eof(AFILE) && eof(BFILE));
exit 1;
