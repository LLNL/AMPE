#! /usr/bin/perl
#########################################################################
##
## This file is part of the SAMRAI distribution.  For full copyright 
## information, see COPYRIGHT and COPYING.LESSER. 
##
## Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
## Description:   shell script to update 1d, 2d, and 3d sed files 
##
#########################################################################

#
# Convert all files of the form "*X.*.sed" into their dimension-dependent
# forms in directories 1d, 2d, and 3d.
#

use File::Basename;
use File::Find;
use Cwd;

#
# Disallow running from certain directories.
#

my $pwd = cwd;
die basename($0) . " should not be run from your current directory"
    if $pwd =~ m:\b(source/test|source/scripts)(/|$):;


my $debug = 0;


#
# Convert all files of the form "*[123]d.m4" into their dimension-dependent
# forms in directories 1d, 2d, and 3d.
#

#
# Find the m4 files to convert.
#
@allfiles = ();
sub selectm4file {
    # This subroutine selects the dimensional m4 files in a find command,
    # excluding certain directories.
    # Do not run in applications directory.
    if ( $File::Find::name =~ m!/(test|scripts|CVS|\.hg|[123]d)$!o ) {
	$File::Find::prune = true;
    }
    elsif ( -f && m/.*[123]d\.m4$/o ) {
	push @allfiles, $File::Find::name;
	$allfiles[$#allfiles] =~ s|^\./||o;
    }
}
print "Scanning...\n" if ( $debug > 0 );
find( \&selectm4file, '.' );
print "Done.\n" if ( $debug > 0 );
# for (@allfiles) { print "$_\n"; }

#
# Recurse through all directories with sed files and update 1d, 2d, and 3d
#

my $savepwd = `pwd`; chop $savepwd;
for $dim (1,2,3) {
    @dimfiles = grep /${dim}d\.m4$/, @allfiles;
    print "dimfiles: @dimfiles\n" if ( $debug > 0 );

    for $mfile (@dimfiles) {
	print "Checking $mfile\n";
	$mdir = dirname $mfile;
	$mbase = basename $mfile;
	( $dfile = $mbase ) =~ s/d\.m4$/d.f/o;
	$tpath = "$mdir/$dfile.tmp"; # Temporary file.
	print "dfile: $dfile\n" if ( $debug > 0 );
	print "tpath: $tpath\n" if ( $debug > 0 );
	print "creating $tpath\n" if ( $debug > 0 );
	chdir $mdir || die "Cannot cd to $mdir";
	open MF, "m4 $mbase |" || die "Cannot launch 'm4 $mbase'";
	chdir($savepwd)
	    || die "Cannot cd back to $savepwd from $mdir due to $!";
	open TF, "> $tpath" || die "Cannot open file $tpath";
	while ( $str = <MF> ) {
	    # Remove CVS-substituted string.
	    $str =~ s/\$((Id|Name|Revision|Date)[^\$]*)\$//go;
	    print TF $str;
	}
	close MF || die "Cannot close process 'm4 $mfile'";
	close TF || die "Cannot close file $tpath";
	if ( &cmpfiles( "$mdir/$dfile", $tpath ) ) {
	    print "renaming $tpath to $dfile due to difference\n"
		if ( $debug > 0 );
	    rename( $tpath, "$mdir/$dfile" )
		|| die "Cannot rename $tpath to $dfile";
	} else {
	    print "Removing $tpath due to no difference\n"
		if ( $debug > 0 );
	    unlink $tpath || die "Cannot remove $tpath";
	}
    }
}



#
# Subroutine to check if two files are the same.
#
sub cmpfiles {

($ANAME,$BNAME) = @_;

return 1 unless ( -r $ANAME && -r $BNAME );
open(AFILE, "$ANAME") || die "Cannot open input file $ANAME...";
open(BFILE, "$BNAME") || die "Cannot open input file $BNAME...";

while (!eof(AFILE) && !eof(BFILE)) {
   $ALINE = <AFILE>;
   $BLINE = <BFILE>;
   if ($ALINE ne $BLINE) {
       close AFILE;
       close BFILE;
       return 1;
   }
}

if (eof(AFILE) && eof(BFILE)) {
    $rvalue = 0;
} else {
    $rvalue = 1;
}

close AFILE;
close BFILE;
return $rvalue;

}



