#! /usr/bin/perl
##
## File:        depend.pl
## Package:     SAMRAI scripts
## Copyright:   (c) 1997-2004 The Regents of the University of California
## Release:     $Name:  $
## Revision:    $Revision: 1.1.1.1 $
## Modified:    $Date: 2006/02/01 00:33:28 $
## Description: perl script to generate dependencies for SAMRAI files
##

##
## Usage: depend.pl <depend directory> <include directory> <.C files>
##

#
# Set up global directory information and output depend file
#

use File::Basename;
$scriptname = basename($0);

$debug = 0;

$this_command = "$scriptname @ARGV";
$src_dir  = shift; $src_dir = '.' if "$src_dir" eq '';
$inc_dir  = shift; $inc_dir = '-' if "$inc_dir" eq '';
@FILES   = @ARGV;
if ( ! @FILES ) {
    opendir( SRCDIR, $src_dir ) || die "Cannot open directory $src_dir";
    @FILES = grep( /.*\.[fCc]$/, sort(readdir SRCDIR) );
    closedir SRCDIR;
}
$DEPEND  = "Makefile.depend.tmp";

# By convention, directories should not end with '/'
# and relative paths should not start with './'.
# This is required to replace file searches with simpler string searches.
for ($src_dir, @FILES) { s|/+$||; s|^./+||o; }

if ( "$inc_dir" eq '-' ) {
    # Find the include directory by climbing back up the path,
    # past the source directory.
    $inc_dir = `pwd`; chop $inc_dir;
    $inc_dir =~ s|^.*/?source/?|source/|o;	# This gets past dir source.
    $inc_dir =~ s|[^/]+|..|go;		# Change each generation to '..'.
    $inc_dir = "$inc_dir/include";	# Append include.
}

@INCPATH = ($src_dir, $inc_dir);
print "$src_dir\n" if $debug;
print "$inc_dir\n" if $debug;
print "@FILES\n" if $debug;

$TABLEN  = 8;
$LINLEN  = 72;

#
# For each of the specified files, get dependency information and write to file
#

open(OUTFILE, ">$DEPEND") || die "Cannot open output file $DEPEND...";

$EMPTY="";

print OUTFILE <<__EOM__;
##
## File:\tMakefile.depend
## Package:\tSAMRAI
## Copyright:\t(c) 1997-2004 The Regents of the University of California
## Description:\tmakefile dependencies
##\n

## This file is automatically generated by $scriptname.


__EOM__

for $cfile (@FILES) {

    print "Checking file: $cfile\n" if $debug;

    undef %dset;	# The set of files that cfile depends on.
    @depfiles = ($cfile);	# List of files an object file compiled
    				# from $cfile would depend on.

    while (@depfiles) {
	$depfile = shift @depfiles;
	next if ( $depfile eq '' );
	if ( defined $dset{$depfile} ) {
	    print "$depfile is already in dependency set.\n" if $debug;
	    next;
	}
	# This file is not part of the dependency set.
	print "$depfile is being added the dependency set\n" if $debug;
	$dset{$depfile} = 1;	# Make current file a part of dependency set.
	# See what files $depfile depends on (and cache that info in
	# the variable filedeps).
	if ( ! defined $filedeps{$depfile} ) {
	    print "Finding filedeps for $depfile\n" if $debug;
	    $filedeps{$depfile} = [ &getMoreDeps( &getFullPath($depfile) ) ];
	}
	# See what other files $depfile depends on.
	for $maybedepfile ( @{$filedeps{$depfile}} ) {
	    $maybedepfile = &getFullPath($maybedepfile);
	    if ( $maybedepfile ne '' ) {
		print "Will also check file: $maybedepfile\n" if $debug;
		# Do not process the newly found include lines here.
		# Just add them to @depfiles to be processed by the
		# while loop.
		push( @depfiles, $maybedepfile )
		    if ! defined $dset{$maybedepfile}
	    }
	}
    }

    @deps = sort(keys %dset);
    for (@deps) { $_ = &fixName($_); }
    print "$cfile depends on @deps\n" if $debug;
    # Add SAMRAI_config.h because everything should depend on it,
    # even though it is ignored for the purpose of finding dependencies.
    # (It is ignored because it is generated at configure time.)
    unshift @deps, '$(SAMRAI)/include/SAMRAI/SAMRAI_config.h';
    printDependencies( $cfile, @deps );

}


sub getMoreDeps {
    my $file = @_[0];
    @deps = ();
    if (open(DEPFILE, $file)) {
	while ( <DEPFILE> ) {
	    if ( s/^\s*\#\s*include\s*\"([^\"]+)\"\s*/\1/o
		 && /[^\s]/o
		 && ! /\.f$/o
		 ) {
		push( @deps, $_ )
		}
	}
	close DEPFILE;
	if ( $debug ) {
	    print "$file is recursively dependent on\n";
	    for ( @deps ) { print "$_\n"; }
	}
    }
    else {
	print "Skipping recursion on file: $file\n" if $debug;
    }
    return @deps;
}


sub getFullPath {
   my $FILE    = shift(@_);

   for (@INCPATH) {
      if (-r "$_/$FILE") {
         return( $_ ne '.' ? "$_/$FILE" : $FILE );
      }
   }

   return("");
}


#
# Print out data dependencies in a pleasing manner
#

sub printDependencies {
   my $FILE = shift(@_);
   my @DEPS = @_;

   # print "FILE: $FILE\n" if $debug;
   my $LIBLINE = fixName($FILE);
   # print "LIBLINE: $LIBLINE\n" if $debug;
   $LIBLINE =~ s/^(.*)\.[Cfc]/$1.o:/o;
   # print "LIBLINE: $LIBLINE\n" if $debug;
   print OUTFILE "$LIBLINE";
   $NTAB = ($LINLEN-length($LIBLINE))/$TABLEN;
   for ($i = 0; $i < $NTAB; $i++) {
      print OUTFILE "\t";
   }
   print OUTFILE "\\\n";

   print OUTFILE "\t";
   $CURLEN = $TABLEN;
   for (@DEPS) {
      my $DEP = $_;
      if (length($DEP)+$CURLEN >= $LINLEN) {
         $NTAB = ($LINLEN-$CURLEN)/$TABLEN;
         for ($i = 0; $i < $NTAB; $i++) {
            print OUTFILE "\t";
         }
         print OUTFILE "\\\n";
         print OUTFILE "\t";
         $CURLEN = $TABLEN;
      }
      if ($CURLEN == $TABLEN) {
         print OUTFILE "$DEP";
         $CURLEN += length($DEP);
      } else {
         print OUTFILE " $DEP";
         $CURLEN += length($DEP)+1;
      }
   }
   print OUTFILE "\n";
}

#
# Convert the filename to a print name - remove include path or ./ prefix
#

sub fixName {
   $_ = shift(@_);
   # print "Fixing name $_\n" if $debug;
   if ( m|include/([^/]*)$|o ) {
      return("\$(INCLUDE_SAM)/$1");
      # return $_;
   } elsif ( m|/([^/]*)$|o ) {
      return($1);
   } else {
      return($_);
   }
}
