#!/usr/bin/perl

#########################################################################
##
## This file is part of the SAMRAI distribution.  For full copyright 
## information, see COPYRIGHT and COPYING.LESSER. 
##
## Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
## Description:    Perl script for running the vizamrai_cat.sh script
##                 on multiple timesteps.  
##
#########################################################################

#
#  Determine path to vizamrai_cat.sh script
#
    # below line just sets the dirname where this perl file is located
    ($dirname = $0) =~ s:[^/]+$::; $dirname = '.' unless "$dirname";
    # the vizamrai_cat.sh file is presumed to be located in the same directory
    $vizamrai_cat = "$dirname/vizamrai_cat.sh";

    if (!((-e $vizamrai_cat) && (-r $vizamrai_cat))) {
       die "ERROR: cannot find $vizamrai_cat";
    }
#
#  Read file basename
#
    print "File basename (e.g. euler.<step>.vis.<proc>, enter 'euler'):";
    $basename = <STDIN>;
    chop($basename);
    print "Start Timestep (include zeros - e.g. 00000-00025, enter '00000'): ";
    $start = <STDIN>;
    chop($start);
    print "End Timestep: ";
    $end = <STDIN>;
    chop($end);

    $length = length($start);
#
#  Loop over files
#
   $inc = 1;
   for ($tstep = $start; $tstep <= $end; $tstep = $tstep + $inc) { 
      $newlength = length($tstep);
      while ($newlength < $length) {
         $tstep = "0$tstep";
      } continue {
         $newlength = length($tstep);
      }
      $firstfile = "$basename.$tstep.vis.00000";
      if (-e $firstfile) {
         system("/bin/sh $vizamrai_cat $basename.$tstep");
         print "Created $basename.$tstep ...\n";
      } else {
         print "$basename.$tstep.vis.<nproc> - files don't exist\n";
      }
   }
   print "Done. \n";
