#!/bin/sh 
#########################################################################
##
## This file is part of the SAMRAI distribution.  For full copyright 
## information, see COPYRIGHT and COPYING.LESSER. 
##
## Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
## Description:   
##
#########################################################################
if [ $1 = "-h" ]
then
       echo "usage: $0 filename"
       echo
       echo "Where filename is the rootname of the .vis file collection."
       echo "     example  $0 euler.00000"
       echo "will collect euler.00000.vis.00000 euler.00000.vis.00001 into euler.vis"
       echo "for a 2 processor run."
       echo
	
       exit 0
fi

# First cat all the files together
rm -f $1.vis
FILES=`ls $1.vis.[0-9]*`
for i in $FILES
do
	cat $i >> $1.vis
done
echo

# Remove files
for i in $FILES
do
	rm -f $i
done



