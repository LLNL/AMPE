c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for constant patchdata transfer routines.
c
define(coarsen_index,`dnl
         if ($1.lt.0) then
            $2=($1+1)/$3-1
         else
            $2=$1/$3
         endif
')dnl
define(coarsen_face_index,`dnl
         it=2*$1+$3
         if (it.le.0) then
            $2=it/(2*$3)-1
         else
            $2=(it-1)/(2*$3)
         endif
')dnl
