c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: FORTRAN routines for spatial coarsening of 1d patch data
c               on a regular Cartesian mesh.
c
include(pdat_m4concoarsenops1d.i)dnl
c
c***********************************************************************
c Constant coarsening for 1d node-centered double data
c***********************************************************************
c
      subroutine conavgnodedoub1d(
conavg_op_node_1d(`double precision')dnl
c
c***********************************************************************
c Constant coarsening for 1d node-centered float data 
c***********************************************************************
c
      subroutine conavgnodeflot1d(
conavg_op_node_1d(`real')dnl
c
c***********************************************************************
c Constant coarsening for 1d node-centered complex data
c***********************************************************************
c
      subroutine conavgnodecplx1d(
conavg_op_node_1d(`double complex')dnl
c
c***********************************************************************
c Constant coarsening for 1d node-centered integer data
c***********************************************************************
c
      subroutine conavgnodeintg1d(
conavg_op_node_1d(`integer')dnl
c

c***********************************************************************
c Constant coarsening for 1d outernode-centered double data
c***********************************************************************
c
      subroutine conavgouternodedoub1d(
conavg_op_outernode_1d(`double precision')dnl

c***********************************************************************
c Constant coarsening for 1d outernode-centered float data
c***********************************************************************
c
      subroutine conavgouternodeflot1d0(
conavg_op_outernode_1d(`real')dnl

c
c***********************************************************************
c Constant coarsening for 1d outernode-centered complex data
c***********************************************************************
c
      subroutine conavgouternodecplx1d0(
conavg_op_outernode_1d(`complex')dnl

c
c***********************************************************************
c Constant coarsening for 1d outernode-centered integer data
c***********************************************************************
c
      subroutine conavgouternodeint1d0(
conavg_op_outernode_1d(`integer')dnl



