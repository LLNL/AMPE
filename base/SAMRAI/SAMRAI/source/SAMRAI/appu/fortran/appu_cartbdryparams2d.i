c
c  File:        $URL$
c  Package:     SAMRAI application utilities
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for 2d boundary constant common blocks
c
      common/cartbdrylocparams2d/
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,
     &  X0Y0,X1Y0,X0Y1,X1Y1
      integer
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,
     &  X0Y0,X1Y0,X0Y1,X1Y1
c
c
      common/cartbdrycondparams2d/
     &  FLOW,XFLOW,YFLOW,
     &  REFLECT,XREFLECT,YREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN
      integer
     &  FLOW,XFLOW,YFLOW,
     &  REFLECT,XREFLECT,YREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN
