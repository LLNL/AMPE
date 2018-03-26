c
c  File:        $URL$
c  Package:     SAMRAI 
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Release:     
c  Revision:    
c  Modified:    
c  Description:    F77 routines for setting embedded boundary conditions.
c
define(SAMRAI_FORTDIR,../../pdat/fortran)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c
c***********************************************************************
c Convert node inout values to cell flags.
c***********************************************************************
c
      subroutine node2cellflag3d(
     &  ifirst0,ifirst1,ifirst2,
     &  ilast0,ilast1,ilast2,
     &  ngc0,ngc1,ngc2,
     &  cgc0,cgc1,cgc2,
     &  node_flag,
     &  cell_flag,
     &  cell_vol)
c***********************************************************************
      implicit none
include(ebparams.i)dnl
c***********************************************************************
c
      integer
     &  ifirst0,ifirst1,ifirst2,
     &  ilast0,ilast1,ilast2,
     &  ngc0,ngc1,ngc2,
     &  cgc0,cgc1,cgc2
      integer
     &  node_flag(NODE3dVECG(ifirst,ilast,ngc))
      integer
     &  cell_flag(CELL3dVECG(ifirst,ilast,cgc))
      double precision
     &  cell_vol(CELL3dVECG(ifirst,ilast,cgc))
      integer ic0,ic1,ic2
c
c  The cell flag is equal to the sum of the node flags around it.
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1  
            do ic0 = ifirst0, ilast0
               cell_flag(ic0,ic1,ic2) = 
     &            node_flag(ic0,ic1,ic2) + 
     &            node_flag(ic0+1,ic1,ic2) +
     &            node_flag(ic0,ic1+1,ic2) + 
     &            node_flag(ic0+1,ic1+1,ic2) +
     &            node_flag(ic0,ic1,ic2+1) + 
     &            node_flag(ic0+1,ic1,ic2+1) +
     &            node_flag(ic0,ic1+1,ic2+1) + 
     &            node_flag(ic0+1,ic1+1,ic2+1) 
            end do
         end do
      end do
c
c  If the cell flag is:
c      flag = 0          - FLOW cell
c      flag = 2^NDIM     - SOLID cell
c      0 < flag < 2^NDIM - CUT cell
c
c  Set volume of cut cell to -1.0 to trip an error if we don't reset it.
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1  
            do ic0 = ifirst0, ilast0
               if (cell_flag(ic0,ic1,ic2).eq.0) then
                  cell_flag(ic0,ic1,ic2) = FLOW_EB
                  cell_vol(ic0,ic1,ic2) = 1.0
               else if (cell_flag(ic0,ic1,ic2).eq.2**3) then
                  cell_flag(ic0,ic1,ic2) = SOLID_EB
                  cell_vol(ic0,ic1,ic2) = 0.0
               else
                  cell_flag(ic0,ic1,ic2) = CUT_EB
                  cell_vol(ic0,ic1,ic2) = -1.0
               endif
            end do
         end do
      end do

      return
      end
               
