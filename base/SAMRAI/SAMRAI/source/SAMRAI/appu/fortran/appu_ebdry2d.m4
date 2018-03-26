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
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
c***********************************************************************
c Set common definitions - integers passed from C++ code that define
c constants.
c***********************************************************************
c
      subroutine setebparams(
     &  SOLIDin,CUTin,BORDERin,FLOWin,
     &  OUTSIDEin,INSIDEin,BOUNDARYin,ONBOUNDARYin)
c***********************************************************************
      implicit none
      integer
     &  SOLIDin,CUTin,BORDERin,FLOWin,
     &  OUTSIDEin,INSIDEin,BOUNDARYin,ONBOUNDARYin
include(ebparams.i)dnl
c***********************************************************************
c
      SOLID_EB      = SOLIDin
      CUT_EB        = CUTin
      BORDER_EB     = BORDERin
      FLOW_EB       = FLOWin
      OUTSIDE_EB    = OUTSIDEin
      INSIDE_EB     = INSIDEin
      BOUNDARY_EB   = BOUNDARYin
      ONBOUNDARY_EB = ONBOUNDARYin

      return
      end
     
c
c***********************************************************************
c Convert node inout values to cell flags.
c***********************************************************************
c
      subroutine node2cellflag2d(
     &  ifirst0,ifirst1,
     &  ilast0,ilast1,
     &  ngc0,ngc1,
     &  cgc0,cgc1,
     &  node_flag,
     &  cell_flag,
     &  cell_vol)
c***********************************************************************
      implicit none
include(ebparams.i)dnl
c***********************************************************************
c
      integer
     &  ifirst0,ifirst1,
     &  ilast0,ilast1,
     &  ngc0,ngc1,
     &  cgc0,cgc1
      integer
     &  node_flag(NODE2dVECG(ifirst,ilast,ngc))
      integer
     &  cell_flag(CELL2dVECG(ifirst,ilast,cgc))
      double precision
     &  cell_vol(CELL2dVECG(ifirst,ilast,cgc))
      integer ic0,ic1
c
c  The cell flag is equal to the sum of the node flags around it.
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0
            cell_flag(ic0,ic1) = 
     &         node_flag(ic0,ic1) + node_flag(ic0+1,ic1) +
     &         node_flag(ic0,ic1+1) + node_flag(ic0+1,ic1+1)
         end do
      end do
c
c  If the cell flag is:
c      flag = 0             - FLOW cell
c      flag = 2^NDIM        - SOLID cell
c      flag = anything else - CUT cell
c
c  Set volume of cut cell to -1.0 to trip an error if we don't reset it.
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0
            if (cell_flag(ic0,ic1).eq.0) then
               cell_flag(ic0,ic1) = FLOW_EB
               cell_vol(ic0,ic1) = 1.0
            else if (cell_flag(ic0,ic1).eq.2**2) then
               cell_flag(ic0,ic1) = SOLID_EB
               cell_vol(ic0,ic1) = 0.0
            else
               cell_flag(ic0,ic1) = CUT_EB
               cell_vol(ic0,ic1) = -1.0
            endif
         end do
      end do

      return
      end
               

c
c***********************************************************************
c Set boundary condtions for embedded boundary
c***********************************************************************
c
      subroutine setebnode2d(
     &  ifirst0,ilast0,
     &  ifirst1,ilast1,
     &  ibeg0,iend0,
     &  ibeg1,iend1,
     &  ngc0, ngc1,
     &  bdry_loc,
     &  bdry_cond,
     &  cell_flag,
     &  cell_vol)
c***********************************************************************
      implicit none
include(appu_cartbdryparams2d.i)dnl
c***********************************************************************
c
      integer
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  ibeg0,iend0,ibeg1,iend1,
     &  ngc0,ngc1,
     &  bdry_loc, bdry_cond
      integer
     &  cell_flag(CELL2dVECG(ifirst,ilast,ngc))
      double precision
     &  cell_vol(CELL2dVECG(ifirst,ilast,ngc))
      integer ic0,ic1
c
c***********************************************************************
c
c
c Set the cell volume and tags appropriately
c
c    DIRICHLET,NEUMANN, - cell_tag(bdry) = cell_tag(interior)
c    REFLECT              cell_vol(bdry) = cell_vol(interior)
c


      if ((bdry_cond.eq.XDIRICHLET).or.
     &    (bdry_cond.eq.XREFLECT).or.
     &    (bdry_cond.eq.XNEUMANN)) then

         if (bdry_loc.eq.X0Y0) then        
            do ic1=ifirst1-ngc1,ifirst1-1     
               do ic0=ifirst0-ngc0,ifirst0-1  
                  cell_flag(ic0,ic1) = cell_flag(ifirst0,ic1)
                  cell_vol(ic0,ic1)  = cell_vol(ifirst0,ic1)
               end do
            end do
         else if (bdry_loc.eq.X0Y1) then
           do ic1=ilast1+1,ilast1+ngc1     
               do ic0=ifirst0-ngc0,ifirst0-1  
                  cell_flag(ic0,ic1) = cell_flag(ifirst0,ic1)
                  cell_vol(ic0,ic1)  = cell_vol(ifirst0,ic1)
               end do
            end do
         else if (bdry_loc.eq.X1Y0) then
           do ic1=ifirst1-ngc1,ifirst1-1     
               do ic0=ilast0+1,ilast0+ngc0  
                  cell_flag(ic0,ic1) = cell_flag(ilast0,ic1)
                  cell_vol(ic0,ic1)  = cell_vol(ilast0,ic1)
               end do
            end do
         else if (bdry_loc.eq.X1Y1) then
           do ic1=ilast1+1,ilast1+ngc1     
               do ic0=ilast0+1,ilast0+ngc0  
                  cell_flag(ic0,ic1) = cell_flag(ilast0,ic1)
                  cell_vol(ic0,ic1)  = cell_vol(ilast0,ic1)
               end do
            end do
         endif
   
      endif

      if ((bdry_cond.eq.YDIRICHLET).or.
     &    (bdry_cond.eq.YREFLECT).or.
     &    (bdry_cond.eq.YNEUMANN)) then

         if (bdry_loc.eq.X0Y0) then        
            do ic1=ifirst1-ngc1,ifirst1-1     
               do ic0=ifirst0-ngc0,ifirst0-1  
                  cell_flag(ic0,ic1) = cell_flag(ic0,ifirst1)
                  cell_vol(ic0,ic1)  = cell_vol(ic0,ifirst1)
               end do
            end do
         else if (bdry_loc.eq.X0Y1) then
           do ic1=ilast1+1,ilast1+ngc1     
               do ic0=ifirst0-ngc0,ifirst0-1  
                  cell_flag(ic0,ic1) = cell_flag(ic0,ilast1)
                  cell_vol(ic0,ic1)  = cell_vol(ic0,ilast1)
               end do
            end do
         else if (bdry_loc.eq.X1Y0) then
           do ic1=ifirst1-ngc1,ifirst1-1     
               do ic0=ilast0+1,ilast0+ngc0  
                  cell_flag(ic0,ic1) = cell_flag(ic0,ifirst1)
                  cell_vol(ic0,ic1)  = cell_vol(ic0,ifirst1)
               end do
            end do
         else if (bdry_loc.eq.X1Y1) then
           do ic1=ilast1+1,ilast1+ngc1     
               do ic0=ilast0+1,ilast0+ngc0  
                  cell_flag(ic0,ic1) = cell_flag(ic0,ilast1)
                  cell_vol(ic0,ic1)  = cell_vol(ic0,ilast1)
               end do
            end do
         endif
   
      endif

c
c    FLOW               - cell_tag(bdry) = -2
c                         cell_vol(bdry) = 0. 
c
      if ((bdry_cond.eq.XFLOW).or.
     &    (bdry_cond.eq.YFLOW)) then

         if (bdry_loc.eq.X0Y0) then        
            do ic1=ifirst1-ngc1,ifirst1-1     
               do ic0=ifirst0-ngc0,ifirst0-1  
                  cell_flag(ic0,ic1) = -2
                  cell_vol(ic0,ic1)  = 0.
               end do
            end do
         else if (bdry_loc.eq.X0Y1) then
            do ic1=ilast1+1,ilast1+ngc1     
               do ic0=ifirst0-ngc0,ifirst0-1       
                  cell_flag(ic0,ic1) = -2
                  cell_vol(ic0,ic1)  = 0.
               end do
            end do
         else if (bdry_loc.eq.X1Y0) then
            do ic1=ifirst1-ngc1,ifirst1-1     
               do ic0=ilast0+1,ilast0+ngc0
                  cell_flag(ic0,ic1) = -2
                  cell_vol(ic0,ic1)  = 0.
               end do
            end do
         else if (bdry_loc.eq.X1Y1) then
            do ic1=ilast1+1,ilast1+ngc1     
               do ic0=ilast0+1,ilast0+ngc0  
                 cell_flag(ic0,ic1) = -2
                  cell_vol(ic0,ic1)  = 0.
               end do
            end do
         endif
      endif
  

      return
      end
c
