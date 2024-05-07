c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
c 
define(NDIM,2)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
c Convert side data to cell data by averaging over cell sides
c
      subroutine side2cell(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   sdata0, sdata1, ngs,
     &   cdata, ngc)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer ngs, ngc
c
c variables in 2d cell indexed
      double precision cdata(CELL2d(ifirst,ilast,ngc))
      double precision sdata0(SIDE2d0(ifirst,ilast,ngs))
      double precision sdata1(SIDE2d1(ifirst,ilast,ngs))
c
c***********************************************************************
c
      integer ic0, ic1
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0

            cdata(ic0,ic1) = 0.25d0 * (
     &         sdata0(ic0,ic1) + sdata0(ic0+1,ic1)
     &       + sdata1(ic0,ic1) + sdata1(ic0,ic1+1) )

         end do
      end do
c
      return
      end
