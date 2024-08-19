c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
c 
define(NDIM,3)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c
c Convert side data to cell data by averaging over cell sides
c
      subroutine side2cell(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   sdata0, sdata1, sdata2, ngs,
     &   cdata, ngc)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngs, ngc
c
c variables in 3d cell indexed
      double precision cdata(CELL3d(ifirst,ilast,ngc))
      double precision sdata0(SIDE3d0(ifirst,ilast,ngs))
      double precision sdata1(SIDE3d1(ifirst,ilast,ngs))
      double precision sdata2(SIDE3d2(ifirst,ilast,ngs))

c
c***********************************************************************
c
      integer ic0, ic1, ic2
      double precision onesixth
c
      onesixth = 1.d0/6.d0
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               cdata(ic0,ic1,ic2) = onesixth * (
     &            sdata0(ic0,ic1,ic2) + sdata0(ic0+1,ic1,ic2)
     &          + sdata1(ic0,ic1,ic2) + sdata1(ic0,ic1+1,ic2)
     &          + sdata2(ic0,ic1,ic2) + sdata2(ic0,ic1,ic2+1) )

            end do
         end do
      end do
c
      return
      end
