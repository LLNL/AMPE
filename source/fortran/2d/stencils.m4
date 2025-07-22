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

      subroutine stencil5pts(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   diag, offdiagx, offdiagy,
     &   field, ngfield,
     &   rhs, ngrhs )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1

      double precision diag, offdiagx, offdiagy
      integer ngfield, ngrhs
c
c variables in 2d cell indexed
      double precision rhs(CELL2d(ifirst,ilast,ngrhs))
      double precision field(CELL2d(ifirst,ilast,ngfield))
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0
            rhs(ic0,ic1) = diag * field(ic0,ic1)
     &               + offdiagx * (field(ic0-1,ic1)+field(ic0+1,ic1))
     &               + offdiagy * (field(ic0,ic1-1)+field(ic0,ic1+1))

         enddo
      enddo

      return
      end

