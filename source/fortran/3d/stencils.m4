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

      subroutine stencil7pts(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   diag, offdiagx, offdiagy, offdiagz,
     &   field, ngfield,
     &   rhs, ngrhs )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2

      double precision diag, offdiagx, offdiagy, offdiagz
      integer ngfield, ngrhs
c
c variables in 3d cell indexed
      double precision rhs(CELL3d(ifirst,ilast,ngrhs))
      double precision field(CELL3d(ifirst,ilast,ngfield))
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1, ic2
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               rhs(ic0,ic1,ic2) = diag * field(ic0,ic1,ic2)
     &            + offdiagx*(field(ic0-1,ic1,ic2)+field(ic0+1,ic1,ic2))
     &            + offdiagy*(field(ic0,ic1-1,ic2)+field(ic0,ic1+1,ic2))
     &            + offdiagz*(field(ic0,ic1,ic2-1)+field(ic0,ic1,ic2+1))
            enddo
         enddo
      enddo

      return
      end
