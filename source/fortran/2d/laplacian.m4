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

      subroutine laplacian(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   dx,
     &   coeff,
     &   field, ngfield,
     &   rhs, ngrhs )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1

      double precision dx(0:1)
      double precision coeff
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
      double precision dxinv2, dyinv2
      double precision diff_term_x, diff_term_y, diff_term

      dxinv2 = 1.d0/(dx(0)*dx(0))
      dyinv2 = 1.d0/(dx(1)*dx(1))
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0
            diff_term_x =
     &           (field(ic0-1,ic1)-2.d0*field(ic0,ic1)+field(ic0+1,ic1))
            diff_term_y =
     &           (field(ic0,ic1-1)-2.d0*field(ic0,ic1)+field(ic0,ic1+1))

            diff_term = diff_term_x*dxinv2 + diff_term_y*dyinv2

            rhs(ic0,ic1) = coeff * diff_term

         enddo
      enddo

      return
      end

