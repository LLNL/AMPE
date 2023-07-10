c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
define(NDIM,3)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine grad_cell(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   diff_x, diff_y, diff_z, ngdiff,
     &   h,
     &   grad_x, grad_y, grad_z, nggrad
     &   )

      implicit none

      integer
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   ngdiff, nggrad

      double precision
     &   diff_x(SIDE3d0(ifirst,ilast,ngdiff)),
     &   diff_y(SIDE3d1(ifirst,ilast,ngdiff)),
     &   diff_z(SIDE3d2(ifirst,ilast,ngdiff)),
     &   h(3),
     &   grad_x(CELL3d(ifirst,ilast,nggrad)),
     &   grad_y(CELL3d(ifirst,ilast,nggrad)),
     &   grad_z(CELL3d(ifirst,ilast,nggrad))

      integer i, j, k
      double precision p5dxinv, p5dyinv, p5dzinv

      p5dxinv = 0.5d0 / h(1)
      p5dyinv = 0.5d0 / h(2)
      p5dzinv = 0.5d0 / h(3)

      do k = ifirst2, ilast2
         do j = ifirst1, ilast1
            do i = ifirst0, ilast0
               grad_x(i,j,k) =
     &            ( diff_x(i+1,j,k) + diff_x(i,j,k) ) * p5dxinv
            enddo
         enddo
      enddo
      
      do k = ifirst2, ilast2
         do j = ifirst1, ilast1
            do i = ifirst0, ilast0
               grad_y(i,j,k) =
     &            ( diff_y(i,j+1,k) + diff_y(i,j,k) ) * p5dyinv
            enddo
         enddo
      enddo

      do k = ifirst2, ilast2
         do j = ifirst1, ilast1
            do i = ifirst0, ilast0
               grad_z(i,j,k) =
     &            ( diff_z(i,j,k+1) + diff_z(i,j,k) ) * p5dzinv
            enddo
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine grad_side(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   diff_x, diff_y, diff_z, ngdiff,
     &   h,
     &   grad_x_xside, grad_y_xside, grad_z_xside,
     &   grad_x_yside, grad_y_yside, grad_z_yside,
     &   grad_x_zside, grad_y_zside, grad_z_zside,
     &   nggrad
     &   )

      implicit none

      integer
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   ngdiff, nggrad

      double precision
     &   diff_x(SIDE3d0(ifirst,ilast,ngdiff)),
     &   diff_y(SIDE3d1(ifirst,ilast,ngdiff)),
     &   diff_z(SIDE3d2(ifirst,ilast,ngdiff)),
     &   h(3),
     &   grad_x_xside(SIDE3d0(ifirst,ilast,nggrad)),
     &   grad_y_xside(SIDE3d0(ifirst,ilast,nggrad)),
     &   grad_z_xside(SIDE3d0(ifirst,ilast,nggrad)),
     &   grad_x_yside(SIDE3d1(ifirst,ilast,nggrad)),
     &   grad_y_yside(SIDE3d1(ifirst,ilast,nggrad)),
     &   grad_z_yside(SIDE3d1(ifirst,ilast,nggrad)),
     &   grad_x_zside(SIDE3d2(ifirst,ilast,nggrad)),
     &   grad_y_zside(SIDE3d2(ifirst,ilast,nggrad)),
     &   grad_z_zside(SIDE3d2(ifirst,ilast,nggrad))

c        local variables:
      integer i, j, k
      double precision dxinv, dyinv, dzinv
      double precision p25_dxinv, p25_dyinv, p25_dzinv

      dxinv = 1.d0 / h(1)
      dyinv = 1.d0 / h(2)
      dzinv = 1.d0 / h(3)
      p25_dxinv = 0.25d0 * dxinv
      p25_dyinv = 0.25d0 * dyinv
      p25_dzinv = 0.25d0 * dzinv

c        Compute gradients on the "X" side of the cell
      do k = ifirst2, ilast2
         do j = ifirst1, ilast1
            do i = ifirst0, ilast0+1

c                 X component of gradient at "X" side
               grad_x_xside(i,j,k) = dxinv * diff_x(i,j,k)

c                 Y component of gradient at "X" side
               grad_y_xside(i,j,k) = p25_dyinv * (
     &            diff_y(i-1,j+1,k) + diff_y(i-1,j,k) +
     &            diff_y(i,j+1,k) + diff_y(i,j,k)
     &            ) 

c                 Z component of gradient at "X" side
               grad_z_xside(i,j,k) = p25_dzinv * (
     &            diff_z(i-1,j,k+1) + diff_z(i-1,j,k) +
     &            diff_z(i,j,k+1) + diff_z(i,j,k)
     &            ) 

            enddo
         enddo
      enddo

c        Compute gradients on the "Y" side of the cell
      do k = ifirst2, ilast2
         do j = ifirst1, ilast1+1
            do i = ifirst0, ilast0

c                 X component of gradient at "Y" side
               grad_x_yside(i,j,k) = p25_dxinv * (
     &            diff_x(i+1,j-1,k) + diff_x(i,j-1,k) +
     &            diff_x(i+1,j,k) + diff_x(i,j,k)
     &            ) 

c                 Y component of gradient at "Y" side
               grad_y_yside(i,j,k) = dyinv * diff_y(i,j,k)

c                 Z component of gradient at "Y" side
               grad_z_yside(i,j,k) = p25_dzinv * (
     &            diff_z(i,j-1,k+1) + diff_z(i,j-1,k) +
     &            diff_z(i,j,k+1) + diff_z(i,j,k)
     &            ) 

            enddo
         enddo
      enddo

c        Compute gradients on the "Z" side of the cell
      do k = ifirst2, ilast2+1
         do j = ifirst1, ilast1
            do i = ifirst0, ilast0
               
c                 X component of gradient at "Z" side
               grad_x_zside(i,j,k) = p25_dxinv * (
     &            diff_x(i+1,j,k-1) + diff_x(i,j,k-1) +
     &            diff_x(i+1,j,k) + diff_x(i,j,k)
     &            ) 

c                 Y component of gradient at "Z" side
               grad_y_zside(i,j,k) = p25_dyinv * (
     &            diff_y(i,j+1,k-1) + diff_y(i,j,k-1) +
     &            diff_y(i,j+1,k) + diff_y(i,j,k)
     &            ) 

c                 Z component of gradient at "Z" side
               grad_z_zside(i,j,k) = dzinv * diff_z(i,j,k)
               
            enddo
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------


      subroutine velocity(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   threshold,
     &   grad_x, grad_y, grad_z, 
     &   phi_dot, 
     &   vel )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2

      double precision
     &   grad_x(lo0:hi0,lo1:hi1,lo2:hi2),
     &   grad_y(lo0:hi0,lo1:hi1,lo2:hi2),
     &   grad_z(lo0:hi0,lo1:hi1,lo2:hi2),
     &   phi_dot(lo0:hi0,lo1:hi1,lo2:hi2),
     &   vel(lo0:hi0,lo1:hi1,lo2:hi2)

      integer i, j, k
      double precision threshold, threshold2
      double precision invthreshold, grad2

      threshold2 = threshold*threshold
      invthreshold = 1./threshold
      
      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0

               grad2 =
     &            grad_x(i,j,k) * grad_x(i,j,k) +
     &            grad_y(i,j,k) * grad_y(i,j,k) +
     &            grad_z(i,j,k) * grad_z(i,j,k)
            
c               if( grad2>threshold2 )then
                vel(i,j,k) = phi_dot(i,j,k)/dsqrt(grad2+threshold2)
c                  vel(i,j,k) = phi_dot(i,j,k)
c               else
c                  vel(i,j,k) = phi_dot(i,j,k)*invthreshold
c                  vel(i,j,k) = phi_dot(i,j,k)
c               endif

            enddo
         enddo
      enddo

      return
      end
