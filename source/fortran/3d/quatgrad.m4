c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
define(NDIM,3)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine quatgrad_cell(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth, h,
     &   diff_x, diff_y, diff_z, ngdiff,
     &   grad_x, grad_y, grad_z, nggrad
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth, ngdiff, nggrad
 
      double precision diff_x(SIDE3d0(lo,hi,ngdiff),depth)
      double precision diff_y(SIDE3d1(lo,hi,ngdiff),depth)
      double precision diff_z(SIDE3d2(lo,hi,ngdiff),depth)
      double precision grad_x(CELL3d(lo,hi,nggrad),depth)
      double precision grad_y(CELL3d(lo,hi,nggrad),depth)
      double precision grad_z(CELL3d(lo,hi,nggrad),depth)
      double precision h(3)

c        local variables:
      integer i, j, k, m
      double precision p5dxinv, p5dyinv, p5dzinv

      p5dxinv = 0.5d0 / h(1)
      p5dyinv = 0.5d0 / h(2)
      p5dzinv = 0.5d0 / h(3)

c        Loop over the quaternion components
      do m = 1, depth

         do k = lo2, hi2
            do j = lo1, hi1
               do i = lo0, hi0
                  grad_x(i,j,k,m) =
     &               ( diff_x(i+1,j,k,m) + diff_x(i,j,k,m) ) * p5dxinv
               enddo
            enddo
         enddo

         do k = lo2, hi2
            do j = lo1, hi1
               do i = lo0, hi0
                  grad_y(i,j,k,m) =
     &               ( diff_y(i,j+1,k,m) + diff_y(i,j,k,m) ) * p5dyinv
               enddo
            enddo
         enddo

         do k = lo2, hi2
            do j = lo1, hi1
               do i = lo0, hi0
                  grad_z(i,j,k,m) =
     &               ( diff_z(i,j,k+1,m) + diff_z(i,j,k,m) ) * p5dzinv
               enddo
            enddo
         enddo

      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine quatgrad_cell_symm(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth, h,
     &   diff_x, diff_y, diff_z, ngdiff,
     &   grad_x, grad_y, grad_z, nggrad,
     &   iqrot_x, iqrot_y, iqrot_z, ngiq
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth, ngdiff, nggrad, ngiq
      double precision h(3)
 
      double precision diff_x(SIDE3d0(lo,hi,ngdiff),depth)
      double precision diff_y(SIDE3d1(lo,hi,ngdiff),depth)
      double precision diff_z(SIDE3d2(lo,hi,ngdiff),depth)
      double precision grad_x(CELL3d(lo,hi,nggrad),depth)
      double precision grad_y(CELL3d(lo,hi,nggrad),depth)
      double precision grad_z(CELL3d(lo,hi,nggrad),depth)

      integer iqrot_x(SIDE3d0(lo,hi,ngiq))
      integer iqrot_y(SIDE3d1(lo,hi,ngiq))
      integer iqrot_z(SIDE3d2(lo,hi,ngiq))

c        local variables:
      integer i, j, k, m, iq
      double precision p5dxinv, p5dyinv, p5dzinv
      double precision dtmp(depth)
      double precision dprime(depth)

      p5dxinv = 0.5d0 / h(1)
      p5dyinv = 0.5d0 / h(2)
      p5dzinv = 0.5d0 / h(3)

      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0

               if ( depth > 1 ) then
c                    Diffs on right side need to have the conjugate
c                    of the symmetry rotation applied before using, but
c                    only for the vector/quaternion orientation.

                  do m = 1, depth
                     dtmp(m) = diff_x(i+1,j,k,m)
                  enddo

                  iq = -1 * iqrot_x(i+1,j,k)
                  call quatsymmrotate( dtmp, iq, dprime, depth )

               else

                  dprime(1) = diff_x(i+1,j,k,1)

               endif

               do m = 1, depth
                  grad_x(i,j,k,m) =
     &               ( dprime(m) + diff_x(i,j,k,m) ) * p5dxinv
               enddo

            enddo
         enddo
      enddo

      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0

               if ( depth > 1 ) then

c                    Diffs on top side need to have the conjugate
c                    of the symmetry rotation applied before using, but
c                    only for the vector/quaternion orientation.

                  do m = 1, depth
                     dtmp(m) = diff_y(i,j+1,k,m)
                  enddo

                  iq = -1 * iqrot_y(i,j+1,k)
                  call quatsymmrotate( dtmp, iq, dprime, depth )

               else

                  dprime(1) = diff_y(i,j+1,k,1)

               endif

               do m = 1, depth
                  grad_y(i,j,k,m) =
     &               ( dprime(m) + diff_y(i,j,k,m) ) * p5dyinv
               enddo

            enddo
         enddo
      enddo

      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0

               if ( depth > 1 ) then

c                    Diffs on back side need to have the conjugate
c                    of the symmetry rotation applied before using, but
c                    only for the vector/quaternion orientation.

                  do m = 1, depth
                     dtmp(m) = diff_z(i,j,k+1,m)
                  enddo

                  iq = -1 * iqrot_z(i,j,k+1)
                  call quatsymmrotate( dtmp, iq, dprime, depth )

               else

                  dprime(1) = diff_z(i,j,k+1,1)

               endif

               do m = 1, depth
                  grad_z(i,j,k,m) =
     &               ( dprime(m) + diff_z(i,j,k,m) ) * p5dzinv
               enddo

            enddo
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine quatgrad_side(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth, h,
     &   diff_x, diff_y, diff_z, ngdiff,
     &   grad_x_xside, grad_y_xside, grad_z_xside,
     &   grad_x_yside, grad_y_yside, grad_z_yside,
     &   grad_x_zside, grad_y_zside, grad_z_zside,
     &   nggrad
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth, ngdiff, nggrad

      double precision diff_x(SIDE3d0(lo,hi,ngdiff),depth)
      double precision diff_y(SIDE3d1(lo,hi,ngdiff),depth)
      double precision diff_z(SIDE3d2(lo,hi,ngdiff),depth)
      double precision grad_x_xside(SIDE3d0(lo,hi,nggrad),depth)
      double precision grad_x_yside(SIDE3d1(lo,hi,nggrad),depth)
      double precision grad_x_zside(SIDE3d2(lo,hi,nggrad),depth)
      double precision grad_y_xside(SIDE3d0(lo,hi,nggrad),depth)
      double precision grad_y_yside(SIDE3d1(lo,hi,nggrad),depth)
      double precision grad_y_zside(SIDE3d2(lo,hi,nggrad),depth)
      double precision grad_z_xside(SIDE3d0(lo,hi,nggrad),depth)
      double precision grad_z_yside(SIDE3d1(lo,hi,nggrad),depth)
      double precision grad_z_zside(SIDE3d2(lo,hi,nggrad),depth)
      double precision h(3)

c        local variables:
      integer i, j, k, m
      double precision dxinv, dyinv, dzinv
      double precision p25_dxinv, p25_dyinv, p25_dzinv

      dxinv = 1.d0 / h(1)
      dyinv = 1.d0 / h(2)
      dzinv = 1.d0 / h(3)
      p25_dxinv = 0.25d0 * dxinv
      p25_dyinv = 0.25d0 * dyinv
      p25_dzinv = 0.25d0 * dzinv

c        Loop over the quaternion components
      do m = 1, depth

c           Compute gradients on the "X" side of the cell
         do k = lo2, hi2
            do j = lo1, hi1
               do i = lo0, hi0+1

c                    X component of gradient at "X" side
                  grad_x_xside(i,j,k,m) = dxinv * diff_x(i,j,k,m)

c                    Y component of gradient at "X" side
                  grad_y_xside(i,j,k,m) = p25_dyinv * (
     &               diff_y(i-1,j+1,k,m) + diff_y(i-1,j,k,m) +
     &               diff_y(i,j+1,k,m) + diff_y(i,j,k,m)
     &               ) 

c                    Z component of gradient at "X" side
                  grad_z_xside(i,j,k,m) = p25_dzinv * (
     &               diff_z(i-1,j,k+1,m) + diff_z(i-1,j,k,m) +
     &               diff_z(i,j,k+1,m) + diff_z(i,j,k,m)
     &               ) 

               enddo
            enddo
         enddo

c           Compute gradients on the "Y" side of the cell
         do k = lo2, hi2
            do j = lo1, hi1+1
               do i = lo0, hi0

c                    X component of gradient at "Y" side
                  grad_x_yside(i,j,k,m) = p25_dxinv * (
     &               diff_x(i+1,j-1,k,m) + diff_x(i,j-1,k,m) +
     &               diff_x(i+1,j,k,m) + diff_x(i,j,k,m)
     &               ) 

c                    Y component of gradient at "Y" side
                  grad_y_yside(i,j,k,m) = dyinv * diff_y(i,j,k,m)

c                    Z component of gradient at "Y" side
                  grad_z_yside(i,j,k,m) = p25_dzinv * (
     &               diff_z(i,j-1,k+1,m) + diff_z(i,j-1,k,m) +
     &               diff_z(i,j,k+1,m) + diff_z(i,j,k,m)
     &               ) 

               enddo
            enddo
         enddo

c           Compute gradients on the "Z" side of the cell
         do k = lo2, hi2+1
            do j = lo1, hi1
               do i = lo0, hi0

c                    X component of gradient at "Z" side
                  grad_x_zside(i,j,k,m) = p25_dxinv * (
     &               diff_x(i+1,j,k-1,m) + diff_x(i,j,k-1,m) +
     &               diff_x(i+1,j,k,m) + diff_x(i,j,k,m)
     &               ) 

c                    Y component of gradient at "Z" side
                  grad_y_zside(i,j,k,m) = p25_dyinv * (
     &               diff_y(i,j+1,k-1,m) + diff_y(i,j,k-1,m) +
     &               diff_y(i,j+1,k,m) + diff_y(i,j,k,m)
     &               ) 

c                    Z component of gradient at "Z" side
                  grad_z_zside(i,j,k,m) = dzinv * diff_z(i,j,k,m)

               enddo
            enddo
         enddo

      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine quatgrad_side_isotropic(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth, h,
     &   diff_x, diff_y, diff_z, ngdiff,
     &   grad_x_xside, grad_y_xside, grad_z_xside,
     &   grad_x_yside, grad_y_yside, grad_z_yside,
     &   grad_x_zside, grad_y_zside, grad_z_zside,
     &   nggrad
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth, ngdiff, nggrad

      double precision diff_x(SIDE3d0(lo,hi,ngdiff),depth)
      double precision diff_y(SIDE3d1(lo,hi,ngdiff),depth)
      double precision diff_z(SIDE3d2(lo,hi,ngdiff),depth)
      double precision grad_x_xside(SIDE3d0(lo,hi,nggrad),depth)
      double precision grad_x_yside(SIDE3d1(lo,hi,nggrad),depth)
      double precision grad_x_zside(SIDE3d2(lo,hi,nggrad),depth)
      double precision grad_y_xside(SIDE3d0(lo,hi,nggrad),depth)
      double precision grad_y_yside(SIDE3d1(lo,hi,nggrad),depth)
      double precision grad_y_zside(SIDE3d2(lo,hi,nggrad),depth)
      double precision grad_z_xside(SIDE3d0(lo,hi,nggrad),depth)
      double precision grad_z_yside(SIDE3d1(lo,hi,nggrad),depth)
      double precision grad_z_zside(SIDE3d2(lo,hi,nggrad),depth)
      double precision h(3)

      print*,'ERROR in function quatgrad_side_isotropic()'
      print*,'Not implemented in 3D!!!'
      stop
      
      return
      end
      
c-----------------------------------------------------------------------

      subroutine quatgrad_side_symm(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth, h,
     &   diff_x, diff_y, diff_z, ngdiff,
     &   grad_x_xside, grad_y_xside, grad_z_xside,
     &   grad_x_yside, grad_y_yside, grad_z_yside,
     &   grad_x_zside, grad_y_zside, grad_z_zside,
     &   nggrad,
     &   iqrot_x, iqrot_y, iqrot_z, ngiq
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth, ngdiff, nggrad, ngiq
      double precision h(3)

      double precision diff_x(SIDE3d0(lo,hi,ngdiff),depth)
      double precision diff_y(SIDE3d1(lo,hi,ngdiff),depth)
      double precision diff_z(SIDE3d2(lo,hi,ngdiff),depth)
      double precision grad_x_xside(SIDE3d0(lo,hi,nggrad),depth)
      double precision grad_x_yside(SIDE3d1(lo,hi,nggrad),depth)
      double precision grad_x_zside(SIDE3d2(lo,hi,nggrad),depth)
      double precision grad_y_xside(SIDE3d0(lo,hi,nggrad),depth)
      double precision grad_y_yside(SIDE3d1(lo,hi,nggrad),depth)
      double precision grad_y_zside(SIDE3d2(lo,hi,nggrad),depth)
      double precision grad_z_xside(SIDE3d0(lo,hi,nggrad),depth)
      double precision grad_z_yside(SIDE3d1(lo,hi,nggrad),depth)
      double precision grad_z_zside(SIDE3d2(lo,hi,nggrad),depth)

      integer iqrot_x(SIDE3d0(lo,hi,ngiq))
      integer iqrot_y(SIDE3d1(lo,hi,ngiq))
      integer iqrot_z(SIDE3d2(lo,hi,ngiq))

c        local variables:
      integer i, j, k, m, iq
      double precision dxinv, dyinv, dzinv
      double precision p25_dxinv, p25_dyinv, p25_dzinv
      double precision d1(depth), d1prime(depth)
      double precision d2(depth), d2prime(depth)
      double precision d3(depth), d3prime(depth)
      double precision d4(depth), d4prime(depth)

      dxinv = 1.d0 / h(1)
      dyinv = 1.d0 / h(2)
      dzinv = 1.d0 / h(3)
      p25_dxinv = 0.25d0 * dxinv
      p25_dyinv = 0.25d0 * dyinv
      p25_dzinv = 0.25d0 * dzinv

c        All diffs must be rotated to be consistent with the
c        local quaternion.

c        Compute gradients on the "X" side of the cell
      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0+1

c                 X component of gradient at "X" side
               do m = 1, depth
                  grad_x_xside(i,j,k,m) = dxinv * diff_x(i,j,k,m)
               enddo

               if ( depth > 1 ) then

                  do m = 1, depth
                     d1(m) = diff_y(i,  j+1,k,m)
                     d2(m) = diff_y(i-1,j+1,k,m)
                     d3(m) = diff_y(i-1,j,  k,m)
                  enddo

                  iq = -1 * iqrot_y(i,j+1,k)
                  call quatsymmrotate( d1, iq, d1prime, depth )

                  iq = -1 * iqrot_y(i-1,j+1,k)
                  call quatsymmrotate( d2, iq, d2prime, depth )

                  do m = 1, depth
                     d4(m) = d2prime(m) + d3(m)
                  enddo

                  iq = iqrot_x(i,j,k)
                  call quatsymmrotate( d4, iq, d4prime, depth )

               else

                  d1prime(1) = diff_y(i,j+1,k,1)
                  d4prime(1) = diff_y(i-1,j+1,k,1) + diff_y(i-1,j,k,1)

               endif

c                 Y component of gradient at "X" side
               do m = 1, depth
                  grad_y_xside(i,j,k,m) = p25_dyinv * (
     &               d4prime(m) + d1prime(m) + diff_y(i,j,k,m) )
               enddo

               if ( depth > 1 ) then

                  do m = 1, depth
                     d1(m) = diff_z(i,  j,k+1,m)
                     d2(m) = diff_z(i-1,j,k+1,m)
                     d3(m) = diff_z(i-1,j,k,  m)
                  enddo

                  iq = -1 * iqrot_z(i,j,k+1)
                  call quatsymmrotate( d1, iq, d1prime, depth )

                  iq = -1 * iqrot_z(i-1,j,k+1)
                  call quatsymmrotate( d2, iq, d2prime, depth )

                  do m = 1, depth
                     d4(m) = d2prime(m) + d3(m)
                  enddo

                  iq = iqrot_x(i,j,k)
                  call quatsymmrotate( d4, iq, d4prime, depth )

               else

                  d1prime(1) = diff_z(i,  j,k+1,1)
                  d4prime(1) = diff_z(i-1,j,k+1,1) + diff_z(i-1,j,k,1)

               endif

c                 Z component of gradient at "X" side
               do m = 1, depth
                  grad_z_xside(i,j,k,m) = p25_dzinv * (
     &               d4prime(m) + d1prime(m) + diff_z(i,j,k,m) )
               enddo

            enddo
         enddo
      enddo

c        Compute gradients on the "Y" side of the cell
      do k = lo2, hi2
         do j = lo1, hi1+1
            do i = lo0, hi0

c                 Y component of gradient at "Y" side
               do m = 1, depth
                  grad_y_yside(i,j,k,m) = dyinv * diff_y(i,j,k,m)
               enddo

               if ( depth > 1 ) then

                  do m = 1, depth
                     d1(m) = diff_x(i+1,j,  k,m)
                     d2(m) = diff_x(i+1,j-1,k,m)
                     d3(m) = diff_x(i,  j-1,k,m)
                  enddo

                  iq = -1 * iqrot_x(i+1,j,k)
                  call quatsymmrotate( d1, iq, d1prime, depth )

                  iq = -1 * iqrot_x(i+1,j-1,k)
                  call quatsymmrotate( d2, iq, d2prime, depth )

                  do m = 1, depth
                     d4(m) = d2prime(m) + d3(m)
                  enddo

                  iq = iqrot_y(i,j,k)
                  call quatsymmrotate( d4, iq, d4prime, depth )

               else

                  d1prime(1) = diff_x(i+1,j,  k,1)
                  d4prime(1) = diff_x(i+1,j-1,k,1) + diff_x(i,j-1,k,1)

               endif

c                 X component of gradient at "Y" side
               do m = 1, depth
                  grad_x_yside(i,j,k,m) = p25_dxinv * (
     &               d4prime(m) + d1prime(m) + diff_x(i,j,k,m) )
               enddo

               if ( depth > 1 ) then

                  do m = 1, depth
                     d1(m) = diff_z(i,j,  k+1,m)
                     d2(m) = diff_z(i,j-1,k+1,m)
                     d3(m) = diff_z(i,j-1,k,  m)
                  enddo

                  iq = -1 * iqrot_z(i,j,k+1)
                  call quatsymmrotate( d1, iq, d1prime, depth )

                  iq = -1 * iqrot_z(i,j-1,k+1)
                  call quatsymmrotate( d2, iq, d2prime, depth )

                  do m = 1, depth
                     d4(m) = d2prime(m) + d3(m)
                  enddo

                  iq = iqrot_y(i,j,k)
                  call quatsymmrotate( d4, iq, d4prime, depth )

               else

                  d1prime(1) = diff_z(i,j,k+1,1)
                  d4prime(1) = diff_z(i,j-1,k+1,1) + diff_z(i,j-1,k,1)

               endif

c                 Z component of gradient at "Y" side
               do m = 1, depth
                  grad_z_yside(i,j,k,m) = p25_dzinv * (
     &               d4prime(m) + d1prime(m) + diff_z(i,j,k,m) )
               enddo

            enddo
         enddo
      enddo

c        Compute gradients on the "Z" side of the cell
      do k = lo2, hi2+1
         do j = lo1, hi1
            do i = lo0, hi0

c                 Z component of gradient at "Z" side
               do m = 1, depth
                  grad_z_zside(i,j,k,m) = dzinv * diff_z(i,j,k,m)
               enddo

               if ( depth > 1 ) then

                  do m = 1, depth
                     d1(m) = diff_x(i+1,j,k,  m)
                     d2(m) = diff_x(i+1,j,k-1,m)
                     d3(m) = diff_x(i,  j,k-1,m)
                  enddo

                  iq = -1 * iqrot_x(i+1,j,k)
                  call quatsymmrotate( d1, iq, d1prime, depth )

                  iq = -1 * iqrot_x(i+1,j,k-1)
                  call quatsymmrotate( d2, iq, d2prime, depth )

                  do m = 1, depth
                     d4(m) = d2prime(m) + d3(m)
                  enddo

                  iq = iqrot_z(i,j,k)
                  call quatsymmrotate( d4, iq, d4prime, depth )

               else

                  d1prime(1) = diff_x(i+1,j,k,1)
                  d4prime(1) = diff_x(i+1,j,k-1,1) + diff_x(i,j,k-1,1)

               endif

c                 X component of gradient at "Z" side
               do m = 1, depth
                  grad_x_zside(i,j,k,m) = p25_dxinv * (
     &               d4prime(m) + d1prime(m) + diff_x(i,j,k,m) )
               enddo

               if ( depth > 1 ) then

                  do m = 1, depth
                     d1(m) = diff_y(i,j+1,k,  m)
                     d2(m) = diff_y(i,j+1,k-1,m)
                     d3(m) = diff_y(i,j,  k-1,m)
                  enddo

                  iq = -1 * iqrot_y(i,j+1,k)
                  call quatsymmrotate( d1, iq, d1prime, depth )

                  iq = -1 * iqrot_y(i,j+1,k-1)
                  call quatsymmrotate( d2, iq, d2prime, depth )

                  do m = 1, depth
                     d4(m) = d2prime(m) + d3(m)
                  enddo

                  iq = iqrot_z(i,j,k)
                  call quatsymmrotate( d4, iq, d4prime, depth )

               else

                  d1prime(1) = diff_y(i,j+1,k,  1)
                  d4prime(1) = diff_y(i,j+1,k-1,1) + diff_y(i,j,k-1,1)

               endif

c                 Y component of gradient at "Z" side
               do m = 1, depth
                  grad_y_zside(i,j,k,m) = p25_dyinv * (
     &               d4prime(m) + d1prime(m) + diff_y(i,j,k,m) )
               enddo

            enddo
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine quatgrad_side_symm_isotropic(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth, h,
     &   diff_x, diff_y, diff_z, ngdiff,
     &   grad_x_xside, grad_y_xside, grad_z_xside,
     &   grad_x_yside, grad_y_yside, grad_z_yside,
     &   grad_x_zside, grad_y_zside, grad_z_zside,
     &   nggrad,
     &   iqrot_x, iqrot_y, iqrot_z, ngiq
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth, ngdiff, nggrad, ngiq
      double precision h(3)

      double precision diff_x(SIDE3d0(lo,hi,ngdiff),depth)
      double precision diff_y(SIDE3d1(lo,hi,ngdiff),depth)
      double precision diff_z(SIDE3d2(lo,hi,ngdiff),depth)
      double precision grad_x_xside(SIDE3d0(lo,hi,nggrad),depth)
      double precision grad_x_yside(SIDE3d1(lo,hi,nggrad),depth)
      double precision grad_x_zside(SIDE3d2(lo,hi,nggrad),depth)
      double precision grad_y_xside(SIDE3d0(lo,hi,nggrad),depth)
      double precision grad_y_yside(SIDE3d1(lo,hi,nggrad),depth)
      double precision grad_y_zside(SIDE3d2(lo,hi,nggrad),depth)
      double precision grad_z_xside(SIDE3d0(lo,hi,nggrad),depth)
      double precision grad_z_yside(SIDE3d1(lo,hi,nggrad),depth)
      double precision grad_z_zside(SIDE3d2(lo,hi,nggrad),depth)

      integer iqrot_x(SIDE3d0(lo,hi,ngiq))
      integer iqrot_y(SIDE3d1(lo,hi,ngiq))
      integer iqrot_z(SIDE3d2(lo,hi,ngiq))

      print*,'ERROR in function quatgrad_side_symm_isotropic()'
      print*,'Not implemented in 3D!!!'
      stop
      
      return
      end

c-----------------------------------------------------------------------

      subroutine quatgrad_modulus_from_sides_compact(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth,
     &   grad_x, grad_y, grad_z, nggq,
     &   grad_mod, ngm)

      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth, nggq, ngm

      double precision
     &   grad_x(SIDE3d0(lo,hi,nggq),depth,NDIM),
     &   grad_y(SIDE3d1(lo,hi,nggq),depth,NDIM),
     &   grad_z(SIDE3d2(lo,hi,nggq),depth,NDIM),
     &   grad_mod(CELL3d(lo,hi,ngm))

      double precision grad_modulus

      integer i, j, k, m, n
      
      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0

               grad_mod(i,j,k) = 0.0d0

c we average contributions from 4 sides
            
               grad_modulus = 0.d0
               do m = 1, depth
                  grad_modulus = grad_modulus + 
     &               grad_x(i,j,k,m,1)*grad_x(i,j,k,m,1)
               enddo
               do m = 1, depth
                  grad_modulus = grad_modulus + 
     &               grad_x(i+1,j,k,m,1)*grad_x(i+1,j,k,m,1)
               enddo
               do m = 1, depth
                  grad_modulus = grad_modulus + 
     &               grad_y(i,j,k,m,2)*grad_y(i,j,k,m,2)
               enddo
               do m = 1, depth
                  grad_modulus = grad_modulus + 
     &               grad_y(i,j+1,k,m,2)*grad_y(i,j+1,k,m,2)
               enddo
               do m = 1, depth
                  grad_modulus = grad_modulus + 
     &               grad_z(i,j,k,m,3)*grad_z(i,j,k,m,3)
               enddo
               do m = 1, depth
                  grad_modulus = grad_modulus + 
     &               grad_z(i,j,k+1,m,3)*grad_z(i,j,k+1,m,3)
               enddo

               grad_mod(i,j,k) = dsqrt(0.5d0*grad_modulus)
            
            enddo
         enddo
      enddo

      return
      end


c-----------------------------------------------------------------------

      subroutine quatgrad_modulus(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth,
     &   grad_x, grad_y, grad_z, nggq,
     &   grad_mod, ngm)

      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth, nggq, ngm

      double precision
     &   grad_x(CELL3d(lo,hi,nggq),depth),
     &   grad_y(CELL3d(lo,hi,nggq),depth),
     &   grad_z(CELL3d(lo,hi,nggq),depth),
     &   grad_mod(CELL3d(lo,hi,ngm))

      integer i, j, k, m

      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0

               grad_mod(i,j,k) = 0.0d0

               do m = 1, depth
                  grad_mod(i,j,k) = grad_mod(i,j,k) +
     &               grad_x(i,j,k,m) * grad_x(i,j,k,m) +
     &               grad_y(i,j,k,m) * grad_y(i,j,k,m) +
     &               grad_z(i,j,k,m) * grad_z(i,j,k,m)
               enddo
            
               grad_mod(i,j,k) = dsqrt( grad_mod(i,j,k) )

            enddo
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
