c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
c Redistribution and use in source and binary forms, with or without 
c modification, are permitted provided that the following conditions are met:
c - Redistributions of source code must retain the above copyright notice,
c   this list of conditions and the disclaimer below.
c - Redistributions in binary form must reproduce the above copyright notice,
c   this list of conditions and the disclaimer (as noted below) in the
c   documentation and/or other materials provided with the distribution.
c - Neither the name of the LLNS/LLNL nor the names of its contributors may be
c   used to endorse or promote products derived from this software without
c   specific prior written permission.
c
c THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
c AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
c IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
c ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
c LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
c DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
c DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
c OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
c HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
c STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
c IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
c POSSIBILITY OF SUCH DAMAGE.
c 
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

      subroutine quatgrad_side_symm_wide(
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

      print*,'function quatgrad_side_symm_wide not implemented in 3D!!!'
      
      return
      end

c-----------------------------------------------------------------------

      subroutine quatgrad_modulus_from_sides_compact(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth,
     &   grad_x, grad_y, 
     &   glo0, ghi0, glo1, ghi1, glo2, ghi2,
     &   grad_mod,
     &   mlo0, mhi0, mlo1, mhi1, mlo2, mhi2,
     &   floor_type,
     &   grad_floor )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth,
     &   glo0, ghi0, glo1, ghi1, glo2, ghi2,
     &   mlo0, mhi0, mlo1, mhi1, mlo2, mhi2

      double precision
     &   grad_x(glo0:ghi0+1,glo1:ghi1,  glo2:ghi2,  depth,NDIM),
     &   grad_y(glo0:ghi0,  glo1:ghi1+1,glo2:ghi2,  depth,NDIM),
     &   grad_z(glo0:ghi0,  glo1:ghi1,  glo2:ghi2+1,depth,NDIM),
     &   grad_mod(mlo0:mhi0,mlo1:mhi1,mlo2:mhi2)

      double precision grad_floor, beta, grad_modulus

      character*(*) floor_type

      integer i, j, k, m, n
      
      if( floor_type(1:1) .eq. 's' )then
         beta = grad_floor*grad_floor
      else
         beta = 0.d0
      endif

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

               grad_mod(i,j,k) = dsqrt(0.5d0*grad_modulus+beta)
            
            enddo
         enddo
      enddo

      return
      end


c-----------------------------------------------------------------------

      subroutine quatgrad_modulus(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth,
     &   grad_x, grad_y, grad_z,
     &   glo0, ghi0, glo1, ghi1, glo2, ghi2,
     &   grad_mod,
     &   mlo0, mhi0, mlo1, mhi1, mlo2, mhi2,
     &   floor_type,
     &   grad_floor )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth,
     &   glo0, ghi0, glo1, ghi1, glo2, ghi2,
     &   mlo0, mhi0, mlo1, mhi1, mlo2, mhi2

      double precision
     &   grad_x(glo0:ghi0,glo1:ghi1,glo2:ghi2,depth),
     &   grad_y(glo0:ghi0,glo1:ghi1,glo2:ghi2,depth),
     &   grad_z(glo0:ghi0,glo1:ghi1,glo2:ghi2,depth),
     &   grad_mod(mlo0:mhi0,mlo1:mhi1,mlo2:mhi2)

      double precision grad_floor

      character*(*) floor_type

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
            
               if( floor_type(1:1) .eq. 's' )then
                  grad_mod(i,j,k) = grad_mod(i,j,k) + grad_floor
               endif
               grad_mod(i,j,k) = dsqrt( grad_mod(i,j,k) )

            enddo
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
