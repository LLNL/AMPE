c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
define(NDIM,2)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl


      subroutine quatgrad_cell(
     &   lo0, hi0, lo1, hi1,
     &   depth, h,
     &   diff_x, diff_y, ngdiff,
     &   grad_x, grad_y, nggrad
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1,
     &   depth, ngdiff, nggrad

      double precision diff_x(SIDE2d0(lo,hi,ngdiff),depth)
      double precision diff_y(SIDE2d1(lo,hi,ngdiff),depth)
      double precision grad_x(CELL2d(lo,hi,nggrad),depth)
      double precision grad_y(CELL2d(lo,hi,nggrad),depth)
      double precision h(2)

      integer i, j, m
      double precision p5dxinv, p5dyinv

      p5dxinv = 0.5d0 / h(1)
      p5dyinv = 0.5d0 / h(2)
      
      do m = 1, depth

         do j = lo1, hi1
            do i = lo0, hi0
               grad_x(i,j,m) =
     &            ( diff_x(i+1,j,m) + diff_x(i,j,m) ) * p5dxinv
            enddo
         enddo

         do j = lo1, hi1
            do i = lo0, hi0
               grad_y(i,j,m) =
     &            ( diff_y(i,j+1,m) + diff_y(i,j,m) ) * p5dyinv
            enddo
         enddo

      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine quatgrad_cell_symm(
     &   lo0, hi0, lo1, hi1,
     &   depth, h,
     &   diff_x, diff_y, ngdiff,
     &   grad_x, grad_y, nggrad,
     &   iqrot_x, iqrot_y, ngiq
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1,
     &   depth, ngdiff, nggrad, ngiq
      double precision h(2)

      double precision diff_x(SIDE2d0(lo,hi,ngdiff),depth)
      double precision diff_y(SIDE2d1(lo,hi,ngdiff),depth)
      double precision grad_x(CELL2d(lo,hi,nggrad),depth)
      double precision grad_y(CELL2d(lo,hi,nggrad),depth)

      integer iqrot_x(SIDE2d0(lo,hi,ngiq))
      integer iqrot_y(SIDE2d1(lo,hi,ngiq))

      integer i, j, m, iq
      double precision p5dxinv, p5dyinv
      double precision dtmp(depth)
      double precision dprime(depth)

      p5dxinv = 0.5d0 / h(1)
      p5dyinv = 0.5d0 / h(2)

      do j = lo1, hi1
         do i = lo0, hi0

            if ( depth > 1 ) then

c                 Diffs on right side need to have the conjugate
c                 of the symmetry rotation applied before using, but
c                 only for the vector/quaternion orientation.
c                 This is required to have "compatible" diffs
c                 one can average.
               do m = 1, depth
                  dtmp(m) = diff_x(i+1,j,m)
               enddo
               
               iq = -1 * iqrot_x(i+1,j)
               call quatsymmrotate( dtmp, iq, dprime, depth )

            else

               dprime(1) = diff_x(i+1,j,1)

            endif

            do m = 1, depth
               grad_x(i,j,m) =
     &            ( dprime(m) + diff_x(i,j,m) ) * p5dxinv
            enddo

         enddo
      enddo

      do j = lo1, hi1
         do i = lo0, hi0

            if ( depth > 1 ) then

c                 Diffs on top side need to have the conjugate
c                 of the symmetry rotation applied before using, but
c                 only for the vector/quaternion orientation.

               do m = 1, depth
                  dtmp(m) = diff_y(i,j+1,m)
               enddo
               
               iq = -1 * iqrot_y(i,j+1)
               call quatsymmrotate( dtmp, iq, dprime, depth )

            else

               dprime(1) = diff_y(i,j+1,1)

            endif

            do m = 1, depth
               grad_y(i,j,m) =
     &            ( dprime(m) + diff_y(i,j,m) ) * p5dyinv
            enddo

         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine quatgrad_side(
     &   lo0, hi0, lo1, hi1,
     &   depth, h,
     &   diff_x, diff_y, ngdiff,
     &   grad_x_xside, grad_y_xside,
     &   grad_x_yside, grad_y_yside,
     &   nggrad
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1,
     &   depth, ngdiff, nggrad

      double precision diff_x(SIDE2d0(lo,hi,ngdiff),depth)
      double precision diff_y(SIDE2d1(lo,hi,ngdiff),depth)
      double precision grad_x_xside(SIDE2d0(lo,hi,nggrad),depth)
      double precision grad_x_yside(SIDE2d1(lo,hi,nggrad),depth)
      double precision grad_y_xside(SIDE2d0(lo,hi,nggrad),depth)
      double precision grad_y_yside(SIDE2d1(lo,hi,nggrad),depth)
      double precision h(2)

c        local variables:
      integer i, j, m
      double precision dxinv, dyinv
      double precision p25_dxinv, p25_dyinv

      dxinv = 1.d0 / h(1)
      dyinv = 1.d0 / h(2)
      p25_dxinv = 0.25d0 * dxinv
      p25_dyinv = 0.25d0 * dyinv

c        Loop over the quaternion components
      do m = 1, depth

c           Compute gradients on the "X" side of the cell
         do j = lo1, hi1
            do i = lo0, hi0+1

c                 X component of gradient at "X" side
               grad_x_xside(i,j,m) = dxinv * diff_x(i,j,m)

c                 Y component of gradient at "X" side
               grad_y_xside(i,j,m) = p25_dyinv * (
     &            diff_y(i-1,j+1,m) + diff_y(i-1,j,m) +
     &            diff_y(i,  j+1,m) + diff_y(i,  j,m)
     &            ) 

            enddo
         enddo

c           Compute gradients on the "Y" side of the cell
         do j = lo1, hi1+1
            do i = lo0, hi0

c                 X component of gradient at "Y" side
               grad_x_yside(i,j,m) = p25_dxinv * (
     &            diff_x(i+1,j-1,m) + diff_x(i,j-1,m) +
     &            diff_x(i+1,j,  m) + diff_x(i,j,  m)
     &            ) 

c                 Y component of gradient at "Y" side
               grad_y_yside(i,j,m) = dyinv * diff_y(i,j,m)

            enddo
         enddo

      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine quatgrad_side_isotropic(
     &   lo0, hi0, lo1, hi1,
     &   depth, h,
     &   diff_x, diff_y, ngdiff,
     &   grad_x_xside, grad_y_xside,
     &   grad_x_yside, grad_y_yside,
     &   nggrad
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1,
     &   depth, ngdiff, nggrad

      double precision diff_x(SIDE2d0(lo,hi,ngdiff),depth)
      double precision diff_y(SIDE2d1(lo,hi,ngdiff),depth)
      double precision grad_x_xside(SIDE2d0(lo,hi,nggrad),depth)
      double precision grad_x_yside(SIDE2d1(lo,hi,nggrad),depth)
      double precision grad_y_xside(SIDE2d0(lo,hi,nggrad),depth)
      double precision grad_y_yside(SIDE2d1(lo,hi,nggrad),depth)
      double precision h(2)

c        local variables:
      integer i, j, m
      double precision dxinv, dyinv
      double precision p25_dxinv, p25_dyinv, p12_dxinv, p12_dyinv

      dxinv = 1.d0 / h(1)
      dyinv = 1.d0 / h(2)
      p25_dxinv = 0.25d0 * dxinv
      p12_dxinv = dxinv/12.d0
      p12_dyinv = dyinv/12.d0
      p25_dyinv = 0.25d0 * dyinv

c        Loop over the quaternion components
      do m = 1, depth

c           Compute gradients on the "X" side of the cell
         do j = lo1, hi1
            do i = lo0, hi0+1

c                 X component of gradient at "X" side
               grad_x_xside(i,j,m) =p12_dxinv 
     &            *(10.d0*diff_x(i,j,m)+diff_x(i,j-1,m)+diff_x(i,j+1,m))

c                 Y component of gradient at "X" side
               grad_y_xside(i,j,m) = p25_dyinv * (
     &            diff_y(i-1,j+1,m) + diff_y(i-1,j,m) +
     &            diff_y(i,  j+1,m) + diff_y(i,  j,m)
     &            ) 

            enddo
         enddo

c           Compute gradients on the "Y" side of the cell
         do j = lo1, hi1+1
            do i = lo0, hi0

c                 X component of gradient at "Y" side
               grad_x_yside(i,j,m) = p25_dxinv * (
     &            diff_x(i+1,j-1,m) + diff_x(i,j-1,m) +
     &            diff_x(i+1,j,  m) + diff_x(i,j,  m)
     &            ) 

c                 Y component of gradient at "Y" side
               grad_y_yside(i,j,m) = p12_dyinv 
     &            *(10.d0*diff_y(i,j,m)+diff_y(i-1,j,m)+diff_y(i+1,j,m))

            enddo
         enddo

      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine quatgrad_side_symm(
     &   lo0, hi0, lo1, hi1,
     &   depth, h,
     &   diff_x, diff_y, ngdiff,
     &   grad_x_xside, grad_y_xside,
     &   grad_x_yside, grad_y_yside,
     &   nggrad,
     &   iqrot_x, iqrot_y, ngiq
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1,
     &   depth, ngdiff, nggrad, ngiq
      double precision h(2)

      double precision diff_x(SIDE2d0(lo,hi,ngdiff),depth)
      double precision diff_y(SIDE2d1(lo,hi,ngdiff),depth)
      double precision grad_x_xside(SIDE2d0(lo,hi,nggrad),depth)
      double precision grad_x_yside(SIDE2d1(lo,hi,nggrad),depth)
      double precision grad_y_xside(SIDE2d0(lo,hi,nggrad),depth)
      double precision grad_y_yside(SIDE2d1(lo,hi,nggrad),depth)

      integer iqrot_x(SIDE2d0(lo,hi,ngiq))
      integer iqrot_y(SIDE2d1(lo,hi,ngiq))

c        local variables:
      integer i, j, m, iq
      double precision dxinv, dyinv
      double precision p25_dxinv, p25_dyinv
      double precision d1(depth), d1prime(depth)
      double precision d2(depth), d2prime(depth)
      double precision d3(depth)
      double precision d4(depth), d4prime(depth)

      dxinv = 1.d0 / h(1)
      dyinv = 1.d0 / h(2)
      p25_dxinv = 0.25d0 * dxinv
      p25_dyinv = 0.25d0 * dyinv

c        All diffs must be rotated to be consistent with the
c        local quaternion.

c        Compute gradients on the "X" side of the cell
      do j = lo1, hi1
         do i = lo0, hi0+1

            if ( depth > 1 ) then

               do m = 1, depth
                  d1(m) = diff_y(i,  j+1,m)
                  d2(m) = diff_y(i-1,j+1,m)
                  d3(m) = diff_y(i-1,j,  m)
               enddo

               iq = -1 * iqrot_y(i,j+1)
               call quatsymmrotate( d1, iq, d1prime, depth )

               iq = -1 * iqrot_y(i-1,j+1)
               call quatsymmrotate( d2, iq, d2prime, depth )

               do m = 1, depth
                  d4(m) = d2prime(m) + d3(m)
               enddo

               iq = iqrot_x(i,j)
               call quatsymmrotate( d4, iq, d4prime, depth )
               
            else

               d1prime(1) = diff_y(i,  j+1,1)
               d4prime(1) = diff_y(i-1,j+1,1) + diff_y(i-1,j,1)

            endif

            do m = 1, depth

c                 Y component of gradient at "X" side
               grad_y_xside(i,j,m) = p25_dyinv * (
     &            d4prime(m) + d1prime(m) + diff_y(i,j,m) )

c                 X component of gradient at "X" side
               grad_x_xside(i,j,m) = dxinv * diff_x(i,j,m)

            enddo

         enddo
      enddo

c        Compute gradients on the "Y" side of the cell
      do j = lo1, hi1+1
         do i = lo0, hi0

            if ( depth > 1 ) then

               do m = 1, depth
                  d1(m) = diff_x(i+1,j,  m)
                  d2(m) = diff_x(i+1,j-1,m)
                  d3(m) = diff_x(i,  j-1,m)
               enddo

               iq = -1 * iqrot_x(i+1,j)
               call quatsymmrotate( d1, iq, d1prime, depth )

               iq = -1 * iqrot_x(i+1,j-1)
               call quatsymmrotate( d2, iq, d2prime, depth )

               do m = 1, depth
                  d4(m) = d2prime(m) + d3(m)
               enddo

               iq = iqrot_y(i,j)
               call quatsymmrotate( d4, iq, d4prime, depth )
               
            else

               d1prime(1) = diff_x(i+1,j,  1)
               d4prime(1) = diff_x(i+1,j-1,1) + diff_x(i,j-1,1)

            endif

            do m = 1, depth

c                 X component of gradient at "Y" side
               grad_x_yside(i,j,m) = p25_dxinv * (
     &            d4prime(m) + d1prime(m) + diff_x(i,j,m) )

c                 Y component of gradient at "Y" side
               grad_y_yside(i,j,m) = dyinv * diff_y(i,j,m)

            enddo

         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c Isotropic scheme for dx on face x and dy on face y
c see Shukla and Giri, J. Comput. Phys. 276 (2014), p.259
c
      subroutine quatgrad_side_symm_isotropic(
     &   lo0, hi0, lo1, hi1,
     &   depth, h,
     &   diff_x, diff_y, ngdiff,
     &   grad_x_xside, grad_y_xside,
     &   grad_x_yside, grad_y_yside,
     &   nggrad,
     &   iqrot_x, iqrot_y, ngiq
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1,
     &   depth, ngdiff, nggrad, ngiq
      double precision h(2)

      double precision diff_x(SIDE2d0(lo,hi,ngdiff),depth)
      double precision diff_y(SIDE2d1(lo,hi,ngdiff),depth)
      double precision grad_x_xside(SIDE2d0(lo,hi,nggrad),depth)
      double precision grad_x_yside(SIDE2d1(lo,hi,nggrad),depth)
      double precision grad_y_xside(SIDE2d0(lo,hi,nggrad),depth)
      double precision grad_y_yside(SIDE2d1(lo,hi,nggrad),depth)

      integer iqrot_x(SIDE2d0(lo,hi,ngiq))
      integer iqrot_y(SIDE2d1(lo,hi,ngiq))

c        local variables:
      integer i, j, m, iq
      double precision dxinv, dyinv
      double precision p25_dxinv, p25_dyinv, p12_dxinv, p12_dyinv
      double precision d1(depth), d1prime(depth)
      double precision d2(depth), d2prime(depth)
      double precision d3(depth)
      double precision d4(depth), d4prime(depth)

c      print*,'call quatgrad_side_symm_wide'

      dxinv = 1.d0 / h(1)
      dyinv = 1.d0 / h(2)
      p25_dxinv = 0.25d0 * dxinv
      p25_dyinv = 0.25d0 * dyinv
      p12_dxinv = dxinv/12.d0
      p12_dyinv = dyinv/12.d0

c        All diffs must be rotated to be consistent with the
c        local quaternion.

c        Compute gradients on the "X" side of the cell
      do j = lo1, hi1
         do i = lo0, hi0+1

            if ( depth > 1 ) then

               do m = 1, depth
                  d1(m) = diff_y(i,  j+1,m)
                  d2(m) = diff_y(i-1,j+1,m)
                  d3(m) = diff_y(i-1,j,  m)
               enddo

               iq = -1 * iqrot_y(i,j+1)
               call quatsymmrotate( d1, iq, d1prime, depth )

               iq = -1 * iqrot_y(i-1,j+1)
               call quatsymmrotate( d2, iq, d2prime, depth )

               do m = 1, depth
                  d4(m) = d2prime(m) + d3(m)
               enddo

               iq = iqrot_x(i,j)
               call quatsymmrotate( d4, iq, d4prime, depth )
               
            else

               d1prime(1) = diff_y(i,  j+1,1)
               d4prime(1) = diff_y(i-1,j+1,1) + diff_y(i-1,j,1)

            endif

            do m = 1, depth

c                 Y component of gradient at "X" side
               grad_y_xside(i,j,m) = p25_dyinv * (
     &            d4prime(m) + d1prime(m) + diff_y(i,j,m) )
            enddo


            if ( depth > 1 ) then

               do m = 1, depth
                  d1(m) = diff_x(i,j-1,m)
                  d2(m) = diff_x(i,j+1,m)
               enddo

               iq = iqrot_y(i,j)
               call quatsymmrotate( d1, iq, d1prime, depth )

               iq = -1 * iqrot_y(i,j+1)
               call quatsymmrotate( d2, iq, d2prime, depth )
               
            else

               d1prime(1) = diff_x(i,j-1,m)
               d2prime(1) = diff_x(i,j+1,m)

            endif

            do m = 1, depth

c                 X component of gradient at "X" side
               grad_x_xside(i,j,m) = p12_dxinv 
     &           * (10.d0*diff_x(i,j,m)+d1prime(m)+d2prime(m))

            enddo

         enddo
      enddo

c        Compute gradients on the "Y" side of the cell
      do j = lo1, hi1+1
         do i = lo0, hi0

            if ( depth > 1 ) then

               do m = 1, depth
                  d1(m) = diff_x(i+1,j,  m)
                  d2(m) = diff_x(i+1,j-1,m)
                  d3(m) = diff_x(i,  j-1,m)
               enddo

               iq = -1 * iqrot_x(i+1,j)
               call quatsymmrotate( d1, iq, d1prime, depth )

               iq = -1 * iqrot_x(i+1,j-1)
               call quatsymmrotate( d2, iq, d2prime, depth )

               do m = 1, depth
                  d4(m) = d2prime(m) + d3(m)
               enddo

               iq = iqrot_y(i,j)
               call quatsymmrotate( d4, iq, d4prime, depth )
               
            else

               d1prime(1) = diff_x(i+1,j,  1)
               d4prime(1) = diff_x(i+1,j-1,1) + diff_x(i,j-1,1)

            endif

            do m = 1, depth

c                 X component of gradient at "Y" side
               grad_x_yside(i,j,m) = p25_dxinv * (
     &            d4prime(m) + d1prime(m) + diff_x(i,j,m) )

            enddo

            if ( depth > 1 ) then

               do m = 1, depth
                  d1(m) = diff_y(i-1,j,m)
                  d2(m) = diff_y(i+1,j,m)
               enddo

               iq = iqrot_x(i,j)
               call quatsymmrotate( d1, iq, d1prime, depth )

               iq = -1 * iqrot_x(i+1,j)
               call quatsymmrotate( d2, iq, d2prime, depth )

            else

               d1prime(1) = diff_y(i-1,j,m)
               d2prime(1) = diff_y(i+1,j,m)

            endif

            do m = 1, depth

c                 Y component of gradient at "Y" side
               grad_y_yside(i,j,m) = p12_dyinv 
     &           * (10.d0*diff_y(i,j,m)+d1prime(m)+d2prime(m))

            enddo

         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine quatgrad_modulus_from_sides(
     &   lo0, hi0, lo1, hi1,
     &   depth,
     &   grad_x, grad_y, nggq,
     &   grad_mod, ngm)

      implicit none

      integer
     &   lo0, hi0, lo1, hi1,
     &   depth, nggq, ngm

      double precision
     &   grad_x(SIDE2d0(lo,hi,nggq),depth,NDIM),
     &   grad_y(SIDE2d1(lo,hi,nggq),depth,NDIM),
     &   grad_mod(CELL2d(lo,hi,ngm))

      double precision grad_modulus

      integer i, j, q, d
      
      do j = lo1, hi1
         do i = lo0, hi0

            grad_mod(i,j) = 0.0d0

c we average contributions from 4 sides
            
            grad_modulus = 0.d0
            do q = 1, depth
c loop over components of nabla q_m
               do d = 1, NDIM
                  grad_modulus = grad_modulus + 
     &               grad_x(i,j,q,d)*grad_x(i,j,q,d)
               enddo
            enddo
            grad_mod(i,j) = grad_mod(i,j) 
     &                    + 0.25*dsqrt(grad_modulus)

            grad_modulus = 0.d0
            do q = 1, depth
c loop over components of nabla q_m
               do d = 1, NDIM
                  grad_modulus = grad_modulus + 
     &               grad_x(i+1,j,q,d)*grad_x(i+1,j,q,d)
               enddo
            enddo
            grad_mod(i,j) = grad_mod(i,j) 
     &                    + 0.25*dsqrt(grad_modulus)

            grad_modulus = 0.d0
            do q = 1, depth
c loop over components of nabla q_m
               do d = 1, NDIM
                  grad_modulus = grad_modulus + 
     &               grad_y(i,j,q,d)*grad_y(i,j,q,d)
               enddo
            enddo
            grad_mod(i,j) = grad_mod(i,j) 
     &                    + 0.25*dsqrt(grad_modulus)

            grad_modulus = 0.d0
            do q = 1, depth
c loop over components of nabla q_m
               do d = 1, NDIM
                  grad_modulus = grad_modulus + 
     &               grad_y(i,j+1,q,d)*grad_y(i,j+1,q,d)
               enddo
            enddo
            grad_mod(i,j) = grad_mod(i,j) 
     &                    + 0.25*dsqrt(grad_modulus)
            
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine quatgrad_modulus_from_sides_harmonic(
     &   lo0, hi0, lo1, hi1,
     &   depth,
     &   grad_x, grad_y, 
     &   glo0, ghi0, glo1, ghi1,
     &   grad_mod,
     &   mlo0, mhi0, mlo1, mhi1,
     &   floor_type,
     &   grad_floor )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1,
     &   depth,
     &   glo0, ghi0, glo1, ghi1,
     &   mlo0, mhi0, mlo1, mhi1

      double precision
     &   grad_x(glo0:ghi0+1,glo1:ghi1,  depth,NDIM),
     &   grad_y(glo0:ghi0,  glo1:ghi1+1,depth,NDIM),
     &   grad_mod(mlo0:mhi0,mlo1:mhi1)

      double precision grad_floor, beta, grad_modulus

      character*(*) floor_type

      integer i, j, q, d
      
      if( floor_type(1:1) .eq. 's' )then
         beta = grad_floor*grad_floor
      else
         beta = 0.d0
      endif

      do j = lo1, hi1
         do i = lo0, hi0

            grad_mod(i,j) = 0.0d0

c we average contributions from 4 sides
            
            grad_modulus = 0.d0
            do q = 1, depth
c loop over components of nabla q_m
               do d = 1, NDIM
                  grad_modulus = grad_modulus + 
     &               grad_x(i,j,q,d)*grad_x(i,j,q,d)
               enddo
            enddo
            grad_mod(i,j) = grad_mod(i,j) 
     &                    + 1./dsqrt(grad_modulus+beta)

            grad_modulus = 0.d0
            do q = 1, depth
c loop over components of nabla q_m
               do d = 1, NDIM
                  grad_modulus = grad_modulus + 
     &               grad_x(i+1,j,q,d)*grad_x(i+1,j,q,d)
               enddo
            enddo
            grad_mod(i,j) = grad_mod(i,j) 
     &                    + 1./dsqrt(grad_modulus+beta)

            grad_modulus = 0.d0
            do q = 1, depth
c loop over components of nabla q_m
               do d = 1, NDIM
                  grad_modulus = grad_modulus + 
     &               grad_y(i,j,q,d)*grad_y(i,j,q,d)
               enddo
            enddo
            grad_mod(i,j) = grad_mod(i,j) 
     &                    + 1./dsqrt(grad_modulus+beta)

            grad_modulus = 0.d0
            do q = 1, depth
c loop over components of nabla q_m
               do d = 1, NDIM
                  grad_modulus = grad_modulus + 
     &               grad_y(i,j+1,q,d)*grad_y(i,j+1,q,d)
               enddo
            enddo
            grad_mod(i,j) = grad_mod(i,j) 
     &                    + 1./dsqrt(grad_modulus+beta)

            grad_mod(i,j) = 4./grad_mod(i,j)
            
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine quatgrad_modulus_from_sides_compact(
     &   lo0, hi0, lo1, hi1,
     &   depth,
     &   grad_x, grad_y, nggq,
     &   grad_mod, ngm)

      implicit none

      integer
     &   lo0, hi0, lo1, hi1,
     &   depth, nggq, ngm

      double precision
     &   grad_x(SIDE2d0(lo,hi,nggq),depth,NDIM),
     &   grad_y(SIDE2d1(lo,hi,nggq),depth,NDIM),
     &   grad_mod(CELL2d(lo,hi,ngm))

      double precision grad_modulus

      integer i, j, q
      
      do j = lo1, hi1
         do i = lo0, hi0

            grad_mod(i,j) = 0.0d0

c we average contributions from 4 sides
            
            grad_modulus = 0.d0
            do q = 1, depth
               grad_modulus = grad_modulus + 
     &            grad_x(i,j,q,1)*grad_x(i,j,q,1)
            enddo
            do q = 1, depth
               grad_modulus = grad_modulus + 
     &            grad_x(i+1,j,q,1)*grad_x(i+1,j,q,1)
            enddo
            do q = 1, depth
               grad_modulus = grad_modulus + 
     &            grad_y(i,j,q,2)*grad_y(i,j,q,2)
            enddo
            do q = 1, depth
               grad_modulus = grad_modulus + 
     &            grad_y(i,j+1,q,2)*grad_y(i,j+1,q,2)
            enddo

            grad_mod(i,j) = dsqrt(0.5d0*grad_modulus)
            
         enddo
      enddo

      return
      end


c-----------------------------------------------------------------------

      subroutine quatgrad_modulus(
     &   lo0, hi0, lo1, hi1,
     &   depth,
     &   grad_x, grad_y, nggq,
     &   grad_mod, ngm)

      implicit none

      integer
     &   lo0, hi0, lo1, hi1,
     &   depth, nggq, ngm

      double precision
     &   grad_x(CELL2d(lo,hi,nggq),depth),
     &   grad_y(CELL2d(lo,hi,nggq),depth),
     &   grad_mod(CELL2d(lo,hi,ngm))

      integer i, j, m
      
      do j = lo1, hi1
         do i = lo0, hi0

            grad_mod(i,j) = 0.0d0

            do m = 1, depth
               grad_mod(i,j) = grad_mod(i,j) +
     &            grad_x(i,j,m) * grad_x(i,j,m) +
     &            grad_y(i,j,m) * grad_y(i,j,m)
            enddo
            
            grad_mod(i,j) = dsqrt( grad_mod(i,j) )

         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

