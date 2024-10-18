c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
define(NDIM,3)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl

C v is CellData with at least 1 ghost layer
C diff_x, diff_y, and diff_z are SideData with at least 1 ghost layer

      subroutine diffs(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   v, ng,
     &   diff_x, diff_y, diff_z, ngd
     &   )
        
      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   ng, ngd

      double precision
     &   v(CELL3d(lo,hi,ng)),
     &   diff_x(SIDE3d0(lo,hi,ngd)),
     &   diff_y(SIDE3d1(lo,hi,ngd)),
     &   diff_z(SIDE3d2(lo,hi,ngd))

c        local variables:
      integer i, j, k

c        x component
      do k = lo2-1, hi2+1
         do j = lo1-1, hi1+1
            do i = lo0, hi0+1
               diff_x(i,j,k) = v(i,j,k) - v(i-1,j,k)
            enddo
         enddo
      enddo

c        y component
      do k = lo2-1, hi2+1
         do j = lo1, hi1+1
            do i = lo0-1, hi0+1
               diff_y(i,j,k) = v(i,j,k) - v(i,j-1,k)
            enddo
         enddo
      enddo

c        z component
      do k = lo2, hi2+1
         do j = lo1-1, hi1+1
            do i = lo0-1, hi0+1
               diff_z(i,j,k) = v(i,j,k) - v(i,j,k-1)
            enddo
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

C        q is CellData with ghosts
C        diff_x, diff_y, and diff_z are SideData WITH ghosts

      subroutine quatdiffs(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth,
     &   q, ngq,
     &   diff_x, diff_y, diff_z, ngdiff
     &   )
        
      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth, ngq, ngdiff

      double precision q(CELL3d(lo,hi,ngq),depth)
      double precision diff_x(SIDE3d0(lo,hi,ngdiff),depth)
      double precision diff_y(SIDE3d1(lo,hi,ngdiff),depth)
      double precision diff_z(SIDE3d2(lo,hi,ngdiff),depth)

c        local variables:
      integer i, j, k, m

c        Loop over the quaterion components
      do m = 1, depth

c           x component
         do k = lo2-1, hi2+1
            do j = lo1-1, hi1+1
               do i = lo0, hi0+1
                  diff_x(i,j,k,m) = q(i,j,k,m) - q(i-1,j,k,m)
               enddo
            enddo
         enddo

c           y component
         do k = lo2-1, hi2+1
            do j = lo1, hi1+1
               do i = lo0-1, hi0+1
                  diff_y(i,j,k,m) = q(i,j,k,m) - q(i,j-1,k,m)
               enddo
            enddo
         enddo

c           z component
         do k = lo2, hi2+1
            do j = lo1-1, hi1+1
               do i = lo0-1, hi0+1
                  diff_z(i,j,k,m) = q(i,j,k,m) - q(i,j,k-1,m)
               enddo
            enddo
         enddo

      enddo

      return
      end

c-----------------------------------------------------------------------

C        q is CellData with ghosts
C        diff_x, diff_y, and diff_z are SideData WITH ghosts

      subroutine quatdiffs_symm(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth,
     &   q, ngq,
     &   diff_x, diff_y, diff_z, ngdiff,
     &   iqrot_x, iqrot_y, iqrot_z, ngiq
     &   )
        
      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth, ngq, ngdiff, ngiq

      double precision q(CELL3d(lo,hi,ngq),depth)
      double precision diff_x(SIDE3d0(lo,hi,ngdiff),depth)
      double precision diff_y(SIDE3d1(lo,hi,ngdiff),depth)
      double precision diff_z(SIDE3d2(lo,hi,ngdiff),depth)
      integer iqrot_x(SIDE3d0(lo,hi,ngiq))
      integer iqrot_y(SIDE3d1(lo,hi,ngiq))
      integer iqrot_z(SIDE3d2(lo,hi,ngiq))

c        local variables:
      integer i, j, k, m
      double precision q2_prime(depth), q2(depth)

c        X component
      do k = lo2-1, hi2+1
         do j = lo1-1, hi1+1
            do i = lo0, hi0+1

               do m = 1, depth
                  q2(m) = q(i-1,j,k,m)
               enddo

               call quatsymmrotate(
     &            q2, iqrot_x(i,j,k), q2_prime, depth )

               do m = 1, depth
                  diff_x(i,j,k,m) = q(i,j,k,m) - q2_prime(m)
               enddo

            enddo
         enddo
      enddo

c        Y component
      do k = lo2-1, hi2+1
         do j = lo1, hi1+1
            do i = lo0-1, hi0+1

               do m = 1, depth
                  q2(m) = q(i,j-1,k,m)
               enddo

               call quatsymmrotate(
     &            q2, iqrot_y(i,j,k), q2_prime, depth )

               do m = 1, depth
                  diff_y(i,j,k,m) = q(i,j,k,m) - q2_prime(m)
               enddo

            enddo
         enddo
      enddo

c        Z component
      do k = lo2, hi2+1
         do j = lo1-1, hi1+1
            do i = lo0-1, hi0+1

               do m = 1, depth
                  q2(m) = q(i,j,k-1,m)
               enddo

               call quatsymmrotate( 
     &            q2, iqrot_z(i,j,k), q2_prime, depth )

               do m = 1, depth
                  diff_z(i,j,k,m) = q(i,j,k,m) - q2_prime(m)
               enddo

            enddo
         enddo
      enddo

      return
      end
