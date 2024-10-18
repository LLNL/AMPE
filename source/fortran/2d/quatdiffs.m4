c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
define(NDIM,2)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl

C        v is CellData with at least 1 ghost layer
C        diff_x and diff_y are SideData with at least 1 ghost layer

      subroutine diffs(
     &   lo0, hi0, lo1, hi1,
     &   v, ng,
     &   diff_x, diff_y, ngd
     &   )
        
      implicit none

      integer
     &   lo0, hi0, lo1, hi1,
     &   ng, ngd

      double precision
     &   v(CELL2d(lo,hi,ng)),
     &   diff_x(SIDE2d0(lo,hi,ngd)),
     &   diff_y(SIDE2d1(lo,hi,ngd))

c        local variables:
      integer i, j

c        x component
      do j = lo1-1, hi1+1
         do i = lo0, hi0+1
            diff_x(i,j) = v(i,j) - v(i-1,j)
         enddo
      enddo

c        y component
      do j = lo1, hi1+1
         do i = lo0-1, hi0+1
            diff_y(i,j) = v(i,j) - v(i,j-1)
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

C        q is CellData with ghosts
C        diff_x and diff_y are SideData WITH GHOSTS

      subroutine quatdiffs(
     &   lo0, hi0, lo1, hi1,
     &   depth,
     &   q, ngq,
     &   diff_x, diff_y, ngdiff
     &   )
        
      implicit none

      integer
     &   lo0, hi0, lo1, hi1,
     &   depth, ngq, ngdiff

      double precision q(CELL2d(lo,hi,ngq),depth)
      double precision diff_x(SIDE2d0(lo,hi,ngdiff),depth)
      double precision diff_y(SIDE2d1(lo,hi,ngdiff),depth)

c        local variables:
      integer i, j, m

c        Loop over the quaternion components
      do m = 1, depth

c           x component
         do j = lo1-1, hi1+1
            do i = lo0, hi0+1
               diff_x(i,j,m) = q(i,j,m) - q(i-1,j,m)
            enddo
         enddo

c           y component
         do j = lo1, hi1+1
            do i = lo0-1, hi0+1
               diff_y(i,j,m) = q(i,j,m) - q(i,j-1,m)
            enddo
         enddo

      enddo

      return
      end

c-----------------------------------------------------------------------

C        q is CellData with ghosts
C        diff_x and diff_y are SideData WITH GHOSTS

      subroutine quatdiffs_symm(
     &   lo0, hi0, lo1, hi1,
     &   depth,
     &   q, ngq,
     &   diff_x, diff_y, ngdiff,
     &   iqrot_x, iqrot_y, ngiq
     &   )
        
      implicit none

      integer
     &   lo0, hi0, lo1, hi1,
     &   depth, ngq, ngdiff, ngiq

      double precision q(CELL2d(lo,hi,ngq),depth)
      double precision diff_x(SIDE2d0(lo,hi,ngdiff),depth)
      double precision diff_y(SIDE2d1(lo,hi,ngdiff),depth)
      integer iqrot_x(SIDE2d0(lo,hi,ngiq))
      integer iqrot_y(SIDE2d1(lo,hi,ngiq))

c        local variables:
      integer i, j, m
      double precision q2_prime(depth), q2(depth)

c        X component
      do j = lo1-1, hi1+1
         do i = lo0, hi0+1

            do m = 1, depth
               q2(m) = q(i-1,j,m)
            enddo

            call quatsymmrotate( q2, iqrot_x(i,j), q2_prime, depth )

            do m = 1, depth
               diff_x(i,j,m) = q(i,j,m) - q2_prime(m)
            enddo

         enddo
      enddo

c        Y component
      do j = lo1, hi1+1
         do i = lo0-1, hi0+1

            do m = 1, depth
               q2(m) = q(i,j-1,m)
            enddo

            call quatsymmrotate( q2, iqrot_y(i,j), q2_prime, depth )

            do m = 1, depth
               diff_y(i,j,m) =  q(i,j,m) - q2_prime(m)
            enddo

         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

C        q is CellData with ghosts
C        diff_x and diff_y are SideData WITH GHOSTS
C
C Test code to see if numerical error in rotation matter.
C jlf July 2011
C
      subroutine quatdiffs_symm_new(
     &   lo0, hi0, lo1, hi1,
     &   depth,
     &   q, ngq,
     &   diff_x, diff_y, ngdiff,
     &   iqrot_x, iqrot_y, ngiq
     &   )
        
      implicit none

      integer
     &   lo0, hi0, lo1, hi1,
     &   depth, ngq, ngdiff, ngiq

      double precision q(CELL2d(lo,hi,ngq),depth)
      double precision diff_x(SIDE2d0(lo,hi,ngdiff),depth)
      double precision diff_y(SIDE2d1(lo,hi,ngdiff),depth)
      integer iqrot_x(SIDE2d0(lo,hi,ngiq))
      integer iqrot_y(SIDE2d1(lo,hi,ngiq))

c        local variables:
      integer i, j, m, iq
      double precision q2_prime(depth), q2(depth)
      double precision diffx_prime(depth), diffy_prime(depth)
      double precision diffx(depth), diffy(depth)

c        X component
      do j = lo1-1, hi1+1
         do i = lo0, hi0+1

            do m = 1, depth
               q2(m) = q(i,j,m)
            enddo

            iq = -1*iqrot_x(i,j)
            call quatsymmrotate( q2, iq, q2_prime, depth)
            do m = 1, depth
               diffx(m) = q2_prime(m) - q(i-1,j,m)
            enddo
            iq = iqrot_x(i,j)
            call quatsymmrotate( diffx,iq,diffx_prime,depth)

            do m = 1, depth
               q2(m) = q(i-1,j,m)
            enddo

            iq = iqrot_x(i,j)
            call quatsymmrotate( q2, iq, q2_prime, depth )

            do m = 1, depth
               diff_x(i,j,m) = 0.5d0*( q(i,j,m) - q2_prime(m) 
     &           + diffx_prime(m) )
            enddo

         enddo
      enddo

c        Y component
      do j = lo1, hi1+1
         do i = lo0-1, hi0+1

            do m = 1, depth
               q2(m) = q(i,j,m)
            enddo

            iq = -1*iqrot_y(i,j)
            call quatsymmrotate( q2, iq, q2_prime, depth )
            do m = 1, depth
               diffy(m) = q2_prime(m) - q(i,j-1,m)
            enddo
            iq = iqrot_y(i,j)
            call quatsymmrotate( diffy,iq,diffy_prime,depth)

            do m = 1, depth
               q2(m) = q(i,j-1,m)
            enddo

            iq = iqrot_y(i,j)
            call quatsymmrotate( q2, iq, q2_prime, depth )

            do m = 1, depth
               diff_y(i,j,m) =  0.5d0*( q(i,j,m) - q2_prime(m)
     &           + diffy_prime(m) )
            enddo

         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine quatavg(
     &   lo0, hi0, lo1, hi1,
     &   depth,
     &   q,
     &   qlo0, qhi0, qlo1, qhi1,
     &   avg, 
     &   alo0, ahi0, alo1, ahi1
     &   )
        
      implicit none

      integer
     &   lo0, hi0, lo1, hi1,
     &   depth,
     &   qlo0, qhi0, qlo1, qhi1,
     &   alo0, ahi0, alo1, ahi1

      double precision
     &   q(qlo0:qhi0,qlo1:qhi1,depth),
     &   avg(alo0:ahi0,alo1:ahi1,depth)

c        local variables:
      integer i, j, m

c        Loop over the quaternion components
      do m = 1, depth

         do j = lo1, hi1
            do i = lo0, hi0
               avg(i,j,m) = 0.5d0*q(i,j,m) 
     &                    + 0.125d0*( q(i-1,j,m)+q(i+1,j,m)
     &                              + q(i,j-1,m)+q(i,j+1,m) )
            enddo
         enddo

      enddo

      return
      end

