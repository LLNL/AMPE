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
c
      subroutine quat_symm_rotation(
     &   lo0, hi0, lo1, hi1,
     &   q, ngq,
     &   depth,
     &   rot_x, rot_y, ngrot
     &   )
        
      implicit none

      integer
     &   lo0, hi0, lo1, hi1,
     &   depth, ngq, ngrot

      double precision
     &   q(CELL2d(lo,hi,ngq),depth)

      integer
     &   rot_x(SIDE2d0(lo,hi,ngrot)),
     &   rot_y(SIDE2d1(lo,hi,ngrot))

c        local variables:
      integer i, j, m
      double precision q2_prime(depth), q1(depth), q2(depth)

c        X component
      do j = lo1-1, hi1+1
         do i = lo0, hi0+1

            do m = 1, depth
               q1(m) = q(i,j,m)
               q2(m) = q(i-1,j,m)
            enddo

            call quatfindsymm( q1, q2, rot_x(i,j), q2_prime, depth )

         enddo
      enddo

c        Y component
      do j = lo1, hi1+1
         do i = lo0-1, hi0+1

            do m = 1, depth
               q1(m) = q(i,j,m)
               q2(m) = q(i,j-1,m)
            enddo

            call quatfindsymm( q1, q2, rot_y(i,j), q2_prime, depth )

         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------

      subroutine quat_fundamental(
     &   lo0, hi0, lo1, hi1,
     &   quat,
     &   qlo0, qhi0, qlo1, qhi1,
     &   depth
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1,
     &   qlo0, qhi0, qlo1, qhi1,
     &   depth

      double precision
     &   quat(qlo0:qhi0,qlo1:qhi1,depth)

      integer i, j, m, iq
      double precision q1(depth)
      double precision q2(depth)
      double precision q2_prime(depth)

      double precision ZERO, ONE
      parameter( ZERO=0.0d0 )
      parameter( ONE=1.0d0 )

      if ( depth .eq. 2 ) then
         call quatset2( q1, ONE, ZERO )
      else if ( depth .eq. 4 ) then
         call quatset( q1, ONE, ZERO, ZERO, ZERO )
      else
         q1(1) = ZERO
      endif

      do j = lo1, hi1
         do i = lo0, hi0

            do m = 1, depth
               q2(m) = quat(i,j,m)
            enddo

            iq = 1
            call quatfindsymm( q1, q2, iq, q2_prime, depth )

            do m = 1, depth
               quat(i,j,m) = q2_prime(m)
            enddo

         enddo
      enddo

      return
      end
