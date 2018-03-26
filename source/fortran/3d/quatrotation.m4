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
c
      subroutine quat_symm_rotation(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   q, ngq,
     &   depth,
     &   rot_x, rot_y, rot_z, ngrot
     &   )
        
      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth, ngq, ngrot

      double precision
     &   q(CELL3d(lo,hi,ngq),depth)

      integer
     &   rot_x(SIDE3d0(lo,hi,ngrot)),
     &   rot_y(SIDE3d1(lo,hi,ngrot)),
     &   rot_z(SIDE3d2(lo,hi,ngrot))

c        local variables:
      integer i, j, k, m, iq
      double precision q2_prime(depth), q1(depth), q2(depth)

c        X component
      do k = lo2-1, hi2+1
         do j = lo1-1, hi1+1
            do i = lo0, hi0+1

               do m = 1, depth
                  q1(m) = q(i,j,k,m)
                  q2(m) = q(i-1,j,k,m)
               enddo

               call quatfindsymm( 
     &            q1, q2, rot_x(i,j,k), q2_prime, depth )

            enddo
         enddo
      enddo

c        Y component
      do k = lo2-1, hi2+1
         do j = lo1, hi1+1
            do i = lo0-1, hi0+1

               do m = 1, depth
                  q1(m) = q(i,j,k,m)
                  q2(m) = q(i,j-1,k,m)
               enddo

               call quatfindsymm(
     &            q1, q2, rot_y(i,j,k), q2_prime, depth )

            enddo
         enddo
      enddo

c        Z component
      do k = lo2, hi2+1
         do j = lo1-1, hi1+1
            do i = lo0-1, hi0+1

               do m = 1, depth
                  q1(m) = q(i,j,k,m)
                  q2(m) = q(i,j,k-1,m)
               enddo

               call quatfindsymm(
     &            q1, q2, rot_z(i,j,k), q2_prime, depth )

            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------

      subroutine quat_fundamental(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   quat,
     &   qlo0, qhi0, qlo1, qhi1, qlo2, qhi2,
     &   depth
     &   )
        
      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   qlo0, qhi0, qlo1, qhi1, qlo2, qhi2,
     &   depth

      double precision
     &   quat(qlo0:qhi0,qlo1:qhi1,qlo2:qhi2,depth)

      integer i, j, k, m, iq
      double precision q1(depth)
      double precision q2(depth)
      double precision q2_prime(depth)

      double precision ZERO, ONE, PI, PI_OVER4
      parameter( ZERO=0.0d0 )
      parameter( ONE=1.0d0 )

      if ( depth .eq. 2 ) then
         call quatset2( q1, ONE, ZERO )
      else if ( depth .eq. 4 ) then
         call quatset( q1, ONE, ZERO, ZERO, ZERO )
      else 
         q1(1) = ZERO
      endif

      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0

               do m = 1, depth
                  q2(m) = quat(i,j,k,m)
               enddo

               iq = 1
               call quatfindsymm( q1, q2, iq, q2_prime, depth )

               do m = 1, depth
                  quat(i,j,k,m) = q2_prime(m)
               enddo

            enddo
         enddo
      enddo

      return
      end
