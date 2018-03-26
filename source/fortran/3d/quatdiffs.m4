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

C        v is CellData with ghosts
C        diff_x, diff_y, and diff_z are SideData WITH GHOSTS

      subroutine diffs(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   v,
     &   vlo0, vhi0, vlo1, vhi1, vlo2, vhi2,
     &   diff_x, diff_y, diff_z,
     &   dlo0, dhi0, dlo1, dhi1, dlo2, dhi2
     &   )
        
      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   vlo0, vhi0, vlo1, vhi1, vlo2, vhi2,
     &   dlo0, dhi0, dlo1, dhi1, dlo2, dhi2

      double precision
     &   v(vlo0:vhi0,vlo1:vhi1,vlo2:vhi2),
     &   diff_x(dlo0:dhi0+1,dlo1:dhi1,dlo2:dhi2),
     &   diff_y(dlo0:dhi0,dlo1:dhi1+1,dlo2:dhi2),
     &   diff_z(dlo0:dhi0,dlo1:dhi1,dlo2:dhi2+1)

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
      integer i, j, k, m, iq
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
