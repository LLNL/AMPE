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

      subroutine grad_cell(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   diff_x, diff_y, diff_z,
     &   dlo0, dhi0, dlo1, dhi1, dlo2, dhi2,
     &   h,
     &   grad_x, grad_y, grad_z,
     &   glo0, ghi0, glo1, ghi1, glo2, ghi2
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   dlo0, dhi0, dlo1, dhi1, dlo2, dhi2,
     &   glo0, ghi0, glo1, ghi1, glo2, ghi2

      double precision
     &   diff_x(dlo0:dhi0+1,dlo1:dhi1,dlo2:dhi2),
     &   diff_y(dlo0:dhi0,dlo1:dhi1+1,dlo2:dhi2),
     &   diff_z(dlo0:dhi0,dlo1:dhi1,dlo2:dhi2+1),
     &   h(3),
     &   grad_x(glo0:ghi0,glo1:ghi1,glo2:ghi2),
     &   grad_y(glo0:ghi0,glo1:ghi1,glo2:ghi2),
     &   grad_z(glo0:ghi0,glo1:ghi1,glo2:ghi2)

      integer i, j, k
      double precision p5dxinv, p5dyinv, p5dzinv

      p5dxinv = 0.5d0 / h(1)
      p5dyinv = 0.5d0 / h(2)
      p5dzinv = 0.5d0 / h(3)

      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0
               grad_x(i,j,k) =
     &            ( diff_x(i+1,j,k) + diff_x(i,j,k) ) * p5dxinv
            enddo
         enddo
      enddo
      
      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0
               grad_y(i,j,k) =
     &            ( diff_y(i,j+1,k) + diff_y(i,j,k) ) * p5dyinv
            enddo
         enddo
      enddo

      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0
               grad_z(i,j,k) =
     &            ( diff_z(i,j,k+1) + diff_z(i,j,k) ) * p5dzinv
            enddo
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine grad_side(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   diff_x, diff_y, diff_z,
     &   dlo0, dhi0, dlo1, dhi1, dlo2, dhi2,
     &   h,
     &   grad_x_xside, grad_y_xside, grad_z_xside,
     &   grad_x_yside, grad_y_yside, grad_z_yside,
     &   grad_x_zside, grad_y_zside, grad_z_zside,
     &   glo0, ghi0, glo1, ghi1, glo2, ghi2
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   dlo0, dhi0, dlo1, dhi1, dlo2, dhi2,
     &   glo0, ghi0, glo1, ghi1, glo2, ghi2

      double precision
     &   diff_x(dlo0:dhi0+1,dlo1:dhi1,dlo2:dhi2),
     &   diff_y(dlo0:dhi0,dlo1:dhi1+1,dlo2:dhi2),
     &   diff_z(dlo0:dhi0,dlo1:dhi1,dlo2:dhi2+1),
     &   h(3),
     &   grad_x_xside(glo0:ghi0+1,glo1:ghi1,glo2:ghi2),
     &   grad_y_xside(glo0:ghi0+1,glo1:ghi1,glo2:ghi2),
     &   grad_z_xside(glo0:ghi0+1,glo1:ghi1,glo2:ghi2),
     &   grad_x_yside(glo0:ghi0,glo1:ghi1+1,glo2:ghi2),
     &   grad_y_yside(glo0:ghi0,glo1:ghi1+1,glo2:ghi2),
     &   grad_z_yside(glo0:ghi0,glo1:ghi1+1,glo2:ghi2),
     &   grad_x_zside(glo0:ghi0,glo1:ghi1,glo2:ghi2+1),
     &   grad_y_zside(glo0:ghi0,glo1:ghi1,glo2:ghi2+1),
     &   grad_z_zside(glo0:ghi0,glo1:ghi1,glo2:ghi2+1)

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
      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0+1

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
      do k = lo2, hi2
         do j = lo1, hi1+1
            do i = lo0, hi0

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
      do k = lo2, hi2+1
         do j = lo1, hi1
            do i = lo0, hi0
               
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
            
               if( grad2>threshold2 )then
                  vel(i,j,k) = phi_dot(i,j,k)/dsqrt(grad2)
c                  vel(i,j,k) = phi_dot(i,j,k)
               else
                  vel(i,j,k) = phi_dot(i,j,k)*invthreshold
c                  vel(i,j,k) = phi_dot(i,j,k)
               endif

            enddo
         enddo
      enddo

      return
      end
