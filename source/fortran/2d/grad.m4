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
define(NDIM,2)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine grad_cell(
     &   lo0, hi0, lo1, hi1,
     &   diff_x, diff_y,
     &   dlo0, dhi0, dlo1, dhi1,
     &   h,
     &   grad_x, grad_y, 
     &   glo0, ghi0, glo1, ghi1
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1,
     &   dlo0, dhi0, dlo1, dhi1,
     &   glo0, ghi0, glo1, ghi1

      double precision
     &   diff_x(dlo0:dhi0+1,dlo1:dhi1),
     &   diff_y(dlo0:dhi0,dlo1:dhi1+1),
     &   h(2),
     &   grad_x(glo0:ghi0,glo1:ghi1),
     &   grad_y(glo0:ghi0,glo1:ghi1)

      integer i, j
      double precision p5dxinv, p5dyinv

      p5dxinv = 0.5d0 / h(1)
      p5dyinv = 0.5d0 / h(2)

      do j = lo1, hi1
         do i = lo0, hi0
            grad_x(i,j) =
     &         ( diff_x(i+1,j) + diff_x(i,j) ) * p5dxinv
         enddo
      enddo
      
      do j = lo1, hi1
         do i = lo0, hi0
            grad_y(i,j) =
     &         ( diff_y(i,j+1) + diff_y(i,j) ) * p5dyinv
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine grad_side(
     &   lo0, hi0, lo1, hi1,
     &   diff_x, diff_y,
     &   dlo0, dhi0, dlo1, dhi1,
     &   h,
     &   grad_x_xside, grad_y_xside,
     &   grad_x_yside, grad_y_yside,
     &   glo0, ghi0, glo1, ghi1
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1,
     &   dlo0, dhi0, dlo1, dhi1,
     &   glo0, ghi0, glo1, ghi1

      double precision
     &   diff_x(dlo0:dhi0+1,dlo1:dhi1),
     &   diff_y(dlo0:dhi0,dlo1:dhi1+1),
     &   h(2),
     &   grad_x_xside(glo0:ghi0+1,glo1:ghi1),
     &   grad_y_xside(glo0:ghi0+1,glo1:ghi1),
     &   grad_x_yside(glo0:ghi0,glo1:ghi1+1),
     &   grad_y_yside(glo0:ghi0,glo1:ghi1+1)

c        local variables:
      integer i, j
      double precision dxinv, dyinv
      double precision p25_dxinv, p25_dyinv

      dxinv = 1.d0 / h(1)
      dyinv = 1.d0 / h(2)
      p25_dxinv = 0.25d0 * dxinv
      p25_dyinv = 0.25d0 * dyinv

c        Compute gradients on the "X" side of the cell
      do j = lo1, hi1
         do i = lo0, hi0+1

c              X component of gradient at "X" side
            grad_x_xside(i,j) = dxinv * diff_x(i,j)

c              Y component of gradient at "X" side
            grad_y_xside(i,j) = p25_dyinv * (
     &         diff_y(i-1,j+1) + diff_y(i-1,j) +
     &         diff_y(i,  j+1) + diff_y(i,  j)
     &         ) 

         enddo
      enddo

c        Compute gradients on the "Y" side of the cell
      do j = lo1, hi1+1
         do i = lo0, hi0

c              X component of gradient at "Y" side
            grad_x_yside(i,j) = p25_dxinv * (
     &         diff_x(i+1,j-1) + diff_x(i,j-1) +
     &         diff_x(i+1,j  ) + diff_x(i,j  )
     &         ) 

c              Y component of gradient at "Y" side
            grad_y_yside(i,j) = dyinv * diff_y(i,j)

         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine velocity(
     &   lo0, hi0, lo1, hi1,
     &   threshold,
     &   grad_x, grad_y, 
     &   phi_dot, 
     &   vel )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1

      double precision
     &   grad_x(lo0:hi0,lo1:hi1),
     &   grad_y(lo0:hi0,lo1:hi1),
     &   phi_dot(lo0:hi0,lo1:hi1),
     &   vel(lo0:hi0,lo1:hi1)

      integer i, j
      double precision threshold, threshold2
      double precision invthreshold, grad2

      threshold2 = threshold*threshold
      invthreshold = 1./threshold
      
      do j = lo1, hi1
         do i = lo0, hi0

            grad2 =
     &            grad_x(i,j) * grad_x(i,j) +
     &            grad_y(i,j) * grad_y(i,j)
            if( isnan(grad2) )stop '"grad2" is a NaN'
            
            if( grad2>threshold2 )then
               vel(i,j) = phi_dot(i,j)/dsqrt(grad2)
c               vel(i,j) = phi_dot(i,j)
            else
               vel(i,j) = phi_dot(i,j)*invthreshold
c               if( abs(vel(i,j)) .gt. 10. )then
c                  print*,'vel=',vel(i,j),', phi_dot(i,j)=',phi_dot(i,j)
c               endif
c               vel(i,j) = phi_dot(i,j)
            endif
            if( isnan(vel(i,j)) )then
               print*,'phi_dot(i,j)=',phi_dot(i,j)
               print*,'vel=',vel(i,j)
               print*,'grad2=',grad2
               stop '"vel(i,j)" is a NaN'
            endif
         enddo
      enddo

      return
      end

