c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
define(NDIM,2)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine grad_cell(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   diff_x, diff_y, ngdiff,
     &   h,
     &   grad_x, grad_y, nggrad
     &   )

      implicit none

      integer
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   ngdiff, nggrad

      double precision
     &   diff_x(SIDE2d0(ifirst,ilast,ngdiff)),
     &   diff_y(SIDE2d1(ifirst,ilast,ngdiff)),
     &   h(2),
     &   grad_x(CELL2d(ifirst,ilast,nggrad)),
     &   grad_y(CELL2d(ifirst,ilast,nggrad))

      integer i, j
      double precision p5dxinv, p5dyinv

      p5dxinv = 0.5d0 / h(1)
      p5dyinv = 0.5d0 / h(2)

      do j = ifirst1, ilast1
         do i = ifirst0, ilast0
            grad_x(i,j) =
     &         ( diff_x(i+1,j) + diff_x(i,j) ) * p5dxinv
         enddo
      enddo
      
      do j = ifirst1, ilast1
         do i = ifirst0, ilast0
            grad_y(i,j) =
     &         ( diff_y(i,j+1) + diff_y(i,j) ) * p5dyinv
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine grad_side(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   diff_x, diff_y, ngdiff,
     &   h,
     &   grad_x_xside, grad_y_xside,
     &   grad_x_yside, grad_y_yside,
     &   nggrad
     &   )

      implicit none

      integer
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   ngdiff, nggrad

      double precision
     &   diff_x(SIDE2d0(ifirst,ilast,ngdiff)),
     &   diff_y(SIDE2d1(ifirst,ilast,ngdiff)),
     &   h(2),
     &   grad_x_xside(SIDE2d0(ifirst,ilast,nggrad)),
     &   grad_y_xside(SIDE2d0(ifirst,ilast,nggrad)),
     &   grad_x_yside(SIDE2d1(ifirst,ilast,nggrad)),
     &   grad_y_yside(SIDE2d1(ifirst,ilast,nggrad))

c        local variables:
      integer i, j
      double precision dxinv, dyinv
      double precision p25_dxinv, p25_dyinv

      dxinv = 1.d0 / h(1)
      dyinv = 1.d0 / h(2)
      p25_dxinv = 0.25d0 * dxinv
      p25_dyinv = 0.25d0 * dyinv

c        Compute gradients on the "X" side of the cell
      do j = ifirst1, ilast1
         do i = ifirst0, ilast0+1

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
      do j = ifirst1, ilast1+1
         do i = ifirst0, ilast0

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
c see Shukla and Giri, J. Comput. Phys. 276 (2014), p.259
c
      subroutine grad_side_isotropic(
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
      double precision p25_dxinv, p25_dyinv, p12_dxinv, p12_dyinv

      dxinv = 1.d0 / h(1)
      dyinv = 1.d0 / h(2)
      p25_dxinv = 0.25d0 * dxinv
      p25_dyinv = 0.25d0 * dyinv
      p12_dxinv = dxinv/12.d0
      p12_dyinv = dyinv/12.d0

c        Compute gradients on the "X" side of the cell
      do j = lo1, hi1
         do i = lo0, hi0+1

c              X component of gradient at "X" side
            grad_x_xside(i,j) = p12_dxinv * (
     &          diff_x(i,j-1) +
     &          diff_x(i,j)*10.d0 +
     &          diff_x(i,j+1) )

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
            grad_y_yside(i,j) = p12_dyinv * (
     &          diff_y(i,j-1) +
     &          diff_y(i,j)*10.d0 +
     &          diff_y(i,j+1) )

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

      logical isnan

      threshold2 = threshold*threshold
      invthreshold = 1./threshold
      
      do j = lo1, hi1
         do i = lo0, hi0

            grad2 =
     &            grad_x(i,j) * grad_x(i,j) +
     &            grad_y(i,j) * grad_y(i,j)
c            if( isnan(grad2) )stop '"grad2" is a NaN'
            
c            if( grad2>threshold2 )then
             vel(i,j) = phi_dot(i,j)/dsqrt(grad2+threshold2)
c               vel(i,j) = phi_dot(i,j)
c            else
c               vel(i,j) = phi_dot(i,j)*invthreshold
c               if( abs(vel(i,j)) .gt. 10. )then
c                  print*,'vel=',vel(i,j),', phi_dot(i,j)=',phi_dot(i,j)
c               endif
c               vel(i,j) = phi_dot(i,j)
c            endif
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

