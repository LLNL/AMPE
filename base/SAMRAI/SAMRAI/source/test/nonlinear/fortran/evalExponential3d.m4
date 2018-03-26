define(NDIM,3)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine evalExponential3d(
     & lo0, hi0, lo1, hi1, lo2, hi2,
     & u,
     & lambda,
     & lexpu )

c  Evaluate exponential term in modified Bratu problem.

      implicit none

      integer lo0
      integer lo1
      integer lo2

      integer hi0
      integer hi1
      integer hi2

      double precision lambda

      double precision u(CELL3d(lo,hi,0))
      double precision lexpu(CELL3d(lo,hi,0))

      integer i
      integer j
      integer k

      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0
               lexpu(i,j,k) = lambda*exp(u(i,j,k))
            end do
         end do
      end do

      return
      end
