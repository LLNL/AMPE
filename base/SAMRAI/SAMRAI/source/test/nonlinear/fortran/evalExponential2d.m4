define(NDIM,2)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine evalExponential2d(
     & lo0, hi0, lo1, hi1, 
     & u,
     & lambda,
     & lexpu )

c  Evaluate exponential term in modified Bratu problem.

      implicit none

      integer lo0
      integer lo1

      integer hi0
      integer hi1

      double precision lambda

      double precision u(CELL2d(lo,hi,0))
      double precision lexpu(CELL2d(lo,hi,0))

      integer i
      integer j

      do j = lo1, hi1
         do i = lo0, hi0
            lexpu(i,j) = lambda*exp(u(i,j))
         end do
      end do

      return
      end
