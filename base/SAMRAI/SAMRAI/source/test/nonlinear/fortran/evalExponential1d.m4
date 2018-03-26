define(NDIM,1)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim1d.i)dnl

      subroutine evalExponential1d(
     & lo0, hi0, 
     & u,
     & lambda,
     & lexpu )

c  Evaluate exponential term in modified Bratu problem.

      implicit none

      integer lo0

      integer hi0

      double precision lambda

      double precision u(CELL1d(lo,hi,0))
      double precision lexpu(CELL1d(lo,hi,0))

      integer i

      do i = lo0, hi0
         lexpu(i) = lambda*exp(u(i))
      end do

      return
      end
