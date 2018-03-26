define(NDIM,2)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine evalF2d(
     & lo0, hi0, lo1, hi1,
     & dx,
     & unew, ucur, 
     & f )

c Evaluate rhs of modified bratu problem at u by subtracting off the
c terms that are due to BE time discretization.

      implicit none

      integer lo0
      integer lo1

      integer hi0
      integer hi1

      double precision dx(0:NDIM-1)

      double precision f(CELL2d(lo,hi,0))
      double precision unew(CELL2d(lo,hi,0))
      double precision ucur(CELL2d(lo,hi,0))
     
      integer i
      integer j

      double precision vol

      vol = dx(0)*dx(1)

      do j = lo1, hi1
         do i = lo0, hi0
            f(i,j) = f(i,j) - vol*(unew(i,j) - ucur(i,j))
         end do
      end do

      return
      end
