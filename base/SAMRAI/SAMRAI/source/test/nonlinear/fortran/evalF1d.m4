define(NDIM,1)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim1d.i)dnl

      subroutine evalF1d(
     & lo0, hi0,
     & dx,
     & unew, ucur, 
     & f )

c Evaluate rhs of modified bratu problem at u by subtracting off the
c terms that are due to BE time discretization.

      implicit none

      integer lo0
      integer hi0

      double precision dx(0:NDIM-1)

      double precision f(CELL1d(lo,hi,0))
      double precision unew(CELL1d(lo,hi,0))
      double precision ucur(CELL1d(lo,hi,0))
     
      integer i

      double precision vol

      vol = dx(0)
      do i = lo0, hi0
         f(i) = f(i) - vol*(unew(i) - ucur(i))
      end do

      return
      end
