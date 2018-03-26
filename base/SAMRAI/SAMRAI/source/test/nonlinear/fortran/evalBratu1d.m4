define(NDIM,1)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim1d.i)dnl

      subroutine evalbratu1d(
     & lo0, hi0, ghostcells,
     & gew, 
     & f, lexpu,
     & v,
     & u,
     & dx, dt,
     & r )

c  Evaluate modified bratu problem at u by assembling fluxes
c  (in gew), sources (in f and lexpu), and previous
c  time step (in v).

      implicit none

      integer lo0

      integer hi0

      integer ghostcells

      double precision dt

      double precision dx(0:NDIM-1)

      double precision gew(SIDE1d(lo,hi,0))
      double precision f(CELL1d(lo,hi,0))
      double precision lexpu(CELL1d(lo,hi,0))
      double precision v(CELL1d(lo,hi,0))
      double precision u(CELL1d(lo,hi,ghostcells))
      double precision r(CELL1d(lo,hi,0))
     
      integer i

      do i = lo0, hi0
         r(i) = u(i) - v(i)
     &            - dt*((gew(i+1) - gew(i))/dx(0) +
     &                  lexpu(i) + f(i))
      end do

      return
      end
