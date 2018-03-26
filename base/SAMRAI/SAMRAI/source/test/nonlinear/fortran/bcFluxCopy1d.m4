define(NDIM,1)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim1d.i)dnl

      subroutine ewfluxcopy1d(
     & lo0, hi0, 
     & ewflux,
     & sideflux,
     & side )

c  Copy face fluxes in the x-direction.

      implicit none

      integer lo0

      integer hi0

      integer side

      double precision ewflux(SIDE1d(lo,hi,0))
      double precision sideflux(OUTERSIDE1d(lo,hi,0))

      integer ihi
      integer ilo

      if (side .eq. 0) then
         ilo = lo0
         sideflux(1) = ewflux(ilo)
      else if (side .eq. 1) then
         ihi = hi0+1
         sideflux(1) = ewflux(ihi)
      end if

      return
      end
