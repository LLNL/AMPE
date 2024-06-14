c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
c 
c
c Based on Moelans 2011 multi-order parameters
c
c Output:
c   Integral of forces over patch added to forces
c
define(NDIM,2)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine phiphi_forces(
     &   lo0, hi0, lo1, hi1,
     &   dx,
     &   phi, pghosts, nphases,
     &   dgamma, m,
     &   weight, forces
     &   )

      implicit none

      integer lo0, hi0, lo1, hi1
      integer pghosts, nphases
      double precision dgamma, m
      double precision phi(CELL2d(lo,hi,pghosts), nphases)
      double precision weight(CELL2d(lo,hi,0))
      double precision dx(0:NDIM-1)
      double precision forces(NDIM*nphases*nphases)

      integer i, j, p1, p2, offset
      double precision dxinv, dyinv
      double precision diffx, diffy
      double precision f, fx, fy, factor
c
      dxinv = 1.d0 / dx(0)
      dyinv = 1.d0 / dx(1)
c
      factor = 0.5d0 * dgamma * m
c
      do p1 = 1, nphases
         do p2 = 1, nphases
            if( p1 .ne. p2 )then
              offset = NDIM*nphases*(p1-1)+p2*2-1
c integral over domain
               do j = lo1, hi1
                  do i = lo0, hi0
                     diffx = dxinv * (
     &                  phi(i+1,j,p2) - phi(i-1,j,p2))
                     diffy = dyinv * (
     &                  phi(i,j+1,p2) - phi(i,j-1,p2))

                     f = phi(i,j,p1)*phi(i,j,p1)
     &                 * phi(i,j,p2)* factor
                     fx = f*diffx
                     fy = f*diffy
c add force density to integral over domain
                     forces(offset)   = forces(offset)
     &                                + fx*weight(i,j)
                     forces(offset+1) = forces(offset+1)
     &                                + fy*weight(i,j)
                  enddo
               enddo
            endif
         enddo
      enddo

      return
      end
c
c
c
      subroutine wang_forces(
     &   lo0, hi0, lo1, hi1,
     &   dx,
     &   phi, pghosts, nphases,
     &   conc, ngc, nc,
     &   rhoe, cthreshold,
     &   weight, forces
     &   )

      implicit none

      integer lo0, hi0, lo1, hi1
      integer pghosts, nphases, ngc, nc
      double precision rhoe, cthreshold
      double precision phi(CELL2d(lo,hi,pghosts), nphases)
      double precision weight(CELL2d(lo,hi,0))
      double precision conc(CELL2d(lo,hi,ngc), nc)
      double precision dx(0:NDIM-1)
      double precision forces(NDIM*nphases*nphases)

      integer i, j, p1, p2, offset, ic
      double precision dxinv, dyinv
      double precision diffx1, diffy1, diffx2, diffy2
      double precision fx, fy
      double precision pp, cdiff, cutoff_slope, cfactor
c
      double precision interp_func
c
      cutoff_slope = 10.d0
c
      dxinv = 0.5d0 / dx(0)
      dyinv = 0.5d0 / dx(1)
c
      do p1 = 1, nphases
         do p2 = 1, nphases
            if( p1 .ne. p2 )then
              offset = NDIM*nphases*(p1-1)+p2*2-1
c integral over domain
               do j = lo1, hi1
                  do i = lo0, hi0
                     diffx1 = dxinv * (
     &                  phi(i+1,j,p1) - phi(i-1,j,p1))
                     diffy1 = dyinv * (
     &                  phi(i,j+1,p1) - phi(i,j-1,p1))
                     diffx2 = dxinv * (
     &                  phi(i+1,j,p2) - phi(i-1,j,p2))
                     diffy2 = dyinv * (
     &                  phi(i,j+1,p2) - phi(i,j-1,p2))
c smooth cutoff
                     pp = 0.5d0+(phi(i,j,p1)*phi(i,j,p2)-cthreshold)
     &                         *cutoff_slope
                     cfactor = interp_func( pp, 'p')
                     cdiff = 0.d0
                     do ic = 1, nc
                        cdiff = cdiff + conc(i,j,ic)
                     enddo
                     cdiff = cdiff-rhoe
                     fx = cdiff*cfactor*(diffx1-diffx2)
                     fy = cdiff*cfactor*(diffy1-diffx2)
c add force density to integral over domain
                     forces(offset)   = forces(offset)
     &                                + fx*weight(i,j)
                     forces(offset+1) = forces(offset+1)
     &                                + fy*weight(i,j)
                  enddo
               enddo
            endif
         enddo
      enddo

      return
      end
