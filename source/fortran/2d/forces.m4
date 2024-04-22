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
      integer eval_per_cell
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
                     forces(offset)   = forces(offset) + fx
                     forces(offset+1) = forces(offset+1) + fy
                  enddo
               enddo
            endif
         enddo
      enddo

      return
      end
