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
define(NDIM,3)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine phiphi_forces(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   dx,
     &   phi, pghosts, nphases,
     &   dgamma, m,
     &   weight, forces
     &   )

      implicit none

      integer lo0, hi0, lo1, hi1, lo2, hi2
      integer pghosts, nphases
      integer eval_per_cell
      double precision dgamma, m
      double precision phi(CELL3d(lo,hi,pghosts), nphases)
      double precision weight(CELL3d(lo,hi,0))
      double precision dx(0:NDIM-1)
      double precision forces(NDIM*nphases*nphases)

      integer i, j, k, p1, p2, offset
      double precision dxinv, dyinv, dzinv
      double precision diffx, diffy, diffz
      double precision f, fx, fy, fz, factor
c
      dxinv = 1.d0 / dx(0)
      dyinv = 1.d0 / dx(1)
      dzinv = 1.d0 / dx(2)
c
      factor = 0.5d0 * dgamma * m
c
      do p1 = 1, nphases
         do p2 = 1, nphases
            if( p1 .ne. p2 )then
              offset = NDIM*nphases*(p1-1)+p2*2-1
c integral over domain
               do k = lo2, hi2
                  do j = lo1, hi1
                     do i = lo0, hi0
                        diffx = dxinv * (
     &                     phi(i+1,j,k,p2) - phi(i-1,j,k,p2))
                        diffy = dyinv * (
     &                     phi(i,j+1,k,p2) - phi(i,j-1,k,p2))
                        diffz = dzinv * (
     &                     phi(i,j,k+1,p2) - phi(i,j,k-1,p2))

                        f = phi(i,j,k,p1)*phi(i,j,k,p1)
     &                    * phi(i,j,k,p2)* factor
                        fx = f*diffx
                        fy = f*diffy
                        fz = f*diffz
                        forces(offset)   = forces(offset) + fx
                        forces(offset+1) = forces(offset+1) + fy
                        forces(offset+2) = forces(offset+2) + fz
                     enddo
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
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   dx,
     &   phi, pghosts, nphases,
     &   conc, ngc, nc,
     &   rhoe, cthreshold,
     &   weight, forces
     &   )

      implicit none

      integer lo0, hi0, lo1, hi1, lo2, hi2
      integer pghosts, nphases, ngc, nc
      double precision rhoe, cthreshold
      double precision phi(CELL3d(lo,hi,pghosts), nphases)
      double precision weight(CELL3d(lo,hi,0))
      double precision conc(CELL3d(lo,hi,ngc), nc)
      double precision dx(0:NDIM-1)
      double precision forces(NDIM*nphases*nphases)

      integer i, j, k, p1, p2, ic
      integer offset1, offset2
      double precision dxinv, dyinv, dzinv
      double precision diffx1, diffy1, diffz1
      double precision diffx2, diffy2, diffz2
      double precision fx, fy, fz
      double precision pp, cdiff, cutoff_slope, cfactor
c
      double precision interp_func
c
      cutoff_slope = 10.d0
c
      dxinv = 0.5d0 / dx(0)
      dyinv = 0.5d0 / dx(1)
      dzinv = 0.5d0 / dx(2)
c
      do p1 = 1, nphases
         do p2 = 1, p1-1
           offset1 = NDIM*(nphases*(p1-1)+(p2-1))
           offset2 = NDIM*(nphases*(p2-1)+(p1-1))
c integral over domain
            do k = lo2, hi2
               do j = lo1, hi1
                  do i = lo0, hi0
c smooth cutoff
                     pp = 0.5d0
     &                  + (phi(i,j,k,p1)*phi(i,j,k,p2)-cthreshold)
     &                    *cutoff_slope
c screen out zero contributions
c if pp<0., cfactor=0.
                     if( pp .gt. 0d0 )then
                        cfactor = interp_func( pp, 'p')
                        cdiff = 0.d0
                        do ic = 1, nc
                           cdiff = cdiff + conc(i,j,k,ic)
                        enddo
                        cdiff = cdiff-rhoe

                        diffx1 = dxinv * (
     &                     phi(i+1,j,k,p1) - phi(i-1,j,k,p1))
                        diffy1 = dyinv * (
     &                     phi(i,j+1,k,p1) - phi(i,j-1,k,p1))
                        diffz1 = dzinv * (
     &                     phi(i,j,k+1,p1) - phi(i,j,k-1,p1))

                        diffx2 = dxinv * (
     &                     phi(i+1,j,k,p2) - phi(i-1,j,k,p2))
                        diffy2 = dyinv * (
     &                     phi(i,j+1,k,p2) - phi(i,j-1,k,p2))
                        diffz2 = dzinv * (
     &                     phi(i,j,k+1,p2) - phi(i,j,k-1,p2))

                        fx = cdiff*cfactor*(diffx1-diffx2)
                        fy = cdiff*cfactor*(diffy1-diffy2)
                        fz = cdiff*cfactor*(diffz1-diffz2)
c add force density to integral over domain
                        forces(offset1+1) = forces(offset1+1)
     &                                    + fx*weight(i,j,k)
                        forces(offset1+2) = forces(offset1+2)
     &                                    + fy*weight(i,j,k)
                        forces(offset1+3) = forces(offset1+3)
     &                                    + fz*weight(i,j,k)
                        forces(offset2+1) = forces(offset2+1)
     &                                    - fx*weight(i,j,k)
                        forces(offset2+2) = forces(offset2+2)
     &                                    - fy*weight(i,j,k)
                        forces(offset2+3) = forces(offset2+3)
     &                                    - fz*weight(i,j,k)
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end
