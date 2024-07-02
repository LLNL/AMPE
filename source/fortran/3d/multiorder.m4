c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE.
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
c
define(NDIM,3)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl

c***********************************************************************
c
c compute r.h.s. for multi order parameters
c
      subroutine computerhsmultiorder(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   norder, dx,
     &   flux0,
     &   flux1,
     &   flux2,
     &   ngflux,
     &   dgamma, m,
     &   phi, ngphi,
     &   phi2sum,
     &   rhs, ngrhs)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2, norder

      double precision dx(3)
      double precision dgamma, m
      integer ngflux, ngphi, ngrhs
c
c variables in 3d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi),norder)
      double precision rhs(CELL3d(ifirst,ilast,ngrhs),norder)
      double precision phi2sum(CELL3d(ifirst,ilast,0))

c variables in 3d side indexed
      double precision
     &     flux0(SIDE3d0(ifirst,ilast,ngflux),norder),
     &     flux1(SIDE3d1(ifirst,ilast,ngflux),norder),
     &     flux2(SIDE3d2(ifirst,ilast,ngflux),norder)
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1, ic2
      integer ip, jp
      double precision diff_term_x, diff_term_y, diff_term_z
      double precision diff_term

      double precision dsum, g_prime
      double precision deriv_interp_func
      double precision deriv_triple_well_func
      double precision dxinv, dyinv, dzinv, fac

c
      dxinv = 1.d0 / dx(1)
      dyinv = 1.d0 / dx(2)
      dzinv = 1.d0 / dx(3)
      fac = 1.d0 / norder
c
c precompute sum of phi**2 at each cell
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0
               phi2sum(ic0,ic1,ic2) =
     &            phi(ic0,ic1,ic2,1)*phi(ic0,ic1,ic2,1)
            enddo
         enddo
      enddo

      do ip = 2, norder
         do ic2 = ifirst2, ilast2
            do ic1 = ifirst1, ilast1
               do ic0 = ifirst0, ilast0

                  phi2sum(ic0,ic1,ic2) = phi2sum(ic0,ic1,ic2)
     &               + phi(ic0,ic1,ic2,ip)*phi(ic0,ic1,ic2,ip)
               enddo
            enddo
         enddo
      enddo
c
      do ip = 1, norder
         do ic2 = ifirst2, ilast2
            do ic1 = ifirst1, ilast1
               do ic0 = ifirst0, ilast0

                  diff_term_x = dxinv *
     &              (flux0(ic0+1,ic1,ic2,ip) - flux0(ic0,ic1,ic2,ip))
                  diff_term_y = dyinv *
     &              (flux1(ic0,ic1+1,ic2,ip) - flux1(ic0,ic1,ic2,ip))
                  diff_term_z = dzinv *
     &              (flux2(ic0,ic1,ic2+1,ip) - flux2(ic0,ic1,ic2,ip))

                  diff_term = diff_term_x + diff_term_y + diff_term_z

                  rhs(ic0,ic1,ic2,ip) = diff_term

c  Phase energy well

                  g_prime =
     &               phi(ic0,ic1,ic2,ip)
     &               * (phi(ic0,ic1,ic2,ip)*phi(ic0,ic1,ic2,ip)-1.d0)
                  dsum = phi2sum(ic0,ic1,ic2) - phi(ic0,ic1,ic2,ip)
                  g_prime = g_prime + 2.*phi(ic0,ic1,ic2,ip)*dgamma*dsum

                  rhs(ic0,ic1,ic2,ip) = rhs(ic0,ic1,ic2,ip) - m*g_prime

               enddo
            enddo
         enddo
      enddo

      return
      end
