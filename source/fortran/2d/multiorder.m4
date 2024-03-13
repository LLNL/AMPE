c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE.
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
c
define(NDIM,2)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl

c***********************************************************************
c
c compute r.h.s. for multi order parameters
c
      subroutine computerhsmultiorder(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   norder, dx,
     &   flux0,
     &   flux1,
     &   ngflux,
     &   dgamma, m,
     &   phi, ngphi,
     &   rhs, ngrhs)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, norder

      double precision dx(2)
      double precision dgamma, m
      integer ngflux, ngphi, ngrhs
c
c variables in 2d cell indexed
      double precision phi(CELL2d(ifirst,ilast,ngphi),norder)
      double precision rhs(CELL2d(ifirst,ilast,ngrhs),norder)

c variables in 2d side indexed
      double precision
     &     flux0(SIDE2d0(ifirst,ilast,ngflux),norder),
     &     flux1(SIDE2d1(ifirst,ilast,ngflux),norder)
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1
      integer ip, jp
      double precision diff_term_x, diff_term_y, diff_term

      double precision dsum, g_prime
      double precision deriv_interp_func
      double precision deriv_triple_well_func
      double precision dxinv, dyinv, fac

c
      dxinv = 1.d0 / dx(1)
      dyinv = 1.d0 / dx(2)
      fac = 1.d0 / norder
c
      do ip = 1, norder
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               diff_term_x =
     &              (flux0(ic0+1,ic1,ip) - flux0(ic0,ic1,ip)) * dxinv
               diff_term_y =
     &              (flux1(ic0,ic1+1,ip) - flux1(ic0,ic1,ip)) * dyinv

               diff_term = diff_term_x + diff_term_y

               rhs(ic0,ic1,ip) = diff_term

c  Phase energy well

               g_prime =
     &            phi(ic0,ic1,ip)*phi(ic0,ic1,ip)*phi(ic0,ic1,ip)
     &            - phi(ic0,ic1,ip)
               dsum = 0.d0
               do jp = 1, norder
                 dsum = dsum +  phi(ic0,ic1,jp)*phi(ic0,ic1,jp)
               enddo
               dsum = dsum - phi(ic0,ic1,ip)
               g_prime = g_prime + 2.*phi(ic0,ic1,ip)*dgamma*dsum

               rhs(ic0,ic1,ip) = rhs(ic0,ic1,ip) - m * g_prime

            enddo
         enddo
      enddo

      return
      end
