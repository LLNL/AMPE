c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE.
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
define(NDIM,3)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl

c***********************************************************************
c
c compute the concentration flux
c
      subroutine add_flux(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   dx,
     &   conc, ngconc,
     &   ncomp,
     &   diffconc0,  diffconc1, diffconc2,  ngdiff,
     &   flux0, flux1, flux2, ngflux )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngflux
      double precision dx(0:2)
      integer ncomp
      integer ngconc
      integer ngdiff
c
c variables in 3d cell indexed
      double precision
     &     flux0(SIDE3d0(ifirst,ilast,ngflux),ncomp),
     &     flux1(SIDE3d1(ifirst,ilast,ngflux),ncomp),
     &     flux2(SIDE3d2(ifirst,ilast,ngflux),ncomp)
      double precision conc(CELL3d(ifirst,ilast,ngconc),ncomp)
      double precision diffconc0(SIDE3d0(ifirst,ilast,ngdiff),
     &                           ncomp*ncomp)
      double precision diffconc1(SIDE3d1(ifirst,ilast,ngdiff),
     &                           ncomp*ncomp)
      double precision diffconc2(SIDE3d2(ifirst,ilast,ngdiff),
     &                           ncomp*ncomp)
c
      double precision dxinv, dyinv, dzinv
      integer          ic0, ic1, ic2, ic, jc, ijc

      dxinv = 1.d0 / dx(0)
      dyinv = 1.d0 / dx(1)
      dzinv = 1.d0 / dx(2)

      do ic = 1, ncomp
         do jc = 1, ncomp
            ijc=ic+(jc-1)*ncomp
            do ic2 = ifirst2, ilast2
               do ic1 = ifirst1, ilast1
                  do ic0 = ifirst0, ilast0+1
                     flux0(ic0,ic1,ic2,ic) = flux0(ic0,ic1,ic2,ic)
     &                  +dxinv * (
     &                   diffconc0(ic0,ic1,ic2,ijc) *
     &                   (conc(ic0,ic1,ic2,jc)-conc(ic0-1,ic1,ic2,jc))
     &                   )
                  enddo
               enddo
            enddo

            do ic2 = ifirst2, ilast2
               do ic1 = ifirst1, ilast1+1
                  do ic0 = ifirst0, ilast0
                     flux1(ic0,ic1,ic2,ic) = flux1(ic0,ic1,ic2,ic)
     &                  +dyinv * (
     &                   diffconc1(ic0,ic1,ic2,ijc) *
     &                   (conc(ic0,ic1,ic2,jc)-conc(ic0,ic1-1,ic2,jc))
     &                   )
                  enddo
               enddo
            enddo

            do ic2 = ifirst2, ilast2+1
               do ic1 = ifirst1, ilast1
                  do ic0 = ifirst0, ilast0
                     flux2(ic0,ic1,ic2,ic) = flux2(ic0,ic1,ic2,ic)
     &                  +dzinv * (
     &                   diffconc2(ic0,ic1,ic2,ijc) *
     &                   (conc(ic0,ic1,ic2,jc)-conc(ic0,ic1,ic2-1,jc))
     &                   )
                  enddo
               enddo
            enddo
         enddo

      enddo

      return
      end

c***********************************************************************
c
c compute the concentration flux using "isotropic" discretization
c see Shukla and Giri, JCP 276 (2014)
c
      subroutine add_flux_iso(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   dx,
     &   conc, ngconc, ncomp,
     &   d0,  d1, d2,  ngdiff,
     &   flux0, flux1, flux2, ngflux )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngflux
      double precision dx(0:2)
      integer ncomp
      integer ngconc, ngdiff
c
c variables in 3d cell indexed
      double precision
     &     flux0(SIDE3d0(ifirst,ilast,ngflux),ncomp),
     &     flux1(SIDE3d1(ifirst,ilast,ngflux),ncomp),
     &     flux2(SIDE3d2(ifirst,ilast,ngflux),ncomp)
      double precision conc(CELL3d(ifirst,ilast,ngconc),ncomp)
      double precision d0(SIDE3d0(ifirst,ilast,ngdiff),
     &                           ncomp*ncomp)
      double precision d1(SIDE3d1(ifirst,ilast,ngdiff),
     &                           ncomp*ncomp)
      double precision d2(SIDE3d2(ifirst,ilast,ngdiff),
     &                           ncomp*ncomp)
c
      double precision dx1, dx8, dy1, dy8, dz1, dz8
      integer          ic0, ic1, ic2, ic, jc, ijc

      dx1 = 1.d0 / (12.d0*dx(0))
      dx8 = 8.d0 / (12.d0*dx(0))
      dy1 = 1.d0 / (12.d0*dx(1))
      dy8 = 8.d0 / (12.d0*dx(1))
      dz1 = 1.d0 / (12.d0*dx(2))
      dz8 = 8.d0 / (12.d0*dx(2))

      do ic = 1, ncomp
      do jc = 1, ncomp
         ijc=ic+(jc-1)*ncomp
         do ic2 = ifirst2, ilast2
            do ic1 = ifirst1, ilast1
               do ic0 = ifirst0, ilast0+1
                  flux0(ic0,ic1,ic2,ic) = flux0(ic0,ic1,ic2,ic)
     &               +dx1 * (
     &                d0(ic0,ic1-1,ic2,ijc) *
     &                (conc(ic0,ic1-1,ic2,jc)-conc(ic0-1,ic1-1,ic2,jc))
     &                )
     &               +dx1 * (
     &                d0(ic0,ic1,ic2-1,ijc) *
     &                (conc(ic0,ic1,ic2-1,jc)-conc(ic0-1,ic1,ic2-1,jc))
     &                )
     &               +dx1 * (
     &                d0(ic0,ic1+1,ic2,ijc) *
     &                (conc(ic0,ic1+1,ic2,jc)-conc(ic0-1,ic1+1,ic2,jc))
     &                )
     &               +dx1 * (
     &                d0(ic0,ic1,ic2+1,ijc) *
     &                (conc(ic0,ic1,ic2+1,jc)-conc(ic0-1,ic1,ic2+1,jc))
     &                )
     &               +dx8 * (
     &                d0(ic0,ic1,ic2,ijc) *
     &                (conc(ic0,ic1,ic2,jc)-conc(ic0-1,ic1,ic2,jc))
     &                )
               enddo
            enddo
         enddo

         do ic2 = ifirst2, ilast2
            do ic1 = ifirst1, ilast1+1
               do ic0 = ifirst0, ilast0
                  flux1(ic0,ic1,ic2,ic) = flux1(ic0,ic1,ic2,ic)
     &               +dy1 * (
     &                d1(ic0-1,ic1,ic2,ijc) *
     &                (conc(ic0-1,ic1,ic2,jc)-conc(ic0-1,ic1-1,ic2,jc))
     &                )
     &               +dy1 * (
     &                d1(ic0+1,ic1,ic2,ijc) *
     &                (conc(ic0+1,ic1,ic2,jc)-conc(ic0+1,ic1-1,ic2,jc))
     &                )
     &               +dy1 * (
     &                d1(ic0,ic1,ic2-1,ijc) *
     &                (conc(ic0,ic1,ic2-1,jc)-conc(ic0,ic1-1,ic2-1,jc))
     &                )
     &               +dy1 * (
     &                d1(ic0,ic1,ic2+1,ijc) *
     &                (conc(ic0,ic1,ic2+1,jc)-conc(ic0,ic1-1,ic2+1,jc))
     &                )
     &               +dy8 * (
     &                d1(ic0,ic1,ic2,ijc) *
     &                (conc(ic0,ic1,ic2,jc)-conc(ic0,ic1-1,ic2,jc))
     &                )
               enddo
            enddo
         enddo

         do ic2 = ifirst2, ilast2+1
            do ic1 = ifirst1, ilast1
               do ic0 = ifirst0, ilast0
                  flux2(ic0,ic1,ic2,ic) = flux2(ic0,ic1,ic2,ic)
     &               +dz1 * (
     &                d2(ic0-1,ic1,ic2,ijc) *
     &                (conc(ic0-1,ic1,ic2,jc)-conc(ic0-1,ic1,ic2-1,jc))
     &                )
     &               +dz1 * (
     &                d2(ic0+1,ic1,ic2,ijc) *
     &                (conc(ic0+1,ic1,ic2,jc)-conc(ic0+1,ic1,ic2-1,jc))
     &                )
     &               +dz1 * (
     &                d2(ic0,ic1-1,ic2,ijc) *
     &                (conc(ic0,ic1-1,ic2,jc)-conc(ic0,ic1-1,ic2-1,jc))
     &                )
     &               +dz1 * (
     &                d2(ic0,ic1+1,ic2,ijc) *
     &                (conc(ic0,ic1+1,ic2,jc)-conc(ic0,ic1+1,ic2-1,jc))
     &                )
     &               +dz8 * (
     &                d2(ic0,ic1,ic2,ijc) *
     &                (conc(ic0,ic1,ic2,jc)-conc(ic0,ic1,ic2-1,jc))
     &                )
               enddo
            enddo
         enddo

      enddo
      enddo

      return
      end

