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
c
c***********************************************************************
c
c compute flux using fourth-order scheme
c see Zhang, Johansen, Colella (2012), eq. (2.5)
c requires two ghost cells for data
c array physb tells if sides touch physical boundaries
c
      subroutine add_flux_4th(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   dx,
     &   data, ng, depth,
     &   d0, d1, d2, ngd,
     &   flux0, flux1, flux2, ngf, physb )
c***********************************************************************
      implicit none
c***********************************************************************
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      double precision dx(0:2)
      integer depth
      integer ng, ngd, ngf
c
      double precision
     &     flux0(SIDE3d0(ifirst,ilast,ngf),depth),
     &     flux1(SIDE3d1(ifirst,ilast,ngf),depth),
     &     flux2(SIDE3d2(ifirst,ilast,ngf),depth)
      double precision data(CELL3d(ifirst,ilast,ng),depth)
      double precision d0(SIDE3d0(ifirst,ilast,ngd),
     &                            depth*depth)
      double precision d1(SIDE3d1(ifirst,ilast,ngd),
     &                            depth*depth)
      double precision d2(SIDE3d2(ifirst,ilast,ngd),
     &                            depth*depth)
      integer physb(6)
c
      double precision dx15, dx1, dy15, dy1, dz15, dz1
      double precision dxi, dyi, dzi
      integer          i0, i1, i2, i, j, ij
      integer          if0, il0, if1, il1, if2, il2

      dxi = 1.d0  / dx(0)
      dyi = 1.d0  / dx(1)
      dzi = 1.d0  / dx(2)

      dx15 = 15.d0 / (12.d0*dx(0))
      dx1  = 1.d0  / (12.d0*dx(0))
      dy15 = 15.d0 / (12.d0*dx(1))
      dy1  = 1.d0  / (12.d0*dx(1))
      dz15 = 15.d0 / (12.d0*dx(2))
      dz1  = 1.d0  / (12.d0*dx(2))

      do i = 1, depth
         do j = 1, depth
            ij=i+(j-1)*depth
            do i2 = ifirst2, ilast2
               do i1 = ifirst1, ilast1
                  if( physb(1) .eq. 1 ) then
                     i0 = ifirst0
                     flux0(i0,i1,i2,i) = flux0(i0,i1,i2,i)
     &                  + d0(i0,i1,i2,ij) *
     &                    dxi *( data(i0,i1,i2,j) - data(i0-1,i1,i2,j) )
                     if0 = ifirst0+1
                  else
                     if0 = ifirst0
                  endif
                  if( physb(2) .eq. 1 ) then
                     i0 = ilast0+1
                     flux0(i0,i1,i2,i) = flux0(i0,i1,i2,i)
     &                 + d0(i0,i1,i2,ij) *
     &                   dxi *( data(i0,i1,i2,j) - data(i0-1,i1,i2,j) )
                     il0 = ilast0
                  else
                     il0 = ilast0+1
                  endif

                  do i0 = if0, il0
                     flux0(i0,i1,i2,i) = flux0(i0,i1,i2,i)
     &                  + d0(i0,i1,i2,ij) * (
     &                   dx15 *
     &                   (data(i0,i1,i2,j)-data(i0-1,i1,i2,j))
     &                 - dx1 *
     &                   (data(i0+1,i1,i2,j)-data(i0-2,i1,i2,j))
     &                   )
                  enddo
               enddo
            enddo
c y flux
            do i2 = ifirst2, ilast2
c first special case for physical boundaries in y-direction
               if( physb(3) .eq. 1 ) then
                  if1 = ifirst1+1
                  i1 = ifirst1
                  do i0 = ifirst0, ilast0
                     flux1(i0,i1,i2,i) = flux1(i0,i1,i2,i)
     &                 + d1(i0,i1,i2,ij) * dyi *( data(i0,i1,i2,j)
     &                                    - data(i0,i1-1,i2,j) )
                  enddo
               else
                 if1 = ifirst1
               endif
               if( physb(4) .eq. 1 ) then
                  il1 = ilast1
                  i1 = ilast1+1
                  do i0 = ifirst0, ilast0
                     flux1(i0,i1,i2,i) = flux1(i0,i1,i2,i)
     &                 + d1(i0,i1,i2,ij) *dyi *( data(i0,i1,i2,j)
     &                                      - data(i0,i1-1,i2,j) )
                  enddo
               else
                  il1 = ilast1+1
               endif

               do i1 = if1, il1
                  do i0 = ifirst0, ilast0
                     flux1(i0,i1,i2,i) = flux1(i0,i1,i2,i)
     &                  +d1(i0,i1,i2,ij) * (
     &                   dy15 *
     &                   (data(i0,i1,i2,j)-data(i0,i1-1,i2,j))
     &                  - dy1 *
     &                   (data(i0,i1+1,i2,j)-data(i0,i1-2,i2,j))
     &                   )
                  enddo
               enddo
            enddo
c z flux
c first special case for physical boundaries in z-direction
            if( physb(5) .eq. 1 ) then
               if2 = ifirst2+1
               i2 = ifirst2
               do i1 = ifirst1, ilast1
                  do i0 = ifirst0, ilast0
                     flux2(i0,i1,i2,i) = flux2(i0,i1,i2,i)
     &                  +  d2(i0,i1,i2,ij) * (
     &                   dzi *
     &                   (data(i0,i1,i2,j)-data(i0,i1,i2-1,j))
     &                   )
                  enddo
               enddo
            else
              if2 = ifirst2
            endif
            if( physb(6) .eq. 1 ) then
               il2 = ilast2
               i2 = ilast2+1
               do i1 = ifirst1, ilast1
                  do i0 = ifirst0, ilast0
                     flux2(i0,i1,i2,i) = flux2(i0,i1,i2,i)
     &                  +  d2(i0,i1,i2,ij) * (
     &                   dzi *
     &                   (data(i0,i1,i2,j)-data(i0,i1,i2-1,j))
     &                   )
                  enddo
               enddo
            else
               il2 = ilast2+1
            endif
            do i2 = if2, il2
               do i1 = ifirst1, ilast1
                  do i0 = ifirst0, ilast0
                     flux2(i0,i1,i2,i) = flux2(i0,i1,i2,i)
     &                  +  d2(i0,i1,i2,ij) * (
     &                   dz15 *
     &                   (data(i0,i1,i2,j)-data(i0,i1,i2-1,j))
     &                  - dz1 *
     &                   (data(i0,i1,i2+1,j)-data(i0,i1,i2-2,j))
     &                   )
                  enddo
               enddo
            enddo
         enddo

      enddo

      return
      end

