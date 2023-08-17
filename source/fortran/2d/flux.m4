c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE.
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
define(NDIM,2)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine add_flux(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   dx,
     &   data, ng, depth,
     &   d0, d1, ngd,
     &   flux0, flux1, ngf )
c***********************************************************************
      implicit none
c***********************************************************************
      integer ifirst0, ilast0, ifirst1, ilast1
      double precision dx(0:1)
      integer depth
      integer ng, ngd, ngf
c
      double precision
     &     flux0(SIDE2d0(ifirst,ilast,ngf),depth),
     &     flux1(SIDE2d1(ifirst,ilast,ngf),depth)
      double precision data(CELL2d(ifirst,ilast,ng),depth)
      double precision d0(SIDE2d0(ifirst,ilast,ngd),
     &                            depth*depth)
      double precision d1(SIDE2d1(ifirst,ilast,ngd),
     &                            depth*depth)
c
      double precision dxinv, dyinv
      integer          ic0, ic1, i, j, ij

      dxinv = 1.d0 / dx(0)
      dyinv = 1.d0 / dx(1)

      do i = 1, depth
         do j = 1, depth
            ij=i+(j-1)*depth
            do ic1 = ifirst1, ilast1
               do ic0 = ifirst0, ilast0+1
                  flux0(ic0,ic1,i) = flux0(ic0,ic1,i) + dxinv *
     &              d0(ic0,ic1,ij) *
     &              ( data(ic0,ic1,j) - data(ic0-1,ic1,j) )
               enddo
            enddo

            do ic1 = ifirst1, ilast1+1
               do ic0 = ifirst0, ilast0
                  flux1(ic0,ic1,i) = flux1(ic0,ic1,i) + dyinv *
     &              d1(ic0,ic1,ij) *
     &              ( data(ic0,ic1,j) - data(ic0,ic1-1,j) )
               enddo
            enddo
         enddo
      enddo

      return
      end
c
c***********************************************************************
c
c compute dataentration flux using "isotropic" discretization
c see Shukla and Giri, JCP 276 (2014)
c
      subroutine add_flux_iso(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   dx,
     &   data, ng, depth,
     &   d0, d1, ngd,
     &   flux0, flux1, ngf )
c***********************************************************************
      implicit none
c***********************************************************************
      integer ifirst0, ilast0, ifirst1, ilast1
      double precision dx(0:1)
      integer depth
      integer ng, ngd, ngf
c
      double precision
     &     flux0(SIDE2d0(ifirst,ilast,ngf),depth),
     &     flux1(SIDE2d1(ifirst,ilast,ngf),depth)
      double precision data(CELL2d(ifirst,ilast,ng),depth)
      double precision d0(SIDE2d0(ifirst,ilast,ngd),
     &                            depth*depth)
      double precision d1(SIDE2d1(ifirst,ilast,ngd),
     &                            depth*depth)
c
      double precision dx10, dx1, dy10, dy1
      integer          i0, i1, i, j, ij

      dx10 = 10.d0 / (12.d0*dx(0))
      dx1  = 1.d0  / (12.d0*dx(0))
      dy10 = 10.d0 / (12.d0*dx(1))
      dy1  = 1.d0  / (12.d0*dx(1))

      do i = 1, depth
         do j = 1, depth
            ij=i+(j-1)*depth
            do i1 = ifirst1, ilast1
               do i0 = ifirst0, ilast0+1
                  flux0(i0,i1,i) = flux0(i0,i1,i)
     &              + dx1 * d0(i0,i1+1,ij) *
     &              ( data(i0,i1+1,j) - data(i0-1,i1+1,j) )
     &              + dx10 * d0(i0,i1,ij) *
     &              ( data(i0,i1,j) - data(i0-1,i1,j) )
     &              + dx1 * d0(i0,i1-1,ij) *
     &              ( data(i0,i1-1,j) - data(i0-1,i1-1,j) )
               enddo
            enddo

            do i1 = ifirst1, ilast1+1
               do i0 = ifirst0, ilast0
                  flux1(i0,i1,i) = flux1(i0,i1,i)
     &              + dy1 *d1(i0+1,i1,ij) *
     &              ( data(i0+1,i1,j) - data(i0+1,i1-1,j) )
     &              + dy10 *d1(i0,i1,ij) *
     &              ( data(i0,i1,j) - data(i0,i1-1,j) )
     &              + dy1 *d1(i0-1,i1,ij) *
     &              ( data(i0-1,i1,j) - data(i0-1,i1-1,j) )
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
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   dx,
     &   data, ng, depth,
     &   d0, d1, ngd,
     &   flux0, flux1, ngf, physb )
c***********************************************************************
      implicit none
c***********************************************************************
      integer ifirst0, ilast0, ifirst1, ilast1
      double precision dx(0:1)
      integer depth
      integer ng, ngd, ngf
c
      double precision
     &     flux0(SIDE2d0(ifirst,ilast,ngf),depth),
     &     flux1(SIDE2d1(ifirst,ilast,ngf),depth)
      double precision data(CELL2d(ifirst,ilast,ng),depth)
      double precision d0(SIDE2d0(ifirst,ilast,ngd),
     &                            depth*depth)
      double precision d1(SIDE2d1(ifirst,ilast,ngd),
     &                            depth*depth)
      integer physb(4)
c
      double precision dx15, dx1, dy15, dy1, dxi, dyi
      integer          i0, i1, i, j, ij
      integer          if0, il0, if1, il1

      dxi = 1.d0  / dx(0)
      dyi = 1.d0  / dx(1)

      dx15 = 15.d0 / (12.d0*dx(0))
      dx1  = 1.d0  / (12.d0*dx(0))
      dy15 = 15.d0 / (12.d0*dx(1))
      dy1  = 1.d0  / (12.d0*dx(1))

      do i = 1, depth
         do j = 1, depth
            ij=i+(j-1)*depth
c x flux
            do i1 = ifirst1, ilast1
               if( physb(1) .eq. 1 ) then
                 i0 = ifirst0
                 flux0(i0,i1,i) = flux0(i0,i1,i)
     &              + d0(i0,i1,ij) *
     &                dxi *( data(i0,i1,j) - data(i0-1,i1,j) )
                 if0 = ifirst0+1
               else
                 if0 = ifirst0
               endif
               if( physb(2) .eq. 1 ) then
                 i0 = ilast0+1
                 flux0(i0,i1,i) = flux0(i0,i1,i)
     &              + d0(i0,i1,ij) *
     &                dxi *( data(i0,i1,j) - data(i0-1,i1,j) )
                 il0 = ilast0
               else
                 il0 = ilast0+1
               endif

               do i0 = if0, il0
                  flux0(i0,i1,i) = flux0(i0,i1,i)
     &              + d0(i0,i1,ij) * (
     &                dx15 *( data(i0,i1,j) - data(i0-1,i1,j) )
     &              - dx1 *( data(i0+1,i1,j) - data(i0-2,i1,j) ) )
               enddo
            enddo
c y flux
            if( physb(3) .eq. 1 ) then
               if1 =  ifirst1+1
               i1 = ifirst1
               do i0 = ifirst0, ilast0
                  flux1(i0,i1,i) = flux1(i0,i1,i)
     &              + d1(i0,i1,ij) * dyi *( data(i0,i1,j)
     &                                    - data(i0,i1-1,j) )
               enddo
            else
              if1 = ifirst1
            endif
            if( physb(4) .eq. 1 ) then
               il1 =  ilast1
               i1 = ilast1+1
               do i0 = ifirst0, ilast0
                  flux1(i0,i1,i) = flux1(i0,i1,i)
     &              + d1(i0,i1,ij) *dyi *( data(i0,i1,j)
     &                                   - data(i0,i1-1,j) )
               enddo
            else
               il1 = ilast1+1
            endif

            do i1 = if1, il1
               do i0 = ifirst0, ilast0
                  flux1(i0,i1,i) = flux1(i0,i1,i)
     &              + d1(i0,i1,ij) * (
     &                dy15 *( data(i0,i1,j) - data(i0,i1-1,j) )
     &              - dy1 *( data(i0,i1+1,j) - data(i0,i1-2,j) ) )
               enddo
            enddo
         enddo
      enddo

      return
      end

