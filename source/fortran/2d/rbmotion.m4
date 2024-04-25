c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE.
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c***********************************************************************
c
      subroutine addrbmotion(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   dx,
     &   phase, ngphase,
     &   vel,
     &   rhs, ngrhs)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1

      double precision dx(2), vel(2)
      integer ngphase, ngrhs
c
c variables in 2d cell indexed
      double precision phase(CELL2d(ifirst,ilast,ngphase))
      double precision rhs(CELL2d(ifirst,ilast,ngrhs))
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1
      double precision dxinv, dyinv
      double precision diff_x, diff_y

      dxinv = 0.5d0 * vel(1) / dx(1)
      dyinv = 0.5d0 * vel(2) / dx(2)
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0

            diff_x = dxinv *
     &           (phase(ic0+1,ic1) - phase(ic0-1,ic1))
            diff_y = dyinv *
     &           (phase(ic0,ic1+1) - phase(ic0,ic1-1))

            rhs(ic0,ic1) = rhs(ic0,ic1)
     &                   - diff_x - diff_y

         enddo
      enddo

      end

