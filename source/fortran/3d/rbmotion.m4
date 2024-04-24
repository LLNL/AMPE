c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE.
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
define(NDIM,3)dnl
define(REAL,`double precision')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c***********************************************************************
c
      subroutine addrbmotion(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   dx,
     &   phase, ngphase,
     &   vel,
     &   rhs, ngrhs)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2

      double precision dx(3), vel(3)
      integer ngphase, ngrhs
c
c variables in 3d cell indexed
      double precision phase(CELL3d(ifirst,ilast,ngphase))
      double precision rhs(CELL3d(ifirst,ilast,ngrhs))
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1, ic2
      double precision dxinv, dyinv, dzinv
      double precision diff_x, diff_y, diff_z

      dxinv = 0.5d0 * vel(1) / dx(1)
      dyinv = 0.5d0 * vel(2) / dx(2)
      dzinv = 0.5d0 * vel(3) / dx(3)
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               diff_x = dxinv *
     &            (phase(ic0+1,ic1,ic2) - phase(ic0-1,ic1,ic2))
               diff_y = dyinv *
     &            (phase(ic0,ic1+1,ic2) - phase(ic0,ic1-1,ic2))
               diff_z = dzinv *
     &            (phase(ic0,ic1,ic2+1) - phase(ic0,ic1,ic2-1))

               rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2)
     &                          - diff_x - diff_y - diff_z

            enddo
         enddo
      enddo

      end
