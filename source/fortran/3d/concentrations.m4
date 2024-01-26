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
c
c Compute c=h(phi)*cs+(1-h(phi))*cl
c
      subroutine compute_concentration_from_phase_concentrations(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phi, ngphi,
     &   cl, cs, ngcphase,
     &   c, ngc,
     &   interp_type )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngphi, ngcphase, ngc
      character*(*) interp_type
c
c variables in 2d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi))
      double precision cl(CELL3d(ifirst,ilast,ngcphase))
      double precision cs(CELL3d(ifirst,ilast,ngcphase))
      double precision c(CELL3d(ifirst,ilast,ngc))
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0, ic1, ic2
      double precision hphi
      double precision interp_func
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               hphi = interp_func( phi(ic0,ic1,ic2), interp_type )

               c(ic0,ic1,ic2) =
     &            ( 1.0d0 - hphi ) * cl(ic0,ic1,ic2) +
     &                        hphi * cs(ic0,ic1,ic2)

           end do
         end do
      end do
c
      return
      end
