c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
define(NDIM,2)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
c Compute c=h(phi)*cs+(1-h(phi))*cl
c
      subroutine compute_concentration_from_phase_concentrations(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phi, ngphi,
     &   cl, cs, ngcphase,
     &   c, ngc,
     &   interp_type )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer ngphi, ngcphase, ngc
      character*(*) interp_type
c
c variables in 2d cell indexed
      double precision phi(CELL2d(ifirst,ilast,ngphi))
      double precision cl(CELL2d(ifirst,ilast,ngcphase))
      double precision cs(CELL2d(ifirst,ilast,ngcphase))
      double precision c(CELL2d(ifirst,ilast,ngc))
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0, ic1
      double precision hphi
      double precision interp_func
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0

            hphi = interp_func( phi(ic0,ic1), interp_type )

            c(ic0,ic1) =
     &         ( 1.0d0 - hphi ) * cl(ic0,ic1) +
     &                     hphi * cs(ic0,ic1)

         end do
      end do
c
      return
      end
