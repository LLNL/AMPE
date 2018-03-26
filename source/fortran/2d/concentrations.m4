c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
c Redistribution and use in source and binary forms, with or without 
c modification, are permitted provided that the following conditions are met:
c - Redistributions of source code must retain the above copyright notice,
c   this list of conditions and the disclaimer below.
c - Redistributions in binary form must reproduce the above copyright notice,
c   this list of conditions and the disclaimer (as noted below) in the
c   documentation and/or other materials provided with the distribution.
c - Neither the name of the LLNS/LLNL nor the names of its contributors may be
c   used to endorse or promote products derived from this software without
c   specific prior written permission.
c
c THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
c AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
c IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
c ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
c LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
c DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
c DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
c OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
c HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
c STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
c IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
c POSSIBILITY OF SUCH DAMAGE.
c 
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
