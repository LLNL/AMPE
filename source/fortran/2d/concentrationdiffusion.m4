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
c Coefficient [h(phi)*d_solid+(1-h(phi))*d_liquid] from S.Y. Hu et al.
c
      subroutine concentrationdiffusion0(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phi, ngphi,
     &   eta, ngeta,
     &   diff0, diff1, ngdiff,
     &   temp, ngtemp,
     &   d_liquid, q0_liquid,
     &   d_solid_A, q0_solid_A,
     &   d_solid_B, q0_solid_B,
     &   gas_constant_R,
     &   interp_type,
     &   avg_type,
     &   three_phase )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer ngphi, ngeta, ngdiff, ngtemp
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solid_A, d_solid_B
      double precision q0_liquid, q0_solid_A, q0_solid_B
      double precision gas_constant_R
      integer three_phase
c
c variables in 2d cell indexed
      double precision phi(CELL2d(ifirst,ilast,ngphi))
      double precision eta(CELL2d(ifirst,ilast,ngeta))
      double precision temp(CELL2d(ifirst,ilast,ngtemp))
      double precision diff0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE2d1(ifirst,ilast,ngdiff))
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0, ic1
      double precision vphi, veta, hphi, heta, invT
      double precision q0_liquid_invR, q0_solid_A_invR, q0_solid_B_invR
      double precision diff_liquid, diff_solid_A, diff_solid_B
      double precision interp_func
      double precision average_func
c
      q0_liquid_invR = q0_liquid / gas_constant_R
      q0_solid_A_invR = q0_solid_A / gas_constant_R
      q0_solid_B_invR = q0_solid_B / gas_constant_R
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0+1

            vphi = average_func(
     &         phi(ic0-1,ic1), phi(ic0,ic1), avg_type )
            
            hphi = interp_func( vphi, interp_type )

            invT = 2.0d0 / ( temp(ic0-1,ic1) + temp(ic0,ic1) )

            diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
            diff_solid_A = d_solid_A * exp( -q0_solid_A_invR * invT )

            if ( three_phase /= 0 ) then
               veta = average_func(
     &            eta(ic0-1,ic1), eta(ic0,ic1), avg_type )

               heta = interp_func( veta, interp_type )
               diff_solid_B =
     &            d_solid_B * exp( -q0_solid_B_invR * invT )
            else
               heta = 0.0d0
               diff_solid_B = 0.0d0
            endif

            diff0(ic0,ic1) =
     &         ( 1.0d0 - hphi ) * diff_liquid +
     &         hphi * (
     &            ( 1.0d0 - heta ) * diff_solid_A +
     &            heta * diff_solid_B 
     &         )

         end do
      end do
c
      do ic1 = ifirst1, ilast1+1
         do ic0 = ifirst0, ilast0

            vphi = average_func(
     &         phi(ic0,ic1-1), phi(ic0,ic1), avg_type )

            hphi = interp_func( vphi, interp_type )

            invT = 2.0d0 / ( temp(ic0,ic1-1) + temp(ic0,ic1) )

            diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
            diff_solid_A = d_solid_A * exp( -q0_solid_A_invR * invT )

            if ( three_phase /= 0 ) then
               veta = average_func(
     &            eta(ic0,ic1-1), eta(ic0,ic1), avg_type )

               heta = interp_func( veta, interp_type )
               diff_solid_B =
     &            d_solid_B * exp( -q0_solid_B_invR * invT )
            else
               heta = 0.0d0
               diff_solid_B = 0.0d0
            endif

            diff1(ic0,ic1) =
     &         ( 1.0d0 - hphi ) * diff_liquid +
     &         hphi * (
     &            ( 1.0d0 - heta ) * diff_solid_A +
     &            heta * diff_solid_B 
     &         )

         end do
      end do
c
      return
      end

c
c Coefficient \tilde D from Acharya, Sharon, Starolesky (2016)
c
      subroutine concentrationdiffusion_utrc(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phi, ngphi,
     &   diff0, diff1, ngdiff,
     &   partition_coeff, ngk,
     &   d_liquid,
     &   d_solid,
     &   interp_type,
     &   avg_type )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer ngphi, ngdiff, ngk
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solid
c
c variables in 2d cell indexed
      double precision phi(CELL2d(ifirst,ilast,ngphi))
      double precision partition_coeff(CELL2d(ifirst,ilast,ngk))
      double precision diff0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE2d1(ifirst,ilast,ngdiff))
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0, ic1
      double precision vphi, hphi, k
      double precision diff_liquid, diff_solid
      double precision interp_func
      double precision average_func
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0+1

            vphi = average_func(
     &         phi(ic0-1,ic1), phi(ic0,ic1), avg_type )
            
            hphi = interp_func( vphi, interp_type )
            
            k = partition_coeff(ic0,ic1)

            diff0(ic0,ic1) = d_solid +
     &         (d_liquid-d_solid)*( 1.0d0 - hphi )
     &                   /(1.0d0 - hphi+k*hphi)

         end do
      end do
c
      do ic1 = ifirst1, ilast1+1
         do ic0 = ifirst0, ilast0

            vphi = average_func(
     &         phi(ic0,ic1-1), phi(ic0,ic1), avg_type )

            hphi = interp_func( vphi, interp_type )

            k = partition_coeff(ic0,ic1)

            diff1(ic0,ic1) = d_solid +
     &         (d_liquid-d_solid)*( 1.0d0 - hphi )
     &                   /(1.0d0 - hphi+k*hphi)

         end do
      end do
c
      return
      end
