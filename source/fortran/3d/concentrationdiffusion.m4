c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
c 
define(NDIM,3)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c
c Coefficient [h(phi)*d_solid+(1-h(phi))*d_liquid]
c
      subroutine concentration_pfmdiffusion(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phi, ngphi,
     &   eta, ngeta,
     &   diff0, diff1, diff2, ngdiff,
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
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngphi, ngeta, ngdiff, ngtemp
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solid_A, d_solid_B
      double precision q0_liquid, q0_solid_A, q0_solid_B
      double precision gas_constant_R
      integer three_phase
c
c variables in 3d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi))
      double precision eta(CELL3d(ifirst,ilast,ngeta))
      double precision temp(CELL3d(ifirst,ilast,ngtemp))
      double precision diff0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diff2(SIDE3d2(ifirst,ilast,ngdiff))
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0, ic1, ic2
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
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1

               vphi = average_func(
     &            phi(ic0-1,ic1,ic2), phi(ic0,ic1,ic2), avg_type )

               hphi = interp_func( vphi, interp_type )

               invT = 2.0d0 /
     &            ( temp(ic0-1,ic1,ic2) + temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solid_A = d_solid_A * exp( -q0_solid_A_invR * invT )

               if ( three_phase /= 0 ) then
                  veta = average_func(
     &               eta(ic0-1,ic1,ic2), eta(ic0,ic1,ic2), avg_type )

                  heta = interp_func( veta, interp_type )
                  diff_solid_B =
     &               d_solid_B * exp( -q0_solid_B_invR * invT )
               else
                  heta = 0.0d0
                  diff_solid_B = 0.0d0
               endif

               diff0(ic0,ic1,ic2) =
     &            ( 1.0d0 - hphi ) * diff_liquid +
     &            hphi * (
     &               ( 1.0d0 - heta ) * diff_solid_A +
     &               heta * diff_solid_B 
     &            )

            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0, ilast0

               vphi = average_func(
     &            phi(ic0,ic1-1,ic2), phi(ic0,ic1,ic2), avg_type )

               hphi = interp_func( vphi, interp_type )

               invT = 2.0d0 /
     &            ( temp(ic0,ic1-1,ic2) + temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solid_A = d_solid_A * exp( -q0_solid_A_invR * invT )

               if ( three_phase /= 0 ) then
                  veta = average_func(
     &               eta(ic0,ic1-1,ic2), eta(ic0,ic1,ic2), avg_type )

                  heta = interp_func( veta, interp_type )
                  diff_solid_B =
     &               d_solid_B * exp( -q0_solid_B_invR * invT )
               else
                  heta = 0.0d0
                  diff_solid_B = 0.0d0
               endif

               diff1(ic0,ic1,ic2) =
     &            ( 1.0d0 - hphi ) * diff_liquid +
     &            hphi * (
     &               ( 1.0d0 - heta ) * diff_solid_A +
     &               heta * diff_solid_B 
     &            )

            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               vphi = average_func(
     &            phi(ic0,ic1,ic2-1), phi(ic0,ic1,ic2), avg_type )

               hphi = interp_func( vphi, interp_type )

               invT = 2.0d0 /
     &            ( temp(ic0,ic1,ic2-1) + temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solid_A = d_solid_A * exp( -q0_solid_A_invR * invT )

               if ( three_phase /= 0 ) then
                  veta = average_func(
     &               eta(ic0,ic1,ic2-1), eta(ic0,ic1,ic2), avg_type )

                  heta = interp_func( veta, interp_type )
                  diff_solid_B =
     &               d_solid_B * exp( -q0_solid_B_invR * invT )
               else
                  heta = 0.0d0
                  diff_solid_B = 0.0d0
               endif

               diff2(ic0,ic1,ic2) =
     &            ( 1.0d0 - hphi ) * diff_liquid +
     &            hphi * (
     &               ( 1.0d0 - heta ) * diff_solid_A +
     &               heta * diff_solid_B 
     &            )

            end do
         end do
      end do
c
      return
      end
c
c same as function concentrationdiffusion0, without accumulating
c component into single D
c
      subroutine concentration_pfmdiffusion_of_temperature(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phi, ngphi,
     &   diffL0, diffL1, diffL2, 
     &   diffA0, diffA1, diffA2, ngdiff,
     &   temp, ngtemp,
     &   d_liquid, q0_liquid,
     &   d_solid_A, q0_solid_A,
     &   gas_constant_R,
     &   interp_type,
     &   avg_type)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngphi, ngdiff, ngtemp
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solid_A
      double precision q0_liquid, q0_solid_A
      double precision gas_constant_R
c
c variables in 3d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi))
      double precision temp(CELL3d(ifirst,ilast,ngtemp))
      double precision diffL0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffL1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffL2(SIDE3d2(ifirst,ilast,ngdiff))
      double precision diffA0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffA1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffA2(SIDE3d2(ifirst,ilast,ngdiff))
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1, ic2
      double precision vphi, hphi, invT
      double precision q0_liquid_invR, q0_solid_A_invR
      double precision diff_liquid, diff_solid_A
      double precision interp_func
      double precision average_func
c
      q0_liquid_invR = q0_liquid / gas_constant_R
      q0_solid_A_invR = q0_solid_A / gas_constant_R
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1

               vphi = average_func(
     &            phi(ic0-1,ic1,ic2), phi(ic0,ic1,ic2), avg_type )

               hphi = interp_func( vphi, interp_type )

               invT = 2.0d0 /
     &            ( temp(ic0-1,ic1,ic2) + temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solid_A = d_solid_A * exp( -q0_solid_A_invR * invT )

               diffL0(ic0,ic1,ic2) =
     &            ( 1.0d0 - hphi ) * diff_liquid
               diffA0(ic0,ic1,ic2) =
     &            hphi * diff_solid_A

            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0, ilast0

               vphi = average_func(
     &            phi(ic0,ic1-1,ic2), phi(ic0,ic1,ic2), avg_type )

               hphi = interp_func( vphi, interp_type )

               invT = 2.0d0 /
     &            ( temp(ic0,ic1-1,ic2) + temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solid_A = d_solid_A * exp( -q0_solid_A_invR * invT )

               diffL1(ic0,ic1,ic2) =
     &            ( 1.0d0 - hphi ) * diff_liquid 
               diffA1(ic0,ic1,ic2) =
     &             hphi * diff_solid_A
            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               vphi = average_func(
     &            phi(ic0,ic1,ic2-1), phi(ic0,ic1,ic2), avg_type )

               hphi = interp_func( vphi, interp_type )

               invT = 2.0d0 /
     &            ( temp(ic0,ic1,ic2-1) + temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solid_A = d_solid_A * exp( -q0_solid_A_invR * invT )

               diffL2(ic0,ic1,ic2) =
     &            ( 1.0d0 - hphi ) * diff_liquid
               diffA2(ic0,ic1,ic2) =
     &            hphi * diff_solid_A

            end do
         end do
      end do
c
      return
      end
c
      subroutine concentration_pfmdiffusion_of_temperature_threephases(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phi, nphi, ngphi,
     &   diffL0, diffL1, diffL2,
     &   diffA0, diffA1, diffA2,
     &   diffB0, diffB1, diffB2, ngdiff,
     &   temp, ngtemp,
     &   d_liquid, q0_liquid,
     &   d_solid_A, q0_solid_A,
     &   d_solid_B, q0_solid_B,
     &   gas_constant_R,
     &   interp_type,
     &   avg_type)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer nphi, ngphi, ngdiff, ngtemp, three_phases
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solid_A, d_solid_B
      double precision q0_liquid, q0_solid_A, q0_solid_B
      double precision gas_constant_R
c
c variables in 3d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi),nphi)
      double precision temp(CELL3d(ifirst,ilast,ngtemp))
      double precision diffL0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffL1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffL2(SIDE3d2(ifirst,ilast,ngdiff))
      double precision diffA0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffA1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffA2(SIDE3d2(ifirst,ilast,ngdiff))
      double precision diffB0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffB1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffB2(SIDE3d2(ifirst,ilast,ngdiff))
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1, ic2
      double precision vphi, phi0, phi1, phi2, invT
      double precision q0_liquid_invR, q0_solid_A_invR, q0_solid_B_invR
      double precision diff_liquid, diff_solid_A, diff_solid_B
      double precision interp_func
      double precision average_func
c
      q0_liquid_invR = q0_liquid / gas_constant_R
      q0_solid_A_invR = q0_solid_A / gas_constant_R
      q0_solid_B_invR = q0_solid_B / gas_constant_R
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1

               vphi = average_func(
     &            phi(ic0-1,ic1,ic2,1), phi(ic0,ic1,ic2,1), avg_type )
               phi0 = interp_func( vphi, interp_type )

               vphi = average_func(
     &            phi(ic0-1,ic1,ic2,2), phi(ic0,ic1,ic2,2), avg_type )
               phi1 = interp_func( vphi, interp_type )

               vphi = average_func(
     &            phi(ic0-1,ic1,ic2,3), phi(ic0,ic1,ic2,3), avg_type )
               phi2 = interp_func( vphi, interp_type )

               invT = 2.0d0 / ( temp(ic0-1,ic1,ic2)+temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solid_A = d_solid_A * exp( -q0_solid_A_invR * invT )
               diff_solid_B = d_solid_B * exp( -q0_solid_B_invR * invT )

               diffL0(ic0,ic1,ic2) = phi0 * diff_liquid
               diffA0(ic0,ic1,ic2) = phi1 * diff_solid_A
               diffB0(ic0,ic1,ic2) = phi2 * diff_solid_B

            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0, ilast0

               vphi = average_func(
     &            phi(ic0,ic1-1,ic2,1), phi(ic0,ic1,ic2,1), avg_type )
               phi0 = interp_func( vphi, interp_type )

               vphi = average_func(
     &             phi(ic0,ic1-1,ic2,2), phi(ic0,ic1,ic2,2), avg_type )
               phi1 = interp_func( vphi, interp_type )

               vphi = average_func(
     &            phi(ic0,ic1-1,ic2,3), phi(ic0,ic1,ic2,3), avg_type )
               phi2 = interp_func( vphi, interp_type )

               invT = 2.0d0 / ( temp(ic0,ic1-1,ic2)+temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solid_A = d_solid_A * exp( -q0_solid_A_invR * invT )
               diff_solid_B = d_solid_B * exp( -q0_solid_B_invR * invT )

               diffL1(ic0,ic1,ic2) = phi0 * diff_liquid
               diffA1(ic0,ic1,ic2) = phi1 * diff_solid_A
               diffB1(ic0,ic1,ic2) = phi2 * diff_solid_B

            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               vphi = average_func(
     &            phi(ic0,ic1,ic2-1,1), phi(ic0,ic1,ic2,1), avg_type )
               phi0 = interp_func( vphi, interp_type )

               vphi = average_func(
     &             phi(ic0,ic1,ic2-1,2), phi(ic0,ic1,ic2,2), avg_type )
               phi1 = interp_func( vphi, interp_type )

               vphi = average_func(
     &            phi(ic0,ic1,ic2-1,3), phi(ic0,ic1,ic2,3), avg_type )
               phi2 = interp_func( vphi, interp_type )

               invT = 2.0d0 / ( temp(ic0,ic1,ic2-1)+temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solid_A = d_solid_A * exp( -q0_solid_A_invR * invT )
               diff_solid_B = d_solid_B * exp( -q0_solid_B_invR * invT )

               diffL2(ic0,ic1,ic2) = phi0 * diff_liquid
               diffA2(ic0,ic1,ic2) = phi1 * diff_solid_A
               diffB2(ic0,ic1,ic2) = phi2 * diff_solid_B

            end do
         end do
      end do

      return
      end
c
c Coefficient \tilde D from Beckermann, Diepers, Steinbach, Karma, Tong, 1999
c
      subroutine concentrationdiffusion_beckermann(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phi, ngphi,
     &   diff0, diff1, diff2, ngdiff,
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
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngphi, ngdiff, ngk
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solid
c
c variables in 3d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi))
      double precision partition_coeff(CELL3d(ifirst,ilast,ngk))
      double precision diff0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diff2(SIDE3d2(ifirst,ilast,ngdiff))
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0, ic1, ic2
      double precision vphi, hphi, k
      double precision diff_liquid, diff_solid
      double precision interp_func
      double precision average_func
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1

               vphi = average_func(
     &            phi(ic0-1,ic1,ic2), phi(ic0,ic1,ic2), avg_type )
               
               hphi = interp_func( vphi, interp_type )
               
               k = partition_coeff(ic0,ic1,ic2)

               diff0(ic0,ic1,ic2) = d_solid +
     &            (d_liquid-d_solid)*( 1.0d0 - hphi )
     &                      /(1.0d0 - hphi+k*hphi)

            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0, ilast0

               vphi = average_func(
     &            phi(ic0,ic1-1,ic2), phi(ic0,ic1,ic2), avg_type )

               hphi = interp_func( vphi, interp_type )

               k = partition_coeff(ic0,ic1,ic2)

               diff1(ic0,ic1,ic2) = d_solid +
     &            (d_liquid-d_solid)*( 1.0d0 - hphi )
     &                      /(1.0d0 - hphi+k*hphi)

            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               vphi = average_func(
     &            phi(ic0,ic1,ic2-1), phi(ic0,ic1,ic2), avg_type )

               hphi = interp_func( vphi, interp_type )

               k = partition_coeff(ic0,ic1,ic2)

               diff2(ic0,ic1,ic2) = d_solid +
     &            (d_liquid-d_solid)*( 1.0d0 - hphi )
     &                      /(1.0d0 - hphi+k*hphi)

            end do
         end do
      end do
c
      return
      end

      subroutine concentration_diffcoeff_of_temperature(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   diffL0, diffL1, diffL2,
     &   diffA0, diffA1, diffA2, ngdiff,
     &   temp, ngtemp,
     &   d_liquid, q0_liquid,
     &   d_solid_A, q0_solid_A,
     &   gas_constant_R)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngdiff, ngtemp
      double precision d_liquid, d_solid_A
      double precision q0_liquid, q0_solid_A
      double precision gas_constant_R
c variables in 3d cell indexed
      double precision temp(CELL3d(ifirst,ilast,ngtemp))
      double precision diffL0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffL1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffL2(SIDE3d2(ifirst,ilast,ngdiff))
      double precision diffA0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffA1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffA2(SIDE3d2(ifirst,ilast,ngdiff))
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1, ic2
      double precision invT
      double precision q0_liquid_invR, q0_solid_A_invR
c
      q0_liquid_invR = q0_liquid / gas_constant_R
      q0_solid_A_invR = q0_solid_A / gas_constant_R
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1

               invT = 2.0d0 /
     &                ( temp(ic0-1,ic1, ic2) + temp(ic0,ic1,ic2) )

               diffL0(ic0,ic1,ic2) = d_liquid *
     &                               exp( -q0_liquid_invR * invT )
               diffA0(ic0,ic1,ic2) = d_solid_A *
     &                               exp( -q0_solid_A_invR * invT )

            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0, ilast0

               invT = 2.0d0 /
     &               ( temp(ic0,ic1-1,ic2) + temp(ic0,ic1,ic2) )

               diffL1(ic0,ic1,ic2) = d_liquid *
     &                               exp( -q0_liquid_invR * invT )
               diffA1(ic0,ic1,ic2) = d_solid_A *
     &                               exp( -q0_solid_A_invR * invT )

            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               invT = 2.0d0 /
     &                ( temp(ic0,ic1,ic2-1) + temp(ic0,ic1,ic2) )

               diffL2(ic0,ic1,ic2) = d_liquid *
     &                               exp( -q0_liquid_invR * invT )
               diffA2(ic0,ic1,ic2) = d_solid_A *
     &                               exp( -q0_solid_A_invR * invT )

            end do
         end do
      end do
c
      return
      end
