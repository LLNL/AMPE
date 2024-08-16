c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
c 
define(NDIM,2)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
c
c Coefficient [h(phi)*d_solid+(1-h(phi))*d_liquid]
c
      subroutine concentration_pfmdiffusion(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phi, ngphi,
     &   diff0, diff1, ngdiff,
     &   temp, ngtemp,
     &   d_liquid, q0_liquid,
     &   d_solidA, q0_solidA,
     &   gas_constant_R,
     &   interp_type,
     &   avg_type )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer ngphi, ngdiff, ngtemp
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solidA
      double precision q0_liquid, q0_solidA
      double precision gas_constant_R
c
c variables in 2d cell indexed
      double precision phi(CELL2d(ifirst,ilast,ngphi))
      double precision temp(CELL2d(ifirst,ilast,ngtemp))
      double precision diff0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE2d1(ifirst,ilast,ngdiff))
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1
      double precision vphi, hphi, invT
      double precision q0_liquid_invR, q0_solidA_invR
      double precision diff_liquid, diff_solidA
      double precision interp_func
      double precision average_func
c
      q0_liquid_invR = q0_liquid / gas_constant_R
      q0_solidA_invR = q0_solidA / gas_constant_R
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0+1

            vphi = average_func(
     &         phi(ic0-1,ic1), phi(ic0,ic1), avg_type )
            hphi = interp_func( vphi, interp_type )

            invT = 2.0d0 / ( temp(ic0-1,ic1) + temp(ic0,ic1) )

            diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
            diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )

            diff0(ic0,ic1) =
     &         ( 1.0d0 - hphi ) * diff_liquid +
     &         hphi * diff_solidA

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
            diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )

            diff1(ic0,ic1) =
     &         ( 1.0d0 - hphi ) * diff_liquid +
     &         hphi * diff_solidA

         end do
      end do
c
      return
      end


      subroutine concentration_pfmdiffusion_threephases(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phi, ngphi,
     &   diff0, diff1, ngdiff,
     &   temp, ngtemp,
     &   d_liquid, q0_liquid,
     &   d_solidA, q0_solidA,
     &   d_solidB, q0_solidB,
     &   gas_constant_R,
     &   interp_type,
     &   avg_type)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer ngphi, ngdiff, ngtemp
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solidA, d_solidB
      double precision q0_liquid, q0_solidA, q0_solidB
      double precision gas_constant_R
      integer three_phase
c
c variables in 2d cell indexed
      double precision phi(CELL2d(ifirst,ilast,ngphi),3)
      double precision temp(CELL2d(ifirst,ilast,ngtemp))
      double precision diff0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE2d1(ifirst,ilast,ngdiff))
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1
      double precision vphi, phi1, phi2, phi3, invT
      double precision q0_liquid_invR, q0_solidA_invR, q0_solidB_invR
      double precision diff_liquid, diff_solidA, diff_solidB
      double precision interp_func
      double precision average_func
c
      q0_liquid_invR = q0_liquid / gas_constant_R
      q0_solidA_invR = q0_solidA / gas_constant_R
      q0_solidB_invR = q0_solidB / gas_constant_R
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0+1

            vphi = average_func(
     &         phi(ic0-1,ic1,1), phi(ic0,ic1,1), avg_type )

            phi1 = interp_func( vphi, interp_type )

            vphi = average_func(
     &         phi(ic0-1,ic1,2), phi(ic0,ic1,2), avg_type )

            phi2 = interp_func( vphi, interp_type )

            vphi = average_func(
     &         phi(ic0-1,ic1,3), phi(ic0,ic1,3), avg_type )

            phi3 = interp_func( vphi, interp_type )

            invT = 2.0d0 / ( temp(ic0-1,ic1) + temp(ic0,ic1) )

            diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
            diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )
            diff_solidB = d_solidB * exp( -q0_solidB_invR * invT )

            diff0(ic0,ic1) = phi1 * diff_liquid
     &                     + phi2 * diff_solidA
     &                     + phi3 * diff_solidB

         end do
      end do
c
      do ic1 = ifirst1, ilast1+1
         do ic0 = ifirst0, ilast0

            vphi = average_func(
     &         phi(ic0,ic1-1,1), phi(ic0,ic1,1), avg_type )
            phi1 = interp_func( vphi, interp_type )

            vphi = average_func(
     &         phi(ic0,ic1-1,2), phi(ic0,ic1,2), avg_type )
            phi2 = interp_func( vphi, interp_type )

            vphi = average_func(
     &         phi(ic0,ic1-1,3), phi(ic0,ic1,3), avg_type )
            phi3 = interp_func( vphi, interp_type )

            invT = 2.0d0 / ( temp(ic0,ic1-1) + temp(ic0,ic1) )

            diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
            diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )
            diff_solidB = d_solidB * exp( -q0_solidB_invR * invT )

            diff1(ic0,ic1) = phi1 * diff_liquid
     &                     + phi2 * diff_solidA
     &                     + phi3 * diff_solidB

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
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phi, ngphi,
     &   diffL0, diffL1,
     &   diffA0, diffA1, ngdiff,
     &   temp, ngtemp,
     &   d_liquid, q0_liquid,
     &   d_solidA, q0_solidA,
     &   gas_constant_R,
     &   interp_type,
     &   avg_type)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer ngphi, ngdiff, ngtemp
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solidA
      double precision q0_liquid, q0_solidA
      double precision gas_constant_R
c
c variables in 2d cell indexed
      double precision phi(CELL2d(ifirst,ilast,ngphi))
      double precision temp(CELL2d(ifirst,ilast,ngtemp))
      double precision diffL0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diffL1(SIDE2d1(ifirst,ilast,ngdiff))
      double precision diffA0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diffA1(SIDE2d1(ifirst,ilast,ngdiff))
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1
      double precision vphi, hphi, invT
      double precision q0_liquid_invR, q0_solidA_invR
      double precision diff_liquid, diff_solidA
      double precision interp_func
      double precision average_func
c
      q0_liquid_invR = q0_liquid / gas_constant_R
      q0_solidA_invR = q0_solidA / gas_constant_R
c
      do ic1 = ifirst1-ngdiff, ilast1+ngdiff
         do ic0 = ifirst0, ilast0+1

            vphi = average_func(
     &         phi(ic0-1,ic1), phi(ic0,ic1), avg_type )

            hphi = interp_func( vphi, interp_type )

            invT = 2.0d0 / ( temp(ic0-1,ic1) + temp(ic0,ic1) )

            diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
            diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )

            diffL0(ic0,ic1) =
     &         ( 1.0d0 - hphi ) * diff_liquid
            diffA0(ic0,ic1) =
     &         hphi * diff_solidA

         end do
      end do
c
      do ic1 = ifirst1, ilast1+1
         do ic0 = ifirst0-ngdiff, ilast0+ngdiff

            vphi = average_func(
     &         phi(ic0,ic1-1), phi(ic0,ic1), avg_type )

            hphi = interp_func( vphi, interp_type )

            invT = 2.0d0 / ( temp(ic0,ic1-1) + temp(ic0,ic1) )

            diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
            diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )

            diffL1(ic0,ic1) =
     &         ( 1.0d0 - hphi ) * diff_liquid
            diffA1(ic0,ic1) =
     &         hphi * diff_solidA
         end do
      end do
c
      return
      end
c
c
      subroutine concentration_pfmdiffusion_of_temperature_multiphases(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phi, nphi, ngphi,
     &   diffL0, diffL1,
     &   diffA0, diffA1, ngdiff,
     &   temp, ngtemp,
     &   d_liquid, q0_liquid,
     &   d_solid, q0_solid,
     &   gas_constant_R,
     &   interp_type,
     &   avg_type)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer nphi, ngphi, ngdiff, ngtemp
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solid
      double precision q0_liquid, q0_solid
      double precision gas_constant_R
c
c variables in 2d cell indexed
      double precision phi(CELL2d(ifirst,ilast,ngphi),nphi)
      double precision temp(CELL2d(ifirst,ilast,ngtemp))
      double precision diffL0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diffL1(SIDE2d1(ifirst,ilast,ngdiff))
      double precision diffA0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diffA1(SIDE2d1(ifirst,ilast,ngdiff))
c
      integer ic0, ic1, ip
      double precision vphi, phil, phis, invT
      double precision q0_liquid_invR, q0_solid_invR
      double precision diff_liquid, diff_solid
      double precision interp_func
      double precision average_func

      q0_liquid_invR = q0_liquid / gas_constant_R
      q0_solid_invR  = q0_solid / gas_constant_R
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0+1
            vphi = 0.d0
c assuming the first nphi-1 order parameters are solid phase
            do ip = 1, nphi-1
               vphi = vphi + average_func(
     &            phi(ic0-1,ic1,ip), phi(ic0,ic1,ip), avg_type )
            enddo
            phis = interp_func( vphi, interp_type )

            vphi = average_func(
     &         phi(ic0-1,ic1,nphi), phi(ic0,ic1,nphi), avg_type )
            phil = interp_func( vphi, interp_type )

            invT = 2.0d0 / ( temp(ic0-1,ic1) + temp(ic0,ic1) )

            diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
            diff_solid  = d_solid * exp( -q0_solid_invR * invT )

            diffL0(ic0,ic1) = phil * diff_liquid
            diffA0(ic0,ic1) = phis * diff_solid
         end do
      end do
c
      do ic1 = ifirst1, ilast1+1
         do ic0 = ifirst0, ilast0
            vphi = 0.d0
            do ip = 1, nphi-1
               vphi = vphi + average_func(
     &            phi(ic0,ic1-1,ip), phi(ic0,ic1,ip), avg_type )
            enddo
            phis = interp_func( vphi, interp_type )

            vphi = average_func(
     &         phi(ic0,ic1-1,nphi), phi(ic0,ic1-1,nphi), avg_type )
            phil = interp_func( vphi, interp_type )

            invT = 2.0d0 / ( temp(ic0,ic1-1) + temp(ic0,ic1) )

            diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
            diff_solid  = d_solid * exp( -q0_solid_invR * invT )

            diffL1(ic0,ic1) = phil * diff_liquid
            diffA1(ic0,ic1) = phis * diff_solid
         end do
      end do
c
      return
      end

c
c
      subroutine concentration_pfmdiffusion_of_temperature_threephases(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phiL, nphiL, phiA, nphiA, phiB, nphiB, ngphi,
     &   diffL0, diffL1,
     &   diffA0, diffA1,
     &   diffB0, diffB1, ngdiff,
     &   temp, ngtemp,
     &   d_liquid, q0_liquid,
     &   d_solidA, q0_solidA,
     &   d_solidB, q0_solidB,
     &   gas_constant_R,
     &   interp_type,
     &   avg_type)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer nphiL, nphiA, nphiB, ngphi, ngdiff, ngtemp
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solidA, d_solidB
      double precision q0_liquid, q0_solidA, q0_solidB
      double precision gas_constant_R
c
c variables in 2d cell indexed
      double precision phiL(CELL2d(ifirst,ilast,ngphi),nphiL)
      double precision phiA(CELL2d(ifirst,ilast,ngphi),nphiA)
      double precision phiB(CELL2d(ifirst,ilast,ngphi),nphiB)

      double precision temp(CELL2d(ifirst,ilast,ngtemp))
      double precision diffL0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diffL1(SIDE2d1(ifirst,ilast,ngdiff))
      double precision diffA0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diffA1(SIDE2d1(ifirst,ilast,ngdiff))
      double precision diffB0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diffB1(SIDE2d1(ifirst,ilast,ngdiff))
c
      integer ic0, ic1, ip
      double precision vphi, vphiL, vphiA, vphiB, invT
      double precision q0_liquid_invR, q0_solidA_invR, q0_solidB_invR
      double precision diff_liquid, diff_solidA, diff_solidB
      double precision interp_func
      double precision average_func
c
      q0_liquid_invR = q0_liquid / gas_constant_R
      q0_solidA_invR = q0_solidA / gas_constant_R
      q0_solidB_invR = q0_solidB / gas_constant_R
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0+1

            vphi = 0.d0
            do ip = 1, nphiA
               vphi = vphi + average_func(
     &            phiA(ic0-1,ic1,ip), phiA(ic0,ic1,ip), avg_type )
            enddo
            vphiA = interp_func( vphi, interp_type )

            vphi = 0.d0
            do ip = 1, nphiB
               vphi = vphi + average_func(
     &            phiB(ic0-1,ic1,ip), phiB(ic0,ic1,ip), avg_type )
            enddo
            vphiB = interp_func( vphi, interp_type )

            vphi = average_func(
     &         phiL(ic0-1,ic1,1), phiL(ic0,ic1,1), avg_type )
            vphiL = interp_func( vphi, interp_type )

            invT = 2.0d0 / ( temp(ic0-1,ic1) + temp(ic0,ic1) )

            diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
            diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )
            diff_solidB = d_solidB * exp( -q0_solidB_invR * invT )

            diffL0(ic0,ic1) = vphiL * diff_liquid
            diffA0(ic0,ic1) = vphiA * diff_solidA
            diffB0(ic0,ic1) = vphiB * diff_solidB
         end do
      end do
c
      do ic1 = ifirst1, ilast1+1
         do ic0 = ifirst0, ilast0

            vphi = 0.d0
            do ip = 1, nphiA
               vphi = vphi + average_func(
     &            phiA(ic0,ic1-1,ip), phiA(ic0,ic1,ip), avg_type )
            enddo
            vphiA = interp_func( vphi, interp_type )

            vphi = 0.d0
            do ip = 1, nphiB
               vphi = vphi + average_func(
     &            phiB(ic0,ic1-1,ip), phiB(ic0,ic1,ip), avg_type )
            enddo
            vphiB = interp_func( vphi, interp_type )

            vphi = average_func(
     &         phiL(ic0,ic1-1,1), phiL(ic0,ic1,1), avg_type )
            vphiL = interp_func( vphi, interp_type )

            invT = 2.0d0 / ( temp(ic0,ic1-1) + temp(ic0,ic1) )

            diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
            diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )
            diff_solidB = d_solidB * exp( -q0_solidB_invR * invT )

            diffL1(ic0,ic1) = vphiL * diff_liquid
            diffA1(ic0,ic1) = vphiA * diff_solidA
            diffB1(ic0,ic1) = vphiB * diff_solidB

         end do
      end do
c
      return
      end

c
c Coefficient \tilde D from Beckermann, Diepers, Steinbach, Karma, Tong, 1999
c
      subroutine concentrationdiffusion_beckermann(
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

      subroutine concentration_diffcoeff_of_temperature(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   diffL0, diffL1,
     &   diffA0, diffA1,
     &   diffB0, diffB1, ngdiff,
     &   temp, ngtemp,
     &   d_liquid, q0_liquid,
     &   d_solidA, q0_solidA,
     &   d_solidB, q0_solidB,
     &   gas_constant_R, three_phases)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer ngdiff, ngtemp, three_phases
      double precision d_liquid, d_solidA, d_solidB
      double precision q0_liquid, q0_solidA, q0_solidB
      double precision gas_constant_R
c
c variables in 2d cell indexed
c
c variables in 2d cell indexed
      double precision temp(CELL2d(ifirst,ilast,ngtemp))
      double precision diffL0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diffL1(SIDE2d1(ifirst,ilast,ngdiff))
      double precision diffA0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diffA1(SIDE2d1(ifirst,ilast,ngdiff))
      double precision diffB0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diffB1(SIDE2d1(ifirst,ilast,ngdiff))
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1
      double precision invT
      double precision q0_liquid_invR, q0_solidA_invR, q0_solidB_invR
c
      q0_liquid_invR = q0_liquid / gas_constant_R
      q0_solidA_invR = q0_solidA / gas_constant_R
      q0_solidB_invR = q0_solidB / gas_constant_R
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0+1

            invT = 2.0d0 / ( temp(ic0-1,ic1) + temp(ic0,ic1) )

            diffL0(ic0,ic1) = d_liquid * exp( -q0_liquid_invR * invT )
            diffA0(ic0,ic1) = d_solidA * exp( -q0_solidA_invR * invT )
            if( three_phases /= 0 ) then
               diffB0(ic0,ic1) =
     &           d_solidB * exp( -q0_solidB_invR * invT )
            endif
         end do
      end do
c
      do ic1 = ifirst1, ilast1+1
         do ic0 = ifirst0, ilast0

            invT = 2.0d0 / ( temp(ic0,ic1-1) + temp(ic0,ic1) )

            diffL1(ic0,ic1) = d_liquid * exp( -q0_liquid_invR * invT )
            diffA1(ic0,ic1) = d_solidA * exp( -q0_solidA_invR * invT )
            if( three_phases /= 0 ) then
               diffB1(ic0,ic1) =
     &            d_solidB * exp( -q0_solidB_invR * invT )
            endif
         end do
      end do
c
      return
      end

c
c add interface diffusion to A and B diffusion
c
      subroutine ab_diffusion_of_temperature(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phia, nphia, phib, nphib, ngphi,
     &   diffA0, diffA1,
     &   diffB0, diffB1, ngdiff,
     &   temp, ngtemp,
     &   d0, q0,
     &   gas_constant_R,
     &   avg_type, dupl)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer nphia, nphib
      integer ngphi, ngdiff, ngtemp, dupl
      character*(*) avg_type
      double precision d0, q0
      double precision gas_constant_R
c
c variables in 2d cell indexed
      double precision phia(CELL2d(ifirst,ilast,ngphi),nphia)
      double precision phib(CELL2d(ifirst,ilast,ngphi),nphib)
      double precision temp(CELL2d(ifirst,ilast,ngtemp))
      double precision diffA0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diffA1(SIDE2d1(ifirst,ilast,ngdiff))
      double precision diffB0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diffB1(SIDE2d1(ifirst,ilast,ngdiff))
c
      integer ic0, ic1, ipa, ipb
      double precision pa, pb, invT, factorT
      double precision q0_invR
      double precision dAB
      double precision average_func
c
      q0_invR = q0 / gas_constant_R
c
      do ipa = 1, nphia
        do ipb = 1, nphib
          if((dupl.eq.0) .or. (ipa.ne.ipb))then
            do ic1 = ifirst1, ilast1
              do ic0 = ifirst0, ilast0+1
                invT = 2.0d0 / ( temp(ic0-1,ic1) + temp(ic0,ic1) )
                factorT = d0*exp(-q0_invR*invT)

                pa = average_func(
     &            phia(ic0-1,ic1,ipa), phia(ic0,ic1,ipa), avg_type )

                pb = average_func(
     &            phib(ic0-1,ic1,ipb), phib(ic0,ic1,ipb), avg_type )

                dAB = 16.d0*pa*pa*pb*pb*factorT

                diffA0(ic0,ic1) = diffA0(ic0,ic1) + dAB
                diffB0(ic0,ic1) = diffB0(ic0,ic1) + dAB
              enddo
            enddo
            do ic1 = ifirst1, ilast1+1
              do ic0 = ifirst0, ilast0
                invT = 2.0d0 / ( temp(ic0,ic1-1) + temp(ic0,ic1) )
                factorT = d0*exp(-q0_invR*invT)

                pa = average_func(
     &            phia(ic0,ic1-1,ipa), phia(ic0,ic1,ipa), avg_type )

                pb = average_func(
     &            phib(ic0,ic1-1,ipb), phib(ic0,ic1,ipb), avg_type )

                dAB = 16.d0*pa*pa*pb*pb*factorT

                diffA1(ic0,ic1) = diffA1(ic0,ic1) + dAB
                diffB1(ic0,ic1) = diffB1(ic0,ic1) + dAB
              end do
            end do
          endif
        end do
      end do
c
      return
      end
