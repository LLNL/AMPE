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
      subroutine pfmdiffusion(
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
c
c Coefficient [h(phi)*d_solid+(1-h(phi))*d_liquid]
c
      subroutine pfmdiffusion_scalar(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phi, ngphi,
     &   diff0, diff1, ngdiff,
     &   d_liquid,
     &   d_solidA,
     &   interp_type,
     &   avg_type )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer ngphi, ngdiff
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solidA
c
c variables in 2d cell indexed
      double precision phi(CELL2d(ifirst,ilast,ngphi))
      double precision diff0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE2d1(ifirst,ilast,ngdiff))
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1
      double precision vphi, hphi
      double precision interp_func
      double precision average_func
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0+1
            vphi = average_func(phi(ic0-1,ic1), phi(ic0,ic1),
     &                          avg_type )
            hphi = interp_func( vphi, interp_type )
            diff0(ic0,ic1) =  (1.d0-hphi) * d_liquid
     &                     + hphi * d_solidA

         end do
      end do

      do ic1 = ifirst1, ilast1+1
         do ic0 = ifirst0, ilast0
            vphi = average_func(phi(ic0,ic1-1), phi(ic0,ic1),
     &                          avg_type )
            hphi = interp_func( vphi, interp_type )
            diff1(ic0,ic1) =  (1.d0-hphi) * d_liquid
     &                     + hphi * d_solidA

         end do
      end do
c
      return
      end
c
c Coefficient [h(phi)*d_solid+(1-h(phi))*d_liquid]
c
      subroutine add_diffusion_scalar_of_temperature(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phi, ngphi,
     &   temperature, ngtemp,
     &   diff0, diff1, ngdiff,
     &   d0_liquid, q0_liquid,
     &   d0_solid, q0_solid,
     &   gas_constant_R,
     &   interp_type)
c***********************************************************************
      implicit none
c***********************************************************************
c input:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer ngphi, ngtemp, ngdiff
      character*(*) interp_type
      double precision d0_liquid, d0_solid, q0_liquid, q0_solid
      double precision gas_constant_R
c
c variables in 2d cell indexed
      double precision phi(CELL2d(ifirst,ilast,ngphi))
      double precision temperature(CELL2d(ifirst,ilast,ngtemp))
      double precision diff0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE2d1(ifirst,ilast,ngdiff))
c
      integer ic0, ic1
      double precision hphi
      double precision interp_func
      double precision invRT, dl, ds, d
c
      do ic1 = ifirst1-1, ilast1+1
        do ic0 = ifirst0-1, ilast0+1
          invRT = 1.0d0 / (gas_constant_R*temperature(ic0,ic1))
          dl = 0.5d0*d0_liquid*invRT*exp(-q0_liquid*invRT)
          ds = 0.5d0*d0_solid*invRT*exp(-q0_solid*invRT)

          hphi = interp_func( phi(ic0,ic1), interp_type )
          d = (1.d0-hphi)*dl+hphi*ds

          diff0(ic0,ic1)   = diff0(ic0,ic1) + d
          diff0(ic0+1,ic1) = diff0(ic0+1,ic1) + d
          diff1(ic0,ic1)   = diff1(ic0,ic1) + d
          diff1(ic0,ic1+1) = diff1(ic0,ic1+1) + d

        end do
      end do
c
      return
      end
c
c Coefficient [h(phi)*d_solid+(1-h(phi))*d_liquid]
c
      subroutine pfmdiffusion_scalar_2phases(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phiL, phiA, nphiA, ngphi,
     &   diff0, diff1, ngdiff,
     &   d_liquid,
     &   d_solidA,
     &   interp_type,
     &   avg_type )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer nphiA, ngphi, ngdiff
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solidA
c
c variables in 2d cell indexed
      double precision phiL(CELL2d(ifirst,ilast,ngphi))
      double precision phiA(CELL2d(ifirst,ilast,ngphi),nphiA)
      double precision diff0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE2d1(ifirst,ilast,ngdiff))
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1, ip
      double precision vphi, hphiL, hphiA
      double precision interp_func
      double precision average_func
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0+1
            vphi = 0.d0
c assuming the first nphi-1 order parameters are solid phase
            do ip = 1, nphiA
               vphi = vphi + average_func(
     &            phiA(ic0-1,ic1,ip), phiA(ic0,ic1,ip), avg_type )
            enddo
            hphiA = interp_func( vphi, interp_type )

            vphi = average_func(phiL(ic0-1,ic1), phiL(ic0,ic1),
     &                          avg_type )
            hphiL = interp_func( vphi, interp_type )
            diff0(ic0,ic1) =  hphiL * d_liquid
     &                     + hphiA * d_solidA

         end do
      end do
c
      do ic1 = ifirst1, ilast1+1
         do ic0 = ifirst0, ilast0
            vphi = 0.d0
c assuming the first nphi-1 order parameters are solid phase
            do ip = 1, nphiA
               vphi = vphi + average_func(
     &            phiA(ic0,ic1-1,ip), phiA(ic0,ic1,ip), avg_type )
            enddo
            hphiA = interp_func( vphi, interp_type )

            vphi = average_func(phiL(ic0,ic1-1), phiL(ic0,ic1),
     &                          avg_type )
            hphiL = interp_func( vphi, interp_type )

            diff1(ic0,ic1) = hphiL * d_liquid
     &                     + hphiA * d_solidA

         end do
      end do
c
      return
      end
c
c
      subroutine pfmdiffusion_of_temperature_folchplapp(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phi, ngphi,
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
      integer ngphi, ngdiff, ngtemp
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solidA, d_solidB
      double precision q0_liquid, q0_solidA, q0_solidB
      double precision gas_constant_R
c
c variables in 2d cell indexed
      double precision phi(CELL2d(ifirst,ilast,ngphi),3)
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

            diffL0(ic0,ic1) = phi1 * diff_liquid
            diffA0(ic0,ic1) = phi2 * diff_solidA
            diffB0(ic0,ic1) = phi3 * diff_solidB

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

            diffL1(ic0,ic1) = phi1 * diff_liquid
            diffA1(ic0,ic1) = phi2 * diff_solidA
            diffB1(ic0,ic1) = phi3 * diff_solidB

         end do
      end do
c
      return
      end

c
c same as function concentrationdiffusion0, without accumulating
c component into single D
c
      subroutine pfmdiffusion_of_temperature(
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
      subroutine pfmdiffusion_of_temperature_multiorder(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phi, nphi, ngphi,
     &   diffL0, diffL1,
     &   diffA0, diffA1, ngdiff,
     &   temp, ngtemp,
     &   d_liquid, q0_liquid,
     &   d_solid, q0_solid,
     &   gas_constant_R,
     &   avg_type)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer nphi, ngphi, ngdiff, ngtemp
      character*(*) avg_type
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
      double precision average_func

      q0_liquid_invR = q0_liquid / gas_constant_R
      q0_solid_invR  = q0_solid / gas_constant_R
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0+1
            phis = 0.d0
c assuming the first nphi-1 order parameters are solid phase
            do ip = 1, nphi-1
               vphi = average_func(
     &            phi(ic0-1,ic1,ip), phi(ic0,ic1,ip), avg_type )
               vphi = max(0.d0, vphi)
               phis = phis + vphi*vphi
            enddo

            vphi = average_func(
     &         phi(ic0-1,ic1,nphi), phi(ic0,ic1,nphi), avg_type )
            vphi = max(0.d0, vphi)
            phil = vphi*vphi

            vphi = 1.d0/(phil+phis)
            phis = phis * vphi
            phil = phil * vphi

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
            phis = 0.d0
            do ip = 1, nphi-1
               vphi = average_func(
     &            phi(ic0,ic1-1,ip), phi(ic0,ic1,ip), avg_type )
               vphi = max(0.d0, vphi)
               phis = phis + vphi*vphi
            enddo

            vphi = average_func(
     &         phi(ic0,ic1-1,nphi), phi(ic0,ic1-1,nphi), avg_type )
            vphi = max(0.d0, vphi)
            phil = vphi*vphi

            vphi = 1.d0/(phil+phis)
            phis = phis * vphi
            phil = phil * vphi

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
      subroutine pfmdiffusion_of_temperature_multiorder_threephases(
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
     &   avg_type)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer nphiL, nphiA, nphiB, ngphi, ngdiff, ngtemp
      character*(*) avg_type
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

            vphiA = 0.d0
            do ip = 1, nphiA
               vphi = average_func(
     &            phiA(ic0-1,ic1,ip), phiA(ic0,ic1,ip), avg_type )
               vphi = max(0.d0, vphi)
               vphiA = vphiA + vphi*vphi
            enddo

            vphiB = 0.d0
            do ip = 1, nphiB
               vphi = average_func(
     &            phiB(ic0-1,ic1,ip), phiB(ic0,ic1,ip), avg_type )
               vphi = max(0.d0, vphi)
               vphiB = vphiB + vphi*vphi
            enddo

            vphi = average_func(
     &         phiL(ic0-1,ic1,1), phiL(ic0,ic1,1), avg_type )
            vphi = max(0.d0, vphi)
            vphiL = vphi*vphi

            vphi = 1.d0/(vphiL+vphiA+vphiB)
            vphiL = vphiL*vphi
            vphiA = vphiA*vphi
            vphiB = vphiB*vphi

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

            vphiA = 0.d0
            do ip = 1, nphiA
               vphi = average_func(
     &            phiA(ic0,ic1-1,ip), phiA(ic0,ic1,ip), avg_type )
               vphi = max(0.d0, vphi)
               vphiA = vphiA + vphi*vphi
            enddo

            vphiB = 0.d0
            do ip = 1, nphiB
               vphi = average_func(
     &            phiB(ic0,ic1-1,ip), phiB(ic0,ic1,ip), avg_type )
               vphi = max(0.d0, vphi)
               vphiB = vphiB + vphi*vphi
            enddo

            vphi = average_func(
     &         phiL(ic0,ic1-1,1), phiL(ic0,ic1,1), avg_type )
            vphi = max(0.d0, vphi)
            vphiL = vphi*vphi

            vphi = 1.d0/(vphiL+vphiA+vphiB)
            vphiL = vphiL*vphi
            vphiA = vphiA*vphi
            vphiB = vphiB*vphi

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
c
      subroutine pfmdiffusion_scalar_multiorder_3phases(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phiL, nphiL, phiA, nphiA, phiB, nphiB, ngphi,
     &   diff0, diff1, ngdiff,
     &   d_liquid,
     &   d_solidA,
     &   d_solidB,
     &   interp_type,
     &   avg_type)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer nphiL, nphiA, nphiB, ngphi, ngdiff
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solidA, d_solidB
c
c variables in 2d cell indexed
      double precision phiL(CELL2d(ifirst,ilast,ngphi),nphiL)
      double precision phiA(CELL2d(ifirst,ilast,ngphi),nphiA)
      double precision phiB(CELL2d(ifirst,ilast,ngphi),nphiB)

      double precision diff0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE2d1(ifirst,ilast,ngdiff))
c
      integer ic0, ic1, ip
      double precision vphi, vphiL, vphiA, vphiB
      double precision interp_func
      double precision average_func
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

            diff0(ic0,ic1) = vphiL * d_liquid
     &                     + vphiA * d_solidA
     &                     + vphiB * d_solidB
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

            diff1(ic0,ic1) = vphiL * d_liquid
     &                     + vphiA * d_solidA
     &                     + vphiB * d_solidB

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
     &   d0c, q0,
     &   gas_constant_R,
     &   same_phase)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer nphia, nphib
      integer ngphi, ngdiff, ngtemp, same_phase
      double precision d0c, q0
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
      integer ic0, ic1, ipa, ipb, ipbmin
      double precision pa, pb, invRT, factorT
      double precision dAB
      double precision threshold, factor
c
      threshold = 1.0d-2
c target value 1 in middle of interface where phi=0.5
      factor = 1.d0/(0.5d0-threshold)
      factor = factor**2
c
      do ic1 = ifirst1, ilast1
        do ic0 = ifirst0, ilast0
          invRT = 1.0d0 / (gas_constant_R*temp(ic0,ic1))
          factorT = d0c*invRT*exp(-q0*invRT)

          do ipa = 1, nphia
            pa =  phia(ic0,ic1,ipa)
            if( pa.gt.threshold )then
              pa = pa - threshold
              ipbmin = 1
              if( same_phase.eq.1 )then
                ipbmin = ipa+1
              endif
              do ipb = ipbmin, nphib
                pb = phib(ic0,ic1,ipb)
                if( pb.gt.threshold )then
                  pb = pb - threshold

c factor 0.5 for two contributions, one from each side
                  dAB = 0.5d0*factor*pa*pb*factorT

c add contribution to four sides of each cell
                  diffA0(ic0,ic1)   = diffA0(ic0,ic1) + dAB
                  diffA0(ic0+1,ic1) = diffA0(ic0+1,ic1) + dAB
                  diffA1(ic0,ic1)   = diffA1(ic0,ic1) + dAB
                  diffA1(ic0,ic1+1) = diffA1(ic0,ic1+1) + dAB
                  if( same_phase.eq.0 )then
                    diffB0(ic0,ic1)   = diffB0(ic0,ic1)
     &                                    + dAB
                    diffB0(ic0+1,ic1) = diffB0(ic0+1,ic1)
     &                                    + dAB
                    diffB1(ic0,ic1)   = diffB1(ic0,ic1)
     &                                    + dAB
                    diffB1(ic0,ic1+1) = diffB1(ic0,ic1+1)
     &                                    + dAB
                  endif
                endif
              end do
            endif
          end do
        end do
      end do
c
      return
      end
c
c
      subroutine add_ab_diffusion_single(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phi, ngphi,
     &   diff0, diff1, ngdiff,
     &   d0c)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer ngphi, ngdiff
      double precision d0c
c
c variables in 2d cell indexed
      double precision phi(CELL2d(ifirst,ilast,ngphi))
      double precision diff0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE2d1(ifirst,ilast,ngdiff))
c
      integer ic0, ic1
      double precision dAB, pa, pb
      double precision factor
c
c factor 0.5 for two contributions, one from each side
      factor = 0.5d0*d0c*(2.d0**4)
c
      do ic1 = ifirst1-1, ilast1+1
        do ic0 = ifirst0-1, ilast0+1
          pa =  phi(ic0,ic1)
          pb = (1.d0-pa)
          dAB = factor*pa*pa*pb*pb
c add contribution to four sides of each cell
          diff0(ic0,ic1)   = diff0(ic0,ic1) + dAB
          diff0(ic0+1,ic1) = diff0(ic0+1,ic1) + dAB
          diff1(ic0,ic1)   = diff1(ic0,ic1) + dAB
          diff1(ic0,ic1+1) = diff1(ic0,ic1+1) + dAB
        end do
      end do
c
      return
      end
c
c
      subroutine add_ab_diffusion_of_temperature_single(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phi, ngphi,
     &   temp, ngtemp,
     &   diff0, diff1, ngdiff,
     &   d0c, q0, gas_constant_R)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer ngphi, ngtemp, ngdiff
      double precision d0c, q0, gas_constant_R
c
c variables in 2d cell indexed
      double precision phi(CELL2d(ifirst,ilast,ngphi))
      double precision temp(CELL2d(ifirst,ilast,ngtemp))
      double precision diff0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE2d1(ifirst,ilast,ngdiff))
c
      integer ic0, ic1
      double precision dAB, pa, invRT
      double precision factor
c
      do ic1 = ifirst1-1, ilast1+1
        do ic0 = ifirst0-1, ilast0+1
          invRT = 1.0d0 / (gas_constant_R*temp(ic0,ic1))
c factor 0.5 for two contributions, one from each side
c factor 4 to get value a in middle of interface
          factor = 2.d0*d0c*exp(-q0*invRT)*invRT

          pa =  phi(ic0,ic1)
          if(pa .lt. 0.d0)then
            pa = 0.d0
          endif
          if(pa .gt. 1.d0)then
            pa = 1.d0
          endif
          dAB = factor*pa*(1.d0-pa)
c add contribution to four sides of each cell
          diff0(ic0,ic1)   = diff0(ic0,ic1) + dAB
          diff0(ic0+1,ic1) = diff0(ic0+1,ic1) + dAB
          diff1(ic0,ic1)   = diff1(ic0,ic1) + dAB
          diff1(ic0,ic1+1) = diff1(ic0,ic1+1) + dAB
        end do
      end do
c
      return
      end
c
c add interface diffusion to A and B diffusion
c
      subroutine add_ab_diffusion(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phia, nphia, phib, nphib, ngphi,
     &   diff0, diff1, ngdiff,
     &   d0c,
     &   same_phase)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer nphia, nphib
      integer ngphi, ngdiff, same_phase
      double precision d0c
c
c variables in 2d cell indexed
      double precision phia(CELL2d(ifirst,ilast,ngphi),nphia)
      double precision phib(CELL2d(ifirst,ilast,ngphi),nphib)
      double precision diff0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE2d1(ifirst,ilast,ngdiff))
c
      integer ic0, ic1, ipa, ipb, ipbmin
      double precision pa, pb
      double precision dAB
      double precision threshold, factor
c
      threshold = 1.0d-2
      factor = 1.d0/(0.5d0-threshold)
c factor 0.5 for two contributions, one from each side
      factor = 0.5d0*d0c*(factor**4)
c
      do ic1 = ifirst1-1, ilast1+1
        do ic0 = ifirst0-1, ilast0+1

          do ipa = 1, nphia
            pa =  phia(ic0,ic1,ipa)
            if( pa.gt.threshold )then
              pa = pa - threshold
              ipbmin = 1
              if( same_phase.eq.1 )then
                ipbmin = ipa+1
              endif
              do ipb = ipbmin, nphib
                pb = phib(ic0,ic1,ipb)
                if( pb.gt.threshold )then
                  pb = pb - threshold

                  dAB = factor*pa*pa*pb*pb

c add contribution to four sides of each cell
                  diff0(ic0,ic1)   = diff0(ic0,ic1) + dAB
                  diff0(ic0+1,ic1) = diff0(ic0+1,ic1) + dAB
                  diff1(ic0,ic1)   = diff1(ic0,ic1) + dAB
                  diff1(ic0,ic1+1) = diff1(ic0,ic1+1) + dAB
                endif
              end do
            endif
          end do
        end do
      end do
c
      return
      end
