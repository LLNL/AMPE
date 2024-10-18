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
c Coefficient [h(phi)*d_solid+(1-h(phi))*d_liquid]
c
      subroutine concentration_pfmdiffusion(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phi, ngphi,
     &   diff0, diff1, diff2, ngdiff,
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
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngphi, ngdiff, ngtemp
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solidA
      double precision q0_liquid, q0_solidA
      double precision gas_constant_R
c
c variables in 3d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi))
      double precision temp(CELL3d(ifirst,ilast,ngtemp))
      double precision diff0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diff2(SIDE3d2(ifirst,ilast,ngdiff))
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0, ic1, ic2
      double precision vphi, hphi, invT
      double precision q0_liquid_invR, q0_solidA_invR
      double precision diff_liquid, diff_solidA
      double precision interp_func
      double precision average_func
c
      q0_liquid_invR = q0_liquid / gas_constant_R
      q0_solidA_invR = q0_solidA / gas_constant_R
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
               diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )

               diff0(ic0,ic1,ic2) =
     &            ( 1.0d0 - hphi ) * diff_liquid +
     &            hphi * diff_solidA

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
               diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )

               diff1(ic0,ic1,ic2) =
     &            ( 1.0d0 - hphi ) * diff_liquid +
     &            hphi * diff_solidA

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
               diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )

               diff2(ic0,ic1,ic2) =
     &            ( 1.0d0 - hphi ) * diff_liquid +
     &            hphi * diff_solidA

            end do
         end do
      end do
c
      return
      end
c
c Coefficient [h(phi)*d_solid+(1-h(phi))*d_liquid]
c
      subroutine concentration_pfmdiffusion_scalar(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phi, ngphi,
     &   diff0, diff1, diff2, ngdiff,
     &   d_liquid,
     &   d_solidA,
     &   interp_type,
     &   avg_type )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngphi, ngdiff
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solidA
c
c variables in 2d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi))
      double precision diff0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diff2(SIDE3d2(ifirst,ilast,ngdiff))
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1, ic2
      double precision vphi, hphi
      double precision diff_liquid, diff_solidA
      double precision interp_func
      double precision average_func
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1
               vphi = average_func(phi(ic0-1,ic1,ic2), phi(ic0,ic1,ic2),
     &                          avg_type )
               hphi = interp_func( vphi, interp_type )
               diff0(ic0,ic1,ic2) =  (1.d0-hphi) * d_liquid
     &                        + hphi * d_solidA
            end do
         end do
      end do

      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0, ilast0
               vphi = average_func(phi(ic0,ic1-1,ic2), phi(ic0,ic1,ic2),
     &                             avg_type )
               hphi = interp_func( vphi, interp_type )
               diff1(ic0,ic1,ic2) =  (1.d0-hphi) * d_liquid
     &                     + hphi * d_solidA
            end do
         end do
      end do

      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0
               vphi = average_func(phi(ic0,ic1,ic2-1), phi(ic0,ic1,ic2),
     &                             avg_type )
               hphi = interp_func( vphi, interp_type )
               diff2(ic0,ic1,ic2) =  (1.d0-hphi) * d_liquid
     &                     + hphi * d_solidA
            end do
         end do
      end do
c
      return
      end
c
c Coefficient [h(phi)*d_solid+(1-h(phi))*d_liquid]
c
      subroutine concentration_pfmdiffusion_scalar_2phases(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phiL, phiA, nphiA, ngphi,
     &   diff0, diff1, diff2, ngdiff,
     &   d_liquid,
     &   d_solidA,
     &   interp_type,
     &   avg_type )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer nphiA, ngphi, ngdiff
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solidA
c
c variables in 3d cell indexed
      double precision phiL(CELL3d(ifirst,ilast,ngphi))
      double precision phiA(CELL3d(ifirst,ilast,ngphi),nphiA)
      double precision diff0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diff2(SIDE3d2(ifirst,ilast,ngdiff))
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1, ic2, ip
      double precision vphi, hphiL, hphiA
      double precision diff_liquid, diff_solidA
      double precision interp_func
      double precision average_func
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1
               vphi = 0.d0
c assuming the first nphi-1 order parameters are solid phase
               do ip = 1, nphiA
                  vphi = vphi + average_func(
     &               phiA(ic0-1,ic1,ic2,ip), phiA(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               hphiA = interp_func( vphi, interp_type )

               vphi = average_func(phiL(ic0-1,ic1,ic2),
     &                             phiL(ic0,ic1,ic2),
     &                             avg_type )
               hphiL = interp_func( vphi, interp_type )
               diff0(ic0,ic1,ic2) =  hphiL * d_liquid
     &                             + hphiA * d_solidA
            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0, ilast0
               vphi = 0.d0
c assuming the first nphi-1 order parameters are solid phase
               do ip = 1, nphiA
                  vphi = vphi + average_func(
     &               phiA(ic0,ic1-1,ic2,ip), phiA(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               hphiA = interp_func( vphi, interp_type )

               vphi = average_func(phiL(ic0,ic1-1,ic2),
     &                             phiL(ic0,ic1,ic2),
     &                             avg_type )
               hphiL = interp_func( vphi, interp_type )

               diff1(ic0,ic1,ic2) = hphiL * d_liquid
     &                        + hphiA * d_solidA
            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0
               vphi = 0.d0
c assuming the first nphi-1 order parameters are solid phase
               do ip = 1, nphiA
                  vphi = vphi + average_func(
     &               phiA(ic0,ic1,ic2-1,ip), phiA(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               hphiA = interp_func( vphi, interp_type )

               vphi = average_func(phiL(ic0,ic1,ic2-1),
     &                             phiL(ic0,ic1,ic2),
     &                             avg_type )
               hphiL = interp_func( vphi, interp_type )
               diff2(ic0,ic1,ic2) =  hphiL * d_liquid
     &                             + hphiA * d_solidA
            end do
         end do
      end do
c
      return
      end
c
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
     &   d_solidA, q0_solidA,
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
      double precision d_liquid, d_solidA
      double precision q0_liquid, q0_solidA
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
      double precision q0_liquid_invR, q0_solidA_invR
      double precision diff_liquid, diff_solidA
      double precision interp_func
      double precision average_func
c
      q0_liquid_invR = q0_liquid / gas_constant_R
      q0_solidA_invR = q0_solidA / gas_constant_R
c
      do ic2 = ifirst2-ngdiff, ilast2+ngdiff
         do ic1 = ifirst1-ngdiff, ilast1+ngdiff
            do ic0 = ifirst0, ilast0+1

               vphi = average_func(
     &            phi(ic0-1,ic1,ic2), phi(ic0,ic1,ic2), avg_type )

               hphi = interp_func( vphi, interp_type )

               invT = 2.0d0 /
     &            ( temp(ic0-1,ic1,ic2) + temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )

               diffL0(ic0,ic1,ic2) =
     &            ( 1.0d0 - hphi ) * diff_liquid
               diffA0(ic0,ic1,ic2) =
     &            hphi * diff_solidA

            end do
         end do
      end do
c
      do ic2 = ifirst2-ngdiff, ilast2+ngdiff
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0-ngdiff, ilast0+ngdiff

               vphi = average_func(
     &            phi(ic0,ic1-1,ic2), phi(ic0,ic1,ic2), avg_type )

               hphi = interp_func( vphi, interp_type )

               invT = 2.0d0 /
     &            ( temp(ic0,ic1-1,ic2) + temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )

               diffL1(ic0,ic1,ic2) =
     &            ( 1.0d0 - hphi ) * diff_liquid 
               diffA1(ic0,ic1,ic2) =
     &             hphi * diff_solidA
            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1-ngdiff, ilast1+ngdiff
            do ic0 = ifirst0-ngdiff, ilast0+ngdiff

               vphi = average_func(
     &            phi(ic0,ic1,ic2-1), phi(ic0,ic1,ic2), avg_type )

               hphi = interp_func( vphi, interp_type )

               invT = 2.0d0 /
     &            ( temp(ic0,ic1,ic2-1) + temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )

               diffL2(ic0,ic1,ic2) =
     &            ( 1.0d0 - hphi ) * diff_liquid
               diffA2(ic0,ic1,ic2) =
     &            hphi * diff_solidA

            end do
         end do
      end do
c
      return
      end
c
      subroutine concentration_pfmdiffusion_of_temperature_threephases(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phiL, nphiL, phiA, nphiA, phiB, nphiB, ngphi,
     &   diffL0, diffL1, diffL2,
     &   diffA0, diffA1, diffA2,
     &   diffB0, diffB1, diffB2, ngdiff,
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
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer nphiL, nphiA, nphiB, ngphi, ngdiff, ngtemp
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solidA, d_solidB
      double precision q0_liquid, q0_solidA, q0_solidB
      double precision gas_constant_R
c
c variables in 3d cell indexed
      double precision phiL(CELL3d(ifirst,ilast,ngphi),nphiL)
      double precision phiA(CELL3d(ifirst,ilast,ngphi),nphiA)
      double precision phiB(CELL3d(ifirst,ilast,ngphi),nphiB)
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
      integer ic0, ic1, ic2, ip
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
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1

               vphi = 0.d0
               do ip = 1, nphiA
                  vphi = vphi + average_func(
     &               phiA(ic0-1,ic1,ic2,ip), phiA(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               vphiA = interp_func( vphi, interp_type )

               vphi = 0.d0
               do ip = 1, nphiB
                  vphi = vphi + average_func(
     &               phiB(ic0-1,ic1,ic2,ip), phiB(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               vphiB = interp_func( vphi, interp_type )

               vphi = average_func(
     &            phiL(ic0-1,ic1,ic2,1), phiL(ic0,ic1,ic2,1),
     &            avg_type )
               vphiL = interp_func( vphi, interp_type )

               invT = 2.0d0 / ( temp(ic0-1,ic1,ic2)+temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )
               diff_solidB = d_solidB * exp( -q0_solidB_invR * invT )

               diffL0(ic0,ic1,ic2) = vphiL * diff_liquid
               diffA0(ic0,ic1,ic2) = vphiA * diff_solidA
               diffB0(ic0,ic1,ic2) = vphiB * diff_solidB

            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0, ilast0

               vphi = 0.d0
               do ip = 1, nphiA
                  vphi = vphi + average_func(
     &               phiA(ic0,ic1-1,ic2,ip), phiA(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               vphiA = interp_func( vphi, interp_type )

               vphi = 0.d0
               do ip = 1, nphiB
                  vphi = vphi + average_func(
     &                phiB(ic0,ic1-1,ic2,ip), phiB(ic0,ic1,ic2,ip),
     &                avg_type )
               enddo
               vphiB = interp_func( vphi, interp_type )

               vphi = average_func(
     &            phiL(ic0,ic1-1,ic2,1), phiL(ic0,ic1,ic2,1),
     &            avg_type )
               vphiL = interp_func( vphi, interp_type )

               invT = 2.0d0 / ( temp(ic0,ic1-1,ic2)+temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )
               diff_solidB = d_solidB * exp( -q0_solidB_invR * invT )

               diffL1(ic0,ic1,ic2) = vphiL * diff_liquid
               diffA1(ic0,ic1,ic2) = vphiA * diff_solidA
               diffB1(ic0,ic1,ic2) = vphiB * diff_solidB

            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               vphi = 0.d0
               do ip = 1, nphiA
                  vphi = vphi + average_func(
     &               phiA(ic0,ic1,ic2-1,ip), phiA(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               vphiA = interp_func( vphi, interp_type )

               vphi = 0.d0
               do ip = 1, nphiB
                  vphi = average_func(
     &               phiB(ic0,ic1,ic2-1,ip), phiB(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               vphiB = interp_func( vphi, interp_type )

               vphi = average_func(
     &            phiL(ic0,ic1,ic2-1,1), phiL(ic0,ic1,ic2,1),
     &            avg_type )
               vphiL = interp_func( vphi, interp_type )

               invT = 2.0d0 / ( temp(ic0,ic1,ic2-1)+temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )
               diff_solidB = d_solidB * exp( -q0_solidB_invR * invT )

               diffL2(ic0,ic1,ic2) = vphiL * diff_liquid
               diffA2(ic0,ic1,ic2) = vphiA * diff_solidA
               diffB2(ic0,ic1,ic2) = vphiB * diff_solidB

            end do
         end do
      end do

      return
      end
c
c
      subroutine concentration_pfmdiffusion_scalar_3phases(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phiL, nphiL, phiA, nphiA, phiB, nphiB, ngphi,
     &   diff0, diff1, diff2, ngdiff,
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
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer nphiL, nphiA, nphiB, ngphi, ngdiff
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solidA, d_solidB
c
c variables in 3d cell indexed
      double precision phiL(CELL3d(ifirst,ilast,ngphi),nphiL)
      double precision phiA(CELL3d(ifirst,ilast,ngphi),nphiA)
      double precision phiB(CELL3d(ifirst,ilast,ngphi),nphiB)

      double precision diff0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diff2(SIDE3d2(ifirst,ilast,ngdiff))
c
      integer ic0, ic1, ic2, ip
      double precision vphi, vphiL, vphiA, vphiB
      double precision diff_liquid, diff_solidA, diff_solidB
      double precision interp_func
      double precision average_func
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1
               vphi = 0.d0
               do ip = 1, nphiA
                  vphi = vphi + average_func(
     &               phiA(ic0-1,ic1,ic2,ip), phiA(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               vphiA = interp_func( vphi, interp_type )

               vphi = 0.d0
               do ip = 1, nphiB
                  vphi = vphi + average_func(
     &               phiB(ic0-1,ic1,ic2,ip), phiB(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               vphiB = interp_func( vphi, interp_type )

               vphi = average_func(
     &            phiL(ic0-1,ic1,ic2,1), phiL(ic0,ic1,ic2,1), avg_type )
               vphiL = interp_func( vphi, interp_type )

               diff0(ic0,ic1,ic2) = vphiL * d_liquid
     &                            + vphiA * d_solidA
     &                            + vphiB * d_solidB
            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0, ilast0

               vphi = 0.d0
               do ip = 1, nphiA
                  vphi = vphi + average_func(
     &               phiA(ic0,ic1-1,ic2,ip), phiA(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               vphiA = interp_func( vphi, interp_type )

               vphi = 0.d0
               do ip = 1, nphiB
                  vphi = vphi + average_func(
     &               phiB(ic0,ic1-1,ic2,ip), phiB(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               vphiB = interp_func( vphi, interp_type )

               vphi = average_func(
     &            phiL(ic0,ic1-1,ic2,1), phiL(ic0,ic1,ic2,1),
     &            avg_type )
               vphiL = interp_func( vphi, interp_type )

               diff1(ic0,ic1,ic2) = vphiL * d_liquid
     &                            + vphiA * d_solidA
     &                            + vphiB * d_solidB
            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               vphi = 0.d0
               do ip = 1, nphiA
                  vphi = vphi + average_func(
     &               phiA(ic0,ic1,ic2-1,ip), phiA(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               vphiA = interp_func( vphi, interp_type )

               vphi = 0.d0
               do ip = 1, nphiB
                  vphi = vphi + average_func(
     &               phiB(ic0,ic1,ic2-1,ip), phiB(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               vphiB = interp_func( vphi, interp_type )

               vphi = average_func(
     &            phiL(ic0,ic1,ic2-1,1), phiL(ic0,ic1,ic2,1),
     &            avg_type )
               vphiL = interp_func( vphi, interp_type )

               diff2(ic0,ic1,ic2) = vphiL * d_liquid
     &                            + vphiA * d_solidA
     &                            + vphiB * d_solidB
            end do
         end do
      end do
c
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
c
c
      subroutine concentration_diffcoeff_of_temperature(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   diffL0, diffL1, diffL2,
     &   diffA0, diffA1, diffA2, ngdiff,
     &   temp, ngtemp,
     &   d_liquid, q0_liquid,
     &   d_solidA, q0_solidA,
     &   gas_constant_R)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngdiff, ngtemp
      double precision d_liquid, d_solidA
      double precision q0_liquid, q0_solidA
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
      double precision q0_liquid_invR, q0_solidA_invR
c
      q0_liquid_invR = q0_liquid / gas_constant_R
      q0_solidA_invR = q0_solidA / gas_constant_R
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1

               invT = 2.0d0 /
     &                ( temp(ic0-1,ic1, ic2) + temp(ic0,ic1,ic2) )

               diffL0(ic0,ic1,ic2) = d_liquid *
     &                               exp( -q0_liquid_invR * invT )
               diffA0(ic0,ic1,ic2) = d_solidA *
     &                               exp( -q0_solidA_invR * invT )

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
               diffA1(ic0,ic1,ic2) = d_solidA *
     &                               exp( -q0_solidA_invR * invT )

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
               diffA2(ic0,ic1,ic2) = d_solidA *
     &                               exp( -q0_solidA_invR * invT )

            end do
         end do
      end do
c
      return
      end

      subroutine concentration_pfmdiffusion_of_temperature_multiphases(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phi, nphi, ngphi,
     &   diffL0, diffL1, diffL2,
     &   diffA0, diffA1, diffA2, ngdiff,
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
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer nphi, ngphi, ngdiff, ngtemp
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solid
      double precision q0_liquid, q0_solid
      double precision gas_constant_R
c
c variables in 2d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi),nphi)
      double precision temp(CELL3d(ifirst,ilast,ngtemp))
      double precision diffL0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffL1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffL2(SIDE3d2(ifirst,ilast,ngdiff))
      double precision diffA0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffA1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffA2(SIDE3d2(ifirst,ilast,ngdiff))
c
      integer ic0, ic1, ic2, ip
      double precision vphi, phil, phis, invT
      double precision q0_liquid_invR, q0_solid_invR
      double precision diff_liquid, diff_solid
      double precision interp_func
      double precision average_func

      q0_liquid_invR = q0_liquid / gas_constant_R
      q0_solid_invR  = q0_solid / gas_constant_R
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1
               vphi = 0.d0
c assuming the first nphi-1 order parameters are solid phase
               do ip = 1, nphi-1
                  vphi = vphi + average_func(
     &               phi(ic0-1,ic1,ic2,ip), phi(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               phis = interp_func( vphi, interp_type )

               vphi = average_func(
     &            phi(ic0-1,ic1,ic2,nphi), phi(ic0,ic1,ic2,nphi),
     &            avg_type )
               phil = interp_func( vphi, interp_type )

               invT = 2.0d0 / ( temp(ic0-1,ic1,ic2)+temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solid  = d_solid * exp( -q0_solid_invR * invT )

               diffL0(ic0,ic1,ic2) = phil * diff_liquid
               diffA0(ic0,ic1,ic2) = phis * diff_solid
            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0, ilast0
               vphi = 0.d0
               do ip = 1, nphi-1
                  vphi = vphi + average_func(
     &               phi(ic0,ic1-1,ic2,ip), phi(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               phis = interp_func( vphi, interp_type )

               vphi = average_func(
     &            phi(ic0,ic1-1,ic2,nphi), phi(ic0,ic1-1,ic2,nphi),
     &            avg_type )
               phil = interp_func( vphi, interp_type )

               invT = 2.0d0 / ( temp(ic0,ic1-1,ic2)+temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solid  = d_solid * exp( -q0_solid_invR * invT )

               diffL1(ic0,ic1,ic2) = phil * diff_liquid
               diffA1(ic0,ic1,ic2) = phis * diff_solid
            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0
               vphi = 0.d0
               do ip = 1, nphi-1
                  vphi = vphi + average_func(
     &               phi(ic0,ic1,ic2-1,ip), phi(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               phis = interp_func( vphi, interp_type )

               vphi = average_func(
     &            phi(ic0,ic1,ic2-1,nphi), phi(ic0,ic1-1,ic2,nphi),
     &            avg_type )
               phil = interp_func( vphi, interp_type )

               invT = 2.0d0 / ( temp(ic0,ic1,ic2-1)+temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solid  = d_solid * exp( -q0_solid_invR * invT )

               diffL2(ic0,ic1,ic2) = phil * diff_liquid
               diffA2(ic0,ic1,ic2) = phis * diff_solid
            end do
         end do
      end do

      return
      end
c
c add interface diffusion to A and B diffusion
c
      subroutine ab_diffusion_of_temperature(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phia, nphia, phib, nphib, ngphi,
     &   diffA0, diffA1, diffA2,
     &   diffB0, diffB1, diffB2, ngdiff,
     &   temp, ngtemp,
     &   d0, q0,
     &   gas_constant_R,
     &   same_phase)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer nphia, nphib, same_phase
      integer ngphi, ngdiff, ngtemp
      double precision d0, q0
      double precision gas_constant_R
c
c
c variables in 3d cell indexed
      double precision phia(CELL3d(ifirst,ilast,ngphi),nphia)
      double precision phib(CELL3d(ifirst,ilast,ngphi),nphib)
      double precision temp(CELL3d(ifirst,ilast,ngtemp))
      double precision diffA0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffA1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffA2(SIDE3d2(ifirst,ilast,ngdiff))
      double precision diffB0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffB1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffB2(SIDE3d2(ifirst,ilast,ngdiff))
c
      integer ic0, ic1, ic2, ipa, ipb, ipbmin
      double precision pa, pb, invT, factorT
      double precision q0_invR
      double precision dAB
      double precision threshold, factor
c
      threshold = 1.0d-2
      factor = 1.d0/(0.5d0-threshold)
      factor = factor**4
      q0_invR = q0 / gas_constant_R
c
      do ic2 = ifirst2, ilast2
        do ic1 = ifirst1, ilast1
          do ic0 = ifirst0, ilast0
            invT = 1.0d0 / temp(ic0,ic1,ic2)
            factorT = d0*exp(-q0_invR*invT)

            do ipa = 1, nphia
              pa =  phia(ic0,ic1,ic2,ipa)
              if( pa.gt.threshold )then
                pa = pa - threshold
                ipbmin = 1
                if( same_phase.eq.1 )then
                  ipbmin = ipa+1
                endif
                do ipb = ipbmin, nphib
                  pb = phib(ic0,ic1,ic2,ipb)
                  if( pb.gt.threshold )then
                    pb = pb - threshold

c factor 0.5 for two contributions, one from each side
                    dAB = 0.5d0*factor*pa*pa*pb*pb*factorT

c add contribution to six sides of each cell
                    diffA0(ic0,ic1,ic2)   = diffA0(ic0,ic1,ic2) + dAB
                    diffA0(ic0+1,ic1,ic2) = diffA0(ic0+1,ic1,ic2) + dAB
                    diffA1(ic0,ic1,ic2)   = diffA1(ic0,ic1,ic2) + dAB
                    diffA1(ic0,ic1+1,ic2) = diffA1(ic0,ic1+1,ic2) + dAB
                    diffA2(ic0,ic1,ic2)   = diffA2(ic0,ic1,ic2) + dAB
                    diffA2(ic0,ic1,ic2+1) = diffA2(ic0,ic1,ic2+1) + dAB
                    if( same_phase.eq.0 )then
                      diffB0(ic0,ic1,ic2)   = diffB0(ic0,ic1,ic2)
     &                                      + dAB
                      diffB0(ic0+1,ic1,ic2) = diffB0(ic0+1,ic1,ic2)
     &                                      + dAB
                      diffB1(ic0,ic1,ic2)   = diffB1(ic0,ic1,ic2)
     &                                      + dAB
                      diffB1(ic0,ic1+1,ic2) = diffB1(ic0,ic1+1,ic2)
     &                                      + dAB
                      diffB2(ic0,ic1,ic2)   = diffB2(ic0,ic1,ic2)
     &                                      + dAB
                      diffB2(ic0,ic1,ic2+1) = diffB2(ic0,ic1,ic2+1)
     &                                      + dAB
                    endif
                  endif
                end do
              endif 
            end do
          end do
        end do
      end do
c
      return
      end
c
c
      subroutine add_ab_diffusion_single(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phi, ngphi,
     &   diff0, diff1, diff2, ngdiff,
     &   d0c)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngphi, ngdiff
      double precision d0c
c
c variables in 3d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi))
      double precision diff0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diff2(SIDE3d2(ifirst,ilast,ngdiff))
c
      integer ic0, ic1, ic2
      double precision dAB, pa, pb
      double precision factor
c
c factor 0.5 for two contributions, one from each side
      factor = 0.5d0*d0c*(2.d0**4)
c
      do ic2 = ifirst2-1, ilast2+1
        do ic1 = ifirst1-1, ilast1+1
          do ic0 = ifirst0-1, ilast0+1
            pa =  phi(ic0,ic1,ic2)
            pb = (1.d0-pa)
            dAB = factor*pa*pa*pb*pb
c add contribution to four sides of each cell
            diff0(ic0,ic1,ic2)   = diff0(ic0,ic1,ic2) + dAB
            diff0(ic0+1,ic1,ic2) = diff0(ic0+1,ic1,ic2) + dAB
            diff1(ic0,ic1,ic2)   = diff1(ic0,ic1,ic2) + dAB
            diff1(ic0,ic1+1,ic2) = diff1(ic0,ic1+1,ic2) + dAB
            diff2(ic0,ic1,ic2)   = diff2(ic0,ic1,ic2) + dAB
            diff2(ic0,ic1,ic2+1) = diff2(ic0,ic1,ic2+1) + dAB
          enddo
        enddo
      enddo
c
      return
      end
c
c add interface diffusion to A and B diffusion
c
      subroutine add_ab_diffusion(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phia, nphia, phib, nphib, ngphi,
     &   diff0, diff1, diff2, ngdiff,
     &   d0c,
     &   same_phase)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer nphia, nphib
      integer ngphi, ngdiff, same_phase
      double precision d0c
c
c variables in 3d cell indexed
      double precision phia(CELL3d(ifirst,ilast,ngphi),nphia)
      double precision phib(CELL3d(ifirst,ilast,ngphi),nphib)
      double precision diff0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diff2(SIDE3d2(ifirst,ilast,ngdiff))
c
      integer ic0, ic1, ic2, ipa, ipb, ipbmin
      double precision pa, pb, factorT
      double precision dAB
      double precision threshold, factor
c
c factor 0.5 for two contributions, one from each side
      threshold = 1.0d-2
      factor = 1.d0/(0.5d0-threshold)
      factor = 0.5d0*d0c*(factor**4)
c
      do ic2 = ifirst2-1, ilast2+1
        do ic1 = ifirst1-1, ilast1+1
          do ic0 = ifirst0-1, ilast0+1
            do ipa = 1, nphia
              pa =  phia(ic0,ic1,ic2,ipa)
              if( pa.gt.threshold )then
                pa = pa - threshold
                ipbmin = 1
                if( same_phase.eq.1 )then
                  ipbmin = ipa+1
                endif
                do ipb = ipbmin, nphib
                  pb = phib(ic0,ic1,ic2,ipb)
                  if( pb.gt.threshold )then
                    pb = pb - threshold

                    dAB = factor*pa*pa*pb*pb

c add contribution to four sides of each cell
                    diff0(ic0,ic1,ic2)   = diff0(ic0,ic1,ic2) + dAB
                    diff0(ic0+1,ic1,ic2) = diff0(ic0+1,ic1,ic2) + dAB
                    diff1(ic0,ic1,ic2)   = diff1(ic0,ic1,ic2) + dAB
                    diff1(ic0,ic1+1,ic2) = diff1(ic0,ic1+1,ic2) + dAB
                    diff2(ic0,ic1,ic2)   = diff2(ic0,ic1,ic2) + dAB
                    diff2(ic0,ic1,ic2+1) = diff2(ic0,ic1,ic2+1) + dAB
                  endif
                end do
              endif
            end do
          end do
        end do
      end do
c
      return
      end
