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

      subroutine gradient_flux(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   h, epsilon,
     &   phase, ngphase,
     &   flux0, flux1, ngflux)

      implicit none
      integer ifirst0, ilast0, ifirst1, ilast1,
     &     ngphase, ngflux
      double precision 
     &     phase(CELL2d(ifirst,ilast,ngphase)),
     &     flux0(SIDE2d0(ifirst,ilast,ngflux)),
     &     flux1(SIDE2d1(ifirst,ilast,ngflux)),
     &     h(2)

c     local variables
      integer i, j
      double precision dxinv, dyinv, epsilon

      double precision epsilon2
c
c      print*,'gradient_flux'
c
      epsilon2 = epsilon * epsilon

      dxinv = epsilon2 / h(1)
      dyinv = epsilon2 / h(2)
      
      do j = ifirst1, ilast1
         do i = ifirst0, ilast0+1
            flux0(i,j) = (phase(i,j) - phase(i-1,j)) * dxinv
         enddo
      enddo

      do j = ifirst1, ilast1+1
         do i = ifirst0, ilast0
            flux1(i,j) = (phase(i,j) - phase(i,j-1)) * dyinv
         enddo
      enddo

      return
      end
      
c***********************************************************************
      subroutine gradient_flux_wide(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   h, epsilon,
     &   phase, ngphase,
     &   flux0, flux1, ngflux)

      implicit none
      integer ifirst0, ilast0, ifirst1, ilast1,
     &     ngphase, ngflux
      double precision 
     &     phase(CELL2d(ifirst,ilast,ngphase)),
     &     flux0(SIDE2d0(ifirst,ilast,ngflux)),
     &     flux1(SIDE2d1(ifirst,ilast,ngflux)),
     &     h(2)

c     local variables
      integer i, j
      double precision dxinv, dyinv, epsilon

      double precision epsilon2
c
      epsilon2 = epsilon * epsilon

      dxinv = epsilon2 / h(1)
      dyinv = epsilon2 / h(2)
      
      do j = ifirst1, ilast1
         do i = ifirst0, ilast0+1
            flux0(i,j) = 0.125d0*dxinv*( 
     &           (phase(i,j-1) - phase(i-1,j-1)) 
     &         + (phase(i,j)   - phase(i-1,j))*6.d0 
     &         + (phase(i,j+1) - phase(i-1,j+1)) )
         enddo
      enddo

      do j = ifirst1, ilast1+1
         do i = ifirst0, ilast0
            flux1(i,j) = 0.125d0*dyinv*(
     &           (phase(i-1,j) - phase(i-1,j-1)) 
     &         + (phase(i,  j) - phase(i,  j-1))*6.d0 
     &         + (phase(i+1,j) - phase(i+1,j-1)) ) 
         enddo
      enddo

      return
      end
      
c***********************************************************************
c see Shukla and Giri, J. Comput. Phys. 276 (2014), p.259
c
      subroutine compute_flux_isotropic(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   h, epsilon,
     &   phase, ngphase,
     &   flux0, flux1, ngflux)

      implicit none
      integer ifirst0, ilast0, ifirst1, ilast1,
     &     ngphase, ngflux
      double precision 
     &     phase(CELL2d(ifirst,ilast,ngphase)),
     &     flux0(SIDE2d0(ifirst,ilast,ngflux)),
     &     flux1(SIDE2d1(ifirst,ilast,ngflux)),
     &     h(2)

c     local variables
      integer i, j
      double precision dxinv, dyinv, epsilon

      double precision epsilon2
c
      epsilon2 = epsilon * epsilon

      dxinv = (1.d0/12.d0) * epsilon2 / h(1)
      dyinv = (1.d0/12.d0) * epsilon2 / h(2)
      
      do j = ifirst1, ilast1
         do i = ifirst0, ilast0+1
            flux0(i,j) = dxinv*( 
     &           (phase(i,j-1) - phase(i-1,j-1)) 
     &         + (phase(i,j)   - phase(i-1,j))*10.d0 
     &         + (phase(i,j+1) - phase(i-1,j+1)) )
         enddo
      enddo

      do j = ifirst1, ilast1+1
         do i = ifirst0, ilast0
            flux1(i,j) = dyinv*(
     &           (phase(i-1,j) - phase(i-1,j-1)) 
     &         + (phase(i,  j) - phase(i,  j-1))*10.d0 
     &         + (phase(i+1,j) - phase(i+1,j-1)) ) 
         enddo
      enddo

      return
      end
      
c***********************************************************************
      subroutine anisotropic_gradient_flux(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   h, epsilon, nu, knumber,
     &   phase, ngphase,
     &   quat, ngq, qlen,
     &   flux0, flux1, ngflux)

      implicit none
      integer ifirst0, ilast0, ifirst1, ilast1,
     &     ngphase, ngq, qlen, ngflux, knumber
      double precision 
     &     phase(CELL2d(ifirst,ilast,ngphase)),
     &     quat(CELL2d(ifirst,ilast,ngq),qlen),
     &     flux0(SIDE2d0(ifirst,ilast,ngflux)),
     &     flux1(SIDE2d1(ifirst,ilast,ngflux)),
     &     nu, h(2)

c     local variables
      integer i, j
      double precision dxinv, dyinv, epsilon
      double precision theta, phi
      double precision pi, q
      double precision epsilon2, dphidx, dphidy
      double precision epstheta, depsdtheta
c
      pi = 4.d0*atan(1.d0)
c
      epsilon2 = epsilon * epsilon

      dxinv = 1.d0 / h(1)
      dyinv = 1.d0 / h(2)

c x faces      
      do j = ifirst1, ilast1
         do i = ifirst0, ilast0+1
            dphidx = (phase(i,j) - phase(i-1,j)) * dxinv
            dphidy = 0.25d0*(phase(i-1,j+1) - phase(i-1,j-1)
     &                   + phase(i  ,j+1) - phase(i  ,j-1)) 
     &                   * dyinv 
            if( abs(dphidx)>1.e-12 )then
               theta=atan(dphidy/dphidx)
            else
               theta=0.5d0*pi
            endif

            q=0.5d0*(quat(i-1,j,1)+quat(i,j,1))
c q could be slightly out of the [-1,+1] range due to finite precision
c need to enforce q to be within that range before call to acos
            if( q .gt.  1.d0 )q=1.d0
            if( q .lt. -1.d0 )q=-1.d0

            if( qlen==4 )then
               phi=2.d0*acos(q)
            else
               phi=acos(q)
            endif
            
            epstheta=epsilon*(1.d0+nu*cos(knumber*(theta-phi)))
            depsdtheta=-knumber*epsilon*nu*sin(knumber*(theta-phi))
            
            flux0(i,j) = epstheta*epstheta*dphidx 
     &                 - epstheta*depsdtheta*dphidy
c            if( flux0(i,j) /= flux0(i,j) )then
c               print*,'flux0,',i,j,q,phi,theta,epstheta,
c     &                          depsdtheta,flux0(i,j)
c            endif
         enddo
      enddo

c y faces      
      do j = ifirst1, ilast1+1
         do i = ifirst0, ilast0
            dphidx = 0.25*(phase(i+1,j-1) - phase(i-1,j-1)
     &                   + phase(i+1,j  ) - phase(i-1,j  )) 
     &                   * dxinv 
            dphidy = (phase(i,j) - phase(i,j-1)) * dyinv
         
            if( abs(dphidx)>1.e-12 )then
               theta=atan(dphidy/dphidx)
            else
               theta=0.5*pi
            endif
            
            q=0.5*(quat(i,j-1,1)+quat(i,j,1))
            if( q .gt.  1.d0 )q=1.d0
            if( q .lt. -1.d0 )q=-1.d0

            if( qlen==4 )then
               phi=2.*acos(q)
            else
               phi=acos(q)
            endif
            
            epstheta=epsilon*(1.+nu*cos(knumber*(theta-phi)))
            depsdtheta=-knumber*epsilon*nu*sin(knumber*(theta-phi))
            
            flux1(i,j) = epstheta*epstheta*dphidy 
     &                 + epstheta*depsdtheta*dphidx
         enddo
      enddo

      return
      end
      
c***********************************************************************
c
c compute r.h.s. for phase variable phi (Pusztai et al. model)
c
c   5 point stencil
c
      subroutine computerhspbg(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   dx, 
     &   misorientation_factor,
     &   epsilonq,
     &   flux0,
     &   flux1,
     &   ngflux,
     &   temp, ngtemp,
     &   phi_well_scale,
     &   eta_well_scale,
     &   phi, ngphi,
     &   eta, ngeta,
     &   orient_grad_mod, ngogm,
     &   rhs, ngrhs,
     &   phi_well_type,
     &   eta_well_type,
     &   energy_interp_type,
     &   orient_interp_type1,
     &   orient_interp_type2,
     &   with_orient,
     &   three_phase )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer with_orient, three_phase

      double precision dx(2)
      double precision misorientation_factor, epsilonq
      double precision phi_well_scale, eta_well_scale
      integer ngflux, ngphi, ngeta, ngogm, ngrhs, ngtemp
      character*(*) phi_well_type
      character*(*) eta_well_type
      character*(*) energy_interp_type
      character*(*) orient_interp_type1, orient_interp_type2
c
c variables in 2d cell indexed
      double precision phi(CELL2d(ifirst,ilast,ngphi))
      double precision eta(CELL2d(ifirst,ilast,ngeta))
      double precision orient_grad_mod(CELL2d(ifirst,ilast,ngogm))
      double precision rhs(CELL2d(ifirst,ilast,ngrhs))
      double precision temp(CELL2d(ifirst,ilast,ngtemp))

c variables in 2d side indexed
      double precision 
     &     flux0(SIDE2d0(ifirst,ilast,ngflux)),
     &     flux1(SIDE2d1(ifirst,ilast,ngflux))
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0, ic1
      double precision diff_term_x, diff_term_y, diff_term

      double precision g, g_prime, h_prime, p1_prime, p2_prime
      double precision deriv_interp_func
      double precision well_func
      double precision deriv_well_func
      double precision dxinv, dyinv
      double precision epsilonq2

      dxinv = 1.d0 / dx(1)
      dyinv = 1.d0 / dx(2)
      epsilonq2 = 0.5d0*epsilonq*epsilonq
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0

            diff_term_x = 
     &           (flux0(ic0+1,ic1) - flux0(ic0,ic1)) * dxinv
            diff_term_y = 
     &           (flux1(ic0,ic1+1) - flux1(ic0,ic1)) * dyinv

            diff_term = diff_term_x + diff_term_y

            rhs(ic0,ic1) = diff_term

c  Phase energy well

            g_prime =
     &         deriv_well_func( phi(ic0,ic1), phi_well_type )

            rhs(ic0,ic1) = rhs(ic0,ic1) -
     &         phi_well_scale * g_prime

         enddo
      enddo

c  Eta energy well

      if ( three_phase /= 0 ) then

         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               h_prime =
     &            deriv_interp_func( phi(ic0,ic1), energy_interp_type )
               
               g = well_func( eta(ic0,ic1), eta_well_type )

               rhs(ic0,ic1) = rhs(ic0,ic1) -
     &            eta_well_scale * h_prime * g

            enddo
         enddo

      endif

c  Misorientation energy
c
      if ( with_orient /= 0 ) then

         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               p1_prime =
     &            deriv_interp_func(
     &               phi(ic0,ic1),
     &               orient_interp_type1 )
               p2_prime =
     &            deriv_interp_func(
     &               phi(ic0,ic1),
     &               orient_interp_type2 )

               rhs(ic0,ic1) = rhs(ic0,ic1) -
     &            misorientation_factor *
     &            temp(ic0,ic1) *
     &            p1_prime *
     &            orient_grad_mod(ic0,ic1)
     &          - p2_prime * epsilonq2 *
     &            orient_grad_mod(ic0,ic1) *orient_grad_mod(ic0,ic1)

            enddo
         enddo

      endif

      return
      end

c***********************************************************************
c
c compute r.h.s. for 3 phases
c
      subroutine computerhsthreephases(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   dx,
     &   flux0,
     &   flux1,
     &   ngflux,
     &   phi_well_scale,
     &   phi, ngphi,
     &   rhs, ngrhs)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1

      double precision dx(2)
      double precision phi_well_scale
      integer ngflux, ngphi, ngrhs
c
c variables in 2d cell indexed
      double precision phi(CELL2d(ifirst,ilast,ngphi),3)
      double precision rhs(CELL2d(ifirst,ilast,ngrhs),3)

c variables in 2d side indexed
      double precision
     &     flux0(SIDE2d0(ifirst,ilast,ngflux),3),
     &     flux1(SIDE2d1(ifirst,ilast,ngflux),3)
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1
      integer ip, jp, ip1, ip2
      double precision diff_term_x, diff_term_y

      double precision g_prime
      double precision deriv_triple_well_func
      double precision dxinv, dyinv, fac, dsum
c
      dxinv = 1.d0 / dx(1)
      dyinv = 1.d0 / dx(2)
c
      do ip = 1, 3
         ip1 = MOD(ip,3)+1
         ip2 = MOD(ip+1,3)+1
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               diff_term_x =
     &              (flux0(ic0+1,ic1,ip) - flux0(ic0,ic1,ip)) * dxinv
               diff_term_y =
     &              (flux1(ic0,ic1+1,ip) - flux1(ic0,ic1,ip)) * dyinv

               rhs(ic0,ic1,ip) = diff_term_x + diff_term_y

c  Phase energy well

               g_prime =
     &            deriv_triple_well_func(
     &               phi(ic0,ic1,ip), phi(ic0,ic1,ip1),
     &               phi(ic0,ic1,ip2) )

               rhs(ic0,ic1,ip) = rhs(ic0,ic1,ip) -
     &            phi_well_scale * g_prime

            enddo
         enddo
      enddo

c enforce constraint sum components = 0
      fac = 1.d0/3.d0
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0
            dsum = 0.d0
            do jp = 1, 3
              dsum = dsum +  rhs(ic0,ic1,jp)
            enddo
            dsum = dsum * fac
            do jp = 1, 3
               rhs(ic0,ic1,jp) = rhs(ic0,ic1,jp) - dsum
            enddo
         enddo
      enddo

      return
      end

c***********************************************************************
c
c compute r.h.s. for eta variable phi
c
c   5 point stencil
c
      subroutine computerhseta(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   dx, 
     &   flux0,
     &   flux1,
     &   ngflux,
     &   temp, ngtemp,
     &   eta_well_scale,
     &   phi, ngphi,
     &   eta, ngeta,
     &   rhs, ngrhs,
     &   eta_well_type,
     &   phi_interp_type )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1

      double precision dx(0:1)
      double precision eta_well_scale
      integer ngflux, ngphi, ngeta, ngrhs, ngtemp
      character*(*) eta_well_type
      character*(*) phi_interp_type
c
c variables in 2d cell indexed
      double precision phi(CELL2d(ifirst,ilast,ngphi))
      double precision eta(CELL2d(ifirst,ilast,ngeta))
      double precision rhs(CELL2d(ifirst,ilast,ngrhs))
      double precision temp(CELL2d(ifirst,ilast,ngtemp))

c variables in 2d side indexed
      double precision 
     &     flux0(SIDE2d0(ifirst,ilast,ngflux)),
     &     flux1(SIDE2d1(ifirst,ilast,ngflux))
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0, ic1
      double precision diff_term_x, diff_term_y, diff_term
      double precision dxinv, dyinv

      double precision g_prime, h
      double precision interp_func
      double precision deriv_well_func
c
      dxinv = 1.d0 / dx(0)
      dyinv = 1.d0 / dx(1)
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0

            diff_term_x = 
     &           ( flux0(ic0+1,ic1) - flux0(ic0,ic1) )*dxinv
            diff_term_y = 
     &           ( flux1(ic0,ic1+1) - flux1(ic0,ic1) )*dyinv

            diff_term = diff_term_x + diff_term_y

            rhs(ic0,ic1) = diff_term

c  Energy well

            h = interp_func( phi(ic0,ic1), phi_interp_type )

            g_prime =
     &         deriv_well_func( eta(ic0,ic1), eta_well_type )

            rhs(ic0,ic1) = rhs(ic0,ic1) -
     &         eta_well_scale * h * g_prime

         enddo
      enddo

      return
      end

c***********************************************************************
c
c compute r.h.s. component due to free energy for phase variable phi
c
      subroutine phaserhs_fenergy(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   fl, fa,
     &   phi, ngphi,
     &   rhs, ngrhs,
     &   energy_interp_type)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1

      integer ngphi, ngrhs
      character*(*) energy_interp_type
c
c variables in 2d cell indexed
      double precision fl(CELL2d(ifirst,ilast,0))
      double precision fa(CELL2d(ifirst,ilast,0))
      double precision phi(CELL2d(ifirst,ilast,ngphi))
      double precision rhs(CELL2d(ifirst,ilast,ngrhs))
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0, ic1

      double precision hphi_prime
      double precision deriv_interp_func
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0

            hphi_prime = 
     &         deriv_interp_func( phi(ic0,ic1), energy_interp_type )

            rhs(ic0,ic1) = rhs(ic0,ic1) +
     &         hphi_prime *
     &         ( fl(ic0,ic1) - fa(ic0,ic1) )

         enddo
      enddo

      return
      end

c
c compute r.h.s. component due to free energy for phase variable phi
c
      subroutine phaserhs_fenergy_multiorderp(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   fl, fa,
     &   phi, ngphi, nphi,
     &   rhs, ngrhs)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1

      integer ngphi, ngrhs, nphi
c
c variables in 2d cell indexed
      double precision fl(CELL2d(ifirst,ilast,0))
      double precision fa(CELL2d(ifirst,ilast,0))
      double precision phi(CELL2d(ifirst,ilast,ngphi),nphi)
      double precision rhs(CELL2d(ifirst,ilast,ngrhs),nphi)
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1, ip

      double precision hphis, hphil, suminv2
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0
            hphis = 0.d0
            do ip = 1, nphi-1
               hphis = hphis + phi(ic0,ic1,ip)*phi(ic0,ic1,ip)
            enddo
            hphil = phi(ic0,ic1,nphi)*phi(ic0,ic1,nphi)
            suminv2 = 1.d0/(hphis+hphil)
            hphis = hphis*suminv2
            hphil = hphil*suminv2

            do ip = 1, nphi-1
               rhs(ic0,ic1,ip) = rhs(ic0,ic1,ip) +
     &            2.d0*phi(ic0,ic1,ip)*hphil*suminv2
     &            *( fl(ic0,ic1) - fa(ic0,ic1) )
            enddo
            rhs(ic0,ic1,nphi) = rhs(ic0,ic1,nphi) +
     &            2.d0*phi(ic0,ic1,nphi)*hphis*suminv2
     &            *( fa(ic0,ic1) - fl(ic0,ic1) )

         enddo
      enddo

      return
      end

c***********************************************************************
c
c compute r.h.s. component due to free energy for phase variable phi
c
      subroutine etarhs_fenergy(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   fl, fa, fb,
     &   phi, ngphi,
     &   eta, ngeta,
     &   rhs, ngrhs,
     &   phi_interp_type,
     &   eta_interp_type )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1

      integer ngphi, ngeta, ngrhs
      character*(*) phi_interp_type
      character*(*) eta_interp_type
c
c variables in 2d cell indexed
      double precision phi(CELL2d(ifirst,ilast,ngphi))
      double precision eta(CELL2d(ifirst,ilast,ngeta))
      double precision fl(CELL2d(ifirst,ilast,0))
      double precision fa(CELL2d(ifirst,ilast,0))
      double precision fb(CELL2d(ifirst,ilast,0))
      double precision rhs(CELL2d(ifirst,ilast,ngrhs))
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0, ic1

      double precision hphi, heta_prime
      double precision interp_func
      double precision deriv_interp_func
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0

            hphi = interp_func( phi(ic0,ic1), phi_interp_type )

            heta_prime =
     &         deriv_interp_func( eta(ic0,ic1), eta_interp_type )

            rhs(ic0,ic1) = rhs(ic0,ic1) +
     &         hphi * heta_prime *
     &         ( fa(ic0,ic1) - fb(ic0,ic1) )

         enddo
      enddo

      return
      end

c***********************************************************************
c
      subroutine computerhstemp(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   dx, 
     &   thermal_diffusivity,
     &   latent_heat,
     &   temp, ngtemp,
     &   cp, ngcp,
     &   with_phase,
     &   phi_rhs, ngphi_rhs,
     &   rhs, ngrhs )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1

      double precision dx(0:1)
      double precision thermal_diffusivity, latent_heat
      integer ngphi_rhs, ngtemp, ngrhs, ngcp, with_phase
c
c variables in 2d cell indexed
      double precision phi_rhs(CELL2d(ifirst,ilast,ngphi_rhs))
      double precision rhs(CELL2d(ifirst,ilast,ngrhs))
      double precision temp(CELL2d(ifirst,ilast,ngtemp))
      double precision cp(CELL2d(ifirst,ilast,ngcp))
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0, ic1
      double precision g, dxinv2, dyinv2
      double precision diag, offdiagx, offdiagy
     
c      print*,'thermal diff=',thermal_diffusivity
 
      dxinv2 = 1.d0/(dx(0)*dx(0))
      dyinv2 = 1.d0/(dx(1)*dx(1))
c
      offdiagx = thermal_diffusivity * dxinv2
      offdiagy = thermal_diffusivity * dyinv2
      diag = -2.d0 * (offdiagx + offdiagy)

      call stencil5pts(ifirst0, ilast0, ifirst1, ilast1,
     &                 diag, offdiagx, offdiagy, temp, ngtemp,
     &                 rhs, ngrhs )

      if( with_phase /= 0 )then
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               g = latent_heat/cp(ic0,ic1)
c               print*,'gamma=',g,', diff=',thermal_diffusivity
            
               rhs(ic0,ic1) = rhs(ic0,ic1) +
     &            g * phi_rhs(ic0,ic1)

            enddo
         enddo
      endif

      return
      end

c***********************************************************************
c
      subroutine computerhsbiaswell(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phi, ngphi,
     &   temp, ngtemp,
     &   alpha, gamma, 
     &   te, ngte,
     &   rhs, ngrhs )
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer ngphi, ngtemp, ngrhs, ngte
      double precision alpha, gamma
c
c variables in 2d cell indexed
      double precision phi(CELL2d(ifirst,ilast,ngphi))
      double precision rhs(CELL2d(ifirst,ilast,ngrhs))
      double precision temp(CELL2d(ifirst,ilast,ngtemp))
      double precision te(CELL2d(ifirst,ilast,ngte))
c
      double precision pi, coeff, m
      integer ic0, ic1
c
      pi = 4.*atan(1.)
      coeff = alpha/pi
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0

            m = coeff*atan(gamma*(te(ic0,ic1)-temp(ic0,ic1)))
c this is just the contribution in addition to regular double well
            rhs(ic0,ic1) = rhs(ic0,ic1) +
     &         m*phi(ic0,ic1)*(1.d0-phi(ic0,ic1))

         enddo
      enddo

      return
      end

c***********************************************************************
c
      subroutine computerhsbiaswellbeckermann(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phi, ngphi,
     &   temp, ngtemp,
     &   alpha, 
     &   te, ngte,
     &   rhs, ngrhs )
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer ngphi, ngtemp, ngrhs, ngte
      double precision alpha
c
c variables in 2d cell indexed
      double precision phi(CELL2d(ifirst,ilast,ngphi))
      double precision rhs(CELL2d(ifirst,ilast,ngrhs))
      double precision temp(CELL2d(ifirst,ilast,ngtemp))
      double precision te(CELL2d(ifirst,ilast,ngte))
c
      double precision m
      integer ic0, ic1
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0

            m = alpha*(te(ic0,ic1)-temp(ic0,ic1))
c this is just the contribution in addition to regular double well
            rhs(ic0,ic1) = rhs(ic0,ic1) +
     &         m*phi(ic0,ic1)*(1.d0-phi(ic0,ic1))

         enddo
      enddo

      return
      end

c***********************************************************************
c
      subroutine computerhsdeltatemperature(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phi, ngphi,
     &   temp, ngtemp,
     &   tm, latentheat,
     &   rhs, ngrhs,
     &   energy_interp_type )
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer ngphi, ngtemp, ngrhs
      double precision tm, latentheat
      character*(*) energy_interp_type
c
c variables in 2d cell indexed
      double precision phi(CELL2d(ifirst,ilast,ngphi))
      double precision rhs(CELL2d(ifirst,ilast,ngrhs))
      double precision temp(CELL2d(ifirst,ilast,ngtemp))
c
      integer ic0, ic1
      double precision m, alpha, wtemp, h_prime
      double precision deriv_interp_func
c
c      print*,'latentheat=',latentheat,', tm=',tm
      alpha = latentheat/tm

      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0

            wtemp = 0.75*temp(ic0,ic1)
     &        +0.0625*( temp(ic0-1,ic1)
     &                 +temp(ic0,ic1-1)
     &                 +temp(ic0+1,ic1)
     &                 +temp(ic0,ic1+1) )

            m = alpha*( tm-wtemp )
            h_prime =
     &            deriv_interp_func( phi(ic0,ic1), energy_interp_type )

            rhs(ic0,ic1) = rhs(ic0,ic1) + m*h_prime

         enddo
      enddo

      return
      end

c***********************************************************************
c
      subroutine computedphidtemperaturedeltatemperature(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   phi, ngphi,
     &   tm, latentheat, mobility,
     &   rhs, ngrhs,
     &   energy_interp_type )
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer ngphi, ngrhs
      double precision tm, latentheat, mobility
      character*(*) energy_interp_type
c
c variables in 2d cell indexed
      double precision phi(CELL2d(ifirst,ilast,ngphi))
      double precision rhs(CELL2d(ifirst,ilast,ngrhs))
c
      integer ic0, ic1
      double precision alpha, h_prime
      double precision deriv_interp_func
c
c      print*,'latentheat=',latentheat,', tm=',tm
      alpha = mobility*latentheat/tm

      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0

            h_prime =
     &            deriv_interp_func( phi(ic0,ic1), energy_interp_type )

            rhs(ic0,ic1) = alpha*h_prime

         enddo
      enddo

      return
      end

c***********************************************************************
c
      subroutine addvdphidx(
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

      double precision dx(2), vel
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
      double precision dxinv, diff_term_x

      dxinv = 0.5d0 * vel / dx(1)
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0

            diff_term_x = dxinv *
     &           (phase(ic0+1,ic1) - phase(ic0-1,ic1))

            rhs(ic0,ic1) = rhs(ic0,ic1)
     &                   + diff_term_x

         enddo
      enddo

      end
c***********************************************************************
c
      subroutine addvdphidx_upwind3(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   dx,
     &   phase, ngphase,
     &   vel,
     &   rhs, ngrhs,
     &   physb)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1

      double precision dx(2), vel
      integer ngphase, ngrhs
      integer physb
c
c variables in 2d cell indexed
      double precision phase(CELL2d(ifirst,ilast,ngphase))
      double precision rhs(CELL2d(ifirst,ilast,ngrhs))

c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1, imax0
      double precision dxinv, diff_term_x

      dxinv = vel / (6.d0*dx(1))
c
      do ic1 = ifirst1, ilast1
         if( physb .eq. 1 ) then
c use regular 2nd order at right boundary since 2nd ghost value
c typically not filled with physical boundar value
            imax0 = ilast0-1
            diff_term_x = 3.d0*dxinv *
     &         (phase(ilast0+1,ic1)-phase(ilast0-1,ic1))
            rhs(ilast0,ic1) = rhs(ilast0,ic1) + diff_term_x
         else
            imax0 = ilast0
         endif
         do ic0 = ifirst0, imax0
            diff_term_x = dxinv *
     &           (-2.d0*phase(ic0-1,ic1) -3.d0*phase(ic0,ic1)
     &            +6.d0*phase(ic0+1,ic1) - phase(ic0+2,ic1))

               rhs(ic0,ic1) = rhs(ic0,ic1)
     &                      + diff_term_x
         enddo
      enddo

      end
