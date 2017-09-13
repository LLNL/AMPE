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
      double precision gx, gy, gn, gni, theta, phi
      double precision nx, ny
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
            
            flux0(i,j) = epsilon2*dphidx 
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
            
            flux1(i,j) = epsilon2*dphidy 
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
     &   phi_interp_type,
     &   orient_interp_type, 
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
      double precision misorientation_factor
      double precision phi_well_scale, eta_well_scale
      integer ngflux, ngphi, ngeta, ngogm, ngrhs, ngtemp
      character*(*) phi_well_type
      character*(*) eta_well_type
      character*(*) phi_interp_type
      character*(*) orient_interp_type
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

      double precision g, g_prime, h_prime, p_prime
      double precision deriv_interp_func
      double precision well_func
      double precision deriv_well_func
      double precision t
      double precision dxinv, dyinv
c
      dxinv = 1.d0 / dx(1)
      dyinv = 1.d0 / dx(2)
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
     &            deriv_interp_func( phi(ic0,ic1), phi_interp_type )
               
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

               p_prime =
     &            deriv_interp_func(
     &               phi(ic0,ic1),
     &               orient_interp_type )

               rhs(ic0,ic1) = rhs(ic0,ic1) -
     &            misorientation_factor *
     &            temp(ic0,ic1) *
     &            p_prime * 
     &            orient_grad_mod(ic0,ic1)
c     &            orient_grad_mod(ic0,ic1) *
c     &            ( 1.d0 +
c     &              epsilon_q2 * orient_grad_mod(ic0,ic1) )

            enddo
         enddo

      endif

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
      double precision t
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
     &   fl, fa, fb,
     &   phi, ngphi,
     &   eta, ngeta,
     &   rhs, ngrhs,
     &   phi_interp_type,
     &   eta_interp_type,
     &   three_phase )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer three_phase

      integer ngphi, ngeta, ngrhs
      character*(*) phi_interp_type
      character*(*) eta_interp_type
c
c variables in 2d cell indexed
      double precision fl(CELL2d(ifirst,ilast,0))
      double precision fa(CELL2d(ifirst,ilast,0))
      double precision fb(CELL2d(ifirst,ilast,0))
      double precision phi(CELL2d(ifirst,ilast,ngphi))
      double precision eta(CELL2d(ifirst,ilast,ngeta))
      double precision rhs(CELL2d(ifirst,ilast,ngrhs))
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0, ic1
      double precision f_l, f_a, f_b

      double precision heta, hphi_prime
      double precision interp_func
      double precision deriv_interp_func
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0

            hphi_prime = 
     &         deriv_interp_func( phi(ic0,ic1), phi_interp_type )

            rhs(ic0,ic1) = rhs(ic0,ic1) +
     &         hphi_prime *
     &         ( fl(ic0,ic1) - fa(ic0,ic1) )

         enddo
      enddo

      if ( three_phase /= 0 ) then

         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               hphi_prime =
     &            deriv_interp_func( phi(ic0,ic1), phi_interp_type )

               heta = interp_func( eta(ic0,ic1), eta_interp_type )

               rhs(ic0,ic1) = rhs(ic0,ic1) +
     &            hphi_prime * heta *
     &            ( fa(ic0,ic1) - fb(ic0,ic1) )

            enddo
         enddo

      endif

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
      double precision f_l, f_a, f_b

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
     &   phi_rhs, ngphi_rhs,
     &   rhs, ngrhs )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer nsp

      double precision dx(0:1)
      double precision thermal_diffusivity, latent_heat
      integer ngphi_rhs, ngtemp, ngrhs, ngcp
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
      double precision gamma, dxinv2, dyinv2
      double precision diff_term_x, diff_term_y, diff_term
      
      dxinv2 = 1.d0/(dx(0)*dx(0))
      dyinv2 = 1.d0/(dx(1)*dx(1))
c
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0

            gamma = latent_heat/cp(ic0,ic1)
            
            diff_term_x = 
     &           (temp(ic0-1,ic1)-2.d0*temp(ic0,ic1)+temp(ic0+1,ic1)) 
            diff_term_y = 
     &           (temp(ic0,ic1-1)-2.d0*temp(ic0,ic1)+temp(ic0,ic1+1)) 

            diff_term = diff_term_x*dxinv2 + diff_term_y*dyinv2

            rhs(ic0,ic1) = thermal_diffusivity * diff_term

            rhs(ic0,ic1) = rhs(ic0,ic1) +
     &         gamma * phi_rhs(ic0,ic1)

         enddo
      enddo

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
c***********************************************************************
c***********************************************************************     
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

