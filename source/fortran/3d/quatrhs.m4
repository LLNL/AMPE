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
define(NDIM,3)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine gradient_flux(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   h, epsilon,
     &   phase, ngphase,
     &   flux0, flux1, flux2, ngflux)

      implicit none
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &     ngphase, ngflux
      double precision 
     &     phase(CELL3d(ifirst,ilast,ngphase)),
     &     flux0(SIDE3d0(ifirst,ilast,ngflux)),
     &     flux1(SIDE3d1(ifirst,ilast,ngflux)),
     &     flux2(SIDE3d2(ifirst,ilast,ngflux)),
     &     h(3)

c     local variables
      integer i, j, k
      double precision dxinv, dyinv, dzinv, epsilon

      double precision epsilon2
c
      epsilon2 = epsilon * epsilon

      dxinv = epsilon2 / h(1)
      dyinv = epsilon2 / h(2)
      dzinv = epsilon2 / h(3)

      do k = ifirst2, ilast2
         do j = ifirst1, ilast1
            do i = ifirst0, ilast0+1
               flux0(i,j,k) = (phase(i,j,k) - phase(i-1,j,k)) * dxinv
            enddo
         enddo
      enddo

      do k = ifirst2, ilast2
         do j = ifirst1, ilast1+1
            do i = ifirst0, ilast0
               flux1(i,j,k) = (phase(i,j,k) - phase(i,j-1,k)) * dyinv
            enddo
         enddo
      enddo

      do k = ifirst2, ilast2+1
         do j = ifirst1, ilast1
            do i = ifirst0, ilast0
               flux2(i,j,k) = (phase(i,j,k) - phase(i,j,k-1)) * dzinv
            enddo
         enddo
      enddo

      return
      end

c***********************************************************************
      subroutine compute_flux_isotropic(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   h, epsilon,
     &   phase, ngphase,
     &   flux0, flux1, flux2, ngflux)

      implicit none
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &     ngphase, ngflux
      double precision 
     &     phase(CELL3d(ifirst,ilast,ngphase)),
     &     flux0(SIDE3d0(ifirst,ilast,ngflux)),
     &     flux1(SIDE3d1(ifirst,ilast,ngflux)),
     &     flux2(SIDE3d2(ifirst,ilast,ngflux)),
     &     h(3)

c     local variables
      integer i, j, k
      double precision dxinv, dyinv, dzinv, epsilon

      double precision epsilon2
c
      print*,'compute_flux_isotropic():'
      print*,'implementation incomplete in 3D!!!'
      stop
      
      epsilon2 = epsilon * epsilon

      dxinv = epsilon2 / h(1)
      dyinv = epsilon2 / h(2)
      dzinv = epsilon2 / h(3)

      do k = ifirst2, ilast2
         do j = ifirst1, ilast1
            do i = ifirst0, ilast0+1
               flux0(i,j,k) = 0.0625d0*dxinv*( 
     &              (phase(i,j-1,k) - phase(i-1,j-1,k)) 
     &            + (phase(i,j  ,k) - phase(i-1,j  ,k))*6.d0 
     &            + (phase(i,j+1,k) - phase(i-1,j+1,k)) 
     &            + (phase(i,j,k-1) - phase(i-1,j,k-1))
     &            + (phase(i,j,k)   - phase(i-1,j,k  ))*6.d0 
     &            + (phase(i,j,k+1) - phase(i-1,j,k+1)) )
            enddo
         enddo
      enddo

      do k = ifirst2, ilast2
         do j = ifirst1, ilast1+1
            do i = ifirst0, ilast0
               flux1(i,j,k) = 0.0625d0*dyinv*(
     &              (phase(i-1,j,k) - phase(i-1,j-1,k)) 
     &            + (phase(i,  j,k) - phase(i,  j-1,k))*6.d0 
     &            + (phase(i+1,j,k) - phase(i+1,j-1,k))
     &            + (phase(i,j,k-1) - phase(i,j-1,k-1)) 
     &            + (phase(i,j,k  ) - phase(i,j-1,k  ))*6.d0 
     &            + (phase(i,j,k+1) - phase(i,j-1,k+1)) )
            enddo
         enddo
      enddo

      do k = ifirst2, ilast2+1
         do j = ifirst1, ilast1
            do i = ifirst0, ilast0
               flux2(i,j,k) = 0.0625d0*dzinv*(
     &              (phase(i-1,j,k) - phase(i-1,j,k-1)) 
     &            + (phase(i,  j,k) - phase(i,  j,k-1))*6.d0 
     &            + (phase(i+1,j,k) - phase(i+1,j,k-1))
     &            + (phase(i,j-1,k) - phase(i,j-1,k-1)) 
     &            + (phase(i,j,  k) - phase(i,j,  k-1))*6.d0 
     &            + (phase(i,j+1,k) - phase(i,j+1,k-1)) )
            enddo
         enddo
      enddo

      return
      end

c***********************************************************************
      subroutine anisotropic_gradient_flux(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   h, epsilon, e4, knumber,
     &   phase, ngphase,
     &   quat, ngq, qlen,
     &   flux0, flux1, flux2, ngflux)

      implicit none
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &     ngphase, ngq, qlen, ngflux, knumber
      double precision 
     &     phase(CELL3d(ifirst,ilast,ngphase)),
     &     quat(CELL3d(ifirst,ilast,ngq),qlen),
     &     flux0(SIDE3d0(ifirst,ilast,ngflux)),
     &     flux1(SIDE3d1(ifirst,ilast,ngflux)),
     &     flux2(SIDE3d2(ifirst,ilast,ngflux)),
     &     e4, h(3)

c     local variables
      integer i, j, k
      double precision dxinv, dyinv, dzinv, epsilon
      double precision pi, q(4), qp(4), n(4), qtmp(4), np(4)
      double precision epsilon2, dphidx, dphidy, dphidz
      double precision gphi2, nni, gamma
c
      pi = 4.*atan(1.)
c
      epsilon2 = epsilon * epsilon

      dxinv = 1. / h(1)
      dyinv = 1. / h(2)
      dzinv = 1. / h(3)

c do averaging across face when taking derivative parallel to face

c x faces      
      do k = ifirst2, ilast2
         do j = ifirst1, ilast1
            do i = ifirst0, ilast0+1
               dphidx = (phase(i,j,k) - phase(i-1,j,k)) * dxinv
               dphidy = 0.25*(phase(i-1,j+1,k) - phase(i-1,j-1,k)
     &                      + phase(i  ,j+1,k) - phase(i  ,j-1,k)) 
     &                       * dyinv 
               dphidz = 0.25*(phase(i-1,j,k+1) - phase(i-1,j,k-1)
     &                      + phase(i  ,j,k+1) - phase(i  ,j,k-1)) 
     &                       * dzinv 
         
               gphi2=dphidx*dphidx+dphidy*dphidy+dphidz*dphidz
               
               n(1)=0.
               if( abs(gphi2)>1.e-12 )then
                  nni=1./sqrt(gphi2)
                  n(2)=dphidx*nni;
                  n(3)=dphidy*nni;
                  n(4)=dphidz*nni;
               else
                  n(2)=0.;
                  n(3)=0.;
                  n(4)=1.;
               endif
               
               q(1)=0.5*(quat(i-1,j,k,1)+quat(i,j,k,1))
               q(2)=0.5*(quat(i-1,j,k,2)+quat(i,j,k,2))
               q(3)=0.5*(quat(i-1,j,k,3)+quat(i,j,k,3))
               q(4)=0.5*(quat(i-1,j,k,4)+quat(i,j,k,4))

               call quatconj(q,qp)
               call quatmult4(n,qp,qtmp)
               call quatmult4(q,qtmp,np)

               gamma=(1.-3.*e4)*(1.+4.*e4*(np(2)**4+np(3)**4+np(4)**4)
     &                          /(1.-3.*e4))
               
               flux0(i,j,k) = epsilon2*dphidx 
     &            + gphi2*epsilon2*gamma*e4*16.*np(2)*np(2)*np(2)
            enddo
         enddo
      enddo

c y faces      
      do k = ifirst2, ilast2
         do j = ifirst1, ilast1+1
            do i = ifirst0, ilast0
               dphidx = 0.25*(phase(i+1,j-1,k) - phase(i-1,j-1,k)
     &                      + phase(i+1,j  ,k) - phase(i-1,j  ,k)) 
     &                      * dxinv 
               dphidy = (phase(i,j,k) - phase(i,j-1,k)) * dyinv
               dphidz = 0.25*(phase(i,j-1,k+1) - phase(i,j-1,k-1)
     &                      + phase(i,j  ,k+1) - phase(i,j  ,k-1)) 
     &                      * dzinv 
            
               gphi2=dphidx*dphidx+dphidy*dphidy+dphidz*dphidz
               
               n(1)=0.
               if( abs(gphi2)>1.e-12 )then
                  nni=1./sqrt(gphi2)
                  n(2)=dphidx*nni;
                  n(3)=dphidy*nni;
                  n(4)=dphidz*nni;
               else
                  n(2)=0.;
                  n(3)=0.;
                  n(4)=1.;
               endif
               
               q(1)=0.5*(quat(i-1,j,k,1)+quat(i,j,k,1))
               q(2)=0.5*(quat(i-1,j,k,2)+quat(i,j,k,2))
               q(3)=0.5*(quat(i-1,j,k,3)+quat(i,j,k,3))
               q(4)=0.5*(quat(i-1,j,k,4)+quat(i,j,k,4))

               call quatconj(q,qp)
               call quatmult4(n,qp,qtmp)
               call quatmult4(q,qtmp,np)

               gamma=(1.-3.*e4)*(1.+4.*e4*(np(2)**4+np(3)**4+np(4)**4)
     &                          /(1.-3.*e4))
               
               flux1(i,j,k) = epsilon2*dphidy 
     &              + gphi2*epsilon2*gamma*e4*16.*np(3)*np(3)*np(3)
            enddo
         enddo
      enddo

c z faces      
      do k = ifirst2, ilast2+1
         do j = ifirst1, ilast1
            do i = ifirst0, ilast0
               dphidx = 0.25*(phase(i+1,j,k-1) - phase(i-1,j,k-1)
     &                      + phase(i+1,j,k  ) - phase(i-1,j,k  )) 
     &                      * dxinv 
               dphidy = 0.25*(phase(i-1,j+1,k-1) - phase(i,j-1,k-1)
     &                      + phase(i  ,j+1,k  ) - phase(i,j-1,k  )) 
     &                      * dyinv             
               dphidz = (phase(i,j,k) - phase(i,j,k-1)) * dzinv

               gphi2=dphidx*dphidx+dphidy*dphidy+dphidz*dphidz
               
               n(1)=0.
               if( abs(gphi2)>1.e-12 )then
                  nni=1./sqrt(gphi2)
                  n(2)=dphidx*nni;
                  n(3)=dphidy*nni;
                  n(4)=dphidz*nni;
               else
                  n(2)=0.;
                  n(3)=0.;
                  n(4)=1.;
               endif
               
               q(1)=0.5*(quat(i-1,j,k,1)+quat(i,j,k,1))
               q(2)=0.5*(quat(i-1,j,k,2)+quat(i,j,k,2))
               q(3)=0.5*(quat(i-1,j,k,3)+quat(i,j,k,3))
               q(4)=0.5*(quat(i-1,j,k,4)+quat(i,j,k,4))

               call quatconj(q,qp)
               call quatmult4(n,qp,qtmp)
               call quatmult4(q,qtmp,np)

               gamma=(1.-3.*e4)*(1.+4.*e4*(np(2)**4+np(3)**4+np(4)**4)
     &                          /(1.-3.*e4))
               
               flux2(i,j,k) = epsilon2*dphidz
     &              + gphi2*epsilon2*gamma*e4*16.*np(4)*np(4)*np(4)

            enddo
         enddo
      enddo

      return
      end
      
c***********************************************************************
c
c compute r.h.s. for phase variable phi (Pusztai et al. model)
*
c   7 point stencil (3d extension of 5pt)
c
      subroutine computerhspbg(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   dx, 
     &   misorientation_factor,
     &   flux0,
     &   flux1,
     &   flux2,
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
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer with_orient, three_phase

      double precision dx(3)
      double precision misorientation_factor
      double precision phi_well_scale, eta_well_scale
      integer ngflux, ngphi, ngeta, ngogm, ngrhs, ngtemp
      character*(*) phi_well_type
      character*(*) eta_well_type
      character*(*) phi_interp_type
      character*(*) orient_interp_type
c
c variables in 3d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi))
      double precision eta(CELL3d(ifirst,ilast,ngeta))
      double precision orient_grad_mod(CELL3d(ifirst,ilast,ngogm))
      double precision rhs(CELL3d(ifirst,ilast,ngrhs))
      double precision temp(CELL3d(ifirst,ilast,ngtemp))

c variables in 3d side indexed
      double precision 
     &     flux0(SIDE3d0(ifirst,ilast,ngflux)),
     &     flux1(SIDE3d1(ifirst,ilast,ngflux)),
     &     flux2(SIDE3d2(ifirst,ilast,ngflux))
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0, ic1, ic2
      double precision diff_term_x, diff_term_y, diff_term_z, diff_term

      double precision g, g_prime, h_prime, p_prime
      double precision deriv_interp_func
      double precision well_func
      double precision deriv_well_func
      double precision t
      double precision dxinv, dyinv, dzinv
c
      dxinv = 1.d0 / dx(1)
      dyinv = 1.d0 / dx(2)
      dzinv = 1.d0 / dx(3)
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               diff_term_x = 
     &              (flux0(ic0+1,ic1,ic2) - flux0(ic0,ic1,ic2)) * dxinv
               diff_term_y = 
     &              (flux1(ic0,ic1+1,ic2) - flux1(ic0,ic1,ic2)) * dyinv
               diff_term_z = 
     &              (flux2(ic0,ic1,ic2+1) - flux2(ic0,ic1,ic2)) * dzinv

               diff_term = diff_term_x + diff_term_y + diff_term_z

               rhs(ic0,ic1,ic2) = diff_term

c  Phase energy well

               g_prime =
     &            deriv_well_func( phi(ic0,ic1,ic2), phi_well_type )

               rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2) -
     &            phi_well_scale * g_prime

            enddo
         enddo
      enddo

c  Eta energy well

      if ( three_phase /= 0 ) then

         do ic2 = ifirst2, ilast2
            do ic1 = ifirst1, ilast1
               do ic0 = ifirst0, ilast0

                  h_prime =
     &               deriv_interp_func(
     &                  phi(ic0,ic1,ic2),
     &                  phi_interp_type )
               
                  g = well_func( eta(ic0,ic1,ic2), eta_well_type )

                  rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2) -
     &               eta_well_scale * h_prime * g

               enddo
            enddo
         enddo

      endif

c  Misorientation energy
c
      if ( with_orient /= 0 ) then

         do ic2 = ifirst2, ilast2
            do ic1 = ifirst1, ilast1
               do ic0 = ifirst0, ilast0

                  p_prime =
     &               deriv_interp_func(
     &                  phi(ic0,ic1,ic2),
     &                  orient_interp_type )

                  rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2) -
     &               misorientation_factor *
     &               temp(ic0,ic1,ic2) *
     &               p_prime *
     &               orient_grad_mod(ic0,ic1,ic2)

               enddo
            enddo
         enddo

      endif

      return
      end

c***********************************************************************
c
c compute r.h.s. for eta variable phi
c
c   7 point stencil
c
      subroutine computerhseta(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   dx,
     &   flux0,
     &   flux1,
     &   flux2,
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
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2

      double precision dx(0:2)
      double precision eta_well_scale
      integer ngphi, ngeta, ngrhs, ngtemp, ngflux
      character*(*) eta_well_type
      character*(*) phi_interp_type
c
      double precision phi(CELL3d(ifirst,ilast,ngphi))
      double precision eta(CELL3d(ifirst,ilast,ngeta))
      double precision rhs(CELL3d(ifirst,ilast,ngrhs))
      double precision temp(CELL3d(ifirst,ilast,ngtemp))
c variables in 3d side indexed
      double precision 
     &     flux0(SIDE3d0(ifirst,ilast,ngflux)),
     &     flux1(SIDE3d1(ifirst,ilast,ngflux)),
     &     flux2(SIDE3d1(ifirst,ilast,ngflux))
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0, ic1, ic2
      double precision diff_term_x, diff_term_y, diff_term_z, diff_term
      double precision dxinv, dyinv, dzinv

      double precision g_prime, h
      double precision interp_func
      double precision deriv_well_func
      double precision t
c
      dxinv = 1.d0 / dx(0)
      dyinv = 1.d0 / dx(1)
      dzinv = 1.d0 / dx(2)
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               diff_term_x = 
     &           ( flux0(ic0+1,ic1,ic2) - flux0(ic0,ic1,ic2) )*dxinv
               diff_term_y = 
     &           ( flux1(ic0,ic1+1,ic2) - flux1(ic0,ic1,ic2) )*dyinv
               diff_term_z = 
     &           ( flux1(ic0,ic1,ic2+1) - flux1(ic0,ic1,ic2) )*dzinv

               diff_term = diff_term_x + diff_term_y + diff_term_z

               rhs(ic0,ic1,ic2) = diff_term

c  Energy well

               h = interp_func( phi(ic0,ic1,ic2), phi_interp_type )

               g_prime =
     &            deriv_well_func( eta(ic0,ic1,ic2), eta_well_type )

               rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2) -
     &            eta_well_scale * h * g_prime

            enddo
         enddo
      enddo

      return
      end

c***********************************************************************
c
c compute r.h.s. component due to free energy for phase variable phi
c
      subroutine phaserhs_fenergy(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
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
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer three_phase

      integer ngphi, ngeta, ngrhs
      character*(*) phi_interp_type
      character*(*) eta_interp_type
c
c variables in 3d cell indexed
      double precision fl(CELL3d(ifirst,ilast,0))
      double precision fa(CELL3d(ifirst,ilast,0))
      double precision fb(CELL3d(ifirst,ilast,0))
      double precision phi(CELL3d(ifirst,ilast,ngphi))
      double precision eta(CELL3d(ifirst,ilast,ngeta))
      double precision rhs(CELL3d(ifirst,ilast,ngrhs))
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0, ic1, ic2
      double precision f_l, f_a, f_b

      double precision heta, hphi_prime
      double precision interp_func
      double precision deriv_interp_func
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               hphi_prime = 
     &            deriv_interp_func(
     &               phi(ic0,ic1,ic2),
     &               phi_interp_type )

               rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2) +
     &            hphi_prime *
     &            ( fl(ic0,ic1,ic2) - fa(ic0,ic1,ic2) )

            enddo
         enddo
      enddo

      if ( three_phase /= 0 ) then

         do ic2 = ifirst2, ilast2
            do ic1 = ifirst1, ilast1
               do ic0 = ifirst0, ilast0

                  hphi_prime =
     &               deriv_interp_func(
     &                  phi(ic0,ic1,ic2),
     &                  phi_interp_type )

                  heta =
     &               interp_func( eta(ic0,ic1,ic2), eta_interp_type )

                  rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2) +
     &               hphi_prime * heta *
     &               ( fa(ic0,ic1,ic2) - fb(ic0,ic1,ic2) )

               enddo
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
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
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
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2

      integer ngphi, ngeta, ngrhs
      character*(*) phi_interp_type
      character*(*) eta_interp_type
c
c variables in 3d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi))
      double precision eta(CELL3d(ifirst,ilast,ngeta))
      double precision fl(CELL3d(ifirst,ilast,0))
      double precision fa(CELL3d(ifirst,ilast,0))
      double precision fb(CELL3d(ifirst,ilast,0))
      double precision rhs(CELL3d(ifirst,ilast,ngrhs))
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0, ic1, ic2
      double precision f_l, f_a, f_b

      double precision hphi, heta_prime
      double precision interp_func
      double precision deriv_interp_func
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               hphi = interp_func( phi(ic0,ic1,ic2), phi_interp_type )

               heta_prime =
     &            deriv_interp_func(
     &               eta(ic0,ic1,ic2),
     &               eta_interp_type )

               rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2) +
     &            hphi * heta_prime *
     &            ( fa(ic0,ic1,ic2) - fb(ic0,ic1,ic2) )

            enddo
         enddo
      enddo

      return
      end

c***********************************************************************
c
      subroutine computerhstemp(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
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
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer nsp

      double precision dx(0:2)
      double precision thermal_diffusivity, latent_heat
      integer ngphi_rhs, ngtemp, ngrhs, ngcp
c
c variables in 3d cell indexed
      double precision phi_rhs(CELL3d(ifirst,ilast,ngphi_rhs))
      double precision rhs(CELL3d(ifirst,ilast,ngrhs))
      double precision temp(CELL3d(ifirst,ilast,ngtemp))
      double precision cp(CELL3d(ifirst,ilast,ngcp))
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0, ic1, ic2
      double precision gamma, dxinv2, dyinv2, dzinv2
      double precision diff_term_x, diff_term_y, diff_term_z
      double precision diff_term
      
      dxinv2 = 1.d0/(dx(0)*dx(0))
      dyinv2 = 1.d0/(dx(1)*dx(1))
      dzinv2 = 1.d0/(dx(2)*dx(2))
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               gamma = latent_heat/cp(ic0,ic1,ic2)
            
               diff_term_x = 
     &           (temp(ic0-1,ic1,ic2)-2.d0*temp(ic0,ic1,ic2)
     &           +temp(ic0+1,ic1,ic2)) 
               diff_term_y = 
     &           (temp(ic0,ic1-1,ic2)-2.d0*temp(ic0,ic1,ic2)
     &           +temp(ic0,ic1+1,ic2)) 
               diff_term_z = 
     &           (temp(ic0,ic1,ic2-1)-2.d0*temp(ic0,ic1,ic2)
     &           +temp(ic0,ic1,ic2+1)) 

               diff_term = diff_term_x*dxinv2 + diff_term_y*dyinv2
     &                   + diff_term_z*dzinv2

               rhs(ic0,ic1,ic2) = thermal_diffusivity * diff_term

               rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2) +
     &            gamma * phi_rhs(ic0,ic1,ic2)

            enddo
         enddo
      enddo

      return
      end

c***********************************************************************
c
      subroutine computerhsbiaswell(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phi, ngphi,
     &   temp, ngtemp,
     &   alpha, gamma,
     &   te, ngte,
     &   rhs, ngrhs )
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2

      integer ngphi, ngtemp, ngrhs, ngte
      double precision alpha, gamma
c
c variables in 2d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi))
      double precision rhs(CELL3d(ifirst,ilast,ngrhs))
      double precision temp(CELL3d(ifirst,ilast,ngtemp))
      double precision te(CELL3d(ifirst,ilast,ngte))
c
      double precision pi, coeff, m
      integer ic0, ic1, ic2
c
      pi = 4.*atan(1.)
      coeff = alpha/pi
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               m = coeff*atan(gamma*(te(ic0,ic1,ic2)-temp(ic0,ic1,ic2)))
c this is just the contribution in addition to regular double well
               rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2) +
     &         m*phi(ic0,ic1,ic2)*(1.d0-phi(ic0,ic1,ic2))

            enddo
         enddo
      enddo

      return
      end

c***********************************************************************
c
      subroutine computerhsbiaswellbeckermann(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phi, ngphi,
     &   temp, ngtemp,
     &   alpha,
     &   te, ngte,
     &   rhs, ngrhs )
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2

      integer ngphi, ngtemp, ngrhs, ngte
      double precision alpha
c
c variables in 2d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi))
      double precision rhs(CELL3d(ifirst,ilast,ngrhs))
      double precision temp(CELL3d(ifirst,ilast,ngtemp))
      double precision te(CELL3d(ifirst,ilast,ngte))
c
      double precision m
      integer ic0, ic1, ic2
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               m = alpha*(te(ic0,ic1,ic2)-temp(ic0,ic1,ic2))
c this is just the contribution in addition to regular double well
               rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2) +
     &         m*phi(ic0,ic1,ic2)*(1.d0-phi(ic0,ic1,ic2))

            enddo
         enddo
      enddo

      return
      end

c***********************************************************************
c
      subroutine computerhsdeltatemperature(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phi, ngphi,
     &   temp, ngtemp,
     &   tm, latentheat,
     &   rhs, ngrhs )
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngphi, ngtemp, ngrhs
      double precision tm, latentheat
c
c variables in 3d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi))
      double precision rhs(CELL3d(ifirst,ilast,ngrhs))
      double precision temp(CELL3d(ifirst,ilast,ngtemp))
c
      integer ic0, ic1, ic2
      double precision m, alpha
c
      alpha = latentheat/tm
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               m = alpha*(tm-temp(ic0,ic1,ic2))

               rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2) +
     &          m*phi(ic0,ic1,ic2)*(1.d0-phi(ic0,ic1,ic2))

            enddo
         enddo
      enddo

      return
      end

