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
c input: quat,n
c output: dgamma,n4
c
      subroutine compute_dgamma(quat,n,dgamma,n4)

      implicit none
      double precision quat(4),n(4),dgamma(4),n4

c     local variables
      double precision qp(4),qtmp(4),np(4)

      call quatconj(quat,qp)

c rotation applied to np
      call quatmult4(n,qp,qtmp)
      call quatmult4(quat,qtmp,np)
c      print 100,n(2),n(3),n(4),np(1),np(2),np(3),np(4)
c100   format (7F6.3)

      n4=np(2)**4+np(3)**4+np(4)**4

      dgamma(1)=0.
      dgamma(2)=np(2)*(np(2)*np(2)-n4)
      dgamma(3)=np(3)*(np(3)*np(3)-n4)
      dgamma(4)=np(4)*(np(4)*np(4)-n4)

c inverse rotation applied to dgamma
      call quatmult4(dgamma,quat,qtmp)
      call quatmult4(qp,qtmp,dgamma)

      return
      end

c***********************************************************************
      subroutine anisotropic_gradient_flux(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   h, epsilon, eps4, knumber,
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
     &     eps4, h(3)

c     local variables
      integer i, j, k
      double precision dxinv, dyinv, dzinv, epsilon
      double precision q(4), n(4), dgamma(4)
      double precision epsilon2, dphidx, dphidy, dphidz
      double precision gphi2, nni, gamma, factor, n4
      double precision threshold
c
c      print*,'anisotropic_gradient_flux, eps4=',eps4
      epsilon2 = epsilon * epsilon

      dxinv = 1. / h(1)
      dyinv = 1. / h(2)
      dzinv = 1. / h(3)

      factor = 4.*eps4/(1.-3.*eps4)
      threshold = 1.e-12

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
               
               if( abs(gphi2)>threshold )then
                  nni=1./sqrt(gphi2)
                  n(1)=0.
                  n(2)=dphidx*nni
                  n(3)=dphidy*nni
                  n(4)=dphidz*nni

                  q(1)=0.5*(quat(i-1,j,k,1)+quat(i,j,k,1))
                  q(2)=0.5*(quat(i-1,j,k,2)+quat(i,j,k,2))
                  q(3)=0.5*(quat(i-1,j,k,3)+quat(i,j,k,3))
                  q(4)=0.5*(quat(i-1,j,k,4)+quat(i,j,k,4))

                  call compute_dgamma(q,n,dgamma,n4)
               else
                  dgamma(1)=0.
                  dgamma(2)=0.
                  dgamma(3)=0.
                  dgamma(4)=1.
                  n4=0.
               endif
               
               gamma=epsilon*(1.-3.*eps4)*(1.+factor*n4)
               
               flux0(i,j,k) = gamma*gamma*dphidx 
     &            + 16.*epsilon*gamma*eps4*sqrt(gphi2)*dgamma(2)
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
               
               if( abs(gphi2)>threshold )then
                  nni=1./sqrt(gphi2)
                  n(1)=0.
                  n(2)=dphidx*nni
                  n(3)=dphidy*nni
                  n(4)=dphidz*nni

                  q(1)=0.5*(quat(i,j-1,k,1)+quat(i,j,k,1))
                  q(2)=0.5*(quat(i,j-1,k,2)+quat(i,j,k,2))
                  q(3)=0.5*(quat(i,j-1,k,3)+quat(i,j,k,3))
                  q(4)=0.5*(quat(i,j-1,k,4)+quat(i,j,k,4))

                  call compute_dgamma(q,n,dgamma,n4)
               else
                  dgamma(1)=0.
                  dgamma(2)=0.
                  dgamma(3)=0.
                  dgamma(4)=1.
                  n4=0.
               endif

               gamma=epsilon*(1.-3.*eps4)*(1.+factor*n4)
               
               flux1(i,j,k) = gamma*gamma*dphidy 
     &              + 16.*epsilon*gamma*eps4*sqrt(gphi2)*dgamma(3)
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
               dphidy = 0.25*(phase(i,j+1,k-1) - phase(i,j-1,k-1)
     &                      + phase(i,j+1,k  ) - phase(i,j-1,k  )) 
     &                      * dyinv             
               dphidz = (phase(i,j,k) - phase(i,j,k-1)) * dzinv

               gphi2=dphidx*dphidx+dphidy*dphidy+dphidz*dphidz
               
               if( abs(gphi2)>threshold )then
                  nni=1./sqrt(gphi2)
                  n(1)=0.
                  n(2)=dphidx*nni
                  n(3)=dphidy*nni
                  n(4)=dphidz*nni

                  q(1)=0.5*(quat(i,j,k-1,1)+quat(i,j,k,1))
                  q(2)=0.5*(quat(i,j,k-1,2)+quat(i,j,k,2))
                  q(3)=0.5*(quat(i,j,k-1,3)+quat(i,j,k,3))
                  q(4)=0.5*(quat(i,j,k-1,4)+quat(i,j,k,4))

                  call compute_dgamma(q,n,dgamma,n4)
               else
                  dgamma(1)=0.
                  dgamma(2)=0.
                  dgamma(3)=1.
                  dgamma(4)=0.
                  n4=0.
               endif
               
               gamma=epsilon*(1.-3.*eps4)*(1.+factor*n4)
               
               flux2(i,j,k) = gamma*gamma*dphidz
     &              + 16.*epsilon*gamma*eps4*sqrt(gphi2)*dgamma(4)

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
     &   epsilonq,
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
     &   orient_interp_type1,
     &   orient_interp_type2,
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
      double precision misorientation_factor, epsilonq
      double precision phi_well_scale, eta_well_scale
      integer ngflux, ngphi, ngeta, ngogm, ngrhs, ngtemp
      character*(*) phi_well_type
      character*(*) eta_well_type
      character*(*) phi_interp_type
      character*(*) orient_interp_type1, orient_interp_type2
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

      double precision g, g_prime, h_prime, p1_prime, p2_prime
      double precision deriv_interp_func
      double precision well_func
      double precision deriv_well_func
      double precision dxinv, dyinv, dzinv
      double precision epsilonq2
c
      dxinv = 1.d0 / dx(1)
      dyinv = 1.d0 / dx(2)
      dzinv = 1.d0 / dx(3)
      epsilonq2 = 0.5d0*epsilonq*epsilonq
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

                  p1_prime =
     &               deriv_interp_func(
     &                  phi(ic0,ic1,ic2),
     &                  orient_interp_type1 )
                  p2_prime =
     &               deriv_interp_func(
     &                  phi(ic0,ic1,ic2),
     &                  orient_interp_type2 )

                  rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2) -
     &               misorientation_factor *
     &               temp(ic0,ic1,ic2) *
     &               p1_prime *
     &               orient_grad_mod(ic0,ic1,ic2)
     &             - p2_prime * epsilonq2 *
     &               orient_grad_mod(ic0,ic1,ic2) *
     &               orient_grad_mod(ic0,ic1,ic2)

               enddo
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
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   dx,
     &   flux0,
     &   flux1,
     &   flux2,
     &   ngflux,
     &   phi_well_scale,
     &   phi, ngphi,
     &   rhs, ngrhs)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2

      double precision dx(3)
      double precision phi_well_scale
      integer ngflux, ngphi, ngrhs
c
c variables in 2d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi),3)
      double precision rhs(CELL3d(ifirst,ilast,ngrhs),3)

c variables in 2d side indexed
      double precision
     &     flux0(SIDE3d0(ifirst,ilast,ngflux),3),
     &     flux1(SIDE3d1(ifirst,ilast,ngflux),3),
     &     flux2(SIDE3d2(ifirst,ilast,ngflux),3)
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1, ic2
      integer ip, jp, ip1, ip2
      double precision diff_term_x, diff_term_y, diff_term_z
      double precision diff_term

      double precision g_prime
      double precision deriv_triple_well_func
      double precision dxinv, dyinv, dzinv, fac, dsum
c
      dxinv = 1.d0 / dx(1)
      dyinv = 1.d0 / dx(2)
      dzinv = 1.d0 / dx(3)
c
      do ip = 1, 3
         ip1 = MOD(ip,3)+1
         ip2 = MOD(ip+1,3)+1
         do ic2 = ifirst2, ilast2
            do ic1 = ifirst1, ilast1
               do ic0 = ifirst0, ilast0

               diff_term_x =
     &            (flux0(ic0+1,ic1,ic2,ip) - flux0(ic0,ic1,ic2,ip))
     &            * dxinv
               diff_term_y =
     &            (flux1(ic0,ic1+1,ic2,ip) - flux1(ic0,ic1,ic2,ip))
     &            * dyinv
               diff_term_z =
     &            (flux2(ic0,ic1,ic2+1,ip) - flux2(ic0,ic1,ic2,ip))
     &            * dzinv

               diff_term = diff_term_x + diff_term_y + diff_term_z

               rhs(ic0,ic1,ic2,ip) = diff_term

c  Phase energy well

               g_prime =
     &            deriv_triple_well_func(
     &               phi(ic0,ic1,ic2,ip), phi(ic0,ic1,ic2,ip1),
     &               phi(ic0,ic1,ic2,ip2) )

               rhs(ic0,ic1,ic2,ip) = rhs(ic0,ic1,ic2,ip) -
     &            phi_well_scale * g_prime

               enddo
            enddo
         enddo
      enddo

c enforce constraint sum components = 0
      fac = 1.d0/3.d0
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0
               dsum = 0.d0
               do jp = 1, 3
                 dsum = dsum +  rhs(ic0,ic1,ic2,jp)
               enddo
               dsum = dsum * fac
               do jp = 1, 3
                  rhs(ic0,ic1,ic2,jp) = rhs(ic0,ic1,ic2,jp) - dsum
               enddo
            enddo
         enddo
      enddo

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
     &           ( flux2(ic0,ic1,ic2+1) - flux2(ic0,ic1,ic2) )*dzinv

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
     &   fl, fa,
     &   phi, ngphi,
     &   rhs, ngrhs,
     &   phi_interp_type)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2

      integer ngphi, ngrhs
      character*(*) phi_interp_type
c
c variables in 3d cell indexed
      double precision fl(CELL3d(ifirst,ilast,0))
      double precision fa(CELL3d(ifirst,ilast,0))
      double precision phi(CELL3d(ifirst,ilast,ngphi))
      double precision rhs(CELL3d(ifirst,ilast,ngrhs))
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0, ic1, ic2

      double precision hphi_prime
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

      return
      end
c
c compute r.h.s. component due to free energy for phase variable phi
c
      subroutine phaserhs_fenergy_multiorderp(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   fl, fa,
     &   phi, ngphi, nphi,
     &   rhs, ngrhs,
     &   phi_interp_type)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2

      integer ngphi, ngrhs, nphi
      character*(*) phi_interp_type
c
c variables in 3d cell indexed
      double precision fl(CELL3d(ifirst,ilast,0))
      double precision fa(CELL3d(ifirst,ilast,0))
      double precision phi(CELL3d(ifirst,ilast,ngphi),nphi)
      double precision rhs(CELL3d(ifirst,ilast,ngrhs),nphi)
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1, ic2, ip

      double precision hphis, hphil, suminv2
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0
               hphis = 0.d0
               do ip = 1, nphi-1
                  hphis = hphis + phi(ic0,ic1,ic2,ip)
     &                           *phi(ic0,ic1,ic2,ip)
               enddo
               hphil = phi(ic0,ic1,ic2,nphi)*phi(ic0,ic1,ic2,nphi)
               suminv2 = 1.d0/(hphis+hphil)
               hphis = hphis*suminv2
               hphil = hphil*suminv2

               do ip = 1, nphi-1
                  rhs(ic0,ic1,ic2,ip) = rhs(ic0,ic1,ic2,ip) +
     &               2.d0*phi(ic0,ic1,ic2,ip)*hphil*suminv2
     &               *( fl(ic0,ic1,ic2) - fa(ic0,ic1,ic2) )
               enddo
               rhs(ic0,ic1,ic2,nphi) = rhs(ic0,ic1,ic2,nphi) +
     &               2.d0*phi(ic0,ic1,ic2,nphi)*hphis*suminv2
     &               *( fa(ic0,ic1,ic2) - fl(ic0,ic1,ic2) )

            enddo
         enddo
      enddo

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
     &   with_phase,
     &   phi_rhs, ngphi_rhs,
     &   rhs, ngrhs )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2

      double precision dx(0:2)
      double precision thermal_diffusivity, latent_heat
      integer ngphi_rhs, ngtemp, ngrhs, ngcp, with_phase
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
      double precision diag, offdiagx, offdiagy, offdiagz

      dxinv2 = 1.d0/(dx(0)*dx(0))
      dyinv2 = 1.d0/(dx(1)*dx(1))
      dzinv2 = 1.d0/(dx(2)*dx(2))

      offdiagx = thermal_diffusivity * dxinv2
      offdiagy = thermal_diffusivity * dyinv2
      offdiagz = thermal_diffusivity * dzinv2

      diag = -2.d0 * (offdiagx + offdiagy + offdiagz)
c
      call stencil7pts(ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &                 diag,offdiagx,offdiagy,offdiagz,
     &                 temp, ngtemp,rhs, ngrhs )

      if( with_phase /= 0 )then
         do ic2 = ifirst2, ilast2
            do ic1 = ifirst1, ilast1
               do ic0 = ifirst0, ilast0

                  gamma = latent_heat/cp(ic0,ic1,ic2)

                  rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2) +
     &               gamma * phi_rhs(ic0,ic1,ic2)

               enddo
            enddo
         enddo
      endif

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
     &   rhs, ngrhs,
     &   phi_interp_type )
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngphi, ngtemp, ngrhs
      double precision tm, latentheat
      character*(*) phi_interp_type
c
c variables in 3d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi))
      double precision rhs(CELL3d(ifirst,ilast,ngrhs))
      double precision temp(CELL3d(ifirst,ilast,ngtemp))
c
      integer ic0, ic1, ic2
      double precision m, alpha, wtemp, h_prime, woff
      double precision deriv_interp_func
c
      alpha = latentheat/tm
      woff = 0.25/6.
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               wtemp = 0.75*temp(ic0,ic1,ic2)
     &           +woff*( temp(ic0-1,ic1,ic2)
     &                  +temp(ic0,ic1-1,ic2)
     &                  +temp(ic0+1,ic1,ic2)
     &                  +temp(ic0,ic1+1,ic2)
     &                  +temp(ic0,ic1,ic2-1)
     &                  +temp(ic0,ic1,ic2+1) )

               m = alpha*( tm-wtemp )
               h_prime =
     &            deriv_interp_func( phi(ic0,ic1,ic2), phi_interp_type )

               rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2) + m*h_prime

            enddo
         enddo
      enddo

      return
      end

c***********************************************************************
c
      subroutine computedphidtemperaturedeltatemperature(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phi, ngphi,
     &   tm, latentheat, mobility,
     &   rhs, ngrhs,
     &   phi_interp_type )
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngphi, ngrhs
      double precision tm, latentheat, mobility
      character*(*) phi_interp_type
c
c variables in 3d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi))
      double precision rhs(CELL3d(ifirst,ilast,ngrhs))
c
      integer ic0, ic1, ic2
      double precision alpha, h_prime
      double precision deriv_interp_func
c
c      print*,'latentheat=',latentheat,', tm=',tm
      alpha = mobility*latentheat/tm

      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               h_prime =
     &            deriv_interp_func( phi(ic0,ic1,ic2), phi_interp_type )

               rhs(ic0,ic1, ic2) = alpha*h_prime

            enddo
         enddo
      enddo

      return
      end

c***********************************************************************
c
      subroutine addvdphidx(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   dx,
     &   phase, ngphase,
     &   vel,
     &   rhs, ngrhs)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2

      double precision dx(3), vel
      integer ngphase, ngrhs
c
c variables in 3d cell indexed
      double precision phase(CELL3d(ifirst,ilast,ngphase))
      double precision rhs(CELL3d(ifirst,ilast,ngrhs))

c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1, ic2
      double precision dxinv, diff_term_x

      dxinv = 0.5d0 * vel / dx(1)
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               diff_term_x = dxinv *
     &              (phase(ic0+1,ic1,ic2) - phase(ic0-1,ic1,ic2))

               rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2)
     &                      + diff_term_x

            enddo
         enddo
      enddo

      end
c***********************************************************************
c
      subroutine addvdphidx_upwind3(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   dx,
     &   phase, ngphase,
     &   vel,
     &   rhs, ngrhs,
     &   physb)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2

      double precision dx(3), vel
      integer ngphase, ngrhs
      integer physb
c
c variables in 3d cell indexed
      double precision phase(CELL3d(ifirst,ilast,ngphase))
      double precision rhs(CELL3d(ifirst,ilast,ngrhs))

c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1, ic2, imax0
      double precision dxinv, diff_term_x

      dxinv = vel / (6.d0*dx(1))
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            if( physb .eq. 1) then
c use regular 2nd order at right boundary since 2nd ghost value
c typically not filled with physical boundar value
               imax0 = ilast0-1
               diff_term_x = 3.d0*dxinv *
     &            (phase(ilast0+1,ic1,ic2)-phase(ilast0-1,ic1,ic2))
               rhs(ilast0,ic1,ic2) = rhs(ilast0,ic1,ic2) + diff_term_x
            else
               imax0 = ilast0
            endif
            do ic0 = ifirst0, imax0

               diff_term_x = dxinv *
     &           (-2.d0*phase(ic0-1,ic1,ic2) -3.d0*phase(ic0,ic1,ic2)
     &            +6.d0*phase(ic0+1,ic1,ic2) - phase(ic0+2,ic1,ic2))

               rhs(ic0,ic1,ic2) = rhs(ic0,ic1,ic2)
     &                      + diff_term_x

            enddo
         enddo
      enddo

      end
