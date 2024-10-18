c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
define(NDIM,2)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl

c
c Evaluate coefficient p2(phi)*eps_q^2+D_q(phi)/|nabla q|
c
      subroutine compute_face_coef2d(
     &     lo0, hi0, lo1, hi1,
     &     depth,
     &     eps_q,
     &     phi, ngp,
     &     temp, ngt,
     &     misorientation_factor,
     &     gqx, gqy, nggq,
     &     fcx, fcy, ngf,
     &     gradient_floor, floor_type,
     &     interp_type1, interp_type2,
     &     avg_type
     &     )
c
      implicit none
c
      integer lo0, hi0, lo1, hi1,
     &        depth
c variables in 2d cell indexed
      integer ngp, ngt, nggq, ngf
      double precision phi(CELL2d(lo,hi,ngp))
      double precision temp(CELL2d(lo,hi,ngt))
      double precision gqx(SIDE2d0(lo,hi,nggq),depth,NDIM)
      double precision gqy(SIDE2d1(lo,hi,nggq),depth,NDIM)
      double precision fcx(SIDE2d0(lo,hi,ngf))
      double precision fcy(SIDE2d1(lo,hi,ngf))

      character*(*) interp_type1, interp_type2
      character*(*) floor_type
      character*(*) avg_type

      double precision
     &           gradient_floor, eps_q

      double precision eval_grad_normi
      double precision misorientation_factor

c     local variables:
      integer i, j, m, n
      double precision grad_norm2, grad_normi,
     &     floor_grad_norm2,
     &     max_grad_normi, eps2, hphi2
      double precision phia, tempa, diff
      double precision interp_func
      double precision average_func

      floor_grad_norm2 = gradient_floor**2
      eps2 = eps_q*eps_q

      max_grad_normi = 1.d0 / gradient_floor

c     x faces

      do j = lo1, hi1
         do i = lo0, hi0+1

            phia = average_func(phi(i-1,j), phi(i,j),
     &                          avg_type)
            tempa = 0.5d0*(temp(i-1,j) + temp(i,j))
            diff = misorientation_factor * tempa 
     &           * interp_func(phia, interp_type1)
            hphi2 = interp_func(phia, interp_type2)

c           compute reciprocal of gradient L2 norm on this face
            grad_norm2 = 0.d0
c loop over components of nabla q_m
            do n = 1, NDIM
c loop over quaternion components
               do m = 1, depth
                  grad_norm2 = grad_norm2 + gqx(i,j,m,n)**2
               enddo
            enddo
            grad_normi = eval_grad_normi(grad_norm2, floor_type, 
     &                                   floor_grad_norm2, 
     &                                   max_grad_normi)
            fcx(i,j) = - grad_normi * diff - eps2 * hphi2

         enddo
      enddo

c     y faces

      do j = lo1, hi1+1
         do i = lo0, hi0

            phia = average_func(phi(i,j-1), phi(i,j),
     &                          avg_type)
            tempa = 0.5d0*(temp(i,j-1) + temp(i,j))
            diff = misorientation_factor * tempa
     &           * interp_func(phia, interp_type1)
            hphi2 = interp_func(phia, interp_type2)

c           compute reciprocal of gradient L2 norm on this face
            grad_norm2 = 0.d0
            do n = 1, NDIM
               do m = 1, depth
                  grad_norm2 = grad_norm2 + gqy(i,j,m,n)**2
               enddo
            enddo
            grad_normi = eval_grad_normi(grad_norm2, floor_type, 
     &                                   floor_grad_norm2, 
     &                                   max_grad_normi)
            fcy(i,j) = - grad_normi * diff - eps2 * hphi2

         enddo
      enddo

      return
      end

      subroutine compute_dquatdphi_face_coef2d(
     &     lo0, hi0, lo1, hi1,
     &     depth,
     &     dprimex, dpxlo0, dpxhi0, dpxlo1, dpxhi1,
     &     dprimey, dpylo0, dpyhi0, dpylo1, dpyhi1,
     &     phi, plo0, phi0, plo1, phi1,
     &     fcx, fcxlo0, fcxhi0, fcxlo1, fcxhi1,
     &     fcy, fcylo0, fcyhi0, fcylo1, fcyhi1
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1,
     &        depth,
     &        dpxlo0, dpxhi0, dpxlo1, dpxhi1,
     &        dpylo0, dpyhi0, dpylo1, dpyhi1,
     &        plo0, phi0, plo1, phi1,
     &        fcxlo0, fcxhi0, fcxlo1, fcxhi1,
     &        fcylo0, fcyhi0, fcylo1, fcyhi1
      double precision
     &           dprimex(dpxlo0:dpxhi0,dpxlo1:dpxhi1,2),
     &           dprimey(dpylo0:dpyhi0,dpylo1:dpyhi1,2),
     &           phi(plo0:phi0,plo1:phi1),
     &           fcx(fcxlo0:fcxhi0,fcxlo1:fcxhi1),
     &           fcy(fcylo0:fcyhi0,fcylo1:fcyhi1)

c     local variables:
      integer i, j

c     x faces

      do j = lo1, hi1
         do i = lo0, hi0+1
            fcx(i,j) = - dprimex(i,j,1) * phi(i-1,j  )
     &                 - dprimex(i,j,2) * phi(i  ,j  )
         enddo
      enddo

c     y faces

      do j = lo1, hi1+1
         do i = lo0, hi0
            fcy(i,j) = - dprimey(i,j,1) * phi(i  ,j-1)
     &                 - dprimey(i,j,2) * phi(i  ,j  )
         enddo
      enddo

      return
      end

      subroutine compute_flux2d(
     &     lo0, hi0, lo1, hi1,
     &     depth,
     &     fcx, fcxlo0, fcxhi0, fcxlo1, fcxhi1,
     &     fcy, fcylo0, fcyhi0, fcylo1, fcyhi1,
     &     q, qlo0, qhi0, qlo1, qhi1,
     &     h,
     &     fx, fxlo0, fxhi0, fxlo1, fxhi1,
     &     fy, fylo0, fyhi0, fylo1, fyhi1
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1,
     &        depth,
     &        fcxlo0, fcxhi0, fcxlo1, fcxhi1,
     &        fcylo0, fcyhi0, fcylo1, fcyhi1,
     &        qlo0, qhi0, qlo1, qhi1,
     &        fxlo0, fxhi0, fxlo1, fxhi1,
     &        fylo0, fyhi0, fylo1, fyhi1
      double precision fcx(fcxlo0:fcxhi0,fcxlo1:fcxhi1),
     &                 fcy(fcylo0:fcyhi0,fcylo1:fcyhi1),
     &                 q(qlo0:qhi0,qlo1:qhi1,depth),
     &                 h(2),
     &                 fx(fxlo0:fxhi0,fxlo1:fxhi1,depth),
     &                 fy(fylo0:fyhi0,fylo1:fyhi1,depth)

c     local variables:
      integer i, j, m
      double precision hinv

c     x faces

      hinv = 1.d0 / h(1)

      do m = 1, depth
         do j = lo1, hi1
            do i = lo0, hi0+1
               fx(i,j,m) = fcx(i,j) * hinv *
     &              (q(i,j,m) - q(i-1,j,m))
            enddo
         enddo
      enddo

c     y faces

      hinv = 1.d0 / h(2)

      do m = 1, depth
         do j = lo1, hi1+1
            do i = lo0, hi0
               fy(i,j,m) = fcy(i,j) * hinv *
     &              (q(i,j,m) - q(i,j-1,m))
            enddo
         enddo
      enddo

      return
      end

      subroutine compute_flux2d_from_gradq(
     &     lo0, hi0, lo1, hi1,
     &     depth,
     &     fcx, fcxlo0, fcxhi0, fcxlo1, fcxhi1,
     &     fcy, fcylo0, fcyhi0, fcylo1, fcyhi1,
     &     grad_x_xside,
     &     grad_y_yside,
     &     fx, fxlo0, fxhi0, fxlo1, fxhi1,
     &     fy, fylo0, fyhi0, fylo1, fyhi1
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1,
     &        depth,
     &        fcxlo0, fcxhi0, fcxlo1, fcxhi1,
     &        fcylo0, fcyhi0, fcylo1, fcyhi1,
     &        fxlo0, fxhi0, fxlo1, fxhi1,
     &        fylo0, fyhi0, fylo1, fyhi1
      double precision fcx(fcxlo0:fcxhi0,fcxlo1:fcxhi1),
     &                 fcy(fcylo0:fcyhi0,fcylo1:fcyhi1),
     &                 fx(fxlo0:fxhi0,fxlo1:fxhi1,depth),
     &                 fy(fylo0:fyhi0,fylo1:fyhi1,depth)

      double precision grad_x_xside(SIDE2d0(lo,hi,0),depth)
      double precision grad_y_yside(SIDE2d1(lo,hi,0),depth)

c     local variables:
      integer i, j, m

c     x faces

      do m = 1, depth
         do j = lo1, hi1
            do i = lo0, hi0+1
               fx(i,j,m) = fcx(i,j) * grad_x_xside(i,j,m)
            enddo
         enddo
      enddo

c     y faces

      do m = 1, depth
         do j = lo1, hi1+1
            do i = lo0, hi0
               fy(i,j,m) = fcy(i,j) * grad_y_yside(i,j,m)
            enddo
         enddo
      enddo

      return
      end

      subroutine compute_sym_flux2d(
     &     lo0, hi0, lo1, hi1,
     &     depth,
     &     fcx, fcxlo0, fcxhi0, fcxlo1, fcxhi1,
     &     fcy, fcylo0, fcyhi0, fcylo1, fcyhi1,
     &     sqrt_m, mlo0, mhi0, mlo1, mhi1,
     &     q, qlo0, qhi0, qlo1, qhi1,
     &     h,
     &     fx, fxlo0, fxhi0, fxlo1, fxhi1,
     &     fy, fylo0, fyhi0, fylo1, fyhi1
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1,
     &        depth,
     &        fcxlo0, fcxhi0, fcxlo1, fcxhi1,
     &        fcylo0, fcyhi0, fcylo1, fcyhi1,
     &        mlo0, mhi0, mlo1, mhi1,
     &        qlo0, qhi0, qlo1, qhi1,
     &        fxlo0, fxhi0, fxlo1, fxhi1,
     &        fylo0, fyhi0, fylo1, fyhi1
      double precision fcx(fcxlo0:fcxhi0,fcxlo1:fcxhi1),
     &                 fcy(fcylo0:fcyhi0,fcylo1:fcyhi1),
     &                 sqrt_m(mlo0:mhi0,mlo1:mhi1),
     &                 q(qlo0:qhi0,qlo1:qhi1,depth),
     &                 h(2),
     &                 fx(fxlo0:fxhi0,fxlo1:fxhi1,depth),
     &                 fy(fylo0:fyhi0,fylo1:fyhi1,depth)

c     local variables:
      integer i, j, m
      double precision hinv

c     x faces

      hinv = 1.d0 / h(1)

      do j = lo1, hi1
         do i = lo0, hi0+1
            do m = 1, depth
               fx(i,j,m) = fcx(i,j) * hinv *
     &              (sqrt_m(i,j)*q(i,j,m) - sqrt_m(i-1,j)*q(i-1,j,m))
            enddo
         enddo
      enddo

c     y faces

      hinv = 1.d0 / h(2)

      do j = lo1, hi1+1
         do i = lo0, hi0
            do m = 1, depth
               fy(i,j,m) = fcy(i,j) * hinv *
     &              (sqrt_m(i,j)*q(i,j,m) - sqrt_m(i,j-1)*q(i,j-1,m))
            enddo
         enddo
      enddo

      return
      end

      subroutine compute_q_residual2d(
     &     lo0, hi0, lo1, hi1,
     &     depth,
     &     sqrt_m, mlo0, mhi0, mlo1, mhi1,
     &     fx, fxlo0, fxhi0, fxlo1, fxhi1,
     &     fy, fylo0, fyhi0, fylo1, fyhi1,
     &     q, qlo0, qhi0, qlo1, qhi1,
     &     h, gamma,
     &     rhs, rhlo0, rhhi0, rhlo1, rhhi1,
     &     residual, rlo0, rhi0, rlo1, rhi1
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1,
     &        depth,
     &        mlo0, mhi0, mlo1, mhi1,
     &        fxlo0, fxhi0, fxlo1, fxhi1,
     &        fylo0, fyhi0, fylo1, fyhi1,
     &        qlo0, qhi0, qlo1, qhi1,
     &        rhlo0, rhhi0, rhlo1, rhhi1,
     &        rlo0, rhi0, rlo1, rhi1
      double precision 
     &        sqrt_m(mlo0:mhi0,mlo1:mhi1),
     &        fx(fxlo0:fxhi0,fxlo1:fxhi1,depth),
     &        fy(fylo0:fyhi0,fylo1:fyhi1,depth),
     &        q(qlo0:qhi0,qlo1:qhi1,depth),
     &        h(2), gamma,
     &        rhs(rhlo0:rhhi0,rhlo1:rhhi1,depth),
     &        residual(rlo0:rhi0,rlo1:rhi1,depth)

c     local variables:
      integer i, j, m
      double precision dxinv, dyinv, right, left, up, down,
     &                 divergence

c      print*,'compute_q_residual2d'

      dxinv = 1.d0 / h(1)
      dyinv = 1.d0 / h(2)

      do m = 1, depth
         do j = lo1, hi1
            do i = lo0, hi0
                  
               right = fx(i+1,j,m)
               left  = fx(i  ,j,m)
               up    = fy(i,j+1,m)
               down  = fy(i,j  ,m)

               divergence = (right - left) * dxinv
     &              + (up - down)    * dyinv

               residual(i,j,m) = rhs(i,j,m) - q(i,j,m)
     &              - gamma * sqrt_m(i,j) * divergence

            enddo
         enddo
      enddo

      return
      end

      subroutine compute_q_residual2d_symm(
     &     lo0, hi0, lo1, hi1,
     &     depth,
     &     sqrt_m, mlo0, mhi0, mlo1, mhi1,
     &     fx, fxlo0, fxhi0, fxlo1, fxhi1,
     &     fy, fylo0, fyhi0, fylo1, fyhi1,
     &     q, qlo0, qhi0, qlo1, qhi1,
     &     h, gamma,
     &     rhs, rhlo0, rhhi0, rhlo1, rhhi1,
     &     residual, rlo0, rhi0, rlo1, rhi1,
     &     iqrot_x, iqrot_y, ngiq
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1,
     &        depth, ngiq,
     &        mlo0, mhi0, mlo1, mhi1,
     &        fxlo0, fxhi0, fxlo1, fxhi1,
     &        fylo0, fyhi0, fylo1, fyhi1,
     &        qlo0, qhi0, qlo1, qhi1,
     &        rhlo0, rhhi0, rhlo1, rhhi1,
     &        rlo0, rhi0, rlo1, rhi1
      double precision 
     &        sqrt_m(mlo0:mhi0,mlo1:mhi1),
     &        fx(fxlo0:fxhi0,fxlo1:fxhi1,depth),
     &        fy(fylo0:fyhi0,fylo1:fyhi1,depth),
     &        q(qlo0:qhi0,qlo1:qhi1,depth),
     &        h(2), gamma,
     &        rhs(rhlo0:rhhi0,rhlo1:rhhi1,depth),
     &        residual(rlo0:rhi0,rlo1:rhi1,depth)
      integer iqrot_x(SIDE2d0(lo,hi,ngiq))
      integer iqrot_y(SIDE2d1(lo,hi,ngiq))

c     local variables:
      integer i, j, m, iq
      double precision dxinv, dyinv, 
     &                   right(depth), left(depth), 
     &                   up(depth), down(depth),
     &                   divergence(depth)
      double precision dtmp(depth)

c      print*,'compute_q_residual2d_symm'

      dxinv = 1.d0 / h(1)
      dyinv = 1.d0 / h(2)

      do j = lo1, hi1
         do i = lo0, hi0
                  
            do m = 1, depth
               dtmp(m) = fx(i+1,j,m)
               left(m) = fx(i  ,j,m)
            enddo
            iq = -1 * iqrot_x(i+1,j)
            call quatsymmrotate( dtmp, iq, right, depth )

            do m = 1, depth
               dtmp(m) = fy(i,j+1,m)
               down(m) = fy(i,j  ,m)
            enddo
            iq = -1 * iqrot_y(i,j+1)
            call quatsymmrotate( dtmp, iq, up, depth )
            
            do m = 1, depth
               divergence(m) = (right(m) - left(m)) * dxinv
     &                       + (up(m)    - down(m)) * dyinv
            enddo

            do m = 1, depth
               residual(i,j,m) = rhs(i,j,m) - q(i,j,m)
     &              - gamma * sqrt_m(i,j) * divergence(m)

            enddo
         enddo
      enddo

      return
      end

      subroutine add_quat_op2d(
     &     lo0, hi0, lo1, hi1,
     &     depth,
     &     mobility, mlo0, mhi0, mlo1, mhi1,
     &     fx, fxlo0, fxhi0, fxlo1, fxhi1,
     &     fy, fylo0, fyhi0, fylo1, fyhi1,
     &     h,
     &     rhs, rhlo0, rhhi0, rhlo1, rhhi1
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1,
     &        depth,
     &        mlo0, mhi0, mlo1, mhi1,
     &        fxlo0, fxhi0, fxlo1, fxhi1,
     &        fylo0, fyhi0, fylo1, fyhi1,
     &        rhlo0, rhhi0, rhlo1, rhhi1
      double precision 
     &        mobility(mlo0:mhi0,mlo1:mhi1),
     &        fx(fxlo0:fxhi0,fxlo1:fxhi1,depth),
     &        fy(fylo0:fyhi0,fylo1:fyhi1,depth),
     &        h(2),
     &        rhs(rhlo0:rhhi0,rhlo1:rhhi1,depth)

c     local variables:
      integer i, j, m
      double precision dxinv, dyinv, right, left, up, down,
     &                 divergence

c      print*,'add_quat_op2d'

      dxinv = 1.d0 / h(1)
      dyinv = 1.d0 / h(2)

      do m = 1, depth
         do j = lo1, hi1
            do i = lo0, hi0

               right = fx(i+1,j,m)
               left  = fx(i  ,j,m)
               up    = fy(i,j+1,m)
               down  = fy(i,j  ,m)

               divergence = (right - left) * dxinv
     &                    + (up - down)    * dyinv

               rhs(i,j,m) = rhs(i,j,m)
     &              - mobility(i,j) * divergence
            enddo
         enddo
      enddo

      return
      end

      subroutine add_quat_proj_op2d(
     &     lo0, hi0, lo1, hi1,
     &     depth,
     &     mobility, mlo0, mhi0, mlo1, mhi1,
     &     fx, fxlo0, fxhi0, fxlo1, fxhi1,
     &     fy, fylo0, fyhi0, fylo1, fyhi1,
     &     q, qlo0, qhi0, qlo1, qhi1,
     &     lambda, llo0, lhi0, llo1, lhi1,
     &     h,
     &     rhs, rhlo0, rhhi0, rhlo1, rhhi1
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1,
     &        depth,
     &        mlo0, mhi0, mlo1, mhi1,
     &        fxlo0, fxhi0, fxlo1, fxhi1,
     &        fylo0, fyhi0, fylo1, fyhi1,
     &        qlo0, qhi0, qlo1, qhi1,
     &        llo0, lhi0, llo1, lhi1,
     &        rhlo0, rhhi0, rhlo1, rhhi1
      double precision 
     &        mobility(mlo0:mhi0,mlo1:mhi1),
     &        fx(fxlo0:fxhi0,fxlo1:fxhi1,depth),
     &        fy(fylo0:fyhi0,fylo1:fyhi1,depth),
     &        q(qlo0:qhi0,qlo1:qhi1,depth),
     &        lambda(llo0:lhi0,llo1:lhi1),
     &        h(2),
     &        rhs(rhlo0:rhhi0,rhlo1:rhhi1,depth)

c     local variables:
      integer i, j, m
      double precision dxinv, dyinv, right, left, up, down,
     &                 divergence

c      print*,'add_quat_proj_op2d'

      dxinv = 1.d0 / h(1)
      dyinv = 1.d0 / h(2)

      do m = 1, depth
         do j = lo1, hi1
            do i = lo0, hi0

               right = fx(i+1,j,m)
               left  = fx(i  ,j,m)
               up    = fy(i,j+1,m)
               down  = fy(i,j  ,m)

               divergence = (right - left) * dxinv
     &                    + (up - down)    * dyinv

               rhs(i,j,m) = rhs(i,j,m)
     &              - mobility(i,j) * (divergence
     &              + 2.d0 * q(i,j,m) * lambda(i,j))
            enddo
         enddo
      enddo

      return
      end

       subroutine add_quat_proj_op2d_symm(
     &     lo0, hi0, lo1, hi1,
     &     depth,
     &     mobility, mlo0, mhi0, mlo1, mhi1,
     &     fx, fxlo0, fxhi0, fxlo1, fxhi1,
     &     fy, fylo0, fyhi0, fylo1, fyhi1,
     &     q, qlo0, qhi0, qlo1, qhi1,
     &     lambda, llo0, lhi0, llo1, lhi1,
     &     h,
     &     rhs, rhlo0, rhhi0, rhlo1, rhhi1,
     &     iqrot_x, iqrot_y, ngiq
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1,
     &        depth, ngiq,
     &        mlo0, mhi0, mlo1, mhi1,
     &        fxlo0, fxhi0, fxlo1, fxhi1,
     &        fylo0, fyhi0, fylo1, fyhi1,
     &        qlo0, qhi0, qlo1, qhi1,
     &        llo0, lhi0, llo1, lhi1,
     &        rhlo0, rhhi0, rhlo1, rhhi1
      double precision 
     &        mobility(mlo0:mhi0,mlo1:mhi1),
     &        fx(fxlo0:fxhi0,fxlo1:fxhi1,depth),
     &        fy(fylo0:fyhi0,fylo1:fyhi1,depth),
     &        q(qlo0:qhi0,qlo1:qhi1,depth),
     &        lambda(llo0:lhi0,llo1:lhi1),
     &        h(2),
     &        rhs(rhlo0:rhhi0,rhlo1:rhhi1,depth)
      integer iqrot_x(SIDE2d0(lo,hi,ngiq))
      integer iqrot_y(SIDE2d1(lo,hi,ngiq))

c     local variables:
      integer i, j, m, iq
      double precision dxinv, dyinv, 
     &                 right(depth), left(depth), 
     &                 up(depth), down(depth),
     &                 divergence(depth)
      double precision dtmp(depth)

c      print*,'add_quat_proj_op2d_symm'

      dxinv = 1.d0 / h(1)
      dyinv = 1.d0 / h(2)

      do j = lo1, hi1
         do i = lo0, hi0

            do m = 1, depth
               dtmp(m) = fx(i+1,j,m)
               left(m) = fx(i  ,j,m)
            enddo
            iq = -1 * iqrot_x(i+1,j)
            call quatsymmrotate( dtmp, iq, right, depth )

            do m = 1, depth
               dtmp(m) = fy(i,j+1,m)
               down(m) = fy(i,j  ,m)
            enddo
            iq = -1 * iqrot_y(i,j+1)
            call quatsymmrotate( dtmp, iq, up, depth )
            
            do m = 1, depth
               divergence(m) = (right(m) - left(m)) * dxinv
     &                       + (up(m)    - down(m)) * dyinv
            enddo

            do m = 1, depth
               rhs(i,j,m) = rhs(i,j,m)
     &              - mobility(i,j) * (divergence(m)
     &              + 2.d0 * q(i,j,m) * lambda(i,j))
            enddo
         enddo
      enddo

      return
      end

      subroutine compute_lambda_flux2d(
     &     lo0, hi0, lo1, hi1,
     &     depth,
     &     fx, fxlo0, fxhi0, fxlo1, fxhi1,
     &     fy, fylo0, fyhi0, fylo1, fyhi1,
     &     q, qlo0, qhi0, qlo1, qhi1,
     &     h,
     &     lambda, llo0, lhi0, llo1, lhi1
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1,
     &        depth,
     &        fxlo0, fxhi0, fxlo1, fxhi1,
     &        fylo0, fyhi0, fylo1, fyhi1,
     &        qlo0, qhi0, qlo1, qhi1,
     &        llo0, lhi0, llo1, lhi1
      double precision fx(fxlo0:fxhi0,fxlo1:fxhi1,depth),
     &                 fy(fylo0:fyhi0,fylo1:fyhi1,depth),
     &                 q(qlo0:qhi0,qlo1:qhi1,depth),
     &                 h(NDIM),
     &                 lambda(llo0:lhi0,llo1:lhi1)

c     local variables:
      integer i, j, m
      double precision xfac, yfac, left, right, down, up, sumq2

c      print*,'compute_lambda_flux2d'
      
      xfac = 0.5d0 / h(1)
      yfac = 0.5d0 / h(2)

      do j = lo1, hi1
         do i = lo0, hi0

            lambda(i,j) = 0.d0
            sumq2 = 0.d0
            do m = 1, depth

               left  = fx(i,j,m)
               right = fx(i+1,j,m)
               down  = fy(i,j,m)
               up    = fy(i,j+1,m)

               lambda(i,j) = lambda(i,j) - q(i,j,m) *
     &                 ((right - left) * xfac + (up - down) * yfac)
               sumq2 = sumq2 + q(i,j,m)**2
            enddo
            lambda(i,j) = lambda(i,j) / sumq2
         enddo
      enddo

      return
      end

      subroutine compute_lambda_flux2d_symm(
     &     lo0, hi0, lo1, hi1,
     &     depth,
     &     fx, fxlo0, fxhi0, fxlo1, fxhi1,
     &     fy, fylo0, fyhi0, fylo1, fyhi1,
     &     q, qlo0, qhi0, qlo1, qhi1,
     &     h,
     &     lambda, llo0, lhi0, llo1, lhi1,
     &     iqrot_x, iqrot_y, ngiq
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1,
     &        depth, ngiq,
     &        fxlo0, fxhi0, fxlo1, fxhi1,
     &        fylo0, fyhi0, fylo1, fyhi1,
     &        qlo0, qhi0, qlo1, qhi1,
     &        llo0, lhi0, llo1, lhi1
      double precision fx(fxlo0:fxhi0,fxlo1:fxhi1,depth),
     &                 fy(fylo0:fyhi0,fylo1:fyhi1,depth),
     &                 q(qlo0:qhi0,qlo1:qhi1,depth),
     &                 h(NDIM),
     &                 lambda(llo0:lhi0,llo1:lhi1)
      integer iqrot_x(SIDE2d0(lo,hi,ngiq))
      integer iqrot_y(SIDE2d1(lo,hi,ngiq))

c     local variables:
      integer i, j, m, iq
      double precision xfac, yfac, sumq2
      double precision right(depth), left(depth), 
     &                 up(depth), down(depth)
      double precision dtmp(depth)

c      print*,'compute_lambda_flux2d_symm'
      
      xfac = 0.5d0 / h(1)
      yfac = 0.5d0 / h(2)

      do j = lo1, hi1
         do i = lo0, hi0

            do m = 1, depth
               dtmp(m) = fx(i+1,j,m)
               left(m) = fx(i  ,j,m)
            enddo
            iq = -1 * iqrot_x(i+1,j)
            call quatsymmrotate( dtmp, iq, right, depth )

            do m = 1, depth
               dtmp(m) = fy(i,j+1,m)
               down(m) = fy(i,j  ,m)
            enddo
            iq = -1 * iqrot_y(i,j+1)
            call quatsymmrotate( dtmp, iq, up, depth )
            
            lambda(i,j) = 0.d0
            sumq2 = 0.d0
            do m = 1, depth

               lambda(i,j) = lambda(i,j) - q(i,j,m) *
     &                 ((right(m) - left(m)) * xfac 
     &                + (up(m) - down(m)) * yfac)
               sumq2 = sumq2 + q(i,j,m)**2
            enddo
            lambda(i,j) = lambda(i,j) / sumq2
         enddo
      enddo

      return
      end

      subroutine fixflux2d(
     &     xflux, yflux, fluxgi, fluxgj, 
     &     fluxnc,
     &     xdc, ydc, dcgi, dcgj,
     &     dcnc,
     &     soln, solngi, solngj,
     &     solnnc,
     &     ifirst, ilast, jfirst, jlast,
     &     nc,
     &     location_index,
     &     ratio_to_coarser,
     &     blower, bupper,
     &     dx)
c
      implicit none
      integer
     &     ifirst, ilast, jfirst, jlast, nc, fluxnc,
     &     dcgi, dcgj, fluxgi, fluxgj, dcnc,
     &     solngi, solngj, solnnc
      double precision
     &     xflux(ifirst-fluxgi:ilast+1+fluxgi,
     &           jfirst-fluxgj:jlast+fluxgj,fluxnc),
     &     yflux(ifirst-fluxgi:ilast+fluxgi,
     &           jfirst-fluxgj:jlast+1+fluxgj,fluxnc),
     &     xdc(ifirst-dcgi:ilast+1+dcgi,
     &         jfirst-dcgj:jlast+dcgj,dcnc),
     &     ydc(ifirst-dcgi:ilast+dcgi,
     &         jfirst-dcgj:jlast+1+dcgj,dcnc),
     &     soln(ifirst-solngi:ilast+solngi,
     &          jfirst-solngj:jlast+solngj, solnnc),
     &     dx(0:1)
      integer location_index,
     &     ratio_to_coarser(0:1)
c     Lower and upper corners of boundary box
      integer blower(0:1), bupper(0:1)

      double precision h
      integer i, ibeg, iend, igho, j, jbeg, jend, jgho, comp
c     Fine grid indices inside one coarse grid.
      integer ip, jp
c     Fine grid indices for point diametrically opposite from (ip,jp).
      integer iq, jq
c     Weights associated with longtitudinal and transverse
c     (with respect to boundary normal) gradients.
      double precision tranwt, longwt

      if ( location_index .eq. 0 ) then
c        min i edge
         tranwt = 1.0/(1+ratio_to_coarser(0))
         longwt = 2*tranwt
         h = dx(0)
         igho = bupper(0)
         i = igho+1
         jbeg = max(blower(1),jfirst)
         jend = min(bupper(1),jlast)
         do comp = 1, nc
            do j=jbeg,jend,ratio_to_coarser(1)
               do jp=0,ratio_to_coarser(1)-1
                  jq = ratio_to_coarser(1) - jp - 1
                  xflux(i,j+jp,comp)
     &                 = longwt*xflux(i,j+jp,comp)
     &                 + tranwt*xdc(i,j+jp,comp)*(
     &                 soln(i,j+jq,comp) - soln(i,j+jp,comp) )/h
               enddo
            enddo
         enddo
      elseif ( location_index .eq. 1 ) then
c        max i edge
         tranwt = 1.0/(1+ratio_to_coarser(0))
         longwt = 2*tranwt
         h = dx(0)
         igho = blower(0)
         i = igho-1
         jbeg = max(blower(1),jfirst)
         jend = min(bupper(1),jlast)
         do comp = 1, nc
            do j=jbeg,jend,ratio_to_coarser(1)
               do jp=0,ratio_to_coarser(1)-1
                  jq = ratio_to_coarser(1) - jp - 1
                  xflux(igho,j+jp,comp)
     &                 = longwt*xflux(igho,j+jp,comp)
     &                 - tranwt*xdc(igho,j+jp,comp)*(
     &                 soln(i,j+jq,comp) - soln(i,j+jp,comp) )/h
               enddo
            enddo
         enddo
      elseif ( location_index .eq. 2 ) then
c        min j edge
         tranwt = 1.0/(1+ratio_to_coarser(1))
         longwt = 2*tranwt
         h = dx(1)
         jgho = bupper(1)
         j = jgho+1
         ibeg = max(blower(0),ifirst)
         iend = min(bupper(0),ilast)
         do comp = 1, nc
            do i=ibeg,iend,ratio_to_coarser(0)
               do ip=0,ratio_to_coarser(0)-1
                  iq = ratio_to_coarser(0) - ip - 1
                  yflux(i+ip,j,comp)
     &                 = longwt*yflux(i+ip,j,comp)
     &                 + tranwt*ydc(i+ip,j,comp)*(
     &                 soln(i+iq,j,comp) - soln(i+ip,j,comp) )/h
               enddo
            enddo
         enddo
      elseif ( location_index .eq. 3 ) then
c        max j edge
         tranwt = 1.0/(1+ratio_to_coarser(1))
         longwt = 2*tranwt
         h = dx(1)
         jgho = blower(1)
         j = jgho-1
         ibeg = max(blower(0),ifirst)
         iend = min(bupper(0),ilast)
         do comp = 1, nc
            do i=ibeg,iend,ratio_to_coarser(0)
               do ip=0,ratio_to_coarser(0)-1
                  iq = ratio_to_coarser(0) - ip - 1
                  yflux(i+ip,jgho,comp)
     &                 = longwt*yflux(i+ip,jgho,comp)
     &                 - tranwt*ydc(i+ip,jgho,comp)*(
     &                 soln(i+iq,j,comp) - soln(i+ip,j,comp) )/h
               enddo
            enddo
         enddo
      endif

      return
      end

      subroutine project2d(
     &     lo0, hi0, lo1, hi1,
     &     depth,
     &     q, qlo0, qhi0, qlo1, qhi1,
     &     corr, clo0, chi0, clo1, chi1,
     &     err, elo0, ehi0, elo1, ehi1
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1,
     &        depth,
     &        qlo0, qhi0, qlo1, qhi1,
     &        clo0, chi0, clo1, chi1,
     &        elo0, ehi0, elo1, ehi1
      double precision 
     &        q(qlo0:qhi0,qlo1:qhi1,depth),
     &        corr(clo0:chi0,clo1:chi1,depth),
     &        err(elo0:ehi0,elo1:ehi1,depth)

c     local variables:
      double precision fac
      integer i, j, m

      do j = lo1, hi1
         do i = lo0, hi0

            fac = 0.d0
            do m = 1, depth
               fac = fac + q(i,j,m)**2
            enddo
            fac = 1.d0 / dsqrt(fac)

c           Store the q projection in the correction array for now
            do m = 1, depth
               corr(i,j,m) = q(i,j,m) * fac
            enddo

c           Compute the dot product of the error with the projected q
            fac = 0.d0
            do m = 1, depth
               fac = fac + corr(i,j,m) * err(i,j,m)
            enddo

c           Subtract the error component in the q direction
            do m = 1, depth
               err(i,j,m) = err(i,j,m) - corr(i,j,m) * fac
            enddo

c           Finalize the correction: q + corr is on the constraint
            do m = 1, depth
               corr(i,j,m) = corr(i,j,m) - q(i,j,m)
            enddo

         enddo
      enddo

      return
      end

      subroutine take_square_root2d(
     &     lo0, hi0, lo1, hi1,
     &     var, vlo0, vhi0, vlo1, vhi1
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1,
     &        vlo0, vhi0, vlo1, vhi1
      double precision 
     &        var(vlo0:vhi0,vlo1:vhi1)

c     local variables:
      integer i, j

      do j = lo1, hi1
         do i = lo0, hi0
            var(i,j) = dsqrt(var(i,j))
         enddo
      enddo

      return
      end

      subroutine multicomponent_multiply2d(
     &     lo0, hi0, lo1, hi1,
     &     factor, flo0, fhi0, flo1, fhi1,
     &     var, vlo0, vhi0, vlo1, vhi1, vnc
     &     )
c
      implicit none
      integer
     &     lo0, hi0, lo1, hi1,
     &     flo0, fhi0, flo1, fhi1,
     &     vlo0, vhi0, vlo1, vhi1, vnc
      double precision 
     &     factor(flo0:fhi0,flo1:fhi1),
     &     var(vlo0:vhi0,vlo1:vhi1,vnc)

c     local variables:
      integer i, j, n

      do n = 1, vnc
         do j = lo1, hi1
            do i = lo0, hi0
               var(i,j,n) = var(i,j,n) * factor(i,j)
            enddo
         enddo
      enddo

      return
      end

      subroutine multicomponent_divide2d(
     &     lo0, hi0, lo1, hi1,
     &     factor, flo0, fhi0, flo1, fhi1,
     &     var, vlo0, vhi0, vlo1, vhi1, vnc
     &     )
c
      implicit none
      integer
     &     lo0, hi0, lo1, hi1,
     &     flo0, fhi0, flo1, fhi1,
     &     vlo0, vhi0, vlo1, vhi1, vnc
      double precision 
     &     factor(flo0:fhi0,flo1:fhi1),
     &     var(vlo0:vhi0,vlo1:vhi1,vnc)

c     local variables:
      integer i, j, n

      do n = 1, vnc
         do j = lo1, hi1
            do i = lo0, hi0
               var(i,j,n) = var(i,j,n) / factor(i,j)
            enddo
         enddo
      enddo

      return
      end
