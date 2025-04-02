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
c input: quat,n
c output: n4
c
      subroutine compute_n4(quat,n,n4)

      implicit none
      double precision quat(4),n(4),n4

c     local variables
      double precision qp(4),qtmp(4),np(4)

      call quatconj(quat,qp)

c rotation applied to np
      call quatmult4(n,qp,qtmp)
      call quatmult4(quat,qtmp,np)
c      print 100,n(2),n(3),n(4),np(1),np(2),np(3),np(4)
c100   format (7F6.3)

      n4=np(2)**4+np(3)**4+np(4)**4

      return
      end
c
c
c
      subroutine interface_anisotropic_energy(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   dx,
     &   phi, pghosts, nphases,
     &   quat, ngq, qlen,
     &   epsilon_phi,
     &   eps4, knumber,
     &   weight,
     &   phi_e,
     &   energy,
     &   eval_per_cell
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &    pghosts, ngq, qlen, knumber, nphases
      integer eval_per_cell

      double precision phi_e, epsilon_phi
      double precision phi(CELL3d(lo,hi,pghosts),nphases)
      double precision quat(CELL3d(lo,hi,ngq),qlen)
      double precision energy(CELL3d(lo,hi,0))
      double precision weight(CELL3d(lo,hi,0))

      double precision dx(NDIM)
      double precision eps4, nni

      double precision dxinv, dyinv, dzinv, theta
      double precision epstheta, diff_term
      double precision dphidx, dphidy, dphidz
      double precision e, q(4), n(4), n4, factor
      double precision threshold

      integer i, j, k, p
c
      phi_e = 0.d0

      dxinv = 0.5d0 / dx(1)
      dyinv = 0.5d0 / dx(2)
      dzinv = 0.5d0 / dx(3)
c
      factor = 4.d0*eps4/(1.d0-3.d0*eps4)
      threshold = 1.e-12
c
      phi_e = 0.d0
      if ( eval_per_cell /= 0 ) then
         do k = lo2, hi2
            do j = lo1, hi1
               do i = lo0, hi0
                  energy(i,j,k) = 0.d0
               enddo
            enddo
         enddo
      endif
c
      do p = 1, nphases
         do k = lo2, hi2
            do j = lo1, hi1
               do i = lo0, hi0
                  dphidx = (phi(i+1,j,k,p) - phi(i-1,j,k,p)) * dxinv
                  dphidy = (phi(i,j+1,k,p) - phi(i,j-1,k,p)) * dyinv
                  dphidz = (phi(i,j,k+1,p) - phi(i,j,k-1,p)) * dzinv

                  diff_term = dphidx*dphidx + dphidy*dphidy
     &                      + dphidz*dphidz

                  if( abs(diff_term)>threshold )then
                     nni=1./sqrt(diff_term)
                     n(1)=0.
                     n(2)=dphidx*nni
                     n(3)=dphidy*nni
                     n(4)=dphidz*nni

                     q(1)=quat(i,j,k,1)
                     q(2)=quat(i,j,k,2)
                     q(3)=quat(i,j,k,3)
                     q(4)=quat(i,j,k,4)

                     call compute_n4(q,n,n4)
                  else
                     n4 = 0.d0
                  endif

                  epstheta=epsilon_phi*(1.d0-3.d0*eps4)
     &                    *(1.d0+factor*n4)

                  e = 0.5d0 * epstheta * epstheta * diff_term
                  if ( eval_per_cell /= 0 ) then
                     energy(i,j,k) = energy(i,j,k) + e
                  endif

                  e = e * weight(i,j,k)
                  phi_e = phi_e + e
               enddo
            enddo
         enddo
      enddo

      return
      end
c
c
c
      subroutine phi_interface_energy(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   dx,
     &   phi, pghosts, nphases,
     &   epsilon_phi,
     &   weight,
     &   phi_e,
     &   energy,
     &   eval_per_cell
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &    pghosts, nphases
      integer eval_per_cell

      double precision phi_e, epsilon_phi, e
      double precision phi(CELL3d(lo,hi,pghosts),nphases)
      double precision energy(CELL3d(lo,hi,0))
      double precision weight(CELL3d(lo,hi,0))

      double precision dx(0:NDIM-1), dx2inv, dy2inv, dz2inv
      double precision diff_term_x, diff_term_y, diff_term_z
      double precision diff_term
      double precision epsilonphi2

      integer i, j, k, p

      phi_e = 0.d0
      if ( eval_per_cell /= 0 ) then
         do k = lo2, hi2
            do j = lo1, hi1
               do i = lo0, hi0
                   energy(i,j,k) = 0.d0
               enddo
            enddo
         enddo
      endif
c
      dx2inv = 1.d0 / dx(0)**2
      dy2inv = 1.d0 / dx(1)**2
      dz2inv = 1.d0 / dx(2)**2
c
      epsilonphi2 = 0.5d0 * epsilon_phi * epsilon_phi
c
      do p = 1, nphases
         do k = lo2, hi2
            do j = lo1, hi1
               do i = lo0, hi0
                  diff_term_x = dx2inv * (
     &               - phi(i+1,j,k,p)
     &               + 2.d0 * phi(i,j,k,p)
     &               - phi(i-1,j,k,p) )

                  diff_term_y = dy2inv * (
     &               - phi(i,j+1,k,p)
     &               + 2.d0 * phi(i,j,k,p)
     &               - phi(i,j-1,k,p) )

                  diff_term_z = dz2inv * (
     &               - phi(i,j,k+1,p)
     &               + 2.d0 * phi(i,j,k,p)
     &               - phi(i,j,k-1,p) )

                  diff_term = diff_term_x + diff_term_y + diff_term_z

                  e = epsilonphi2 * diff_term * phi(i,j,k,p)
                  if ( eval_per_cell /= 0 ) then
                     energy(i,j,k) = energy(i,j,k) + e
                  endif

                  e = e * weight(i,j,k)
                  phi_e = phi_e + e
               enddo
            enddo
         enddo
      enddo

      return
      end
c
      subroutine quatenergy(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth,
     &   dx,
     &   gqx, gqy, gqz, gqghosts,
     &   phi, pghosts, nphases,
     &   quat, ngq,
     &   epsilon_phi, epsilon_q,
     &   anisotropy, knumber,
     &   misorientation_factor,
     &   temperature, tghosts,
     &   phi_well_scale,
     &   weight,
     &   total_energy,
     &   total_phi_e,
     &   total_orient_e,
     &   total_qint_e,
     &   energy,
     &   eval_per_cell,
     &   phi_interp_type,
     &   phi_well_type,
     &   orient_interp_type1,
     &   orient_interp_type2,
     &   avg_type,
     &   floor_type,
     $   gradient_floor
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth, gqghosts, pghosts, tghosts, ngq,
     &   knumber, nphases

      double precision phi(CELL3d(lo,hi,pghosts),nphases)
      double precision temperature(CELL3d(lo,hi,tghosts))
      
      double precision quat(CELL3d(lo,hi,ngq),depth)
      double precision gqx(SIDE3d0(lo,hi,gqghosts),depth,NDIM)
      double precision gqy(SIDE3d1(lo,hi,gqghosts),depth,NDIM)
      double precision gqz(SIDE3d2(lo,hi,gqghosts),depth,NDIM)
      double precision energy(CELL3d(lo,hi,0))
      
      double precision weight(CELL3d(lo,hi,0))
      character*(*) phi_interp_type
      character*(*) phi_well_type
      character*(*) orient_interp_type1, orient_interp_type2
      character*(*) avg_type
      character*(*) floor_type
      double precision gradient_floor, floor2
      integer eval_per_cell
      integer three_phase

      double precision anisotropy, misorientation_factor
      double precision phi_well_scale
      double precision epsilon_phi, epsilon_q
      double precision epsilonq2
      double precision total_energy, e, o2
      double precision total_phi_e, phi_e, total_orient_e
      double precision total_qint_e, total_free_e
      double precision p_phi
      double precision dx(NDIM), dx2inv, dy2inv, dz2inv
      double precision diff_term_x, diff_term_y, diff_term_z
      double precision diff_term

      double precision g_phi, aphi, h_phi
      double precision interp_func
      double precision well_func
      double precision average_func

      integer i, j, k, m, n, p
c
      dx2inv = 1.d0 / dx(1)**2
      dy2inv = 1.d0 / dx(2)**2
      dz2inv = 1.d0 / dx(3)**2
c
c phi interface energy
c
      if ( anisotropy .gt. 0.d0 )then
         call interface_anisotropic_energy(
     &      lo0, hi0, lo1, hi1, lo2, hi2,
     &      dx, phi, pghosts, nphases, quat, ngq, depth,
     &      epsilon_phi, anisotropy, knumber, weight,
     &      phi_e, energy,
     &      eval_per_cell)
      else
         call phi_interface_energy(lo0, hi0, lo1, hi1, lo2, hi2,
     &      dx, phi, pghosts, nphases,
     &      epsilon_phi, weight, phi_e, energy, eval_per_cell)
      endif

      total_phi_e = total_phi_e + phi_e
      total_energy = total_energy + phi_e

c
c Orientational energy: average over 2 cell faces in each direction
c
      if ( misorientation_factor .gt. 0.d0 ) then
         if ( floor_type(1:1) .eq. 's' )then
            floor2 = gradient_floor*gradient_floor
         else
            floor2 = 0.d0
         endif
         total_orient_e = 0.d0
         do k = lo2, hi2
            do j = lo1, hi1
               do i = lo0, hi0

                  e = 0.d0

                  aphi = average_func(phi(i-1,j,k,1),phi(i,j,k,1),
     &                                   avg_type)
                  p_phi = interp_func( aphi, orient_interp_type1 )
                  o2 = 0.d0
                  do m = 1, depth
                     do n = 1, NDIM
                        o2 = o2 + gqx(i,j,k,m,n) * gqx(i,j,k,m,n)
                     enddo
                  enddo
                  o2 = o2 + floor2
                  e = e + sqrt(o2) * p_phi
c
                  aphi = average_func(phi(i,j,k,1),phi(i+1,j,k,1),
     &                                   avg_type)
                  p_phi = interp_func( aphi, orient_interp_type1 )
                  o2 = 0.d0
                  do m = 1, depth
                     do n = 1, NDIM
                        o2 = o2 + gqx(i+1,j,k,m,n) * gqx(i+1,j,k,m,n)
                     enddo
                  enddo
                  o2 = o2 + floor2
                  e = e + sqrt(o2) * p_phi
c
                  aphi = average_func(phi(i,j-1,k,1),phi(i,j,k,1),
     &                                   avg_type)
                  p_phi = interp_func( aphi, orient_interp_type1 )
                  o2 = 0.d0
                  do m = 1, depth
                     do n = 1, NDIM
                        o2 = o2 + gqy(i,j,k,m,n) * gqy(i,j,k,m,n)
                     enddo
                  enddo
                  o2 = o2 + floor2
                  e = e + sqrt(o2) * p_phi
c
                  aphi = average_func(phi(i,j,k,1),phi(i,j+1,k,1),
     &                                   avg_type)
                  p_phi = interp_func( aphi, orient_interp_type1 )
                  o2 = 0.d0
                  do m = 1, depth
                     do n = 1, NDIM
                        o2 = o2 + gqy(i,j+1,k,m,n) * gqy(i,j+1,k,m,n)
                     enddo
                  enddo
                  e = e + sqrt(o2) * p_phi
                  o2 = o2 + floor2
c
                  aphi = average_func(phi(i,j,k,1),phi(i,j,k-1,1),
     &                                   avg_type)
                  p_phi = interp_func( aphi, orient_interp_type1 )
                  o2 = 0.d0
                  do m = 1, depth
                     do n = 1, NDIM
                        o2 = o2 + gqz(i,j,k,m,n) * gqz(i,j,k,m,n)
                     enddo
                  enddo
                  o2 = o2 + floor2
                  e = e + sqrt(o2) * p_phi
c
                  aphi = average_func(phi(i,j,k,1),phi(i,j,k+1,1),
     &                                   avg_type)
                  p_phi = interp_func( aphi, orient_interp_type1 )
                  o2 = 0.d0
                  do m = 1, depth
                     do n = 1, NDIM
                        o2 = o2 + gqz(i,j,k+1,m,n) * gqz(i,j,k+1,m,n)
                     enddo
                  enddo
                  o2 = o2 + floor2
                  e = e + sqrt(o2) * p_phi

                  e = e * temperature(i,j,k)
c        factor 1/6 because of cell averaging
                  e = e * misorientation_factor / 6.d0
                  if ( eval_per_cell /= 0 ) then
                     energy(i,j,k) = energy(i,j,k) + e
                  endif
                  e = e * weight(i,j,k)

                  total_energy = total_energy + e
                  total_orient_e = total_orient_e + e
               enddo
            enddo
         enddo
c
c q interface energy
c
         total_qint_e = 0.d0
         epsilonq2 = 0.5d0 * epsilon_q * epsilon_q
         do k = lo2, hi2
            do j = lo1, hi1
               do i = lo0, hi0
                  p_phi = interp_func( phi(i,j,k,1),
     &                                 orient_interp_type2 )
                  e = 0.d0
                  do m = 1, depth
                     do n = 1, NDIM
                        e = e + (
     &                     gqx(i,  j,k  ,m,n) * gqx(i,  j,k  ,m,n)
     &                   + gqx(i+1,j,k  ,m,n) * gqx(i+1,j,k  ,m,n)
     &                   + gqy(i,j  ,k  ,m,n) * gqy(i,j  ,k  ,m,n)
     &                   + gqy(i,j+1,k  ,m,n) * gqy(i,j+1,k  ,m,n)
     &                   + gqz(i,j  ,k  ,m,n) * gqz(i,j  ,k  ,m,n)
     &                   + gqz(i,j  ,k+1,m,n) * gqz(i,j  ,k+1,m,n)
     &                   )
                     enddo
                  enddo
c                 factor 1/6 because of cell average
                  e = e* p_phi * epsilonq2 / 6.d0
                  if ( eval_per_cell /= 0 ) then
                     energy(i,j,k) = energy(i,j,k) + e
                  endif
                  e = e* weight(i,j,k)

                  total_energy = total_energy + e
                  total_qint_e = total_qint_e + e
               enddo
            enddo
         enddo
      endif
c
c double well energy
c
      do p = 1, nphases
         do k = lo2, hi2
            do j = lo1, hi1
               do i = lo0, hi0

                  g_phi = well_func( phi(i,j,k,p), phi_well_type )
                  e = phi_well_scale * g_phi

                  if ( eval_per_cell /= 0 ) then
                     energy(i,j,k) = energy(i,j,k) + e
                  endif
                  e = e * weight(i,j,k)

                  total_energy = total_energy + e
                  total_phi_e = total_phi_e + e
               enddo
            enddo
         enddo
      enddo

      return
      end

c
c bulk energy
c
      subroutine bulkenergy(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   phi, pghosts,
     &   fl, fa,
     &   weight,
     &   total_energy,
     &   bulk_energy,
     &   energy, eval_per_cell,
     &   phi_interp_type
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2, pghosts
      integer eval_per_cell
      character*(*) phi_interp_type
      double precision total_energy, bulk_energy

      double precision phi(CELL3d(lo,hi,pghosts))

      double precision fl(CELL3d(lo,hi,0))
      double precision fa(CELL3d(lo,hi,0))
      double precision energy(CELL3d(lo,hi,0))
      double precision weight(CELL3d(lo,hi,0))

      double precision f_l, f_a, e, h_phi
      integer i, j, k
      double precision interp_func

      bulk_energy = 0.d0

       do k = lo2, hi2
          do j = lo1, hi1
             do i = lo0, hi0

                f_l = fl(i,j,k)
                f_a = fa(i,j,k)

                h_phi = interp_func( phi(i,j,k),
     &                               phi_interp_type )

                e = ( 1.0d0 - h_phi ) * f_l + h_phi * f_a

                if ( eval_per_cell /= 0 ) then
                   energy(i,j,k) = energy(i,j,k) + e
                endif
                e = e * weight(i,j,k)

                total_energy = total_energy + e
                bulk_energy  = bulk_energy + e
            enddo
         enddo
      enddo

      return
      end

c
c
c
      subroutine temperature_energy(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   temperature, tghosts,
     &   fs,
     &   T_M, L_A
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2, i, j, k,
     &   tghosts

      double precision temperature(CELL3d(lo,hi,tghosts))
      double precision fs(CELL3d(lo,hi,0))

      double precision L_A, T_M

      double precision factor

      factor = L_A/T_M

      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0

                fs(i,j,k) = factor*(T_M-temperature(i,j,k))

            enddo
         enddo
      enddo

      return
      end
c
c see Moelans 2011
c
      subroutine phiphi_interfacialenergy(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   phi, pghosts, nphases,
     &   gamma, m,
     &   weight,
     &   edensity,
     &   total_energy,
     &   eval_per_cell
     &   )

      implicit none

      integer lo0, hi0, lo1, hi1, lo2, hi2
      integer pghosts, nphases
      integer eval_per_cell
      double precision gamma, m, total_energy(nphases*nphases)
      double precision phi(CELL3d(lo,hi,pghosts), nphases)
      double precision edensity(CELL3d(lo,hi,0))
      double precision weight(CELL3d(lo,hi,0))

      integer i, j, k, p1, p2
      double precision e, factor
c
      factor = 0.5d0 * gamma * m
      if ( eval_per_cell /= 0 ) then
         do k = lo2, hi2
            do j = lo1, hi1
               do i = lo0, hi0
                  edensity(i,j,k) = 0.d0
               enddo
            enddo
         enddo
      endif

c
      do p1 = 1, nphases
         do p2 = 1, nphases
            total_energy(nphases*(p1-1)+p2) = 0.d0
            do k = lo2, hi2
               do j = lo1, hi1
                  do i = lo0, hi0
                     e = phi(i,j,k,p1)*phi(i,j,k,p1)
     &                 * phi(i,j,k,p2)*phi(i,j,k,p2)
                     e = e * factor
                     if ( eval_per_cell /= 0 ) then
                        edensity(i,j,k) = e
                     endif
                     e = e * weight(i,j,k)

                     total_energy(nphases*(p1-1)+p2) =
     &                  total_energy(nphases*(p1-1)+p2) + e
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end
