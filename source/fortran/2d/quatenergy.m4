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
      subroutine interface_anisotropic_energy(
     &   lo0, hi0, lo1, hi1,
     &   dx,
     &   phi, pghosts, nphases,
     &   quat, ngq, qlen,
     &   epsilon,
     &   anisotropy, knumber,
     &   weight,
     &   phi_e,
     &   energy,
     &   eval_per_cell
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1,
     &    pghosts, ngq, qlen, knumber, nphases
      integer eval_per_cell

      double precision phi_e, epsilon
      double precision phi(CELL2d(lo,hi,pghosts),nphases)
      double precision quat(CELL2d(lo,hi,ngq),qlen)
      double precision energy(CELL2d(lo,hi,0))
      double precision weight(CELL2d(lo,hi,0))

      double precision dx(NDIM)
      double precision anisotropy

      double precision dxinv, dyinv, theta, q
      double precision epstheta, diff_term, dphidx, dphidy
      double precision pi, e, refangle

      integer i, j, p
c
      pi = 4.d0*atan(1.d0)
c
      dxinv = 0.5d0 / dx(1)
      dyinv = 0.5d0 / dx(2)
c
      phi_e = 0.d0
      if ( eval_per_cell /= 0 ) then
         do j = lo1, hi1
            do i = lo0, hi0
               energy(i,j) = 0.d0
            enddo
         enddo
      endif
c
      do p = 1, nphases
         do j = lo1, hi1
            do i = lo0, hi0
               dphidx = (phi(i+1,j,p) - phi(i-1,j,p)) * dxinv
               dphidy = (phi(i,j+1,p) - phi(i,j-1,p)) * dyinv
               if( abs(dphidx)>1.e-12 )then
                  theta=atan(dphidy/dphidx)
               else
                  theta=0.5d0*pi
               endif

               q=quat(i,j,1)
c q could be slightly out of the [-1,+1] range due to finite precision
c need to enforce q to be within that range before call to acos
               if( q .gt.  1.d0 )q=1.d0
               if( q .lt. -1.d0 )q=-1.d0

               if( qlen==4 )then
                  refangle=2.d0*acos(q)
               else
                  refangle=acos(q)
               endif

               epstheta=epsilon*(1.d0
     &                  +anisotropy*cos(knumber*(theta-refangle)))

               diff_term = dphidx*dphidx + dphidy*dphidy

               e = 0.5d0 * epstheta * epstheta * diff_term
               if ( eval_per_cell /= 0 ) then
                  energy(i,j) = energy(i,j) + e
               endif

               e = e * weight(i,j)
               phi_e = phi_e + e
            enddo
         enddo
      enddo

      return
      end
c
c
c
      subroutine phi_interface_energy(
     &   lo0, hi0, lo1, hi1,
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
     &   lo0, hi0, lo1, hi1,
     &    pghosts, nphases
      integer eval_per_cell

      double precision phi_e, epsilon_phi, e
      double precision phi(CELL2d(lo,hi,pghosts), nphases)
      double precision energy(CELL2d(lo,hi,0))
      double precision weight(CELL2d(lo,hi,0))

      double precision dx(0:NDIM-1), dx2inv, dy2inv
      double precision diff_term_x, diff_term_y
      double precision diff_term
      double precision epsilonphi2

      integer i, j, p
c
      dx2inv = 1.d0 / dx(0)**2
      dy2inv = 1.d0 / dx(1)**2
c
      epsilonphi2 = 0.5d0 * epsilon_phi * epsilon_phi
c
      phi_e = 0.d0
      if ( eval_per_cell /= 0 ) then
         do j = lo1, hi1
            do i = lo0, hi0
               energy(i,j) = 0.d0
            enddo
         enddo
      endif
c
      do p = 1, nphases
         do j = lo1, hi1
            do i = lo0, hi0
               diff_term_x = dx2inv * (
     &            - phi(i+1,j,p)
     &            + 2.d0 * phi(i,j,p)
     &            - phi(i-1,j,p) )

               diff_term_y = dy2inv * (
     &            - phi(i,j+1,p)
     &            + 2.d0 * phi(i,j,p)
     &            - phi(i,j-1,p) )

               diff_term = diff_term_x + diff_term_y

               e = epsilonphi2 * diff_term * phi(i,j,p)
               if ( eval_per_cell /= 0 ) then
                  energy(i,j) = energy(i,j) + e
               endif

               e = e * weight(i,j)
               phi_e = phi_e + e
            enddo
         enddo
      enddo

      return
      end
c
      subroutine quatenergy(
     &   lo0, hi0, lo1, hi1,
     &   depth,
     &   dx,
     &   gqx,gqy,gqghosts,
     &   phi, pghosts, nphases,
     &   quat, ngq,
     &   epsilon_phi, epsilon_q,
     &   anisotropy, knumber,
     &   misorientation_factor,
     &   temperature, tghosts,
     &   phi_well_scale,
     &   weight,
     &   total_energy,
     &   total_interface_e,
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
     &   lo0, hi0, lo1, hi1, nphases,
     &   depth, gqghosts, pghosts, tghosts, ngq,
     &   knumber

      double precision phi(CELL2d(lo,hi,pghosts),nphases)
      double precision temperature(CELL2d(lo,hi,tghosts))

      double precision quat(CELL2d(lo,hi,ngq),depth)
      double precision gqx(SIDE2d0(lo,hi,gqghosts),depth,NDIM)
      double precision gqy(SIDE2d1(lo,hi,gqghosts),depth,NDIM)
      double precision energy(CELL2d(lo,hi,0))
      
      double precision weight(CELL2d(lo,hi,0))
      character*(*) phi_interp_type
      character*(*) phi_well_type
      character*(*) orient_interp_type1, orient_interp_type2
      character*(*) avg_type
      character*(*) floor_type
      double precision gradient_floor, floor2
      integer eval_per_cell

      double precision anisotropy, misorientation_factor
      double precision phi_well_scale
      double precision epsilon_phi, epsilon_q
      double precision epsilonq2
      double precision total_energy, e, o2
      double precision total_interface_e, phi_e, total_orient_e
      double precision total_qint_e
      double precision p_phi
      double precision dx(NDIM), dx2inv, dy2inv

      double precision g_phi, aphi
      double precision interp_func
      double precision well_func
      double precision average_func

      integer i, j, m, n, p
c
      dx2inv = 1.d0 / dx(1)**2
      dy2inv = 1.d0 / dx(2)**2
c
c phi interface energy
c
      if ( anisotropy .gt. 0.d0 )then
         call interface_anisotropic_energy(lo0, hi0, lo1, hi1,
     &      dx, phi, pghosts, nphases, quat, ngq, depth,
     &      epsilon_phi, anisotropy, knumber, weight,
     &      phi_e, energy,
     &      eval_per_cell)
      else
         call phi_interface_energy(lo0, hi0, lo1, hi1, dx, phi, pghosts,
     &      nphases,
     &      epsilon_phi, weight, phi_e, energy, eval_per_cell)
      endif

      total_interface_e = total_interface_e + phi_e
      total_energy = total_energy + phi_e

c
c Orientational energy: average over 2 cell sides in each direction
c
      if ( misorientation_factor .gt. 0.d0 ) then
         if ( floor_type(1:1) .eq. 's' )then
            floor2 = gradient_floor*gradient_floor
         else
            floor2 = 0.d0
         endif
         do p = 1, nphases
            do j = lo1, hi1
               do i = lo0, hi0

                  e = 0.d0

                  aphi = average_func( phi(i-1,j,p), phi(i,j,p),
     &                                 avg_type)
                  p_phi = interp_func( aphi, orient_interp_type1 )
                  o2 = 0.d0               
                  do n = 1, NDIM
                     do m = 1, depth
                        o2 = o2 + gqx(i,j,m,n) * gqx(i,j,m,n)
                     enddo
                  enddo
                  o2 = o2 + floor2
                  e = e + sqrt(o2) * p_phi
c
                  aphi = average_func( phi(i,j,p), phi(i+1,j,p),
     &                                 avg_type)
                  p_phi = interp_func( aphi, orient_interp_type1 )
                  o2 = 0.d0
                  do n = 1, NDIM
                     do m = 1, depth
                        o2 = o2 + gqx(i+1,j,m,n) * gqx(i+1,j,m,n)
                     enddo
                  enddo
                  o2 = o2 + floor2
                  e = e + sqrt(o2) * p_phi
c
                  aphi = average_func( phi(i,j-1,p), phi(i,j,p),
     &                                 avg_type)
                  p_phi = interp_func( aphi, orient_interp_type1 )
                  o2 = 0.d0
                  do n = 1, NDIM
                     do m = 1, depth
                        o2 = o2 + gqy(i,j,m,n) * gqy(i,j,m,n)
                     enddo
                  enddo
                  o2 = o2 + floor2
                  e = e + sqrt(o2) * p_phi
c
                  aphi = average_func( phi(i,j,p), phi(i,j+1,p),
     &                                 avg_type)
                  p_phi = interp_func( aphi, orient_interp_type1 )
                  o2 = 0.d0
                  do n = 1, NDIM
                     do m = 1, depth
                        o2 = o2 + gqy(i,j+1,m,n) * gqy(i,j+1,m,n)
                     enddo
                  enddo
                  o2 = o2 + floor2
                  e = e + sqrt(o2) * p_phi
               
                  e = e * temperature(i,j)
c        factor 0.25 because of cell averaging and double counting
                  e = e * 0.25d0 * misorientation_factor
                  if ( eval_per_cell /= 0 ) then
                     energy(i,j) = energy(i,j) + e
                  endif
                  e = e * weight(i,j)

                  total_energy   = total_energy   + e
                  total_orient_e = total_orient_e + e
               enddo
            enddo
c
c q interface energy
c
            epsilonq2 = 0.5d0 * epsilon_q * epsilon_q
            do j = lo1, hi1
               do i = lo0, hi0
                  p_phi = interp_func( phi(i,j,p), orient_interp_type2 )
                  e = 0.d0
                  do n = 1, NDIM
                     do m = 1, depth
                        e = e + (
     &                     gqx(i,  j,m,n) * gqx(i,  j,m,n)
     &                   + gqx(i+1,j,m,n) * gqx(i+1,j,m,n)
     &                   + gqy(i,j  ,m,n) * gqy(i,j  ,m,n)
     &                   + gqy(i,j+1,m,n) * gqy(i,j+1,m,n)
     &                   )
                     enddo
                  enddo
c                 factor 0.25 because of cell average and double counting
                  e = e * 0.25d0 * epsilonq2 * p_phi
                  if ( eval_per_cell /= 0 ) then
                     energy(i,j) = energy(i,j) + e
                  endif
                  e = e * weight(i,j)

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
         do j = lo1, hi1
            do i = lo0, hi0

               g_phi = well_func( phi(i,j,p), phi_well_type )

               e = phi_well_scale * g_phi

               if ( eval_per_cell /= 0 ) then
                  energy(i,j) = energy(i,j) + e
               endif
               e = e * weight(i,j)

               total_energy = total_energy + e
               total_interface_e = total_interface_e + e
            enddo
         enddo
      enddo

      return
      end

c
c bulk energy
c
      subroutine bulkenergy(
     &   lo0, hi0, lo1, hi1,
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
     &   lo0, hi0, lo1, hi1, pghosts
      integer eval_per_cell
      character*(*) phi_interp_type
      double precision total_energy, bulk_energy

      double precision phi(CELL2d(lo,hi,pghosts))

      double precision fl(CELL2d(lo,hi,0))
      double precision fa(CELL2d(lo,hi,0))
      double precision energy(CELL2d(lo,hi,0))
      double precision weight(CELL2d(lo,hi,0))

      double precision f_l, f_a, e, h_phi
      integer i, j
      double precision interp_func

      do j = lo1, hi1
         do i = lo0, hi0

            f_l = fl(i,j)
            f_a = fa(i,j)

            h_phi = interp_func( phi(i,j), phi_interp_type )

            e = ( 1.0d0 - h_phi ) * f_l + h_phi * f_a

            if ( eval_per_cell /= 0 ) then
               energy(i,j) = energy(i,j) + e
            endif
            e = e * weight(i,j)

            total_energy = total_energy + e
            bulk_energy = bulk_energy + e
         enddo
      enddo

      return
      end

c
c
c
      subroutine temperature_energy(
     &   lo0, hi0, lo1, hi1,
     &   temperature, tghosts,
     &   fl,
     &   T_M, L_A
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1, i, j,
     &   tghosts

      double precision temperature(CELL2d(lo,hi,tghosts))
      double precision fl(CELL2d(lo,hi,0))

      double precision L_A, T_M

      double precision factor

      factor = L_A/T_M

      do j = lo1, hi1
         do i = lo0, hi0

             fl(i,j) = factor*(T_M-temperature(i,j))

         enddo
      enddo

      return
      end
c
c see Moelans 2011
c
      subroutine phiphi_interfacialenergy(
     &   lo0, hi0, lo1, hi1,
     &   phi, pghosts, nphases,
     &   gamma, m,
     &   weight,
     &   edensity,
     &   total_energy,
     &   eval_per_cell
     &   )

      implicit none

      integer lo0, hi0, lo1, hi1
      integer pghosts, nphases
      integer eval_per_cell
      double precision gamma, m, total_energy(nphases*nphases)
      double precision phi(CELL2d(lo,hi,pghosts), nphases)
      double precision edensity(CELL2d(lo,hi,0))
      double precision weight(CELL2d(lo,hi,0))

      integer i, j, p1, p2
      double precision e, factor
c
      factor = 0.5d0 * gamma * m
      if ( eval_per_cell /= 0 ) then
         do j = lo1, hi1
            do i = lo0, hi0
               edensity(i,j) = 0.d0
            enddo
         enddo
      endif

c
      do p1 = 1, nphases
         do p2 = 1, nphases
            total_energy(nphases*(p1-1)+p2) = 0.d0
            do j = lo1, hi1
               do i = lo0, hi0
                  e = phi(i,j,p1)*phi(i,j,p1)
     &              * phi(i,j,p2)*phi(i,j,p2)
                  e = e * factor
                  if ( eval_per_cell /= 0 ) then
                     edensity(i,j) = e
                  endif
                  e = e * weight(i,j)

                  total_energy(nphases*(p1-1)+p2) =
     &               total_energy(nphases*(p1-1)+p2) + e
               enddo
            enddo
         enddo
      enddo

      return
      end
