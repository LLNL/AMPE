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
c
      subroutine quatenergy(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth,
     &   dx,
     &   gqx, gqy, gqz, gqghosts,
     &   phi, pghosts,
     &   eta, ngeta,
     &   epsilon_phi, epsilon_eta, epsilon_q,
     &   misorientation_factor,
     &   temperature, tghosts,
     &   phi_well_scale,
     &   eta_well_scale,
     &   fl, fa, fb,
     &   three_phase,
     &   weight,
     &   total_energy,
     &   total_phi_e,
     &   total_eta_e,
     &   total_orient_e,
     &   total_qint_e,
     &   total_well_e,
     &   total_free_e,
     &   energy,
     &   eval_per_cell,
     &   phi_interp_type,
     &   eta_interp_type,
     &   phi_well_type,
     &   eta_well_type,
     &   orient_interp_type,
     &   avg_type,
     &   floor_type,
     $   gradient_floor
     &   )

      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth, gqghosts, pghosts, tghosts, ngeta

      double precision phi(CELL3d(lo,hi,pghosts))
      double precision eta(CELL3d(lo,hi,ngeta))
      double precision temperature(CELL3d(lo,hi,tghosts))
      
      double precision gqx(SIDE3d0(lo,hi,gqghosts),depth,NDIM)
      double precision gqy(SIDE3d1(lo,hi,gqghosts),depth,NDIM)
      double precision gqz(SIDE3d2(lo,hi,gqghosts),depth,NDIM)
      double precision fl(CELL3d(lo,hi,0))
      double precision fa(CELL3d(lo,hi,0))
      double precision fb(CELL3d(lo,hi,0))
      double precision energy(CELL3d(lo,hi,0))
      
      double precision weight(CELL3d(lo,hi,0))
      character*(*) phi_interp_type
      character*(*) eta_interp_type
      character*(*) phi_well_type
      character*(*) eta_well_type
      character*(*) orient_interp_type
      character*(*) avg_type
      character*(*) floor_type
      double precision gradient_floor, floor2
      integer eval_per_cell
      integer three_phase

      double precision misorientation_factor
      double precision phi_well_scale, eta_well_scale
      double precision f_l, f_a, f_b
      double precision epsilon_phi, epsilon_q, epsilon_eta
      double precision epsilonphi2, epsilonq2, epsiloneta2
      double precision total_energy, e, o2
      double precision total_phi_e, total_orient_e
      double precision total_qint_e, total_well_e, total_free_e
      double precision total_eta_e
      double precision p_phi
      double precision dx(0:NDIM-1), dx2inv, dy2inv, dz2inv
      double precision diff_term_x, diff_term_y, diff_term_z
      double precision diff_term

      double precision g_phi, aphi, g_eta, h_phi, h_eta
      double precision interp_func
      double precision well_func
      double precision average_func

      integer i, j, k, m, n
c
      dx2inv = 1.d0 / dx(0)**2
      dy2inv = 1.d0 / dx(1)**2
      dz2inv = 1.d0 / dx(2)**2
c
      if ( floor_type(1:1) .eq. 's' )then
         floor2 = gradient_floor*gradient_floor
      else
         floor2 = 0.d0
      endif
      epsilonphi2 = 0.5d0 * epsilon_phi * epsilon_phi
      if ( three_phase /= 0 ) then
         epsiloneta2 = 0.5d0 * epsilon_eta * epsilon_eta
      endif

c
c
c phi interface energy
c      
      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0
               diff_term_x = dx2inv * (
     &            - phi(i+1,j,k)
     &            + 2.d0 * phi(i,j,k)
     &            - phi(i-1,j,k) )

               diff_term_y = dy2inv * (
     &            - phi(i,j+1,k) 
     &            + 2.d0 * phi(i,j,k) 
     &            - phi(i,j-1,k) )

               diff_term_z = dz2inv * (
     &            - phi(i,j,k+1) 
     &            + 2.d0 * phi(i,j,k) 
     &            - phi(i,j,k-1) )

               diff_term = diff_term_x + diff_term_y + diff_term_z

               e = epsilonphi2 * diff_term * phi(i,j,k)
               if ( eval_per_cell /= 0 ) then
                  energy(i,j,k) = e
               endif
               e = e * weight(i,j,k)

               total_energy = total_energy + e
               total_phi_e = total_phi_e + e
            enddo
         enddo
      enddo

c
c eta interface energy
c    
      if ( three_phase /= 0 ) then
         do k = lo2, hi2
            do j = lo1, hi1
               do i = lo0, hi0
                  diff_term_x = dx2inv * (
     &               - eta(i+1,j,k)
     &               + 2.d0 * eta(i,j,k) 
     &               - eta(i-1,j,k) )

                  diff_term_y = dy2inv * (
     &               - eta(i,j+1,k) 
     &               + 2.d0 * eta(i,j,k) 
     &               - eta(i,j-1,k) )

                  diff_term_z = dz2inv * (
     &               - eta(i,j,k+1) 
     &               + 2.d0 * eta(i,j,k) 
     &               - eta(i,j,k-1) )

                  diff_term = diff_term_x + diff_term_y + diff_term_z

                  e = epsiloneta2 * diff_term * eta(i,j,k)
                  if ( eval_per_cell /= 0 ) then
                     energy(i,j,k) = energy(i,j,k) + e
                  endif
                  e = e * weight(i,j,k)

                  total_energy = total_energy + e
                  total_eta_e = total_eta_e + e
               enddo
            enddo
         enddo
      endif

c
c Orientational energy: average over 2 cell faces in each direction
c
      if ( misorientation_factor .gt. 0.d0 ) then
         do k = lo2, hi2
            do j = lo1, hi1
               do i = lo0, hi0

                  e = 0.d0

                  aphi = average_func(phi(i-1,j,k),phi(i,j,k),avg_type)
                  p_phi = interp_func( aphi, orient_interp_type )
                  o2 = 0.d0
                  do m = 1, depth
                     do n = 1, NDIM
                        o2 = o2 + gqx(i,j,k,m,n) * gqx(i,j,k,m,n)
                     enddo
                  enddo
                  o2 = o2 + floor2
                  e = e + sqrt(o2) * p_phi
c
                  aphi = average_func(phi(i,j,k),phi(i+1,j,k),avg_type)
                  p_phi = interp_func( aphi, orient_interp_type )
                  o2 = 0.d0
                  do m = 1, depth
                     do n = 1, NDIM
                        o2 = o2 + gqx(i+1,j,k,m,n) * gqx(i+1,j,k,m,n)
                     enddo
                  enddo
                  o2 = o2 + floor2
                  e = e + sqrt(o2) * p_phi
c
                  aphi = average_func(phi(i,j-1,k),phi(i,j,k),avg_type)
                  p_phi = interp_func( aphi, orient_interp_type )
                  o2 = 0.d0
                  do m = 1, depth
                     do n = 1, NDIM
                        o2 = o2 + gqy(i,j,k,m,n) * gqy(i,j,k,m,n)
                     enddo
                  enddo
                  o2 = o2 + floor2
                  e = e + sqrt(o2) * p_phi
c
                  aphi = average_func(phi(i,j,k),phi(i,j+1,k),avg_type)
                  p_phi = interp_func( aphi, orient_interp_type )
                  o2 = 0.d0
                  do m = 1, depth
                     do n = 1, NDIM
                        o2 = o2 + gqy(i,j+1,k,m,n) * gqy(i,j+1,k,m,n)
                     enddo
                  enddo
                  e = e + sqrt(o2) * p_phi
                  o2 = o2 + floor2
c
                  aphi = average_func(phi(i,j,k),phi(i,j,k-1),avg_type)
                  p_phi = interp_func( aphi, orient_interp_type )
                  o2 = 0.d0
                  do m = 1, depth
                     do n = 1, NDIM
                        o2 = o2 + gqz(i,j,k,m,n) * gqz(i,j,k,m,n)
                     enddo
                  enddo
                  o2 = o2 + floor2
                  e = e + sqrt(o2) * p_phi
c
                  aphi = average_func(phi(i,j,k),phi(i,j,k+1),avg_type)
                  p_phi = interp_func( aphi, orient_interp_type )
                  o2 = 0.d0
                  do m = 1, depth
                     do n = 1, NDIM
                        o2 = o2 + gqz(i,j,k+1,n,m) * gqz(i,j,k+1,n,m)
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
         epsilonq2 = 0.5d0 * epsilon_q * epsilon_q
         do k = lo2, hi2
            do j = lo1, hi1
               do i = lo0, hi0
                  e = 0.d0
                  do m = 1, depth
                     do n = 1, NDIM
                        e = e + (
     &                     gqx(i,  j,k  ,n,m) * gqx(i,  j,k  ,n,m)
     &                   + gqx(i+1,j,k  ,n,m) * gqx(i+1,j,k  ,n,m)
     &                   + gqy(i,j  ,k  ,n,m) * gqy(i,j  ,k  ,n,m)
     &                   + gqy(i,j+1,k  ,n,m) * gqy(i,j+1,k  ,n,m)
     &                   + gqz(i,j  ,k  ,n,m) * gqz(i,j  ,k  ,n,m)
     &                   + gqz(i,j  ,k+1,n,m) * gqz(i,j  ,k+1,n,m)
     &                   )
                     enddo
                  enddo
c                 factor 1/6 because of cell average
                  e = e* epsilonq2/6.d0
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
      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0

               g_phi = well_func( phi(i,j,k), phi_well_type )

               if ( three_phase /= 0 ) then
                  g_eta = well_func( eta, eta_well_type )
               else
                  g_eta = 0.0d0
               endif
            
               e =(
     &            phi_well_scale * g_phi
     &            )

               if ( three_phase /= 0 ) then
                  h_phi = interp_func( phi(i,j,k), phi_interp_type )
                  e = e + eta_well_scale * h_phi * g_eta
               endif
               if ( eval_per_cell /= 0 ) then
                  energy(i,j,k) = energy(i,j,k) + e
               endif
               e = e * weight(i,j,k)

               total_energy = total_energy + e
               total_well_e = total_well_e + e
            enddo
         enddo
      enddo

c
c free energy
c
      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0

               f_l = fl(i,j,k)
               f_a = fa(i,j,k)

               h_phi = interp_func( phi(i,j,k), phi_interp_type )

               if ( three_phase /= 0 ) then
                  f_b = fb(i,j,k)
                  h_eta = interp_func( eta, eta_interp_type )
               else
                  f_b = 0.0d0
                  h_eta = 0.0d0
               endif
            
               e = (
     &            ( 1.0d0 - h_phi ) * f_l +
     &            h_phi * ( ( 1.0d0 - h_eta ) * f_a + h_eta * f_b )
     &            )
               if ( eval_per_cell /= 0 ) then
                  energy(i,j,k) = energy(i,j,k) + e
               endif
               e = e * weight(i,j,k)

               total_energy = total_energy + e
               total_free_e = total_free_e + e
            enddo
         enddo
      enddo

      return
      end
