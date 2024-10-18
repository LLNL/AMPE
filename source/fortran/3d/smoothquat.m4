c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
define(NDIM,3)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c        
      subroutine smoothquat(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth,
     &   quat,  ngq,
     &   quat_out, ngq_out,
     &   symm_diff_x, symm_diff_y, symm_diff_z,
     &   ngdiff,
     &   phase, ngphi,
     &   phase_threshold
     &   )

      implicit none

      integer
     &   lo0, lo1, hi0, hi1, lo2, hi2,
     &   depth, ngdiff, ngq, ngq_out, ngphi
     
      double precision phase_threshold

      double precision symm_diff_x(SIDE3d0(lo,hi,ngdiff),depth)
      double precision symm_diff_y(SIDE3d1(lo,hi,ngdiff),depth)
      double precision symm_diff_z(SIDE3d2(lo,hi,ngdiff),depth)

      double precision quat(CELL3d(lo,hi,ngq),depth)
      double precision quat_out(CELL3d(lo,hi,ngq_out),depth)

      double precision phase(CELL3d(lo,hi,ngphi))
      
      integer i, j, k, m
      double precision fac
      
      fac=1.d0/6.d0
      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0         
               
               if( phase(i,j,k)<phase_threshold )then
                  do m = 1, depth
                     quat_out(i,j,k,m) = quat(i,j,k,  m) 
     &                  - fac*symm_diff_x(i,  j,  k,  m)
     &                  + fac*symm_diff_x(i+1,j,  k,  m)
     &                  - fac*symm_diff_y(i,  j,  k,  m)
     &                  + fac*symm_diff_y(i,  j+1,k,  m)
     &                  - fac*symm_diff_z(i,  j,  k,  m)
     &                  + fac*symm_diff_z(i,  j,  k+1,m)
                  enddo

               endif

            enddo
         enddo
      enddo

      return
      end
