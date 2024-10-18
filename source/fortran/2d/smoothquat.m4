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
      subroutine smoothquat(
     &   lo0, hi0, lo1, hi1,
     &   depth,
     &   quat,  ngq,
     &   quat_out, ngq_out,
     &   symm_diff_x, symm_diff_y,
     &   ngdiff,
     &   phase, ngphi,
     &   phase_threshold
     &   )

      implicit none

      integer
     &   lo0, lo1, hi0, hi1,
     &   depth, ngdiff, ngq, ngq_out, ngphi
     
      double precision phase_threshold

      double precision symm_diff_x(SIDE2d0(lo,hi,ngdiff),depth)
      double precision symm_diff_y(SIDE2d1(lo,hi,ngdiff),depth)

      double precision quat(CELL2d(lo,hi,ngq),depth)
      double precision quat_out(CELL2d(lo,hi,ngq_out),depth)

      double precision phase(CELL2d(lo,hi,ngphi))
      
      integer i, j, m
      do j = lo1, hi1
         do i = lo0, hi0         
               
            if( phase(i,j)<phase_threshold )then
               do m = 1, depth
                  quat_out(i,j,m) = quat(i,j,m) 
     &               - 0.125*symm_diff_x(i,  j,  m)
     &               + 0.125*symm_diff_x(i+1,j,  m)
     &               - 0.125*symm_diff_y(i,  j,  m)
     &               + 0.125*symm_diff_y(i,  j+1,m)
               enddo

            endif

         enddo
      enddo

      return
      end
