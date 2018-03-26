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
      
      integer i, j, m, iq
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
