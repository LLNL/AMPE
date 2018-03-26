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
      subroutine tagfromgrads(
     &   lo0, hi0, lo1, hi1,
     &   gradx, grady,
     &   glo0, ghi0, glo1, ghi1, 
     &   h,
     &   tags,
     &   tlo0, thi0, tlo1, thi1, 
     &   refine_tag_val,
     &   tolerance )

      implicit none

c        input arrays:
      integer lo0, hi0, lo1, hi1
      integer glo0, ghi0, glo1, ghi1
      integer tlo0, thi0, tlo1, thi1
      double precision tolerance
      integer refine_tag_val

c        variables in 2d cell indexed         
      integer tags(tlo0:thi0,tlo1:thi1)
      double precision
     &   gradx(glo0:ghi0,glo1:ghi1),
     &   grady(glo0:ghi0,glo1:ghi1),
     &   h(2)

      integer ii, jj
      double precision gx, gy, g, tol_sq

      tol_sq = tolerance * tolerance

      do jj = lo1, hi1
         do ii = lo0, hi0

            gx = gradx(ii,jj) * h(1)
            gy = grady(ii,jj) * h(2)

            g = gx*gx + gy*gy

            if ( g .ge. tol_sq )
     &         tags(ii,jj) = refine_tag_val

         enddo
      enddo

      return
      end

c=======================================================================

      subroutine tagfromquatgrads(
     &   lo0, hi0, lo1, hi1,
     &   qgradx, qgrady,
     &   glo0, ghi0, glo1, ghi1, 
     &   depth,
     &   h,
     &   tags,
     &   tlo0, thi0, tlo1, thi1, 
     &   refine_tag_val,
     &   tolerance )

      implicit none

c        input arrays:
      integer lo0, hi0, lo1, hi1
      integer glo0, ghi0, glo1, ghi1
      integer tlo0, thi0, tlo1, thi1
      integer depth
      double precision tolerance
      integer refine_tag_val

c        variables in 2d cell indexed         
      integer tags(tlo0:thi0,tlo1:thi1)
      double precision
     &   qgradx(glo0:ghi0,glo1:ghi1,depth),
     &   qgrady(glo0:ghi0,glo1:ghi1,depth),
     &   h(2)

      integer ii, jj, m
      double precision gx, gy, g, tol_sq

      tol_sq = tolerance * tolerance

      do jj = lo1, hi1
         do ii = lo0, hi0

            g = 0.0d0

            do m = 1, depth
               gx = qgradx(ii,jj,m) * h(1)
               g = g + gx*gx

               gy = qgrady(ii,jj,m) * h(2)
               g = g + gy*gy
            enddo

            if ( g .ge. tol_sq )
     &         tags(ii,jj) = refine_tag_val

         enddo
      enddo

      return
      end
