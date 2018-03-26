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
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   gradx, grady, gradz,
     &   dlo0, dhi0, dlo1, dhi1, dlo2, dhi2,
     &   h,
     &   tags,
     &   tlo0, thi0, tlo1, thi1, tlo2, thi2,
     &   refine_tag_val,
     &   tolerance )

      implicit none

c        input arrays:
      integer lo0, hi0, lo1, hi1, lo2, hi2
      integer dlo0, dhi0, dlo1, dhi1, dlo2, dhi2
      integer tlo0, thi0, tlo1, thi1, tlo2, thi2
      double precision tolerance
      integer refine_tag_val

c        variables in 3d cell indexed         
      integer tags(tlo0:thi0,tlo1:thi1,tlo2:thi2)
      double precision
     &   gradx(dlo0:dhi0,dlo1:dhi1,dlo2:dhi2),
     &   grady(dlo0:dhi0,dlo1:dhi1,dlo2:dhi2),
     &   gradz(dlo0:dhi0,dlo1:dhi1,dlo2:dhi2),
     &   h(3)

      integer ii, jj, kk
      double precision gx, gy, gz, g, tol_sq

      tol_sq = tolerance * tolerance

      do kk = lo2, hi2
         do jj = lo1, hi1
            do ii = lo0, hi0

               gx = gradx(ii,jj,kk) * h(1)
               gy = grady(ii,jj,kk) * h(2)
               gz = gradz(ii,jj,kk) * h(3)

               g = gx*gx + gy*gy + gz*gz

               if ( g .ge. tol_sq )
     &            tags(ii,jj,kk) = refine_tag_val

            enddo
         enddo
      enddo

      return
      end

c=======================================================================

      subroutine tagfromquatgrads(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   qgradx, qgrady, qgradz,
     &   dlo0, dhi0, dlo1, dhi1, dlo2, dhi2, 
     &   depth,
     &   h,
     &   tags,
     &   tlo0, thi0, tlo1, thi1, tlo2, thi2,
     &   refine_tag_val,
     &   tolerance )

      implicit none

c        input arrays:
      integer lo0, hi0, lo1, hi1, lo2, hi2
      integer dlo0, dhi0, dlo1, dhi1, dlo2, dhi2
      integer tlo0, thi0, tlo1, thi1, tlo2, thi2
      integer depth
      double precision tolerance
      integer refine_tag_val

c        variables in 3d cell indexed         
      integer tags(tlo0:thi0,tlo1:thi1,tlo2:thi2)
      double precision
     &   qgradx(dlo0:dhi0,dlo1:dhi1,dlo2:dhi2,depth),
     &   qgrady(dlo0:dhi0,dlo1:dhi1,dlo2:dhi2,depth),
     &   qgradz(dlo0:dhi0,dlo1:dhi1,dlo2:dhi2,depth),
     &   h(3)

      integer ii, jj, kk, m
      double precision gx, gy, gz, g, tol_sq

      tol_sq = tolerance * tolerance

      do kk = lo2, hi2
         do jj = lo1, hi1
            do ii = lo0, hi0

               g = 0.0d0

               do m = 1, depth
                  gx = qgradx(ii,jj,kk,m) * h(1)
                  g = g + gx*gx

                  gy = qgrady(ii,jj,kk,m) * h(2)
                  g = g + gy*gy

                  gz = qgradz(ii,jj,kk,m) * h(3)
                  g = g + gz*gz
               enddo

               if ( g .ge. tol_sq )
     &            tags(ii,jj,kk) = refine_tag_val

            enddo
         enddo
      enddo

      return
      end
