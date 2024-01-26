c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
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
