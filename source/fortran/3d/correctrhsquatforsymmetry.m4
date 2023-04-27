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
c        Correct r.h.s. for quaternion implicit time evolution equation
c        Generalization of method proposed by Warren, Kobayashi et al., 
c        Act. Mater.51 (2003)
c        
      subroutine correctrhsquatforsymmetry(
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   depth, dx,
     &   nonsymm_diff_x, nonsymm_diff_y, nonsymm_diff_z,
     &   symm_diff_x, symm_diff_y, symm_diff_z,
     &   ngdiff,
     &   rhs, ngrhs,
     &   quat, ngq,
     &   facecoeffx, facecoeffy, facecoeffz, ngfacecoeff,
     &   mobility, ngmob,
     &   iqrot_x, iqrot_y, iqrot_z, ngiq
     &   )

      implicit none

      integer
     &   lo0, lo1, lo2, hi0, hi1, hi2,
     &   depth, ngdiff, ngrhs, ngq, ngfacecoeff, ngmob, ngiq

      double precision dx(NDIM)

      double precision nonsymm_diff_x(SIDE3d0(lo,hi,ngdiff),depth)
      double precision nonsymm_diff_y(SIDE3d1(lo,hi,ngdiff),depth)
      double precision nonsymm_diff_z(SIDE3d2(lo,hi,ngdiff),depth)

      double precision symm_diff_x(SIDE3d0(lo,hi,ngdiff),depth)
      double precision symm_diff_y(SIDE3d1(lo,hi,ngdiff),depth)
      double precision symm_diff_z(SIDE3d2(lo,hi,ngdiff),depth)

      double precision rhs(CELL3d(lo,hi,ngrhs),depth)

      double precision quat(CELL3d(lo,hi,ngq),depth)

      double precision facecoeffx(SIDE3d0(lo,hi,ngfacecoeff))
      double precision facecoeffy(SIDE3d1(lo,hi,ngfacecoeff))
      double precision facecoeffz(SIDE3d2(lo,hi,ngfacecoeff))
      
      double precision mobility(CELL3d(lo,hi,ngmob))
      
      integer iqrot_x(SIDE3d0(lo,hi,ngiq))
      integer iqrot_y(SIDE3d1(lo,hi,ngiq))
      integer iqrot_z(SIDE3d2(lo,hi,ngiq))

      integer i, j, k, m, iq
      double precision beta, lambda
      double precision invdx2(NDIM)
      double precision tmp(depth)
      double precision dtmp_x(depth), dtmp_y(depth), dtmp_z(depth)
      double precision dprime_x(depth), dprime_y(depth), dprime_z(depth)
      
      invdx2(1) = 1.d0 / ( dx(1) * dx(1) )
      invdx2(2) = 1.d0 / ( dx(2) * dx(2) )
      invdx2(3) = 1.d0 / ( dx(3) * dx(3) )

      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0         
               
               if ( depth > 1 ) then

                  do m = 1, depth
                     dtmp_x(m) = symm_diff_x(i+1,j,k,m)
                     dtmp_y(m) = symm_diff_y(i,j+1,k,m)
                     dtmp_z(m) = symm_diff_z(i,j,k+1,m)
                  enddo
                  
                  iq = -1 * iqrot_x(i+1,j,k)
                  call quatsymmrotate( dtmp_x, iq, dprime_x, depth )

                  iq = -1 * iqrot_y(i,j+1,k)
                  call quatsymmrotate( dtmp_y, iq, dprime_y, depth )

                  iq = -1 * iqrot_z(i,j,k+1)
                  call quatsymmrotate( dtmp_z, iq, dprime_z, depth )

               else

                  dprime_x(1) = symm_diff_x(i+1,j,k,1)
                  dprime_y(1) = symm_diff_y(i,j+1,k,1)
                  dprime_z(1) = symm_diff_z(i,j,k+1,1)

               endif

               do m = 1, depth

                  tmp(m) = 
     &               invdx2(1) * (
     &                  facecoeffx(i+1,j,k) *
     &                     ( nonsymm_diff_x(i+1,j,k,m) -
     &                       dprime_x(m) )
     &                  -
     &                  facecoeffx(i,j,k) * 
     &                     ( nonsymm_diff_x(i,j,k,m) -
     &                       symm_diff_x(i,j,k,m) )
     &               )
     &               +
     &               invdx2(2) * (
     &                  facecoeffy(i,j+1,k) *
     &                  ( nonsymm_diff_y(i,j+1,k,m) - 
     &                    dprime_y(m) )
     &                  -
     &                  facecoeffy(i,j,k) *
     &                  ( nonsymm_diff_y(i,j,k,m) - 
     &                    symm_diff_y(i,j,k,m) )
     &               )
     &               +
     &               invdx2(3) * (
     &                  facecoeffz(i,j,k+1) *
     &                  ( nonsymm_diff_z(i,j,k+1,m) - 
     &                    dprime_z(m) )
     &                  -
     &                  facecoeffz(i,j,k) *
     &                  ( nonsymm_diff_z(i,j,k,m) - 
     &                    symm_diff_z(i,j,k,m) )
     &               )

               enddo

c subtract projection of tmp along quat
               if ( depth > 1 ) then

                  beta = 0.d0
                  lambda = 0.d0

                  do m = 1, depth
                     beta   = beta   + quat(i,j,k,m) * quat(i,j,k,m)
                     lambda = lambda + quat(i,j,k,m) * tmp(m)
                  enddo

                  lambda = lambda / beta

                  do m = 1, depth
                     rhs(i,j,k,m) = rhs(i,j,k,m) +
     &                  mobility(i,j,k) *
     &                  ( tmp(m) - lambda * quat(i,j,k,m) )
                  enddo

               else

                  do m = 1, depth
                     rhs(i,j,k,m) = rhs(i,j,k,m) +
     &                  mobility(i,j,k) * tmp(m)
                  enddo

               endif

            enddo
         enddo
      enddo

      return
      end
