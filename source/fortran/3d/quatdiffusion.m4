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
      subroutine quatdiffusionderiv(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  misorientation_factor,
     &  temperature, tghosts,
     &  var, ngvar,
     &  depth,
     &  gradq0, gradq1, gradq2, nggradq,
     &  diff0, diff1, diff2, ngdiff,
     &  gradient_floor, floor_type,
     &  diff_type,
     &  avg_type )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      double precision    misorientation_factor, gradient_floor
      integer ngvar, nggradq, ngdiff, depth, tghosts
      character*(*) floor_type
c
c variables in 3d cell indexed
      double precision var(CELL3d(ifirst,ilast,ngvar))
      double precision temperature(CELL3d(ifirst,ilast,tghosts))
      double precision gradq0(SIDE3d0(ifirst,ilast,nggradq),depth,NDIM)
      double precision gradq1(SIDE3d1(ifirst,ilast,nggradq),depth,NDIM)
      double precision gradq2(SIDE3d2(ifirst,ilast,nggradq),depth,NDIM)
      double precision diff0(SIDE3d0(ifirst,ilast,ngdiff),2)
      double precision diff1(SIDE3d1(ifirst,ilast,ngdiff),2)
      double precision diff2(SIDE3d2(ifirst,ilast,ngdiff),2)
      character*(*) diff_type
      character*(*) avg_type
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0,ic1,ic2,m,n
      double precision phi, dphi, jlf_threshold, d_deriv,
     &     fac, norm_gradq,
     &     floor_grad_norm2, max_grad_normi,
     &     grad_norm2, grad_normi
      parameter( jlf_threshold=1.0d-16 )
      double precision deriv_interp_func
      double precision average_func
      double precision deriv_average_func
      double precision t 

      double precision eval_grad_normi

      floor_grad_norm2 = gradient_floor**2
c
      max_grad_normi = 1.d0 / gradient_floor
c
      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            do ic0=ifirst0,ilast0+1
               if (var(ic0-1,ic1,ic2) .lt. jlf_threshold .or.
     &             var(ic0  ,ic1,ic2) .lt. jlf_threshold ) then
                  do m = 1, depth
                     diff0(ic0,ic1,ic2,1) = 0.d0
                     diff0(ic0,ic1,ic2,2) = 0.d0
                  end do
               else
                  phi = average_func(var(ic0-1,ic1,ic2),
     &                               var(ic0,  ic1,ic2),
     &                               avg_type)

                  t = 0.5 * (
     &               temperature(ic0-1,ic1,ic2) +
     &               temperature(ic0,ic1,ic2) )

                  d_deriv = misorientation_factor * t *
     &               deriv_interp_func( phi, diff_type )

                  grad_norm2 = 0.d0
                  do n = 1, NDIM
                     do m = 1, depth
                        grad_norm2 = grad_norm2 + 
     &                       gradq0(ic0,ic1,ic2,m,n)**2
                     enddo
                  enddo

                  grad_normi = eval_grad_normi(grad_norm2, floor_type, 
     &                                   floor_grad_norm2, 
     &                                   max_grad_normi)

                  fac = grad_normi * d_deriv

                  dphi = deriv_average_func( phi, var(ic0-1,ic1,ic2),
     &                                       avg_type)
                  diff0(ic0,ic1,ic2,1) = fac * dphi

                  dphi = deriv_average_func( phi, var(ic0,ic1,ic2),
     &                                       avg_type)
                  diff0(ic0,ic1,ic2,2) = fac * dphi
               end if
            end do
         end do
      end do
c
      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1+1
            do ic0=ifirst0,ilast0
               if (var(ic0,ic1-1,ic2) .lt. jlf_threshold .or.
     &             var(ic0,ic1  ,ic2) .lt. jlf_threshold ) then
                  do m = 1, depth
                     diff1(ic0,ic1,ic2,1) = 0.d0
                     diff1(ic0,ic1,ic2,2) = 0.d0
                  end do
               else
                  phi = average_func(var(ic0,ic1-1,ic2),
     &                               var(ic0,ic1  ,ic2),
     &                               avg_type)
            
                  t = 0.5 * (
     &               temperature(ic0,ic1-1,ic2) +
     &               temperature(ic0,ic1,ic2) )

                  d_deriv = misorientation_factor * t *
     &               deriv_interp_func( phi, diff_type )

                  grad_norm2 = 0.d0
                  do n = 1, NDIM
                     do m = 1, depth
                        grad_norm2 = grad_norm2 +
     &                       gradq1(ic0,ic1,ic2,m,n)**2
                     enddo
                  enddo

                  grad_normi = eval_grad_normi(grad_norm2, floor_type, 
     &                                   floor_grad_norm2, 
     &                                   max_grad_normi)

                  fac = grad_normi * d_deriv

                  dphi = deriv_average_func( phi, var(ic0,ic1-1,ic2),
     &                                       avg_type)
                  diff1(ic0,ic1,ic2,1) = fac * dphi

                  dphi = deriv_average_func( phi, var(ic0,ic1,ic2),
     &                                       avg_type)
                  diff1(ic0,ic1,ic2,2) = fac * dphi

               end if
            end do
         end do
      end do
c
      do ic2=ifirst2,ilast2+1
         do ic1=ifirst1,ilast1
            do ic0=ifirst0,ilast0
               if (var(ic0,ic1,ic2-1) .lt. jlf_threshold .or.
     &             var(ic0,ic1,ic2  ) .lt. jlf_threshold ) then
                  do m = 1, depth
                     diff2(ic0,ic1,ic2,1) = 0.d0
                     diff2(ic0,ic1,ic2,2) = 0.d0
                  end do
               else

                  phi = average_func(var(ic0,ic1,ic2-1),
     &                               var(ic0,ic1,ic2  ),
     &                               avg_type)
            
                  t = 0.5 * (
     &               temperature(ic0,ic1,ic2-1) +
     &               temperature(ic0,ic1,ic2) )

                  d_deriv = misorientation_factor * t *
     &               deriv_interp_func( phi, diff_type )

                  grad_norm2 = 0.d0
                  do n = 1, NDIM
                     do m = 1, depth
                        grad_norm2 = grad_norm2 +
     &                       gradq2(ic0,ic1,ic2,m,n)**2
                     enddo
                  enddo
                  
                  grad_normi = eval_grad_normi(grad_norm2, floor_type,
     &                                   floor_grad_norm2,
     &                                   max_grad_normi)

                  fac = grad_normi * d_deriv

                  dphi = deriv_average_func( phi, var(ic0,ic1,ic2-1),
     &                                       avg_type)
                  diff2(ic0,ic1,ic2,1) = fac * dphi

                  dphi = deriv_average_func( phi, var(ic0,ic1,ic2),
     &                                       avg_type)
                  diff2(ic0,ic1,ic2,2) = fac * dphi

               end if
            end do
         end do
      end do
c
      return
      end
