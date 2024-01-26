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
c
      subroutine quatmobility(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  phase,
     &  ngphase,
     &  mobility,
     &  ngmobility,
     &  scale_mobility,
     &  min_mobility,
     &  func_type,
     &  alt_scale_factor )
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer ngphase,ngmobility
      double precision phase(CELL2d(ifirst,ilast,ngphase))
c ouput arrays:
      double precision mobility(CELL2d(ifirst,ilast,ngmobility))
      double precision scale_mobility, min_mobility
      double precision alt_scale_factor
      character*(*) func_type
c
c***********************************************************************
c
      integer i0,i1
      double precision phi, qfunc
      double precision c
c
      do i1=ifirst1-ngmobility,ilast1+ngmobility
         do i0=ifirst0-ngmobility,ilast0+ngmobility
            
            phi = phase(i0,i1)

            if ( func_type(1:1) == 'p' .or.
     &           func_type(1:1) == 'P' ) then

               phi = max( 0.d0, min( 1.d0, phi ) )
            
               qfunc = 
     &            phi * phi * phi *
     &            ( 10.d0 - 15.d0 * phi +
     &            6.d0 * phi * phi )

               qfunc = 1.d0 - qfunc

            elseif ( func_type(1:1) == 'e' .or.
     &               func_type(1:1) == 'E' ) then

c               c = 10.d0
               c = alt_scale_factor

               phi = max( 0.d0, min( 1.d0, phi ) )
            
               qfunc = ( 1.d0 - exp( c * phi ) ) / ( 1.d0 - exp( c ) )

               qfunc = 1.d0 - qfunc

            elseif ( func_type(1:1) == 'i' .or.
     &               func_type(1:1) == 'I' ) then

               phi = max( 1.d-6, min( 1.d0, phi ) )
               
               qfunc = max( 0.0d0, ( 1.d0 - phi ) / ( phi * phi ) )

               qfunc = min( qfunc, alt_scale_factor )

            else

               print *, "Error in quatmobility: unknown function type"
               stop

            endif

            mobility(i0,i1) = min_mobility +
     &         (scale_mobility-min_mobility) * qfunc

         enddo
      enddo
c
      return
      end
c
c
      subroutine quatmobilityderiv(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  phase,
     &  ngphase,
     &  dmobility,
     &  ngmobility,
     &  scale_mobility,
     &  min_mobility,
     &  func_type,
     &  alt_scale_factor )
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer ngphase,ngmobility
      double precision phase(CELL2d(ifirst,ilast,ngphase))
c ouput arrays:
      double precision dmobility(CELL2d(ifirst,ilast,ngmobility))
      double precision scale_mobility, min_mobility
      double precision alt_scale_factor
      character*(*) func_type
c
c***********************************************************************
c
      integer i0,i1
      double precision phi, dqfunc
      double precision c
c
      do i1=ifirst1,ilast1
         do i0=ifirst0,ilast0
            
            phi = phase(i0,i1)

            if ( func_type(1:1) == 'p' .or.
     &           func_type(1:1) == 'P' ) then

               phi = max( 0.d0, min( 1.d0, phi ) )
            
               dqfunc = 
     &            -30.d0 * phi * phi *
     &            ( 1.d0 - phi  ) * ( 1.d0 - phi )

            elseif ( func_type(1:1) == 'e' .or.
     &               func_type(1:1) == 'E' ) then

*               c = 10.0
               c = alt_scale_factor

               phi = max( 0.d0, min( 1.d0, phi ) )
            
               dqfunc = ( c * exp( c * phi ) ) / ( 1. - exp( c ) )

            elseif ( func_type(1:1) == 'i' .or.
     &               func_type(1:1) == 'I' ) then

               phi = max( 1.d-6, min( 1.d0, phi ) )
               
               dqfunc = ( phi - 2.d0 ) / ( phi * phi * phi )

            else

               print *,
     &            "Error in quatmobilityderiv: unknown function type"
               stop

            endif

            dmobility(i0,i1) = 
     &         (scale_mobility-min_mobility) * dqfunc

         enddo
      enddo

      return
      end
