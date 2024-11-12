c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
c 
c        These routines are written for simplicity, flexibility and ease
c        of maintenance for now, but in the future they may need to be
c        reworked completely for performance reasons.

c=======================================================================

c Valid values for type are "quadratic", "pbg", "harmonic",...

      function interp_func( phi, type )

      implicit none

      double precision interp_func
      double precision phi, phit
      character*(*) type

      double precision gamma
      gamma = 10.d0

      if ( type(1:1) .eq. 'q' ) then

         phit = max( 0.d0, phi )

         interp_func = phit * phit

      else if ( type(1:1) .eq. 'p') then 

         phit = max( 0.d0, min( 1.d0, phi ) )
         
         interp_func = 
     &      phit * phit * phit *
     &      ( 10.d0 - 15.d0 * phit + 6.d0 * phit * phit )

      else if ( type(1:1) .eq. 'h') then

         phit = max( 0.d0, min( 1.d0, phi ) )
         
         interp_func =
     &      phit * phit *
     &      ( 3.d0 - 2.d0 * phit )

      else if ( type(1:1) .eq. 'w') then

         phit = max( 0.d0, phi )

         interp_func = phit * phit * ( 2.d0 - phit )

      else if ( type(1:1) .eq. 'l' ) then

         phit = max( 0.d0, phi )
         phit = min( 1.d0, phit )
         interp_func = phit

      else if ( type(1:1) .eq. 'm' ) then
c Moleans 2011
         phit = max( 0.d0, min( 1.d0, phi ) )

         interp_func =
     &      phit * phit / ( 2.d0 * phit * (phit -1.d0) + 1.d0 )

      else if ( type(1:1) .eq. '3'  ) then

         phit = max( 0.d0, phi )
         interp_func = phit*phit*phit

      else if ( type(1:1) .eq. 's' ) then

         phit = max( 0.d0, phi )
         interp_func = log ( cosh (gamma*phi) )/log(cosh(gamma))

      else if ( type(1:1) .eq. 'c' ) then

         interp_func = 1.d0

      else

         print *, "Error in interp_func: type unknown"
         stop

      endif

      return
      end

c-----------------------------------------------------------------------

      function deriv_interp_func( phi, type )

      implicit none

      double precision deriv_interp_func
      double precision phi, phit
      
      character*(*) type

      double precision gamma, tmp
      gamma = 10.d0
      
      if ( type(1:1) .eq. 'q' ) then

         phit = max( 0.d0, phi )

         deriv_interp_func = 2.d0 * phit

      else if ( type(1:1) .eq. 'p' ) then

         phit = max( 0.d0, min( 1.d0, phi ) )
         
         deriv_interp_func = 
     &      30.d0 * phit * phit *
     &      ( 1.d0 - phit  ) * ( 1.d0 - phit )

      else if ( type(1:1) .eq. 'h' ) then

         phit = max( 0.d0, min( 1.d0, phi ) )

         deriv_interp_func =
     &      6.d0 * phit * ( 1.d0 - phit )

      else if ( type(1:1) .eq. 'w' ) then

         phit = max( 0.d0, phi )

         deriv_interp_func = phit * ( 4.d0 - 3.d0 * phit )

      else if ( type(1:1) .eq. 'l' ) then

         if( phi .gt. 0.d0 .or.
     &       phi .lt. 1.d0 )then
            deriv_interp_func = 1.d0
         else
             deriv_interp_func = 0.d0
         endif

      else if ( type(1:1) .eq. 'm' ) then

         phit = max( 0.d0, min( 1.d0, phi ) )
         tmp =  2.d0*phit*(phit-1.d0)+1.d0
         deriv_interp_func = 2.d0*phit*(1.d0-phit)
     &      / (tmp*tmp)

      else if ( type(1:1) .eq. '3' ) then

         phit = max( 0.d0, phi )
         deriv_interp_func = 3.d0*phit*phit

      else if ( type(1:1) .eq. 's' ) then

         phit = max( 0.d0, phi )
         deriv_interp_func = gamma*tanh(gamma*phit)/log(cosh(gamma))

      else if ( type(1:1) .eq. 'c' ) then

         deriv_interp_func = 0.d0

      else

         print *, "Error in deriv_interp_func: unknown type ",type(1:1)
         stop

      endif

      return
      end

c=======================================================================

      function second_deriv_interp_func( phi, type )

      implicit none

      double precision second_deriv_interp_func
      double precision phi, phit
      character*(*) type

      if ( type(1:1) .eq. 'q' ) then

         if( phi .ge. 0.d0 )then
            second_deriv_interp_func = 2.d0
         else
            second_deriv_interp_func = 0.d0
         endif

      else if ( type(1:1) .eq. 'p' ) then

         phit = max( 0.d0, min( 1.d0, phi ) )
         
         second_deriv_interp_func = 
     &      60.d0 * phit *
     &      ( 1.d0 - 3.d0 * phit + 2.d0 * phit * phit )

      else if ( type(1:1) .eq. 'h' ) then

         phit = max( 0.d0, min( 1.d0, phi ) )

         second_deriv_interp_func =
     &      6.d0 * ( 1.d0 - 2.d0 * phit )

      else if ( type(1:1) .eq. 'w' ) then

         phit = max( 0.d0, phi )

         second_deriv_interp_func = 4.d0 - 6.d0 * phit

      else if ( type(1:1) .eq. 'l' ) then

         second_deriv_interp_func = 0.d0

      else if ( type(1:1) .eq. 'a' ) then

         if( phi < 0.08333333333333333d0 )then
            phit = max( 0.d0, phi )
            second_deriv_interp_func = 324.d0*phit
         else if( phi>0.9166666666666666d0 )then
            phit = min( 1.d0, phi )
            phit = 1.d0-phit
            second_deriv_interp_func = -324.d0*phit
         else
            second_deriv_interp_func = 0.d0
         endif

      else if ( type(1:1) .eq. 'c' ) then

         second_deriv_interp_func = 0.d0

      else

         print *, "Error in deriv_interp_func: type unknown"
         stop

      endif

      return
      end

c=======================================================================

      function well_func( phi, type )

      implicit none

      double precision well_func
      double precision phi
      character*(*) type

      if ( type(1:1) .eq. 'd' ) then

         well_func =
     &      16.d0 * phi * phi * ( 1.d0 - phi ) * ( 1.d0 - phi )

      else if ( type(1:1) .eq. 's' ) then

         well_func = ( 1.d0 - phi ) * ( 1.d0 - phi )

      else

         print *, "Error in well_func: type unknown"
         stop

      endif

      return
      end

c-----------------------------------------------------------------------

      function deriv_well_func( phi, type )

      implicit none

      double precision deriv_well_func
      double precision phi
      character*(*) type

      if ( type(1:1) .eq. 'd' ) then

         deriv_well_func =
     &      32.d0 * phi * ( 1.d0 - phi ) * ( 1.d0 - 2.d0 * phi )

      else if ( type(1:1) .eq. 's' ) then

         deriv_well_func = 2.d0 * ( phi - 1.d0 )

      else

         print *, "Error in deriv_well_func: type unknown"
         stop

      endif

      return
      end

c-----------------------------------------------------------------------

      function second_deriv_well_func( phi, type )

      implicit none

      double precision second_deriv_well_func
      double precision phi
      character*(*) type

      if ( type(1:1) .eq. 'd' ) then

         second_deriv_well_func =
     &      32.d0 * ( 1.d0 + 6.d0 * phi * ( phi - 1.d0 ) )

      else if ( type(1:1) .eq. 's' ) then

         second_deriv_well_func = 2.d0

      else

         print *, "Error in second_deriv_well_func: type unknown"
         stop

      endif

      return
      end

c=======================================================================

      function average_func( phi1, phi2, avg_type )

      implicit none

      double precision average_func
      double precision phi1,phi2
      double precision threshold
      character*(*) avg_type
      parameter( threshold=1.0d-16 )
      
      if ( avg_type(1:1) .eq. 'a' ) then

c     arithmetic average
         average_func = 0.5d0 * (phi1+phi2)

      else if ( avg_type(1:1) .eq. 'h') then

c     harmonic average
         if (phi1 .lt. threshold .or.
     &       phi2 .lt. threshold ) then
            average_func = 0.d0
         else
            average_func = 2.d0 / (1.d0/phi1 +
     &                             1.d0/phi2)
         endif

      else

         print *, "Error in average_func: type unknown"
         stop

      endif
      
      return
      end

c-----------------------------------------------------------------------

      function deriv_average_func( avg_phi, next_phi, avg_type )

      implicit none

      double precision deriv_average_func
      double precision avg_phi,next_phi
      double precision threshold
      character*(*) avg_type
      parameter( threshold=1.0d-16 )
      
      if ( avg_type(1:1) .eq. 'a' ) then

         deriv_average_func = 0.5d0

      else if ( avg_type(1:1) .eq. 'h') then

         if ( avg_phi .lt. threshold ) then
            deriv_average_func = 0.d0
         else
            deriv_average_func = 0.5d0 * next_phi * next_phi
     &                                 / (avg_phi * avg_phi)
         endif

      else

         print *, "Error in deriv_average_func: type unknown"
         stop

      endif
      
      return
      end

c-----------------------------------------------------------------------
c mapping from "double index" to "single index" in stiffness matrices
c 11 -> 1, 22 -> 2, 12 -> 3
c (Voigt notation extended to 2D)
      function single_index_stiffness2d( i, j )

      implicit none
      
      integer single_index_stiffness2d
      integer i,j
      
      if( i .eq. j ) then
         single_index_stiffness2d = i
      else
         single_index_stiffness2d = 3
      endif

      return
      end

c-----------------------------------------------------------------------
c mapping from "double index" to "single index" in stiffness matrices
c 11 -> 1, 22 -> 2, 33 -> 3
c 23,32 -> 4, 31,13 -> 5, 12,21 -> 6
c (Voigt notation)
      function single_index_stiffness3d( i, j )

      implicit none
      
      integer single_index_stiffness3d
      integer i,j
      
      if( i .eq. j ) then
         single_index_stiffness3d = i
      else
         single_index_stiffness3d = -i -j + 9
      endif

      return
      end

c-----------------------------------------------------------------------
c get storage location for element (i,j) in compact storage
c of symmetric stiffness matrix
c 1<=i,j<=3
c
      function storage_index_stiffness2d( i, j )

      implicit none
      
      integer storage_index_stiffness2d
      integer i,j
      integer ii,jj,n
      
      if( i .ge. j ) then
        ii=j
        jj=i
      else
        ii=i
        jj=j
      endif
      
      n=4-ii
      storage_index_stiffness2d = 7 - (n*(n+1))/2 + (jj-ii)

      return
      end

c-----------------------------------------------------------------------
c get storage location for element (i,j) in compact storage
c of symmetric stiffness matrix
c 1<=i,j<=6
c
      function storage_index_stiffness3d( i, j )

      implicit none
      
      integer storage_index_stiffness3d
      integer i,j
      integer ii,jj,n
      
      if( i .ge. j ) then
        ii=j
        jj=i
      else
        ii=i
        jj=j
      endif
      
      n=7-ii
      storage_index_stiffness3d = 22 - (n*(n+1))/2 + (jj-ii)

      return
      end

c-----------------------------------------------------------------------

      function limitedinverse( val )

      implicit none
      
      double precision limitedinverse
      double precision val
      
      double precision largeval
      double precision threshold   
      parameter( threshold=1.0d-4 )
      parameter( largeval=1.d4 )
      
      if(abs(val)>threshold)then
         limitedinverse=1.d0/val
      else
         limitedinverse=largeval
      endif
      
      return
      end

c-----------------------------------------------------------------------

      function interp_ratio_func(phi, type1, type2)

      implicit none

      double precision interp_ratio_func
      double precision phi, phit
      character*(*) type1, type2

      if ( type1(1:1) .eq. type2(1:1) )then

         interp_ratio_func = 1.d0

      else if ( type1(1:1) .eq. 'p' )then

         if ( type2(1:1) .eq. 'l' )then

            phit = max( 0.d0, min( 1.d0, phi ) )

            interp_ratio_func = phit * phit * 
     &                     ( 10.d0 - 15.d0 * phit + 6.d0 * phit * phit )

         else

            print *, "Error, interp_ratio: unknown/incompatible types"
            stop

         endif
      else

         print *, "Error, interp_ratio: unknown/incompatible types"
         stop

      endif

      return
      end

c-----------------------------------------------------------------------

      function compl_interp_ratio_func(phi, type1, type2)

      implicit none

      double precision compl_interp_ratio_func
      double precision phi, phit
      character*(*) type1, type2

      if ( type1(1:1) .eq. type2(1:1) )then

         compl_interp_ratio_func = 1.d0

      else if ( type1(1:1) .eq. 'p' )then

         if ( type2(1:1) .eq. 'p' )then

            compl_interp_ratio_func = 1.d0

         else if ( type2(1:1) .eq. 'l' )then

            phit = max( 0.d0, min( 1.d0, phi ) )

            compl_interp_ratio_func = (1.d0-phit) * (1.d0-phit) *
     &                     ( 1.d0 + 3.d0 * phit + 6.d0 * phit * phit )

         else

            print *, "Error:"
            print *, "compl_interp_ratio: unknown/incompatible types"
            stop

         endif
      else

         print *, "Error:"
         print *, "compl_interp_ratio: unknown/incompatible types"
         stop

      endif

      return
      end

c-----------------------------------------------------------------------
      function triple_well_func(p0, p1, p2)

      implicit none

      double precision triple_well_func
      double precision p0, p1, p2

      triple_well_func = 16.d0 * (p0*p0*(1.d0-p0)*(1.d0-p0)
     &                         +  p1*p1*(1.d0-p1)*(1.d0-p1)
     &                         +  p2*p2*(1.d0-p2)*(1.d0-p2))

      return
      end

c-----------------------------------------------------------------------
c dF/dp0 - 1/3 * sum_j dF/dpj
c
      function deriv_triple_well_func(p0, p1, p2)

      implicit none

      double precision deriv_triple_well_func
      double precision p0, p1, p2

      deriv_triple_well_func = 32.d0 * (1.d0 / 3.d0) *
     &    (2.d0 * p0 * (1.d0 - p0) * (1.d0 - 2.d0 * p0)
     &     - p1 * (1.d0 - p1) * (1.d0 - 2.d0 * p1)
     &     - p2 * (1.d0 - p2) * (1.d0 - 2.d0 * p2));

      return
      end

