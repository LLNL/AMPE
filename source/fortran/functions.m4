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
c        These routines are written for simplicity, flexibility and ease
c        of maintenance for now, but in the future they may need to be
c        reworked completely for performance reasons.

c=======================================================================

c        Valid values for type are "quadratic" and "pbg"

      function interp_func( phi, type )

      implicit none

      double precision interp_func
      double precision phi, phit
      character*(*) type

      double precision gamma
      gamma = 10.d0
      
      if ( type(1:1) .eq. 'q' .or.
     &     type(1:1) .eq. 'Q' ) then

         phit = max( 0.d0, phi )

         interp_func = phit * phit

      else if ( type(1:1) .eq. 'p' .or.
     &          type(1:1) .eq. 'P') then

         phit = max( 0.d0, min( 1.d0, phi ) )
         
         interp_func = 
     &      phit * phit * phit *
     &      ( 10.d0 - 15.d0 * phit + 6.d0 * phit * phit )

      else if ( type(1:1) .eq. 't' .or.
     &          type(1:1) .eq. 'T') then

         phit = max( 0.05d0, min( 1.d0, phi ) )
         
         interp_func = 
     &      phit * phit * phit *
     &      ( 10.d0 - 15.d0 * phit + 6.d0 * phit * phit )

      else if ( type(1:1) .eq. 'h' .or.
     &          type(1:1) .eq. 'H') then

         phit = max( 0.d0, min( 1.d0, phi ) )
         
         interp_func =
     &      phit * phit *
     &      ( 3.d0 - 2.d0 * phit )

      else if ( type(1:1) .eq. 'w' .or.
     &          type(1:1) .eq. 'W') then

         phit = max( 0.d0, phi )

         interp_func = phit * phit * ( 2.d0 - phit )

      else if ( type(1:1) .eq. 'l' .or.
     &          type(1:1) .eq. 'L' ) then

         phit = max( 0.d0, phi )
         phit = min( 1.d0, phit )
         interp_func = phit

      else if ( type(1:1) .eq. '3'  ) then

         phit = max( 0.d0, phi )
         interp_func = phit*phit*phit

      else if ( type(1:1) .eq. 's' .or.
     &          type(1:1) .eq. 'S' ) then

         phit = max( 0.d0, phi )
         interp_func = log ( cosh (gamma*phi) )/log(cosh(gamma))

      else if ( type(1:1) .eq. 'c' .or.
     &          type(1:1) .eq. 'C' ) then

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

      double precision gamma
      gamma = 10.d0
      
      if ( type(1:1) .eq. 'q' .or.
     &     type(1:1) .eq. 'Q' ) then

         phit = max( 0.d0, phi )

         deriv_interp_func = 2.d0 * phit

      else if ( type(1:1) .eq. 'p' .or.
     &          type(1:1) .eq. 'P') then

         phit = max( 0.d0, min( 1.d0, phi ) )
         
         deriv_interp_func = 
     &      30.d0 * phit * phit *
     &      ( 1.d0 - phit  ) * ( 1.d0 - phit )

      else if ( type(1:1) .eq. 't' .or.
     &          type(1:1) .eq. 'T') then

         phit = max( 0.05d0, min( 1.d0, phi ) )
         
         deriv_interp_func = 
     &      30.d0 * phit * phit *
     &      ( 1.d0 - phit  ) * ( 1.d0 - phit )

      else if ( type(1:1) .eq. 'h' .or.
     &          type(1:1) .eq. 'H') then

         phit = max( 0.d0, min( 1.d0, phi ) )

         deriv_interp_func =
     &      6.d0 * phit * ( 1.d0 - phit )

      else if ( type(1:1) .eq. 'w' .or.
     &          type(1:1) .eq. 'W') then

         phit = max( 0.d0, phi )

         deriv_interp_func = phit * ( 4.d0 - 3.d0 * phit )

c      else if ( type(1:1) .eq. 'l' .or.
c     &          type(1:1) .eq. 'L' ) then
c
c         deriv_interp_func = 1.d0
c
      else if ( type(1:1) .eq. '3' ) then

         phit = max( 0.d0, phi )
         deriv_interp_func = 3.d0*phit*phit

      else if ( type(1:1) .eq. 's' .or.
     &          type(1:1) .eq. 'S' ) then

         phit = max( 0.d0, phi )
         deriv_interp_func = gamma*tanh(gamma*phit)/log(cosh(gamma))

      else if ( type(1:1) .eq. 'c' .or.
     &          type(1:1) .eq. 'C' ) then

         deriv_interp_func = 0.d0

      else

         print *, "Error in deriv_interp_func: type unknown"
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

      if ( type(1:1) .eq. 'q' .or.
     &     type(1:1) .eq. 'Q' ) then

         second_deriv_interp_func = 2.d0

      else if ( type(1:1) .eq. 'p' .or.
     &          type(1:1) .eq. 'P') then

         phit = max( 0.d0, min( 1.d0, phi ) )
         
         second_deriv_interp_func = 
     &      60.d0 * phit *
     &      ( 1.d0 - 3.d0 * phit + 2.d0 * phit * phit )

      else if ( type(1:1) .eq. 'h' .or.
     &          type(1:1) .eq. 'H') then

         phit = max( 0.d0, min( 1.d0, phi ) )

         second_deriv_interp_func =
     &      6.d0 * ( 1.d0 - 2.d0 * phit )

      else if ( type(1:1) .eq. 'w' .or.
     &          type(1:1) .eq. 'W') then

         phit = max( 0.d0, phi )

         second_deriv_interp_func = 4.d0 - 6.d0 * phit

      else if ( type(1:1) .eq. 'l' .or.
     &          type(1:1) .eq. 'L' ) then

         second_deriv_interp_func = 0.d0

      else if ( type(1:1) .eq. 'c' .or.
     &          type(1:1) .eq. 'C' ) then

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

      if ( type(1:1) .eq. 's' .or.
     &     type(1:1) .eq. 'S' ) then

         well_func = ( 1.d0 - phi ) * ( 1.d0 - phi )

      else if ( type(1:1) .eq. 'd' .or.
     &          type(1:1) .eq. 'D') then

         well_func =
     &      16.d0 * phi * phi * ( 1.d0 - phi ) * ( 1.d0 - phi )

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

      if ( type(1:1) .eq. 's' .or.
     &     type(1:1) .eq. 'S' ) then

         deriv_well_func = 2.d0 * ( phi - 1.d0 )

      else if ( type(1:1) .eq. 'd' .or.
     &          type(1:1) .eq. 'D') then

         deriv_well_func =
     &      32.d0 * phi * ( 1.d0 - phi ) * ( 1.d0 - 2.d0 * phi )

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

      if ( type(1:1) .eq. 's' .or.
     &     type(1:1) .eq. 'S' ) then

         second_deriv_well_func = 2.d0

      else if ( type(1:1) .eq. 'd' .or.
     &          type(1:1) .eq. 'D') then

         second_deriv_well_func =
     &      32.d0 * ( 1.d0 + 6.d0 * phi * ( phi - 1.d0 ) )

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
