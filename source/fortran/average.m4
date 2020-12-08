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

      end

