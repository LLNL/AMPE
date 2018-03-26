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
define(REAL,`double precision')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine source_temperature(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   conc, ngc,
     &   rhs, ngrhs,
     &   cp, ncp,
     &   source,
     &   nsp)

      implicit none
      integer ifirst0, ilast0, ifirst1, ilast1,
     &        ngc, ngrhs, nsp, ncp
      REAL 
     &     conc(CELL2d(ifirst,ilast,ngc)),
     &     rhs(CELL2d(ifirst,ilast,ngrhs)),
     &     cp(CELL2d(ifirst,ilast,ncp)),
     &     source(nsp), c

c     local variables
      integer i, j

      do j = ifirst1, ilast1
         do i = ifirst0, ilast0
            c=conc(i,j)
            rhs(i,j) = ( source(1)*c
     &                  +source(2)*(1.d0-c) )
     &                 /cp(i,j)
         enddo
      enddo
      
      return
      end

c***********************************************************************
      
      subroutine initgaussian(dx,xlo,xhi,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,
     &  var,
     &  center, ll,
     &  standard_dev, temperature_base, temperature_peak)
c***********************************************************************
      implicit none
c***********************************************************************     
c***********************************************************************     
      REAL half
      parameter (half=.5d0)
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0,gcw1
      REAL standard_dev,temperature_base, temperature_peak
      REAL radius2,factor,deltaT
      REAL center(0:NDIM-1),ll(0:NDIM-1)
      REAL dx(0:NDIM-1),xlo(0:NDIM-1),xhi(0:NDIM-1)
c
c variables in 2d cell indexed         
      REAL var(CELL2dVECG(ifirst,ilast,gcw))
c
c***********************************************************************     
c
      integer ic0,ic1
      REAL xc(0:NDIM-1),x0,x1
      
      deltaT = temperature_peak-temperature_base
      
      factor = 1./(2.*standard_dev*standard_dev)

      do ic1=ifirst1,ilast1
        xc(1) = xlo(1)+dx(1)*(dble(ic1-ifirst1)+half)
        x1 = abs(xc(1)-center(1))
        if( ll(1).gt.0 )then
           x1 = min (x1,abs(xc(1)-center(1)+ll(1)))
           x1 = min (x1,abs(xc(1)-center(1)-ll(1)))
        endif
        do ic0=ifirst0,ilast0
           xc(0) = xlo(0)+dx(0)*(dble(ic0-ifirst0)+half)
           x0 = abs(xc(0)-center(0))
           if( ll(0).gt.0 )then
              x0 = min (x0,abs(xc(0)-center(0)+ll(0)))
              x0 = min (x0,abs(xc(0)-center(0)-ll(0)))
           endif

           radius2=x0**2+x1**2
           var(ic0,ic1) = temperature_base
     &                  + deltaT*exp(-radius2*factor)

         enddo
      enddo
c
      return
      end

c***********************************************************************

      subroutine heat_capacity_nkr(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   conc, ngc,
     &   temp, ngt,
     &   cp, ncp,
     &   cp_powers, npowers,
     &   cpcoeffs, ncoeffs)

      implicit none
      integer ifirst0, ilast0, ifirst1, ilast1,
     &        ngc, ngt, ncp, npowers, ncoeffs
      integer cp_powers(npowers)
      REAL 
     &     conc(CELL2d(ifirst,ilast,ngc)),
     &     temp(CELL2d(ifirst,ilast,ngt)),
     &     cp(CELL2d(ifirst,ilast,ncp)),
     &     cpcoeffs(ncoeffs)

c     local variables
      integer i, j, p
      REAL c, t, cpvals(2)

      do j = ifirst1, ilast1
         do i = ifirst0, ilast0
            c=conc(i,j)
            t=temp(i,j)
            
            cpvals(1)=0.
            do p=1, npowers
              cpvals(1)=cpvals(1)+cpcoeffs(p)*T**cp_powers(p)
            enddo
            cpvals(2)=0.
            do p=1, npowers
              cpvals(2)=cpvals(2)+cpcoeffs(npowers+p)*T**cp_powers(p)
            enddo
            cp(i,j) = ( cpvals(1)*c
     &                 +cpvals(2)*(1.d0-c) )
         enddo
      enddo
      
      return
      end

c***********************************************************************
      
      subroutine initgaussiansource(dx,xlo,xhi,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  rhs, ngrhs,
     &  cp, ngcp,
     &  center, ll,
     &  standard_dev, heat_peak)
c***********************************************************************
      implicit none
c***********************************************************************     
c***********************************************************************     
      REAL half
      parameter (half=.5d0)
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer ngcp, ngrhs
      REAL standard_dev, heat_peak
      REAL radius2,factor
      REAL center(0:NDIM-1),ll(0:NDIM-1)
      REAL dx(0:NDIM-1),xlo(0:NDIM-1),xhi(0:NDIM-1)
c
c variables in 2d cell indexed         
      REAL cp(CELL2d(ifirst,ilast,ngcp))
      REAL rhs(CELL2d(ifirst,ilast,ngrhs))
c
c***********************************************************************     
c
      integer ic0,ic1
      REAL xc(0:NDIM-1),x0,x1
      
      factor = 1./(2.*standard_dev*standard_dev)

      do ic1=ifirst1,ilast1
        xc(1) = xlo(1)+dx(1)*(dble(ic1-ifirst1)+half)
        x1 = abs(xc(1)-center(1))
        if( ll(1).gt.0 )then
           x1 = min (x1,abs(xc(1)-center(1)+ll(1)))
           x1 = min (x1,abs(xc(1)-center(1)-ll(1)))
        endif
        do ic0=ifirst0,ilast0
           xc(0) = xlo(0)+dx(0)*(dble(ic0-ifirst0)+half)
           x0 = abs(xc(0)-center(0))
           if( ll(0).gt.0 )then
              x0 = min (x0,abs(xc(0)-center(0)+ll(0)))
              x0 = min (x0,abs(xc(0)-center(0)-ll(0)))
           endif

           radius2=x0**2+x1**2
           rhs(ic0,ic1) = heat_peak*exp(-radius2*factor)/cp(ic0,ic1)

         enddo
      enddo
c
      return
      end

c***********************************************************************
      
      subroutine linearmeltingline(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  temp, ngt,
     &  conc, ngc,
     &  tref, cref, slope)
c***********************************************************************
      implicit none
c***********************************************************************     
      integer ifirst0,ilast0,ifirst1,ilast1
      integer ngt, ngc
      REAL tref,cref,slope
c
c variables in 2d cell indexed         
      REAL conc(CELL2d(ifirst,ilast,ngc))
      REAL temp(CELL2d(ifirst,ilast,ngt))
c
c***********************************************************************     
c
      integer i,j
      
      do j = ifirst1, ilast1
         do i = ifirst0, ilast0
            
            temp(i,j)=tref+slope*(conc(i,j)-cref)

         enddo
      enddo
c
      return
      end

