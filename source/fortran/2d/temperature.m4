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
     &   conc, ngc, nc,
     &   temp, ngt,
     &   cp, ncp,
     &   cp_powers, npow,
     &   cpcoeffs, ncoeffs)

      implicit none
      integer ifirst0, ilast0, ifirst1, ilast1,
     &        ngc, nc, ngt, ncp, npow, ncoeffs
      integer cp_powers(npow)
      REAL 
     &     conc(CELL2d(ifirst,ilast,ngc),nc),
     &     temp(CELL2d(ifirst,ilast,ngt)),
     &     cp(CELL2d(ifirst,ilast,ncp)),
     &     cpcoeffs(ncoeffs)

c     local variables
      integer i, j, p, m
      REAL mc, te, cpval
 
      do j = ifirst1, ilast1
         do i = ifirst0, ilast0
            te=temp(i,j)
            
            mc = 1.d0
            cp(i,j) = 0.d0

c loop over species
            do m = 1, nc
               cpval=0.d0
               do p=1, npow
                  cpval=cpval+cpcoeffs((m-1)*npow+p)*te**cp_powers(p)
               enddo
               mc = mc - conc(i,j,m)
               cp(i,j) = cp(i,j) + cpval*conc(i,j,m)
            enddo

            cpval=0.d0
            do p=1, npow
               cpval=cpval+cpcoeffs(nc*npow+p)*te**cp_powers(p)
            enddo
            cp(i,j) = cp(i,j) + cpval*mc

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
c            temp(i,j)=tref

         enddo
      enddo
c
      return
      end

c***********************************************************************

      subroutine initgradient(dx,xlo,xhi,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,
     &  var,
     &  center,temperature_center,
     &  gradient)
c***********************************************************************
      implicit none
c***********************************************************************
      REAL half
      parameter (half=.5d0)

      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0,gcw1
      REAL center(0:NDIM-1)
      REAL temperature_center
      REAL dx(0:NDIM-1),xlo(0:NDIM-1),xhi(0:NDIM-1)
      REAL gradient(0:NDIM-1)
c variables in 2d cell indexed
      REAL var(CELL2dVECG(ifirst,ilast,gcw))

c***********************************************************************
      integer ic0,ic1
      REAL xc(0:NDIM-1),x0,x1

      do ic1=ifirst1-gcw1,ilast1+gcw1
        xc(1) = xlo(1)+dx(1)*(dble(ic1-ifirst1)+half)
        x1 = xc(1)-center(1)
        do ic0=ifirst0-gcw0,ilast0+gcw0
           xc(0) = xlo(0)+dx(0)*(dble(ic0-ifirst0)+half)
           x0 = xc(0)-center(0)

           var(ic0,ic1) = temperature_center
     &                  + x0*gradient(0)
     &                  + x1*gradient(1)

         enddo
      enddo

      return
      end

