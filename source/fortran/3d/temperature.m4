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
define(NDIM,3)dnl
define(REAL,`double precision')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine source_temperature(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   conc, ngc,
     &   rhs, ngrhs,
     &   cp, ncp,
     &   source,
     &   nsp)

      implicit none
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &        ngc, ngrhs, nsp, ncp
      double precision 
     &     conc(CELL3d(ifirst,ilast,ngc)),
     &     rhs(CELL3d(ifirst,ilast,ngrhs)),
     &     cp(CELL3d(ifirst,ilast,ncp)),
     &     source(nsp), c

c     local variables
      integer i, j, k
      REAL cpl

      do k = ifirst2, ilast2
         do j = ifirst1, ilast1
            do i = ifirst0, ilast0
               c=conc(i,j,k)
               cpl=cp(i,j,k)
               rhs(i,j,k) = ( source(1)*c
     &                       +source(2)*(1.d0-c) )
     &                      /cpl
            enddo
         enddo
      enddo
      
      return
      end
      
      subroutine initgaussian(dx,xlo,xhi,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  gcw0,gcw1,gcw2,
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
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer gcw0,gcw1,gcw2
      REAL standard_dev,temperature_base, temperature_peak
      REAL radius2,factor,deltaT
      REAL center(0:NDIM-1),ll(0:NDIM-1)
      REAL dx(0:NDIM-1),xlo(0:NDIM-1),xhi(0:NDIM-1)
c
c variables in 2d cell indexed         
      REAL var(CELL3dVECG(ifirst,ilast,gcw))
c
c***********************************************************************     
c
      integer ic0,ic1,ic2
      REAL xc(0:NDIM-1),x0,x1,x2
      
      deltaT = temperature_peak-temperature_base
      
      factor = 1./(2.*standard_dev*standard_dev)

      do ic2=ifirst2,ilast2
        xc(2) = xlo(2)+dx(2)*(dble(ic2-ifirst2)+half)
        x2 = abs(xc(2)-center(2))
        if( ll(2).gt.0 )then
           x2 = min (x2,abs(xc(2)-center(2)+ll(2)))
           x2 = min (x2,abs(xc(2)-center(2)-ll(2)))
        endif
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

             radius2=x0**2+x1**2+x2**2
             var(ic0,ic1,ic2) = temperature_base
     &                        + deltaT*exp(-radius2*factor)

           enddo
        enddo
      enddo
c
      return
      end

c***********************************************************************

      subroutine heat_capacity_nkr(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   conc, ngc,
     &   temp, ngt,
     &   cp, ncp,
     &   cp_powers, npowers,
     &   cpcoeffs, ncoeffs)

      implicit none
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &        ngc, ngt, ncp, npowers, ncoeffs
      integer cp_powers(npowers)
      REAL 
     &     conc(CELL3d(ifirst,ilast,ngc)),
     &     temp(CELL3d(ifirst,ilast,ngt)),
     &     cp(CELL3d(ifirst,ilast,ncp)),
     &     cpcoeffs(ncoeffs)

c     local variables
      integer i, j, k, p
      REAL c, t, cpvals(2)

c      print*,'heat_capacity_nkr, npowers=',npowers,', ncoeffs=',ncoeffs
c      print*,'cpcoeffs(1)=',cpcoeffs(1)
c      print*,'cpcoeffs(2)=',cpcoeffs(2)
c      print*,'cp_powers(1)=',cp_powers(1)

      do k = ifirst2, ilast2
         do j = ifirst1, ilast1
            do i = ifirst0, ilast0
               c=conc(i,j,k)
               t=temp(i,j,k)
            
               cpvals(1)=0.
               do p=1, npowers
                 cpvals(1)=cpvals(1)+cpcoeffs(p)*T**cp_powers(p)
               enddo
               cpvals(2)=0.
               do p=1, npowers
                 cpvals(2)=cpvals(2)+cpcoeffs(npowers+p)*T**cp_powers(p)
               enddo
               cp(i,j,k) = ( cpvals(1)*c
     &                 +cpvals(2)*(1.d0-c) )
c               if(cp(i,j,k).lt.1.E-8)then
c                  print*,'cp=',cp(i,j,k),', c=',c,', T=',t
c                  print*,'cpvals(1)=',cpvals(1),', cpvals(2)=',cpvals(2)
c                  stop
c               endif
            enddo
         enddo
      enddo
      
      return
      end

c***********************************************************************
      
      subroutine initgaussiansource(dx,xlo,xhi,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
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
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ngcp, ngrhs
      REAL standard_dev, heat_peak
      REAL radius2,factor
      REAL center(0:NDIM-1),ll(0:NDIM-1)
      REAL dx(0:NDIM-1),xlo(0:NDIM-1),xhi(0:NDIM-1)
c
c variables in 2d cell indexed         
      REAL cp(CELL3d(ifirst,ilast,ngcp))
      REAL rhs(CELL3d(ifirst,ilast,ngrhs))
c
c***********************************************************************     
c
      integer ic0,ic1,ic2
      REAL xc(0:NDIM-1),x0,x1,x2
      
      factor = 1./(2.*standard_dev*standard_dev)

      do ic2=ifirst2,ilast2
        xc(2) = xlo(2)+dx(2)*(dble(ic2-ifirst2)+half)
        x2 = abs(xc(2)-center(2))
        if( ll(2).gt.0 )then
           x2 = min (x2,abs(xc(2)-center(2)+ll(2)))
           x2 = min (x2,abs(xc(2)-center(2)-ll(2)))
        endif
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

            radius2=x0**2+x1**2+x2**2
            rhs(ic0,ic1,ic2) = 
     &        heat_peak*exp(-radius2*factor)/cp(ic0,ic1,ic2)
          enddo
        enddo
      enddo
c
      return
      end

c***********************************************************************
      
      subroutine linearmeltingline(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  temp, ngt,
     &  conc, ngc,
     &  tref, cref, slope)
c***********************************************************************
      implicit none
c***********************************************************************     
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ngt, ngc
      REAL tref,cref,slope
c
c variables in 3d cell indexed         
      REAL conc(CELL3d(ifirst,ilast,ngc))
      REAL temp(CELL3d(ifirst,ilast,ngt))
c
c***********************************************************************     
c
      integer i,j,k
      
      do k = ifirst2, ilast2
         do j = ifirst1, ilast1
            do i = ifirst0, ilast0
            
               temp(i,j,k)=tref+slope*(conc(i,j,k)-cref)

            enddo
         enddo
      enddo
c
      return
      end

c***********************************************************************

      subroutine initgradient(dx,xlo,xhi,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  gcw0,gcw1,gcw2,
     &  var,
     &  center,temperature_center,
     &  gradient)
c***********************************************************************
      implicit none
c***********************************************************************
      REAL half
      parameter (half=.5d0)

      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer gcw0,gcw1,gcw2
      REAL center(0:NDIM-1)
      REAL temperature_center
      REAL dx(0:NDIM-1),xlo(0:NDIM-1),xhi(0:NDIM-1)
      REAL gradient(0:NDIM-1)
c variables in 3d cell indexed
      REAL var(CELL3dVECG(ifirst,ilast,gcw))

c***********************************************************************
      integer ic0,ic1,ic2
      REAL xc(0:NDIM-1),x0,x1,x2

      do ic2=ifirst2,ilast2
         xc(2) = xlo(2)+dx(2)*(dble(ic2-ifirst2)+half)
         x2 = xc(2)-center(2)
         do ic1=ifirst1,ilast1
            xc(1) = xlo(1)+dx(1)*(dble(ic1-ifirst1)+half)
            x1 = xc(1)-center(1)
            do ic0=ifirst0,ilast0
               xc(0) = xlo(0)+dx(0)*(dble(ic0-ifirst0)+half)
               x0 = xc(0)-center(0)

               var(ic0,ic1,ic2) = temperature_center
     &                          + x0*gradient(0)
     &                          + x1*gradient(1)
     &                          + x2*gradient(2)

            enddo
         enddo
      enddo

      return
      end

