c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
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
      
      subroutine initgaussian(dx,xlo,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  var, ng,
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
      integer ng
      REAL standard_dev,temperature_base, temperature_peak
      REAL radius2,factor,deltaT
      REAL center(0:NDIM-1),ll(0:NDIM-1)
      REAL dx(0:NDIM-1),xlo(0:NDIM-1)
c
c variables in 2d cell indexed
      REAL var(CELL2d(ifirst,ilast,ng))
c
c***********************************************************************     
c
      integer ic0,ic1
      REAL xc(0:NDIM-1),x0,x1
      
      deltaT = temperature_peak-temperature_base
      
      factor = 1./(2.*standard_dev*standard_dev)

      do ic1=ifirst1-ng,ilast1+ng
        xc(1) = xlo(1)+dx(1)*(dble(ic1-ifirst1)+half)
        x1 = abs(xc(1)-center(1))
        if( ll(1).gt.0 )then
           x1 = min (x1,abs(xc(1)-center(1)+ll(1)))
           x1 = min (x1,abs(xc(1)-center(1)-ll(1)))
        endif
        do ic0=ifirst0-ng,ilast0+ng
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
      
      subroutine initgaussiansource(dx,xlo,
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
      REAL dx(0:NDIM-1),xlo(0:NDIM-1)
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
c Input:
c   ref_position: position of global node (0,0)
c   ref_temperature: temperature at ref_position
c
      subroutine initgradient(dx,xlo,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  var, ng,
     &  ref_position,ref_temperature,
     &  gradient)
c***********************************************************************
      implicit none
c***********************************************************************
      REAL half
      parameter (half=.5d0)

      integer ifirst0,ilast0,ifirst1,ilast1, ng
      REAL ref_position(0:NDIM-1)
      REAL ref_temperature
      REAL dx(0:NDIM-1),xlo(0:NDIM-1)
      REAL gradient(0:NDIM-1)
c variables in 2d cell indexed
      REAL var(CELL2d(ifirst,ilast,ng))

c***********************************************************************
      integer ic0,ic1
      REAL xc(0:NDIM-1),x0,x1

      do ic1=ifirst1-ng,ilast1+ng
        xc(1) = xlo(1)+dx(1)*(dble(ic1-ifirst1)+half)
        x1 = xc(1)-ref_position(1)
        do ic0=ifirst0-ng,ilast0+ng
           xc(0) = xlo(0)+dx(0)*(dble(ic0-ifirst0)+half)
           x0 = xc(0)-ref_position(0)

           var(ic0,ic1) = ref_temperature
     &                  + x0*gradient(0)
     &                  + x1*gradient(1)

         enddo
      enddo

      return
      end
