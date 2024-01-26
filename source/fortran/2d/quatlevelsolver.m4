c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
c 
      subroutine set_j_ij2d(
     &     lo0, hi0, lo1, hi1,
     &     fcx, fcxlo0, fcxhi0, fcxlo1, fcxhi1,
     &     fcy, fcylo0, fcyhi0, fcylo1, fcyhi1,
     &     h,
     &     diag, dlo0, dhi0, dlo1, dhi1,
     &     offdiagx, odxlo0, odxhi0, odxlo1, odxhi1,
     &     offdiagy, odylo0, odyhi0, odylo1, odyhi1
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1,
     &        fcxlo0, fcxhi0, fcxlo1, fcxhi1,
     &        fcylo0, fcyhi0, fcylo1, fcyhi1,
     &        dlo0, dhi0, dlo1, dhi1,
     &        odxlo0, odxhi0, odxlo1, odxhi1,
     &        odylo0, odyhi0, odylo1, odyhi1
      double precision 
     &        fcx(fcxlo0:fcxhi0,fcxlo1:fcxhi1),
     &        fcy(fcylo0:fcyhi0,fcylo1:fcyhi1),
     &        diag(dlo0:dhi0,dlo1:dhi1),
     &        offdiagx(odxlo0:odxhi0,odxlo1:odxhi1),
     &        offdiagy(odylo0:odyhi0,odylo1:odyhi1),
     &        h(2)

c     local variables:
      double precision fac
      integer i, j

c     Handle x faces
      fac = 1.d0 / h(1)**2
      do j = lo1, hi1
         do i = lo0, hi0+1
            offdiagx(i,j) = fcx(i,j) * fac
         enddo
      enddo

c     Handle y faces
      fac = 1.d0 / h(2)**2
      do j = lo1, hi1+1
         do i = lo0, hi0
            offdiagy(i,j) = fcy(i,j) * fac
         enddo
      enddo

c     Set the diagonal contributions from the cell boundaries
      do j = lo1, hi1
         do i = lo0, hi0
            diag(i,j) = - offdiagx(i,j) - offdiagx(i+1,j)
     &                  - offdiagy(i,j) - offdiagy(i,j+1)
         enddo
      enddo

      return
      end

      subroutine set_stencil2d(
     &     lo0, hi0, lo1, hi1,
     &     gamma,
     &     sqrt_m, mlo0, mhi0, mlo1, mhi1,
     &     diag, dlo0, dhi0, dlo1, dhi1,
     &     offdiagx, odxlo0, odxhi0, odxlo1, odxhi1,
     &     offdiagy, odylo0, odyhi0, odylo1, odyhi1,
     &     values
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1,
     &        mlo0, mhi0, mlo1, mhi1,
     &        dlo0, dhi0, dlo1, dhi1,
     &        odxlo0, odxhi0, odxlo1, odxhi1,
     &        odylo0, odyhi0, odylo1, odyhi1
      double precision gamma,
     &        sqrt_m(mlo0:mhi0,mlo1:mhi1),
     &        diag(dlo0:dhi0,dlo1:dhi1),
     &        offdiagx(odxlo0:odxhi0,odxlo1:odxhi1),
     &        offdiagy(odylo0:odyhi0,odylo1:odyhi1),
     &        values(1)

c     local variables:
      double precision fac
      integer i, j, index

      index = -4
      do j = lo1, hi1
         do i = lo0, hi0
                  
            fac = gamma * sqrt_m(i,j)

            index = index + 5

c           center
            values(index)   = fac * diag(i,j) * sqrt_m(i,j) + 1.d0

c           west
            values(index+1) = fac * offdiagx(i,j) * sqrt_m(i-1,j)

c           east
            values(index+2) = fac * offdiagx(i+1,j) * sqrt_m(i+1,j)

c           south
            values(index+3) = fac * offdiagy(i,j) * sqrt_m(i,j-1)

c           north
            values(index+4) = fac * offdiagy(i,j+1) * sqrt_m(i,j+1)

         enddo
      enddo
      
      return
      end

      subroutine set_symmetric_stencil2d(
     &     lo0, hi0, lo1, hi1,
     &     gamma,
     &     sqrt_m, mlo0, mhi0, mlo1, mhi1,
     &     diag, dlo0, dhi0, dlo1, dhi1,
     &     offdiagx, odxlo0, odxhi0, odxlo1, odxhi1,
     &     offdiagy, odylo0, odyhi0, odylo1, odyhi1,
     &     values
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1,
     &        mlo0, mhi0, mlo1, mhi1,
     &        dlo0, dhi0, dlo1, dhi1,
     &        odxlo0, odxhi0, odxlo1, odxhi1,
     &        odylo0, odyhi0, odylo1, odyhi1
      double precision gamma,
     &        sqrt_m(mlo0:mhi0,mlo1:mhi1),
     &        diag(dlo0:dhi0,dlo1:dhi1),
     &        offdiagx(odxlo0:odxhi0,odxlo1:odxhi1),
     &        offdiagy(odylo0:odyhi0,odylo1:odyhi1),
     &        values(1)

c     local variables:
      double precision fac
      integer i, j, index

      index = -2
      do j = lo1, hi1
         do i = lo0, hi0
                  
            fac = gamma * sqrt_m(i,j)

            index = index + 3

c           center
            values(index)   = fac * diag(i,j) * sqrt_m(i,j) + 1.d0

c           west
            values(index+1) = fac * offdiagx(i,j) * sqrt_m(i-1,j)

c           south
            values(index+2) = fac * offdiagy(i,j) * sqrt_m(i,j-1)

         enddo
      enddo
      
      return
      end

      subroutine copy2d(
     &     lo0, hi0, lo1, hi1,
     &     src, slo0, shi0, slo1, shi1,
     &     dst, dlo0, dhi0, dlo1, dhi1
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1,
     &        slo0, shi0, slo1, shi1,
     &        dlo0, dhi0, dlo1, dhi1
      double precision src(slo0:shi0,slo1:shi1),
     &                 dst(dlo0:dhi0,dlo1:dhi1)

c     local variables:
      integer i, j

      do j = lo1, hi1
         do i = lo0, hi0
            dst(i,j) = src(i,j)
         enddo
      enddo

      return
      end

      subroutine adjust_bdry2d(
     &     gamma,
     &     sqrt_m, mlo0, mhi0, mlo1, mhi1,
     &     diag, offdiagi, offdiagj,
     &     pifirst, pilast, pjfirst, pjlast,
     &     acoef, bcoef,
     &     aifirst, ailast, ajfirst, ajlast,
     &     auk0,
     &     kifirst, kilast, kjfirst, kjlast,
     &     lower, upper,
     &     location, h )

      implicit none
      double precision gamma
      integer mlo0, mhi0, mlo1, mhi1
      double precision sqrt_m(mlo0:mhi0,mlo1:mhi1)
      integer pifirst, pilast, pjfirst, pjlast
      double precision diag(pifirst:pilast,pjfirst:pjlast)
      double precision offdiagi(pifirst:pilast+1,pjfirst:pjlast)
      double precision offdiagj(pifirst:pilast,pjfirst:pjlast+1)
      integer aifirst, ailast, ajfirst, ajlast
      double precision acoef(aifirst:ailast,ajfirst:ajlast)
      double precision bcoef(aifirst:ailast,ajfirst:ajlast)
      integer kifirst, kilast, kjfirst, kjlast
      double precision auk0(kifirst:kilast,kjfirst:kjlast)
      integer lower(0:1), upper(0:1)
      integer location
      double precision h(0:1), hh
      integer igho, iint, ifac
      integer jgho, jint, jfac
      double precision uk0, k1

c     Nomenclature for indices: gho=ghost, int=interior, fac=surface

      integer i, j

      hh = h(location/2)

      if ( location .eq. 0 ) then
c        low i edge
         igho = upper(0)
         iint = igho + 1
         ifac = igho + 1
         do j=lower(1),upper(1)
            uk0 = (hh)
     &         / (1-acoef(ifac,j)*(1-0.5*hh))
            k1 = (1-acoef(ifac,j)*(1+0.5*hh))
     &         / (1-acoef(ifac,j)*(1-0.5*hh))
            diag(iint,j) = diag(iint,j) 
     &           + k1*offdiagi(ifac,j)
            auk0(ifac,j) = uk0*offdiagi(ifac,j)*sqrt_m(igho,j)
     &                     *gamma*sqrt_m(iint,j)
            offdiagi(ifac,j) = 0.0
         enddo
      elseif ( location .eq. 1 ) then
c        high i edge
         igho = lower(0)
         iint = igho - 1
         ifac = igho
         do j=lower(1),upper(1)
            uk0 = (hh)
     &         / (1-acoef(ifac,j)*(1-0.5*hh))
            k1 = (1-acoef(ifac,j)*(1+0.5*hh))
     &         / (1-acoef(ifac,j)*(1-0.5*hh))
            diag(iint,j) = diag(iint,j)
     &           + k1*offdiagi(ifac,j)
            auk0(ifac,j) = uk0*offdiagi(ifac,j)*sqrt_m(igho,j)
     &                     *gamma*sqrt_m(iint,j)
            offdiagi(ifac,j) = 0.0
         enddo
      elseif ( location .eq. 2 ) then
c        low j edge
         jgho = upper(1)
         jint = jgho + 1
         jfac = jgho + 1
         do i=lower(0),upper(0)
            uk0 = (hh)
     &         / (1-acoef(i,jfac)*(1-0.5*hh))
            k1 = (1-acoef(i,jfac)*(1+0.5*hh))
     &         / (1-acoef(i,jfac)*(1-0.5*hh))
            diag(i,jint) = diag(i,jint)
     &           + k1*offdiagj(i,jfac)
            auk0(i,jfac) = uk0*offdiagj(i,jfac)*sqrt_m(i,jgho)
     &                     *gamma*sqrt_m(i,jint)
            offdiagj(i,jfac) = 0.0
         enddo
      elseif ( location .eq. 3 ) then
c        high j edge
         jgho = lower(1)
         jint = jgho - 1
         jfac = jgho
         do i=lower(0),upper(0)
            uk0 = (hh)
     &         / (1-acoef(i,jfac)*(1-0.5*hh))
            k1 = (1-acoef(i,jfac)*(1+0.5*hh))
     &         / (1-acoef(i,jfac)*(1-0.5*hh))
            diag(i,jint) = diag(i,jint)
     &           + k1*offdiagj(i,jfac)
            auk0(i,jfac) = uk0*offdiagj(i,jfac)*sqrt_m(i,jgho)
     &                     *gamma*sqrt_m(i,jint)
            offdiagj(i,jfac) = 0.0
         enddo
      endif

      return
      end

      subroutine adjust_qrhs2d(
     &     rhs, rifirst, rilast, rjfirst, rjlast,
     &     auk0, kifirst, kilast, kjfirst, kjlast,
     &     gcoef, aifirst, ailast, ajfirst, ajlast,
     &     lower, upper,
     &     location )

      implicit none
      integer rifirst, rilast, rjfirst, rjlast,
     &        kifirst, kilast, kjfirst, kjlast,
     &        aifirst, ailast, ajfirst, ajlast,
     &        lower(0:1), upper(0:1), location
      double precision rhs(rifirst:rilast,rjfirst:rjlast),
     &        auk0(kifirst:kilast,kjfirst:kjlast),
     &        gcoef(aifirst:ailast,ajfirst:ajlast)

c     local data:
      integer igho, jgho, iint, jint, ifac, jfac, i, j

c     Nomenclature for indices: cel=first-cell, gho=ghost,
c     beg=beginning, end=ending.

      if ( location .eq. 0 ) then
c        min i edge
         igho = upper(0)
         ifac = igho + 1
         iint = igho + 1
         do j = lower(1), upper(1)
            rhs(iint,j) = rhs(iint,j) - auk0(ifac,j)*gcoef(ifac,j)
         enddo
      elseif ( location .eq. 1 ) then
c        max i edge
         igho = lower(0)
         ifac = igho
         iint = igho - 1
         do j = lower(1), upper(1)
            rhs(iint,j) = rhs(iint,j) - auk0(ifac,j)*gcoef(ifac,j)
            enddo
      elseif ( location .eq. 2 ) then
c        min j edge
         jgho = upper(1)
         jfac = jgho + 1
         jint = jgho + 1
         do i = lower(0), upper(0)
            rhs(i,jint) = rhs(i,jint) - auk0(i,jfac)*gcoef(i,jfac)
         enddo
      elseif ( location .eq. 3 ) then
c        max j edge
         jgho = lower(1)
         jfac = jgho
         jint = jgho - 1
         do i = lower(0), upper(0)
            rhs(i,jint) = rhs(i,jint) - auk0(i,jfac)*gcoef(i,jfac)
         enddo
      endif

      return
      end
