c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
c 
      subroutine set_j_ij3d(
     &     lo0, hi0, lo1, hi1, lo2, hi2,
     &     fcx, fcxlo0, fcxhi0, fcxlo1, fcxhi1, fcxlo2, fcxhi2,
     &     fcy, fcylo0, fcyhi0, fcylo1, fcyhi1, fcylo2, fcyhi2,
     &     fcz, fczlo0, fczhi0, fczlo1, fczhi1, fczlo2, fczhi2,
     &     h,
     &     diag, dlo0, dhi0, dlo1, dhi1, dlo2, dhi2,
     &     offdiagx, odxlo0, odxhi0, odxlo1, odxhi1, odxlo2, odxhi2,
     &     offdiagy, odylo0, odyhi0, odylo1, odyhi1, odylo2, odyhi2,
     &     offdiagz, odzlo0, odzhi0, odzlo1, odzhi1, odzlo2, odzhi2
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1, lo2, hi2,
     &        fcxlo0, fcxhi0, fcxlo1, fcxhi1, fcxlo2, fcxhi2,
     &        fcylo0, fcyhi0, fcylo1, fcyhi1, fcylo2, fcyhi2,
     &        fczlo0, fczhi0, fczlo1, fczhi1, fczlo2, fczhi2,
     &        dlo0, dhi0, dlo1, dhi1, dlo2, dhi2,
     &        odxlo0, odxhi0, odxlo1, odxhi1, odxlo2, odxhi2,
     &        odylo0, odyhi0, odylo1, odyhi1, odylo2, odyhi2,
     &        odzlo0, odzhi0, odzlo1, odzhi1, odzlo2, odzhi2
      double precision 
     &        fcx(fcxlo0:fcxhi0,fcxlo1:fcxhi1,fcxlo2:fcxhi2),
     &        fcy(fcylo0:fcyhi0,fcylo1:fcyhi1,fcylo2:fcyhi2),
     &        fcz(fczlo0:fczhi0,fczlo1:fczhi1,fczlo2:fczhi2),
     &        diag(dlo0:dhi0,dlo1:dhi1,dlo2:dhi2),
     &        offdiagx(odxlo0:odxhi0,odxlo1:odxhi1,odxlo2:odxhi2),
     &        offdiagy(odylo0:odyhi0,odylo1:odyhi1,odylo2:odyhi2),
     &        offdiagz(odzlo0:odzhi0,odzlo1:odzhi1,odzlo2:odzhi2),
     &        h(3)

c     local variables:
      double precision fac
      integer i, j, k

c     Handle x faces
      fac = 1.d0 / h(1)**2
      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0+1
               offdiagx(i,j,k) = fcx(i,j,k) * fac
            enddo
         enddo
      enddo

c     Handle y faces
      fac = 1.d0 / h(2)**2
      do k = lo2, hi2
         do j = lo1, hi1+1
            do i = lo0, hi0
               offdiagy(i,j,k) = fcy(i,j,k) * fac
            enddo
         enddo
      enddo

c     Handle z faces
      fac = 1.d0 / h(3)**2
      do k = lo2, hi2+1
         do j = lo1, hi1
            do i = lo0, hi0
               offdiagz(i,j,k) = fcz(i,j,k) * fac
            enddo
         enddo
      enddo

c     Set the diagonal contributions from the cell boundaries
      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0
               diag(i,j,k) = - offdiagx(i,j,k) - offdiagx(i+1,j,k)
     &                       - offdiagy(i,j,k) - offdiagy(i,j+1,k)
     &                       - offdiagz(i,j,k) - offdiagz(i,j,k+1)
            enddo
         enddo
      enddo

      return
      end

      subroutine set_stencil3d(
     &     lo0, hi0, lo1, hi1, lo2, hi2,
     &     gamma,
     &     sqrt_m, mlo0, mhi0, mlo1, mhi1, mlo2, mhi2,
     &     diag, dlo0, dhi0, dlo1, dhi1, dlo2, dhi2,
     &     offdiagx, odxlo0, odxhi0, odxlo1, odxhi1, odxlo2, odxhi2,
     &     offdiagy, odylo0, odyhi0, odylo1, odyhi1, odylo2, odyhi2,
     &     offdiagz, odzlo0, odzhi0, odzlo1, odzhi1, odzlo2, odzhi2,
     &     values
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1, lo2, hi2,
     &        mlo0, mhi0, mlo1, mhi1, mlo2, mhi2,
     &        dlo0, dhi0, dlo1, dhi1, dlo2, dhi2,
     &        odxlo0, odxhi0, odxlo1, odxhi1, odxlo2, odxhi2,
     &        odylo0, odyhi0, odylo1, odyhi1, odylo2, odyhi2,
     &        odzlo0, odzhi0, odzlo1, odzhi1, odzlo2, odzhi2
      double precision gamma,
     &           sqrt_m(mlo0:mhi0,mlo1:mhi1,mlo2:mhi2),
     &           diag(dlo0:dhi0,dlo1:dhi1,dlo2:dhi2),
     &           offdiagx(odxlo0:odxhi0,odxlo1:odxhi1,odxlo2:odxhi2),
     &           offdiagy(odylo0:odyhi0,odylo1:odyhi1,odylo2:odyhi2),
     &           offdiagz(odzlo0:odzhi0,odzlo1:odzhi1,odzlo2:odzhi2),
     &           values(1)

c     local variables:
      double precision fac
      integer i, j, k, index

      index = -6
      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0

               fac = gamma * sqrt_m(i,j,k)

               index = index + 7

c              center
               values(index)   = fac*diag(i,j,k)*sqrt_m(i,j,k) + 1.d0

c              west
               values(index+1) = fac*offdiagx(i,j,k)*sqrt_m(i-1,j,k)

c              east
               values(index+2) = fac*offdiagx(i+1,j,k)*sqrt_m(i+1,j,k)

c              south
               values(index+3) = fac*offdiagy(i,j,k)*sqrt_m(i,j-1,k)

c              north
               values(index+4) = fac*offdiagy(i,j+1,k)*sqrt_m(i,j+1,k)

c              front
               values(index+5) = fac*offdiagz(i,j,k)*sqrt_m(i,j,k-1)

c              back
               values(index+6) = fac*offdiagz(i,j,k+1)*sqrt_m(i,j,k+1)

            enddo
         enddo
      enddo
      
      return
      end

      subroutine set_symmetric_stencil3d(
     &     lo0, hi0, lo1, hi1, lo2, hi2,
     &     gamma,
     &     sqrt_m, mlo0, mhi0, mlo1, mhi1, mlo2, mhi2,
     &     diag, dlo0, dhi0, dlo1, dhi1, dlo2, dhi2,
     &     offdiagx, odxlo0, odxhi0, odxlo1, odxhi1, odxlo2, odxhi2,
     &     offdiagy, odylo0, odyhi0, odylo1, odyhi1, odylo2, odyhi2,
     &     offdiagz, odzlo0, odzhi0, odzlo1, odzhi1, odzlo2, odzhi2,
     &     values
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1, lo2, hi2,
     &        mlo0, mhi0, mlo1, mhi1, mlo2, mhi2,
     &        dlo0, dhi0, dlo1, dhi1, dlo2, dhi2,
     &        odxlo0, odxhi0, odxlo1, odxhi1, odxlo2, odxhi2,
     &        odylo0, odyhi0, odylo1, odyhi1, odylo2, odyhi2,
     &        odzlo0, odzhi0, odzlo1, odzhi1, odzlo2, odzhi2
      double precision gamma,
     &           sqrt_m(mlo0:mhi0,mlo1:mhi1,mlo2:mhi2),
     &           diag(dlo0:dhi0,dlo1:dhi1,dlo2:dhi2),
     &           offdiagx(odxlo0:odxhi0,odxlo1:odxhi1,odxlo2:odxhi2),
     &           offdiagy(odylo0:odyhi0,odylo1:odyhi1,odylo2:odyhi2),
     &           offdiagz(odzlo0:odzhi0,odzlo1:odzhi1,odzlo2:odzhi2),
     &           values(1)

c     local variables:
      double precision fac
      integer i, j, k, index

      index = -3
      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0

               fac = gamma * sqrt_m(i,j,k)

               index = index + 4

c              center
               values(index)   = fac*diag(i,j,k)*sqrt_m(i,j,k) + 1.d0

c              west
               values(index+1) = fac*offdiagx(i,j,k)*sqrt_m(i-1,j,k)

c              south
               values(index+2) = fac*offdiagy(i,j,k)*sqrt_m(i,j-1,k)

c              front
               values(index+3) = fac*offdiagz(i,j,k)*sqrt_m(i,j,k-1)

            enddo
         enddo
      enddo
      
      return
      end

      subroutine copy3d(
     &     lo0, hi0, lo1, hi1, lo2, hi2,
     &     src, slo0, shi0, slo1, shi1, slo2, shi2,
     &     dst, dlo0, dhi0, dlo1, dhi1, dlo2, dhi2
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1, lo2, hi2,
     &        slo0, shi0, slo1, shi1, slo2, shi2,
     &        dlo0, dhi0, dlo1, dhi1, dlo2, dhi2
      double precision src(slo0:shi0,slo1:shi1,slo2:shi2),
     &                 dst(dlo0:dhi0,dlo1:dhi1,dlo2:dhi2)

c     local variables:
      integer i, j, k

      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0
               dst(i,j,k) = src(i,j,k)
            enddo
         enddo
      enddo

      return
      end

      subroutine adjust_bdry3d(
     &     gamma,
     &     sqrt_m, mlo0, mhi0, mlo1, mhi1, mlo2, mhi2,
     &     diag, offdiagi, offdiagj, offdiagk,
     &     pifirst, pilast, pjfirst, pjlast, pkfirst, pklast,
     &     acoef,
     &     bcoef,
     &     aifirst, ailast, ajfirst, ajlast, akfirst, aklast,
     &     auk0,
     &     kifirst, kilast, kjfirst, kjlast, kkfirst, kklast,
     &     lower, upper,
     &     location, h )

      implicit none
      double precision gamma
      integer mlo0, mhi0, mlo1, mhi1, mlo2, mhi2
      double precision sqrt_m(mlo0:mhi0,mlo1:mhi1,mlo2:mhi2)
      integer pifirst, pilast, pjfirst, pjlast, pkfirst, pklast
      double precision 
     &     diag(pifirst:pilast,pjfirst:pjlast,pkfirst:pklast)
      double precision
     &     offdiagi(pifirst:pilast+1,pjfirst:pjlast,pkfirst:pklast),
     &     offdiagj(pifirst:pilast,pjfirst:pjlast+1,pkfirst:pklast),
     &     offdiagk(pifirst:pilast,pjfirst:pjlast,pkfirst:pklast+1)
      integer aifirst, ailast, ajfirst, ajlast, akfirst, aklast
      double precision 
     &     acoef(aifirst:ailast,ajfirst:ajlast,akfirst:aklast)
      double precision 
     &     bcoef(aifirst:ailast,ajfirst:ajlast,akfirst:aklast)
      integer kifirst, kilast, kjfirst, kjlast, kkfirst, kklast
      double precision
     &     auk0(kifirst:kilast,kjfirst:kjlast,kkfirst:kklast)
      integer lower(0:2), upper(0:2)
      integer location
      double precision h(0:2), hh
      integer igho, iint, ifac
      integer jgho, jint, jfac
      integer kgho, kint, kfac
      double precision uk0, k1
c     Nomenclature for indices: gho=ghost, int=interior,
c     fac=surface, beg=beginning, end=ending.
      integer i, j, k
      hh = h(location/2)
      if ( location .eq. 0 ) then
c     min i side
         igho = upper(0)
         ifac = igho + 1
         iint = igho + 1
         do k=lower(2),upper(2)
            do j=lower(1),upper(1)
               uk0 = (hh)
     &            / (1-acoef(ifac,j,k)*(1-0.5*hh))
               k1 = (1-acoef(ifac,j,k)*(1+0.5*hh))
     &            / (1-acoef(ifac,j,k)*(1-0.5*hh))
               diag(iint,j,k) = diag(iint,j,k)
     &                              + k1*offdiagi(ifac,j,k)
               auk0(ifac,j,k) 
     &              = uk0*offdiagi(ifac,j,k)*sqrt_m(igho,j,k)
     &                *gamma*sqrt_m(iint,j,k)
               offdiagi(ifac,j,k) = 0.0
            enddo
         enddo
      elseif ( location .eq. 1 ) then
c     max i side
         igho = lower(0)
         ifac = igho
         iint = igho - 1
         do k=lower(2),upper(2)
            do j=lower(1),upper(1)
               uk0 = (hh)
     &            / (1-acoef(ifac,j,k)*(1-0.5*hh))
               k1 = (1-acoef(ifac,j,k)*(1+0.5*hh))
     &            / (1-acoef(ifac,j,k)*(1-0.5*hh))
               diag(iint,j,k) = diag(iint,j,k)
     &                              + k1*offdiagi(ifac,j,k)
               auk0(ifac,j,k) 
     &              = uk0*offdiagi(ifac,j,k)*sqrt_m(igho,j,k)
     &                *gamma*sqrt_m(iint,j,k)
               offdiagi(ifac,j,k) = 0.0
            enddo
         enddo
      elseif ( location .eq. 2 ) then
c     min j side
         jgho = upper(1)
         jfac = jgho + 1
         jint = jgho + 1
         do k=lower(2),upper(2)
            do i=lower(0),upper(0)
               uk0 = (hh)
     &            / (1-acoef(i,jfac,k)*(1-0.5*hh))
               k1 = (1-acoef(i,jfac,k)*(1+0.5*hh))
     &            / (1-acoef(i,jfac,k)*(1-0.5*hh))
               diag(i,jint,k) = diag(i,jint,k)
     &                              + k1*offdiagj(i,jfac,k)
               auk0(i,jfac,k) 
     &              = uk0*offdiagj(i,jfac,k)*sqrt_m(i,jgho,k)
     &                *gamma*sqrt_m(i,jint,k)
               offdiagj(i,jfac,k) = 0.0
            enddo
         enddo
      elseif ( location .eq. 3 ) then
c     max j side
         jgho = lower(1)
         jfac = jgho
         jint = jgho - 1
         do k=lower(2),upper(2)
            do i=lower(0),upper(0)
               uk0 = (hh)
     &            / (1-acoef(i,jfac,k)*(1-0.5*hh))
               k1 = (1-acoef(i,jfac,k)*(1+0.5*hh))
     &            / (1-acoef(i,jfac,k)*(1-0.5*hh))
               diag(i,jint,k) = diag(i,jint,k)
     &                              + k1*offdiagj(i,jfac,k)
               auk0(i,jfac,k) 
     &              = uk0*offdiagj(i,jfac,k)*sqrt_m(i,jgho,k)
     &                *gamma*sqrt_m(i,jint,k)
               offdiagj(i,jfac,k) = 0.0
            enddo
         enddo
      elseif ( location .eq. 4 ) then
c     min k side
         kgho = upper(2)
         kfac = kgho + 1
         kint = kgho + 1
         do j=lower(1),upper(1)
            do i=lower(0),upper(0)
               uk0 = (hh)
     &            / (1-acoef(i,j,kfac)*(1-0.5*hh))
               k1 = (1-acoef(i,j,kfac)*(1+0.5*hh))
     &            / (1-acoef(i,j,kfac)*(1-0.5*hh))
               diag(i,j,kint) = diag(i,j,kint)
     &                              + k1*offdiagk(i,j,kfac)
               auk0(i,j,kfac) 
     &              = uk0*offdiagk(i,j,kfac)*sqrt_m(i,j,kgho)
     &                *gamma*sqrt_m(i,j,kint)
               offdiagk(i,j,kfac) = 0.0
            enddo
         enddo
      elseif ( location .eq. 5 ) then
c     max k side
         kgho = lower(2)
         kfac = kgho
         kint = kgho - 1
         do j=lower(1),upper(1)
            do i=lower(0),upper(0)
               uk0 = (hh)
     &            / (1-acoef(i,j,kfac)*(1-0.5*hh))
               k1 = (1-acoef(i,j,kfac)*(1+0.5*hh))
     &            / (1-acoef(i,j,kfac)*(1-0.5*hh))
               diag(i,j,kint) = diag(i,j,kint)
     &                              + k1*offdiagk(i,j,kfac)
               auk0(i,j,kfac) 
     &              = uk0*offdiagk(i,j,kfac)*sqrt_m(i,j,kgho)
     &                *gamma*sqrt_m(i,j,kint)
               offdiagk(i,j,kfac) = 0.0
            enddo
         enddo
      endif
      return
      end

      subroutine adjust_qrhs3d(
     &     rhs, rifirst, rilast, rjfirst, rjlast, rkfirst, rklast,
     &     auk0, kifirst, kilast, kjfirst, kjlast, kkfirst, kklast,
     &     gcoef, aifirst, ailast, ajfirst, ajlast, akfirst, aklast,
     &     lower, upper,
     &     location)

      implicit none
      integer rifirst, rilast, rjfirst, rjlast, rkfirst, rklast,
     &        kifirst, kilast, kjfirst, kjlast, kkfirst, kklast,
     &     aifirst, ailast, ajfirst, ajlast, akfirst, aklast,
     &     lower(0:2), upper(0:2), location
      double precision
     &     rhs(rifirst:rilast,rjfirst:rjlast,rkfirst:rklast),
     &     auk0(kifirst:kilast,kjfirst:kjlast,kkfirst:kklast),
     &     gcoef(aifirst:ailast,ajfirst:ajlast,akfirst:aklast)

c     local data:
      integer igho, jgho, kgho, icel, jcel, kcel, ifac, jfac, kfac,
     &     i, j, k

c     Nomenclature for indices: cel=first-cell, gho=ghost,
c     beg=beginning, end=ending.

      if ( location .eq. 0 ) then
c        min i side
         igho = upper(0)
         ifac = igho + 1
         icel = igho + 1
         do j = lower(1), upper(1)
            do k = lower(2), upper(2)
               rhs(icel,j,k) = rhs(icel,j,k) 
     &              - auk0(ifac,j,k)*gcoef(ifac,j,k)
            enddo
         enddo
      elseif ( location .eq. 1 ) then
c        max i side
         igho = lower(0)
         ifac = igho
         icel = igho - 1
         do j = lower(1), upper(1)
            do k = lower(2), upper(2)
               rhs(icel,j,k) = rhs(icel,j,k)
     &              - auk0(ifac,j,k)*gcoef(ifac,j,k)
            enddo
         enddo
      elseif ( location .eq. 2 ) then
c        min j side
         jgho = upper(1)
         jfac = jgho + 1
         jcel = jgho + 1
         do k = lower(2), upper(2)
            do i = lower(0), upper(0)
               rhs(i,jcel,k) = rhs(i,jcel,k)
     &              - auk0(i,jfac,k)*gcoef(i,jfac,k)
            enddo
         enddo
      elseif ( location .eq. 3 ) then
c        max j side
         jgho = lower(1)
         jfac = jgho
         jcel = jgho - 1
         do k = lower(2), upper(2)
            do i = lower(0), upper(0)
               rhs(i,jcel,k) = rhs(i,jcel,k)
     &              - auk0(i,jfac,k)*gcoef(i,jfac,k)
            enddo
         enddo
      elseif ( location .eq. 4 ) then
c     min k side
         kgho = upper(2)
         kfac = kgho + 1
         kcel = kgho + 1
         do j = lower(1), upper(1)
            do i = lower(0), upper(0)
               rhs(i,j,kcel) = rhs(i,j,kcel)
     &              - auk0(i,j,kfac)*gcoef(i,j,kfac)
            enddo
         enddo
      elseif ( location .eq. 5 ) then
c     max k side
         kgho = lower(2)
         kfac = kgho
         kcel = kgho - 1
         do j = lower(1), upper(1)
            do i = lower(0), upper(0)
               rhs(i,j,kcel) = rhs(i,j,kcel)
     &              - auk0(i,j,kfac)*gcoef(i,j,kfac)
            enddo
         enddo
      endif

      return
      end
