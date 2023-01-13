c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
define(NDIM,3)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c***********************************************************************
c
c compute the concentration flux
c
      subroutine concentrationflux(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   dx,
     &   conc, ngconc,
     &   phi, ngphi,
     &   eta, ngeta,
     &   diffconc0,  diffconc1,  diffconc2, ngdiffconc,
     &   dphicoupl0, dphicoupl1, dphicoupl2, ngdphicoupl,
     &   detacoupl0, detacoupl1, detacoupl2, ngdetacoupl,
     &   flux0, flux1, flux2, ngflux,
     &   three_phase )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2, ngflux
      integer three_phase
      double precision 
     &     flux0(SIDE3d0(ifirst,ilast,ngflux)),
     &     flux1(SIDE3d1(ifirst,ilast,ngflux)),
     &     flux2(SIDE3d2(ifirst,ilast,ngflux))
      double precision dx(0:2)
      integer ngconc, ngphi, ngeta
      integer ngdiffconc, ngdphicoupl, ngdetacoupl
c
c variables in 3d cell indexed
      double precision conc(CELL3d(ifirst,ilast,ngconc))
      double precision phi(CELL3d(ifirst,ilast,ngphi))
      double precision eta(CELL3d(ifirst,ilast,ngphi))
      double precision diffconc0(SIDE3d0(ifirst,ilast,ngdiffconc))
      double precision diffconc1(SIDE3d1(ifirst,ilast,ngdiffconc))
      double precision diffconc2(SIDE3d2(ifirst,ilast,ngdiffconc))

c variables in 3d side indexed
      double precision
     &     dphicoupl0(SIDE3d0(ifirst,ilast,ngdphicoupl))
      double precision
     &     dphicoupl1(SIDE3d1(ifirst,ilast,ngdphicoupl))
      double precision
     &     dphicoupl2(SIDE3d2(ifirst,ilast,ngdphicoupl))
      double precision
     &     detacoupl0(SIDE3d0(ifirst,ilast,ngdetacoupl))
      double precision
     &     detacoupl1(SIDE3d1(ifirst,ilast,ngdetacoupl))
      double precision
     &     detacoupl2(SIDE3d2(ifirst,ilast,ngdetacoupl))
c
      double precision dxinv, dyinv, dzinv
      integer          ic0, ic1, ic2

      dxinv = 1.d0 / dx(0)

      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1
               flux0(ic0,ic1,ic2) = dxinv * (
     &              diffconc0(ic0,ic1,ic2) *
     &              ( conc(ic0,ic1,ic2) - conc(ic0-1,ic1,ic2) )
     &            + dphicoupl0(ic0,ic1,ic2) *
     &              ( phi(ic0,ic1,ic2) - phi(ic0-1,ic1,ic2) )
     &              )
            enddo
         enddo
      enddo

      dyinv = 1.d0 / dx(1)

      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0, ilast0
               flux1(ic0,ic1,ic2) = dyinv * (
     &              diffconc1(ic0,ic1,ic2) *
     &              ( conc(ic0,ic1,ic2) - conc(ic0,ic1-1,ic2) )
     &            + dphicoupl1(ic0,ic1,ic2) *
     &              ( phi(ic0,ic1,ic2) - phi(ic0,ic1-1,ic2) )
     &              )
            enddo
         enddo
      enddo

      dzinv = 1.d0 / dx(2)

      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0
               flux2(ic0,ic1,ic2) = dzinv * (
     &              diffconc2(ic0,ic1,ic2) *
     &              ( conc(ic0,ic1,ic2) - conc(ic0,ic1,ic2-1) )
     &            + dphicoupl2(ic0,ic1,ic2) *
     &              ( phi(ic0,ic1,ic2) - phi(ic0,ic1,ic2-1) )
     &              )
            enddo
         enddo
      enddo

      if ( three_phase /= 0 ) then

         do ic2 = ifirst2, ilast2
            do ic1 = ifirst1, ilast1
               do ic0 = ifirst0, ilast0+1
                  flux0(ic0,ic1,ic2) = flux0(ic0,ic1,ic2)
     &                 + dxinv * detacoupl0(ic0,ic1,ic2) *
     &                 ( eta(ic0,ic1,ic2) - eta(ic0-1,ic1,ic2) )
               enddo
            enddo
         enddo

         do ic2 = ifirst2, ilast2
            do ic1 = ifirst1, ilast1+1
               do ic0 = ifirst0, ilast0
                  flux1(ic0,ic1,ic2) = flux1(ic0,ic1,ic2)
     &                 + dyinv * detacoupl1(ic0,ic1,ic2) *
     &                 ( eta(ic0,ic1,ic2) - eta(ic0,ic1-1,ic2) )
               enddo
            enddo
         enddo

         do ic2 = ifirst2, ilast2+1
            do ic1 = ifirst1, ilast1
               do ic0 = ifirst0, ilast0
                  flux2(ic0,ic1,ic2) = flux2(ic0,ic1,ic2)
     &                 + dzinv * detacoupl2(ic0,ic1,ic2) *
     &                 ( eta(ic0,ic1,ic2) - eta(ic0,ic1,ic2-1) )
               enddo
            enddo
         enddo

      endif

      return
      end

c***********************************************************************
c
c compute the concentration flux
c
      subroutine add_concentrationflux_ebs(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   dx,
     &   conc, ngconc,
     &   ncomp,
     &   diffconc0,  diffconc1, diffconc2,  ngdiff,
     &   flux0, flux1, flux2, ngflux )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngflux
      double precision dx(0:2)
      integer ncomp
      integer ngconc
      integer ngdiff
c
c variables in 3d cell indexed
      double precision 
     &     flux0(SIDE3d0(ifirst,ilast,ngflux),ncomp),
     &     flux1(SIDE3d1(ifirst,ilast,ngflux),ncomp),
     &     flux2(SIDE3d2(ifirst,ilast,ngflux),ncomp)
      double precision conc(CELL3d(ifirst,ilast,ngconc),ncomp)
      double precision diffconc0(SIDE3d0(ifirst,ilast,ngdiff),
     &                           ncomp*ncomp)
      double precision diffconc1(SIDE3d1(ifirst,ilast,ngdiff),
     &                           ncomp*ncomp)
      double precision diffconc2(SIDE3d2(ifirst,ilast,ngdiff),
     &                           ncomp*ncomp)
c
      double precision dxinv, dyinv, dzinv
      integer          ic0, ic1, ic2, ic, jc, ijc

      dxinv = 1.d0 / dx(0)
      dyinv = 1.d0 / dx(1)
      dzinv = 1.d0 / dx(2)

      do ic = 1, ncomp
         do jc = 1, ncomp
            ijc=ic+(jc-1)*ncomp
            do ic2 = ifirst2, ilast2
               do ic1 = ifirst1, ilast1
                  do ic0 = ifirst0, ilast0+1
                     flux0(ic0,ic1,ic2,ic) = flux0(ic0,ic1,ic2,ic)
     &                  +dxinv * (
     &                   diffconc0(ic0,ic1,ic2,ijc) *
     &                   (conc(ic0,ic1,ic2,jc)-conc(ic0-1,ic1,ic2,jc))
     &                   )
                  enddo
               enddo
            enddo

            do ic2 = ifirst2, ilast2
               do ic1 = ifirst1, ilast1+1
                  do ic0 = ifirst0, ilast0
                     flux1(ic0,ic1,ic2,ic) = flux1(ic0,ic1,ic2,ic)
     &                  +dyinv * (
     &                   diffconc1(ic0,ic1,ic2,ijc) *
     &                   (conc(ic0,ic1,ic2,jc)-conc(ic0,ic1-1,ic2,jc))
     &                   )
                  enddo
               enddo
            enddo

            do ic2 = ifirst2, ilast2+1
               do ic1 = ifirst1, ilast1
                  do ic0 = ifirst0, ilast0
                     flux2(ic0,ic1,ic2,ic) = flux2(ic0,ic1,ic2,ic)
     &                  +dzinv * (
     &                   diffconc2(ic0,ic1,ic2,ijc) *
     &                   (conc(ic0,ic1,ic2,jc)-conc(ic0,ic1,ic2-1,jc))
     &                   )
                  enddo
               enddo
            enddo
         enddo

      enddo

      return
      end

c
c compute the concentration flux
c
      subroutine concentrationflux_spinodal(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   dx,
     &   conc, ngconc,
     &   ncomp,
     &   conca, ngconca,
     &   concb, ngconcb,
     &   diffconc0,  diffconc1, diffconc2, ngdiff,
     &   eta, ngeta,
     $   kappa,
     &   flux0, flux1, flux2, ngflux )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
      integer maxncomp
      parameter ( maxncomp=1 )
c max matrix depth (should be ncomp*ncomp)
      integer mdepth
      parameter ( mdepth=1 )
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngflux
      double precision dx(0:2)
      double precision kappa
      integer ncomp
      integer three_phase
      integer ngconc, ngconca, ngconcb
      integer ngeta, ngdiff
c
c variables in 2d cell indexed
      double precision 
     &     flux0(SIDE3d0(ifirst,ilast,ngflux),maxncomp),
     &     flux1(SIDE3d1(ifirst,ilast,ngflux),maxncomp),
     &     flux2(SIDE3d2(ifirst,ilast,ngflux),maxncomp)
      double precision conc(CELL3d(ifirst,ilast,ngconc),maxncomp)
      double precision conca(CELL3d(ifirst,ilast,ngconca),maxncomp)
      double precision concb(CELL3d(ifirst,ilast,ngconcb),maxncomp)
      double precision eta(CELL3d(ifirst,ilast,ngeta))
      double precision diffconc0(SIDE3d0(ifirst,ilast,ngdiff),
     &                            mdepth)
      double precision diffconc1(SIDE3d1(ifirst,ilast,ngdiff),
     &                            mdepth)
      double precision diffconc2(SIDE3d2(ifirst,ilast,ngdiff),
     &                            mdepth)
c
      double precision dxinv, dyinv, dzinv, sideca, sidecb
      integer          ic0, ic1, ic2, ic, jc, ijc

      dxinv = 1.d0 / dx(0)
      dyinv = 1.d0 / dx(1)
      dzinv = 1.d0 / dx(2)
c      print*,kappa

      do ic = 1, ncomp
         do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0+1
            flux0(ic0,ic1,ic2,ic) = 0.d0
            do jc = 1, ncomp
               ijc=ic+(jc-1)*ncomp
               sideca = 0.5*(
     &            conca(ic0-1,ic1,ic2,jc)+conca(ic0,ic1,ic2,jc) )
               sidecb = 0.5*(
     &            concb(ic0-1,ic1,ic2,jc)+concb(ic0,ic1,ic2,jc) )
               flux0(ic0,ic1,ic2,ic) = flux0(ic0,ic1,ic2,ic)
     &         + dxinv*(
     &           diffconc0(ic0,ic1,ic2,ijc) *
     &           ( conc(ic0,ic1,ic2,jc) - conc(ic0-1,ic1,ic2,jc) )
     &         - kappa * (sideca-sidecb)
     &           *( eta(ic0,ic1,ic2) - eta(ic0-1,ic1,ic2) )
     &           )
            enddo
         enddo
         enddo
         enddo

         do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1+1
         do ic0 = ifirst0, ilast0
            flux1(ic0,ic1,ic2,ic) = 0.d0
            do jc = 1, ncomp
               ijc=ic+(jc-1)*ncomp
               sideca = 0.5*(
     &            conca(ic0,ic1-1,ic2,jc)+conca(ic0,ic1,ic2,jc) )
               sidecb = 0.5*(
     &            concb(ic0,ic1-1,ic2,jc)+concb(ic0,ic1,ic2,jc) )
               flux1(ic0,ic1,ic2,ic) = flux1(ic0,ic1,ic2,ic)
     &         + dyinv*(
     &           diffconc1(ic0,ic1,ic2,ijc) *
     &           ( conc(ic0,ic1,ic2,jc) - conc(ic0,ic1-1,ic2,jc) )
     &         - kappa * (sideca-sidecb)
     &           *( eta(ic0,ic1,ic2) - eta(ic0,ic1-1,ic2) )
     &           )
            enddo
         enddo
         enddo
         enddo

         do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0
            flux1(ic0,ic1,ic2,ic) = 0.d0
            do jc = 1, ncomp
               ijc=ic+(jc-1)*ncomp
               sideca = 0.5*(
     &            conca(ic0,ic1,ic2-1,jc)+conca(ic0,ic1,ic2,jc) )
               sidecb = 0.5*(
     &            concb(ic0,ic1,ic2-1,jc)+concb(ic0,ic1,ic2,jc) )
               flux1(ic0,ic1,ic2,ic) = flux1(ic0,ic1,ic2,ic)
     &         + dyinv*(
     &           diffconc1(ic0,ic1,ic2,ijc) *
     &           ( conc(ic0,ic1,ic2,jc) - conc(ic0,ic1,ic2-1,jc) )
     &         - kappa * (sideca-sidecb)
     &           *( eta(ic0,ic1,ic2) - eta(ic0,ic1,ic2-1) )
     &           )
            enddo
         enddo
         enddo
         enddo

      enddo

      return
      end

c***********************************************************************
c
c compute r.h.s. for concentration variable
c
      subroutine computerhsconcentration(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   dx,
     &   flux0, flux1, flux2, ngflux,
     &   mobility,
     &   rhs, ngrhs)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &     ngflux, ngrhs
      double precision
     &     dx(0:2),
     &     flux0(SIDE3d0(ifirst,ilast,ngflux)),
     &     flux1(SIDE3d1(ifirst,ilast,ngflux)),
     &     flux2(SIDE3d2(ifirst,ilast,ngflux)),
     &     mobility,
     &     rhs(CELL3d(ifirst,ilast,ngrhs))

c local variables:
      double precision dxinv, dyinv, dzinv
      integer          ic0, ic1, ic2

      dxinv = 1.d0 / dx(0)
      dyinv = 1.d0 / dx(1)
      dzinv = 1.d0 / dx(2)

      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0
               rhs(ic0,ic1,ic2) = mobility * (
     &              dxinv * (flux0(ic0+1,ic1,ic2) - flux0(ic0,ic1,ic2))
     &            + dyinv * (flux1(ic0,ic1+1,ic2) - flux1(ic0,ic1,ic2))
     &            + dzinv * (flux2(ic0,ic1,ic2+1) - flux2(ic0,ic1,ic2))
     &              )
            enddo
         enddo
      enddo

      return
      end

c***********************************************************************
c
c compute the concentration flux
c
      subroutine addconcentrationfluxfromgradt(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   dx,
     &   temperature, ngt,
     &   mq0, mq1, mq2, ngmq,
     &   flux0, flux1, flux2, ngflux,
     &   avg_type )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngflux, ngt, ngmq
      double precision 
     &     flux0(SIDE3d0(ifirst,ilast,ngflux)),
     &     flux1(SIDE3d1(ifirst,ilast,ngflux)),
     &     flux2(SIDE3d2(ifirst,ilast,ngflux))
      double precision 
     &     mq0(SIDE3d0(ifirst,ilast,ngmq)),
     &     mq1(SIDE3d1(ifirst,ilast,ngmq)),
     &     mq2(SIDE3d2(ifirst,ilast,ngmq))
      double precision dx(0:2)
      character*(*) avg_type

      double precision temperature(CELL3d(ifirst,ilast,ngt))

c
      double precision dxinv, dyinv, dzinv, sideT
      integer          ic0, ic1, ic2
      double precision average_func

      dxinv = 1.d0 / dx(0)
      dyinv = 1.d0 / dx(1)
      dzinv = 1.d0 / dx(2)

      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1
               sideT = average_func(
     &            temperature(ic0-1,ic1,ic2), temperature(ic0,ic1,ic2), 
     &            avg_type )
            
               flux0(ic0,ic1,ic2) = flux0(ic0,ic1,ic2) -
     &           dxinv * mq0(ic0,ic1,ic2) * (
     &           ( temperature(ic0,  ic1,ic2) 
     &           - temperature(ic0-1,ic1,ic2) )
     &           / sideT
     &           )
            enddo
         enddo
      enddo

      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0, ilast0
               sideT = average_func(
     &            temperature(ic0,ic1-1,ic2), temperature(ic0,ic1,ic2), 
     &            avg_type )

            flux1(ic0,ic1,ic2) = flux1(ic0,ic1,ic2) -
     &           dyinv * mq1(ic0,ic1,ic2) * (
     &           ( temperature(ic0,ic1,  ic2) 
     &           - temperature(ic0,ic1-1,ic2) )
     &           / sideT
     &           )
            enddo
         enddo
      enddo

      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0
               sideT = average_func(
     &            temperature(ic0,ic1,ic2-1), temperature(ic0,ic1,ic2), 
     &            avg_type )

            flux2(ic0,ic1,ic2) = flux2(ic0,ic1,ic2) -
     &           dzinv * mq2(ic0,ic1,ic2) * (
     &           ( temperature(ic0,ic1,ic2)
     &           - temperature(ic0,ic1,ic2-1) )
     &           / sideT
     &           )
            enddo
         enddo
      enddo

      return
      end

c***********************************************************************
c
c compute the concentration flux
c 0.25 coeff is 0.5 for average over 2 gradients times 0.5
c for 1./2dx
c
      subroutine addconcentrationfluxfromantitrapping(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   dx,
     &   phase, ngp,
     &   cl, ca, ngc,
     &   ncomp,
     &   dphidt, ngd,
     &   alpha,
     &   flux0, flux1, flux2, ngflux )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngflux, ngd, ngp, ngc, ncomp
      double precision 
     &     flux0(SIDE3d0(ifirst,ilast,ngflux),ncomp),
     &     flux1(SIDE3d1(ifirst,ilast,ngflux),ncomp),
     &     flux2(SIDE3d2(ifirst,ilast,ngflux),ncomp)
      double precision phase(CELL3d(ifirst,ilast,ngp))
      double precision cl(CELL3d(ifirst,ilast,ngc),ncomp)
      double precision ca(CELL3d(ifirst,ilast,ngc),ncomp)
      double precision dphidt(CELL3d(ifirst,ilast,ngd))
      double precision dx(0:2)
      double precision alpha
c
      double precision dxinv, dyinv, dzinv
      double precision dphix, dphiy, dphiz, dphi2, dphin
      integer          ic, ic0, ic1, ic2
      double precision tol, tol2, dpdt, acl, aca

      tol = 1.e-8
      tol2 = tol*tol

      dxinv = 1.d0 / dx(0)
      dyinv = 1.d0 / dx(1)
      dzinv = 1.d0 / dx(2)

      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1
               dphix = dxinv 
     &               * ( phase(ic0,ic1,ic2)-phase(ic0-1,ic1,ic2) )
               dphiy = dyinv*0.25d0*(
     &               phase(ic0,  ic1+1,ic2)-phase(ic0,  ic1-1,ic2)
     &             + phase(ic0-1,ic1+1,ic2)-phase(ic0-1,ic1-1,ic2)
     &                )
               dphiz = dzinv*0.25d0*(
     &               phase(ic0,  ic1,ic2+1)-phase(ic0,  ic1,ic2-1)
     &             + phase(ic0-1,ic1,ic2+1)-phase(ic0-1,ic1,ic2-1)
     &              )

               dphi2 = dphix*dphix+dphiy*dphiy+dphiz*dphiz
               if( abs(dphi2) .gt. tol2 ) then
                  dphin = sqrt(dphi2)
                  dpdt=0.5d0*(dphidt(ic0-1,ic1,ic2)+dphidt(ic0,ic1,ic2))

                  do ic = 1, ncomp
                     acl=0.5d0*(cl(ic0-1,ic1,ic2,ic)
     &                         +cl(ic0,ic1,ic2,ic))
                     aca=0.5d0*(ca(ic0-1,ic1,ic2,ic)
     &                         +ca(ic0,ic1,ic2,ic))

                     flux0(ic0,ic1,ic2,ic) = flux0(ic0,ic1,ic2,ic) +
     &                  alpha*(dphix/dphin)*(acl-aca)*dpdt 
                  enddo
               endif
            enddo
         enddo
      enddo

      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0, ilast0
               dphiy = dyinv 
     &            * ( phase(ic0,ic1,ic2)-phase(ic0,ic1-1,ic2) )
               dphix = dxinv*0.25d0*(
     &            phase(ic0+1,ic1-1,ic2)-phase(ic0-1,ic1-1,ic2)
     &          + phase(ic0+1,ic1  ,ic2)-phase(ic0-1,ic1  ,ic2)
     &              )
               dphiz = dzinv*0.25d0*(
     &            phase(ic0,ic1-1,ic2+1)-phase(ic0,ic1-1,ic2-1)
     &          + phase(ic0,ic1,  ic2+1)-phase(ic0,ic1,  ic2-1)
     &              )

               dphi2 = dphix*dphix+dphiy*dphiy+dphiz*dphiz
               if( abs(dphi2) .gt. tol2 ) then
                  dphin = sqrt(dphi2)
                  dpdt=0.5d0*(dphidt(ic0,ic1-1,ic2)+dphidt(ic0,ic1,ic2))
                  do ic = 1, ncomp
                     acl=0.5d0*(cl(ic0,ic1-1,ic2,ic)
     &                         +cl(ic0,ic1,ic2,ic))
                     aca=0.5d0*(ca(ic0,ic1-1,ic2,ic)
     &                         +ca(ic0,ic1,ic2,ic))
                     flux1(ic0,ic1,ic2,ic) = flux1(ic0,ic1,ic2,ic) +
     &                  alpha*(dphiy/dphin)*(acl-aca)*dpdt
                  enddo
               endif
            enddo
         enddo
      enddo

      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0
               dphiz = dzinv 
     &            * ( phase(ic0,ic1,ic2)-phase(ic0,ic1,ic2-1) )
               dphix = dxinv*0.25d0*(
     &            phase(ic0+1,ic1,ic2-1)-phase(ic0-1,ic1,ic2-1)
     &          + phase(ic0+1,ic1,ic2  )-phase(ic0-1,ic1,ic2  )
     &              )
               dphiy = dyinv*0.25d0*(
     &               phase(ic0,ic1+1,ic2-1)-phase(ic0,ic1-1,ic2-1)
     &             + phase(ic0,ic1+1,ic2  )-phase(ic0,ic1-1,ic2  )
     &                )

               dphi2 = dphix*dphix+dphiy*dphiy+dphiz*dphiz
               if( abs(dphi2) .gt. tol2 ) then
                  dphin = sqrt(dphi2)
                  dpdt=0.5d0*(dphidt(ic0,ic1,ic2-1)+dphidt(ic0,ic1,ic2))
               
                  do ic = 1, ncomp
                     acl=0.5d0*(cl(ic0,ic1,ic2-1,ic)
     &                         +cl(ic0,ic1,ic2,ic))
                     aca=0.5d0*(ca(ic0,ic1,ic2-1,ic)
     &                         +ca(ic0,ic1,ic2,ic))
                     flux2(ic0,ic1,ic2,ic) = flux2(ic0,ic1,ic2,ic) +
     &                  alpha*(dphiz/dphin)*(acl-aca)*dpdt
                  enddo
               endif
            enddo
         enddo
      enddo

      return
      end

c***********************************************************************
c
c compute the concentration flux
c 0.25 coeff is 0.5 for average over 2 gradients times 0.5
c for 1./2dx
c
      subroutine addconcentrationfluxfromantitrapping3phases(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   dx,
     &   phase, ngp,
     &   cl, ca, cb, ngc,
     &   dphidt, ngd,
     &   alpha,
     &   flux0, flux1, flux2, ngflux )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngflux, ngd, ngp, ngc, ncomp
      double precision
     &     flux0(SIDE3d0(ifirst,ilast,ngflux)),
     &     flux1(SIDE3d1(ifirst,ilast,ngflux)),
     &     flux2(SIDE3d2(ifirst,ilast,ngflux))
      double precision phase(CELL3d(ifirst,ilast,ngp),3)
      double precision cl(CELL3d(ifirst,ilast,ngc))
      double precision ca(CELL3d(ifirst,ilast,ngc))
      double precision cb(CELL3d(ifirst,ilast,ngc))
      double precision dphidt(CELL3d(ifirst,ilast,ngd),3)
      double precision dx(0:2)
      double precision alpha
c
      double precision dxinv, dyinv, dzinv
      double precision dphix, dphiy, dphiz, dphi2, dphin
      double precision dphipx, dphipy, dphipz
      double precision dphipn, dphip2
      integer          ic0, ic1, ic2, ip
      double precision tol, tol2, dpdt
      double precision ac(3)
      double precision factor

      tol = 1.e-8
      tol2 = tol*tol

      dxinv = 1.d0 / dx(0)
      dyinv = 1.d0 / dx(1)
      dzinv = 1.d0 / dx(2)

c x-faces
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1
c compute gradient liquid phase first
               dphix = dxinv
     &               * ( phase(ic0,ic1,ic2,1)-phase(ic0-1,ic1,ic2,1) )
               dphiy = dyinv*0.25d0*(
     &               phase(ic0,  ic1+1,ic2,1)- phase(ic0,  ic1-1,ic2,1)
     &             + phase(ic0-1,ic1+1,ic2,1)- phase(ic0-1,ic1-1,ic2,1)
     &                )
               dphiz = dzinv*0.25d0*(
     &               phase(ic0,  ic1,ic2+1,1)-phase(ic0,  ic1,ic2-1,1)
     &             + phase(ic0-1,ic1,ic2+1,1)-phase(ic0-1,ic1,ic2-1,1)
     &              )
               dphi2 = dphix*dphix+dphiy*dphiy+dphiz*dphiz
               if( abs(dphi2) .gt. tol2 ) then
                  dphin = sqrt(dphi2)

                  ac(1)=0.5d0*(cl(ic0-1,ic1,ic2)+cl(ic0,ic1,ic2))
                  ac(2)=0.5d0*(ca(ic0-1,ic1,ic2)+ca(ic0,ic1,ic2))
                  ac(3)=0.5d0*(cb(ic0-1,ic1,ic2)+cb(ic0,ic1,ic2))

c loop over two solid phases
                  do ip = 2, 3
                     dpdt=0.5d0*(dphidt(ic0-1,ic1,ic2,ip)
     &                          +dphidt(ic0,ic1,ic2,ip))

                     dphipx = dxinv * ( phase(ic0,ic1,ic2,ip)
     &                                - phase(ic0-1,ic1,ic2,ip) )
                     dphipy = dyinv*0.25d0
     &                      * (phase(ic0-1,ic1+1,ic2,ip)
     &                        -phase(ic0-1,ic1-1,ic2,ip)
     &                        +phase(ic0  ,ic1+1,ic2,ip)
     &                        -phase(ic0  ,ic1-1,ic2,ip))
                     dphipz = dzinv*0.25d0
     &                      * (phase(ic0-1,ic1,ic2+1,ip)
     &                        -phase(ic0-1,ic1,ic2-1,ip)
     &                        +phase(ic0  ,ic1,ic2+1,ip)
     &                        -phase(ic0  ,ic1,ic2-1,ip))

                     dphip2 = dphipx*dphipx+dphipy*dphipy+dphipz*dphipz

                     if( abs(dphip2) .gt. tol2 ) then
                        dphipn = sqrt(dphip2)
c factor should be one when only two phases are present
                        factor=-1.d0*(dphipx*dphix+dphipy*dphiy
     &                               +dphipz*dphiz)
     &                        /(dphipn*dphin)
                        flux0(ic0,ic1,ic2) = flux0(ic0,ic1,ic2)
     &                     + alpha*factor*(dphix/dphin)
     &                       *(ac(1)-ac(ip))*dpdt
                     endif
                  enddo
               endif
            enddo
         enddo
      enddo

c y-faces
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0, ilast0
c compute gradient liquid phase first
               dphiy = dyinv
     &            * ( phase(ic0,ic1,ic2,1)-phase(ic0,ic1-1,ic2,1) )
               dphix = dxinv*0.25d0*(
     &            phase(ic0+1,ic1-1,ic2,1)-phase(ic0-1,ic1-1,ic2,1)
     &          + phase(ic0+1,ic1  ,ic2,1)-phase(ic0-1,ic1  ,ic2,1)
     &              )
               dphiz = dzinv*0.25d0*(
     &            phase(ic0,ic1-1,ic2+1,1)-phase(ic0,ic1-1,ic2-1,1)
     &          + phase(ic0,ic1,  ic2+1,1)-phase(ic0,ic1,  ic2-1,1)
     &              )
               dphi2 = dphix*dphix+dphiy*dphiy+dphiz*dphiz
               if( abs(dphi2) .gt. tol2 ) then
                  dphin = sqrt(dphi2)
                  ac(1)=0.5d0*(cl(ic0,ic1-1,ic2)+cl(ic0,ic1,ic2))
                  ac(2)=0.5d0*(ca(ic0,ic1-1,ic2)+ca(ic0,ic1,ic2))
                  ac(3)=0.5d0*(cb(ic0,ic1-1,ic2)+cb(ic0,ic1,ic2))

c loop over two solid phases
                  do ip = 2, 3
                     dpdt=0.5d0*(dphidt(ic0,ic1-1,ic2,ip)
     &                          +dphidt(ic0,ic1,ic2,ip))
                     dphipy = dyinv
     &                  * ( phase(ic0,ic1,ic2,ip)
     &                     -phase(ic0,ic1-1,ic2,ip) )
                     dphipx = dxinv*0.25d0*(
     &                  phase(ic0+1,ic1-1,ic2,ip)
     &                - phase(ic0-1,ic1-1,ic2,ip)
     &                + phase(ic0+1,ic1  ,ic2,ip)
     &                - phase(ic0-1,ic1  ,ic2,ip)
     &                    )
                     dphipz = dzinv*0.25d0*(
     &                  phase(ic0,ic1-1,ic2+1,ip)
     &                - phase(ic0,ic1-1,ic2-1,ip)
     &                + phase(ic0,ic1,  ic2+1,ip)
     &                - phase(ic0,ic1,  ic2-1,ip)
     &                    )
                     dphip2 = dphipx*dphipx+dphipy*dphipy+dphipz*dphipz

                     if( abs(dphip2) .gt. tol2 ) then
                        dphipn = sqrt(dphip2)
c factor should be one when only two phases are present
                        factor=-1.d0*(dphipx*dphix+dphipy*dphiy
     &                               +dphipz*dphiz)
     &                        /(dphipn*dphin)
                        flux1(ic0,ic1,ic2) = flux1(ic0,ic1,ic2)
     &                     + alpha*factor*(dphiy/dphin)
     &                       *(ac(1)-ac(ip))*dpdt
                     endif
                  enddo
               endif
            enddo
         enddo
      enddo

c z-face
      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0
c compute gradient liquid phase first
               dphiz = dzinv
     &            * ( phase(ic0,ic1,ic2,1)-phase(ic0,ic1,ic2-1,1) )
               dphix = dxinv*0.25d0*(
     &            phase(ic0+1,ic1,ic2-1,1)-phase(ic0-1,ic1,ic2-1,1)
     &          + phase(ic0+1,ic1,ic2,1  )-phase(ic0-1,ic1,ic2,1  )
     &              )
               dphiy = dyinv*0.25d0*(
     &               phase(ic0,ic1+1,ic2-1,1)-phase(ic0,ic1-1,ic2-1,1)
     &             + phase(ic0,ic1+1,ic2,1  )-phase(ic0,ic1-1,ic2,1  )
     &                )
               dphi2 = dphix*dphix+dphiy*dphiy+dphiz*dphiz
               if( abs(dphi2) .gt. tol2 ) then
                  dphin = sqrt(dphi2)
                  ac(1)=0.5d0*(cl(ic0,ic1,ic2-1)+cl(ic0,ic1,ic2))
                  ac(2)=0.5d0*(ca(ic0,ic1,ic2-1)+ca(ic0,ic1,ic2))
                  ac(3)=0.5d0*(cb(ic0,ic1,ic2-1)+cb(ic0,ic1,ic2))

c loop over two solid phases
                  do ip = 2, 3
                     dpdt=0.5d0*(dphidt(ic0,ic1,ic2-1,ip)
     &                          +dphidt(ic0,ic1,ic2,ip))
                     dphipz = dzinv
     &                  * ( phase(ic0,ic1,ic2,ip)
     &                     -phase(ic0,ic1,ic2-1,ip) )
                     dphipx = dxinv*0.25d0*(
     &                  phase(ic0+1,ic1,ic2-1,ip)
     &                - phase(ic0-1,ic1,ic2-1,ip)
     &                + phase(ic0+1,ic1,ic2,ip)
     &                - phase(ic0-1,ic1,ic2,ip)
     &                    )
                     dphipy = dyinv*0.25d0*(
     &                  phase(ic0,ic1+1,ic2-1,ip)
     &                - phase(ic0,ic1-1,ic2-1,ip)
     &                + phase(ic0,ic1+1,ic2,ip)
     &                - phase(ic0,ic1-1,ic2,ip)
     &                    )
                     dphip2 = dphipx*dphipx+dphipy*dphipy+dphipz*dphipz

                     if( abs(dphip2) .gt. tol2 ) then
                        dphipn = sqrt(dphip2)
c factor should be one when only two phases are present
                        factor=-1.d0*(dphipx*dphix+dphipy*dphiy
     &                               +dphipz*dphiz)
     &                        /(dphipn*dphin)
                        flux2(ic0,ic1,ic2) = flux2(ic0,ic1,ic2)
     &                    + alpha*factor*(dphiz/dphin)
     &                      *(ac(1)-ac(ip))*dpdt
                     endif
                  enddo
               endif
            enddo
         enddo
      enddo

      return
      end
