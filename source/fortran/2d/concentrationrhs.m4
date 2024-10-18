c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
define(NDIM,2)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl

c***********************************************************************
c
c compute the concentration flux
c
      subroutine concentrationflux(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   dx,
     &   conc, ngconc,
     &   phi, ngphi,
     &   diffconc0,  diffconc1,  ngdiffconc,
     &   dphicoupl0, dphicoupl1, ngdphicoupl,
     &   flux0, flux1, ngflux)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ngflux
      double precision 
     &     flux0(SIDE2d0(ifirst,ilast,ngflux)),
     &     flux1(SIDE2d1(ifirst,ilast,ngflux))
      double precision dx(0:1)
c
      integer ngconc, ngphi
      integer ngdiffconc, ngdphicoupl
c
c variables in 2d cell indexed
      double precision conc(CELL2d(ifirst,ilast,ngconc))
      double precision phi(CELL2d(ifirst,ilast,ngphi))
      double precision diffconc0(SIDE2d0(ifirst,ilast,ngdiffconc))
      double precision diffconc1(SIDE2d1(ifirst,ilast,ngdiffconc))

c variables in 2d side indexed
      double precision
     &   dphicoupl0(SIDE2d0(ifirst,ilast,ngdphicoupl))
      double precision
     &   dphicoupl1(SIDE2d1(ifirst,ilast,ngdphicoupl))
c
      double precision dxinv, dyinv
      integer          ic0, ic1

      dxinv = 1.d0 / dx(0)

      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0+1
            flux0(ic0,ic1) = dxinv * (
     &           diffconc0(ic0,ic1) *
     &           ( conc(ic0,ic1) - conc(ic0-1,ic1) )
     &         + dphicoupl0(ic0,ic1) *
     &           ( phi(ic0,ic1) - phi(ic0-1,ic1) )
     &           )
         enddo
      enddo

      dyinv = 1.d0 / dx(1)

      do ic1 = ifirst1, ilast1+1
         do ic0 = ifirst0, ilast0
            flux1(ic0,ic1) = dyinv * (
     &           diffconc1(ic0,ic1) *
     &           ( conc(ic0,ic1) - conc(ic0,ic1-1) )
     &         + dphicoupl1(ic0,ic1) *
     &           ( phi(ic0,ic1) - phi(ic0,ic1-1) )
     &           )
         enddo
      enddo

      return
      end

c***********************************************************************
c
c Cahn-Hilliard double well flux
c
      subroutine add_cahnhilliarddoublewell_flux(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   dx,
     &   conc, ngconc,
     &   mobility,
     &   ca, cb, well_scale, kappa,
     &   flux0, flux1, ngflux)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ngflux
      double precision conc(CELL2d(ifirst,ilast,ngconc))
      double precision
     &     flux0(SIDE2d0(ifirst,ilast,ngflux)),
     &     flux1(SIDE2d1(ifirst,ilast,ngflux))
      double precision dx(0:1)
c
      integer ngconc
      double precision mobility, ca, cb, well_scale, kappa
c
      double precision dxinv, dyinv, dxinv2, dyinv2
      double precision c, lap, mu
      integer          ic0, ic1

      dxinv = 1.d0 / dx(0)
      dyinv = 1.d0 / dx(1)

      dxinv2 = dxinv*dxinv
      dyinv2 = dyinv*dyinv

      do ic1 = ifirst1-1, ilast1+1
         do ic0 = ifirst0-1, ilast0+1
            lap = dxinv2*(-2.d0*conc(ic0,ic1)
     &                  + conc(ic0-1,ic1) +conc(ic0+1,ic1))
     &          + dyinv2*(-2.d0*conc(ic0,ic1)
     &                  + conc(ic0,ic1-1) +conc(ic0,ic1+1))
            c = conc(ic0,ic1)
            mu = 2.d0*well_scale*(c-ca)*(cb-c)*(cb+ca-2.d0*c)-kappa*lap
            flux0(ic0,ic1)   = flux0(ic0,ic1)
     &           + mobility*dxinv * mu
            flux0(ic0+1,ic1) = flux0(ic0+1,ic1)
     &           - mobility*dxinv * mu
            flux1(ic0,ic1)   = flux1(ic0,ic1)
     &           + mobility*dyinv * mu
            flux1(ic0,ic1+1) = flux1(ic0,ic1+1)
     &           - mobility*dyinv * mu
         enddo
      enddo

      return
      end

c***********************************************************************
c
c Wang sintering double well flux
c
      subroutine add_wang_sintering_flux(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   dx,
     &   conc, ngconc,
     &   diff0, diff1, ngdiff,
     &   phi, ngphi, norder,
     &   mobility,
     &   parameter_a, parameter_b, kappa,
     &   flux0, flux1, ngflux,
     &   phi2sum, phi3sum)
c***********************************************************************
      implicit none
c***********************************************************************
      integer ngconc, ngphi, ngflux, ngdiff
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      double precision conc(CELL2d(ifirst,ilast,ngconc))
      double precision diff0(SIDE2d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE2d1(ifirst,ilast,ngdiff))
      double precision phi(CELL2d(ifirst,ilast,ngphi),norder)
      double precision flux0(SIDE2d0(ifirst,ilast,ngflux))
      double precision flux1(SIDE2d1(ifirst,ilast,ngflux))
      double precision phi2sum(CELL2d(ifirst,ilast,ngphi))
      double precision phi3sum(CELL2d(ifirst,ilast,ngphi))
      double precision dx(0:1)
c
      double precision mobility, parameter_a, parameter_b, kappa
c
      double precision dxinv, dyinv, dxinv2, dyinv2
      double precision c, lap, mu
      integer          ic0, ic1, ip, norder

      dxinv = 1.d0 / dx(0)
      dyinv = 1.d0 / dx(1)
      dxinv2 = dxinv*dxinv
      dyinv2 = dyinv*dyinv
c
c precompute sum of phi**2 at each cell
c precompute sum of phi**3 at each cell
      do ic1 = ifirst1-1, ilast1+1
         do ic0 = ifirst0-1, ilast0+1
            phi2sum(ic0,ic1) =
     &         phi(ic0,ic1,1)*phi(ic0,ic1,1)
            phi3sum(ic0,ic1) =
     &         phi(ic0,ic1,1)*phi(ic0,ic1,1)*phi(ic0,ic1,1)
         enddo
      enddo
c
      do ip = 2, norder
         do ic1 = ifirst1-1, ilast1+1
            do ic0 = ifirst0-1, ilast0+1
               phi2sum(ic0,ic1) = phi2sum(ic0,ic1)
     &            + phi(ic0,ic1,ip)*phi(ic0,ic1,ip)
               phi3sum(ic0,ic1) = phi3sum(ic0,ic1)
     &            + phi(ic0,ic1,ip)*phi(ic0,ic1,ip)*phi(ic0,ic1,ip)
            enddo
         enddo
      enddo

      do ic1 = ifirst1-1, ilast1+1
         do ic0 = ifirst0-1, ilast0+1
            lap = dxinv2*(-2.d0*conc(ic0,ic1)
     &                  + conc(ic0-1,ic1) +conc(ic0+1,ic1))
     &          + dyinv2*(-2.d0*conc(ic0,ic1)
     &                  + conc(ic0,ic1-1) +conc(ic0,ic1+1))
            c = conc(ic0,ic1)
            mu = 2.d0*parameter_a*c*(1.d0-c)*(1.d0-2.d0*c)
     &           +2.d0*parameter_b*(
     &           c-3.d0*phi2sum(ic0,ic1)+2.d0*phi3sum(ic0,ic1))
     &           -kappa*lap
            flux0(ic0,ic1)   = flux0(ic0,ic1)
     &           + mobility*diff0(ic0,ic1)*dxinv*mu
            flux0(ic0+1,ic1) = flux0(ic0+1,ic1)
     &           - mobility*diff0(ic0+1,ic1)*dxinv*mu
            flux1(ic0,ic1)   = flux1(ic0,ic1)
     &           + mobility*diff1(ic0,ic1)*dyinv*mu
            flux1(ic0,ic1+1) = flux1(ic0,ic1+1)
     &           - mobility*diff1(ic0,ic1+1)*dyinv*mu
         enddo
      enddo

      return
      end

c***********************************************************************
c
c compute the concentration flux
c
      subroutine concentrationflux_spinodal(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   dx,
     &   conc, ngconc,
     &   ncomp,
     &   conca, ngconca,
     &   concb, ngconcb,
     &   diffconc0,  diffconc1,  ngdiff,
     &   eta, ngeta,
     $   kappa,
     &   flux0, flux1, ngflux )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
      integer maxncomp
      parameter ( maxncomp=1 )
c max matrix depth (should be ncomp*ncomp)
      integer mdepth
      parameter ( mdepth=1 )
      integer ifirst0, ilast0, ifirst1, ilast1, ngflux
      double precision dx(0:1)
      double precision kappa
      integer ncomp
      integer ngconc, ngconca, ngconcb
      integer ngeta, ngdiff
c
c variables in 2d cell indexed
      double precision 
     &     flux0(SIDE2d0(ifirst,ilast,ngflux),maxncomp),
     &     flux1(SIDE2d1(ifirst,ilast,ngflux),maxncomp)
      double precision conc(CELL2d(ifirst,ilast,ngconc),maxncomp)
      double precision conca(CELL2d(ifirst,ilast,ngconca),maxncomp)
      double precision concb(CELL2d(ifirst,ilast,ngconcb),maxncomp)
      double precision eta(CELL2d(ifirst,ilast,ngeta))
      double precision diffconc0(SIDE2d0(ifirst,ilast,ngdiff),
     &                            mdepth)
      double precision diffconc1(SIDE2d1(ifirst,ilast,ngdiff),
     &                            mdepth)
c
      double precision dxinv, dyinv, sideca, sidecb
      integer          ic0, ic1, ic, jc, ijc

      dxinv = 1.d0 / dx(0)
      dyinv = 1.d0 / dx(1)
c      print*,kappa

      do ic = 1, ncomp
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1
               flux0(ic0,ic1,ic) = 0.d0
               do jc = 1, ncomp
                  ijc=ic+(jc-1)*ncomp
                  sideca = 0.5*(
     &               conca(ic0-1,ic1,jc)+conca(ic0,ic1,jc) )
                  sidecb = 0.5*(
     &               concb(ic0-1,ic1,jc)+concb(ic0,ic1,jc) )
                  flux0(ic0,ic1,ic) = flux0(ic0,ic1,ic) + dxinv * (
     &              diffconc0(ic0,ic1,ijc) *
     &              ( conc(ic0,ic1,jc) - conc(ic0-1,ic1,jc) )
     &            - kappa * (sideca-sidecb)
     &              *( eta(ic0,ic1) - eta(ic0-1,ic1) )
     &              )
               enddo
            enddo
         enddo

         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0, ilast0
               flux1(ic0,ic1,ic) = 0.d0
               do jc = 1, ncomp
                  ijc=ic+(jc-1)*ncomp
                  sideca = 0.5*(
     &               conca(ic0,ic1-1,jc)+conca(ic0,ic1,jc) )
                  sidecb = 0.5*(
     &               concb(ic0,ic1-1,jc)+concb(ic0,ic1,jc) )
                  flux1(ic0,ic1,ic) = flux1(ic0,ic1,ic) + dyinv * (
     &              diffconc1(ic0,ic1,ijc) *
     &              ( conc(ic0,ic1,jc) - conc(ic0,ic1-1,jc) )
     &            - kappa * (sideca-sidecb)
     &              *( eta(ic0,ic1) - eta(ic0,ic1-1) )
     &              )
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
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   dx,
     &   flux0, flux1, ngflux,
     &   mobility,
     &   rhs, ngrhs)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1,
     &     ngflux, ngrhs
      double precision
     &     dx(0:1),
     &     flux0(SIDE2d0(ifirst,ilast,ngflux)),
     &     flux1(SIDE2d1(ifirst,ilast,ngflux)),
     &     mobility,
     &     rhs(CELL2d(ifirst,ilast,ngrhs))

c local variables:
      double precision dxinv, dyinv
      integer          ic0, ic1

      dxinv = 1.d0 / dx(0)
      dyinv = 1.d0 / dx(1)

      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0
            rhs(ic0,ic1) = mobility * (
     &           dxinv * (flux0(ic0+1,ic1) - flux0(ic0,ic1))
     &         + dyinv * (flux1(ic0,ic1+1) - flux1(ic0,ic1))
     &           )
         enddo
      enddo

      return
      end

c***********************************************************************
c
c compute the concentration flux
c
      subroutine addconcentrationfluxfromgradt(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   dx,
     &   temperature, ngt,
     &   mq0, mq1, ngmq,
     &   flux0, flux1, ngflux,
     &   avg_type )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer ngflux, ngt, ngmq
      double precision 
     &     flux0(SIDE2d0(ifirst,ilast,ngflux)),
     &     flux1(SIDE2d1(ifirst,ilast,ngflux))
      double precision 
     &     mq0(SIDE2d0(ifirst,ilast,ngmq)),
     &     mq1(SIDE2d1(ifirst,ilast,ngmq))
      double precision dx(0:1)
      character*(*) avg_type

      double precision temperature(CELL2d(ifirst,ilast,ngt))

c
      double precision dxinv, dyinv, sideT
      integer          ic0, ic1
      double precision average_func

      dxinv = 1.d0 / dx(0)
      dyinv = 1.d0 / dx(1)

      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0+1
            sideT = average_func(
     &         temperature(ic0-1,ic1), temperature(ic0,ic1), 
     &         avg_type )
            
            flux0(ic0,ic1) = flux0(ic0,ic1) -
     &           dxinv * mq0(ic0,ic1) * (
     &           ( temperature(ic0,ic1) - temperature(ic0-1,ic1) )
     &           / sideT
     &           )
         enddo
      enddo

      do ic1 = ifirst1, ilast1+1
         do ic0 = ifirst0, ilast0
            sideT = average_func(
     &         temperature(ic0,ic1-1), temperature(ic0,ic1), 
     &         avg_type )

            flux1(ic0,ic1) = flux1(ic0,ic1) -
     &           dyinv * mq1(ic0,ic1) * (
     &           ( temperature(ic0,ic1) - temperature(ic0,ic1-1) )
     &           / sideT
     &           )
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
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   dx,
     &   phase, ngp,
     &   cl, ca, ngc,
     &   ncomp,
     &   dphidt, ngd,
     &   alpha,
     &   flux0, flux1, ngflux )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer ngflux, ngd, ngp, ngc, ncomp
      double precision 
     &     flux0(SIDE2d0(ifirst,ilast,ngflux),ncomp),
     &     flux1(SIDE2d1(ifirst,ilast,ngflux),ncomp)
      double precision phase(CELL2d(ifirst,ilast,ngp))
      double precision cl(CELL2d(ifirst,ilast,ngc),ncomp)
      double precision ca(CELL2d(ifirst,ilast,ngc),ncomp)
      double precision dphidt(CELL2d(ifirst,ilast,ngd))
      double precision dx(0:1)
      double precision alpha
c
      double precision dxinv, dyinv
      double precision dphix, dphiy, dphi2, dphin
      integer          ic, ic0, ic1
      double precision tol, tol2, dpdt, acl, aca

      tol = 1.e-8
      tol2 = tol*tol

      dxinv = 1.d0 / dx(0)
      dyinv = 1.d0 / dx(1)

      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0+1
            dphix = dxinv * ( phase(ic0,ic1) - phase(ic0-1,ic1) )
            dphiy = dyinv*0.25d0*(phase(ic0-1,ic1+1)-phase(ic0-1,ic1-1)
     &                          + phase(ic0  ,ic1+1)-phase(ic0  ,ic1-1)
     &                          )
            dphi2 = dphix*dphix+dphiy*dphiy
            if( abs(dphi2) .gt. tol2 ) then
               dphin = sqrt(dphi2)
               dpdt=0.5d0*( dphidt(ic0-1,ic1)+dphidt(ic0,ic1) )
               
               do ic = 1, ncomp
                  acl=0.5d0*(cl(ic0-1,ic1,ic)+cl(ic0,ic1,ic))
                  aca=0.5d0*(ca(ic0-1,ic1,ic)+ca(ic0,ic1,ic))

                  flux0(ic0,ic1,ic) = flux0(ic0,ic1,ic) +
     &               alpha*(dphix/dphin)*(acl-aca)*dpdt
               enddo
            endif
         enddo
      enddo

      do ic1 = ifirst1, ilast1+1
         do ic0 = ifirst0, ilast0
            dphiy = dyinv * ( phase(ic0,ic1) - phase(ic0,ic1-1) )
            dphix = dxinv*0.25d0*(phase(ic0+1,ic1-1)-phase(ic0-1,ic1-1)
     &                          + phase(ic0+1,ic1  )-phase(ic0-1,ic1  )
     &                           )
            dphi2 = dphix*dphix+dphiy*dphiy
            if( abs(dphi2) .gt. tol2 ) then
               dphin = sqrt(dphi2)
               dpdt=0.5d0*( dphidt(ic0,ic1-1)+dphidt(ic0,ic1) )

               do ic = 1, ncomp
                  acl=0.5d0*(cl(ic0,ic1-1,ic)+cl(ic0,ic1,ic))
                  aca=0.5d0*(ca(ic0,ic1-1,ic)+ca(ic0,ic1,ic))

                  flux1(ic0,ic1,ic) = flux1(ic0,ic1,ic) +
     &               alpha*(dphiy/dphin)*(acl-aca)*dpdt
               enddo
            endif
         enddo
      enddo

      return
      end
c
c***********************************************************************
c
c compute the concentration flux
c 0.25 coeff is 0.5 for average over 2 gradients times 0.5
c for 1./2dx
c
      subroutine addconcentrationfluxfromantitrapping3phases(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   dx,
     &   phase, ngp,
     &   cl, ca, cb, ngc,
     &   dphidt, ngd,
     &   alpha,
     &   flux0, flux1, ngflux )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer ngflux, ngd, ngp, ngc
      double precision
     &     flux0(SIDE2d0(ifirst,ilast,ngflux)),
     &     flux1(SIDE2d1(ifirst,ilast,ngflux))
      double precision phase(CELL2d(ifirst,ilast,ngp),3)
      double precision cl(CELL2d(ifirst,ilast,ngc))
      double precision ca(CELL2d(ifirst,ilast,ngc))
      double precision cb(CELL2d(ifirst,ilast,ngc))
      double precision dphidt(CELL2d(ifirst,ilast,ngd),3)
      double precision dx(0:1)
      double precision alpha
c
      double precision dxinv, dyinv
      double precision dphix, dphiy, dphi2, dphin
      double precision dphipx, dphipy, dphipn, dphip2
      integer          ic0, ic1, ip
      double precision tol, tol2, dpdt
      double precision ac(3)
      double precision factor

      tol = 1.e-8
      tol2 = tol*tol

      dxinv = 1.d0 / dx(0)
      dyinv = 1.d0 / dx(1)

c x-faces
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0+1
c compute gradient liquid phase first
            dphix = dxinv * ( phase(ic0,ic1,1)
     &                      - phase(ic0-1,ic1,1) )
            dphiy = dyinv*0.25d0*(phase(ic0-1,ic1+1,1)
     &                           -phase(ic0-1,ic1-1,1)
     &                           +phase(ic0  ,ic1+1,1)
     &                           -phase(ic0  ,ic1-1,1)
     &                          )
            dphi2 = dphix*dphix+dphiy*dphiy
            if( abs(dphi2) .gt. tol2 ) then
               dphin = sqrt(dphi2)

c average compositions to get values on x-face
               ac(1) = 0.5d0*(cl(ic0-1,ic1)+cl(ic0,ic1))
               ac(2) = 0.5d0*(ca(ic0-1,ic1)+ca(ic0,ic1))
               ac(3) = 0.5d0*(cb(ic0-1,ic1)+cb(ic0,ic1))

c loop over two solid phases
               do ip = 2, 3
                  dphipx = dxinv * ( phase(ic0,ic1,ip)
     &                             - phase(ic0-1,ic1,ip) )
                  dphipy = dyinv*0.25d0
     &                   * (phase(ic0-1,ic1+1,ip)
     &                     -phase(ic0-1,ic1-1,ip)
     &                     +phase(ic0  ,ic1+1,ip)
     &                     -phase(ic0  ,ic1-1,ip))
                  dphip2 = dphipx*dphipx+dphipy*dphipy

                  if( abs(dphip2) .gt. tol2 ) then
                     dpdt=0.5d0*( dphidt(ic0-1,ic1,ip)
     &                           +dphidt(ic0,ic1,ip) )
                     dphipn = sqrt(dphip2)
c factor should be one when only two phases are present
                     factor=-1.d0*(dphipx*dphix+dphipy*dphiy)
     &                     /(dphipn*dphin)

                     flux0(ic0,ic1) = flux0(ic0,ic1) +
     &                  alpha*factor*(dphix/dphin)*(ac(1)-ac(ip))*dpdt
                  endif
               enddo
            endif
         enddo
      enddo

      do ic1 = ifirst1, ilast1+1
         do ic0 = ifirst0, ilast0
c compute gradient liquid phase first
            dphiy = dyinv * ( phase(ic0,ic1,1)
     &                      - phase(ic0,ic1-1,1) )
            dphix = dxinv*0.25d0*(phase(ic0+1,ic1-1,1)
     &                           -phase(ic0-1,ic1-1,1)
     &                           +phase(ic0+1,ic1  ,1)
     &                           -phase(ic0-1,ic1  ,1)
     &                           )
            dphi2 = dphix*dphix+dphiy*dphiy
            if( abs(dphi2) .gt. tol2 ) then
               dphin = sqrt(dphi2)

c average compositions to get values on y-face
               ac(1) = 0.5d0*(cl(ic0,ic1-1)+cl(ic0,ic1))
               ac(2) = 0.5d0*(ca(ic0,ic1-1)+ca(ic0,ic1))
               ac(3) = 0.5d0*(cb(ic0,ic1-1)+cb(ic0,ic1))

c loop over two solid phases
               do ip = 2, 3
                  dphipy = dyinv * ( phase(ic0,ic1,ip)
     &                             - phase(ic0,ic1-1,ip) )
                  dphipx = dxinv*0.25d0
     &                   * (phase(ic0+1,ic1-1,ip)
     &                     -phase(ic0-1,ic1-1,ip)
     &                     +phase(ic0+1,ic1  ,ip)
     &                     -phase(ic0-1,ic1  ,ip))
                  dphip2 = dphipx*dphipx+dphipy*dphipy

                  if( abs(dphip2) .gt. tol2 ) then
                     dpdt=0.5d0*( dphidt(ic0,ic1-1,ip)
     &                           +dphidt(ic0,ic1,ip) )
                     dphipn = sqrt(dphip2)
c factor should be one when only two phases are present
                     factor=-1.d0*(dphipx*dphix+dphipy*dphiy)
     &                     /(dphipn*dphin)

                     flux1(ic0,ic1) = flux1(ic0,ic1) +
     &                  alpha*(factor*dphiy/dphin)*(ac(1)-ac(ip))*dpdt
                  endif
               enddo
            endif
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
      subroutine addconcentrationfluxfromantitrappingmultiorderp(
     &   ifirst0, ilast0, ifirst1, ilast1,
     &   dx,
     &   phase, ngp, norderp,
     &   cl, ca, ngc,
     &   dphidt, ngd,
     &   alpha,
     &   flux0, flux1, ngflux )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1
      integer ngflux, ngd, ngp, ngc, norderp
      double precision
     &     flux0(SIDE2d0(ifirst,ilast,ngflux)),
     &     flux1(SIDE2d1(ifirst,ilast,ngflux))
      double precision phase(CELL2d(ifirst,ilast,ngp),norderp)
      double precision cl(CELL2d(ifirst,ilast,ngc))
      double precision ca(CELL2d(ifirst,ilast,ngc))
      double precision dphidt(CELL2d(ifirst,ilast,ngd),norderp)
      double precision dx(0:1)
      double precision alpha
c
      double precision dxinv, dyinv
      double precision dphix, dphiy, dphi2, dphin
      double precision dphipx, dphipy, dphipn, dphip2
      integer          ic0, ic1, ip
      double precision tol, tol2, dpdt
      double precision ac(2)

      tol = 1.e-8
      tol2 = tol*tol

      dxinv = 1.d0 / dx(0)
      dyinv = 1.d0 / dx(1)

c x-faces
      do ic1 = ifirst1, ilast1
         do ic0 = ifirst0, ilast0+1
c compute gradient liquid phase first
            dphix = dxinv * ( phase(ic0,ic1,norderp)
     &                      - phase(ic0-1,ic1,norderp) )
            dphiy = dyinv*0.25d0*(phase(ic0-1,ic1+1,norderp)
     &                           -phase(ic0-1,ic1-1,norderp)
     &                           +phase(ic0  ,ic1+1,norderp)
     &                           -phase(ic0  ,ic1-1,norderp)
     &                          )
            dphi2 = dphix*dphix+dphiy*dphiy
            if( dphi2 .gt. tol2 ) then
               dphin = sqrt(dphi2)

c average compositions to get values on x-face
               ac(1) = 0.5d0*(cl(ic0-1,ic1)+cl(ic0,ic1))
               ac(2) = 0.5d0*(ca(ic0-1,ic1)+ca(ic0,ic1))

c loop over solid order parameters
               dphipx = 0.d0
               dphipy = 0.d0
               dpdt   = 0.d0
               do ip = 1, norderp-1
                  dphipx = dphipx + dxinv * ( phase(ic0,ic1,ip)
     &                                      - phase(ic0-1,ic1,ip) )
                  dphipy = dphipy + dyinv*0.25d0
     &                   * (phase(ic0-1,ic1+1,ip)
     &                     -phase(ic0-1,ic1-1,ip)
     &                     +phase(ic0  ,ic1+1,ip)
     &                     -phase(ic0  ,ic1-1,ip))
                  dpdt = dpdt + 0.5d0*( dphidt(ic0-1,ic1,ip)
     &                                 +dphidt(ic0,ic1,ip) )
               enddo
               dphip2 = dphipx*dphipx+dphipy*dphipy

               if( dphip2 .gt. tol2 ) then
                  dphipn = sqrt(dphip2)

                  flux0(ic0,ic1) = flux0(ic0,ic1) -
     &               alpha*(dphix/dphin)*(ac(1)-ac(2))*dpdt
               endif
            endif
         enddo
      enddo

      do ic1 = ifirst1, ilast1+1
         do ic0 = ifirst0, ilast0
c compute gradient liquid phase first
            dphiy = dyinv * ( phase(ic0,ic1,norderp)
     &                      - phase(ic0,ic1-1,norderp) )
            dphix = dxinv*0.25d0*(phase(ic0+1,ic1-1,norderp)
     &                           -phase(ic0-1,ic1-1,norderp)
     &                           +phase(ic0+1,ic1  ,norderp)
     &                           -phase(ic0-1,ic1  ,norderp)
     &                           )
            dphi2 = dphix*dphix+dphiy*dphiy
            if( dphi2 .gt. tol2 ) then
               dphin = sqrt(dphi2)

c average compositions to get values on y-face
               ac(1) = 0.5d0*(cl(ic0,ic1-1)+cl(ic0,ic1))
               ac(2) = 0.5d0*(ca(ic0,ic1-1)+ca(ic0,ic1))

c loop over solid order parameters
               dphipx = 0.d0
               dphipy = 0.d0
               dpdt   = 0.d0
               do ip = 1, norderp-1
                  dphipy = dphipy + dyinv * ( phase(ic0,ic1,ip)
     &                                      - phase(ic0,ic1-1,ip) )
                  dphipx = dphipx + dxinv*0.25d0
     &                   * (phase(ic0+1,ic1-1,ip)
     &                     -phase(ic0-1,ic1-1,ip)
     &                     +phase(ic0+1,ic1  ,ip)
     &                     -phase(ic0-1,ic1  ,ip))
                  dpdt =  dpdt + 0.5d0*( dphidt(ic0,ic1-1,ip)
     &                           +dphidt(ic0,ic1,ip) )
               enddo
               dphip2 = dphipx*dphipx+dphipy*dphipy

               if( dphip2 .gt. tol2 ) then
                  dphipn = sqrt(dphip2)

                  flux1(ic0,ic1) = flux1(ic0,ic1) -
     &               alpha*(dphiy/dphin)*(ac(1)-ac(2))*dpdt
               endif
            endif
         enddo
      enddo

      return
      end

