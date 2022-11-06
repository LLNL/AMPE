define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine phasesetexactandrhs2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  exact,rhs,dx,xlower,
     &  m,delta,omega,gamma)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
c variables in 1d axis indexed
c
      REAL
     &     dx(0:NDIM-1),
     &     xlower(0:NDIM-1)
c variables in 2d cell indexed
      REAL
     &     exact(CELL2d(ifirst,ilast,1)),
     &     rhs(CELL2d(ifirst,ilast,0))
      REAL m,delta,omega,gamma
c
c***********************************************************************
c
      integer ic0,ic1
      REAL x, y, tanhx, inv2d, factor1, factor2

      inv2d=0.5/delta
c r.h.s. resulting from diffusion term
      factor1=32.*gamma*m*omega
c r.h.s. resulting from "C*u" term
      factor2=32.*gamma*m*omega
c      print*,gamma,m,omega,delta
      do ic1=ifirst1,ilast1
        do ic0=ifirst0,ilast0
          x = xlower(0) + dx(0)*(ic0-ifirst0+0.5)
          tanhx=0.5*(1.+tanh(inv2d*x))
          exact(ic0,ic1) = tanhx
          rhs(ic0,ic1) = -factor1*tanhx*(1.-tanhx)*(1.-2.*tanhx)
     &                 + tanhx
     &                 + factor2*(1.+6.*tanhx*(tanhx-1.))*tanhx
        enddo
      enddo

      return
      end
c
      subroutine projectphi2d(
     &     lo0, hi0, lo1, hi1,
     &     depth,
     &     phi, plo0, phi0, plo1, phi1,
     &     corr, clo0, chi0, clo1, chi1,
     &     err, elo0, ehi0, elo1, ehi1
     &     )
c
      implicit none
      integer lo0, hi0, lo1, hi1,
     &        depth,
     &        plo0, phi0, plo1, phi1,
     &        clo0, chi0, clo1, chi1,
     &        elo0, ehi0, elo1, ehi1
      double precision
     &        phi(plo0:phi0,plo1:phi1,depth),
     &        corr(clo0:chi0,clo1:chi1,depth),
     &        err(elo0:ehi0,elo1:ehi1,depth)

c     local variables:
      double precision fac
      integer i, j, m

      do j = lo1, hi1
         do i = lo0, hi0

c           Projection into plane sum_m phi_m=1
c           consists in adding a correction in the direction (1,1,1)
            fac = 0.d0
            do m = 1, depth
               fac = fac + phi(i,j,m)
            enddo
            fac = (fac - 1.d0)/3.d0

c           Store the projection in the correction array for now
            do m = 1, depth
               corr(i,j,m) = phi(i,j,m) - fac
            enddo

c
            fac = 0.d0
            do m = 1, depth
               fac = fac + err(i,j,m)
            enddo
            fac = fac / depth;

c           Subtract the error component in the (1,1,1) direction
            do m = 1, depth
               err(i,j,m) = err(i,j,m) - fac
            enddo

c           Finalize the correction: phi + corr is on the constraint
            do m = 1, depth
               corr(i,j,m) = corr(i,j,m) - phi(i,j,m)
            enddo

         enddo
      enddo

      return
      end

