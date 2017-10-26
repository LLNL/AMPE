define(NDIM,3)dnl
define(REAL,`double precision')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine phasesetexactandrhs3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  exact,rhs,dx,xlower,
     &  m,delta,omega,gamma)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
c variables in 1d axis indexed
c
      REAL
     &     dx(0:NDIM-1),
     &     xlower(0:NDIM-1)
c variables in 3d cell indexed
      REAL
     &     exact(CELL3d(ifirst,ilast,1)),
     &     rhs(CELL3d(ifirst,ilast,0))
      REAL m,delta,omega,gamma
c
c***********************************************************************
c
      integer ic0,ic1,ic2
      REAL x, y, tanhx, inv2d, factor1, factor2

      inv2d=0.5/delta
c r.h.s. resulting from diffusion term
      factor1=32.*gamma*m*omega
c r.h.s. resulting from "C*u" term
      factor2=32.*gamma*m*omega

      do ic2=ifirst2,ilast2
        do ic1=ifirst1,ilast1
          do ic0=ifirst0,ilast0
            x = xlower(0) + dx(0)*(ic0-ifirst0+0.5)
            tanhx=0.5*(1.+tanh(inv2d*x))
            exact(ic0,ic1,ic2) = tanhx
            rhs(ic0,ic1,ic2) = -factor1*tanhx*(1.-tanhx)*(1.-2.*tanhx)
     &                 + tanhx
     &                 + factor2*(1.+6.*tanhx*(tanhx-1.))*tanhx
          enddo
        enddo
      enddo

      return
      end
c

