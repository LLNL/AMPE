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
      REAL x, tanhx, inv2d, factor1, factor2

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
      subroutine projectphi3d(
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     depth,
     &     phi, ng,
     &     corr, ngc,
     &     err, nge
     &     )
c
      implicit none
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer depth, ng, ngc, nge

      double precision phi(CELL3d(ifirst,ilast,ng),depth)
      double precision corr(CELL3d(ifirst,ilast,ngc),depth)
      double precision err(CELL3d(ifirst,ilast,nge),depth)

c     local variables:
      double precision fac, avg
      integer i, j, k, m

      avg = 1.d0/depth

      do k = ifirst2, ilast2
         do j = ifirst1, ilast1
            do i = ifirst0, ilast0

c           Projection into plane sum_m phi_m=1
c           consists in adding a correction in the direction (1,1,1)
               fac = 0.d0
               do m = 1, depth
                  fac = fac + phi(i,j,k,m)
               enddo
               fac = (fac - 1.d0) * avg

c              Store the projection in the correction array for now
               do m = 1, depth
                  corr(i,j,k,m) = phi(i,j,k,m) - fac
               enddo
c
               fac = 0.d0
               do m = 1, depth
                  fac = fac + err(i,j,k,m)
               enddo
               fac = fac * avg

c              Subtract the error component in the (1,1,1) direction
               do m = 1, depth
                  err(i,j,k,m) = err(i,j,k,m) - fac
               enddo

c              Finalize the correction: phi + corr is on the constraint
               do m = 1, depth
                  corr(i,j,k,m) = corr(i,j,k,m) - phi(i,j,k,m)
               enddo

            enddo
         enddo
      enddo

      return
      end
