define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim2d.i)dnl
      subroutine multiplyoffdiagbym(
     &  lo0, hi0, lo1, hi1,
     &  offdiagx, offdiagy, ngoffdiag,
     &  sqrtm, ngsqrtm)

      implicit none

      integer
     &   lo0, hi0, lo1, hi1,
     &   ngoffdiag, ngsqrtm

      double precision offdiagx(SIDE2d0(lo,hi,ngoffdiag))
      double precision offdiagy(SIDE2d1(lo,hi,ngoffdiag))
      double precision sqrtm(CELL2d(lo,hi,ngsqrtm))
c local variables:
      integer i, j
c x component
      do j = lo1, hi1
         do i = lo0, hi0+1
            offdiagx(i,j) = offdiagx(i,j)*sqrtm(i-1,j)*sqrtm(i,j)
         enddo
      enddo

c y component
      do j = lo1, hi1+1
         do i = lo0, hi0
            offdiagy(i,j) = offdiagy(i,j)*sqrtm(i,j-1)*sqrtm(i,j)
         enddo
      enddo

      return
      end
c

c Adapted from SAMRAI test suite

      subroutine setexactandrhs2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  exact,rhs,dx,xlower)
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
c
c***********************************************************************     
c
      integer ic0,ic1
      REAL x, y, sinsin, pi

      pi=3.141592654

      do ic1=ifirst1,ilast1
        y = xlower(1) + dx(1)*(ic1-ifirst1+0.5)
        do ic0=ifirst0,ilast0
          x = xlower(0) + dx(0)*(ic0-ifirst0+0.5)
          sinsin = sin(pi*x) * sin(pi*y)
          exact(ic0,ic1) = 1.+sinsin
          rhs(ic0,ic1) = 5.*exact(ic0,ic1)+NDIM*pi*pi*sinsin
        enddo
      enddo

      return
      end   
c

