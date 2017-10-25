define(NDIM,3)dnl
define(REAL,`double precision')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
      subroutine multiplyoffdiagbym(
     &  lo0, hi0, lo1, hi1, lo2, hi2,
     &  offdiagx, offdiagy, offdiagz, ngoffdiag,
     &  sqrtm, ngsqrtm)

      implicit none

      integer
     &   lo0, hi0, lo1, hi1, lo2, hi2,
     &   ngoffdiag, ngsqrtm

      double precision offdiagx(SIDE3d0(lo,hi,ngoffdiag))
      double precision offdiagy(SIDE3d1(lo,hi,ngoffdiag))
      double precision offdiagz(SIDE3d2(lo,hi,ngoffdiag))

      double precision sqrtm(CELL3d(lo,hi,ngsqrtm))
c local variables:
      integer i, j, k
c x component
      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0+1
               offdiagx(i,j,k) = offdiagx(i,j,k)
     &                         *sqrtm(i-1,j,k)*sqrtm(i,j,k)
            enddo
         enddo
      enddo

c y component
      do k = lo2, hi2
         do j = lo1, hi1+1
            do i = lo0, hi0
               offdiagy(i,j,k) = offdiagy(i,j,k)
     &                         *sqrtm(i,j-1,k)*sqrtm(i,j,k)
            enddo
         enddo
      enddo

c z component
      do k = lo2, hi2+1
         do j = lo1, hi1
            do i = lo0, hi0
               offdiagz(i,j,k) = offdiagz(i,j,k)
     &                         *sqrtm(i,j,k-1)*sqrtm(i,j,k)
            enddo
         enddo
      enddo

      return
      end
c


c Adapted from SAMRAI test suite

      subroutine setexactandrhs3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  exact,rhs,dx,xlower)
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
c
c***********************************************************************     
c
      integer ic0,ic1,ic2
      REAL x, y, z, sinsin, pi

      pi=3.141592654

      do ic2=ifirst2,ilast2
         z = xlower(2) + dx(2)*(ic2-ifirst2+0.5)
         do ic1=ifirst1,ilast1
            y = xlower(1) + dx(1)*(ic1-ifirst1+0.5)
            do ic0=ifirst0,ilast0
               x = xlower(0) + dx(0)*(ic0-ifirst0+0.5)
               sinsin = sin(pi*x) * sin(pi*y) * sin(pi*z)
               exact(ic0,ic1,ic2) = 1.+sinsin
               rhs(ic0,ic1,ic2) = 5.*exact(ic0,ic1,ic2)
     &                            +NDIM*pi*pi*sinsin
            enddo
         enddo
      enddo

      return
      end   
c
