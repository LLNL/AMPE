c
c  File:        $URL$
c  Package:     SAMRAI geometry
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: FORTRAN routines for spatial refining of 2d patch data
c               on a regular Cartesian mesh.
c
c
c  File:        $URL$
c  Package:     SAMRAI geometry
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for 2d Cartesian refine operators
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for dimensioning 2d arrays in FORTRAN routines.
c
c
c  File:        $URL$
c  Package:     SAMRAI geometry
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for cartesian geometry transfer routines.
c




c
c
c
c
c
c
c
c
c
c
c
c***********************************************************************
c Linear interpolation for 2d cell-centered double data
c***********************************************************************
c
      subroutine cartlinrefcelldoub2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1)
      double precision
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1)
      double precision deltax(0:15,0:2-1),x,y
      integer ic0,ic1,if0,if1,ir0,ir1
c
c***********************************************************************
c
      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         ir1=if1-ic1*ratio(1)
         y=deltax(ir1,1)/dxc(1)
         if( y .lt. 0.d0 ) then
            ic1 = ic1-1
            y = y + one
         endif
         do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            ir0=if0-ic0*ratio(0)
            x=deltax(ir0,0)/dxc(0)
            if( x .lt. 0.d0 ) then
               ic0 = ic0-1
               x = x + one
            endif
            arrayf(if0,if1)=
     &      (arrayc(ic0,ic1)+(arrayc(ic0+1,ic1)-arrayc(ic0,ic1))*x)
     &        *(one-y)
     &      +(arrayc(ic0,ic1+1)
     &        +(arrayc(ic0+1,ic1+1)-arrayc(ic0,ic1+1))*x)*y
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear interpolation for 2d cell-centered float data
c***********************************************************************
c
      subroutine cartlinrefcellflot2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1)
      real
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1)
      double precision deltax(0:15,0:2-1),x,y
      integer ic0,ic1,if0,if1,ir0,ir1
c
c***********************************************************************
c
      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         ir1=if1-ic1*ratio(1)
         y=deltax(ir1,1)/dxc(1)
         if( y .lt. 0.d0 ) then
            ic1 = ic1-1
            y = y + one
         endif
         do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            ir0=if0-ic0*ratio(0)
            x=deltax(ir0,0)/dxc(0)
            if( x .lt. 0.d0 ) then
               ic0 = ic0-1
               x = x + one
            endif
            arrayf(if0,if1)=
     &      (arrayc(ic0,ic1)+(arrayc(ic0+1,ic1)-arrayc(ic0,ic1))*x)
     &        *(one-y)
     &      +(arrayc(ic0,ic1+1)
     &        +(arrayc(ic0+1,ic1+1)-arrayc(ic0,ic1+1))*x)*y
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear interpolation for 2d cell-centered complex data
c***********************************************************************
c
      subroutine cartlinrefcellcplx2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1)
      double complex
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1)
      double precision deltax(0:15,0:2-1),x,y
      integer ic0,ic1,if0,if1,ir0,ir1
c
c***********************************************************************
c
      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         ir1=if1-ic1*ratio(1)
         y=deltax(ir1,1)/dxc(1)
         if( y .lt. 0.d0 ) then
            ic1 = ic1-1
            y = y + one
         endif
         do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            ir0=if0-ic0*ratio(0)
            x=deltax(ir0,0)/dxc(0)
            if( x .lt. 0.d0 ) then
               ic0 = ic0-1
               x = x + one
            endif
            arrayf(if0,if1)=
     &      (arrayc(ic0,ic1)+(arrayc(ic0+1,ic1)-arrayc(ic0,ic1))*x)
     &        *(one-y)
     &      +(arrayc(ic0,ic1+1)
     &        +(arrayc(ic0+1,ic1+1)-arrayc(ic0,ic1+1))*x)*y
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 2d cell-centered double data
c***********************************************************************
c
      subroutine cartclinrefcelldoub2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1),
     &  deltax(0:15,0:2-1)
      double precision
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0,
     &          cilo1:cihi1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0,
     &          cilo1:cihi1)
      integer ic0,ic1,ie0,ie1,if0,if1,ir0,ir1
      double precision
     &  coef2,bound
      double precision deltax1
c
c***********************************************************************
c
      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ic1=ifirstc1,ilastc1
      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ie0,ic1)
     &               -arrayc(ie0-1,ic1)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic0,ic1)=zero
         endif
      enddo
      enddo

      do ic0=ifirstc0,ilastc0
      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ic0,ie1)
     &               -arrayc(ic0,ie1-1)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic0,ic1)=zero
         endif
      enddo
      enddo

      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         ir1=if1-ic1*ratio(1)
         deltax1=deltax(ir1,1)
         do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            ir0=if0-ic0*ratio(0)
            arrayf(if0,if1)=arrayc(ic0,ic1)
     &                      +slope0(ic0,ic1)*deltax(ir0,0)
     &                      +slope1(ic0,ic1)*deltax1
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 2d cell-centered float data
c***********************************************************************
c
      subroutine cartclinrefcellflot2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1),
     &  deltax(0:15,0:2-1)
      real
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0,
     &          cilo1:cihi1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0,
     &          cilo1:cihi1)
      integer ic0,ic1,ie0,ie1,if0,if1,ir0,ir1
      real
     &  coef2,bound
      double precision deltax1
c
c***********************************************************************
c
      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ic1=ifirstc1,ilastc1
      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ie0,ic1)
     &               -arrayc(ie0-1,ic1)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic0,ic1)=zero
         endif
      enddo
      enddo

      do ic0=ifirstc0,ilastc0
      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ic0,ie1)
     &               -arrayc(ic0,ie1-1)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic0,ic1)=zero
         endif
      enddo
      enddo

      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         ir1=if1-ic1*ratio(1)
         deltax1=deltax(ir1,1)
         do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            ir0=if0-ic0*ratio(0)
            arrayf(if0,if1)=arrayc(ic0,ic1)
     &                      +slope0(ic0,ic1)*deltax(ir0,0)
     &                      +slope1(ic0,ic1)*deltax1
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 2d cell-centered complex data
c***********************************************************************
c
      subroutine cartclinrefcellcplx2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1),
     &  deltax(0:15,0:2-1)
      double complex
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0,
     &          cilo1:cihi1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0,
     &          cilo1:cihi1)
      integer ic0,ic1,ie0,ie1,if0,if1,ir0,ir1
      double precision
     &  coef2real,coef2imag,boundreal,boundimag,
     &  diff0real,diff0imag,diff1real,diff1imag,
     &  slopereal,slopeimag
      double precision deltax1
c
c***********************************************************************
c
      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ic1=ifirstc1,ilastc1
      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ie0,ic1)
     &               -arrayc(ie0-1,ic1)
      enddo
      do ic0=ifirstc0,ilastc0
         diff0real=dble(diff0(ic0))
         diff0imag=imag(diff0(ic0))
         diff1real=dble(diff0(ic0+1))
         diff1imag=imag(diff0(ic0+1))
         coef2real=half*(diff0real+diff1real)
         coef2imag=half*(diff0imag+diff1imag)
         boundreal=two*min(abs(diff1real),abs(diff0real))
         boundimag=two*min(abs(diff1imag),abs(diff0imag))
         if (diff0real*diff1real.gt.zero) then
            slopereal=sign(min(abs(coef2real),boundreal),coef2real)
     &                /dxc(0)
         else
            slopereal=zero
         endif
         if (diff0imag*diff1imag.gt.zero) then
            slopeimag=sign(min(abs(coef2imag),boundimag),coef2imag) 
     &                /dxc(0)
         else
            slopeimag=zero
         endif
         slope0(ic0,ic1) = slopereal + cmplx(zero,one)*slopeimag
      enddo
      enddo

      do ic0=ifirstc0,ilastc0
      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ic0,ie1)
     &               -arrayc(ic0,ie1-1)
      enddo
      do ic1=ifirstc1,ilastc1
         diff0real=dble(diff1(ic1))
         diff0imag=imag(diff1(ic1))
         diff1real=dble(diff1(ic1+1))
         diff1imag=imag(diff1(ic1+1))
         coef2real=half*(diff0real+diff1real)
         coef2imag=half*(diff0imag+diff1imag)
         boundreal=two*min(abs(diff1real),abs(diff0real))
         boundimag=two*min(abs(diff1imag),abs(diff0imag))
         if (diff0real*diff1real.gt.zero) then
            slopereal=sign(min(abs(coef2real),boundreal),coef2real)
     &                /dxc(1)
         else
            slopereal=zero
         endif
         if (diff0imag*diff1imag.gt.zero) then
            slopeimag=sign(min(abs(coef2imag),boundimag),coef2imag) 
     &                /dxc(1)
         else
            slopeimag=zero
         endif
         slope1(ic0,ic1) = slopereal + cmplx(zero,one)*slopeimag
      enddo
      enddo

      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         ir1=if1-ic1*ratio(1)
         deltax1=deltax(ir1,1)
         do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            ir0=if0-ic0*ratio(0)
            arrayf(if0,if1)=arrayc(ic0,ic1)
     &                      +slope0(ic0,ic1)*deltax(ir0,0)
     &                      +slope1(ic0,ic1)*deltax1
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 2d edge-centered double data
c***********************************************************************
c
      subroutine cartclinrefedgedoub2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1),
     &  deltax(0:15,0:2-1)
      double precision
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0,
     &          cilo1:cihi1+1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0,
     &          cilo1:cihi1+1)
      integer ic0,ic1,ie0,ie1,if0,if1,ir0,ir1
      double precision
     &  coef2,bound
      double precision deltax1
c
c***********************************************************************
c


      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo


      do ir1=0,ratio(1)-1
         deltax(ir1,1)=dble(ir1)*dxf(1)
      enddo

      do ic1=ifirstc1,ilastc1+1

      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ie0,ic1)
     &               -arrayc(ie0-1,ic1)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic0,ic1)=zero
         endif
      enddo
      enddo

      do ic0=ifirstc0,ilastc0+0

      do ic1=ifirstc1-1,ilastc1+1
         diff1(ic1)=arrayc(ic0,ic1+1)
     &                -arrayc(ic0,ic1)
      enddo
      do ie1=ifirstc1,ilastc1+1
         coef2=half*(diff1(ie1-1)+diff1(ie1))
         bound=two*min(abs(diff1(ie1-1)),abs(diff1(ie1)))
         if (diff1(ie1)*diff1(ie1-1).gt.zero) then
           slope1(ic0,ie1)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(1) 
         else
            slope1(ic0,ie1)=zero
         endif
      enddo
      enddo

      do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         ir1=if1-ic1*ratio(1)
         deltax1=deltax(ir1,1)
         do if0=ifirstf0,ilastf0+0
         if (if0.lt.0) then
            ie0=(if0+1)/ratio(0)-1
         else
            ie0=if0/ratio(0)
         endif
            ir0=if0-ie0*ratio(0)
            arrayf(if0,if1)=arrayc(ie0,ic1)
     &                +slope0(ie0,ic1)*deltax(ir0,0)
     &                +slope1(ie0,ic1)*deltax1
         enddo
      enddo
c
      return
      end
c
      subroutine cartclinrefedgedoub2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff1,slope1,diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1),
     &  deltax(0:15,0:2-1)
      double precision
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0+1,
     &          cilo1:cihi1)
      integer ic0,ic1,ie0,ie1,if0,if1,ir0,ir1
      double precision
     &  coef2,bound
      double precision deltax1
c
c***********************************************************************
c


      do ir0=0,ratio(0)-1
         deltax(ir0,0)=dble(ir0)*dxf(0)
      enddo


      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ic1=ifirstc1,ilastc1+0

      do ic0=ifirstc0-1,ilastc0+1
         diff0(ic0)=arrayc(ic0+1,ic1)
     &                -arrayc(ic0,ic1)
      enddo
      do ie0=ifirstc0,ilastc0+1
         coef2=half*(diff0(ie0-1)+diff0(ie0))
         bound=two*min(abs(diff0(ie0-1)),abs(diff0(ie0)))
         if (diff0(ie0)*diff0(ie0-1).gt.zero) then
           slope0(ie0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(0) 
         else
            slope0(ie0,ic1)=zero
         endif
      enddo
      enddo

      do ic0=ifirstc0,ilastc0+1

      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ic0,ie1)
     &               -arrayc(ic0,ie1-1)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic0,ic1)=zero
         endif
      enddo
      enddo

      do if1=ifirstf1,ilastf1+0
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         ir1=if1-ic1*ratio(1)
         deltax1=deltax(ir1,1)
         do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ie0=(if0+1)/ratio(0)-1
         else
            ie0=if0/ratio(0)
         endif
            ir0=if0-ie0*ratio(0)
            arrayf(if0,if1)=arrayc(ie0,ic1)
     &                +slope0(ie0,ic1)*deltax(ir0,0)
     &                +slope1(ie0,ic1)*deltax1
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 2d edge-centered float data
c***********************************************************************
c
      subroutine cartclinrefedgeflot2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1),
     &  deltax(0:15,0:2-1)
      real
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0,
     &          cilo1:cihi1+1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0,
     &          cilo1:cihi1+1)
      integer ic0,ic1,ie0,ie1,if0,if1,ir0,ir1
      real
     &  coef2,bound
      double precision deltax1
c
c***********************************************************************
c


      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo


      do ir1=0,ratio(1)-1
         deltax(ir1,1)=dble(ir1)*dxf(1)
      enddo

      do ic1=ifirstc1,ilastc1+1

      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ie0,ic1)
     &               -arrayc(ie0-1,ic1)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic0,ic1)=zero
         endif
      enddo
      enddo

      do ic0=ifirstc0,ilastc0+0

      do ic1=ifirstc1-1,ilastc1+1
         diff1(ic1)=arrayc(ic0,ic1+1)
     &                -arrayc(ic0,ic1)
      enddo
      do ie1=ifirstc1,ilastc1+1
         coef2=half*(diff1(ie1-1)+diff1(ie1))
         bound=two*min(abs(diff1(ie1-1)),abs(diff1(ie1)))
         if (diff1(ie1)*diff1(ie1-1).gt.zero) then
           slope1(ic0,ie1)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(1) 
         else
            slope1(ic0,ie1)=zero
         endif
      enddo
      enddo

      do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         ir1=if1-ic1*ratio(1)
         deltax1=deltax(ir1,1)
         do if0=ifirstf0,ilastf0+0
         if (if0.lt.0) then
            ie0=(if0+1)/ratio(0)-1
         else
            ie0=if0/ratio(0)
         endif
            ir0=if0-ie0*ratio(0)
            arrayf(if0,if1)=arrayc(ie0,ic1)
     &                +slope0(ie0,ic1)*deltax(ir0,0)
     &                +slope1(ie0,ic1)*deltax1
         enddo
      enddo
c
      return
      end
c
      subroutine cartclinrefedgeflot2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff1,slope1,diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1),
     &  deltax(0:15,0:2-1)
      real
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0+1,
     &          cilo1:cihi1)
      integer ic0,ic1,ie0,ie1,if0,if1,ir0,ir1
      real
     &  coef2,bound
      double precision deltax1
c
c***********************************************************************
c


      do ir0=0,ratio(0)-1
         deltax(ir0,0)=dble(ir0)*dxf(0)
      enddo


      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ic1=ifirstc1,ilastc1+0

      do ic0=ifirstc0-1,ilastc0+1
         diff0(ic0)=arrayc(ic0+1,ic1)
     &                -arrayc(ic0,ic1)
      enddo
      do ie0=ifirstc0,ilastc0+1
         coef2=half*(diff0(ie0-1)+diff0(ie0))
         bound=two*min(abs(diff0(ie0-1)),abs(diff0(ie0)))
         if (diff0(ie0)*diff0(ie0-1).gt.zero) then
           slope0(ie0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(0) 
         else
            slope0(ie0,ic1)=zero
         endif
      enddo
      enddo

      do ic0=ifirstc0,ilastc0+1

      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ic0,ie1)
     &               -arrayc(ic0,ie1-1)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic0,ic1)=zero
         endif
      enddo
      enddo

      do if1=ifirstf1,ilastf1+0
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         ir1=if1-ic1*ratio(1)
         deltax1=deltax(ir1,1)
         do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ie0=(if0+1)/ratio(0)-1
         else
            ie0=if0/ratio(0)
         endif
            ir0=if0-ie0*ratio(0)
            arrayf(if0,if1)=arrayc(ie0,ic1)
     &                +slope0(ie0,ic1)*deltax(ir0,0)
     &                +slope1(ie0,ic1)*deltax1
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 2d face-centered double data
c***********************************************************************
c
      subroutine cartclinreffacedoub2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1),
     &  deltax(0:15,0:2-1)
      double precision
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0+1,
     &          cilo1:cihi1)
      integer ic0,ic1,ie0,ie1,if0,if1,ir0,ir1
      double precision
     &  coef2,bound
      double precision deltax1
c
c***********************************************************************
c

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ir0=0,ratio(0)-1
         deltax(ir0,0)=dble(ir0)*dxf(0)
      enddo

      do ic1=ifirstc1,ilastc1
      do ic0=ifirstc0-1,ilastc0+1
         diff0(ic0)=arrayc(ic0+1,ic1)
     &                -arrayc(ic0,ic1)
      enddo
      do ie0=ifirstc0,ilastc0+1
         coef2=half*(diff0(ie0-1)+diff0(ie0))
         bound=two*min(abs(diff0(ie0-1)),abs(diff0(ie0)))
         if (diff0(ie0)*diff0(ie0-1).gt.zero) then
           slope0(ie0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(0) 
         else
            slope0(ie0,ic1)=zero
         endif
      enddo
      enddo

      do ie0=ifirstc0,ilastc0+1
      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ie0,ie1)
     &               -arrayc(ie0,ie1-1)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ie0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ie0,ic1)=zero
         endif
      enddo
      enddo

      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         ir1=if1-ic1*ratio(1)
         deltax1=deltax(ir1,1)
         do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ie0=(if0+1)/ratio(0)-1
         else
            ie0=if0/ratio(0)
         endif
            ir0=if0-ie0*ratio(0)
            arrayf(if0,if1)=arrayc(ie0,ic1)
     &                +slope0(ie0,ic1)*deltax(ir0,0)
     &                +slope1(ie0,ic1)*deltax1
         enddo
      enddo
c
      return
      end
c
      subroutine cartclinreffacedoub2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff1,slope1,diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1),
     &  deltax(0:15,0:2-1)
      double precision
     &  arrayc(cilo1:cihi1+1,
     &          cilo0:cihi0),
     &  arrayf(filo1:fihi1+1,
     &          filo0:fihi0),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo1:cihi1+1,
     &          cilo0:cihi0),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo1:cihi1+1,
     &          cilo0:cihi0)
      integer ic1,ic0,ie1,ie0,if1,if0,ir1,ir0
      double precision
     &  coef2,bound
      double precision deltax0
c
c***********************************************************************
c

      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=dble(ir1)*dxf(1)
      enddo

      do ic0=ifirstc0,ilastc0
      do ic1=ifirstc1-1,ilastc1+1
         diff1(ic1)=arrayc(ic1+1,ic0)
     &                -arrayc(ic1,ic0)
      enddo
      do ie1=ifirstc1,ilastc1+1
         coef2=half*(diff1(ie1-1)+diff1(ie1))
         bound=two*min(abs(diff1(ie1-1)),abs(diff1(ie1)))
         if (diff1(ie1)*diff1(ie1-1).gt.zero) then
           slope1(ie1,ic0)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(1) 
         else
            slope1(ie1,ic0)=zero
         endif
      enddo
      enddo

      do ie1=ifirstc1,ilastc1+1
      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ie1,ie0)
     &               -arrayc(ie1,ie0-1)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ie1,ic0)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ie1,ic0)=zero
         endif
      enddo
      enddo

      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         ir0=if0-ic0*ratio(0)
         deltax0=deltax(ir0,0)
         do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ie1=(if1+1)/ratio(1)-1
         else
            ie1=if1/ratio(1)
         endif
            ir1=if1-ie1*ratio(1)
            arrayf(if1,if0)=arrayc(ie1,ic0)
     &                +slope1(ie1,ic0)*deltax(ir1,1)
     &                +slope0(ie1,ic0)*deltax0
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 2d face-centered float data
c***********************************************************************
c
      subroutine cartclinreffaceflot2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1),
     &  deltax(0:15,0:2-1)
      real
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0+1,
     &          cilo1:cihi1)
      integer ic0,ic1,ie0,ie1,if0,if1,ir0,ir1
      real
     &  coef2,bound
      double precision deltax1
c
c***********************************************************************
c

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ir0=0,ratio(0)-1
         deltax(ir0,0)=dble(ir0)*dxf(0)
      enddo

      do ic1=ifirstc1,ilastc1
      do ic0=ifirstc0-1,ilastc0+1
         diff0(ic0)=arrayc(ic0+1,ic1)
     &                -arrayc(ic0,ic1)
      enddo
      do ie0=ifirstc0,ilastc0+1
         coef2=half*(diff0(ie0-1)+diff0(ie0))
         bound=two*min(abs(diff0(ie0-1)),abs(diff0(ie0)))
         if (diff0(ie0)*diff0(ie0-1).gt.zero) then
           slope0(ie0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(0) 
         else
            slope0(ie0,ic1)=zero
         endif
      enddo
      enddo

      do ie0=ifirstc0,ilastc0+1
      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ie0,ie1)
     &               -arrayc(ie0,ie1-1)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ie0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ie0,ic1)=zero
         endif
      enddo
      enddo

      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         ir1=if1-ic1*ratio(1)
         deltax1=deltax(ir1,1)
         do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ie0=(if0+1)/ratio(0)-1
         else
            ie0=if0/ratio(0)
         endif
            ir0=if0-ie0*ratio(0)
            arrayf(if0,if1)=arrayc(ie0,ic1)
     &                +slope0(ie0,ic1)*deltax(ir0,0)
     &                +slope1(ie0,ic1)*deltax1
         enddo
      enddo
c
      return
      end
c
      subroutine cartclinreffaceflot2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff1,slope1,diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1),
     &  deltax(0:15,0:2-1)
      real
     &  arrayc(cilo1:cihi1+1,
     &          cilo0:cihi0),
     &  arrayf(filo1:fihi1+1,
     &          filo0:fihi0),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo1:cihi1+1,
     &          cilo0:cihi0),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo1:cihi1+1,
     &          cilo0:cihi0)
      integer ic1,ic0,ie1,ie0,if1,if0,ir1,ir0
      real
     &  coef2,bound
      double precision deltax0
c
c***********************************************************************
c

      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=dble(ir1)*dxf(1)
      enddo

      do ic0=ifirstc0,ilastc0
      do ic1=ifirstc1-1,ilastc1+1
         diff1(ic1)=arrayc(ic1+1,ic0)
     &                -arrayc(ic1,ic0)
      enddo
      do ie1=ifirstc1,ilastc1+1
         coef2=half*(diff1(ie1-1)+diff1(ie1))
         bound=two*min(abs(diff1(ie1-1)),abs(diff1(ie1)))
         if (diff1(ie1)*diff1(ie1-1).gt.zero) then
           slope1(ie1,ic0)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(1) 
         else
            slope1(ie1,ic0)=zero
         endif
      enddo
      enddo

      do ie1=ifirstc1,ilastc1+1
      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ie1,ie0)
     &               -arrayc(ie1,ie0-1)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ie1,ic0)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ie1,ic0)=zero
         endif
      enddo
      enddo

      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         ir0=if0-ic0*ratio(0)
         deltax0=deltax(ir0,0)
         do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ie1=(if1+1)/ratio(1)-1
         else
            ie1=if1/ratio(1)
         endif
            ir1=if1-ie1*ratio(1)
            arrayf(if1,if0)=arrayc(ie1,ic0)
     &                +slope1(ie1,ic0)*deltax(ir1,1)
     &                +slope0(ie1,ic0)*deltax0
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 2d face-centered complex data
c***********************************************************************
c
c      subroutine cartclinreffacecplx2d0(
ccart_clinref_op_face_2d(double complex,0,1)cc
c      subroutine cartclinreffacecplx2d1(
ccart_clinref_op_face_2d(double complex,1,0)c
c***********************************************************************
c Linear interpolation for 2d node-centered double data
c***********************************************************************
c
       subroutine cartlinrefnodedoub2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1)
      double precision
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1)
      double precision x,y,realrat0,realrat1
      integer ic0,ic1,if0,if1,ie0,ie1,ir0,ir1,i,j
c
c***********************************************************************
c
      realrat0=one/dble(ratio(0))
      realrat1=one/dble(ratio(1))

      do ic1=ifirstc1,ilastc1
         if1=ic1*ratio(1)
         do ir1=0,ratio(1)
            ie1=if1+ir1
            if ((ie1.ge.filo1).and.(ie1.le.(fihi1+1))) then
      do ic0=ifirstc0,ilastc0
         if0=ic0*ratio(0)
         do ir0=0,ratio(0)
            ie0=if0+ir0
            if ((ie0.ge.filo0).and.(ie0.le.(fihi0+1))) then
               x = dble(ir0)*realrat0
               y = dble(ir1)*realrat1
               arrayf(ie0,ie1)=
     &              (arrayc(ic0,ic1)*(one-x) + 
     &               arrayc(ic0+1,ic1)*x)*(one-y) +
     &              (arrayc(ic0,ic1+1)*(one-x) + 
     &               arrayc(ic0+1,ic1+1)*x)*y
           endif
         end do
      end do
           endif
         end do
      end do
c
      return
      end
c
c***********************************************************************
c Linear interpolation for 2d node-centered float data
c***********************************************************************
c
       subroutine cartlinrefnodeflot2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1)
      real
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1)
      double precision x,y,realrat0,realrat1
      integer ic0,ic1,if0,if1,ie0,ie1,ir0,ir1,i,j
c
c***********************************************************************
c
      realrat0=one/dble(ratio(0))
      realrat1=one/dble(ratio(1))

      do ic1=ifirstc1,ilastc1
         if1=ic1*ratio(1)
         do ir1=0,ratio(1)
            ie1=if1+ir1
            if ((ie1.ge.filo1).and.(ie1.le.(fihi1+1))) then
      do ic0=ifirstc0,ilastc0
         if0=ic0*ratio(0)
         do ir0=0,ratio(0)
            ie0=if0+ir0
            if ((ie0.ge.filo0).and.(ie0.le.(fihi0+1))) then
               x = dble(ir0)*realrat0
               y = dble(ir1)*realrat1
               arrayf(ie0,ie1)=
     &              (arrayc(ic0,ic1)*(one-x) + 
     &               arrayc(ic0+1,ic1)*x)*(one-y) +
     &              (arrayc(ic0,ic1+1)*(one-x) + 
     &               arrayc(ic0+1,ic1+1)*x)*y
           endif
         end do
      end do
           endif
         end do
      end do
c
      return
      end
c
c***********************************************************************
c Linear interpolation for 2d node-centered complex data
c***********************************************************************
c
       subroutine cartlinrefnodecplx2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1)
      double complex
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1)
      double precision x,y,realrat0,realrat1
      integer ic0,ic1,if0,if1,ie0,ie1,ir0,ir1,i,j
c
c***********************************************************************
c
      realrat0=one/dble(ratio(0))
      realrat1=one/dble(ratio(1))

      do ic1=ifirstc1,ilastc1
         if1=ic1*ratio(1)
         do ir1=0,ratio(1)
            ie1=if1+ir1
            if ((ie1.ge.filo1).and.(ie1.le.(fihi1+1))) then
      do ic0=ifirstc0,ilastc0
         if0=ic0*ratio(0)
         do ir0=0,ratio(0)
            ie0=if0+ir0
            if ((ie0.ge.filo0).and.(ie0.le.(fihi0+1))) then
               x = dble(ir0)*realrat0
               y = dble(ir1)*realrat1
               arrayf(ie0,ie1)=
     &              (arrayc(ic0,ic1)*(one-x) + 
     &               arrayc(ic0+1,ic1)*x)*(one-y) +
     &              (arrayc(ic0,ic1+1)*(one-x) + 
     &               arrayc(ic0+1,ic1+1)*x)*y
           endif
         end do
      end do
           endif
         end do
      end do
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 2d outerface double data
c***********************************************************************
c
      subroutine cartclinrefoutfacedoub2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1),
     &  deltax(0:15,0:2-1)
      double precision
     &  arrayc(cilo1:cihi1),
     &  arrayf(filo1:fihi1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo1:cihi1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo1:cihi1)
      integer ic1,ie1,if1,ir1
      double precision
     &  coef2,bound
c
c***********************************************************************
c

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ie1)
     &               -arrayc(ie1-1)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic1)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic1)=zero
         endif
      enddo

      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         ir1=if1-ic1*ratio(1)
         arrayf(if1)=arrayc(ic1)
     &                +slope1(ic1)*deltax(ir1,1)
      enddo
c
      return
      end
c
      subroutine cartclinrefoutfacedoub2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff1,slope1,diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1),
     &  deltax(0:15,0:2-1)
      double precision
     &  arrayc(cilo0:cihi0),
     &  arrayf(filo0:fihi0),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0)
      integer ic0,ie0,if0,ir0
      double precision
     &  coef2,bound
c
c***********************************************************************
c

      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ie0)
     &               -arrayc(ie0-1)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic0)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic0)=zero
         endif
      enddo

      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         ir0=if0-ic0*ratio(0)
         arrayf(if0)=arrayc(ic0)
     &                +slope0(ic0)*deltax(ir0,0)
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 2d outerface float data
c***********************************************************************
c
      subroutine cartclinrefoutfaceflot2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1),
     &  deltax(0:15,0:2-1)
      real
     &  arrayc(cilo1:cihi1),
     &  arrayf(filo1:fihi1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo1:cihi1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo1:cihi1)
      integer ic1,ie1,if1,ir1
      real
     &  coef2,bound
c
c***********************************************************************
c

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ie1)
     &               -arrayc(ie1-1)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic1)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic1)=zero
         endif
      enddo

      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         ir1=if1-ic1*ratio(1)
         arrayf(if1)=arrayc(ic1)
     &                +slope1(ic1)*deltax(ir1,1)
      enddo
c
      return
      end
c
      subroutine cartclinrefoutfaceflot2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff1,slope1,diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1),
     &  deltax(0:15,0:2-1)
      real
     &  arrayc(cilo0:cihi0),
     &  arrayf(filo0:fihi0),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0)
      integer ic0,ie0,if0,ir0
      real
     &  coef2,bound
c
c***********************************************************************
c

      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ie0)
     &               -arrayc(ie0-1)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic0)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic0)=zero
         endif
      enddo

      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         ir0=if0-ic0*ratio(0)
         arrayf(if0)=arrayc(ic0)
     &                +slope0(ic0)*deltax(ir0,0)
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 2d outerface complex data
c***********************************************************************
c
c      subroutine cartclinrefoutfacecplx2d0(
ccart_clinref_op_outerface_2d(double complex,0,1)cc
c      subroutine cartclinrefoutfacecplx2d1(
ccart_clinref_op_outerface_2d(double complex,1,0)c
c***********************************************************************
c Conservative linear interpolation for 2d side-centered double data
c***********************************************************************
c
      subroutine cartclinrefsidedoub2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1),
     &  deltax(0:15,0:2-1)
      double precision
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0+1,
     &          cilo1:cihi1)
      integer ic0,ic1,ie0,ie1,if0,if1,ir0,ir1
      double precision
     &  coef2,bound
      double precision deltax1
c
c***********************************************************************
c


      do ir0=0,ratio(0)-1
         deltax(ir0,0)=dble(ir0)*dxf(0)
      enddo


      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ic1=ifirstc1,ilastc1+0

      do ic0=ifirstc0-1,ilastc0+1
         diff0(ic0)=arrayc(ic0+1,ic1)
     &                -arrayc(ic0,ic1)
      enddo
      do ie0=ifirstc0,ilastc0+1
         coef2=half*(diff0(ie0-1)+diff0(ie0))
         bound=two*min(abs(diff0(ie0-1)),abs(diff0(ie0)))
         if (diff0(ie0)*diff0(ie0-1).gt.zero) then
           slope0(ie0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(0) 
         else
            slope0(ie0,ic1)=zero
         endif
      enddo
      enddo

      do ic0=ifirstc0,ilastc0+1

      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ic0,ie1)
     &               -arrayc(ic0,ie1-1)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic0,ic1)=zero
         endif
      enddo
      enddo

      do if1=ifirstf1,ilastf1+0
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         ir1=if1-ic1*ratio(1)
         deltax1=deltax(ir1,1)
         do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            ir0=if0-ic0*ratio(0)
            arrayf(if0,if1)=arrayc(ic0,ic1)
     &                +slope0(ic0,ic1)*deltax(ir0,0)
     &                +slope1(ic0,ic1)*deltax1
         enddo
      enddo
c
      return
      end
c
      subroutine cartclinrefsidedoub2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff1,slope1,diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1),
     &  deltax(0:15,0:2-1)
      double precision
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0,
     &          cilo1:cihi1+1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0,
     &          cilo1:cihi1+1)
      integer ic0,ic1,ie0,ie1,if0,if1,ir0,ir1
      double precision
     &  coef2,bound
      double precision deltax1
c
c***********************************************************************
c


      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo


      do ir1=0,ratio(1)-1
         deltax(ir1,1)=dble(ir1)*dxf(1)
      enddo

      do ic1=ifirstc1,ilastc1+1

      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ie0,ic1)
     &               -arrayc(ie0-1,ic1)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic0,ic1)=zero
         endif
      enddo
      enddo

      do ic0=ifirstc0,ilastc0+0

      do ic1=ifirstc1-1,ilastc1+1
         diff1(ic1)=arrayc(ic0,ic1+1)
     &                -arrayc(ic0,ic1)
      enddo
      do ie1=ifirstc1,ilastc1+1
         coef2=half*(diff1(ie1-1)+diff1(ie1))
         bound=two*min(abs(diff1(ie1-1)),abs(diff1(ie1)))
         if (diff1(ie1)*diff1(ie1-1).gt.zero) then
           slope1(ic0,ie1)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(1) 
         else
            slope1(ic0,ie1)=zero
         endif
      enddo
      enddo

      do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         ir1=if1-ic1*ratio(1)
         deltax1=deltax(ir1,1)
         do if0=ifirstf0,ilastf0+0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            ir0=if0-ic0*ratio(0)
            arrayf(if0,if1)=arrayc(ic0,ic1)
     &                +slope0(ic0,ic1)*deltax(ir0,0)
     &                +slope1(ic0,ic1)*deltax1
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 2d side-centered float data
c***********************************************************************
c
      subroutine cartclinrefsideflot2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1),
     &  deltax(0:15,0:2-1)
      real
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0+1,
     &          cilo1:cihi1)
      integer ic0,ic1,ie0,ie1,if0,if1,ir0,ir1
      real
     &  coef2,bound
      double precision deltax1
c
c***********************************************************************
c


      do ir0=0,ratio(0)-1
         deltax(ir0,0)=dble(ir0)*dxf(0)
      enddo


      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ic1=ifirstc1,ilastc1+0

      do ic0=ifirstc0-1,ilastc0+1
         diff0(ic0)=arrayc(ic0+1,ic1)
     &                -arrayc(ic0,ic1)
      enddo
      do ie0=ifirstc0,ilastc0+1
         coef2=half*(diff0(ie0-1)+diff0(ie0))
         bound=two*min(abs(diff0(ie0-1)),abs(diff0(ie0)))
         if (diff0(ie0)*diff0(ie0-1).gt.zero) then
           slope0(ie0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(0) 
         else
            slope0(ie0,ic1)=zero
         endif
      enddo
      enddo

      do ic0=ifirstc0,ilastc0+1

      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ic0,ie1)
     &               -arrayc(ic0,ie1-1)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic0,ic1)=zero
         endif
      enddo
      enddo

      do if1=ifirstf1,ilastf1+0
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         ir1=if1-ic1*ratio(1)
         deltax1=deltax(ir1,1)
         do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            ir0=if0-ic0*ratio(0)
            arrayf(if0,if1)=arrayc(ic0,ic1)
     &                +slope0(ic0,ic1)*deltax(ir0,0)
     &                +slope1(ic0,ic1)*deltax1
         enddo
      enddo
c
      return
      end
c
      subroutine cartclinrefsideflot2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff1,slope1,diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  dxc(0:2-1),
     &  dxf(0:2-1),
     &  deltax(0:15,0:2-1)
      real
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0,
     &          cilo1:cihi1+1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0,
     &          cilo1:cihi1+1)
      integer ic0,ic1,ie0,ie1,if0,if1,ir0,ir1
      real
     &  coef2,bound
      double precision deltax1
c
c***********************************************************************
c


      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo


      do ir1=0,ratio(1)-1
         deltax(ir1,1)=dble(ir1)*dxf(1)
      enddo

      do ic1=ifirstc1,ilastc1+1

      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ie0,ic1)
     &               -arrayc(ie0-1,ic1)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic0,ic1)=zero
         endif
      enddo
      enddo

      do ic0=ifirstc0,ilastc0+0

      do ic1=ifirstc1-1,ilastc1+1
         diff1(ic1)=arrayc(ic0,ic1+1)
     &                -arrayc(ic0,ic1)
      enddo
      do ie1=ifirstc1,ilastc1+1
         coef2=half*(diff1(ie1-1)+diff1(ie1))
         bound=two*min(abs(diff1(ie1-1)),abs(diff1(ie1)))
         if (diff1(ie1)*diff1(ie1-1).gt.zero) then
           slope1(ic0,ie1)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(1) 
         else
            slope1(ic0,ie1)=zero
         endif
      enddo
      enddo

      do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         ir1=if1-ic1*ratio(1)
         deltax1=deltax(ir1,1)
         do if0=ifirstf0,ilastf0+0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            ir0=if0-ic0*ratio(0)
            arrayf(if0,if1)=arrayc(ic0,ic1)
     &                +slope0(ic0,ic1)*deltax(ir0,0)
     &                +slope1(ic0,ic1)*deltax1
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 2d side-centered complex data
c***********************************************************************
c
c      subroutine cartclinrefsidecplx2d0(
ccart_clinref_op_side_2d(double complex,0,1)cc
c      subroutine cartclinrefsidecplx2d1(
ccart_clinref_op_side_2d(double complex,1,0)c
