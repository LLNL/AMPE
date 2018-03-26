c
c  File:        $URL$
c  Package:     SAMRAI geometry
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: FORTRAN routines for spatial refining of 1d patch data
c               on a regular Cartesian mesh.
c
c
c  File:        $URL$
c  Package:     SAMRAI geometry
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for 1d Cartesian refine operators
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for dimensioning 1d arrays in FORTRAN routines.
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
c***********************************************************************
c Linear interpolation for 1d cell-centered double data
c***********************************************************************
c
      subroutine cartlinrefcelldoub1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0
      integer ratio(0:1-1)
      double precision
     &  dxc(0:1-1),
     &  dxf(0:1-1)
      double precision
     &  arrayc(cilo0:cihi0),
     &  arrayf(filo0:fihi0)
      double precision deltax(0:15,0:1-1),x
      integer ic0,if0,ir0
c
c***********************************************************************
c
      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

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
        arrayf(if0)=arrayc(ic0)+(arrayc(ic0+1)-arrayc(ic0))*x
      enddo
c
      return
      end
c
c***********************************************************************
c Linear interpolation for 1d cell-centered float data
c***********************************************************************
c
      subroutine cartlinrefcellflot1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0
      integer ratio(0:1-1)
      double precision
     &  dxc(0:1-1),
     &  dxf(0:1-1)
      real
     &  arrayc(cilo0:cihi0),
     &  arrayf(filo0:fihi0)
      double precision deltax(0:15,0:1-1),x
      integer ic0,if0,ir0
c
c***********************************************************************
c
      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

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
        arrayf(if0)=arrayc(ic0)+(arrayc(ic0+1)-arrayc(ic0))*x
      enddo
c
      return
      end
c
c***********************************************************************
c Linear interpolation for 1d cell-centered complex data
c***********************************************************************
c
      subroutine cartlinrefcellcplx1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0
      integer ratio(0:1-1)
      double precision
     &  dxc(0:1-1),
     &  dxf(0:1-1)
      double complex
     &  arrayc(cilo0:cihi0),
     &  arrayf(filo0:fihi0)
      double precision deltax(0:15,0:1-1),x
      integer ic0,if0,ir0
c
c***********************************************************************
c
      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

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
        arrayf(if0)=arrayc(ic0)+(arrayc(ic0+1)-arrayc(ic0))*x
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 1d cell-centered double data
c***********************************************************************
c
      subroutine cartclinrefcelldoub1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0
      integer ratio(0:1-1)
      double precision
     &  dxc(0:1-1),
     &  dxf(0:1-1),
     &  deltax(0:15,0:1-1)
      double precision
     &  arrayc(cilo0:cihi0),
     &  arrayf(filo0:fihi0),
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
         arrayf(if0)=arrayc(ic0)+slope0(ic0)*deltax(ir0,0)
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 1d cell-centered float data
c***********************************************************************
c
      subroutine cartclinrefcellflot1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0
      integer ratio(0:1-1)
      double precision
     &  dxc(0:1-1),
     &  dxf(0:1-1),
     &  deltax(0:15,0:1-1)
      real
     &  arrayc(cilo0:cihi0),
     &  arrayf(filo0:fihi0),
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
         arrayf(if0)=arrayc(ic0)+slope0(ic0)*deltax(ir0,0)
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 1d cell-centered complex data
c***********************************************************************
c
      subroutine cartclinrefcellcplx1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0
      integer ratio(0:1-1)
      double precision
     &  dxc(0:1-1),
     &  dxf(0:1-1),
     &  deltax(0:15,0:1-1)
      double complex
     &  arrayc(cilo0:cihi0),
     &  arrayf(filo0:fihi0),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0)
      integer ic0,ie0,if0,ir0
      double precision
     &  coef2real,coef2imag,boundreal,boundimag,
     &  diff0real,diff0imag,diff1real,diff1imag,
     &  slopereal,slopeimag
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
         slope0(ic0) = slopereal + cmplx(zero,one)*slopeimag
      enddo

      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         ir0=if0-ic0*ratio(0)
         arrayf(if0)=arrayc(ic0)+slope0(ic0)*deltax(ir0,0)
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 1d edge-centered double data
c***********************************************************************
c
      subroutine cartclinrefedgedoub1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0
      integer ratio(0:1-1)
      double precision
     &  dxc(0:1-1),
     &  dxf(0:1-1),
     &  deltax(0:15,0:1-1)
      double precision
     &  arrayc(cilo0:cihi0),
     &  arrayf(filo0:fihi0),
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
         arrayf(if0)=arrayc(ic0)+slope0(ic0)*deltax(ir0,0)
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 1d edge-centered float data
c***********************************************************************
c
      subroutine cartclinrefedgeflot1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0
      integer ratio(0:1-1)
      double precision
     &  dxc(0:1-1),
     &  dxf(0:1-1),
     &  deltax(0:15,0:1-1)
      real
     &  arrayc(cilo0:cihi0),
     &  arrayf(filo0:fihi0),
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
         arrayf(if0)=arrayc(ic0)+slope0(ic0)*deltax(ir0,0)
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 1d face-centered double data
c***********************************************************************
c
      subroutine cartclinreffacedoub1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0
      integer ratio(0:1-1)
      double precision
     &  dxc(0:1-1),
     &  dxf(0:1-1),
     &  deltax(0:15,0:1-1)
      double precision
     &  arrayc(cilo0:cihi0+1),
     &  arrayf(filo0:fihi0+1),
     &  diff0(cilo0:cihi0),
     &  slope0(cilo0:cihi0+1)
      integer ic0,ie0,if0,ir0
      double precision
     &  coef2,bound
c
c***********************************************************************
c

      do ir0=0,ratio(0)-1
         deltax(ir0,0)=dble(ir0)*dxf(0)
      enddo

      do ic0=ifirstc0-1,ilastc0+1
         diff0(ic0)=arrayc(ic0+1)
     &                -arrayc(ic0)
      enddo
      do ie0=ifirstc0,ilastc0+1
         coef2=half*(diff0(ie0-1)+diff0(ie0))
         bound=two*min(abs(diff0(ie0-1)),abs(diff0(ie0)))
         if (diff0(ie0)*diff0(ie0-1).gt.zero) then
           slope0(ie0)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(0) 
         else
            slope0(ie0)=zero
         endif
      enddo

      do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         ir0=if0-ic0*ratio(0)
         arrayf(if0)=arrayc(ic0)+slope0(ic0)*deltax(ir0,0)
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 1d face-centered float data
c***********************************************************************
c
      subroutine cartclinreffaceflot1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0
      integer ratio(0:1-1)
      double precision
     &  dxc(0:1-1),
     &  dxf(0:1-1),
     &  deltax(0:15,0:1-1)
      real
     &  arrayc(cilo0:cihi0+1),
     &  arrayf(filo0:fihi0+1),
     &  diff0(cilo0:cihi0),
     &  slope0(cilo0:cihi0+1)
      integer ic0,ie0,if0,ir0
      real
     &  coef2,bound
c
c***********************************************************************
c

      do ir0=0,ratio(0)-1
         deltax(ir0,0)=dble(ir0)*dxf(0)
      enddo

      do ic0=ifirstc0-1,ilastc0+1
         diff0(ic0)=arrayc(ic0+1)
     &                -arrayc(ic0)
      enddo
      do ie0=ifirstc0,ilastc0+1
         coef2=half*(diff0(ie0-1)+diff0(ie0))
         bound=two*min(abs(diff0(ie0-1)),abs(diff0(ie0)))
         if (diff0(ie0)*diff0(ie0-1).gt.zero) then
           slope0(ie0)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(0) 
         else
            slope0(ie0)=zero
         endif
      enddo

      do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         ir0=if0-ic0*ratio(0)
         arrayf(if0)=arrayc(ic0)+slope0(ic0)*deltax(ir0,0)
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 1d face-centered complex data
c***********************************************************************
c
c      subroutine cartclinreffacecplx1d(
ccart_clinref_op_face_1d(double complex)c
c***********************************************************************
c Linear interpolation for 1d node-centered double data
c***********************************************************************
c
       subroutine cartlinrefnodedoub1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0
      integer ratio(0:1-1)
      double precision
     &  dxc(0:1-1),
     &  dxf(0:1-1)
      double precision
     &  arrayc(cilo0:cihi0+1),
     &  arrayf(filo0:fihi0+1)
      double precision realrat,x
      integer i,ic0,if0,ie0,ir0
c
c***********************************************************************
c
      realrat=one/dble(ratio(0))

      do ic0=ifirstc0,ilastc0
         if0=ic0*ratio(0)
         if (if0.ge.filo0.and.if0.le.fihi0+1) then 
            do ir0=0,ratio(0)-1
               ie0=if0+ir0
               x = dble(ir0)*realrat
               if (ie0.ge.filo0.and.ie0.le.fihi0+1) then
                  arrayf(ie0) = arrayc(ic0)*(one-x)
     &                              + arrayc(ic0+1)*x
               endif
            enddo
         endif
      enddo
c
      ic0 = ilastc0+1
      if0 = ic0*ratio(0)
      if (if0.ge.filo0.and.if0.le.fihi0+1) then
         arrayf(if0) = arrayc(ic0)
      endif
c
      return
      end
c
c***********************************************************************
c Linear interpolation for 1d node-centered float data
c***********************************************************************
c
       subroutine cartlinrefnodeflot1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0
      integer ratio(0:1-1)
      double precision
     &  dxc(0:1-1),
     &  dxf(0:1-1)
      real
     &  arrayc(cilo0:cihi0+1),
     &  arrayf(filo0:fihi0+1)
      double precision realrat,x
      integer i,ic0,if0,ie0,ir0
c
c***********************************************************************
c
      realrat=one/dble(ratio(0))

      do ic0=ifirstc0,ilastc0
         if0=ic0*ratio(0)
         if (if0.ge.filo0.and.if0.le.fihi0+1) then 
            do ir0=0,ratio(0)-1
               ie0=if0+ir0
               x = dble(ir0)*realrat
               if (ie0.ge.filo0.and.ie0.le.fihi0+1) then
                  arrayf(ie0) = arrayc(ic0)*(one-x)
     &                              + arrayc(ic0+1)*x
               endif
            enddo
         endif
      enddo
c
      ic0 = ilastc0+1
      if0 = ic0*ratio(0)
      if (if0.ge.filo0.and.if0.le.fihi0+1) then
         arrayf(if0) = arrayc(ic0)
      endif
c
      return
      end
c
c***********************************************************************
c Linear interpolation for 1d node-centered complex data
c***********************************************************************
c
       subroutine cartlinrefnodecplx1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0
      integer ratio(0:1-1)
      double precision
     &  dxc(0:1-1),
     &  dxf(0:1-1)
      double complex
     &  arrayc(cilo0:cihi0+1),
     &  arrayf(filo0:fihi0+1)
      double precision realrat,x
      integer i,ic0,if0,ie0,ir0
c
c***********************************************************************
c
      realrat=one/dble(ratio(0))

      do ic0=ifirstc0,ilastc0
         if0=ic0*ratio(0)
         if (if0.ge.filo0.and.if0.le.fihi0+1) then 
            do ir0=0,ratio(0)-1
               ie0=if0+ir0
               x = dble(ir0)*realrat
               if (ie0.ge.filo0.and.ie0.le.fihi0+1) then
                  arrayf(ie0) = arrayc(ic0)*(one-x)
     &                              + arrayc(ic0+1)*x
               endif
            enddo
         endif
      enddo
c
      ic0 = ilastc0+1
      if0 = ic0*ratio(0)
      if (if0.ge.filo0.and.if0.le.fihi0+1) then
         arrayf(if0) = arrayc(ic0)
      endif
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 1d outerface double data
c***********************************************************************
c
      subroutine cartclinrefoutfacedoub1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0
      integer ratio(0:1-1)
      double precision
     &  dxc(0:1-1),
     &  dxf(0:1-1)
      double precision
     &  arrayc(1),
     &  arrayf(1)
c
c***********************************************************************
c
      arrayf(1)=arrayc(1)
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 1d outerface float data
c***********************************************************************
c
      subroutine cartclinrefoutfaceflot1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0
      integer ratio(0:1-1)
      double precision
     &  dxc(0:1-1),
     &  dxf(0:1-1),
     &  deltax(0:15,0:1-1)
      real
     &  arrayc(cilo0:cihi0+1),
     &  arrayf(filo0:fihi0+1),
     &  diff0(cilo0:cihi0),
     &  slope0(cilo0:cihi0+1)
      integer ic0,ie0,if0,ir0
      real
     &  coef2,bound
c
c***********************************************************************
c

      do ir0=0,ratio(0)-1
         deltax(ir0,0)=dble(ir0)*dxf(0)
      enddo

      do ic0=ifirstc0-1,ilastc0+1
         diff0(ic0)=arrayc(ic0+1)
     &                -arrayc(ic0)
      enddo
      do ie0=ifirstc0,ilastc0+1
         coef2=half*(diff0(ie0-1)+diff0(ie0))
         bound=two*min(abs(diff0(ie0-1)),abs(diff0(ie0)))
         if (diff0(ie0)*diff0(ie0-1).gt.zero) then
           slope0(ie0)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(0) 
         else
            slope0(ie0)=zero
         endif
      enddo

      do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         ir0=if0-ic0*ratio(0)
         arrayf(if0)=arrayc(ic0)+slope0(ic0)*deltax(ir0,0)
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 1d outerface complex data
c***********************************************************************
c
c      subroutine cartclinrefoutfacecplx1d(
ccart_clinref_op_face_1d(double complex)c
c***********************************************************************
c Conservative linear interpolation for 1d side-centered double data
c***********************************************************************
c
      subroutine cartclinrefsidedoub1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0
      integer ratio(0:1-1)
      double precision
     &  dxc(0:1-1),
     &  dxf(0:1-1),
     &  deltax(0:15,0:1-1)
      double precision
     &  arrayc(cilo0:cihi0+1),
     &  arrayf(filo0:fihi0+1),
     &  diff0(cilo0:cihi0),
     &  slope0(cilo0:cihi0+1)
      integer ic0,ie0,if0,ir0
      double precision
     &  coef2,bound
c
c***********************************************************************
c

      do ir0=0,ratio(0)-1
         deltax(ir0,0)=dble(ir0)*dxf(0)
      enddo

      do ic0=ifirstc0-1,ilastc0+1
         diff0(ic0)=arrayc(ic0+1)
     &                -arrayc(ic0)
      enddo
      do ie0=ifirstc0,ilastc0+1
         coef2=half*(diff0(ie0-1)+diff0(ie0))
         bound=two*min(abs(diff0(ie0-1)),abs(diff0(ie0)))
         if (diff0(ie0)*diff0(ie0-1).gt.zero) then
           slope0(ie0)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(0) 
         else
            slope0(ie0)=zero
         endif
      enddo

      do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         ir0=if0-ic0*ratio(0)
         arrayf(if0)=arrayc(ic0)+slope0(ic0)*deltax(ir0,0)
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 1d side-centered float data
c***********************************************************************
c
      subroutine cartclinrefsideflot1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0
      integer ratio(0:1-1)
      double precision
     &  dxc(0:1-1),
     &  dxf(0:1-1),
     &  deltax(0:15,0:1-1)
      real
     &  arrayc(cilo0:cihi0+1),
     &  arrayf(filo0:fihi0+1),
     &  diff0(cilo0:cihi0),
     &  slope0(cilo0:cihi0+1)
      integer ic0,ie0,if0,ir0
      real
     &  coef2,bound
c
c***********************************************************************
c

      do ir0=0,ratio(0)-1
         deltax(ir0,0)=dble(ir0)*dxf(0)
      enddo

      do ic0=ifirstc0-1,ilastc0+1
         diff0(ic0)=arrayc(ic0+1)
     &                -arrayc(ic0)
      enddo
      do ie0=ifirstc0,ilastc0+1
         coef2=half*(diff0(ie0-1)+diff0(ie0))
         bound=two*min(abs(diff0(ie0-1)),abs(diff0(ie0)))
         if (diff0(ie0)*diff0(ie0-1).gt.zero) then
           slope0(ie0)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(0) 
         else
            slope0(ie0)=zero
         endif
      enddo

      do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         ir0=if0-ic0*ratio(0)
         arrayf(if0)=arrayc(ic0)+slope0(ic0)*deltax(ir0,0)
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 1d side-centered complex data
c***********************************************************************
c
c      subroutine cartclinrefsidecplx1d(
ccart_clinref_op_side_1d(double complex)c
