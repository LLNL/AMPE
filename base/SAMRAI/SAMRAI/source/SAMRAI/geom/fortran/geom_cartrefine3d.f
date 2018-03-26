c
c  File:        $URL$
c  Package:     SAMRAI geometry
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: FORTRAN routines for spatial refining of 3d patch data
c               on a regular Cartesian mesh.
c
c
c  File:        $URL$
c  Package:     SAMRAI geometry
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for 3d Cartesian refine operators
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for dimensioning 3d arrays in FORTRAN routines.
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
c Linear interpolation for 3d cell-centered double data
c***********************************************************************
c
      subroutine cartlinrefcelldoub3d(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1)
      double precision
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2)
      double precision deltax(0:15,0:3-1),x,y,z
      integer ic0,ic1,ic2,if0,if1,if2,ir0,ir1,ir2
c
c***********************************************************************
c
      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ir2=0,ratio(2)-1
         deltax(ir2,2)=(dble(ir2)+half)*dxf(2)-dxc(2)*half
      enddo

      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         z=deltax(ir2,2)/dxc(2)
         if( z .lt. 0.d0 ) then
           ic2 = ic2-1
           z = z + one
         endif
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
               arrayf(if0,if1,if2)=
     &          ( (arrayc(ic0,ic1,ic2)
     &      +(arrayc(ic0+1,ic1,ic2)-arrayc(ic0,ic1,ic2))*x)*(one-y)
     &      +(arrayc(ic0,ic1+1,ic2)
     &      +(arrayc(ic0+1,ic1+1,ic2)-arrayc(ic0,ic1+1,ic2))*x)*y )
     &       *(one-z)
     &         +( (arrayc(ic0,ic1,ic2+1)
     &      +(arrayc(ic0+1,ic1,ic2+1)-arrayc(ic0,ic1,ic2+1))*x)
     &       *(one-y)
     &      +(arrayc(ic0,ic1+1,ic2+1)
     &      +(arrayc(ic0+1,ic1+1,ic2+1)-arrayc(ic0,ic1+1,ic2+1))*x)*y)
     &       *z
          enddo
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear interpolation for 3d cell-centered float data
c***********************************************************************
c
      subroutine cartlinrefcellflot3d(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1)
      real
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2)
      double precision deltax(0:15,0:3-1),x,y,z
      integer ic0,ic1,ic2,if0,if1,if2,ir0,ir1,ir2
c
c***********************************************************************
c
      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ir2=0,ratio(2)-1
         deltax(ir2,2)=(dble(ir2)+half)*dxf(2)-dxc(2)*half
      enddo

      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         z=deltax(ir2,2)/dxc(2)
         if( z .lt. 0.d0 ) then
           ic2 = ic2-1
           z = z + one
         endif
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
               arrayf(if0,if1,if2)=
     &          ( (arrayc(ic0,ic1,ic2)
     &      +(arrayc(ic0+1,ic1,ic2)-arrayc(ic0,ic1,ic2))*x)*(one-y)
     &      +(arrayc(ic0,ic1+1,ic2)
     &      +(arrayc(ic0+1,ic1+1,ic2)-arrayc(ic0,ic1+1,ic2))*x)*y )
     &       *(one-z)
     &         +( (arrayc(ic0,ic1,ic2+1)
     &      +(arrayc(ic0+1,ic1,ic2+1)-arrayc(ic0,ic1,ic2+1))*x)
     &       *(one-y)
     &      +(arrayc(ic0,ic1+1,ic2+1)
     &      +(arrayc(ic0+1,ic1+1,ic2+1)-arrayc(ic0,ic1+1,ic2+1))*x)*y)
     &       *z
          enddo
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear interpolation for 3d cell-centered complex data
c***********************************************************************
c
      subroutine cartlinrefcellcplx3d(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1)
      double complex
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2)
      double precision deltax(0:15,0:3-1),x,y,z
      integer ic0,ic1,ic2,if0,if1,if2,ir0,ir1,ir2
c
c***********************************************************************
c
      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ir2=0,ratio(2)-1
         deltax(ir2,2)=(dble(ir2)+half)*dxf(2)-dxc(2)*half
      enddo

      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         z=deltax(ir2,2)/dxc(2)
         if( z .lt. 0.d0 ) then
           ic2 = ic2-1
           z = z + one
         endif
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
               arrayf(if0,if1,if2)=
     &          ( (arrayc(ic0,ic1,ic2)
     &      +(arrayc(ic0+1,ic1,ic2)-arrayc(ic0,ic1,ic2))*x)*(one-y)
     &      +(arrayc(ic0,ic1+1,ic2)
     &      +(arrayc(ic0+1,ic1+1,ic2)-arrayc(ic0,ic1+1,ic2))*x)*y )
     &       *(one-z)
     &         +( (arrayc(ic0,ic1,ic2+1)
     &      +(arrayc(ic0+1,ic1,ic2+1)-arrayc(ic0,ic1,ic2+1))*x)
     &       *(one-y)
     &      +(arrayc(ic0,ic1+1,ic2+1)
     &      +(arrayc(ic0+1,ic1+1,ic2+1)-arrayc(ic0,ic1+1,ic2+1))*x)*y)
     &       *z
          enddo
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 3d cell-centered double data
c***********************************************************************
c
      subroutine cartclinrefcelldoub3d(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1,diff2,slope2)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      double precision
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2)
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,ir0,ir1,ir2
      double precision
     &  coef2,bound
      double precision deltax1,deltax2
c
c***********************************************************************
c
      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ir2=0,ratio(2)-1
         deltax(ir2,2)=(dble(ir2)+half)*dxf(2)-dxc(2)*half
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ie0,ic1,ic2)
     &               -arrayc(ie0-1,ic1,ic2)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic0=ifirstc0,ilastc0
      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ic0,ie1,ic2)
     &               -arrayc(ic0,ie1-1,ic2)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
      do ie2=ifirstc2,ilastc2+1
         diff2(ie2)=arrayc(ic0,ic1,ie2)
     &               -arrayc(ic0,ic1,ie2-1)
      enddo
      do ic2=ifirstc2,ilastc2
         coef2=half*(diff2(ic2+1)+diff2(ic2))
         bound=two*min(abs(diff2(ic2+1)),abs(diff2(ic2)))
         if (diff2(ic2)*diff2(ic2+1).gt.zero) then
            slope2(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(2)
         else
            slope2(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo

      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)
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
               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
     &                     +slope0(ic0,ic1,ic2)*deltax(ir0,0)
     &                     +slope1(ic0,ic1,ic2)*deltax1
     &                     +slope2(ic0,ic1,ic2)*deltax2
          enddo
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 3d cell-centered float data
c***********************************************************************
c
      subroutine cartclinrefcellflot3d(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1,diff2,slope2)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      real
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2)
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,ir0,ir1,ir2
      real
     &  coef2,bound
      double precision deltax1,deltax2
c
c***********************************************************************
c
      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ir2=0,ratio(2)-1
         deltax(ir2,2)=(dble(ir2)+half)*dxf(2)-dxc(2)*half
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ie0,ic1,ic2)
     &               -arrayc(ie0-1,ic1,ic2)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic0=ifirstc0,ilastc0
      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ic0,ie1,ic2)
     &               -arrayc(ic0,ie1-1,ic2)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
      do ie2=ifirstc2,ilastc2+1
         diff2(ie2)=arrayc(ic0,ic1,ie2)
     &               -arrayc(ic0,ic1,ie2-1)
      enddo
      do ic2=ifirstc2,ilastc2
         coef2=half*(diff2(ic2+1)+diff2(ic2))
         bound=two*min(abs(diff2(ic2+1)),abs(diff2(ic2)))
         if (diff2(ic2)*diff2(ic2+1).gt.zero) then
            slope2(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(2)
         else
            slope2(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo

      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)
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
               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
     &                     +slope0(ic0,ic1,ic2)*deltax(ir0,0)
     &                     +slope1(ic0,ic1,ic2)*deltax1
     &                     +slope2(ic0,ic1,ic2)*deltax2
          enddo
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 3d cell-centered complex data
c***********************************************************************
c
      subroutine cartclinrefcellcplx3d(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1,diff2,slope2)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      double complex
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2)
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,ir0,ir1,ir2
      double precision
     &  coef2real,coef2imag,boundreal,boundimag,
     &  diff0real,diff0imag,diff1real,diff1imag,
     &  slopereal,slopeimag
      double precision deltax1,deltax2
c
c***********************************************************************
c
      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ir2=0,ratio(2)-1
         deltax(ir2,2)=(dble(ir2)+half)*dxf(2)-dxc(2)*half
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ie0,ic1,ic2)
     &               -arrayc(ie0-1,ic1,ic2)
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
         slope0(ic0,ic1,ic2) = slopereal + cmplx(zero,one)*slopeimag
      enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic0=ifirstc0,ilastc0
      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ic0,ie1,ic2)
     &               -arrayc(ic0,ie1-1,ic2)
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
         slope1(ic0,ic1,ic2) = slopereal + cmplx(zero,one)*slopeimag
      enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
      do ie2=ifirstc2,ilastc2+1
         diff2(ie2)=arrayc(ic0,ic1,ie2)
     &               -arrayc(ic0,ic1,ie2-1)
      enddo
      do ic2=ifirstc2,ilastc2
         diff0real=dble(diff2(ic2))
         diff0imag=imag(diff2(ic2))
         diff1real=dble(diff2(ic2+1))
         diff1imag=imag(diff2(ic2+1))
         coef2real=half*(diff0real+diff1real)
         coef2imag=half*(diff0imag+diff1imag)
         boundreal=two*min(abs(diff1real),abs(diff0real))
         boundimag=two*min(abs(diff1imag),abs(diff0imag))
         if (diff0real*diff1real.gt.zero) then
            slopereal=sign(min(abs(coef2real),boundreal),coef2real)
     &                /dxc(2)
         else
            slopereal=zero
         endif
         if (diff0imag*diff1imag.gt.zero) then
            slopeimag=sign(min(abs(coef2imag),boundimag),coef2imag) 
     &                /dxc(2)
         else
            slopeimag=zero
         endif
         slope2(ic0,ic1,ic2) = slopereal + cmplx(zero,one)*slopeimag
      enddo
         enddo
      enddo

      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)
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
               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
     &                     +slope0(ic0,ic1,ic2)*deltax(ir0,0)
     &                     +slope1(ic0,ic1,ic2)*deltax1
     &                     +slope2(ic0,ic1,ic2)*deltax2
          enddo
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 3d edge-centered double data
c***********************************************************************
c
      subroutine cartclinrefedgedoub3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1,diff2,slope2)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      double precision
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1,
     &          filo2:fihi2+1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1)
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,
     &        ir0,ir1,ir2
      double precision
     &  coef2,bound
      double precision deltax1,deltax2
c
c***********************************************************************
c


      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo


      do ir1=0,ratio(1)-1
         deltax(ir1,1)=dble(ir1)*dxf(1)
      enddo


      do ir2=0,ratio(2)-1
         deltax(ir2,2)=dble(ir2)*dxf(2)
      enddo


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1+1

      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ie0,ic1,ic2)
     &               -arrayc(ie0-1,ic1,ic2)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1

            do ic0=ifirstc0,ilastc0


      do ic1=ifirstc1-1,ilastc1+1
         diff1(ic1)=arrayc(ic0,ic1+1,ic2)
     &                -arrayc(ic0,ic1,ic2)
      enddo
      do ie1=ifirstc1,ilastc1+1
         coef2=half*(diff1(ie1-1)+diff1(ie1))
         bound=two*min(abs(diff1(ie1-1)),abs(diff1(ie1)))
         if (diff1(ie1)*diff1(ie1-1).gt.zero) then
           slope1(ic0,ie1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(1) 
         else
            slope1(ic0,ie1,ic2)=zero
         endif
      enddo
         enddo
      enddo


         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0

      do ic2=ifirstc2-1,ilastc2+1
         diff2(ic2)=arrayc(ic0,ic1,ic2+1)
     &                -arrayc(ic0,ic1,ic2)
      enddo
      do ie2=ifirstc2,ilastc2+1
         coef2=half*(diff2(ie2-1)+diff2(ie2))
         bound=two*min(abs(diff2(ie2-1)),abs(diff2(ie2)))
         if (diff2(ie2)*diff2(ie2-1).gt.zero) then
           slope2(ic0,ic1,ie2)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(2) 
         else
            slope2(ic0,ic1,ie2)=zero
         endif
      enddo
         enddo
      enddo


      do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)


         do if1=ifirstf1,ilastf1+1
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
               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
     &                +slope0(ic0,ic1,ic2)*deltax(ir0,0)
     &                +slope1(ic0,ic1,ic2)*deltax1
     &                +slope2(ic0,ic1,ic2)*deltax2
          enddo
        enddo
      enddo
c
      return
      end
c
      subroutine cartclinrefedgedoub3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff1,slope1,diff2,slope2,diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      double precision
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2+1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1)
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,
     &        ir0,ir1,ir2
      double precision
     &  coef2,bound
      double precision deltax1,deltax2
c
c***********************************************************************
c


      do ir0=0,ratio(0)-1
         deltax(ir0,0)=dble(ir0)*dxf(0)
      enddo


      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo


      do ir2=0,ratio(2)-1
         deltax(ir2,2)=dble(ir2)*dxf(2)
      enddo


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1

      do ic0=ifirstc0-1,ilastc0+1
         diff0(ic0)=arrayc(ic0+1,ic1,ic2)
     &                -arrayc(ic0,ic1,ic2)
      enddo
      do ie0=ifirstc0,ilastc0+1
         coef2=half*(diff0(ie0-1)+diff0(ie0))
         bound=two*min(abs(diff0(ie0-1)),abs(diff0(ie0)))
         if (diff0(ie0)*diff0(ie0-1).gt.zero) then
           slope0(ie0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(0) 
         else
            slope0(ie0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1

            do ic0=ifirstc0,ilastc0+1


      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ic0,ie1,ic2)
     &               -arrayc(ic0,ie1-1,ic2)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


         do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0+1

      do ic2=ifirstc2-1,ilastc2+1
         diff2(ic2)=arrayc(ic0,ic1,ic2+1)
     &                -arrayc(ic0,ic1,ic2)
      enddo
      do ie2=ifirstc2,ilastc2+1
         coef2=half*(diff2(ie2-1)+diff2(ie2))
         bound=two*min(abs(diff2(ie2-1)),abs(diff2(ie2)))
         if (diff2(ie2)*diff2(ie2-1).gt.zero) then
           slope2(ic0,ic1,ie2)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(2) 
         else
            slope2(ic0,ic1,ie2)=zero
         endif
      enddo
         enddo
      enddo


      do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)


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
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
               ir0=if0-ic0*ratio(0)
               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
     &                +slope0(ic0,ic1,ic2)*deltax(ir0,0)
     &                +slope1(ic0,ic1,ic2)*deltax1
     &                +slope2(ic0,ic1,ic2)*deltax2
          enddo
        enddo
      enddo
c
      return
      end
c
      subroutine cartclinrefedgedoub3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff2,slope2,diff0,slope0,diff1,slope1)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      double precision
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1,
     &          filo2:fihi2),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2)
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,
     &        ir0,ir1,ir2
      double precision
     &  coef2,bound
      double precision deltax1,deltax2
c
c***********************************************************************
c


      do ir0=0,ratio(0)-1
         deltax(ir0,0)=dble(ir0)*dxf(0)
      enddo


      do ir1=0,ratio(1)-1
         deltax(ir1,1)=dble(ir1)*dxf(1)
      enddo


      do ir2=0,ratio(2)-1
         deltax(ir2,2)=(dble(ir2)+half)*dxf(2)-dxc(2)*half
      enddo


      do ic2=ifirstc2,ilastc2

         do ic1=ifirstc1,ilastc1+1

      do ic0=ifirstc0-1,ilastc0+1
         diff0(ic0)=arrayc(ic0+1,ic1,ic2)
     &                -arrayc(ic0,ic1,ic2)
      enddo
      do ie0=ifirstc0,ilastc0+1
         coef2=half*(diff0(ie0-1)+diff0(ie0))
         bound=two*min(abs(diff0(ie0-1)),abs(diff0(ie0)))
         if (diff0(ie0)*diff0(ie0-1).gt.zero) then
           slope0(ie0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(0) 
         else
            slope0(ie0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2

            do ic0=ifirstc0,ilastc0+1


      do ic1=ifirstc1-1,ilastc1+1
         diff1(ic1)=arrayc(ic0,ic1+1,ic2)
     &                -arrayc(ic0,ic1,ic2)
      enddo
      do ie1=ifirstc1,ilastc1+1
         coef2=half*(diff1(ie1-1)+diff1(ie1))
         bound=two*min(abs(diff1(ie1-1)),abs(diff1(ie1)))
         if (diff1(ie1)*diff1(ie1-1).gt.zero) then
           slope1(ic0,ie1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(1) 
         else
            slope1(ic0,ie1,ic2)=zero
         endif
      enddo
         enddo
      enddo


         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0+1

      do ie2=ifirstc2,ilastc2+1
         diff2(ie2)=arrayc(ic0,ic1,ie2)
     &               -arrayc(ic0,ic1,ie2-1)
      enddo
      do ic2=ifirstc2,ilastc2
         coef2=half*(diff2(ic2+1)+diff2(ic2))
         bound=two*min(abs(diff2(ic2+1)),abs(diff2(ic2)))
         if (diff2(ic2)*diff2(ic2+1).gt.zero) then
            slope2(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(2)
         else
            slope2(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)


         do if1=ifirstf1,ilastf1+1
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
               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
     &                +slope0(ic0,ic1,ic2)*deltax(ir0,0)
     &                +slope1(ic0,ic1,ic2)*deltax1
     &                +slope2(ic0,ic1,ic2)*deltax2
          enddo
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 3d edge-centered float data
c***********************************************************************
c
      subroutine cartclinrefedgeflot3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1,diff2,slope2)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      real
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1,
     &          filo2:fihi2+1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1)
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,
     &        ir0,ir1,ir2
      real
     &  coef2,bound
      double precision deltax1,deltax2
c
c***********************************************************************
c


      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo


      do ir1=0,ratio(1)-1
         deltax(ir1,1)=dble(ir1)*dxf(1)
      enddo


      do ir2=0,ratio(2)-1
         deltax(ir2,2)=dble(ir2)*dxf(2)
      enddo


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1+1

      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ie0,ic1,ic2)
     &               -arrayc(ie0-1,ic1,ic2)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1

            do ic0=ifirstc0,ilastc0


      do ic1=ifirstc1-1,ilastc1+1
         diff1(ic1)=arrayc(ic0,ic1+1,ic2)
     &                -arrayc(ic0,ic1,ic2)
      enddo
      do ie1=ifirstc1,ilastc1+1
         coef2=half*(diff1(ie1-1)+diff1(ie1))
         bound=two*min(abs(diff1(ie1-1)),abs(diff1(ie1)))
         if (diff1(ie1)*diff1(ie1-1).gt.zero) then
           slope1(ic0,ie1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(1) 
         else
            slope1(ic0,ie1,ic2)=zero
         endif
      enddo
         enddo
      enddo


         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0

      do ic2=ifirstc2-1,ilastc2+1
         diff2(ic2)=arrayc(ic0,ic1,ic2+1)
     &                -arrayc(ic0,ic1,ic2)
      enddo
      do ie2=ifirstc2,ilastc2+1
         coef2=half*(diff2(ie2-1)+diff2(ie2))
         bound=two*min(abs(diff2(ie2-1)),abs(diff2(ie2)))
         if (diff2(ie2)*diff2(ie2-1).gt.zero) then
           slope2(ic0,ic1,ie2)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(2) 
         else
            slope2(ic0,ic1,ie2)=zero
         endif
      enddo
         enddo
      enddo


      do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)


         do if1=ifirstf1,ilastf1+1
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
               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
     &                +slope0(ic0,ic1,ic2)*deltax(ir0,0)
     &                +slope1(ic0,ic1,ic2)*deltax1
     &                +slope2(ic0,ic1,ic2)*deltax2
          enddo
        enddo
      enddo
c
      return
      end
c
      subroutine cartclinrefedgeflot3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff1,slope1,diff2,slope2,diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      real
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2+1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1)
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,
     &        ir0,ir1,ir2
      real
     &  coef2,bound
      double precision deltax1,deltax2
c
c***********************************************************************
c


      do ir0=0,ratio(0)-1
         deltax(ir0,0)=dble(ir0)*dxf(0)
      enddo


      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo


      do ir2=0,ratio(2)-1
         deltax(ir2,2)=dble(ir2)*dxf(2)
      enddo


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1

      do ic0=ifirstc0-1,ilastc0+1
         diff0(ic0)=arrayc(ic0+1,ic1,ic2)
     &                -arrayc(ic0,ic1,ic2)
      enddo
      do ie0=ifirstc0,ilastc0+1
         coef2=half*(diff0(ie0-1)+diff0(ie0))
         bound=two*min(abs(diff0(ie0-1)),abs(diff0(ie0)))
         if (diff0(ie0)*diff0(ie0-1).gt.zero) then
           slope0(ie0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(0) 
         else
            slope0(ie0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1

            do ic0=ifirstc0,ilastc0+1


      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ic0,ie1,ic2)
     &               -arrayc(ic0,ie1-1,ic2)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


         do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0+1

      do ic2=ifirstc2-1,ilastc2+1
         diff2(ic2)=arrayc(ic0,ic1,ic2+1)
     &                -arrayc(ic0,ic1,ic2)
      enddo
      do ie2=ifirstc2,ilastc2+1
         coef2=half*(diff2(ie2-1)+diff2(ie2))
         bound=two*min(abs(diff2(ie2-1)),abs(diff2(ie2)))
         if (diff2(ie2)*diff2(ie2-1).gt.zero) then
           slope2(ic0,ic1,ie2)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(2) 
         else
            slope2(ic0,ic1,ie2)=zero
         endif
      enddo
         enddo
      enddo


      do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)


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
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
               ir0=if0-ic0*ratio(0)
               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
     &                +slope0(ic0,ic1,ic2)*deltax(ir0,0)
     &                +slope1(ic0,ic1,ic2)*deltax1
     &                +slope2(ic0,ic1,ic2)*deltax2
          enddo
        enddo
      enddo
c
      return
      end
c
      subroutine cartclinrefedgeflot3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff2,slope2,diff0,slope0,diff1,slope1)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      real
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1,
     &          filo2:fihi2),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2)
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,
     &        ir0,ir1,ir2
      real
     &  coef2,bound
      double precision deltax1,deltax2
c
c***********************************************************************
c


      do ir0=0,ratio(0)-1
         deltax(ir0,0)=dble(ir0)*dxf(0)
      enddo


      do ir1=0,ratio(1)-1
         deltax(ir1,1)=dble(ir1)*dxf(1)
      enddo


      do ir2=0,ratio(2)-1
         deltax(ir2,2)=(dble(ir2)+half)*dxf(2)-dxc(2)*half
      enddo


      do ic2=ifirstc2,ilastc2

         do ic1=ifirstc1,ilastc1+1

      do ic0=ifirstc0-1,ilastc0+1
         diff0(ic0)=arrayc(ic0+1,ic1,ic2)
     &                -arrayc(ic0,ic1,ic2)
      enddo
      do ie0=ifirstc0,ilastc0+1
         coef2=half*(diff0(ie0-1)+diff0(ie0))
         bound=two*min(abs(diff0(ie0-1)),abs(diff0(ie0)))
         if (diff0(ie0)*diff0(ie0-1).gt.zero) then
           slope0(ie0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(0) 
         else
            slope0(ie0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2

            do ic0=ifirstc0,ilastc0+1


      do ic1=ifirstc1-1,ilastc1+1
         diff1(ic1)=arrayc(ic0,ic1+1,ic2)
     &                -arrayc(ic0,ic1,ic2)
      enddo
      do ie1=ifirstc1,ilastc1+1
         coef2=half*(diff1(ie1-1)+diff1(ie1))
         bound=two*min(abs(diff1(ie1-1)),abs(diff1(ie1)))
         if (diff1(ie1)*diff1(ie1-1).gt.zero) then
           slope1(ic0,ie1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(1) 
         else
            slope1(ic0,ie1,ic2)=zero
         endif
      enddo
         enddo
      enddo


         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0+1

      do ie2=ifirstc2,ilastc2+1
         diff2(ie2)=arrayc(ic0,ic1,ie2)
     &               -arrayc(ic0,ic1,ie2-1)
      enddo
      do ic2=ifirstc2,ilastc2
         coef2=half*(diff2(ic2+1)+diff2(ic2))
         bound=two*min(abs(diff2(ic2+1)),abs(diff2(ic2)))
         if (diff2(ic2)*diff2(ic2+1).gt.zero) then
            slope2(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(2)
         else
            slope2(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)


         do if1=ifirstf1,ilastf1+1
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
               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
     &                +slope0(ic0,ic1,ic2)*deltax(ir0,0)
     &                +slope1(ic0,ic1,ic2)*deltax1
     &                +slope2(ic0,ic1,ic2)*deltax2
          enddo
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 3d face-centered double data
c***********************************************************************
c
      subroutine cartclinreffacedoub3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1,diff2,slope2)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      double precision
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2)
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,
     &        ir0,ir1,ir2
      double precision
     &  coef2,bound
      double precision deltax1,deltax2
c
c***********************************************************************
c

      do ir2=0,ratio(2)-1
         deltax(ir2,2)=(dble(ir2)+half)*dxf(2)-dxc(2)*half
      enddo

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ir0=0,ratio(0)-1
         deltax(ir0,0)=dble(ir0)*dxf(0)
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
      do ic0=ifirstc0-1,ilastc0+1
         diff0(ic0)=arrayc(ic0+1,ic1,ic2)
     &                -arrayc(ic0,ic1,ic2)
      enddo
      do ie0=ifirstc0,ilastc0+1
         coef2=half*(diff0(ie0-1)+diff0(ie0))
         bound=two*min(abs(diff0(ie0-1)),abs(diff0(ie0)))
         if (diff0(ie0)*diff0(ie0-1).gt.zero) then
           slope0(ie0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(0) 
         else
            slope0(ie0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic0=ifirstc0,ilastc0+1
      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ic0,ie1,ic2)
     &               -arrayc(ic0,ie1-1,ic2)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0+1
      do ie2=ifirstc2,ilastc2+1
         diff2(ie2)=arrayc(ic0,ic1,ie2)
     &               -arrayc(ic0,ic1,ie2-1)
      enddo
      do ic2=ifirstc2,ilastc2
         coef2=half*(diff2(ic2+1)+diff2(ic2))
         bound=two*min(abs(diff2(ic2+1)),abs(diff2(ic2)))
         if (diff2(ic2)*diff2(ic2+1).gt.zero) then
            slope2(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(2)
         else
            slope2(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo

      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)
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
               arrayf(if0,if1,if2)=arrayc(ie0,ic1,ic2)
     &                +slope0(ie0,ic1,ic2)*deltax(ir0,0)
     &                +slope1(ie0,ic1,ic2)*deltax1
     &                +slope2(ie0,ic1,ic2)*deltax2
          enddo
        enddo
      enddo
c
      return
      end
c
      subroutine cartclinreffacedoub3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff1,slope1,diff2,slope2,diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      double precision
     &  arrayc(cilo1:cihi1+1,
     &          cilo2:cihi2,
     &          cilo0:cihi0),
     &  arrayf(filo1:fihi1+1,
     &          filo2:fihi2,
     &          filo0:fihi0),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo1:cihi1+1,
     &          cilo2:cihi2,
     &          cilo0:cihi0),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo1:cihi1+1,
     &          cilo2:cihi2,
     &          cilo0:cihi0),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo1:cihi1+1,
     &          cilo2:cihi2,
     &          cilo0:cihi0)
      integer ic1,ic2,ic0,ie1,ie2,ie0,if1,if2,if0,
     &        ir1,ir2,ir0
      double precision
     &  coef2,bound
      double precision deltax2,deltax0
c
c***********************************************************************
c

      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ir2=0,ratio(2)-1
         deltax(ir2,2)=(dble(ir2)+half)*dxf(2)-dxc(2)*half
      enddo

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=dble(ir1)*dxf(1)
      enddo

      do ic0=ifirstc0,ilastc0
         do ic2=ifirstc2,ilastc2
      do ic1=ifirstc1-1,ilastc1+1
         diff1(ic1)=arrayc(ic1+1,ic2,ic0)
     &                -arrayc(ic1,ic2,ic0)
      enddo
      do ie1=ifirstc1,ilastc1+1
         coef2=half*(diff1(ie1-1)+diff1(ie1))
         bound=two*min(abs(diff1(ie1-1)),abs(diff1(ie1)))
         if (diff1(ie1)*diff1(ie1-1).gt.zero) then
           slope1(ie1,ic2,ic0)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(1) 
         else
            slope1(ie1,ic2,ic0)=zero
         endif
      enddo
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         do ic1=ifirstc1,ilastc1+1
      do ie2=ifirstc2,ilastc2+1
         diff2(ie2)=arrayc(ic1,ie2,ic0)
     &               -arrayc(ic1,ie2-1,ic0)
      enddo
      do ic2=ifirstc2,ilastc2
         coef2=half*(diff2(ic2+1)+diff2(ic2))
         bound=two*min(abs(diff2(ic2+1)),abs(diff2(ic2)))
         if (diff2(ic2)*diff2(ic2+1).gt.zero) then
            slope2(ic1,ic2,ic0)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(2)
         else
            slope2(ic1,ic2,ic0)=zero
         endif
      enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1+1
      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ic1,ic2,ie0)
     &               -arrayc(ic1,ic2,ie0-1)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic1,ic2,ic0)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic1,ic2,ic0)=zero
         endif
      enddo
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
         do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
            ir2=if2-ic2*ratio(2)
            deltax2=deltax(ir2,2)
            do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ie1=(if1+1)/ratio(1)-1
         else
            ie1=if1/ratio(1)
         endif
               ir1=if1-ie1*ratio(1)
               arrayf(if1,if2,if0)=arrayc(ie1,ic2,ic0)
     &                +slope1(ie1,ic2,ic0)*deltax(ir1,1)
     &                +slope2(ie1,ic2,ic0)*deltax2
     &                +slope0(ie1,ic2,ic0)*deltax0
          enddo
        enddo
      enddo
c
      return
      end
c
      subroutine cartclinreffacedoub3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff2,slope2,diff0,slope0,diff1,slope1)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      double precision
     &  arrayc(cilo2:cihi2+1,
     &          cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo2:fihi2+1,
     &          filo0:fihi0,
     &          filo1:fihi1),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo2:cihi2+1,
     &          cilo0:cihi0,
     &          cilo1:cihi1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo2:cihi2+1,
     &          cilo0:cihi0,
     &          cilo1:cihi1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo2:cihi2+1,
     &          cilo0:cihi0,
     &          cilo1:cihi1)
      integer ic2,ic0,ic1,ie2,ie0,ie1,if2,if0,if1,
     &        ir2,ir0,ir1
      double precision
     &  coef2,bound
      double precision deltax0,deltax1
c
c***********************************************************************
c

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ir2=0,ratio(2)-1
         deltax(ir2,2)=dble(ir2)*dxf(2)
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
      do ic2=ifirstc2-1,ilastc2+1
         diff2(ic2)=arrayc(ic2+1,ic0,ic1)
     &                -arrayc(ic2,ic0,ic1)
      enddo
      do ie2=ifirstc2,ilastc2+1
         coef2=half*(diff2(ie2-1)+diff2(ie2))
         bound=two*min(abs(diff2(ie2-1)),abs(diff2(ie2)))
         if (diff2(ie2)*diff2(ie2-1).gt.zero) then
           slope2(ie2,ic0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(2) 
         else
            slope2(ie2,ic0,ic1)=zero
         endif
      enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic2=ifirstc2,ilastc2+1
      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ic2,ie0,ic1)
     &               -arrayc(ic2,ie0-1,ic1)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic2,ic0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic2,ic0,ic1)=zero
         endif
      enddo
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         do ic2=ifirstc2,ilastc2+1
      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ic2,ic0,ie1)
     &               -arrayc(ic2,ic0,ie1-1)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic2,ic0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic2,ic0,ic1)=zero
         endif
      enddo
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
            deltax0=deltax(ir0,0)
            do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ie2=(if2+1)/ratio(2)-1
         else
            ie2=if2/ratio(2)
         endif
               ir2=if2-ie2*ratio(2)
               arrayf(if2,if0,if1)=arrayc(ie2,ic0,ic1)
     &                +slope2(ie2,ic0,ic1)*deltax(ir2,2)
     &                +slope0(ie2,ic0,ic1)*deltax0
     &                +slope1(ie2,ic0,ic1)*deltax1
          enddo
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 3d face-centered float data
c***********************************************************************
c
      subroutine cartclinreffaceflot3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1,diff2,slope2)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      real
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2)
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,
     &        ir0,ir1,ir2
      real
     &  coef2,bound
      double precision deltax1,deltax2
c
c***********************************************************************
c

      do ir2=0,ratio(2)-1
         deltax(ir2,2)=(dble(ir2)+half)*dxf(2)-dxc(2)*half
      enddo

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ir0=0,ratio(0)-1
         deltax(ir0,0)=dble(ir0)*dxf(0)
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
      do ic0=ifirstc0-1,ilastc0+1
         diff0(ic0)=arrayc(ic0+1,ic1,ic2)
     &                -arrayc(ic0,ic1,ic2)
      enddo
      do ie0=ifirstc0,ilastc0+1
         coef2=half*(diff0(ie0-1)+diff0(ie0))
         bound=two*min(abs(diff0(ie0-1)),abs(diff0(ie0)))
         if (diff0(ie0)*diff0(ie0-1).gt.zero) then
           slope0(ie0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(0) 
         else
            slope0(ie0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic0=ifirstc0,ilastc0+1
      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ic0,ie1,ic2)
     &               -arrayc(ic0,ie1-1,ic2)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0+1
      do ie2=ifirstc2,ilastc2+1
         diff2(ie2)=arrayc(ic0,ic1,ie2)
     &               -arrayc(ic0,ic1,ie2-1)
      enddo
      do ic2=ifirstc2,ilastc2
         coef2=half*(diff2(ic2+1)+diff2(ic2))
         bound=two*min(abs(diff2(ic2+1)),abs(diff2(ic2)))
         if (diff2(ic2)*diff2(ic2+1).gt.zero) then
            slope2(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(2)
         else
            slope2(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo

      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)
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
               arrayf(if0,if1,if2)=arrayc(ie0,ic1,ic2)
     &                +slope0(ie0,ic1,ic2)*deltax(ir0,0)
     &                +slope1(ie0,ic1,ic2)*deltax1
     &                +slope2(ie0,ic1,ic2)*deltax2
          enddo
        enddo
      enddo
c
      return
      end
c
      subroutine cartclinreffaceflot3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff1,slope1,diff2,slope2,diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      real
     &  arrayc(cilo1:cihi1+1,
     &          cilo2:cihi2,
     &          cilo0:cihi0),
     &  arrayf(filo1:fihi1+1,
     &          filo2:fihi2,
     &          filo0:fihi0),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo1:cihi1+1,
     &          cilo2:cihi2,
     &          cilo0:cihi0),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo1:cihi1+1,
     &          cilo2:cihi2,
     &          cilo0:cihi0),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo1:cihi1+1,
     &          cilo2:cihi2,
     &          cilo0:cihi0)
      integer ic1,ic2,ic0,ie1,ie2,ie0,if1,if2,if0,
     &        ir1,ir2,ir0
      real
     &  coef2,bound
      double precision deltax2,deltax0
c
c***********************************************************************
c

      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ir2=0,ratio(2)-1
         deltax(ir2,2)=(dble(ir2)+half)*dxf(2)-dxc(2)*half
      enddo

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=dble(ir1)*dxf(1)
      enddo

      do ic0=ifirstc0,ilastc0
         do ic2=ifirstc2,ilastc2
      do ic1=ifirstc1-1,ilastc1+1
         diff1(ic1)=arrayc(ic1+1,ic2,ic0)
     &                -arrayc(ic1,ic2,ic0)
      enddo
      do ie1=ifirstc1,ilastc1+1
         coef2=half*(diff1(ie1-1)+diff1(ie1))
         bound=two*min(abs(diff1(ie1-1)),abs(diff1(ie1)))
         if (diff1(ie1)*diff1(ie1-1).gt.zero) then
           slope1(ie1,ic2,ic0)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(1) 
         else
            slope1(ie1,ic2,ic0)=zero
         endif
      enddo
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         do ic1=ifirstc1,ilastc1+1
      do ie2=ifirstc2,ilastc2+1
         diff2(ie2)=arrayc(ic1,ie2,ic0)
     &               -arrayc(ic1,ie2-1,ic0)
      enddo
      do ic2=ifirstc2,ilastc2
         coef2=half*(diff2(ic2+1)+diff2(ic2))
         bound=two*min(abs(diff2(ic2+1)),abs(diff2(ic2)))
         if (diff2(ic2)*diff2(ic2+1).gt.zero) then
            slope2(ic1,ic2,ic0)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(2)
         else
            slope2(ic1,ic2,ic0)=zero
         endif
      enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1+1
      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ic1,ic2,ie0)
     &               -arrayc(ic1,ic2,ie0-1)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic1,ic2,ic0)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic1,ic2,ic0)=zero
         endif
      enddo
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
         do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
            ir2=if2-ic2*ratio(2)
            deltax2=deltax(ir2,2)
            do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ie1=(if1+1)/ratio(1)-1
         else
            ie1=if1/ratio(1)
         endif
               ir1=if1-ie1*ratio(1)
               arrayf(if1,if2,if0)=arrayc(ie1,ic2,ic0)
     &                +slope1(ie1,ic2,ic0)*deltax(ir1,1)
     &                +slope2(ie1,ic2,ic0)*deltax2
     &                +slope0(ie1,ic2,ic0)*deltax0
          enddo
        enddo
      enddo
c
      return
      end
c
      subroutine cartclinreffaceflot3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff2,slope2,diff0,slope0,diff1,slope1)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      real
     &  arrayc(cilo2:cihi2+1,
     &          cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo2:fihi2+1,
     &          filo0:fihi0,
     &          filo1:fihi1),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo2:cihi2+1,
     &          cilo0:cihi0,
     &          cilo1:cihi1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo2:cihi2+1,
     &          cilo0:cihi0,
     &          cilo1:cihi1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo2:cihi2+1,
     &          cilo0:cihi0,
     &          cilo1:cihi1)
      integer ic2,ic0,ic1,ie2,ie0,ie1,if2,if0,if1,
     &        ir2,ir0,ir1
      real
     &  coef2,bound
      double precision deltax0,deltax1
c
c***********************************************************************
c

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ir2=0,ratio(2)-1
         deltax(ir2,2)=dble(ir2)*dxf(2)
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
      do ic2=ifirstc2-1,ilastc2+1
         diff2(ic2)=arrayc(ic2+1,ic0,ic1)
     &                -arrayc(ic2,ic0,ic1)
      enddo
      do ie2=ifirstc2,ilastc2+1
         coef2=half*(diff2(ie2-1)+diff2(ie2))
         bound=two*min(abs(diff2(ie2-1)),abs(diff2(ie2)))
         if (diff2(ie2)*diff2(ie2-1).gt.zero) then
           slope2(ie2,ic0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(2) 
         else
            slope2(ie2,ic0,ic1)=zero
         endif
      enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic2=ifirstc2,ilastc2+1
      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ic2,ie0,ic1)
     &               -arrayc(ic2,ie0-1,ic1)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic2,ic0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic2,ic0,ic1)=zero
         endif
      enddo
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         do ic2=ifirstc2,ilastc2+1
      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ic2,ic0,ie1)
     &               -arrayc(ic2,ic0,ie1-1)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic2,ic0,ic1)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic2,ic0,ic1)=zero
         endif
      enddo
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
            deltax0=deltax(ir0,0)
            do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ie2=(if2+1)/ratio(2)-1
         else
            ie2=if2/ratio(2)
         endif
               ir2=if2-ie2*ratio(2)
               arrayf(if2,if0,if1)=arrayc(ie2,ic0,ic1)
     &                +slope2(ie2,ic0,ic1)*deltax(ir2,2)
     &                +slope0(ie2,ic0,ic1)*deltax0
     &                +slope1(ie2,ic0,ic1)*deltax1
          enddo
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 3d face-centered complex data
c***********************************************************************
c
c     subroutine cartclinreffacecplx3d0(
ccart_clinref_op_face_3d(double complex,0,1,2)c
c      subroutine cartclinreffacecplx3d1(
ccart_clinref_op_face_3d(double complex,1,2,0)c
c      subroutine cartclinreffacecplx3d2(
ccart_clinref_op_face_3d(double complex,2,0,1)c
c***********************************************************************
c Linear interpolation for 3d node-centered double data
c***********************************************************************
c
       subroutine cartlinrefnodedoub3d(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1)
      double precision
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1,
     &          filo2:fihi2+1)
      double precision x,y,z,realrat0,realrat1,realrat2
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,ir0,ir1,ir2,i,j,k
c
c***********************************************************************
c
      realrat0=one/dble(ratio(0))
      realrat1=one/dble(ratio(1))
      realrat2=one/dble(ratio(2))

      do ic2=ifirstc2,ilastc2
         if2=ic2*ratio(2)
         do ir2=0,ratio(2)
            ie2=if2+ir2
            if ((ie2.ge.filo2).and.(ie2.le.(fihi2+1))) then
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
               z = dble(ir2)*realrat2
               arrayf(ie0,ie1,ie2)=
     &            ( (arrayc(ic0,ic1,ic2)*(one-x) +
     &               arrayc(ic0+1,ic1,ic2)*x)*(one-y)
     &            + (arrayc(ic0,ic1+1,ic2)*(one-x) +
     &               arrayc(ic0+1,ic1+1,ic2)*x)*y ) * (one-z) + 
     &            ( (arrayc(ic0,ic1,ic2+1)*(one-x) +
     &               arrayc(ic0+1,ic1,ic2+1)*x)*(one-y)
     &            + (arrayc(ic0,ic1+1,ic2+1)*(one-x) +
     &               arrayc(ic0+1,ic1+1,ic2+1)*x)*y ) * z
           endif
         end do
      end do
           endif
         end do
      end do
           endif
         end do
      end do

      return
      end
c
c***********************************************************************
c Linear interpolation for 3d node-centered float data
c***********************************************************************
c
       subroutine cartlinrefnodeflot3d(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1)
      real
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1,
     &          filo2:fihi2+1)
      double precision x,y,z,realrat0,realrat1,realrat2
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,ir0,ir1,ir2,i,j,k
c
c***********************************************************************
c
      realrat0=one/dble(ratio(0))
      realrat1=one/dble(ratio(1))
      realrat2=one/dble(ratio(2))

      do ic2=ifirstc2,ilastc2
         if2=ic2*ratio(2)
         do ir2=0,ratio(2)
            ie2=if2+ir2
            if ((ie2.ge.filo2).and.(ie2.le.(fihi2+1))) then
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
               z = dble(ir2)*realrat2
               arrayf(ie0,ie1,ie2)=
     &            ( (arrayc(ic0,ic1,ic2)*(one-x) +
     &               arrayc(ic0+1,ic1,ic2)*x)*(one-y)
     &            + (arrayc(ic0,ic1+1,ic2)*(one-x) +
     &               arrayc(ic0+1,ic1+1,ic2)*x)*y ) * (one-z) + 
     &            ( (arrayc(ic0,ic1,ic2+1)*(one-x) +
     &               arrayc(ic0+1,ic1,ic2+1)*x)*(one-y)
     &            + (arrayc(ic0,ic1+1,ic2+1)*(one-x) +
     &               arrayc(ic0+1,ic1+1,ic2+1)*x)*y ) * z
           endif
         end do
      end do
           endif
         end do
      end do
           endif
         end do
      end do

      return
      end
c
c***********************************************************************
c Linear interpolation for 3d node-centered complex data
c***********************************************************************
c
       subroutine cartlinrefnodecplx3d(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1)
      double complex
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1,
     &          filo2:fihi2+1)
      double precision x,y,z,realrat0,realrat1,realrat2
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,ir0,ir1,ir2,i,j,k
c
c***********************************************************************
c
      realrat0=one/dble(ratio(0))
      realrat1=one/dble(ratio(1))
      realrat2=one/dble(ratio(2))

      do ic2=ifirstc2,ilastc2
         if2=ic2*ratio(2)
         do ir2=0,ratio(2)
            ie2=if2+ir2
            if ((ie2.ge.filo2).and.(ie2.le.(fihi2+1))) then
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
               z = dble(ir2)*realrat2
               arrayf(ie0,ie1,ie2)=
     &            ( (arrayc(ic0,ic1,ic2)*(one-x) +
     &               arrayc(ic0+1,ic1,ic2)*x)*(one-y)
     &            + (arrayc(ic0,ic1+1,ic2)*(one-x) +
     &               arrayc(ic0+1,ic1+1,ic2)*x)*y ) * (one-z) + 
     &            ( (arrayc(ic0,ic1,ic2+1)*(one-x) +
     &               arrayc(ic0+1,ic1,ic2+1)*x)*(one-y)
     &            + (arrayc(ic0,ic1+1,ic2+1)*(one-x) +
     &               arrayc(ic0+1,ic1+1,ic2+1)*x)*y ) * z
           endif
         end do
      end do
           endif
         end do
      end do
           endif
         end do
      end do

      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 3d outerface double data
c***********************************************************************
c
      subroutine cartclinrefoutfacedoub3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1,diff2,slope2)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      double precision
     &  arrayc(cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo1:fihi1,
     &          filo2:fihi2),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo1:cihi1,
     &          cilo2:cihi2),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo1:cihi1,
     &          cilo2:cihi2),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo1:cihi1,
     &          cilo2:cihi2)
      integer ic1,ic2,ie1,ie2,if1,if2,ir1,ir2
      double precision
     &  coef2,bound
      double precision deltax2
c
c***********************************************************************
c

      do ir2=0,ratio(2)-1
         deltax(ir2,2)=(dble(ir2)+half)*dxf(2)-dxc(2)*half
      enddo

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ic2=ifirstc2,ilastc2
      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ie1,ic2)
     &               -arrayc(ie1-1,ic2)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic1,ic2)=zero
         endif
      enddo
      enddo

      do ic1=ifirstc1,ilastc1
      do ie2=ifirstc2,ilastc2+1
         diff2(ie2)=arrayc(ic1,ie2)
     &               -arrayc(ic1,ie2-1)
      enddo
      do ic2=ifirstc2,ilastc2
         coef2=half*(diff2(ic2+1)+diff2(ic2))
         bound=two*min(abs(diff2(ic2+1)),abs(diff2(ic2)))
         if (diff2(ic2)*diff2(ic2+1).gt.zero) then
            slope2(ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(2)
         else
            slope2(ic1,ic2)=zero
         endif
      enddo
      enddo

      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)
         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
            ir1=if1-ic1*ratio(1)
            arrayf(if1,if2)=arrayc(ic1,ic2)
     &            +slope1(ic1,ic2)*deltax(ir1,1)
     &            +slope2(ic1,ic2)*deltax2
         enddo
      enddo
c
      return
      end
c
      subroutine cartclinrefoutfacedoub3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff1,slope1,diff2,slope2,diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      double precision
     &  arrayc(cilo2:cihi2,
     &          cilo0:cihi0),
     &  arrayf(filo2:fihi2,
     &          filo0:fihi0),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo2:cihi2,
     &          cilo0:cihi0),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo2:cihi2,
     &          cilo0:cihi0),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo2:cihi2,
     &          cilo0:cihi0)
      integer ic2,ic0,ie2,ie0,if2,if0,ir2,ir0
      double precision
     &  coef2,bound
      double precision deltax0
c
c***********************************************************************
c

      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ir2=0,ratio(2)-1
         deltax(ir2,2)=(dble(ir2)+half)*dxf(2)-dxc(2)*half
      enddo

      do ic0=ifirstc0,ilastc0
      do ie2=ifirstc2,ilastc2+1
         diff2(ie2)=arrayc(ie2,ic0)
     &               -arrayc(ie2-1,ic0)
      enddo
      do ic2=ifirstc2,ilastc2
         coef2=half*(diff2(ic2+1)+diff2(ic2))
         bound=two*min(abs(diff2(ic2+1)),abs(diff2(ic2)))
         if (diff2(ic2)*diff2(ic2+1).gt.zero) then
            slope2(ic2,ic0)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(2)
         else
            slope2(ic2,ic0)=zero
         endif
      enddo
      enddo

      do ic2=ifirstc2,ilastc2
      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ic2,ie0)
     &               -arrayc(ic2,ie0-1)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic2,ic0)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic2,ic0)=zero
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
         do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
            ir2=if2-ic2*ratio(2)
            arrayf(if2,if0)=arrayc(ic2,ic0)
     &            +slope2(ic2,ic0)*deltax(ir2,2)
     &            +slope0(ic2,ic0)*deltax0
         enddo
      enddo
c
      return
      end
c
      subroutine cartclinrefoutfacedoub3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff2,slope2,diff0,slope0,diff1,slope1)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      double precision
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo0:cihi0,
     &          cilo1:cihi1),
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

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
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
     &            +slope0(ic0,ic1)*deltax(ir0,0)
     &            +slope1(ic0,ic1)*deltax1
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 3d outerface float data
c***********************************************************************
c
      subroutine cartclinrefoutfaceflot3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1,diff2,slope2)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      real
     &  arrayc(cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo1:fihi1,
     &          filo2:fihi2),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo1:cihi1,
     &          cilo2:cihi2),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo1:cihi1,
     &          cilo2:cihi2),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo1:cihi1,
     &          cilo2:cihi2)
      integer ic1,ic2,ie1,ie2,if1,if2,ir1,ir2
      real
     &  coef2,bound
      double precision deltax2
c
c***********************************************************************
c

      do ir2=0,ratio(2)-1
         deltax(ir2,2)=(dble(ir2)+half)*dxf(2)-dxc(2)*half
      enddo

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ic2=ifirstc2,ilastc2
      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ie1,ic2)
     &               -arrayc(ie1-1,ic2)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic1,ic2)=zero
         endif
      enddo
      enddo

      do ic1=ifirstc1,ilastc1
      do ie2=ifirstc2,ilastc2+1
         diff2(ie2)=arrayc(ic1,ie2)
     &               -arrayc(ic1,ie2-1)
      enddo
      do ic2=ifirstc2,ilastc2
         coef2=half*(diff2(ic2+1)+diff2(ic2))
         bound=two*min(abs(diff2(ic2+1)),abs(diff2(ic2)))
         if (diff2(ic2)*diff2(ic2+1).gt.zero) then
            slope2(ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(2)
         else
            slope2(ic1,ic2)=zero
         endif
      enddo
      enddo

      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)
         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
            ir1=if1-ic1*ratio(1)
            arrayf(if1,if2)=arrayc(ic1,ic2)
     &            +slope1(ic1,ic2)*deltax(ir1,1)
     &            +slope2(ic1,ic2)*deltax2
         enddo
      enddo
c
      return
      end
c
      subroutine cartclinrefoutfaceflot3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff1,slope1,diff2,slope2,diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      real
     &  arrayc(cilo2:cihi2,
     &          cilo0:cihi0),
     &  arrayf(filo2:fihi2,
     &          filo0:fihi0),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo2:cihi2,
     &          cilo0:cihi0),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo2:cihi2,
     &          cilo0:cihi0),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo2:cihi2,
     &          cilo0:cihi0)
      integer ic2,ic0,ie2,ie0,if2,if0,ir2,ir0
      real
     &  coef2,bound
      double precision deltax0
c
c***********************************************************************
c

      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ir2=0,ratio(2)-1
         deltax(ir2,2)=(dble(ir2)+half)*dxf(2)-dxc(2)*half
      enddo

      do ic0=ifirstc0,ilastc0
      do ie2=ifirstc2,ilastc2+1
         diff2(ie2)=arrayc(ie2,ic0)
     &               -arrayc(ie2-1,ic0)
      enddo
      do ic2=ifirstc2,ilastc2
         coef2=half*(diff2(ic2+1)+diff2(ic2))
         bound=two*min(abs(diff2(ic2+1)),abs(diff2(ic2)))
         if (diff2(ic2)*diff2(ic2+1).gt.zero) then
            slope2(ic2,ic0)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(2)
         else
            slope2(ic2,ic0)=zero
         endif
      enddo
      enddo

      do ic2=ifirstc2,ilastc2
      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ic2,ie0)
     &               -arrayc(ic2,ie0-1)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic2,ic0)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic2,ic0)=zero
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
         do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
            ir2=if2-ic2*ratio(2)
            arrayf(if2,if0)=arrayc(ic2,ic0)
     &            +slope2(ic2,ic0)*deltax(ir2,2)
     &            +slope0(ic2,ic0)*deltax0
         enddo
      enddo
c
      return
      end
c
      subroutine cartclinrefoutfaceflot3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff2,slope2,diff0,slope0,diff1,slope1)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      real
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo0:cihi0,
     &          cilo1:cihi1),
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

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
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
     &            +slope0(ic0,ic1)*deltax(ir0,0)
     &            +slope1(ic0,ic1)*deltax1
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 3d outerface complex data
c***********************************************************************
c
c     subroutine cartclinrefoutfacecplx3d0(
ccart_clinref_op_outerface_3d(double complex,0,1,2)c
c      subroutine cartclinrefoutfacecplx3d1(
ccart_clinref_op_outerface_3d(double complex,1,2,0)c
c      subroutine cartclinrefoutfacecplx3d2(
ccart_clinref_op_outerface_3d(double complex,2,0,1)c
c***********************************************************************
c Conservative linear interpolation for 3d side-centered double data
c***********************************************************************
c
      subroutine cartclinrefsidedoub3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1,diff2,slope2)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      double precision
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2)
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,
     &        ir0,ir1,ir2
      double precision
     &  coef2,bound
      double precision deltax2,deltax1
c
c***********************************************************************
c


      do ir0=0,ratio(0)-1
         deltax(ir0,0)=dble(ir0)*dxf(0)
      enddo


      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo


      do ir2=0,ratio(2)-1
         deltax(ir2,2)=(dble(ir2)+half)*dxf(2)-dxc(2)*half
      enddo


      do ic2=ifirstc2,ilastc2

         do ic1=ifirstc1,ilastc1

      do ic0=ifirstc0-1,ilastc0+1
         diff0(ic0)=arrayc(ic0+1,ic1,ic2)
     &                -arrayc(ic0,ic1,ic2)
      enddo
      do ie0=ifirstc0,ilastc0+1
         coef2=half*(diff0(ie0-1)+diff0(ie0))
         bound=two*min(abs(diff0(ie0-1)),abs(diff0(ie0)))
         if (diff0(ie0)*diff0(ie0-1).gt.zero) then
           slope0(ie0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(0) 
         else
            slope0(ie0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


            do ic2=ifirstc2,ilastc2

            do ic0=ifirstc0,ilastc0+1

      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ic0,ie1,ic2)
     &               -arrayc(ic0,ie1-1,ic2)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


            do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0+1

      do ie2=ifirstc2,ilastc2+1
         diff2(ie2)=arrayc(ic0,ic1,ie2)
     &               -arrayc(ic0,ic1,ie2-1)
      enddo
      do ic2=ifirstc2,ilastc2
         coef2=half*(diff2(ic2+1)+diff2(ic2))
         bound=two*min(abs(diff2(ic2+1)),abs(diff2(ic2)))
         if (diff2(ic2)*diff2(ic2+1).gt.zero) then
            slope2(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(2)
         else
            slope2(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)


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
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
               ir0=if0-ic0*ratio(0)

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
     &                +slope0(ic0,ic1,ic2)*deltax(ir0,0)
     &                +slope1(ic0,ic1,ic2)*deltax1
     &                +slope2(ic0,ic1,ic2)*deltax2
          enddo
        enddo
      enddo
c
      return
      end
c
      subroutine cartclinrefsidedoub3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff1,slope1,diff2,slope2,diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      double precision
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1,
     &          filo2:fihi2),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2)
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,
     &        ir0,ir1,ir2
      double precision
     &  coef2,bound
      double precision deltax2,deltax1
c
c***********************************************************************
c


      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo


      do ir1=0,ratio(1)-1
         deltax(ir1,1)=dble(ir1)*dxf(1)
      enddo


      do ir2=0,ratio(2)-1
         deltax(ir2,2)=(dble(ir2)+half)*dxf(2)-dxc(2)*half
      enddo


      do ic2=ifirstc2,ilastc2

         do ic1=ifirstc1,ilastc1+1

      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ie0,ic1,ic2)
     &               -arrayc(ie0-1,ic1,ic2)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


            do ic2=ifirstc2,ilastc2

            do ic0=ifirstc0,ilastc0

      do ic1=ifirstc1-1,ilastc1+1
         diff1(ic1)=arrayc(ic0,ic1+1,ic2)
     &                -arrayc(ic0,ic1,ic2)
      enddo
      do ie1=ifirstc1,ilastc1+1
         coef2=half*(diff1(ie1-1)+diff1(ie1))
         bound=two*min(abs(diff1(ie1-1)),abs(diff1(ie1)))
         if (diff1(ie1)*diff1(ie1-1).gt.zero) then
           slope1(ic0,ie1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(1) 
         else
            slope1(ic0,ie1,ic2)=zero
         endif
      enddo
         enddo
      enddo


            do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0

      do ie2=ifirstc2,ilastc2+1
         diff2(ie2)=arrayc(ic0,ic1,ie2)
     &               -arrayc(ic0,ic1,ie2-1)
      enddo
      do ic2=ifirstc2,ilastc2
         coef2=half*(diff2(ic2+1)+diff2(ic2))
         bound=two*min(abs(diff2(ic2+1)),abs(diff2(ic2)))
         if (diff2(ic2)*diff2(ic2+1).gt.zero) then
            slope2(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(2)
         else
            slope2(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)


         do if1=ifirstf1,ilastf1+1
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

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
     &                +slope0(ic0,ic1,ic2)*deltax(ir0,0)
     &                +slope1(ic0,ic1,ic2)*deltax1
     &                +slope2(ic0,ic1,ic2)*deltax2
          enddo
        enddo
      enddo
c
      return
      end
c
      subroutine cartclinrefsidedoub3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff2,slope2,diff0,slope0,diff1,slope1)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      double precision
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2+1),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1)
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,
     &        ir0,ir1,ir2
      double precision
     &  coef2,bound
      double precision deltax2,deltax1
c
c***********************************************************************
c


      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo


      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo


      do ir2=0,ratio(2)-1
         deltax(ir2,2)=dble(ir2)*dxf(2)
      enddo


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1

      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ie0,ic1,ic2)
     &               -arrayc(ie0-1,ic1,ic2)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


            do ic2=ifirstc2,ilastc2+1

            do ic0=ifirstc0,ilastc0

      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ic0,ie1,ic2)
     &               -arrayc(ic0,ie1-1,ic2)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


            do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0

      do ic2=ifirstc2-1,ilastc2+1
         diff2(ic2)=arrayc(ic0,ic1,ic2+1)
     &                -arrayc(ic0,ic1,ic2)
      enddo
      do ie2=ifirstc2,ilastc2+1
         coef2=half*(diff2(ie2-1)+diff2(ie2))
         bound=two*min(abs(diff2(ie2-1)),abs(diff2(ie2)))
         if (diff2(ie2)*diff2(ie2-1).gt.zero) then
           slope2(ic0,ic1,ie2)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(2) 
         else
            slope2(ic0,ic1,ie2)=zero
         endif
      enddo
         enddo
      enddo


      do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)


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

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
     &                +slope0(ic0,ic1,ic2)*deltax(ir0,0)
     &                +slope1(ic0,ic1,ic2)*deltax1
     &                +slope2(ic0,ic1,ic2)*deltax2
          enddo
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 3d side-centered float data
c***********************************************************************
c
      subroutine cartclinrefsideflot3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0,diff1,slope1,diff2,slope2)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      real
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2)
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,
     &        ir0,ir1,ir2
      real
     &  coef2,bound
      double precision deltax2,deltax1
c
c***********************************************************************
c


      do ir0=0,ratio(0)-1
         deltax(ir0,0)=dble(ir0)*dxf(0)
      enddo


      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo


      do ir2=0,ratio(2)-1
         deltax(ir2,2)=(dble(ir2)+half)*dxf(2)-dxc(2)*half
      enddo


      do ic2=ifirstc2,ilastc2

         do ic1=ifirstc1,ilastc1

      do ic0=ifirstc0-1,ilastc0+1
         diff0(ic0)=arrayc(ic0+1,ic1,ic2)
     &                -arrayc(ic0,ic1,ic2)
      enddo
      do ie0=ifirstc0,ilastc0+1
         coef2=half*(diff0(ie0-1)+diff0(ie0))
         bound=two*min(abs(diff0(ie0-1)),abs(diff0(ie0)))
         if (diff0(ie0)*diff0(ie0-1).gt.zero) then
           slope0(ie0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(0) 
         else
            slope0(ie0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


            do ic2=ifirstc2,ilastc2

            do ic0=ifirstc0,ilastc0+1

      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ic0,ie1,ic2)
     &               -arrayc(ic0,ie1-1,ic2)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


            do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0+1

      do ie2=ifirstc2,ilastc2+1
         diff2(ie2)=arrayc(ic0,ic1,ie2)
     &               -arrayc(ic0,ic1,ie2-1)
      enddo
      do ic2=ifirstc2,ilastc2
         coef2=half*(diff2(ic2+1)+diff2(ic2))
         bound=two*min(abs(diff2(ic2+1)),abs(diff2(ic2)))
         if (diff2(ic2)*diff2(ic2+1).gt.zero) then
            slope2(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(2)
         else
            slope2(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)


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
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
               ir0=if0-ic0*ratio(0)

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
     &                +slope0(ic0,ic1,ic2)*deltax(ir0,0)
     &                +slope1(ic0,ic1,ic2)*deltax1
     &                +slope2(ic0,ic1,ic2)*deltax2
          enddo
        enddo
      enddo
c
      return
      end
c
      subroutine cartclinrefsideflot3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff1,slope1,diff2,slope2,diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      real
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1,
     &          filo2:fihi2),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2)
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,
     &        ir0,ir1,ir2
      real
     &  coef2,bound
      double precision deltax2,deltax1
c
c***********************************************************************
c


      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo


      do ir1=0,ratio(1)-1
         deltax(ir1,1)=dble(ir1)*dxf(1)
      enddo


      do ir2=0,ratio(2)-1
         deltax(ir2,2)=(dble(ir2)+half)*dxf(2)-dxc(2)*half
      enddo


      do ic2=ifirstc2,ilastc2

         do ic1=ifirstc1,ilastc1+1

      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ie0,ic1,ic2)
     &               -arrayc(ie0-1,ic1,ic2)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


            do ic2=ifirstc2,ilastc2

            do ic0=ifirstc0,ilastc0

      do ic1=ifirstc1-1,ilastc1+1
         diff1(ic1)=arrayc(ic0,ic1+1,ic2)
     &                -arrayc(ic0,ic1,ic2)
      enddo
      do ie1=ifirstc1,ilastc1+1
         coef2=half*(diff1(ie1-1)+diff1(ie1))
         bound=two*min(abs(diff1(ie1-1)),abs(diff1(ie1)))
         if (diff1(ie1)*diff1(ie1-1).gt.zero) then
           slope1(ic0,ie1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(1) 
         else
            slope1(ic0,ie1,ic2)=zero
         endif
      enddo
         enddo
      enddo


            do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0

      do ie2=ifirstc2,ilastc2+1
         diff2(ie2)=arrayc(ic0,ic1,ie2)
     &               -arrayc(ic0,ic1,ie2-1)
      enddo
      do ic2=ifirstc2,ilastc2
         coef2=half*(diff2(ic2+1)+diff2(ic2))
         bound=two*min(abs(diff2(ic2+1)),abs(diff2(ic2)))
         if (diff2(ic2)*diff2(ic2+1).gt.zero) then
            slope2(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(2)
         else
            slope2(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)


         do if1=ifirstf1,ilastf1+1
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

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
     &                +slope0(ic0,ic1,ic2)*deltax(ir0,0)
     &                +slope1(ic0,ic1,ic2)*deltax1
     &                +slope2(ic0,ic1,ic2)*deltax2
          enddo
        enddo
      enddo
c
      return
      end
c
      subroutine cartclinrefsideflot3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff2,slope2,diff0,slope0,diff1,slope1)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  dxc(0:3-1),
     &  dxf(0:3-1),
     &  deltax(0:15,0:3-1)
      real
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2+1),
     &  diff2(cilo2:cihi2+1),
     &  slope2(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1),
     &  diff0(cilo0:cihi0+1),
     &  slope0(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1),
     &  diff1(cilo1:cihi1+1),
     &  slope1(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1)
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,
     &        ir0,ir1,ir2
      real
     &  coef2,bound
      double precision deltax2,deltax1
c
c***********************************************************************
c


      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo


      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo


      do ir2=0,ratio(2)-1
         deltax(ir2,2)=dble(ir2)*dxf(2)
      enddo


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1

      do ie0=ifirstc0,ilastc0+1
         diff0(ie0)=arrayc(ie0,ic1,ic2)
     &               -arrayc(ie0-1,ic1,ic2)
      enddo
      do ic0=ifirstc0,ilastc0
         coef2=half*(diff0(ic0+1)+diff0(ic0))
         bound=two*min(abs(diff0(ic0+1)),abs(diff0(ic0)))
         if (diff0(ic0)*diff0(ic0+1).gt.zero) then
            slope0(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(0)
         else
            slope0(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


            do ic2=ifirstc2,ilastc2+1

            do ic0=ifirstc0,ilastc0

      do ie1=ifirstc1,ilastc1+1
         diff1(ie1)=arrayc(ic0,ie1,ic2)
     &               -arrayc(ic0,ie1-1,ic2)
      enddo
      do ic1=ifirstc1,ilastc1
         coef2=half*(diff1(ic1+1)+diff1(ic1))
         bound=two*min(abs(diff1(ic1+1)),abs(diff1(ic1)))
         if (diff1(ic1)*diff1(ic1+1).gt.zero) then
            slope1(ic0,ic1,ic2)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc(1)
         else
            slope1(ic0,ic1,ic2)=zero
         endif
      enddo
         enddo
      enddo


            do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0

      do ic2=ifirstc2-1,ilastc2+1
         diff2(ic2)=arrayc(ic0,ic1,ic2+1)
     &                -arrayc(ic0,ic1,ic2)
      enddo
      do ie2=ifirstc2,ilastc2+1
         coef2=half*(diff2(ie2-1)+diff2(ie2))
         bound=two*min(abs(diff2(ie2-1)),abs(diff2(ie2)))
         if (diff2(ie2)*diff2(ie2-1).gt.zero) then
           slope2(ic0,ic1,ie2)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc(2) 
         else
            slope2(ic0,ic1,ie2)=zero
         endif
      enddo
         enddo
      enddo


      do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)


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

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
     &                +slope0(ic0,ic1,ic2)*deltax(ir0,0)
     &                +slope1(ic0,ic1,ic2)*deltax1
     &                +slope2(ic0,ic1,ic2)*deltax2
          enddo
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Conservative linear interpolation for 3d side-centered complex data
c***********************************************************************
c
c     subroutine cartclinrefsidecplx3d0(
ccart_clinref_op_side_3d(double complex,0,1,2)c
c      subroutine cartclinrefsidecplx3d1(
ccart_clinref_op_side_3d(double complex,1,2,0)c
c      subroutine cartclinrefsidecplx3d2(
ccart_clinref_op_side_3d(double complex,2,0,1)c
