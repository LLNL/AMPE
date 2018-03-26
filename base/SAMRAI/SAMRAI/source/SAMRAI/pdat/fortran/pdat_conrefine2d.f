c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: FORTRAN routines for spatial refining of 2d patch data
c               on a regular Cartesian mesh.
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for 2d constant refine operators
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
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for constant patchdata transfer routines.
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
c Constant interpolation for 2d cell-centered double data
c***********************************************************************
c
      subroutine conrefcelldoub2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 2d cell-centered float data
c***********************************************************************
c
      subroutine conrefcellflot2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      real
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 2d cell-centered complex data
c***********************************************************************
c
      subroutine conrefcellcplx2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      double complex
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 2d cell-centered integer data
c***********************************************************************
c
      subroutine conrefcellintg2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      integer
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 2d edge-centered double data
c***********************************************************************
c
      subroutine conrefedgedoub2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0+0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
      subroutine conrefedgedoub2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1+0
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 2d edge-centered float data
c***********************************************************************
c
      subroutine conrefedgeflot2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      real
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0+0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
      subroutine conrefedgeflot2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      real
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1+0
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 2d edge-centered complex data
c***********************************************************************

      subroutine conrefedgecplx2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      double complex
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0+0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
      subroutine conrefedgecplx2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      double complex
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1+0
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 2d edge-centered integer data
c***********************************************************************
c
      subroutine conrefedgeintg2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      integer
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0+0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
      subroutine conrefedgeintg2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      integer
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1+0
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 2d face-centered double data
c***********************************************************************
c
      subroutine conreffacedoub2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1)
      integer ie0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ie0=(if0+1)/ratio(0)-1
         else
            ie0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ie0,ic1)
         enddo
      enddo
c
      return
      end
c
      subroutine conreffacedoub2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
     &  arrayc(cilo1:cihi1+1,
     &          cilo0:cihi0),
     &  arrayf(filo1:fihi1+1,
     &          filo0:fihi0)
      integer ie1,ic0,if1,if0
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ie1=(if1+1)/ratio(1)-1
         else
            ie1=if1/ratio(1)
         endif
            arrayf(if1,if0)=arrayc(ie1,ic0)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 2d face-centered float data
c***********************************************************************
c
      subroutine conreffaceflot2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      real
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1)
      integer ie0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ie0=(if0+1)/ratio(0)-1
         else
            ie0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ie0,ic1)
         enddo
      enddo
c
      return
      end
c
      subroutine conreffaceflot2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      real
     &  arrayc(cilo1:cihi1+1,
     &          cilo0:cihi0),
     &  arrayf(filo1:fihi1+1,
     &          filo0:fihi0)
      integer ie1,ic0,if1,if0
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ie1=(if1+1)/ratio(1)-1
         else
            ie1=if1/ratio(1)
         endif
            arrayf(if1,if0)=arrayc(ie1,ic0)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 2d face-centered complex data
c***********************************************************************

      subroutine conreffacecplx2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      double complex
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1)
      integer ie0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ie0=(if0+1)/ratio(0)-1
         else
            ie0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ie0,ic1)
         enddo
      enddo
c
      return
      end
c
      subroutine conreffacecplx2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      double complex
     &  arrayc(cilo1:cihi1+1,
     &          cilo0:cihi0),
     &  arrayf(filo1:fihi1+1,
     &          filo0:fihi0)
      integer ie1,ic0,if1,if0
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ie1=(if1+1)/ratio(1)-1
         else
            ie1=if1/ratio(1)
         endif
            arrayf(if1,if0)=arrayc(ie1,ic0)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 2d face-centered integer data
c***********************************************************************
c
      subroutine conreffaceintg2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      integer
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1)
      integer ie0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ie0=(if0+1)/ratio(0)-1
         else
            ie0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ie0,ic1)
         enddo
      enddo
c
      return
      end
c
      subroutine conreffaceintg2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      integer
     &  arrayc(cilo1:cihi1+1,
     &          cilo0:cihi0),
     &  arrayf(filo1:fihi1+1,
     &          filo0:fihi0)
      integer ie1,ic0,if1,if0
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ie1=(if1+1)/ratio(1)-1
         else
            ie1=if1/ratio(1)
         endif
            arrayf(if1,if0)=arrayc(ie1,ic0)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 2d outerface double data
c***********************************************************************
c
      subroutine conrefoutfacedoub2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
     &  arrayc(cilo1:cihi1),
     &  arrayf(filo1:fihi1)
      integer ic1,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         arrayf(if1)=arrayc(ic1)
      enddo
c
      return
      end
c
      subroutine conrefoutfacedoub2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
     &  arrayc(cilo0:cihi0),
     &  arrayf(filo0:fihi0)
      integer ic0,if0
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         arrayf(if0)=arrayc(ic0)
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 2d outerface float data
c***********************************************************************
c
      subroutine conrefoutfaceflot2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      real
     &  arrayc(cilo1:cihi1),
     &  arrayf(filo1:fihi1)
      integer ic1,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         arrayf(if1)=arrayc(ic1)
      enddo
c
      return
      end
c
      subroutine conrefoutfaceflot2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      real
     &  arrayc(cilo0:cihi0),
     &  arrayf(filo0:fihi0)
      integer ic0,if0
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         arrayf(if0)=arrayc(ic0)
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 2d outerface complex data
c***********************************************************************

      subroutine conrefoutfacecplx2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      double complex
     &  arrayc(cilo1:cihi1),
     &  arrayf(filo1:fihi1)
      integer ic1,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         arrayf(if1)=arrayc(ic1)
      enddo
c
      return
      end
c
      subroutine conrefoutfacecplx2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      double complex
     &  arrayc(cilo0:cihi0),
     &  arrayf(filo0:fihi0)
      integer ic0,if0
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         arrayf(if0)=arrayc(ic0)
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 2d outerface integer data
c***********************************************************************
c
      subroutine conrefoutfaceintg2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      integer
     &  arrayc(cilo1:cihi1),
     &  arrayf(filo1:fihi1)
      integer ic1,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         arrayf(if1)=arrayc(ic1)
      enddo
c
      return
      end
c
      subroutine conrefoutfaceintg2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      integer
     &  arrayc(cilo0:cihi0),
     &  arrayf(filo0:fihi0)
      integer ic0,if0
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         arrayf(if0)=arrayc(ic0)
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 2d side-centered double data
c***********************************************************************
c
      subroutine conrefsidedoub2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1+0
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
      subroutine conrefsidedoub2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0+0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 2d side-centered float data
c***********************************************************************
c
      subroutine conrefsideflot2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      real
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1+0
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
      subroutine conrefsideflot2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      real
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0+0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 2d side-centered complex data
c***********************************************************************

      subroutine conrefsidecplx2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      double complex
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1+0
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
      subroutine conrefsidecplx2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      double complex
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0+0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 2d side-centered integer data
c***********************************************************************
c
      subroutine conrefsideintg2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      integer
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1+0
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
      subroutine conrefsideintg2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
      integer
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0+0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
