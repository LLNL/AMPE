c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: FORTRAN routines for spatial refining of 1d patch data
c               on a regular Cartesian mesh.
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
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
c***********************************************************************
c Constant interpolation for 1d cell-centered double data
c***********************************************************************
c
      subroutine conrefcelldoub1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,
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
c Constant interpolation for 1d cell-centered float data
c***********************************************************************
c
      subroutine conrefcellflot1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,
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
c Constant interpolation for 1d cell-centered complex data
c***********************************************************************
c
      subroutine conrefcellcplx1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,
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
c Constant interpolation for 1d cell-centered integer data
c***********************************************************************
c
      subroutine conrefcellintg1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,
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
c Constant interpolation for 1d edge-centered double data
c***********************************************************************
c
      subroutine conrefedgedoub1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,
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
     &  arrayc(cilo0:cihi0),
     &  arrayf(filo0:fihi0)
      integer ie0,if0
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ie0=(if0+1)/ratio(0)-1
         else
            ie0=if0/ratio(0)
         endif
         arrayf(if0)=arrayc(ie0)
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 1d edge-centered float data
c***********************************************************************
c
      subroutine conrefedgeflot1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,
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
      real
     &  arrayc(cilo0:cihi0),
     &  arrayf(filo0:fihi0)
      integer ie0,if0
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ie0=(if0+1)/ratio(0)-1
         else
            ie0=if0/ratio(0)
         endif
         arrayf(if0)=arrayc(ie0)
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 1d edge-centered complex data
c***********************************************************************
c
      subroutine conrefedgecplx1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,
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
      double complex
     &  arrayc(cilo0:cihi0),
     &  arrayf(filo0:fihi0)
      integer ie0,if0
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ie0=(if0+1)/ratio(0)-1
         else
            ie0=if0/ratio(0)
         endif
         arrayf(if0)=arrayc(ie0)
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 1d edge-centered integer data
c***********************************************************************
c
      subroutine conrefedgeintg1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,
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
      integer
     &  arrayc(cilo0:cihi0),
     &  arrayf(filo0:fihi0)
      integer ie0,if0
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ie0=(if0+1)/ratio(0)-1
         else
            ie0=if0/ratio(0)
         endif
         arrayf(if0)=arrayc(ie0)
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 1d face-centered double data
c***********************************************************************
c
      subroutine conreffacedoub1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,
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
     &  arrayc(cilo0:cihi0+1),
     &  arrayf(filo0:fihi0+1)
      integer ie0,if0,it
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0+1
         it=2*if0+ratio(0)
         if (it.le.0) then
            ie0=it/(2*ratio(0))-1
         else
            ie0=(it-1)/(2*ratio(0))
         endif
         arrayf(if0)=arrayc(ie0)
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 1d face-centered float data
c***********************************************************************
c
      subroutine conreffaceflot1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,
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
      real
     &  arrayc(cilo0:cihi0+1),
     &  arrayf(filo0:fihi0+1)
      integer ie0,if0,it
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0+1
         it=2*if0+ratio(0)
         if (it.le.0) then
            ie0=it/(2*ratio(0))-1
         else
            ie0=(it-1)/(2*ratio(0))
         endif
         arrayf(if0)=arrayc(ie0)
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 1d face-centered complex data
c***********************************************************************
c
      subroutine conreffacecplx1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,
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
      double complex
     &  arrayc(cilo0:cihi0+1),
     &  arrayf(filo0:fihi0+1)
      integer ie0,if0,it
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0+1
         it=2*if0+ratio(0)
         if (it.le.0) then
            ie0=it/(2*ratio(0))-1
         else
            ie0=(it-1)/(2*ratio(0))
         endif
         arrayf(if0)=arrayc(ie0)
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 1d face-centered integer data
c***********************************************************************
c
      subroutine conreffaceintg1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,
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
      integer
     &  arrayc(cilo0:cihi0+1),
     &  arrayf(filo0:fihi0+1)
      integer ie0,if0,it
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0+1
         it=2*if0+ratio(0)
         if (it.le.0) then
            ie0=it/(2*ratio(0))-1
         else
            ie0=(it-1)/(2*ratio(0))
         endif
         arrayf(if0)=arrayc(ie0)
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 1d outerface double data
c***********************************************************************
c
      subroutine conrefoutfacedoub1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,
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
c Constant interpolation for 1d outerface float data
c***********************************************************************
c
      subroutine conrefoutfaceflot1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,
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
      real
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
c Constant interpolation for 1d outerface complex data
c***********************************************************************
c
      subroutine conrefoutfacecplx1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,
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
      double complex
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
c Constant interpolation for 1d outerface integer data
c***********************************************************************
c
      subroutine conrefoutfaceintg1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,
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
      integer
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
c Constant interpolation for 1d side-centered double data
c***********************************************************************
c
      subroutine conrefsidedoub1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,
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
     &  arrayc(cilo0:cihi0+1),
     &  arrayf(filo0:fihi0+1)
      integer ie0,if0,it
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0+1
         it=2*if0+ratio(0)
         if (it.le.0) then
            ie0=it/(2*ratio(0))-1
         else
            ie0=(it-1)/(2*ratio(0))
         endif
         arrayf(if0)=arrayc(ie0)
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 1d side-centered float data
c***********************************************************************
c
      subroutine conrefsideflot1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,
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
      real
     &  arrayc(cilo0:cihi0+1),
     &  arrayf(filo0:fihi0+1)
      integer ie0,if0,it
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0+1
         it=2*if0+ratio(0)
         if (it.le.0) then
            ie0=it/(2*ratio(0))-1
         else
            ie0=(it-1)/(2*ratio(0))
         endif
         arrayf(if0)=arrayc(ie0)
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 1d side-centered complex data
c***********************************************************************
c
      subroutine conrefsidecplx1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,
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
      double complex
     &  arrayc(cilo0:cihi0+1),
     &  arrayf(filo0:fihi0+1)
      integer ie0,if0,it
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0+1
         it=2*if0+ratio(0)
         if (it.le.0) then
            ie0=it/(2*ratio(0))-1
         else
            ie0=(it-1)/(2*ratio(0))
         endif
         arrayf(if0)=arrayc(ie0)
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 1d side-centered integer data
c***********************************************************************
c
      subroutine conrefsideintg1d(
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,
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
      integer
     &  arrayc(cilo0:cihi0+1),
     &  arrayf(filo0:fihi0+1)
      integer ie0,if0,it
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0+1
         it=2*if0+ratio(0)
         if (it.le.0) then
            ie0=it/(2*ratio(0))-1
         else
            ie0=(it-1)/(2*ratio(0))
         endif
         arrayf(if0)=arrayc(ie0)
      enddo
c
      return
      end
c
