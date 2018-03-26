c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: FORTRAN routines for standard linear time interpolation 
c               of 2d patch data types.
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for 2d std linear time interpolation operators.
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for dimensioning 2d arrays in FORTRAN routines.
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
c Linear time interpolation for 2d cell-centered double data
c***********************************************************************
c
      subroutine lintimeintcelldoub2d(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0,
     &          oilo1:oihi1),
     &  arraynew(nilo0:nihi0,
     &          nilo1:nihi1),
     &  arraydst(dilo0:dihi0,
     &          dilo1:dihi1)
      integer ic0,ic1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic1=ifirst1,ilast1
         do ic0=ifirst0,ilast0
            arraydst(ic0,ic1)=arrayold(ic0,ic1)*oldfrac
     &                       +arraynew(ic0,ic1)*tfrac
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 2d cell-centered float data
c***********************************************************************
c
      subroutine lintimeintcellfloat2d(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0,
     &          oilo1:oihi1),
     &  arraynew(nilo0:nihi0,
     &          nilo1:nihi1),
     &  arraydst(dilo0:dihi0,
     &          dilo1:dihi1)
      integer ic0,ic1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic1=ifirst1,ilast1
         do ic0=ifirst0,ilast0
            arraydst(ic0,ic1)=arrayold(ic0,ic1)*oldfrac
     &                       +arraynew(ic0,ic1)*tfrac
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 2d cell-centered complex data
c***********************************************************************
c
      subroutine lintimeintcellcmplx2d(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0,
     &          oilo1:oihi1),
     &  arraynew(nilo0:nihi0,
     &          nilo1:nihi1),
     &  arraydst(dilo0:dihi0,
     &          dilo1:dihi1)
      integer ic0,ic1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic1=ifirst1,ilast1
         do ic0=ifirst0,ilast0
            arraydst(ic0,ic1)=arrayold(ic0,ic1)*oldfrac
     &                       +arraynew(ic0,ic1)*tfrac
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 2d edge-centered double data
c***********************************************************************
c
      subroutine lintimeintedgedoub2d0(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0,
     &          oilo1:oihi1+1),
     &  arraynew(nilo0:nihi0,
     &          nilo1:nihi1+1),
     &  arraydst(dilo0:dihi0,
     &          dilo1:dihi1+1)
      integer ie0,ie1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie1=ifirst1,ilast1+1
         do ie0=ifirst0,ilast0+0
            arraydst(ie0,ie1)=arrayold(ie0,ie1)*oldfrac
     &                       +arraynew(ie0,ie1)*tfrac
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintedgedoub2d1(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1)
      integer ie0,ie1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie1=ifirst1,ilast1+0
         do ie0=ifirst0,ilast0+1
            arraydst(ie0,ie1)=arrayold(ie0,ie1)*oldfrac
     &                       +arraynew(ie0,ie1)*tfrac
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 2d edge-centered float data
c***********************************************************************
c
      subroutine lintimeintedgefloat2d0(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0,
     &          oilo1:oihi1+1),
     &  arraynew(nilo0:nihi0,
     &          nilo1:nihi1+1),
     &  arraydst(dilo0:dihi0,
     &          dilo1:dihi1+1)
      integer ie0,ie1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie1=ifirst1,ilast1+1
         do ie0=ifirst0,ilast0+0
            arraydst(ie0,ie1)=arrayold(ie0,ie1)*oldfrac
     &                       +arraynew(ie0,ie1)*tfrac
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintedgefloat2d1(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1)
      integer ie0,ie1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie1=ifirst1,ilast1+0
         do ie0=ifirst0,ilast0+1
            arraydst(ie0,ie1)=arrayold(ie0,ie1)*oldfrac
     &                       +arraynew(ie0,ie1)*tfrac
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 2d edge-centered complex data
c***********************************************************************
c
      subroutine lintimeintedgecmplx2d0(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0,
     &          oilo1:oihi1+1),
     &  arraynew(nilo0:nihi0,
     &          nilo1:nihi1+1),
     &  arraydst(dilo0:dihi0,
     &          dilo1:dihi1+1)
      integer ie0,ie1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie1=ifirst1,ilast1+1
         do ie0=ifirst0,ilast0+0
            arraydst(ie0,ie1)=arrayold(ie0,ie1)*oldfrac
     &                       +arraynew(ie0,ie1)*tfrac
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintedgecmplx2d1(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1)
      integer ie0,ie1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie1=ifirst1,ilast1+0
         do ie0=ifirst0,ilast0+1
            arraydst(ie0,ie1)=arrayold(ie0,ie1)*oldfrac
     &                       +arraynew(ie0,ie1)*tfrac
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 2d face-centered double data
c***********************************************************************
c
      subroutine lintimeintfacedoub2d0(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1)
      integer ie0,ic1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic1=ifirst1,ilast1
         do ie0=ifirst0,ilast0+1
            arraydst(ie0,ic1)=arrayold(ie0,ic1)*oldfrac
     &                       +arraynew(ie0,ic1)*tfrac
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintfacedoub2d1(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo1:oihi1+1,
     &          oilo0:oihi0),
     &  arraynew(nilo1:nihi1+1,
     &          nilo0:nihi0),
     &  arraydst(dilo1:dihi1+1,
     &          dilo0:dihi0)
      integer ie1,ic0
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic0=ifirst0,ilast0
         do ie1=ifirst1,ilast1+1
            arraydst(ie1,ic0)=arrayold(ie1,ic0)*oldfrac
     &                       +arraynew(ie1,ic0)*tfrac
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 2d face-centered float data
c***********************************************************************
c
      subroutine lintimeintfacefloat2d0(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1)
      integer ie0,ic1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic1=ifirst1,ilast1
         do ie0=ifirst0,ilast0+1
            arraydst(ie0,ic1)=arrayold(ie0,ic1)*oldfrac
     &                       +arraynew(ie0,ic1)*tfrac
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintfacefloat2d1(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo1:oihi1+1,
     &          oilo0:oihi0),
     &  arraynew(nilo1:nihi1+1,
     &          nilo0:nihi0),
     &  arraydst(dilo1:dihi1+1,
     &          dilo0:dihi0)
      integer ie1,ic0
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic0=ifirst0,ilast0
         do ie1=ifirst1,ilast1+1
            arraydst(ie1,ic0)=arrayold(ie1,ic0)*oldfrac
     &                       +arraynew(ie1,ic0)*tfrac
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 2d face-centered complex data
c***********************************************************************
c
      subroutine lintimeintfacecmplx2d0(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1)
      integer ie0,ic1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic1=ifirst1,ilast1
         do ie0=ifirst0,ilast0+1
            arraydst(ie0,ic1)=arrayold(ie0,ic1)*oldfrac
     &                       +arraynew(ie0,ic1)*tfrac
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintfacecmplx2d1(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo1:oihi1+1,
     &          oilo0:oihi0),
     &  arraynew(nilo1:nihi1+1,
     &          nilo0:nihi0),
     &  arraydst(dilo1:dihi1+1,
     &          dilo0:dihi0)
      integer ie1,ic0
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic0=ifirst0,ilast0
         do ie1=ifirst1,ilast1+1
            arraydst(ie1,ic0)=arrayold(ie1,ic0)*oldfrac
     &                       +arraynew(ie1,ic0)*tfrac
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 2d node-centered double data
c***********************************************************************
c
      subroutine lintimeintnodedoub2d(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1+1),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1+1),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1+1)
      integer ie0,ie1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie1=ifirst1,ilast1+1
         do ie0=ifirst0,ilast0+1
            arraydst(ie0,ie1)=arrayold(ie0,ie1)*oldfrac
     &                       +arraynew(ie0,ie1)*tfrac
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 2d node-centered float data
c***********************************************************************
c
      subroutine lintimeintnodefloat2d(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1+1),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1+1),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1+1)
      integer ie0,ie1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie1=ifirst1,ilast1+1
         do ie0=ifirst0,ilast0+1
            arraydst(ie0,ie1)=arrayold(ie0,ie1)*oldfrac
     &                       +arraynew(ie0,ie1)*tfrac
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 2d node-centered complex data
c***********************************************************************
c
      subroutine lintimeintnodecmplx2d(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1+1),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1+1),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1+1)
      integer ie0,ie1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie1=ifirst1,ilast1+1
         do ie0=ifirst0,ilast0+1
            arraydst(ie0,ie1)=arrayold(ie0,ie1)*oldfrac
     &                       +arraynew(ie0,ie1)*tfrac
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 2d outerface double data
c***********************************************************************
c
      subroutine lintimeintoutfacedoub2d0(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo1:oihi1),
     &  arraynew(nilo1:nihi1),
     &  arraydst(dilo1:dihi1)
      integer ic1 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic1=ifirst1,ilast1
         arraydst(ic1)=arrayold(ic1)*oldfrac
     &                +arraynew(ic1)*tfrac
      enddo
c
      return
      end
c
      subroutine lintimeintoutfacedoub2d1(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0),
     &  arraynew(nilo0:nihi0),
     &  arraydst(dilo0:dihi0)
      integer ic0 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic0=ifirst0,ilast0
         arraydst(ic0)=arrayold(ic0)*oldfrac
     &                +arraynew(ic0)*tfrac
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 2d outerface float data
c***********************************************************************
c
      subroutine lintimeintoutfacefloat2d0(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo1:oihi1),
     &  arraynew(nilo1:nihi1),
     &  arraydst(dilo1:dihi1)
      integer ic1 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic1=ifirst1,ilast1
         arraydst(ic1)=arrayold(ic1)*oldfrac
     &                +arraynew(ic1)*tfrac
      enddo
c
      return
      end
c
      subroutine lintimeintoutfacefloat2d1(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0),
     &  arraynew(nilo0:nihi0),
     &  arraydst(dilo0:dihi0)
      integer ic0 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic0=ifirst0,ilast0
         arraydst(ic0)=arrayold(ic0)*oldfrac
     &                +arraynew(ic0)*tfrac
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 2d outerface complex data
c***********************************************************************
c
      subroutine lintimeintoutfacecmplx2d0(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo1:oihi1),
     &  arraynew(nilo1:nihi1),
     &  arraydst(dilo1:dihi1)
      integer ic1 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic1=ifirst1,ilast1
         arraydst(ic1)=arrayold(ic1)*oldfrac
     &                +arraynew(ic1)*tfrac
      enddo
c
      return
      end
c
      subroutine lintimeintoutfacecmplx2d1(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0),
     &  arraynew(nilo0:nihi0),
     &  arraydst(dilo0:dihi0)
      integer ic0 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic0=ifirst0,ilast0
         arraydst(ic0)=arrayold(ic0)*oldfrac
     &                +arraynew(ic0)*tfrac
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 2d outerside double data
c***********************************************************************
c
      subroutine lintimeintoutsidedoub2d0(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo1:oihi1),
     &  arraynew(nilo1:nihi1),
     &  arraydst(dilo1:dihi1)
      integer ic1 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic1=ifirst1,ilast1
         arraydst(ic1)=arrayold(ic1)*oldfrac
     &                +arraynew(ic1)*tfrac
      enddo
c
      return
      end
c
      subroutine lintimeintoutsidedoub2d1(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0),
     &  arraynew(nilo0:nihi0),
     &  arraydst(dilo0:dihi0)
      integer ic0 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic0=ifirst0,ilast0
         arraydst(ic0)=arrayold(ic0)*oldfrac
     &                +arraynew(ic0)*tfrac
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 2d outerside float data
c***********************************************************************
c
      subroutine lintimeintoutsidefloat2d0(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo1:oihi1),
     &  arraynew(nilo1:nihi1),
     &  arraydst(dilo1:dihi1)
      integer ic1 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic1=ifirst1,ilast1
         arraydst(ic1)=arrayold(ic1)*oldfrac
     &                +arraynew(ic1)*tfrac
      enddo
c
      return
      end
c
      subroutine lintimeintoutsidefloat2d1(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0),
     &  arraynew(nilo0:nihi0),
     &  arraydst(dilo0:dihi0)
      integer ic0 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic0=ifirst0,ilast0
         arraydst(ic0)=arrayold(ic0)*oldfrac
     &                +arraynew(ic0)*tfrac
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 2d outerside complex data
c***********************************************************************
c
      subroutine lintimeintoutsidecmplx2d0(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo1:oihi1),
     &  arraynew(nilo1:nihi1),
     &  arraydst(dilo1:dihi1)
      integer ic1 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic1=ifirst1,ilast1
         arraydst(ic1)=arrayold(ic1)*oldfrac
     &                +arraynew(ic1)*tfrac
      enddo
c
      return
      end
c
      subroutine lintimeintoutsidecmplx2d1(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0),
     &  arraynew(nilo0:nihi0),
     &  arraydst(dilo0:dihi0)
      integer ic0 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic0=ifirst0,ilast0
         arraydst(ic0)=arrayold(ic0)*oldfrac
     &                +arraynew(ic0)*tfrac
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 2d side-centered double data
c***********************************************************************
c
      subroutine lintimeintsidedoub2d0(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1)
      integer ie0,ie1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie1=ifirst1,ilast1+0
         do ie0=ifirst0,ilast0+1
            arraydst(ie0,ie1)=arrayold(ie0,ie1)*oldfrac
     &                       +arraynew(ie0,ie1)*tfrac
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintsidedoub2d1(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0,
     &          oilo1:oihi1+1),
     &  arraynew(nilo0:nihi0,
     &          nilo1:nihi1+1),
     &  arraydst(dilo0:dihi0,
     &          dilo1:dihi1+1)
      integer ie0,ie1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie1=ifirst1,ilast1+1
         do ie0=ifirst0,ilast0+0
            arraydst(ie0,ie1)=arrayold(ie0,ie1)*oldfrac
     &                       +arraynew(ie0,ie1)*tfrac
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 2d side-centered float data
c***********************************************************************
c
      subroutine lintimeintsidefloat2d0(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1)
      integer ie0,ie1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie1=ifirst1,ilast1+0
         do ie0=ifirst0,ilast0+1
            arraydst(ie0,ie1)=arrayold(ie0,ie1)*oldfrac
     &                       +arraynew(ie0,ie1)*tfrac
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintsidefloat2d1(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0,
     &          oilo1:oihi1+1),
     &  arraynew(nilo0:nihi0,
     &          nilo1:nihi1+1),
     &  arraydst(dilo0:dihi0,
     &          dilo1:dihi1+1)
      integer ie0,ie1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie1=ifirst1,ilast1+1
         do ie0=ifirst0,ilast0+0
            arraydst(ie0,ie1)=arrayold(ie0,ie1)*oldfrac
     &                       +arraynew(ie0,ie1)*tfrac
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 2d side-centered complex data
c***********************************************************************
c
      subroutine lintimeintsidecmplx2d0(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1)
      integer ie0,ie1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie1=ifirst1,ilast1+0
         do ie0=ifirst0,ilast0+1
            arraydst(ie0,ie1)=arrayold(ie0,ie1)*oldfrac
     &                       +arraynew(ie0,ie1)*tfrac
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintsidecmplx2d1(
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  oilo0,oilo1,oihi0,oihi1,
     &  nilo0,nilo1,nihi0,nihi1,
     &  dilo0,dilo1,dihi0,dihi1
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0,
     &          oilo1:oihi1+1),
     &  arraynew(nilo0:nihi0,
     &          nilo1:nihi1+1),
     &  arraydst(dilo0:dihi0,
     &          dilo1:dihi1+1)
      integer ie0,ie1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie1=ifirst1,ilast1+1
         do ie0=ifirst0,ilast0+0
            arraydst(ie0,ie1)=arrayold(ie0,ie1)*oldfrac
     &                       +arraynew(ie0,ie1)*tfrac
         enddo
      enddo
c
      return
      end
c
