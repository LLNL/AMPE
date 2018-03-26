c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: FORTRAN routines for standard linear time interpolation 
c               of 1d patch data types.
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for standard 1d time interpolation operators.
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for dimensioning 1d arrays in FORTRAN routines.
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
c Linear time interpolation for 1d cell-centered double data
c***********************************************************************
c
      subroutine lintimeintcelldoub1d(
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
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
c Linear time interpolation for 1d cell-centered float data
c***********************************************************************
c
      subroutine lintimeintcellfloat1d(
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
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
c Linear time interpolation for 1d cell-centered complex data
c***********************************************************************
c
      subroutine lintimeintcellcmplx1d(
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
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
c Linear time interpolation for 1d edge-centered double data
c***********************************************************************
c
      subroutine lintimeintedgedoub1d(
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0),
     &  arraynew(nilo0:nihi0),
     &  arraydst(dilo0:dihi0)
      integer ie0
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie0=ifirst0,ilast0
         arraydst(ie0)=arrayold(ie0)*oldfrac
     &                +arraynew(ie0)*tfrac
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 1d edge-centered float data
c***********************************************************************
c
      subroutine lintimeintedgefloat1d(
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0),
     &  arraynew(nilo0:nihi0),
     &  arraydst(dilo0:dihi0)
      integer ie0
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie0=ifirst0,ilast0
         arraydst(ie0)=arrayold(ie0)*oldfrac
     &                +arraynew(ie0)*tfrac
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 1d edge-centered complex data
c***********************************************************************
c
      subroutine lintimeintedgecmplx1d(
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0),
     &  arraynew(nilo0:nihi0),
     &  arraydst(dilo0:dihi0)
      integer ie0
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie0=ifirst0,ilast0
         arraydst(ie0)=arrayold(ie0)*oldfrac
     &                +arraynew(ie0)*tfrac
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 1d face-centered double data
c***********************************************************************
c
      subroutine lintimeintfacedoub1d(
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0+1),
     &  arraynew(nilo0:nihi0+1),
     &  arraydst(dilo0:dihi0+1)
      integer ie0
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie0=ifirst0,ilast0+1
         arraydst(ie0)=arrayold(ie0)*oldfrac
     &                +arraynew(ie0)*tfrac
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 1d face-centered float data
c***********************************************************************
c
      subroutine lintimeintfacefloat1d(
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0+1),
     &  arraynew(nilo0:nihi0+1),
     &  arraydst(dilo0:dihi0+1)
      integer ie0
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie0=ifirst0,ilast0+1
         arraydst(ie0)=arrayold(ie0)*oldfrac
     &                +arraynew(ie0)*tfrac
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 1d face-centered complex data
c***********************************************************************
c
      subroutine lintimeintfacecmplx1d(
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0+1),
     &  arraynew(nilo0:nihi0+1),
     &  arraydst(dilo0:dihi0+1)
      integer ie0
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie0=ifirst0,ilast0+1
         arraydst(ie0)=arrayold(ie0)*oldfrac
     &                +arraynew(ie0)*tfrac
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 1d node-centered double data
c***********************************************************************
c
      subroutine lintimeintnodedoub1d(
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0+1),
     &  arraynew(nilo0:nihi0+1),
     &  arraydst(dilo0:dihi0+1)
      integer ie0
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie0=ifirst0,ilast0+1
         arraydst(ie0)=arrayold(ie0)*oldfrac
     &                +arraynew(ie0)*tfrac
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 1d node-centered float data
c***********************************************************************
c
      subroutine lintimeintnodefloat1d(
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0+1),
     &  arraynew(nilo0:nihi0+1),
     &  arraydst(dilo0:dihi0+1)
      integer ie0
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie0=ifirst0,ilast0+1
         arraydst(ie0)=arrayold(ie0)*oldfrac
     &                +arraynew(ie0)*tfrac
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 1d node-centered complex data
c***********************************************************************
c
      subroutine lintimeintnodecmplx1d(
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0+1),
     &  arraynew(nilo0:nihi0+1),
     &  arraydst(dilo0:dihi0+1)
      integer ie0
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie0=ifirst0,ilast0+1
         arraydst(ie0)=arrayold(ie0)*oldfrac
     &                +arraynew(ie0)*tfrac
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 1d outerface double data
c***********************************************************************
c
      subroutine lintimeintoutfacedoub1d(
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(1),
     &  arraynew(1),
     &  arraydst(1)
c
c***********************************************************************
c
      oldfrac=one-tfrac

      arraydst(1)=arrayold(1)*oldfrac
     &           +arraynew(1)*tfrac
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 1d outerface float data
c***********************************************************************
c
      subroutine lintimeintoutfacefloat1d(
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(1),
     &  arraynew(1),
     &  arraydst(1)
c
c***********************************************************************
c
      oldfrac=one-tfrac

      arraydst(1)=arrayold(1)*oldfrac
     &           +arraynew(1)*tfrac
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 1d outerface complex data
c***********************************************************************
c
      subroutine lintimeintoutfacecmplx1d(
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(1),
     &  arraynew(1),
     &  arraydst(1)
c
c***********************************************************************
c
      oldfrac=one-tfrac

      arraydst(1)=arrayold(1)*oldfrac
     &           +arraynew(1)*tfrac
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 1d outerside double data
c***********************************************************************
c
      subroutine lintimeintoutsidedoub1d(
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(1),
     &  arraynew(1),
     &  arraydst(1)
c
c***********************************************************************
c
      oldfrac=one-tfrac

      arraydst(1)=arrayold(1)*oldfrac
     &           +arraynew(1)*tfrac
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 1d outerside float data
c***********************************************************************
c
      subroutine lintimeintoutsidefloat1d(
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(1),
     &  arraynew(1),
     &  arraydst(1)
c
c***********************************************************************
c
      oldfrac=one-tfrac

      arraydst(1)=arrayold(1)*oldfrac
     &           +arraynew(1)*tfrac
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 1d outerside complex data
c***********************************************************************
c
      subroutine lintimeintoutsidecmplx1d(
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(1),
     &  arraynew(1),
     &  arraydst(1)
c
c***********************************************************************
c
      oldfrac=one-tfrac

      arraydst(1)=arrayold(1)*oldfrac
     &           +arraynew(1)*tfrac
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 1d side-centered double data
c***********************************************************************
c
      subroutine lintimeintsidedoub1d(
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0+1),
     &  arraynew(nilo0:nihi0+1),
     &  arraydst(dilo0:dihi0+1)
      integer ie0
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie0=ifirst0,ilast0+1
         arraydst(ie0)=arrayold(ie0)*oldfrac
     &                +arraynew(ie0)*tfrac
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 1d side-centered float data
c***********************************************************************
c
      subroutine lintimeintsidefloat1d(
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0+1),
     &  arraynew(nilo0:nihi0+1),
     &  arraydst(dilo0:dihi0+1)
      integer ie0
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie0=ifirst0,ilast0+1
         arraydst(ie0)=arrayold(ie0)*oldfrac
     &                +arraynew(ie0)*tfrac
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 1d side-centered complex data
c***********************************************************************
c
      subroutine lintimeintsidecmplx1d(
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0+1),
     &  arraynew(nilo0:nihi0+1),
     &  arraydst(dilo0:dihi0+1)
      integer ie0
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie0=ifirst0,ilast0+1
         arraydst(ie0)=arrayold(ie0)*oldfrac
     &                +arraynew(ie0)*tfrac
      enddo
c
      return
      end
c
