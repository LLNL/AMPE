c
c  File:        $URL$
c  Package:     SAMRAI gerometry
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: FORTRAN routines for standard linear time interpolation 
c               of 3d patch data types.
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for 3d std linear time interpolation operators.
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for dimensioning 3d arrays in FORTRAN routines.
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
c Linear time interpolation for 3d cell-centered double data
c***********************************************************************
c
      subroutine lintimeintcelldoub3d(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0,
     &          oilo1:oihi1,
     &          oilo2:oihi2),
     &  arraynew(nilo0:nihi0,
     &          nilo1:nihi1,
     &          nilo2:nihi2),
     &  arraydst(dilo0:dihi0,
     &          dilo1:dihi1,
     &          dilo2:dihi2)
      integer ic0,ic1,ic2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            do ic0=ifirst0,ilast0
               arraydst(ic0,ic1,ic2)=
     &                              +arrayold(ic0,ic1,ic2)*oldfrac
     &                              +arraynew(ic0,ic1,ic2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 3d cell-centered float data
c***********************************************************************
c
      subroutine lintimeintcellfloat3d(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0,
     &          oilo1:oihi1,
     &          oilo2:oihi2),
     &  arraynew(nilo0:nihi0,
     &          nilo1:nihi1,
     &          nilo2:nihi2),
     &  arraydst(dilo0:dihi0,
     &          dilo1:dihi1,
     &          dilo2:dihi2)
      integer ic0,ic1,ic2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            do ic0=ifirst0,ilast0
               arraydst(ic0,ic1,ic2)=
     &                              +arrayold(ic0,ic1,ic2)*oldfrac
     &                              +arraynew(ic0,ic1,ic2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 3d cell-centered complex data
c***********************************************************************
c
      subroutine lintimeintcellcmplx3d(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0,
     &          oilo1:oihi1,
     &          oilo2:oihi2),
     &  arraynew(nilo0:nihi0,
     &          nilo1:nihi1,
     &          nilo2:nihi2),
     &  arraydst(dilo0:dihi0,
     &          dilo1:dihi1,
     &          dilo2:dihi2)
      integer ic0,ic1,ic2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            do ic0=ifirst0,ilast0
               arraydst(ic0,ic1,ic2)=
     &                              +arrayold(ic0,ic1,ic2)*oldfrac
     &                              +arraynew(ic0,ic1,ic2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 3d edge-centered double data
c***********************************************************************
c
      subroutine lintimeintedgedoub3d0(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0,
     &          oilo1:oihi1+1,
     &          oilo2:oihi2+1),
     &  arraynew(nilo0:nihi0,
     &          nilo1:nihi1+1,
     &          nilo2:nihi2+1),
     &  arraydst(dilo0:dihi0,
     &          dilo1:dihi1+1,
     &          dilo2:dihi2+1)
      integer ie0,ie1,ie2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie2=ifirst2,ilast2+1
         do ie1=ifirst1,ilast1+1
            do ie0=ifirst0,ilast0
               arraydst(ie0,ie1,ie2)=
     &                              +arrayold(ie0,ie1,ie2)*oldfrac
     &                              +arraynew(ie0,ie1,ie2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintedgedoub3d1(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1,
     &          oilo2:oihi2+1),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1,
     &          nilo2:nihi2+1),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1,
     &          dilo2:dihi2+1)
      integer ie0,ie1,ie2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie2=ifirst2,ilast2+1
         do ie1=ifirst1,ilast1
            do ie0=ifirst0,ilast0+1
               arraydst(ie0,ie1,ie2)=
     &                              +arrayold(ie0,ie1,ie2)*oldfrac
     &                              +arraynew(ie0,ie1,ie2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintedgedoub3d2(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1+1,
     &          oilo2:oihi2),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1+1,
     &          nilo2:nihi2),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1+1,
     &          dilo2:dihi2)
      integer ie0,ie1,ie2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie2=ifirst2,ilast2
         do ie1=ifirst1,ilast1+1
            do ie0=ifirst0,ilast0+1
               arraydst(ie0,ie1,ie2)=
     &                              +arrayold(ie0,ie1,ie2)*oldfrac
     &                              +arraynew(ie0,ie1,ie2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 3d edge-centered float data
c***********************************************************************
c
      subroutine lintimeintedgefloat3d0(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0,
     &          oilo1:oihi1+1,
     &          oilo2:oihi2+1),
     &  arraynew(nilo0:nihi0,
     &          nilo1:nihi1+1,
     &          nilo2:nihi2+1),
     &  arraydst(dilo0:dihi0,
     &          dilo1:dihi1+1,
     &          dilo2:dihi2+1)
      integer ie0,ie1,ie2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie2=ifirst2,ilast2+1
         do ie1=ifirst1,ilast1+1
            do ie0=ifirst0,ilast0
               arraydst(ie0,ie1,ie2)=
     &                              +arrayold(ie0,ie1,ie2)*oldfrac
     &                              +arraynew(ie0,ie1,ie2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintedgefloat3d1(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1,
     &          oilo2:oihi2+1),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1,
     &          nilo2:nihi2+1),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1,
     &          dilo2:dihi2+1)
      integer ie0,ie1,ie2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie2=ifirst2,ilast2+1
         do ie1=ifirst1,ilast1
            do ie0=ifirst0,ilast0+1
               arraydst(ie0,ie1,ie2)=
     &                              +arrayold(ie0,ie1,ie2)*oldfrac
     &                              +arraynew(ie0,ie1,ie2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintedgefloat3d2(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1+1,
     &          oilo2:oihi2),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1+1,
     &          nilo2:nihi2),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1+1,
     &          dilo2:dihi2)
      integer ie0,ie1,ie2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie2=ifirst2,ilast2
         do ie1=ifirst1,ilast1+1
            do ie0=ifirst0,ilast0+1
               arraydst(ie0,ie1,ie2)=
     &                              +arrayold(ie0,ie1,ie2)*oldfrac
     &                              +arraynew(ie0,ie1,ie2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 3d edge-centered complex data
c***********************************************************************
c
      subroutine lintimeintedgecmplx3d0(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0,
     &          oilo1:oihi1+1,
     &          oilo2:oihi2+1),
     &  arraynew(nilo0:nihi0,
     &          nilo1:nihi1+1,
     &          nilo2:nihi2+1),
     &  arraydst(dilo0:dihi0,
     &          dilo1:dihi1+1,
     &          dilo2:dihi2+1)
      integer ie0,ie1,ie2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie2=ifirst2,ilast2+1
         do ie1=ifirst1,ilast1+1
            do ie0=ifirst0,ilast0
               arraydst(ie0,ie1,ie2)=
     &                              +arrayold(ie0,ie1,ie2)*oldfrac
     &                              +arraynew(ie0,ie1,ie2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintedgecmplx3d1(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1,
     &          oilo2:oihi2+1),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1,
     &          nilo2:nihi2+1),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1,
     &          dilo2:dihi2+1)
      integer ie0,ie1,ie2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie2=ifirst2,ilast2+1
         do ie1=ifirst1,ilast1
            do ie0=ifirst0,ilast0+1
               arraydst(ie0,ie1,ie2)=
     &                              +arrayold(ie0,ie1,ie2)*oldfrac
     &                              +arraynew(ie0,ie1,ie2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintedgecmplx3d2(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1+1,
     &          oilo2:oihi2),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1+1,
     &          nilo2:nihi2),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1+1,
     &          dilo2:dihi2)
      integer ie0,ie1,ie2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie2=ifirst2,ilast2
         do ie1=ifirst1,ilast1+1
            do ie0=ifirst0,ilast0+1
               arraydst(ie0,ie1,ie2)=
     &                              +arrayold(ie0,ie1,ie2)*oldfrac
     &                              +arraynew(ie0,ie1,ie2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 3d face-centered double data
c***********************************************************************
c
      subroutine lintimeintfacedoub3d0(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1,
     &          oilo2:oihi2),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1,
     &          nilo2:nihi2),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1,
     &          dilo2:dihi2)
      integer ie0,ic1,ic2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            do ie0=ifirst0,ilast0+1
               arraydst(ie0,ic1,ic2)=
     &                              +arrayold(ie0,ic1,ic2)*oldfrac
     &                              +arraynew(ie0,ic1,ic2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintfacedoub3d1(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo1:oihi1+1,
     &          oilo2:oihi2,
     &          oilo0:oihi0),
     &  arraynew(nilo1:nihi1+1,
     &          nilo2:nihi2,
     &          nilo0:nihi0),
     &  arraydst(dilo1:dihi1+1,
     &          dilo2:dihi2,
     &          dilo0:dihi0)
      integer ie1,ic2,ic0
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic0=ifirst0,ilast0
         do ic2=ifirst2,ilast2
            do ie1=ifirst1,ilast1+1
               arraydst(ie1,ic2,ic0)=
     &                              +arrayold(ie1,ic2,ic0)*oldfrac
     &                              +arraynew(ie1,ic2,ic0)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintfacedoub3d2(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo2:oihi2+1,
     &          oilo0:oihi0,
     &          oilo1:oihi1),
     &  arraynew(nilo2:nihi2+1,
     &          nilo0:nihi0,
     &          nilo1:nihi1),
     &  arraydst(dilo2:dihi2+1,
     &          dilo0:dihi0,
     &          dilo1:dihi1)
      integer ie2,ic0,ic1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic1=ifirst1,ilast1
         do ic0=ifirst0,ilast0
            do ie2=ifirst2,ilast2+1
               arraydst(ie2,ic0,ic1)=
     &                              +arrayold(ie2,ic0,ic1)*oldfrac
     &                              +arraynew(ie2,ic0,ic1)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 3d face-centered float data
c***********************************************************************
c
      subroutine lintimeintfacefloat3d0(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1,
     &          oilo2:oihi2),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1,
     &          nilo2:nihi2),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1,
     &          dilo2:dihi2)
      integer ie0,ic1,ic2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            do ie0=ifirst0,ilast0+1
               arraydst(ie0,ic1,ic2)=
     &                              +arrayold(ie0,ic1,ic2)*oldfrac
     &                              +arraynew(ie0,ic1,ic2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintfacefloat3d1(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo1:oihi1+1,
     &          oilo2:oihi2,
     &          oilo0:oihi0),
     &  arraynew(nilo1:nihi1+1,
     &          nilo2:nihi2,
     &          nilo0:nihi0),
     &  arraydst(dilo1:dihi1+1,
     &          dilo2:dihi2,
     &          dilo0:dihi0)
      integer ie1,ic2,ic0
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic0=ifirst0,ilast0
         do ic2=ifirst2,ilast2
            do ie1=ifirst1,ilast1+1
               arraydst(ie1,ic2,ic0)=
     &                              +arrayold(ie1,ic2,ic0)*oldfrac
     &                              +arraynew(ie1,ic2,ic0)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintfacefloat3d2(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo2:oihi2+1,
     &          oilo0:oihi0,
     &          oilo1:oihi1),
     &  arraynew(nilo2:nihi2+1,
     &          nilo0:nihi0,
     &          nilo1:nihi1),
     &  arraydst(dilo2:dihi2+1,
     &          dilo0:dihi0,
     &          dilo1:dihi1)
      integer ie2,ic0,ic1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic1=ifirst1,ilast1
         do ic0=ifirst0,ilast0
            do ie2=ifirst2,ilast2+1
               arraydst(ie2,ic0,ic1)=
     &                              +arrayold(ie2,ic0,ic1)*oldfrac
     &                              +arraynew(ie2,ic0,ic1)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 3d face-centered complex data
c***********************************************************************
c
      subroutine lintimeintfacecmplx3d0(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1,
     &          oilo2:oihi2),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1,
     &          nilo2:nihi2),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1,
     &          dilo2:dihi2)
      integer ie0,ic1,ic2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            do ie0=ifirst0,ilast0+1
               arraydst(ie0,ic1,ic2)=
     &                              +arrayold(ie0,ic1,ic2)*oldfrac
     &                              +arraynew(ie0,ic1,ic2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintfacecmplx3d1(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo1:oihi1+1,
     &          oilo2:oihi2,
     &          oilo0:oihi0),
     &  arraynew(nilo1:nihi1+1,
     &          nilo2:nihi2,
     &          nilo0:nihi0),
     &  arraydst(dilo1:dihi1+1,
     &          dilo2:dihi2,
     &          dilo0:dihi0)
      integer ie1,ic2,ic0
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic0=ifirst0,ilast0
         do ic2=ifirst2,ilast2
            do ie1=ifirst1,ilast1+1
               arraydst(ie1,ic2,ic0)=
     &                              +arrayold(ie1,ic2,ic0)*oldfrac
     &                              +arraynew(ie1,ic2,ic0)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintfacecmplx3d2(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo2:oihi2+1,
     &          oilo0:oihi0,
     &          oilo1:oihi1),
     &  arraynew(nilo2:nihi2+1,
     &          nilo0:nihi0,
     &          nilo1:nihi1),
     &  arraydst(dilo2:dihi2+1,
     &          dilo0:dihi0,
     &          dilo1:dihi1)
      integer ie2,ic0,ic1
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic1=ifirst1,ilast1
         do ic0=ifirst0,ilast0
            do ie2=ifirst2,ilast2+1
               arraydst(ie2,ic0,ic1)=
     &                              +arrayold(ie2,ic0,ic1)*oldfrac
     &                              +arraynew(ie2,ic0,ic1)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 3d node-centered double data
c***********************************************************************
c
      subroutine lintimeintnodedoub3d(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1+1,
     &          oilo2:oihi2+1),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1+1,
     &          nilo2:nihi2+1),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1+1,
     &          dilo2:dihi2+1)
      integer ie0,ie1,ie2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie2=ifirst2,ilast2+1
         do ie1=ifirst1,ilast1+1
            do ie0=ifirst0,ilast0+1
               arraydst(ie0,ie1,ie2)=
     &                              +arrayold(ie0,ie1,ie2)*oldfrac
     &                              +arraynew(ie0,ie1,ie2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 3d node-centered float data
c***********************************************************************
c
      subroutine lintimeintnodefloat3d(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1+1,
     &          oilo2:oihi2+1),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1+1,
     &          nilo2:nihi2+1),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1+1,
     &          dilo2:dihi2+1)
      integer ie0,ie1,ie2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie2=ifirst2,ilast2+1
         do ie1=ifirst1,ilast1+1
            do ie0=ifirst0,ilast0+1
               arraydst(ie0,ie1,ie2)=
     &                              +arrayold(ie0,ie1,ie2)*oldfrac
     &                              +arraynew(ie0,ie1,ie2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 3d node-centered complex data
c***********************************************************************
c
      subroutine lintimeintnodecmplx3d(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1+1,
     &          oilo2:oihi2+1),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1+1,
     &          nilo2:nihi2+1),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1+1,
     &          dilo2:dihi2+1)
      integer ie0,ie1,ie2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie2=ifirst2,ilast2+1
         do ie1=ifirst1,ilast1+1
            do ie0=ifirst0,ilast0+1
               arraydst(ie0,ie1,ie2)=
     &                              +arrayold(ie0,ie1,ie2)*oldfrac
     &                              +arraynew(ie0,ie1,ie2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 3d outerface double data
c***********************************************************************
c
      subroutine lintimeintoutfacedoub3d0(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo1:oihi1,
     &          oilo2:oihi2),
     &  arraynew(nilo1:nihi1,
     &          nilo2:nihi2),
     &  arraydst(dilo1:dihi1,
     &          dilo2:dihi2)
      integer ic1,ic2 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            arraydst(ic1,ic2)=arrayold(ic1,ic2)*oldfrac
     &                       +arraynew(ic1,ic2)*tfrac
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintoutfacedoub3d1(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo2:oihi2,
     &          oilo0:oihi0),
     &  arraynew(nilo2:nihi2,
     &          nilo0:nihi0),
     &  arraydst(dilo2:dihi2,
     &          dilo0:dihi0)
      integer ic2,ic0 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic0=ifirst0,ilast0
         do ic2=ifirst2,ilast2
            arraydst(ic2,ic0)=arrayold(ic2,ic0)*oldfrac
     &                       +arraynew(ic2,ic0)*tfrac
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintoutfacedoub3d2(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
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
c Linear time interpolation for 3d outerface float data
c***********************************************************************
c
      subroutine lintimeintoutfacefloat3d0(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo1:oihi1,
     &          oilo2:oihi2),
     &  arraynew(nilo1:nihi1,
     &          nilo2:nihi2),
     &  arraydst(dilo1:dihi1,
     &          dilo2:dihi2)
      integer ic1,ic2 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            arraydst(ic1,ic2)=arrayold(ic1,ic2)*oldfrac
     &                       +arraynew(ic1,ic2)*tfrac
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintoutfacefloat3d1(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo2:oihi2,
     &          oilo0:oihi0),
     &  arraynew(nilo2:nihi2,
     &          nilo0:nihi0),
     &  arraydst(dilo2:dihi2,
     &          dilo0:dihi0)
      integer ic2,ic0 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic0=ifirst0,ilast0
         do ic2=ifirst2,ilast2
            arraydst(ic2,ic0)=arrayold(ic2,ic0)*oldfrac
     &                       +arraynew(ic2,ic0)*tfrac
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintoutfacefloat3d2(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
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
c Linear time interpolation for 3d outerface complex data
c***********************************************************************
c
      subroutine lintimeintoutfacecmplx3d0(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo1:oihi1,
     &          oilo2:oihi2),
     &  arraynew(nilo1:nihi1,
     &          nilo2:nihi2),
     &  arraydst(dilo1:dihi1,
     &          dilo2:dihi2)
      integer ic1,ic2 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            arraydst(ic1,ic2)=arrayold(ic1,ic2)*oldfrac
     &                       +arraynew(ic1,ic2)*tfrac
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintoutfacecmplx3d1(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo2:oihi2,
     &          oilo0:oihi0),
     &  arraynew(nilo2:nihi2,
     &          nilo0:nihi0),
     &  arraydst(dilo2:dihi2,
     &          dilo0:dihi0)
      integer ic2,ic0 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic0=ifirst0,ilast0
         do ic2=ifirst2,ilast2
            arraydst(ic2,ic0)=arrayold(ic2,ic0)*oldfrac
     &                       +arraynew(ic2,ic0)*tfrac
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintoutfacecmplx3d2(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
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
c Linear time interpolation for 3d outerside double data
c***********************************************************************
c
      subroutine lintimeintoutsidedoub3d0(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo1:oihi1,
     &          oilo2:oihi2),
     &  arraynew(nilo1:nihi1,
     &          nilo2:nihi2),
     &  arraydst(dilo1:dihi1,
     &          dilo2:dihi2)
      integer ic1,ic2 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            arraydst(ic1,ic2)=arrayold(ic1,ic2)*oldfrac
     &                       +arraynew(ic1,ic2)*tfrac
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintoutsidedoub3d1(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0,
     &          oilo2:oihi2),
     &  arraynew(nilo0:nihi0,
     &          nilo2:nihi2),
     &  arraydst(dilo0:dihi0,
     &          dilo2:dihi2)
      integer ic0,ic2 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic2=ifirst2,ilast2
         do ic0=ifirst0,ilast0
            arraydst(ic0,ic2)=arrayold(ic0,ic2)*oldfrac
     &                       +arraynew(ic0,ic2)*tfrac
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintoutsidedoub3d2(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
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
c Linear time interpolation for 3d outerside float data
c***********************************************************************
c
      subroutine lintimeintoutsidefloat3d0(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo1:oihi1,
     &          oilo2:oihi2),
     &  arraynew(nilo1:nihi1,
     &          nilo2:nihi2),
     &  arraydst(dilo1:dihi1,
     &          dilo2:dihi2)
      integer ic1,ic2 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            arraydst(ic1,ic2)=arrayold(ic1,ic2)*oldfrac
     &                       +arraynew(ic1,ic2)*tfrac
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintoutsidefloat3d1(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0,
     &          oilo2:oihi2),
     &  arraynew(nilo0:nihi0,
     &          nilo2:nihi2),
     &  arraydst(dilo0:dihi0,
     &          dilo2:dihi2)
      integer ic0,ic2 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic2=ifirst2,ilast2
         do ic0=ifirst0,ilast0
            arraydst(ic0,ic2)=arrayold(ic0,ic2)*oldfrac
     &                       +arraynew(ic0,ic2)*tfrac
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintoutsidefloat3d2(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
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
c Linear time interpolation for 3d outerside complex data
c***********************************************************************
c
      subroutine lintimeintoutsidecmplx3d0(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo1:oihi1,
     &          oilo2:oihi2),
     &  arraynew(nilo1:nihi1,
     &          nilo2:nihi2),
     &  arraydst(dilo1:dihi1,
     &          dilo2:dihi2)
      integer ic1,ic2 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            arraydst(ic1,ic2)=arrayold(ic1,ic2)*oldfrac
     &                       +arraynew(ic1,ic2)*tfrac
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintoutsidecmplx3d1(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0,
     &          oilo2:oihi2),
     &  arraynew(nilo0:nihi0,
     &          nilo2:nihi2),
     &  arraydst(dilo0:dihi0,
     &          dilo2:dihi2)
      integer ic0,ic2 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic2=ifirst2,ilast2
         do ic0=ifirst0,ilast0
            arraydst(ic0,ic2)=arrayold(ic0,ic2)*oldfrac
     &                       +arraynew(ic0,ic2)*tfrac
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintoutsidecmplx3d2(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
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
c Linear time interpolation for 3d side-centered double data
c***********************************************************************
c
      subroutine lintimeintsidedoub3d0(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1,
     &          oilo2:oihi2),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1,
     &          nilo2:nihi2),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1,
     &          dilo2:dihi2)
      integer ie0,ie1,ie2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie2=ifirst2,ilast2
         do ie1=ifirst1,ilast1
            do ie0=ifirst0,ilast0+1
               arraydst(ie0,ie1,ie2)=
     &                              +arrayold(ie0,ie1,ie2)*oldfrac
     &                              +arraynew(ie0,ie1,ie2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintsidedoub3d1(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0,
     &          oilo1:oihi1+1,
     &          oilo2:oihi2),
     &  arraynew(nilo0:nihi0,
     &          nilo1:nihi1+1,
     &          nilo2:nihi2),
     &  arraydst(dilo0:dihi0,
     &          dilo1:dihi1+1,
     &          dilo2:dihi2)
      integer ie0,ie1,ie2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie2=ifirst2,ilast2
         do ie1=ifirst1,ilast1+1
            do ie0=ifirst0,ilast0
               arraydst(ie0,ie1,ie2)=
     &                              +arrayold(ie0,ie1,ie2)*oldfrac
     &                              +arraynew(ie0,ie1,ie2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintsidedoub3d2(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double precision
     &  arrayold(oilo0:oihi0,
     &          oilo1:oihi1,
     &          oilo2:oihi2+1),
     &  arraynew(nilo0:nihi0,
     &          nilo1:nihi1,
     &          nilo2:nihi2+1),
     &  arraydst(dilo0:dihi0,
     &          dilo1:dihi1,
     &          dilo2:dihi2+1)
      integer ie0,ie1,ie2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie2=ifirst2,ilast2+1
         do ie1=ifirst1,ilast1
            do ie0=ifirst0,ilast0
               arraydst(ie0,ie1,ie2)=
     &                              +arrayold(ie0,ie1,ie2)*oldfrac
     &                              +arraynew(ie0,ie1,ie2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 3d side-centered float data
c***********************************************************************
c
      subroutine lintimeintsidefloat3d0(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1,
     &          oilo2:oihi2),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1,
     &          nilo2:nihi2),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1,
     &          dilo2:dihi2)
      integer ie0,ie1,ie2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie2=ifirst2,ilast2
         do ie1=ifirst1,ilast1
            do ie0=ifirst0,ilast0+1
               arraydst(ie0,ie1,ie2)=
     &                              +arrayold(ie0,ie1,ie2)*oldfrac
     &                              +arraynew(ie0,ie1,ie2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintsidefloat3d1(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0,
     &          oilo1:oihi1+1,
     &          oilo2:oihi2),
     &  arraynew(nilo0:nihi0,
     &          nilo1:nihi1+1,
     &          nilo2:nihi2),
     &  arraydst(dilo0:dihi0,
     &          dilo1:dihi1+1,
     &          dilo2:dihi2)
      integer ie0,ie1,ie2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie2=ifirst2,ilast2
         do ie1=ifirst1,ilast1+1
            do ie0=ifirst0,ilast0
               arraydst(ie0,ie1,ie2)=
     &                              +arrayold(ie0,ie1,ie2)*oldfrac
     &                              +arraynew(ie0,ie1,ie2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintsidefloat3d2(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      real
     &  arrayold(oilo0:oihi0,
     &          oilo1:oihi1,
     &          oilo2:oihi2+1),
     &  arraynew(nilo0:nihi0,
     &          nilo1:nihi1,
     &          nilo2:nihi2+1),
     &  arraydst(dilo0:dihi0,
     &          dilo1:dihi1,
     &          dilo2:dihi2+1)
      integer ie0,ie1,ie2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie2=ifirst2,ilast2+1
         do ie1=ifirst1,ilast1
            do ie0=ifirst0,ilast0
               arraydst(ie0,ie1,ie2)=
     &                              +arrayold(ie0,ie1,ie2)*oldfrac
     &                              +arraynew(ie0,ie1,ie2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Linear time interpolation for 3d side-centered complex data
c***********************************************************************
c
      subroutine lintimeintsidecmplx3d0(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0+1,
     &          oilo1:oihi1,
     &          oilo2:oihi2),
     &  arraynew(nilo0:nihi0+1,
     &          nilo1:nihi1,
     &          nilo2:nihi2),
     &  arraydst(dilo0:dihi0+1,
     &          dilo1:dihi1,
     &          dilo2:dihi2)
      integer ie0,ie1,ie2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie2=ifirst2,ilast2
         do ie1=ifirst1,ilast1
            do ie0=ifirst0,ilast0+1
               arraydst(ie0,ie1,ie2)=
     &                              +arrayold(ie0,ie1,ie2)*oldfrac
     &                              +arraynew(ie0,ie1,ie2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintsidecmplx3d1(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0,
     &          oilo1:oihi1+1,
     &          oilo2:oihi2),
     &  arraynew(nilo0:nihi0,
     &          nilo1:nihi1+1,
     &          nilo2:nihi2),
     &  arraydst(dilo0:dihi0,
     &          dilo1:dihi1+1,
     &          dilo2:dihi2)
      integer ie0,ie1,ie2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie2=ifirst2,ilast2
         do ie1=ifirst1,ilast1+1
            do ie0=ifirst0,ilast0
               arraydst(ie0,ie1,ie2)=
     &                              +arrayold(ie0,ie1,ie2)*oldfrac
     &                              +arraynew(ie0,ie1,ie2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine lintimeintsidecmplx3d2(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
      double complex
     &  arrayold(oilo0:oihi0,
     &          oilo1:oihi1,
     &          oilo2:oihi2+1),
     &  arraynew(nilo0:nihi0,
     &          nilo1:nihi1,
     &          nilo2:nihi2+1),
     &  arraydst(dilo0:dihi0,
     &          dilo1:dihi1,
     &          dilo2:dihi2+1)
      integer ie0,ie1,ie2
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ie2=ifirst2,ilast2+1
         do ie1=ifirst1,ilast1
            do ie0=ifirst0,ilast0
               arraydst(ie0,ie1,ie2)=
     &                              +arrayold(ie0,ie1,ie2)*oldfrac
     &                              +arraynew(ie0,ie1,ie2)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
c
