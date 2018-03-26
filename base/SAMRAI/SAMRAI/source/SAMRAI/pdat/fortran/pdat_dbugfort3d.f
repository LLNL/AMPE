c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: F77 routines for debugging 3d patch data types.
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for 3d patchdata debugging routines
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
c***********************************************************************
c Debugging routines for 3d cell-centered data
c***********************************************************************
c
      subroutine dbugcelldoub3d(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double precision
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+ng,
     &          fi2-ng:la2+ng)
      integer ic0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2
         do ic1=ibeg1,iend1
            do ic0=ibeg0,iend0 
               write(6,*) "array[",ic0,",",ic1,",",ic2,"] = ",
     &         array(ic0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugcellflot3d(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      real
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+ng,
     &          fi2-ng:la2+ng)
      integer ic0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2
         do ic1=ibeg1,iend1
            do ic0=ibeg0,iend0 
               write(6,*) "array[",ic0,",",ic1,",",ic2,"] = ",
     &         array(ic0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugcellcplx3d(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double complex
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+ng,
     &          fi2-ng:la2+ng)
      integer ic0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2
         do ic1=ibeg1,iend1
            do ic0=ibeg0,iend0 
               write(6,*) "array[",ic0,",",ic1,",",ic2,"] = ",
     &         array(ic0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
c
c***********************************************************************
c Debugging routines for 3d face-centered data
c***********************************************************************
c
      subroutine dbugfacedoub3d0(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double precision
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+ng,
     &          fi2-ng:la2+ng)
      integer ie0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2
         do ic1=ibeg1,iend1
            do ie0=ibeg0,iend0+1 
               write(6,*) "array[",ie0,",",ic1,",",ic2,"] = ",
     &         array(ie0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugfacedoub3d1(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double precision
     &  array(fi1-ng:la1+1+ng,
     &          fi2-ng:la2+ng,
     &          fi0-ng:la0+ng)
      integer ie1,ic2,ic0
c     =============================================================

      do ic0=ibeg0,iend0
         do ic2=ibeg2,iend2
            do ie1=ibeg1,iend1+1 
               write(6,*) "array[",ie1,",",ic2,",",ic0,"] = ",
     &         array(ie1,ic2,ic0)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugfacedoub3d2(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double precision
     &  array(fi2-ng:la2+1+ng,
     &          fi0-ng:la0+ng,
     &          fi1-ng:la1+ng)
      integer ie2,ic0,ic1
c     =============================================================

      do ic1=ibeg1,iend1
         do ic0=ibeg0,iend0
            do ie2=ibeg2,iend2+1 
               write(6,*) "array[",ie2,",",ic0,",",ic1,"] = ",
     &         array(ie2,ic0,ic1)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugfaceflot3d0(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      real
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+ng,
     &          fi2-ng:la2+ng)
      integer ie0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2
         do ic1=ibeg1,iend1
            do ie0=ibeg0,iend0+1 
               write(6,*) "array[",ie0,",",ic1,",",ic2,"] = ",
     &         array(ie0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugfaceflot3d1(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      real
     &  array(fi1-ng:la1+1+ng,
     &          fi2-ng:la2+ng,
     &          fi0-ng:la0+ng)
      integer ie1,ic2,ic0
c     =============================================================

      do ic0=ibeg0,iend0
         do ic2=ibeg2,iend2
            do ie1=ibeg1,iend1+1 
               write(6,*) "array[",ie1,",",ic2,",",ic0,"] = ",
     &         array(ie1,ic2,ic0)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugfaceflot3d2(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      real
     &  array(fi2-ng:la2+1+ng,
     &          fi0-ng:la0+ng,
     &          fi1-ng:la1+ng)
      integer ie2,ic0,ic1
c     =============================================================

      do ic1=ibeg1,iend1
         do ic0=ibeg0,iend0
            do ie2=ibeg2,iend2+1 
               write(6,*) "array[",ie2,",",ic0,",",ic1,"] = ",
     &         array(ie2,ic0,ic1)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugfacecplx3d0(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double complex
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+ng,
     &          fi2-ng:la2+ng)
      integer ie0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2
         do ic1=ibeg1,iend1
            do ie0=ibeg0,iend0+1 
               write(6,*) "array[",ie0,",",ic1,",",ic2,"] = ",
     &         array(ie0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugfacecplx3d1(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double complex
     &  array(fi1-ng:la1+1+ng,
     &          fi2-ng:la2+ng,
     &          fi0-ng:la0+ng)
      integer ie1,ic2,ic0
c     =============================================================

      do ic0=ibeg0,iend0
         do ic2=ibeg2,iend2
            do ie1=ibeg1,iend1+1 
               write(6,*) "array[",ie1,",",ic2,",",ic0,"] = ",
     &         array(ie1,ic2,ic0)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugfacecplx3d2(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double complex
     &  array(fi2-ng:la2+1+ng,
     &          fi0-ng:la0+ng,
     &          fi1-ng:la1+ng)
      integer ie2,ic0,ic1
c     =============================================================

      do ic1=ibeg1,iend1
         do ic0=ibeg0,iend0
            do ie2=ibeg2,iend2+1 
               write(6,*) "array[",ie2,",",ic0,",",ic1,"] = ",
     &         array(ie2,ic0,ic1)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
c
c***********************************************************************
c Debugging routines for 3d node-centered data
c***********************************************************************
c
      subroutine dbugnodedoub3d(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double precision
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+1+ng,
     &          fi2-ng:la2+1+ng)
      integer ie0,ie1,ie2
c     =============================================================

      do ie2=ibeg2,iend2+1
         do ie1=ibeg1,iend1+1
            do ie0=ibeg0,iend0+1 
               write(6,*) "array[",ie0,",",ie1,",",ie2,"] = ",
     &         array(ie0,ie1,ie2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugnodeflot3d(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      real
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+1+ng,
     &          fi2-ng:la2+1+ng)
      integer ie0,ie1,ie2
c     =============================================================

      do ie2=ibeg2,iend2+1
         do ie1=ibeg1,iend1+1
            do ie0=ibeg0,iend0+1 
               write(6,*) "array[",ie0,",",ie1,",",ie2,"] = ",
     &         array(ie0,ie1,ie2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugnodecplx3d(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double complex
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+1+ng,
     &          fi2-ng:la2+1+ng)
      integer ie0,ie1,ie2
c     =============================================================

      do ie2=ibeg2,iend2+1
         do ie1=ibeg1,iend1+1
            do ie0=ibeg0,iend0+1 
               write(6,*) "array[",ie0,",",ie1,",",ie2,"] = ",
     &         array(ie0,ie1,ie2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
c
c***********************************************************************
c Debugging routines for 3d outerface data
c***********************************************************************
c
      subroutine dbugoutfacedoub3d0(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double precision
     &  array(fi1-ng:la1+ng,
     &          fi2-ng:la2+ng)
      integer ic1,ic2  
c     =============================================================

      do ic2=ibeg2,iend2
         do ic1=ibeg1,iend1
            write(6,*) "array[",ic1,",",ic2,"] = ",
     &      array(ic1,ic2)
         enddo
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutfacedoub3d1(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double precision
     &  array(fi2-ng:la2+ng,
     &          fi0-ng:la0+ng)
      integer ic2,ic0  
c     =============================================================

      do ic0=ibeg0,iend0
         do ic2=ibeg2,iend2
            write(6,*) "array[",ic2,",",ic0,"] = ",
     &      array(ic2,ic0)
         enddo
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutfacedoub3d2(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double precision
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+ng)
      integer ic0,ic1  
c     =============================================================

      do ic1=ibeg1,iend1
         do ic0=ibeg0,iend0
            write(6,*) "array[",ic0,",",ic1,"] = ",
     &      array(ic0,ic1)
         enddo
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutfaceflot3d0(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      real
     &  array(fi1-ng:la1+ng,
     &          fi2-ng:la2+ng)
      integer ic1,ic2  
c     =============================================================

      do ic2=ibeg2,iend2
         do ic1=ibeg1,iend1
            write(6,*) "array[",ic1,",",ic2,"] = ",
     &      array(ic1,ic2)
         enddo
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutfaceflot3d1(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      real
     &  array(fi2-ng:la2+ng,
     &          fi0-ng:la0+ng)
      integer ic2,ic0  
c     =============================================================

      do ic0=ibeg0,iend0
         do ic2=ibeg2,iend2
            write(6,*) "array[",ic2,",",ic0,"] = ",
     &      array(ic2,ic0)
         enddo
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutfaceflot3d2(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      real
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+ng)
      integer ic0,ic1  
c     =============================================================

      do ic1=ibeg1,iend1
         do ic0=ibeg0,iend0
            write(6,*) "array[",ic0,",",ic1,"] = ",
     &      array(ic0,ic1)
         enddo
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutfacecplx3d0(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double complex
     &  array(fi1-ng:la1+ng,
     &          fi2-ng:la2+ng)
      integer ic1,ic2  
c     =============================================================

      do ic2=ibeg2,iend2
         do ic1=ibeg1,iend1
            write(6,*) "array[",ic1,",",ic2,"] = ",
     &      array(ic1,ic2)
         enddo
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutfacecplx3d1(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double complex
     &  array(fi2-ng:la2+ng,
     &          fi0-ng:la0+ng)
      integer ic2,ic0  
c     =============================================================

      do ic0=ibeg0,iend0
         do ic2=ibeg2,iend2
            write(6,*) "array[",ic2,",",ic0,"] = ",
     &      array(ic2,ic0)
         enddo
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutfacecplx3d2(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double complex
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+ng)
      integer ic0,ic1  
c     =============================================================

      do ic1=ibeg1,iend1
         do ic0=ibeg0,iend0
            write(6,*) "array[",ic0,",",ic1,"] = ",
     &      array(ic0,ic1)
         enddo
      enddo
 
      call flush(6)
      return
      end
c
c
c***********************************************************************
c Debugging routines for 3d outerside data
c***********************************************************************
c
      subroutine dbugoutsidedoub3d0(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double precision
     &  array(fi1-ng:la1+ng,
     &          fi2-ng:la2+ng)
      integer ic1,ic2  
c     =============================================================

      do ic2=ibeg2,iend2
         do ic1=ibeg1,iend1
            write(6,*) "array[",ic1,",",ic2,"] = ",
     &      array(ic1,ic2)
         enddo
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutsidedoub3d1(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double precision
     &  array(fi0-ng:la0+ng,
     &          fi2-ng:la2+ng)
      integer ic2,ic0  
c     =============================================================

      do ic0=ibeg0,iend0
         do ic2=ibeg2,iend2
            write(6,*) "array[",ic2,",",ic0,"] = ",
     &      array(ic2,ic0)
         enddo
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutsidedoub3d2(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double precision
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+ng)
      integer ic0,ic1  
c     =============================================================

      do ic1=ibeg1,iend1
         do ic0=ibeg0,iend0
            write(6,*) "array[",ic0,",",ic1,"] = ",
     &      array(ic0,ic1)
         enddo
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutsideflot3d0(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      real
     &  array(fi1-ng:la1+ng,
     &          fi2-ng:la2+ng)
      integer ic1,ic2  
c     =============================================================

      do ic2=ibeg2,iend2
         do ic1=ibeg1,iend1
            write(6,*) "array[",ic1,",",ic2,"] = ",
     &      array(ic1,ic2)
         enddo
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutsideflot3d1(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      real
     &  array(fi0-ng:la0+ng,
     &          fi2-ng:la2+ng)
      integer ic2,ic0  
c     =============================================================

      do ic0=ibeg0,iend0
         do ic2=ibeg2,iend2
            write(6,*) "array[",ic2,",",ic0,"] = ",
     &      array(ic2,ic0)
         enddo
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutsideflot3d2(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      real
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+ng)
      integer ic0,ic1  
c     =============================================================

      do ic1=ibeg1,iend1
         do ic0=ibeg0,iend0
            write(6,*) "array[",ic0,",",ic1,"] = ",
     &      array(ic0,ic1)
         enddo
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutsidecplx3d0(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double complex
     &  array(fi1-ng:la1+ng,
     &          fi2-ng:la2+ng)
      integer ic1,ic2  
c     =============================================================

      do ic2=ibeg2,iend2
         do ic1=ibeg1,iend1
            write(6,*) "array[",ic1,",",ic2,"] = ",
     &      array(ic1,ic2)
         enddo
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutsidecplx3d1(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double complex
     &  array(fi0-ng:la0+ng,
     &          fi2-ng:la2+ng)
      integer ic2,ic0  
c     =============================================================

      do ic0=ibeg0,iend0
         do ic2=ibeg2,iend2
            write(6,*) "array[",ic2,",",ic0,"] = ",
     &      array(ic2,ic0)
         enddo
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutsidecplx3d2(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double complex
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+ng)
      integer ic0,ic1  
c     =============================================================

      do ic1=ibeg1,iend1
         do ic0=ibeg0,iend0
            write(6,*) "array[",ic0,",",ic1,"] = ",
     &      array(ic0,ic1)
         enddo
      enddo
 
      call flush(6)
      return
      end
c
c
c***********************************************************************
c Debugging routines for 3d side-centered data
c***********************************************************************
c
      subroutine dbugsidedoub3d0(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double precision
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+ng,
     &          fi2-ng:la2+ng)
      integer ie0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2
         do ic1=ibeg1,iend1
            do ie0=ibeg0,iend0+1 
               write(6,*) "array[",ie0,",",ic1,",",ic2,"] = ",
     &         array(ie0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugsidedoub3d1(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double precision
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+1+ng,
     &          fi2-ng:la2+ng)
      integer ie0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2
         do ic1=ibeg1,iend1+1
            do ie0=ibeg0,iend0 
               write(6,*) "array[",ie0,",",ic1,",",ic2,"] = ",
     &         array(ie0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugsidedoub3d2(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double precision
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+ng,
     &          fi2-ng:la2+1+ng)
      integer ie0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2+1
         do ic1=ibeg1,iend1
            do ie0=ibeg0,iend0 
               write(6,*) "array[",ie0,",",ic1,",",ic2,"] = ",
     &         array(ie0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugsideflot3d0(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      real
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+ng,
     &          fi2-ng:la2+ng)
      integer ie0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2
         do ic1=ibeg1,iend1
            do ie0=ibeg0,iend0+1 
               write(6,*) "array[",ie0,",",ic1,",",ic2,"] = ",
     &         array(ie0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugsideflot3d1(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      real
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+1+ng,
     &          fi2-ng:la2+ng)
      integer ie0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2
         do ic1=ibeg1,iend1+1
            do ie0=ibeg0,iend0 
               write(6,*) "array[",ie0,",",ic1,",",ic2,"] = ",
     &         array(ie0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugsideflot3d2(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      real
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+ng,
     &          fi2-ng:la2+1+ng)
      integer ie0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2+1
         do ic1=ibeg1,iend1
            do ie0=ibeg0,iend0 
               write(6,*) "array[",ie0,",",ic1,",",ic2,"] = ",
     &         array(ie0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugsidecplx3d0(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double complex
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+ng,
     &          fi2-ng:la2+ng)
      integer ie0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2
         do ic1=ibeg1,iend1
            do ie0=ibeg0,iend0+1 
               write(6,*) "array[",ie0,",",ic1,",",ic2,"] = ",
     &         array(ie0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugsidecplx3d1(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double complex
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+1+ng,
     &          fi2-ng:la2+ng)
      integer ie0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2
         do ic1=ibeg1,iend1+1
            do ie0=ibeg0,iend0 
               write(6,*) "array[",ie0,",",ic1,",",ic2,"] = ",
     &         array(ie0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugsidecplx3d2(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double complex
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+ng,
     &          fi2-ng:la2+1+ng)
      integer ie0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2+1
         do ic1=ibeg1,iend1
            do ie0=ibeg0,iend0 
               write(6,*) "array[",ie0,",",ic1,",",ic2,"] = ",
     &         array(ie0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
c
c***********************************************************************
c Debugging routines for 3d edge-centered data
c***********************************************************************
c
      subroutine dbugedgedoub3d0(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double precision
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+1+ng,
     &          fi2-ng:la2+1+ng)
      integer ie0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2+1
         do ic1=ibeg1,iend1+1
            do ie0=ibeg0,iend0 
               write(6,*) "array[",ie0,",",ic1,",",ic2,"] = ",
     &         array(ie0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugedgedoub3d1(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double precision
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+ng,
     &          fi2-ng:la2+1+ng)
      integer ie0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2+1
         do ic1=ibeg1,iend1
            do ie0=ibeg0,iend0+1 
               write(6,*) "array[",ie0,",",ic1,",",ic2,"] = ",
     &         array(ie0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugedgedoub3d2(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double precision
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+1+ng,
     &          fi2-ng:la2+ng)
      integer ie0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2
         do ic1=ibeg1,iend1+1
            do ie0=ibeg0,iend0+1 
               write(6,*) "array[",ie0,",",ic1,",",ic2,"] = ",
     &         array(ie0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugedgeflot3d0(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      real
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+1+ng,
     &          fi2-ng:la2+1+ng)
      integer ie0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2+1
         do ic1=ibeg1,iend1+1
            do ie0=ibeg0,iend0 
               write(6,*) "array[",ie0,",",ic1,",",ic2,"] = ",
     &         array(ie0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugedgeflot3d1(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      real
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+ng,
     &          fi2-ng:la2+1+ng)
      integer ie0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2+1
         do ic1=ibeg1,iend1
            do ie0=ibeg0,iend0+1 
               write(6,*) "array[",ie0,",",ic1,",",ic2,"] = ",
     &         array(ie0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugedgeflot3d2(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      real
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+1+ng,
     &          fi2-ng:la2+ng)
      integer ie0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2
         do ic1=ibeg1,iend1+1
            do ie0=ibeg0,iend0+1 
               write(6,*) "array[",ie0,",",ic1,",",ic2,"] = ",
     &         array(ie0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugedgecplx3d0(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double complex
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+1+ng,
     &          fi2-ng:la2+1+ng)
      integer ie0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2+1
         do ic1=ibeg1,iend1+1
            do ie0=ibeg0,iend0 
               write(6,*) "array[",ie0,",",ic1,",",ic2,"] = ",
     &         array(ie0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugedgecplx3d1(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double complex
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+ng,
     &          fi2-ng:la2+1+ng)
      integer ie0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2+1
         do ic1=ibeg1,iend1
            do ie0=ibeg0,iend0+1 
               write(6,*) "array[",ie0,",",ic1,",",ic2,"] = ",
     &         array(ie0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugedgecplx3d2(
     &  fi0,la0,fi1,la1,fi2,la2,ng,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,fi2,la2,ng,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      double complex
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+1+ng,
     &          fi2-ng:la2+ng)
      integer ie0,ic1,ic2
c     =============================================================

      do ic2=ibeg2,iend2
         do ic1=ibeg1,iend1+1
            do ie0=ibeg0,iend0+1 
               write(6,*) "array[",ie0,",",ic1,",",ic2,"] = ",
     &         array(ie0,ic1,ic2)
            enddo
         enddo
      enddo

      call flush(6)
      return
      end
c
