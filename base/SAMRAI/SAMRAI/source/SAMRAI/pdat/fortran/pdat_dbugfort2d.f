c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: F77 routines for debugging 2d patch data types.
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for 2d patchdata debugging routines
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
c***********************************************************************
c Debugging routines for 2d cell-centered data
c***********************************************************************
c
      subroutine dbugcelldoub2d(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double precision
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+ng)
      integer ic0,ic1
c     =============================================================

      do ic1=ibeg1,iend1
         do ic0=ibeg0,iend0 
            write(6,*) "array[",ic0,",",ic1,"] = ",array(ic0,ic1)
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugcellflot2d(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      real
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+ng)
      integer ic0,ic1
c     =============================================================

      do ic1=ibeg1,iend1
         do ic0=ibeg0,iend0 
            write(6,*) "array[",ic0,",",ic1,"] = ",array(ic0,ic1)
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugcellcplx2d(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double complex
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+ng)
      integer ic0,ic1
c     =============================================================

      do ic1=ibeg1,iend1
         do ic0=ibeg0,iend0 
            write(6,*) "array[",ic0,",",ic1,"] = ",array(ic0,ic1)
         enddo
      enddo

      call flush(6)
      return
      end
c
c
c***********************************************************************
c Debugging routines for 2d face-centered data
c***********************************************************************
c
      subroutine dbugfacedoub2d0(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double precision
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+ng)
      integer ie0,ic1
c     =============================================================

      do ic1=ibeg1,iend1
         do ie0=ibeg0,iend0+1 
            write(6,*) "array[",ie0,",",ic1,"] = ",array(ie0,ic1)
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugfacedoub2d1(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double precision
     &  array(fi1-ng:la1+1+ng,
     &          fi0-ng:la0+ng)
      integer ie1,ic0
c     =============================================================

      do ic0=ibeg0,iend0
         do ie1=ibeg1,iend1+1 
            write(6,*) "array[",ie1,",",ic0,"] = ",array(ie1,ic0)
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugfaceflot2d0(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      real
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+ng)
      integer ie0,ic1
c     =============================================================

      do ic1=ibeg1,iend1
         do ie0=ibeg0,iend0+1 
            write(6,*) "array[",ie0,",",ic1,"] = ",array(ie0,ic1)
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugfaceflot2d1(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      real
     &  array(fi1-ng:la1+1+ng,
     &          fi0-ng:la0+ng)
      integer ie1,ic0
c     =============================================================

      do ic0=ibeg0,iend0
         do ie1=ibeg1,iend1+1 
            write(6,*) "array[",ie1,",",ic0,"] = ",array(ie1,ic0)
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugfacecplx2d0(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double complex
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+ng)
      integer ie0,ic1
c     =============================================================

      do ic1=ibeg1,iend1
         do ie0=ibeg0,iend0+1 
            write(6,*) "array[",ie0,",",ic1,"] = ",array(ie0,ic1)
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugfacecplx2d1(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double complex
     &  array(fi1-ng:la1+1+ng,
     &          fi0-ng:la0+ng)
      integer ie1,ic0
c     =============================================================

      do ic0=ibeg0,iend0
         do ie1=ibeg1,iend1+1 
            write(6,*) "array[",ie1,",",ic0,"] = ",array(ie1,ic0)
         enddo
      enddo

      call flush(6)
      return
      end
c
c
c***********************************************************************
c Debugging routines for 2d node-centered data
c***********************************************************************
c
      subroutine dbugnodedoub2d(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double precision
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+1+ng)
      integer ie0,ie1
c     =============================================================

      do ie1=ibeg1,iend1+1
         do ie0=ibeg0,iend0+1 
            write(6,*) "array[",ie0,",",ie1,"] = ",array(ie0,ie1)
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugnodeflot2d(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      real
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+1+ng)
      integer ie0,ie1
c     =============================================================

      do ie1=ibeg1,iend1+1
         do ie0=ibeg0,iend0+1 
            write(6,*) "array[",ie0,",",ie1,"] = ",array(ie0,ie1)
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugnodecplx2d(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double complex
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+1+ng)
      integer ie0,ie1
c     =============================================================

      do ie1=ibeg1,iend1+1
         do ie0=ibeg0,iend0+1 
            write(6,*) "array[",ie0,",",ie1,"] = ",array(ie0,ie1)
         enddo
      enddo

      call flush(6)
      return
      end
c
c
c***********************************************************************
c Debugging routines for 2d outerface data
c***********************************************************************
c
      subroutine dbugoutfacedoub2d0(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double precision
     &  array(fi1-ng:la1+ng)
      integer ic1
c     =============================================================

      do ic1=ibeg1,iend1
         write(6,*) "array[",ic1,"] = ",array(ic1)
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutfacedoub2d1(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double precision
     &  array(fi0-ng:la0+ng)
      integer ic0
c     =============================================================

      do ic0=ibeg0,iend0
         write(6,*) "array[",ic0,"] = ",array(ic0)
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutfaceflot2d0(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      real
     &  array(fi1-ng:la1+ng)
      integer ic1
c     =============================================================

      do ic1=ibeg1,iend1
         write(6,*) "array[",ic1,"] = ",array(ic1)
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutfaceflot2d1(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      real
     &  array(fi0-ng:la0+ng)
      integer ic0
c     =============================================================

      do ic0=ibeg0,iend0
         write(6,*) "array[",ic0,"] = ",array(ic0)
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutfacecplx2d0(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double complex
     &  array(fi1-ng:la1+ng)
      integer ic1
c     =============================================================

      do ic1=ibeg1,iend1
         write(6,*) "array[",ic1,"] = ",array(ic1)
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutfacecplx2d1(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double complex
     &  array(fi0-ng:la0+ng)
      integer ic0
c     =============================================================

      do ic0=ibeg0,iend0
         write(6,*) "array[",ic0,"] = ",array(ic0)
      enddo
 
      call flush(6)
      return
      end
c
c
c***********************************************************************
c Debugging routines for 2d outerside data
c***********************************************************************
c
      subroutine dbugoutsidedoub2d0(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double precision
     &  array(fi1-ng:la1+ng)
      integer ic1
c     =============================================================

      do ic1=ibeg1,iend1
         write(6,*) "array[",ic1,"] = ",array(ic1)
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutsidedoub2d1(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double precision
     &  array(fi0-ng:la0+ng)
      integer ic0
c     =============================================================

      do ic0=ibeg0,iend0
         write(6,*) "array[",ic0,"] = ",array(ic0)
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutsideflot2d0(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      real
     &  array(fi1-ng:la1+ng)
      integer ic1
c     =============================================================

      do ic1=ibeg1,iend1
         write(6,*) "array[",ic1,"] = ",array(ic1)
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutsideflot2d1(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      real
     &  array(fi0-ng:la0+ng)
      integer ic0
c     =============================================================

      do ic0=ibeg0,iend0
         write(6,*) "array[",ic0,"] = ",array(ic0)
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutsidecplx2d0(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double complex
     &  array(fi1-ng:la1+ng)
      integer ic1
c     =============================================================

      do ic1=ibeg1,iend1
         write(6,*) "array[",ic1,"] = ",array(ic1)
      enddo
 
      call flush(6)
      return
      end
c
      subroutine dbugoutsidecplx2d1(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double complex
     &  array(fi0-ng:la0+ng)
      integer ic0
c     =============================================================

      do ic0=ibeg0,iend0
         write(6,*) "array[",ic0,"] = ",array(ic0)
      enddo
 
      call flush(6)
      return
      end
c
c
c***********************************************************************
c Debugging routines for 2d side-centered data
c***********************************************************************
c
      subroutine dbugsidedoub2d0(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double precision
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+ng)
      integer ie0,ic1
c     =============================================================

      do ic1=ibeg1,iend1+0
         do ie0=ibeg0,iend0+1 
            write(6,*) "array[",ie0,",",ic1,"] = ",array(ie0,ic1)
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugsidedoub2d1(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double precision
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+1+ng)
      integer ie0,ic1
c     =============================================================

      do ic1=ibeg1,iend1+1
         do ie0=ibeg0,iend0+0 
            write(6,*) "array[",ie0,",",ic1,"] = ",array(ie0,ic1)
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugsideflot2d0(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      real
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+ng)
      integer ie0,ic1
c     =============================================================

      do ic1=ibeg1,iend1+0
         do ie0=ibeg0,iend0+1 
            write(6,*) "array[",ie0,",",ic1,"] = ",array(ie0,ic1)
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugsideflot2d1(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      real
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+1+ng)
      integer ie0,ic1
c     =============================================================

      do ic1=ibeg1,iend1+1
         do ie0=ibeg0,iend0+0 
            write(6,*) "array[",ie0,",",ic1,"] = ",array(ie0,ic1)
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugsidecplx2d0(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double complex
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+ng)
      integer ie0,ic1
c     =============================================================

      do ic1=ibeg1,iend1+0
         do ie0=ibeg0,iend0+1 
            write(6,*) "array[",ie0,",",ic1,"] = ",array(ie0,ic1)
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugsidecplx2d1(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double complex
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+1+ng)
      integer ie0,ic1
c     =============================================================

      do ic1=ibeg1,iend1+1
         do ie0=ibeg0,iend0+0 
            write(6,*) "array[",ie0,",",ic1,"] = ",array(ie0,ic1)
         enddo
      enddo

      call flush(6)
      return
      end
c
c
c***********************************************************************
c Debugging routines for 2d edge-centered data
c***********************************************************************
c
      subroutine dbugedgedoub2d0(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double precision
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+1+ng)
      integer ie0,ic1
c     =============================================================

      do ic1=ibeg1,iend1+1
         do ie0=ibeg0,iend0+0 
            write(6,*) "array[",ie0,",",ic1,"] = ",array(ie0,ic1)
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugedgedoub2d1(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double precision
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+ng)
      integer ie0,ic1
c     =============================================================

      do ic1=ibeg1,iend1+0
         do ie0=ibeg0,iend0+1 
            write(6,*) "array[",ie0,",",ic1,"] = ",array(ie0,ic1)
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugedgeflot2d0(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      real
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+1+ng)
      integer ie0,ic1
c     =============================================================

      do ic1=ibeg1,iend1+1
         do ie0=ibeg0,iend0+0 
            write(6,*) "array[",ie0,",",ic1,"] = ",array(ie0,ic1)
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugedgeflot2d1(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      real
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+ng)
      integer ie0,ic1
c     =============================================================

      do ic1=ibeg1,iend1+0
         do ie0=ibeg0,iend0+1 
            write(6,*) "array[",ie0,",",ic1,"] = ",array(ie0,ic1)
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugedgecplx2d0(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double complex
     &  array(fi0-ng:la0+ng,
     &          fi1-ng:la1+1+ng)
      integer ie0,ic1
c     =============================================================

      do ic1=ibeg1,iend1+1
         do ie0=ibeg0,iend0+0 
            write(6,*) "array[",ie0,",",ic1,"] = ",array(ie0,ic1)
         enddo
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugedgecplx2d1(
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
      double complex
     &  array(fi0-ng:la0+1+ng,
     &          fi1-ng:la1+ng)
      integer ie0,ic1
c     =============================================================

      do ic1=ibeg1,iend1+0
         do ie0=ibeg0,iend0+1 
            write(6,*) "array[",ie0,",",ic1,"] = ",array(ie0,ic1)
         enddo
      enddo

      call flush(6)
      return
      end
c
