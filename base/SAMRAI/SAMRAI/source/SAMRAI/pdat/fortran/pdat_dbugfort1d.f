c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: F77 routines for debugging 1d patch data types.
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for 1d patchdata debugging routines
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
c***********************************************************************
c Debugging routines for 1d cell-centered data
c***********************************************************************
c
      subroutine dbugcelldoub1d(
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
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
      subroutine dbugcellflot1d(
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
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
      subroutine dbugcellcplx1d(
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
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
c Debugging routines for 1d face-centered data
c***********************************************************************
c
      subroutine dbugfacedoub1d(
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
      double precision
     &  array(fi0-ng:la0+1+ng)
      integer ie0
c     =============================================================

      do ie0=ibeg0,iend0+1
         write(6,*) "array[",ie0,"] = ",array(ie0)
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugfaceflot1d(
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
      real
     &  array(fi0-ng:la0+1+ng)
      integer ie0
c     =============================================================

      do ie0=ibeg0,iend0+1
         write(6,*) "array[",ie0,"] = ",array(ie0)
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugfacecplx1d(
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
      double complex
     &  array(fi0-ng:la0+1+ng)
      integer ie0
c     =============================================================

      do ie0=ibeg0,iend0+1
         write(6,*) "array[",ie0,"] = ",array(ie0)
      enddo

      call flush(6)
      return
      end
c
c
c***********************************************************************
c Debugging routines for 1d node-centered data
c***********************************************************************
c
      subroutine dbugnodedoub1d(
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
      double precision
     &  array(fi0-ng:la0+1+ng)
      integer ie0
c     =============================================================

      do ie0=ibeg0,iend0+1
         write(6,*) "array[",ie0,"] = ",array(ie0)
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugnodeflot1d(
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
      real
     &  array(fi0-ng:la0+1+ng)
      integer ie0
c     =============================================================

      do ie0=ibeg0,iend0+1
         write(6,*) "array[",ie0,"] = ",array(ie0)
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugnodecplx1d(
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
      double complex
     &  array(fi0-ng:la0+1+ng)
      integer ie0
c     =============================================================

      do ie0=ibeg0,iend0+1
         write(6,*) "array[",ie0,"] = ",array(ie0)
      enddo

      call flush(6)
      return
      end
c
c
c***********************************************************************
c Debugging routines for 1d outerface data
c***********************************************************************
c
      subroutine dbugoutfacedoub1d(
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
      double precision
     &  array(1)
c     =============================================================

      write(6,*) "array[",1,"] = ",array(1)

      call flush(6)
      return
      end
c
      subroutine dbugoutfaceflot1d(
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
      real
     &  array(1)
c     =============================================================

      write(6,*) "array[",1,"] = ",array(1)

      call flush(6)
      return
      end
c
      subroutine dbugoutfacecplx1d(
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
      double complex
     &  array(1)
c     =============================================================

      write(6,*) "array[",1,"] = ",array(1)

      call flush(6)
      return
      end
c
c
c***********************************************************************
c Debugging routines for 1d outerside data
c***********************************************************************
c
      subroutine dbugoutsidedoub1d(
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
      double precision
     &  array(1)
c     =============================================================

      write(6,*) "array[",1,"] = ",array(1)

      call flush(6)
      return
      end
c
      subroutine dbugoutsideflot1d(
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
      real
     &  array(1)
c     =============================================================

      write(6,*) "array[",1,"] = ",array(1)

      call flush(6)
      return
      end
c
      subroutine dbugoutsidecplx1d(
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
      double complex
     &  array(1)
c     =============================================================

      write(6,*) "array[",1,"] = ",array(1)

      call flush(6)
      return
      end
c
c
c***********************************************************************
c Debugging routines for 1d side-centered data
c***********************************************************************
c
      subroutine dbugsidedoub1d(
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
      double precision
     &  array(fi0-ng:la0+1+ng)
      integer ie0
c     =============================================================

      do ie0=ibeg0,iend0+1
         write(6,*) "array[",ie0,"] = ",array(ie0)
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugsideflot1d(
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
      real
     &  array(fi0-ng:la0+1+ng)
      integer ie0
c     =============================================================

      do ie0=ibeg0,iend0+1
         write(6,*) "array[",ie0,"] = ",array(ie0)
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugsidecplx1d(
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
      double complex
     &  array(fi0-ng:la0+1+ng)
      integer ie0
c     =============================================================

      do ie0=ibeg0,iend0+1
         write(6,*) "array[",ie0,"] = ",array(ie0)
      enddo

      call flush(6)
      return
      end
c
c
c***********************************************************************
c Debugging routines for 1d edge-centered data
c***********************************************************************
c
      subroutine dbugedgedoub1d(
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
      double precision
     &  array(fi0-ng:la0+ng)
      integer ie0
c     =============================================================

      do ie0=ibeg0,iend0
         write(6,*) "array[",ie0,"] = ",array(ie0)
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugedgeflot1d(
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
      real
     &  array(fi0-ng:la0+ng)
      integer ie0
c     =============================================================

      do ie0=ibeg0,iend0
         write(6,*) "array[",ie0,"] = ",array(ie0)
      enddo

      call flush(6)
      return
      end
c
      subroutine dbugedgecplx1d(
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
      double complex
     &  array(fi0-ng:la0+ng)
      integer ie0
c     =============================================================

      do ie0=ibeg0,iend0
         write(6,*) "array[",ie0,"] = ",array(ie0)
      enddo

      call flush(6)
      return
      end
c
