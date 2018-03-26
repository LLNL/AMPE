c
c  File:        $URL$
c  Package:     SAMRAI geometry
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: FORTRAN routines for spatial coarsening of 1d patch data
c               on a regular Cartesian mesh.
c
c
c  File:        $URL$
c  Package:     SAMRAI geometry
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for 1d Cartesian coarsen operators
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
c***********************************************************************
c Weighted averaging for 1d cell-centered double data
c***********************************************************************
c
      subroutine cartwgtavgcelldoub1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0
      integer ratio(0:1-1)
      double precision
     &  dxf(0:1-1),
     &  dxc(0:1-1)
      double precision
     &  arrayf(filo0:fihi0),
     &  arrayc(cilo0:cihi0)
      double precision dVf,dVc
      integer ic0,if0,ir0
c
c***********************************************************************
c
      dVf = dxf(0)
      dVc = dxc(0)

      do ic0=ifirstc0,ilastc0
         arrayc(ic0) = zero
      enddo

      do ir0=0,ratio(0)-1
         do ic0=ifirstc0,ilastc0
            if0=ic0*ratio(0)+ir0
            arrayc(ic0)=arrayc(ic0)+arrayf(if0)*dVf
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         arrayc(ic0)=arrayc(ic0)/dVc
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 1d cell-centered float data
c***********************************************************************
c
      subroutine cartwgtavgcellflot1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0
      integer ratio(0:1-1)
      double precision
     &  dxf(0:1-1),
     &  dxc(0:1-1)
      real
     &  arrayf(filo0:fihi0),
     &  arrayc(cilo0:cihi0)
      double precision dVf,dVc
      integer ic0,if0,ir0
c
c***********************************************************************
c
      dVf = dxf(0)
      dVc = dxc(0)

      do ic0=ifirstc0,ilastc0
         arrayc(ic0) = zero
      enddo

      do ir0=0,ratio(0)-1
         do ic0=ifirstc0,ilastc0
            if0=ic0*ratio(0)+ir0
            arrayc(ic0)=arrayc(ic0)+arrayf(if0)*dVf
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         arrayc(ic0)=arrayc(ic0)/dVc
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 1d cell-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgcellcplx1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0
      integer ratio(0:1-1)
      double precision
     &  dxf(0:1-1),
     &  dxc(0:1-1)
      double complex
     &  arrayf(filo0:fihi0),
     &  arrayc(cilo0:cihi0)
      double precision dVf,dVc
      integer ic0,if0,ir0
c
c***********************************************************************
c
      dVf = dxf(0)
      dVc = dxc(0)

      do ic0=ifirstc0,ilastc0
         arrayc(ic0) = cmplx(zero, zero)
      enddo

      do ir0=0,ratio(0)-1
         do ic0=ifirstc0,ilastc0
            if0=ic0*ratio(0)+ir0
            arrayc(ic0)=arrayc(ic0)+arrayf(if0)*dVf
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         arrayc(ic0)=arrayc(ic0)/dVc
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 1d edge-centered double data
c***********************************************************************
c
      subroutine cartwgtavgedgedoub1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0
      integer ratio(0:1-1)
      double precision
     &  dxf(0:1-1),
     &  dxc(0:1-1)
      double precision
     &  arrayf(filo0:fihi0),
     &  arrayc(cilo0:cihi0)
      double precision dVf,dVc
      integer ic0,if0,ir0
c
c***********************************************************************
c
      dVf = dxf(0)
      dVc = dxc(0)

      do ic0=ifirstc0,ilastc0
         arrayc(ic0) = zero
      enddo

      do ir0=0,ratio(0)-1
         do ic0=ifirstc0,ilastc0
            if0=ic0*ratio(0)+ir0
            arrayc(ic0)=arrayc(ic0)+arrayf(if0)*dVf
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         arrayc(ic0)=arrayc(ic0)/dVc
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 1d edge-centered float data
c***********************************************************************
c
      subroutine cartwgtavgedgeflot1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0
      integer ratio(0:1-1)
      double precision
     &  dxf(0:1-1),
     &  dxc(0:1-1)
      real
     &  arrayf(filo0:fihi0),
     &  arrayc(cilo0:cihi0)
      double precision dVf,dVc
      integer ic0,if0,ir0
c
c***********************************************************************
c
      dVf = dxf(0)
      dVc = dxc(0)

      do ic0=ifirstc0,ilastc0
         arrayc(ic0) = zero
      enddo

      do ir0=0,ratio(0)-1
         do ic0=ifirstc0,ilastc0
            if0=ic0*ratio(0)+ir0
            arrayc(ic0)=arrayc(ic0)+arrayf(if0)*dVf
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         arrayc(ic0)=arrayc(ic0)/dVc
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 1d edge-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgedgecplx1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0
      integer ratio(0:1-1)
      double precision
     &  dxf(0:1-1),
     &  dxc(0:1-1)
      double complex
     &  arrayf(filo0:fihi0),
     &  arrayc(cilo0:cihi0)
      double precision dVf,dVc
      integer ic0,if0,ir0
c
c***********************************************************************
c
      dVf = dxf(0)
      dVc = dxc(0)

      do ic0=ifirstc0,ilastc0
         arrayc(ic0) = cmplx(zero, zero)
      enddo

      do ir0=0,ratio(0)-1
         do ic0=ifirstc0,ilastc0
            if0=ic0*ratio(0)+ir0
            arrayc(ic0)=arrayc(ic0)+arrayf(if0)*dVf
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         arrayc(ic0)=arrayc(ic0)/dVc
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 1d face-centered double data
c***********************************************************************
c
      subroutine cartwgtavgfacedoub1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0
      integer ratio(0:1-1)
      double precision
     &  dxf(0:1-1),
     &  dxc(0:1-1)
      double precision
     &  arrayf(filo0:fihi0+1),
     &  arrayc(cilo0:cihi0+1)
      integer ie0
c
c***********************************************************************
c
      do ie0=ifirstc0,ilastc0+1
         arrayc(ie0)=arrayf(ie0*ratio(0))
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 1d face-centered float data
c***********************************************************************
c
      subroutine cartwgtavgfaceflot1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0
      integer ratio(0:1-1)
      double precision
     &  dxf(0:1-1),
     &  dxc(0:1-1)
      real
     &  arrayf(filo0:fihi0+1),
     &  arrayc(cilo0:cihi0+1)
      integer ie0
c
c***********************************************************************
c
      do ie0=ifirstc0,ilastc0+1
         arrayc(ie0)=arrayf(ie0*ratio(0))
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 1d face-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgfacecplx1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0
      integer ratio(0:1-1)
      double precision
     &  dxf(0:1-1),
     &  dxc(0:1-1)
      double complex
     &  arrayf(filo0:fihi0+1),
     &  arrayc(cilo0:cihi0+1)
      integer ie0
c
c***********************************************************************
c
      do ie0=ifirstc0,ilastc0+1
         arrayc(ie0)=arrayf(ie0*ratio(0))
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 1d outerface double data
c***********************************************************************
c
      subroutine cartwgtavgoutfacedoub1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0
      integer ratio(0:1-1)
      double precision
     &  dxf(0:1-1),
     &  dxc(0:1-1)
      double precision
     &  arrayf(1),
     &  arrayc(1)
c
c***********************************************************************
c
      arrayc(1)=arrayf(1)
c
      return
      end    
c
c***********************************************************************
c Weighted averaging for 1d outerface float data
c***********************************************************************
c
      subroutine cartwgtavgoutfaceflot1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0
      integer ratio(0:1-1)
      double precision
     &  dxf(0:1-1),
     &  dxc(0:1-1)
      real
     &  arrayf(1),
     &  arrayc(1)
c
c***********************************************************************
c
      arrayc(1)=arrayf(1)
c
      return
      end    
c
c***********************************************************************
c Weighted averaging for 1d outerface complex data
c***********************************************************************
c
      subroutine cartwgtavgoutfacecplx1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0
      integer ratio(0:1-1)
      double precision
     &  dxf(0:1-1),
     &  dxc(0:1-1)
      double complex
     &  arrayf(1),
     &  arrayc(1)
c
c***********************************************************************
c
      arrayc(1)=arrayf(1)
c
      return
      end    
c
c***********************************************************************
c Weighted averaging for 1d outerside double data
c***********************************************************************
c
      subroutine cartwgtavgoutsidedoub1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0
      integer ratio(0:1-1)
      double precision
     &  dxf(0:1-1),
     &  dxc(0:1-1)
      double precision
     &  arrayf(1),
     &  arrayc(1)
c
c***********************************************************************
c
      arrayc(1)=arrayf(1)
c
      return
      end    
c
c***********************************************************************
c Weighted averaging for 1d side-centered double data
c***********************************************************************
c
      subroutine cartwgtavgsidedoub1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0
      integer ratio(0:1-1)
      double precision
     &  dxf(0:1-1),
     &  dxc(0:1-1)
      double precision
     &  arrayf(filo0:fihi0+1),
     &  arrayc(cilo0:cihi0+1)
      integer ie0
c
c***********************************************************************
c
      do ie0=ifirstc0,ilastc0+1
         arrayc(ie0)=arrayf(ie0*ratio(0))
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 1d side-centered float data
c***********************************************************************
c
      subroutine cartwgtavgsideflot1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0
      integer ratio(0:1-1)
      double precision
     &  dxf(0:1-1),
     &  dxc(0:1-1)
      real
     &  arrayf(filo0:fihi0+1),
     &  arrayc(cilo0:cihi0+1)
      integer ie0
c
c***********************************************************************
c
      do ie0=ifirstc0,ilastc0+1
         arrayc(ie0)=arrayf(ie0*ratio(0))
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 1d side-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgsidecplx1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0
      integer ratio(0:1-1)
      double precision
     &  dxf(0:1-1),
     &  dxc(0:1-1)
      double complex
     &  arrayf(filo0:fihi0+1),
     &  arrayc(cilo0:cihi0+1)
      integer ie0
c
c***********************************************************************
c
      do ie0=ifirstc0,ilastc0+1
         arrayc(ie0)=arrayf(ie0*ratio(0))
      enddo
c
      return
      end
c
