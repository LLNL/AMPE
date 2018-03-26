c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: FORTRAN routines for spatial coarsening of 1d patch data
c               on a regular Cartesian mesh.
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for 1d constant coarsen operators
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
c***********************************************************************
c Constant coarsening for 1d node-centered double data
c***********************************************************************
c
      subroutine conavgnodedoub1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,
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
c Constant coarsening for 1d node-centered float data 
c***********************************************************************
c
      subroutine conavgnodeflot1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,
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
c Constant coarsening for 1d node-centered complex data
c***********************************************************************
c
      subroutine conavgnodecplx1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,
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
c Constant coarsening for 1d node-centered integer data
c***********************************************************************
c
      subroutine conavgnodeintg1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,
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
      integer
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
c Constant coarsening for 1d outernode-centered double data
c***********************************************************************
c
      subroutine conavgouternodedoub1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,
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
     &  arrayf(1),
     &  arrayc(1)
c
c***********************************************************************
c
      arrayc(1)=arrayf(1)

      return
      end

c***********************************************************************
c Constant coarsening for 1d outernode-centered float data
c***********************************************************************
c
      subroutine conavgouternodeflot1d0(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,
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
      real
     &  arrayf(1),
     &  arrayc(1)
c
c***********************************************************************
c
      arrayc(1)=arrayf(1)

      return
      end

c
c***********************************************************************
c Constant coarsening for 1d outernode-centered complex data
c***********************************************************************
c
      subroutine conavgouternodecplx1d0(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,
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
      complex
     &  arrayf(1),
     &  arrayc(1)
c
c***********************************************************************
c
      arrayc(1)=arrayf(1)

      return
      end

c
c***********************************************************************
c Constant coarsening for 1d outernode-centered integer data
c***********************************************************************
c
      subroutine conavgouternodeint1d0(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,
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
      integer
     &  arrayf(1),
     &  arrayc(1)
c
c***********************************************************************
c
      arrayc(1)=arrayf(1)

      return
      end



