c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: FORTRAN routines for spatial coarsening of 2d patch data
c               on a regular Cartesian mesh.
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for 2d Cartesian coarsen operators
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
c***********************************************************************
c Constant coarsening for 2d node-centered double data
c***********************************************************************
c
      subroutine conavgnodedoub2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1)
      integer ie0,ie1,if1
c
c***********************************************************************
c
      do ie1=ifirstc1,ilastc1+1
         if1=ie1*ratio(1)
         do ie0=ifirstc0,ilastc0+1
            arrayc(ie0,ie1)=arrayf(ie0*ratio(0),if1)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant coarsening for 2d node-centered float data 
c***********************************************************************
c
      subroutine conavgnodeflot2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      real
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1)
      integer ie0,ie1,if1
c
c***********************************************************************
c
      do ie1=ifirstc1,ilastc1+1
         if1=ie1*ratio(1)
         do ie0=ifirstc0,ilastc0+1
            arrayc(ie0,ie1)=arrayf(ie0*ratio(0),if1)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant coarsening for 2d node-centered complex data
c***********************************************************************
c
      subroutine conavgnodecplx2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double complex
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1)
      integer ie0,ie1,if1
c
c***********************************************************************
c
      do ie1=ifirstc1,ilastc1+1
         if1=ie1*ratio(1)
         do ie0=ifirstc0,ilastc0+1
            arrayc(ie0,ie1)=arrayf(ie0*ratio(0),if1)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant coarsening for 2d node-centered integer data
c***********************************************************************
c
      subroutine conavgnodeintg2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      integer
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1)
      integer ie0,ie1,if1
c
c***********************************************************************
c
      do ie1=ifirstc1,ilastc1+1
         if1=ie1*ratio(1)
         do ie0=ifirstc0,ilastc0+1
            arrayc(ie0,ie1)=arrayf(ie0*ratio(0),if1)
         enddo
      enddo
c
      return
      end
c
c

c***********************************************************************
c Constant coarsening for 2d outernode-centered double data
c***********************************************************************
c
      subroutine conavgouternodedoub2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  arrayf(filo1+1:fihi1),
     &  arrayc(cilo1+1:cihi1)
      integer ic1,if1
c
c***********************************************************************
c
      do ic1=ifirstc1+1,ilastc1
         if1=ic1*ratio(1)
         arrayc(ic1)=arrayf(if1)
      enddo

      return
      end

      subroutine conavgouternodedoub2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      double precision
     &  arrayf(filo0:fihi0+1),
     &  arrayc(cilo0:cihi0+1)
      integer ic0,if0
c
c***********************************************************************
c
      do ic0=ifirstc0,ilastc0+1
         if0=ic0*ratio(0)
         arrayc(ic0)=arrayf(if0)
      enddo

      return
      end

c***********************************************************************
c Constant coarsening for 2d outernode-centered float data
c***********************************************************************
c
      subroutine conavgouternodeflot2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      real
     &  arrayf(filo1+1:fihi1),
     &  arrayc(cilo1+1:cihi1)
      integer ic1,if1
c
c***********************************************************************
c
      do ic1=ifirstc1+1,ilastc1
         if1=ic1*ratio(1)
         arrayc(ic1)=arrayf(if1)
      enddo

      return
      end

      subroutine conavgouternodeflot2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      real
     &  arrayf(filo0:fihi0+1),
     &  arrayc(cilo0:cihi0+1)
      integer ic0,if0
c
c***********************************************************************
c
      do ic0=ifirstc0,ilastc0+1
         if0=ic0*ratio(0)
         arrayc(ic0)=arrayf(if0)
      enddo

      return
      end

c
c***********************************************************************
c Constant coarsening for 2d outernode-centered complex data
c***********************************************************************
c
      subroutine conavgouternodecplx2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      complex
     &  arrayf(filo1+1:fihi1),
     &  arrayc(cilo1+1:cihi1)
      integer ic1,if1
c
c***********************************************************************
c
      do ic1=ifirstc1+1,ilastc1
         if1=ic1*ratio(1)
         arrayc(ic1)=arrayf(if1)
      enddo

      return
      end

      subroutine conavgouternodecplx2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      complex
     &  arrayf(filo0:fihi0+1),
     &  arrayc(cilo0:cihi0+1)
      integer ic0,if0
c
c***********************************************************************
c
      do ic0=ifirstc0,ilastc0+1
         if0=ic0*ratio(0)
         arrayc(ic0)=arrayf(if0)
      enddo

      return
      end

c
c***********************************************************************
c Constant coarsening for 2d outernode-centered integer data
c***********************************************************************
c
      subroutine conavgouternodeint2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      integer
     &  arrayf(filo1+1:fihi1),
     &  arrayc(cilo1+1:cihi1)
      integer ic1,if1
c
c***********************************************************************
c
      do ic1=ifirstc1+1,ilastc1
         if1=ic1*ratio(1)
         arrayc(ic1)=arrayf(if1)
      enddo

      return
      end

      subroutine conavgouternodeint2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:2-1)
      integer
     &  arrayf(filo0:fihi0+1),
     &  arrayc(cilo0:cihi0+1)
      integer ic0,if0
c
c***********************************************************************
c
      do ic0=ifirstc0,ilastc0+1
         if0=ic0*ratio(0)
         arrayc(ic0)=arrayf(if0)
      enddo

      return
      end

c

