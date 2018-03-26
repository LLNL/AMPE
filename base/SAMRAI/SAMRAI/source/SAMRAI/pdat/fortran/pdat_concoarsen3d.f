c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: FORTRAN routines for spatial coarsening of 3d patch data
c               on a regular Cartesian mesh.
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for 3d Cartesian coarsen operators
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
c***********************************************************************
c Constant coarsening for 3d node-centered double data
c***********************************************************************
c
      subroutine conavgnodedoub3d(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2
      integer ratio(0:3-1)
      double precision
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1,
     &          filo2:fihi2+1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1)
      integer ie0,ie1,ie2,if1,if2
c
c***********************************************************************
c
      do ie2=ifirstc2,ilastc2+1
         if2=ie2*ratio(2)
         do ie1=ifirstc1,ilastc1+1
            if1=ie1*ratio(1)
            do ie0=ifirstc0,ilastc0+1
               arrayc(ie0,ie1,ie2)=arrayf(ie0*ratio(0),if1,if2)
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant coarsening for 3d node-centered float data 
c***********************************************************************
c
      subroutine conavgnodeflot3d(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2
      integer ratio(0:3-1)
      real
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1,
     &          filo2:fihi2+1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1)
      integer ie0,ie1,ie2,if1,if2
c
c***********************************************************************
c
      do ie2=ifirstc2,ilastc2+1
         if2=ie2*ratio(2)
         do ie1=ifirstc1,ilastc1+1
            if1=ie1*ratio(1)
            do ie0=ifirstc0,ilastc0+1
               arrayc(ie0,ie1,ie2)=arrayf(ie0*ratio(0),if1,if2)
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant coarsening for 3d node-centered complex data
c***********************************************************************
c
      subroutine conavgnodecplx3d(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2
      integer ratio(0:3-1)
      double complex
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1,
     &          filo2:fihi2+1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1)
      integer ie0,ie1,ie2,if1,if2
c
c***********************************************************************
c
      do ie2=ifirstc2,ilastc2+1
         if2=ie2*ratio(2)
         do ie1=ifirstc1,ilastc1+1
            if1=ie1*ratio(1)
            do ie0=ifirstc0,ilastc0+1
               arrayc(ie0,ie1,ie2)=arrayf(ie0*ratio(0),if1,if2)
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant coarsening for 3d node-centered integer data
c***********************************************************************
c
      subroutine conavgnodeintg3d(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2
      integer ratio(0:3-1)
      integer
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1,
     &          filo2:fihi2+1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1)
      integer ie0,ie1,ie2,if1,if2
c
c***********************************************************************
c
      do ie2=ifirstc2,ilastc2+1
         if2=ie2*ratio(2)
         do ie1=ifirstc1,ilastc1+1
            if1=ie1*ratio(1)
            do ie0=ifirstc0,ilastc0+1
               arrayc(ie0,ie1,ie2)=arrayf(ie0*ratio(0),if1,if2)
            enddo
         enddo
      enddo
c
      return
      end
c

c***********************************************************************
c Constant coarsening for 3d outernode-centered double data
c***********************************************************************
c
      subroutine conavgouternodedoub3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2
      integer ratio(0:3-1)
      double precision
     &  arrayf(filo1+1:fihi1,
     &          filo2+1:fihi2),
     &  arrayc(cilo1+1:cihi1,
     &          cilo2+1:cihi2)
      integer ic1,ic2,if1,if2
c
c***********************************************************************
c
      do ic2=ifirstc2+1,ilastc2
         if2=ic2*ratio(2)
         do ic1=ifirstc1+1,ilastc1
            if1=ic1*ratio(1)
            arrayc(ic1,ic2)=arrayf(if1,if2)
         enddo
      enddo

      return
      end

      subroutine conavgouternodedoub3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2
      integer ratio(0:3-1)
      double precision
     &  arrayf(filo0:fihi0+1,
     &          filo2+1:fihi2),
     &  arrayc(cilo0:cihi0+1,
     &          cilo2+1:cihi2)
      integer ic0,ic2,if0,if2
c
c***********************************************************************
c
      do ic2=ifirstc2+1,ilastc2
         if2=ic2*ratio(2)
         do ic0=ifirstc0,ilastc0+1
            if0=ic0*ratio(0)
            arrayc(ic0,ic2)=arrayf(if0,if2)
         enddo
      enddo

      return
      end

      subroutine conavgouternodedoub3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2
      integer ratio(0:3-1)
      double precision
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do ic1=ifirstc1,ilastc1+1
         if1=ic1*ratio(1)
         do ic0=ifirstc0,ilastc0+1
            if0=ic0*ratio(0)
            arrayc(ic0,ic1)=arrayf(if0,if1)
         enddo
      enddo

      return
      end

c***********************************************************************
c Constant coarsening for 3d outernode-centered float data
c***********************************************************************
c
      subroutine conavgouternodeflot3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2
      integer ratio(0:3-1)
      real
     &  arrayf(filo1+1:fihi1,
     &          filo2+1:fihi2),
     &  arrayc(cilo1+1:cihi1,
     &          cilo2+1:cihi2)
      integer ic1,ic2,if1,if2
c
c***********************************************************************
c
      do ic2=ifirstc2+1,ilastc2
         if2=ic2*ratio(2)
         do ic1=ifirstc1+1,ilastc1
            if1=ic1*ratio(1)
            arrayc(ic1,ic2)=arrayf(if1,if2)
         enddo
      enddo

      return
      end

      subroutine conavgouternodeflot3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2
      integer ratio(0:3-1)
      real
     &  arrayf(filo0:fihi0+1,
     &          filo2+1:fihi2),
     &  arrayc(cilo0:cihi0+1,
     &          cilo2+1:cihi2)
      integer ic0,ic2,if0,if2
c
c***********************************************************************
c
      do ic2=ifirstc2+1,ilastc2
         if2=ic2*ratio(2)
         do ic0=ifirstc0,ilastc0+1
            if0=ic0*ratio(0)
            arrayc(ic0,ic2)=arrayf(if0,if2)
         enddo
      enddo

      return
      end

      subroutine conavgouternodeflot3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2
      integer ratio(0:3-1)
      real
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do ic1=ifirstc1,ilastc1+1
         if1=ic1*ratio(1)
         do ic0=ifirstc0,ilastc0+1
            if0=ic0*ratio(0)
            arrayc(ic0,ic1)=arrayf(if0,if1)
         enddo
      enddo

      return
      end

c
c***********************************************************************
c Constant coarsening for 3d outernode-centered complex data
c***********************************************************************
c
      subroutine conavgouternodecplx3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2
      integer ratio(0:3-1)
      complex
     &  arrayf(filo1+1:fihi1,
     &          filo2+1:fihi2),
     &  arrayc(cilo1+1:cihi1,
     &          cilo2+1:cihi2)
      integer ic1,ic2,if1,if2
c
c***********************************************************************
c
      do ic2=ifirstc2+1,ilastc2
         if2=ic2*ratio(2)
         do ic1=ifirstc1+1,ilastc1
            if1=ic1*ratio(1)
            arrayc(ic1,ic2)=arrayf(if1,if2)
         enddo
      enddo

      return
      end

      subroutine conavgouternodecplx3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2
      integer ratio(0:3-1)
      complex
     &  arrayf(filo0:fihi0+1,
     &          filo2+1:fihi2),
     &  arrayc(cilo0:cihi0+1,
     &          cilo2+1:cihi2)
      integer ic0,ic2,if0,if2
c
c***********************************************************************
c
      do ic2=ifirstc2+1,ilastc2
         if2=ic2*ratio(2)
         do ic0=ifirstc0,ilastc0+1
            if0=ic0*ratio(0)
            arrayc(ic0,ic2)=arrayf(if0,if2)
         enddo
      enddo

      return
      end

      subroutine conavgouternodecplx3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2
      integer ratio(0:3-1)
      complex
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do ic1=ifirstc1,ilastc1+1
         if1=ic1*ratio(1)
         do ic0=ifirstc0,ilastc0+1
            if0=ic0*ratio(0)
            arrayc(ic0,ic1)=arrayf(if0,if1)
         enddo
      enddo

      return
      end


c
c***********************************************************************
c Constant coarsening for 3d outernode-centered integer data
c***********************************************************************
c
      subroutine conavgouternodeint3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2
      integer ratio(0:3-1)
      integer
     &  arrayf(filo1+1:fihi1,
     &          filo2+1:fihi2),
     &  arrayc(cilo1+1:cihi1,
     &          cilo2+1:cihi2)
      integer ic1,ic2,if1,if2
c
c***********************************************************************
c
      do ic2=ifirstc2+1,ilastc2
         if2=ic2*ratio(2)
         do ic1=ifirstc1+1,ilastc1
            if1=ic1*ratio(1)
            arrayc(ic1,ic2)=arrayf(if1,if2)
         enddo
      enddo

      return
      end

      subroutine conavgouternodeint3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2
      integer ratio(0:3-1)
      integer
     &  arrayf(filo0:fihi0+1,
     &          filo2+1:fihi2),
     &  arrayc(cilo0:cihi0+1,
     &          cilo2+1:cihi2)
      integer ic0,ic2,if0,if2
c
c***********************************************************************
c
      do ic2=ifirstc2+1,ilastc2
         if2=ic2*ratio(2)
         do ic0=ifirstc0,ilastc0+1
            if0=ic0*ratio(0)
            arrayc(ic0,ic2)=arrayf(if0,if2)
         enddo
      enddo

      return
      end

      subroutine conavgouternodeint3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2
      integer ratio(0:3-1)
      integer
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do ic1=ifirstc1,ilastc1+1
         if1=ic1*ratio(1)
         do ic0=ifirstc0,ilastc0+1
            if0=ic0*ratio(0)
            arrayc(ic0,ic1)=arrayf(if0,if1)
         enddo
      enddo

      return
      end



