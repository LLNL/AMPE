c
c  File:        $URL$
c  Package:     SAMRAI geometry
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: FORTRAN routines for spatial coarsening of 2d patch data
c               on a regular Cartesian mesh.
c
c
c  File:        $URL$
c  Package:     SAMRAI geometry
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
c
c
c
c
c***********************************************************************
c Weighted averaging for 2d cell-centered double data
c***********************************************************************
c
      subroutine cartwgtavgcelldoub2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      double precision
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1)
      double precision dVf,dVc
      integer ic0,ic1,if0,if1,ir0,ir1
c
c***********************************************************************
c
      dVf = dxf(0)*dxf(1)
      dVc = dxc(0)*dxc(1)

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            arrayc(ic0,ic1)=zero
         enddo
      enddo

      do ir1=0,ratio(1)-1
         do ir0=0,ratio(0)-1
            do ic1=ifirstc1,ilastc1
               if1=ic1*ratio(1)+ir1
               do ic0=ifirstc0,ilastc0
                  if0=ic0*ratio(0)+ir0
                  arrayc(ic0,ic1)=arrayc(ic0,ic1)
     &                            +arrayf(if0,if1)*dVf
               enddo
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            arrayc(ic0,ic1)=arrayc(ic0,ic1)/dVc
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 2d cell-centered float data
c***********************************************************************
c
      subroutine cartwgtavgcellflot2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      real
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1)
      double precision dVf,dVc
      integer ic0,ic1,if0,if1,ir0,ir1
c
c***********************************************************************
c
      dVf = dxf(0)*dxf(1)
      dVc = dxc(0)*dxc(1)

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            arrayc(ic0,ic1)=zero
         enddo
      enddo

      do ir1=0,ratio(1)-1
         do ir0=0,ratio(0)-1
            do ic1=ifirstc1,ilastc1
               if1=ic1*ratio(1)+ir1
               do ic0=ifirstc0,ilastc0
                  if0=ic0*ratio(0)+ir0
                  arrayc(ic0,ic1)=arrayc(ic0,ic1)
     &                            +arrayf(if0,if1)*dVf
               enddo
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            arrayc(ic0,ic1)=arrayc(ic0,ic1)/dVc
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 2d cell-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgcellcplx2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      double complex
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1)
      double precision dVf,dVc
      integer ic0,ic1,if0,if1,ir0,ir1
c
c***********************************************************************
c
      dVf = dxf(0)*dxf(1)
      dVc = dxc(0)*dxc(1)

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            arrayc(ic0,ic1)=cmplx(zero,zero)
         enddo
      enddo

      do ir1=0,ratio(1)-1
         do ir0=0,ratio(0)-1
            do ic1=ifirstc1,ilastc1
               if1=ic1*ratio(1)+ir1
               do ic0=ifirstc0,ilastc0
                  if0=ic0*ratio(0)+ir0
                  arrayc(ic0,ic1)=arrayc(ic0,ic1)
     &                            +arrayf(if0,if1)*dVf
               enddo
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            arrayc(ic0,ic1)=arrayc(ic0,ic1)/dVc
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 2d edge-centered double data
c***********************************************************************
c
      subroutine cartwgtavgedgedoub2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      double precision
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1)
      double precision lengthf, lengthc
      integer ic0,ic1,if0,if1,ir
c
c***********************************************************************
c
      lengthf=dxf(0)
      lengthc=dxc(0)

      do ic1=ifirstc1,ilastc1+1
         do ic0=ifirstc0,ilastc0+0
            arrayc(ic0,ic1)=zero
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+1

            if1=ic1*ratio(1)
            do ic0=ifirstc0,ilastc0+0

               do ir=0,ratio(0)-1
                  if0=ic0*ratio(0)+ir
                  arrayc(ic0,ic1)=arrayc(ic0,ic1)
     &                           +arrayf(if0,if1)*lengthf
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+1
         do ic0=ifirstc0,ilastc0+0
            arrayc(ic0,ic1)=arrayc(ic0,ic1)/lengthc
        enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgedgedoub2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      double precision
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1)
      double precision lengthf, lengthc
      integer ic0,ic1,if0,if1,ir
c
c***********************************************************************
c
      lengthf=dxf(1)
      lengthc=dxc(1)

      do ic1=ifirstc1,ilastc1+0
         do ic0=ifirstc0,ilastc0+1
            arrayc(ic0,ic1)=zero
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+0

         do ir=0,ratio(1)-1
            if1=ic1*ratio(1)+ir
            do ic0=ifirstc0,ilastc0+1

                  if0=ic0*ratio(0)
                  arrayc(ic0,ic1)=arrayc(ic0,ic1)
     &                           +arrayf(if0,if1)*lengthf
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+0
         do ic0=ifirstc0,ilastc0+1
            arrayc(ic0,ic1)=arrayc(ic0,ic1)/lengthc
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 2d edge-centered float data
c***********************************************************************
c
      subroutine cartwgtavgedgeflot2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      real
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1)
      double precision lengthf, lengthc
      integer ic0,ic1,if0,if1,ir
c
c***********************************************************************
c
      lengthf=dxf(0)
      lengthc=dxc(0)

      do ic1=ifirstc1,ilastc1+1
         do ic0=ifirstc0,ilastc0+0
            arrayc(ic0,ic1)=zero
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+1

            if1=ic1*ratio(1)
            do ic0=ifirstc0,ilastc0+0

               do ir=0,ratio(0)-1
                  if0=ic0*ratio(0)+ir
                  arrayc(ic0,ic1)=arrayc(ic0,ic1)
     &                           +arrayf(if0,if1)*lengthf
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+1
         do ic0=ifirstc0,ilastc0+0
            arrayc(ic0,ic1)=arrayc(ic0,ic1)/lengthc
        enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgedgeflot2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      real
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1)
      double precision lengthf, lengthc
      integer ic0,ic1,if0,if1,ir
c
c***********************************************************************
c
      lengthf=dxf(1)
      lengthc=dxc(1)

      do ic1=ifirstc1,ilastc1+0
         do ic0=ifirstc0,ilastc0+1
            arrayc(ic0,ic1)=zero
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+0

         do ir=0,ratio(1)-1
            if1=ic1*ratio(1)+ir
            do ic0=ifirstc0,ilastc0+1

                  if0=ic0*ratio(0)
                  arrayc(ic0,ic1)=arrayc(ic0,ic1)
     &                           +arrayf(if0,if1)*lengthf
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+0
         do ic0=ifirstc0,ilastc0+1
            arrayc(ic0,ic1)=arrayc(ic0,ic1)/lengthc
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 2d edge-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgedgecplx2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      double complex
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1)
      double precision lengthf, lengthc
      integer ic0,ic1,if0,if1,ir
c
c***********************************************************************
c
      lengthf=dxf(0)
      lengthc=dxc(0)

      do ic1=ifirstc1,ilastc1+1
         do ic0=ifirstc0,ilastc0+0
            arrayc(ic0,ic1)=cmplx(zero,zero)
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+1

            if1=ic1*ratio(1)
            do ic0=ifirstc0,ilastc0+0

               do ir=0,ratio(0)-1
                  if0=ic0*ratio(0)+ir
                  arrayc(ic0,ic1)=arrayc(ic0,ic1)
     &                           +arrayf(if0,if1)*lengthf
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+1
         do ic0=ifirstc0,ilastc0+0
            arrayc(ic0,ic1)=arrayc(ic0,ic1)/lengthc
        enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgedgecplx2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      double complex
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1)
      double precision lengthf, lengthc
      integer ic0,ic1,if0,if1,ir
c
c***********************************************************************
c
      lengthf=dxf(1)
      lengthc=dxc(1)

      do ic1=ifirstc1,ilastc1+0
         do ic0=ifirstc0,ilastc0+1
            arrayc(ic0,ic1)=cmplx(zero,zero)
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+0

         do ir=0,ratio(1)-1
            if1=ic1*ratio(1)+ir
            do ic0=ifirstc0,ilastc0+1

                  if0=ic0*ratio(0)
                  arrayc(ic0,ic1)=arrayc(ic0,ic1)
     &                           +arrayf(if0,if1)*lengthf
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+0
         do ic0=ifirstc0,ilastc0+1
            arrayc(ic0,ic1)=arrayc(ic0,ic1)/lengthc
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 2d face-centered double data
c***********************************************************************
c
      subroutine cartwgtavgfacedoub2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      double precision
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1)
      double precision lengthf, lengthc
      integer ie0,ic1,if0,if1,ir1
c
c***********************************************************************
c
      lengthf=dxf(1)
      lengthc=dxc(1)

      do ic1=ifirstc1,ilastc1
         do ie0=ifirstc0,ilastc0+1
            arrayc(ie0,ic1)=zero
         enddo
      enddo

      do ir1=0,ratio(1)-1
         do ic1=ifirstc1,ilastc1
            if1=ic1*ratio(1)+ir1
            do ie0=ifirstc0,ilastc0+1
               if0=ie0*ratio(0)
               arrayc(ie0,ic1)=arrayc(ie0,ic1)
     &                           +arrayf(if0,if1)*lengthf
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ie0=ifirstc0,ilastc0+1
            arrayc(ie0,ic1)=arrayc(ie0,ic1)/lengthc
        enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgfacedoub2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      double precision
     &  arrayf(filo1:fihi1+1,
     &          filo0:fihi0),
     &  arrayc(cilo1:cihi1+1,
     &          cilo0:cihi0)
      double precision lengthf, lengthc
      integer ie1,ic0,if1,if0,ir0
c
c***********************************************************************
c
      lengthf=dxf(0)
      lengthc=dxc(0)

      do ic0=ifirstc0,ilastc0
         do ie1=ifirstc1,ilastc1+1
            arrayc(ie1,ic0)=zero
         enddo
      enddo

      do ir0=0,ratio(0)-1
         do ic0=ifirstc0,ilastc0
            if0=ic0*ratio(0)+ir0
            do ie1=ifirstc1,ilastc1+1
               if1=ie1*ratio(1)
               arrayc(ie1,ic0)=arrayc(ie1,ic0)
     &                           +arrayf(if1,if0)*lengthf
            enddo
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         do ie1=ifirstc1,ilastc1+1
            arrayc(ie1,ic0)=arrayc(ie1,ic0)/lengthc
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 2d face-centered float data
c***********************************************************************
c
      subroutine cartwgtavgfaceflot2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      real
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1)
      double precision lengthf, lengthc
      integer ie0,ic1,if0,if1,ir1
c
c***********************************************************************
c
      lengthf=dxf(1)
      lengthc=dxc(1)

      do ic1=ifirstc1,ilastc1
         do ie0=ifirstc0,ilastc0+1
            arrayc(ie0,ic1)=zero
         enddo
      enddo

      do ir1=0,ratio(1)-1
         do ic1=ifirstc1,ilastc1
            if1=ic1*ratio(1)+ir1
            do ie0=ifirstc0,ilastc0+1
               if0=ie0*ratio(0)
               arrayc(ie0,ic1)=arrayc(ie0,ic1)
     &                           +arrayf(if0,if1)*lengthf
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ie0=ifirstc0,ilastc0+1
            arrayc(ie0,ic1)=arrayc(ie0,ic1)/lengthc
        enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgfaceflot2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      real
     &  arrayf(filo1:fihi1+1,
     &          filo0:fihi0),
     &  arrayc(cilo1:cihi1+1,
     &          cilo0:cihi0)
      double precision lengthf, lengthc
      integer ie1,ic0,if1,if0,ir0
c
c***********************************************************************
c
      lengthf=dxf(0)
      lengthc=dxc(0)

      do ic0=ifirstc0,ilastc0
         do ie1=ifirstc1,ilastc1+1
            arrayc(ie1,ic0)=zero
         enddo
      enddo

      do ir0=0,ratio(0)-1
         do ic0=ifirstc0,ilastc0
            if0=ic0*ratio(0)+ir0
            do ie1=ifirstc1,ilastc1+1
               if1=ie1*ratio(1)
               arrayc(ie1,ic0)=arrayc(ie1,ic0)
     &                           +arrayf(if1,if0)*lengthf
            enddo
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         do ie1=ifirstc1,ilastc1+1
            arrayc(ie1,ic0)=arrayc(ie1,ic0)/lengthc
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 2d face-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgfacecplx2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      double complex
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1)
      double precision lengthf, lengthc
      integer ie0,ic1,if0,if1,ir1
c
c***********************************************************************
c
      lengthf=dxf(1)
      lengthc=dxc(1)

      do ic1=ifirstc1,ilastc1
         do ie0=ifirstc0,ilastc0+1
            arrayc(ie0,ic1)=cmplx(zero,zero)
         enddo
      enddo

      do ir1=0,ratio(1)-1
         do ic1=ifirstc1,ilastc1
            if1=ic1*ratio(1)+ir1
            do ie0=ifirstc0,ilastc0+1
               if0=ie0*ratio(0)
               arrayc(ie0,ic1)=arrayc(ie0,ic1)
     &                           +arrayf(if0,if1)*lengthf
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ie0=ifirstc0,ilastc0+1
            arrayc(ie0,ic1)=arrayc(ie0,ic1)/lengthc
        enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgfacecplx2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      double complex
     &  arrayf(filo1:fihi1+1,
     &          filo0:fihi0),
     &  arrayc(cilo1:cihi1+1,
     &          cilo0:cihi0)
      double precision lengthf, lengthc
      integer ie1,ic0,if1,if0,ir0
c
c***********************************************************************
c
      lengthf=dxf(0)
      lengthc=dxc(0)

      do ic0=ifirstc0,ilastc0
         do ie1=ifirstc1,ilastc1+1
            arrayc(ie1,ic0)=cmplx(zero,zero)
         enddo
      enddo

      do ir0=0,ratio(0)-1
         do ic0=ifirstc0,ilastc0
            if0=ic0*ratio(0)+ir0
            do ie1=ifirstc1,ilastc1+1
               if1=ie1*ratio(1)
               arrayc(ie1,ic0)=arrayc(ie1,ic0)
     &                           +arrayf(if1,if0)*lengthf
            enddo
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         do ie1=ifirstc1,ilastc1+1
            arrayc(ie1,ic0)=arrayc(ie1,ic0)/lengthc
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 2d outerface double data
c***********************************************************************
c
      subroutine cartwgtavgoutfacedoub2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      double precision
     &  arrayf(filo1:fihi1),
     &  arrayc(cilo1:cihi1)
      double precision lengthf, lengthc
      integer ic1,if1,ir1
c
c***********************************************************************
c
      lengthf=dxf(1)
      lengthc=dxc(1)

      do ic1=ifirstc1,ilastc1
         arrayc(ic1)=zero
      enddo

      do ir1=0,ratio(1)-1
         do ic1=ifirstc1,ilastc1
            if1=ic1*ratio(1)+ir1
            arrayc(ic1)=arrayc(ic1)+arrayf(if1)*lengthf
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         arrayc(ic1)=arrayc(ic1)/lengthc
      enddo
c
      return
      end
c
      subroutine cartwgtavgoutfacedoub2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      double precision
     &  arrayf(filo0:fihi0),
     &  arrayc(cilo0:cihi0)
      double precision lengthf, lengthc
      integer ic0,if0,ir0
c
c***********************************************************************
c
      lengthf=dxf(0)
      lengthc=dxc(0)

      do ic0=ifirstc0,ilastc0
         arrayc(ic0)=zero
      enddo

      do ir0=0,ratio(0)-1
         do ic0=ifirstc0,ilastc0
            if0=ic0*ratio(0)+ir0
            arrayc(ic0)=arrayc(ic0)+arrayf(if0)*lengthf
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         arrayc(ic0)=arrayc(ic0)/lengthc
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 2d outerface float data
c***********************************************************************
c
      subroutine cartwgtavgoutfaceflot2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      real
     &  arrayf(filo1:fihi1),
     &  arrayc(cilo1:cihi1)
      double precision lengthf, lengthc
      integer ic1,if1,ir1
c
c***********************************************************************
c
      lengthf=dxf(1)
      lengthc=dxc(1)

      do ic1=ifirstc1,ilastc1
         arrayc(ic1)=zero
      enddo

      do ir1=0,ratio(1)-1
         do ic1=ifirstc1,ilastc1
            if1=ic1*ratio(1)+ir1
            arrayc(ic1)=arrayc(ic1)+arrayf(if1)*lengthf
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         arrayc(ic1)=arrayc(ic1)/lengthc
      enddo
c
      return
      end
c
      subroutine cartwgtavgoutfaceflot2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      real
     &  arrayf(filo0:fihi0),
     &  arrayc(cilo0:cihi0)
      double precision lengthf, lengthc
      integer ic0,if0,ir0
c
c***********************************************************************
c
      lengthf=dxf(0)
      lengthc=dxc(0)

      do ic0=ifirstc0,ilastc0
         arrayc(ic0)=zero
      enddo

      do ir0=0,ratio(0)-1
         do ic0=ifirstc0,ilastc0
            if0=ic0*ratio(0)+ir0
            arrayc(ic0)=arrayc(ic0)+arrayf(if0)*lengthf
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         arrayc(ic0)=arrayc(ic0)/lengthc
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 2d outerface complex data
c***********************************************************************
c
      subroutine cartwgtavgoutfacecplx2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      double complex
     &  arrayf(filo1:fihi1),
     &  arrayc(cilo1:cihi1)
      double precision lengthf, lengthc
      integer ic1,if1,ir1
c
c***********************************************************************
c
      lengthf=dxf(1)
      lengthc=dxc(1)

      do ic1=ifirstc1,ilastc1
         arrayc(ic1)=cmplx(zero,zero)
      enddo

      do ir1=0,ratio(1)-1
         do ic1=ifirstc1,ilastc1
            if1=ic1*ratio(1)+ir1
            arrayc(ic1)=arrayc(ic1)+arrayf(if1)*lengthf
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         arrayc(ic1)=arrayc(ic1)/lengthc
      enddo
c
      return
      end
c
      subroutine cartwgtavgoutfacecplx2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      double complex
     &  arrayf(filo0:fihi0),
     &  arrayc(cilo0:cihi0)
      double precision lengthf, lengthc
      integer ic0,if0,ir0
c
c***********************************************************************
c
      lengthf=dxf(0)
      lengthc=dxc(0)

      do ic0=ifirstc0,ilastc0
         arrayc(ic0)=cmplx(zero,zero)
      enddo

      do ir0=0,ratio(0)-1
         do ic0=ifirstc0,ilastc0
            if0=ic0*ratio(0)+ir0
            arrayc(ic0)=arrayc(ic0)+arrayf(if0)*lengthf
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         arrayc(ic0)=arrayc(ic0)/lengthc
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 2d outerside double data
c***********************************************************************
c
      subroutine cartwgtavgoutsidedoub2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      double precision
     &  arrayf(filo1:fihi1),
     &  arrayc(cilo1:cihi1)
      double precision lengthf, lengthc
      integer ic1,if1,ir1
c
c***********************************************************************
c
      lengthf=dxf(1)
      lengthc=dxc(1)

      do ic1=ifirstc1,ilastc1
         arrayc(ic1)=zero
      enddo

      do ir1=0,ratio(1)-1
         do ic1=ifirstc1,ilastc1
            if1=ic1*ratio(1)+ir1
            arrayc(ic1)=arrayc(ic1)+arrayf(if1)*lengthf
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         arrayc(ic1)=arrayc(ic1)/lengthc
      enddo
c
      return
      end
c
      subroutine cartwgtavgoutsidedoub2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      double precision
     &  arrayf(filo0:fihi0),
     &  arrayc(cilo0:cihi0)
      double precision lengthf, lengthc
      integer ic0,if0,ir0
c
c***********************************************************************
c
      lengthf=dxf(0)
      lengthc=dxc(0)

      do ic0=ifirstc0,ilastc0
         arrayc(ic0)=zero
      enddo

      do ir0=0,ratio(0)-1
         do ic0=ifirstc0,ilastc0
            if0=ic0*ratio(0)+ir0
            arrayc(ic0)=arrayc(ic0)+arrayf(if0)*lengthf
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         arrayc(ic0)=arrayc(ic0)/lengthc
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 2d side-centered double data
c***********************************************************************
c
      subroutine cartwgtavgsidedoub2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      double precision
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1)
      double precision lengthf, lengthc
      integer ic0,ic1,if0,if1,ir
c
c***********************************************************************
c
      lengthf=dxf(1)
      lengthc=dxc(1)

      do ic1=ifirstc1,ilastc1+0
         do ic0=ifirstc0,ilastc0+1
            arrayc(ic0,ic1)=zero
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+0

         do ir=0,ratio(1)-1
            if1=ic1*ratio(1)+ir
            do ic0=ifirstc0,ilastc0+1

                  if0=ic0*ratio(0)
               arrayc(ic0,ic1)=arrayc(ic0,ic1)
     &                           +arrayf(if0,if1)*lengthf
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+0
         do ic0=ifirstc0,ilastc0+1
            arrayc(ic0,ic1)=arrayc(ic0,ic1)/lengthc
        enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgsidedoub2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      double precision
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1)
      double precision lengthf, lengthc
      integer ic0,ic1,if0,if1,ir
c
c***********************************************************************
c
      lengthf=dxf(0)
      lengthc=dxc(0)

      do ic1=ifirstc1,ilastc1+1
         do ic0=ifirstc0,ilastc0+0
            arrayc(ic0,ic1)=zero
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+1

            if1=ic1*ratio(1)
            do ic0=ifirstc0,ilastc0+0

               do ir=0,ratio(0)-1
                  if0=ic0*ratio(0)+ir
               arrayc(ic0,ic1)=arrayc(ic0,ic1)
     &                           +arrayf(if0,if1)*lengthf
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+1
         do ic0=ifirstc0,ilastc0+0
            arrayc(ic0,ic1)=arrayc(ic0,ic1)/lengthc
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 2d side-centered float data
c***********************************************************************
c
      subroutine cartwgtavgsideflot2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      real
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1)
      double precision lengthf, lengthc
      integer ic0,ic1,if0,if1,ir
c
c***********************************************************************
c
      lengthf=dxf(1)
      lengthc=dxc(1)

      do ic1=ifirstc1,ilastc1+0
         do ic0=ifirstc0,ilastc0+1
            arrayc(ic0,ic1)=zero
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+0

         do ir=0,ratio(1)-1
            if1=ic1*ratio(1)+ir
            do ic0=ifirstc0,ilastc0+1

                  if0=ic0*ratio(0)
               arrayc(ic0,ic1)=arrayc(ic0,ic1)
     &                           +arrayf(if0,if1)*lengthf
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+0
         do ic0=ifirstc0,ilastc0+1
            arrayc(ic0,ic1)=arrayc(ic0,ic1)/lengthc
        enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgsideflot2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      real
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1)
      double precision lengthf, lengthc
      integer ic0,ic1,if0,if1,ir
c
c***********************************************************************
c
      lengthf=dxf(0)
      lengthc=dxc(0)

      do ic1=ifirstc1,ilastc1+1
         do ic0=ifirstc0,ilastc0+0
            arrayc(ic0,ic1)=zero
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+1

            if1=ic1*ratio(1)
            do ic0=ifirstc0,ilastc0+0

               do ir=0,ratio(0)-1
                  if0=ic0*ratio(0)+ir
               arrayc(ic0,ic1)=arrayc(ic0,ic1)
     &                           +arrayf(if0,if1)*lengthf
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+1
         do ic0=ifirstc0,ilastc0+0
            arrayc(ic0,ic1)=arrayc(ic0,ic1)/lengthc
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 2d side-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgsidecplx2d0(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      double complex
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1)
      double precision lengthf, lengthc
      integer ic0,ic1,if0,if1,ir
c
c***********************************************************************
c
      lengthf=dxf(1)
      lengthc=dxc(1)

      do ic1=ifirstc1,ilastc1+0
         do ic0=ifirstc0,ilastc0+1
            arrayc(ic0,ic1)=cmplx(zero,zero)
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+0

         do ir=0,ratio(1)-1
            if1=ic1*ratio(1)+ir
            do ic0=ifirstc0,ilastc0+1

                  if0=ic0*ratio(0)
               arrayc(ic0,ic1)=arrayc(ic0,ic1)
     &                           +arrayf(if0,if1)*lengthf
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+0
         do ic0=ifirstc0,ilastc0+1
            arrayc(ic0,ic1)=arrayc(ic0,ic1)/lengthc
        enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgsidecplx2d1(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
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
     &  dxf(0:2-1),
     &  dxc(0:2-1)
      double complex
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1)
      double precision lengthf, lengthc
      integer ic0,ic1,if0,if1,ir
c
c***********************************************************************
c
      lengthf=dxf(0)
      lengthc=dxc(0)

      do ic1=ifirstc1,ilastc1+1
         do ic0=ifirstc0,ilastc0+0
            arrayc(ic0,ic1)=cmplx(zero,zero)
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+1

            if1=ic1*ratio(1)
            do ic0=ifirstc0,ilastc0+0

               do ir=0,ratio(0)-1
                  if0=ic0*ratio(0)+ir
               arrayc(ic0,ic1)=arrayc(ic0,ic1)
     &                           +arrayf(if0,if1)*lengthf
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+1
         do ic0=ifirstc0,ilastc0+0
            arrayc(ic0,ic1)=arrayc(ic0,ic1)/lengthc
        enddo
      enddo
c
      return
      end
c
