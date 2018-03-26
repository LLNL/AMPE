c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: FORTRAN routines for spatial refining of 3d patch data
c               on a regular Cartesian mesh.
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for 3d Cartesian refine operators
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for dimensioning 3d arrays in FORTRAN routines.
c
c
c  File:        $URL$
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: m4 include file for constant patchdata transfer routines.
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
c Constant interpolation for 3d cell-centered double data
c***********************************************************************
c
      subroutine conrefcelldoub3d(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c
      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
            do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
          enddo
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 3d cell-centered float data
c***********************************************************************
c
      subroutine conrefcellflot3d(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      real
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c
      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
            do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
          enddo
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 3d cell-centered complex data
c***********************************************************************
c
      subroutine conrefcellcplx3d(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double complex
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c
      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
            do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
          enddo
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 3d cell-centered integer data
c***********************************************************************
c
      subroutine conrefcellintg3d(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      integer
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c
      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
            do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
          enddo
        enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 3d edge-centered double data
c***********************************************************************
c
      subroutine conrefedgedoub3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1,
     &          filo2:fihi2+1)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conrefedgedoub3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2+1)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conrefedgedoub3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1,
     &          filo2:fihi2)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 3d edge-centered float data
c***********************************************************************
c
      subroutine conrefedgeflot3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      real
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1,
     &          filo2:fihi2+1)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conrefedgeflot3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      real
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2+1)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conrefedgeflot3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      real
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1,
     &          filo2:fihi2)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 3d edge-centered complex data
c***********************************************************************
c
      subroutine conrefedgecplx3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double complex
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1,
     &          filo2:fihi2+1)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conrefedgecplx3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double complex
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2+1)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conrefedgecplx3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double complex
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1,
     &          filo2:fihi2)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 3d edge-centered integer data
c***********************************************************************
c
      subroutine conrefedgeintg3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      integer
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1,
     &          filo2:fihi2+1)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conrefedgeintg3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      integer
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2+1)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conrefedgeintg3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      integer
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1,
     &          filo2:fihi2)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 3d face-centered double data
c***********************************************************************
c
      subroutine conreffacedoub3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2)
      integer ie0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c
      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
            do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ie0=(if0+1)/ratio(0)-1
         else
            ie0=if0/ratio(0)
         endif
               arrayf(if0,if1,if2)=arrayc(ie0,ic1,ic2)
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conreffacedoub3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  arrayc(cilo1:cihi1+1,
     &          cilo2:cihi2,
     &          cilo0:cihi0),
     &  arrayf(filo1:fihi1+1,
     &          filo2:fihi2,
     &          filo0:fihi0)
      integer ie1,ic2,ic0,if1,if2,if0
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
            do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ie1=(if1+1)/ratio(1)-1
         else
            ie1=if1/ratio(1)
         endif
               arrayf(if1,if2,if0)=arrayc(ie1,ic2,ic0)
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conreffacedoub3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  arrayc(cilo2:cihi2+1,
     &          cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo2:fihi2+1,
     &          filo0:fihi0,
     &          filo1:fihi1)
      integer ie2,ic0,ic1,if2,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ie2=(if2+1)/ratio(2)-1
         else
            ie2=if2/ratio(2)
         endif
               arrayf(if2,if0,if1)=arrayc(ie2,ic0,ic1)
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 3d face-centered float data
c***********************************************************************
c
      subroutine conreffaceflot3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      real
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2)
      integer ie0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c
      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
            do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ie0=(if0+1)/ratio(0)-1
         else
            ie0=if0/ratio(0)
         endif
               arrayf(if0,if1,if2)=arrayc(ie0,ic1,ic2)
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conreffaceflot3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      real
     &  arrayc(cilo1:cihi1+1,
     &          cilo2:cihi2,
     &          cilo0:cihi0),
     &  arrayf(filo1:fihi1+1,
     &          filo2:fihi2,
     &          filo0:fihi0)
      integer ie1,ic2,ic0,if1,if2,if0
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
            do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ie1=(if1+1)/ratio(1)-1
         else
            ie1=if1/ratio(1)
         endif
               arrayf(if1,if2,if0)=arrayc(ie1,ic2,ic0)
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conreffaceflot3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      real
     &  arrayc(cilo2:cihi2+1,
     &          cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo2:fihi2+1,
     &          filo0:fihi0,
     &          filo1:fihi1)
      integer ie2,ic0,ic1,if2,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ie2=(if2+1)/ratio(2)-1
         else
            ie2=if2/ratio(2)
         endif
               arrayf(if2,if0,if1)=arrayc(ie2,ic0,ic1)
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 3d face-centered complex data
c***********************************************************************
c
      subroutine conreffacecplx3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double complex
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2)
      integer ie0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c
      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
            do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ie0=(if0+1)/ratio(0)-1
         else
            ie0=if0/ratio(0)
         endif
               arrayf(if0,if1,if2)=arrayc(ie0,ic1,ic2)
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conreffacecplx3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double complex
     &  arrayc(cilo1:cihi1+1,
     &          cilo2:cihi2,
     &          cilo0:cihi0),
     &  arrayf(filo1:fihi1+1,
     &          filo2:fihi2,
     &          filo0:fihi0)
      integer ie1,ic2,ic0,if1,if2,if0
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
            do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ie1=(if1+1)/ratio(1)-1
         else
            ie1=if1/ratio(1)
         endif
               arrayf(if1,if2,if0)=arrayc(ie1,ic2,ic0)
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conreffacecplx3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double complex
     &  arrayc(cilo2:cihi2+1,
     &          cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo2:fihi2+1,
     &          filo0:fihi0,
     &          filo1:fihi1)
      integer ie2,ic0,ic1,if2,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ie2=(if2+1)/ratio(2)-1
         else
            ie2=if2/ratio(2)
         endif
               arrayf(if2,if0,if1)=arrayc(ie2,ic0,ic1)
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 3d face-centered integer data
c***********************************************************************
c
      subroutine conreffaceintg3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      integer
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2)
      integer ie0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c
      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
            do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ie0=(if0+1)/ratio(0)-1
         else
            ie0=if0/ratio(0)
         endif
               arrayf(if0,if1,if2)=arrayc(ie0,ic1,ic2)
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conreffaceintg3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      integer
     &  arrayc(cilo1:cihi1+1,
     &          cilo2:cihi2,
     &          cilo0:cihi0),
     &  arrayf(filo1:fihi1+1,
     &          filo2:fihi2,
     &          filo0:fihi0)
      integer ie1,ic2,ic0,if1,if2,if0
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
            do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ie1=(if1+1)/ratio(1)-1
         else
            ie1=if1/ratio(1)
         endif
               arrayf(if1,if2,if0)=arrayc(ie1,ic2,ic0)
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conreffaceintg3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      integer
     &  arrayc(cilo2:cihi2+1,
     &          cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo2:fihi2+1,
     &          filo0:fihi0,
     &          filo1:fihi1)
      integer ie2,ic0,ic1,if2,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ie2=(if2+1)/ratio(2)-1
         else
            ie2=if2/ratio(2)
         endif
               arrayf(if2,if0,if1)=arrayc(ie2,ic0,ic1)
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 3d outerface double data
c***********************************************************************
c
      subroutine conrefoutfacedoub3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  arrayc(cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo1:fihi1,
     &          filo2:fihi2)
      integer ic1,ic2,if1,if2
c
c***********************************************************************
c
      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
            arrayf(if1,if2)=arrayc(ic1,ic2)
         enddo
      enddo
c
      return
      end
c
      subroutine conrefoutfacedoub3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  arrayc(cilo2:cihi2,
     &          cilo0:cihi0),
     &  arrayf(filo2:fihi2,
     &          filo0:fihi0)
      integer ic2,ic0,if2,if0
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
            arrayf(if2,if0)=arrayc(ic2,ic0)
         enddo
      enddo
c
      return
      end
c
      subroutine conrefoutfacedoub3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 3d outerface float data
c***********************************************************************
c
      subroutine conrefoutfaceflot3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      real
     &  arrayc(cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo1:fihi1,
     &          filo2:fihi2)
      integer ic1,ic2,if1,if2
c
c***********************************************************************
c
      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
            arrayf(if1,if2)=arrayc(ic1,ic2)
         enddo
      enddo
c
      return
      end
c
      subroutine conrefoutfaceflot3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      real
     &  arrayc(cilo2:cihi2,
     &          cilo0:cihi0),
     &  arrayf(filo2:fihi2,
     &          filo0:fihi0)
      integer ic2,ic0,if2,if0
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
            arrayf(if2,if0)=arrayc(ic2,ic0)
         enddo
      enddo
c
      return
      end
c
      subroutine conrefoutfaceflot3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      real
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 3d outerface complex data
c***********************************************************************
c
      subroutine conrefoutfacecplx3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double complex
     &  arrayc(cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo1:fihi1,
     &          filo2:fihi2)
      integer ic1,ic2,if1,if2
c
c***********************************************************************
c
      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
            arrayf(if1,if2)=arrayc(ic1,ic2)
         enddo
      enddo
c
      return
      end
c
      subroutine conrefoutfacecplx3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double complex
     &  arrayc(cilo2:cihi2,
     &          cilo0:cihi0),
     &  arrayf(filo2:fihi2,
     &          filo0:fihi0)
      integer ic2,ic0,if2,if0
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
            arrayf(if2,if0)=arrayc(ic2,ic0)
         enddo
      enddo
c
      return
      end
c
      subroutine conrefoutfacecplx3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double complex
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 3d outerface integer data
c***********************************************************************
c
      subroutine conrefoutfaceintg3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      integer
     &  arrayc(cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo1:fihi1,
     &          filo2:fihi2)
      integer ic1,ic2,if1,if2
c
c***********************************************************************
c
      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
            arrayf(if1,if2)=arrayc(ic1,ic2)
         enddo
      enddo
c
      return
      end
c
      subroutine conrefoutfaceintg3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      integer
     &  arrayc(cilo2:cihi2,
     &          cilo0:cihi0),
     &  arrayf(filo2:fihi2,
     &          filo0:fihi0)
      integer ic2,ic0,if2,if0
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
         do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
            arrayf(if2,if0)=arrayc(ic2,ic0)
         enddo
      enddo
c
      return
      end
c
      subroutine conrefoutfaceintg3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      integer
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1)
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif
         do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 3d side-centered double data
c***********************************************************************
c
      subroutine conrefsidedoub3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conrefsidedoub3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1,
     &          filo2:fihi2)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conrefsidedoub3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double precision
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2+1)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 3d side-centered float data
c***********************************************************************
c
      subroutine conrefsideflot3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      real
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conrefsideflot3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      real
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1,
     &          filo2:fihi2)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conrefsideflot3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      real
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2+1)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 3d side-centered complex data
c***********************************************************************
c
      subroutine conrefsidecplx3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double complex
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conrefsidecplx3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double complex
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1,
     &          filo2:fihi2)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conrefsidecplx3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      double complex
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2+1)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Constant interpolation for 3d side-centered integer data
c***********************************************************************
c
      subroutine conrefsideintg3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      integer
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0+1
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conrefsideintg3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      integer
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1,
     &          filo2:fihi2)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1+1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine conrefsideintg3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:3-1)
      integer
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1),
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2+1)
      integer ic0,ic1,ic2,if0,if1,if2
c
c***********************************************************************
c


      do if2=ifirstf2,ilastf2+1
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif


         do if1=ifirstf1,ilastf1
         if (if1.lt.0) then
            ic1=(if1+1)/ratio(1)-1
         else
            ic1=if1/ratio(1)
         endif


            do if0=ifirstf0,ilastf0
         if (if0.lt.0) then
            ic0=(if0+1)/ratio(0)-1
         else
            ic0=if0/ratio(0)
         endif

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
c
