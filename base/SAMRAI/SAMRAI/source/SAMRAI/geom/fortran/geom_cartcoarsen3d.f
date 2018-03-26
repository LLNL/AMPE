c
c  File:        $URL$
c  Package:     SAMRAI geometry
c  Copyright:   (c) 1997-2012 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision$
c  Description: FORTRAN routines for spatial coarsening of 3d patch data
c               on a regular Cartesian mesh.
c
c
c  File:        $URL$
c  Package:     SAMRAI geometry
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
c
c
c
c
c***********************************************************************
c Weighted averaging for 3d cell-centered double data
c***********************************************************************
c
      subroutine cartwgtavgcelldoub3d(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double precision
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2)
      double precision dVf,dVc
      integer ic0,ic1,ic2,if0,if1,if2,ir0,ir1,ir2
c
c***********************************************************************
c
      dVf = dxf(0)*dxf(1)*dxf(2)
      dVc = dxc(0)*dxc(1)*dxc(2)

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ic0=ifirstc0,ilastc0
               arrayc(ic0,ic1,ic2)=zero
            enddo
         enddo
      enddo

      do ir2=0,ratio(2)-1
         do ir1=0,ratio(1)-1
            do ir0=0,ratio(0)-1
               do ic2=ifirstc2,ilastc2
                  if2=ic2*ratio(2)+ir2
                  do ic1=ifirstc1,ilastc1
                     if1=ic1*ratio(1)+ir1
                     do ic0=ifirstc0,ilastc0
                        if0=ic0*ratio(0)+ir0
                        arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                                      +arrayf(if0,if1,if2)*dVf
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ic0=ifirstc0,ilastc0
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/dVc
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 3d cell-centered float data
c***********************************************************************
c
      subroutine cartwgtavgcellflot3d(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      real
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2)
      double precision dVf,dVc
      integer ic0,ic1,ic2,if0,if1,if2,ir0,ir1,ir2
c
c***********************************************************************
c
      dVf = dxf(0)*dxf(1)*dxf(2)
      dVc = dxc(0)*dxc(1)*dxc(2)

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ic0=ifirstc0,ilastc0
               arrayc(ic0,ic1,ic2)=zero
            enddo
         enddo
      enddo

      do ir2=0,ratio(2)-1
         do ir1=0,ratio(1)-1
            do ir0=0,ratio(0)-1
               do ic2=ifirstc2,ilastc2
                  if2=ic2*ratio(2)+ir2
                  do ic1=ifirstc1,ilastc1
                     if1=ic1*ratio(1)+ir1
                     do ic0=ifirstc0,ilastc0
                        if0=ic0*ratio(0)+ir0
                        arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                                      +arrayf(if0,if1,if2)*dVf
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ic0=ifirstc0,ilastc0
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/dVc
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 3d cell-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgcellcplx3d(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double complex
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2)
      double precision dVf,dVc
      integer ic0,ic1,ic2,if0,if1,if2,ir0,ir1,ir2
c
c***********************************************************************
c
      dVf = dxf(0)*dxf(1)*dxf(2)
      dVc = dxc(0)*dxc(1)*dxc(2)

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ic0=ifirstc0,ilastc0
               arrayc(ic0,ic1,ic2)=cmplx(zero,zero)
            enddo
         enddo
      enddo

      do ir2=0,ratio(2)-1
         do ir1=0,ratio(1)-1
            do ir0=0,ratio(0)-1
               do ic2=ifirstc2,ilastc2
                  if2=ic2*ratio(2)+ir2
                  do ic1=ifirstc1,ilastc1
                     if1=ic1*ratio(1)+ir1
                     do ic0=ifirstc0,ilastc0
                        if0=ic0*ratio(0)+ir0
                        arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                                      +arrayf(if0,if1,if2)*dVf
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ic0=ifirstc0,ilastc0
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/dVc
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 3d edge-centered double data
c***********************************************************************
c
      subroutine cartwgtavgedgedoub3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double precision
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1,
     &          filo2:fihi2+1),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1)
      double precision lengthf,lengthc
      integer ic0,ic1,ic2,if0,if1,if2,ir
c
c***********************************************************************
c
      lengthf=dxf(0)
      lengthc=dxc(0)


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0

               arrayc(ic0,ic1,ic2)=zero
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1
            if2=ic2*ratio(2)


         do ic1=ifirstc1,ilastc1+1
               if1=ic1*ratio(1)


            do ic0=ifirstc0,ilastc0
               do ir=0,ratio(0)-1
                  if0=ic0*ratio(0)+ir
                  arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                             +arrayf(if0,if1,if2)*lengthf
               enddo
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/lengthc
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgedgedoub3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double precision
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2+1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1)
      double precision lengthf,lengthc
      integer ic0,ic1,ic2,if0,if1,if2,ir
c
c***********************************************************************
c
      lengthf=dxf(1)
      lengthc=dxc(1)


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0+1

               arrayc(ic0,ic1,ic2)=zero
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1
            if2=ic2*ratio(2)


         do ic1=ifirstc1,ilastc1
            do ir=0,ratio(1)-1
               if1=ic1*ratio(1)+ir


            do ic0=ifirstc0,ilastc0+1
                  if0=ic0*ratio(0)
                  arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                             +arrayf(if0,if1,if2)*lengthf
               enddo
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0+1
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/lengthc
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgedgedoub3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double precision
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1,
     &          filo2:fihi2),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2)
      double precision lengthf,lengthc
      integer ic0,ic1,ic2,if0,if1,if2,ir
c
c***********************************************************************
c
      lengthf=dxf(2)
      lengthc=dxc(2)


      do ic2=ifirstc2,ilastc2 

         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0+1

               arrayc(ic0,ic1,ic2)=zero
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2 
         do ir=0,ratio(2)-1
            if2=ic2*ratio(2)+ir


         do ic1=ifirstc1,ilastc1+1
               if1=ic1*ratio(1)


            do ic0=ifirstc0,ilastc0+1
                  if0=ic0*ratio(0)
                  arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                             +arrayf(if0,if1,if2)*lengthf
               enddo
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2 

         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0+1
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/lengthc
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 3d edge-centered float data
c***********************************************************************
c
      subroutine cartwgtavgedgeflot3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      real
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1,
     &          filo2:fihi2+1),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1)
      double precision lengthf,lengthc
      integer ic0,ic1,ic2,if0,if1,if2,ir
c
c***********************************************************************
c
      lengthf=dxf(0)
      lengthc=dxc(0)


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0

               arrayc(ic0,ic1,ic2)=zero
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1
            if2=ic2*ratio(2)


         do ic1=ifirstc1,ilastc1+1
               if1=ic1*ratio(1)


            do ic0=ifirstc0,ilastc0
               do ir=0,ratio(0)-1
                  if0=ic0*ratio(0)+ir
                  arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                             +arrayf(if0,if1,if2)*lengthf
               enddo
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/lengthc
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgedgeflot3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      real
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2+1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1)
      double precision lengthf,lengthc
      integer ic0,ic1,ic2,if0,if1,if2,ir
c
c***********************************************************************
c
      lengthf=dxf(1)
      lengthc=dxc(1)


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0+1

               arrayc(ic0,ic1,ic2)=zero
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1
            if2=ic2*ratio(2)


         do ic1=ifirstc1,ilastc1
            do ir=0,ratio(1)-1
               if1=ic1*ratio(1)+ir


            do ic0=ifirstc0,ilastc0+1
                  if0=ic0*ratio(0)
                  arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                             +arrayf(if0,if1,if2)*lengthf
               enddo
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0+1
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/lengthc
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgedgeflot3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      real
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1,
     &          filo2:fihi2),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2)
      double precision lengthf,lengthc
      integer ic0,ic1,ic2,if0,if1,if2,ir
c
c***********************************************************************
c
      lengthf=dxf(2)
      lengthc=dxc(2)


      do ic2=ifirstc2,ilastc2 

         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0+1

               arrayc(ic0,ic1,ic2)=zero
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2 
         do ir=0,ratio(2)-1
            if2=ic2*ratio(2)+ir


         do ic1=ifirstc1,ilastc1+1
               if1=ic1*ratio(1)


            do ic0=ifirstc0,ilastc0+1
                  if0=ic0*ratio(0)
                  arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                             +arrayf(if0,if1,if2)*lengthf
               enddo
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2 

         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0+1
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/lengthc
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 3d edge-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgedgecplx3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double complex
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1,
     &          filo2:fihi2+1),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2+1)
      double precision lengthf,lengthc
      integer ic0,ic1,ic2,if0,if1,if2,ir
c
c***********************************************************************
c
      lengthf=dxf(0)
      lengthc=dxc(0)


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0

               arrayc(ic0,ic1,ic2)=cmplx(zero,zero)
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1
            if2=ic2*ratio(2)


         do ic1=ifirstc1,ilastc1+1
               if1=ic1*ratio(1)


            do ic0=ifirstc0,ilastc0
               do ir=0,ratio(0)-1
                  if0=ic0*ratio(0)+ir
                  arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                             +arrayf(if0,if1,if2)*lengthf
               enddo
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/lengthc
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgedgecplx3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double complex
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2+1),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1)
      double precision lengthf,lengthc
      integer ic0,ic1,ic2,if0,if1,if2,ir
c
c***********************************************************************
c
      lengthf=dxf(1)
      lengthc=dxc(1)


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0+1

               arrayc(ic0,ic1,ic2)=cmplx(zero,zero)
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1
            if2=ic2*ratio(2)


         do ic1=ifirstc1,ilastc1
            do ir=0,ratio(1)-1
               if1=ic1*ratio(1)+ir


            do ic0=ifirstc0,ilastc0+1
                  if0=ic0*ratio(0)
                  arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                             +arrayf(if0,if1,if2)*lengthf
               enddo
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0+1
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/lengthc
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgedgecplx3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double complex
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1+1,
     &          filo2:fihi2),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2)
      double precision lengthf,lengthc
      integer ic0,ic1,ic2,if0,if1,if2,ir
c
c***********************************************************************
c
      lengthf=dxf(2)
      lengthc=dxc(2)


      do ic2=ifirstc2,ilastc2 

         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0+1

               arrayc(ic0,ic1,ic2)=cmplx(zero,zero)
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2 
         do ir=0,ratio(2)-1
            if2=ic2*ratio(2)+ir


         do ic1=ifirstc1,ilastc1+1
               if1=ic1*ratio(1)


            do ic0=ifirstc0,ilastc0+1
                  if0=ic0*ratio(0)
                  arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                             +arrayf(if0,if1,if2)*lengthf
               enddo
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2 

         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0+1
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/lengthc
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 3d face-centered double data
c***********************************************************************
c
      subroutine cartwgtavgfacedoub3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double precision
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2)
      double precision areaf,areac
      integer ie0,ic1,ic2,if0,if1,if2,ir1,ir2
c
c***********************************************************************
c
      areaf=dxf(1)*dxf(2)
      areac=dxc(1)*dxc(2)

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ie0=ifirstc0,ilastc0+1
               arrayc(ie0,ic1,ic2)=zero
            enddo
         enddo
      enddo

      do ir2=0,ratio(2)-1
         do ir1=0,ratio(1)-1
            do ic2=ifirstc2,ilastc2
               if2=ic2*ratio(2)+ir2
               do ic1=ifirstc1,ilastc1
                  if1=ic1*ratio(1)+ir1
                  do ie0=ifirstc0,ilastc0+1
                     if0=ie0*ratio(0)
                     arrayc(ie0,ic1,ic2)=arrayc(ie0,ic1,ic2)
     &                                +arrayf(if0,if1,if2)*areaf
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ie0=ifirstc0,ilastc0+1
               arrayc(ie0,ic1,ic2)=arrayc(ie0,ic1,ic2)/areac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgfacedoub3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double precision
     &  arrayf(filo1:fihi1+1,
     &          filo2:fihi2,
     &          filo0:fihi0),
     &  arrayc(cilo1:cihi1+1,
     &          cilo2:cihi2,
     &          cilo0:cihi0)
      double precision areaf,areac
      integer ie1,ic2,ic0,if1,if2,if0,ir2,ir0
c
c***********************************************************************
c
      areaf=dxf(2)*dxf(0)
      areac=dxc(2)*dxc(0)

      do ic0=ifirstc0,ilastc0
         do ic2=ifirstc2,ilastc2
            do ie1=ifirstc1,ilastc1+1
               arrayc(ie1,ic2,ic0)=zero
            enddo
         enddo
      enddo

      do ir0=0,ratio(0)-1
         do ir2=0,ratio(2)-1
            do ic0=ifirstc0,ilastc0
               if0=ic0*ratio(0)+ir0
               do ic2=ifirstc2,ilastc2
                  if2=ic2*ratio(2)+ir2
                  do ie1=ifirstc1,ilastc1+1
                     if1=ie1*ratio(1)
                     arrayc(ie1,ic2,ic0)=arrayc(ie1,ic2,ic0)
     &                                +arrayf(if1,if2,if0)*areaf
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         do ic2=ifirstc2,ilastc2
            do ie1=ifirstc1,ilastc1+1
               arrayc(ie1,ic2,ic0)=arrayc(ie1,ic2,ic0)/areac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgfacedoub3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double precision
     &  arrayf(filo2:fihi2+1,
     &          filo0:fihi0,
     &          filo1:fihi1),
     &  arrayc(cilo2:cihi2+1,
     &          cilo0:cihi0,
     &          cilo1:cihi1)
      double precision areaf,areac
      integer ie2,ic0,ic1,if2,if0,if1,ir0,ir1
c
c***********************************************************************
c
      areaf=dxf(0)*dxf(1)
      areac=dxc(0)*dxc(1)

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            do ie2=ifirstc2,ilastc2+1
               arrayc(ie2,ic0,ic1)=zero
            enddo
         enddo
      enddo

      do ir1=0,ratio(1)-1
         do ir0=0,ratio(0)-1
            do ic1=ifirstc1,ilastc1
               if1=ic1*ratio(1)+ir1
               do ic0=ifirstc0,ilastc0
                  if0=ic0*ratio(0)+ir0
                  do ie2=ifirstc2,ilastc2+1
                     if2=ie2*ratio(2)
                     arrayc(ie2,ic0,ic1)=arrayc(ie2,ic0,ic1)
     &                                +arrayf(if2,if0,if1)*areaf
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            do ie2=ifirstc2,ilastc2+1
               arrayc(ie2,ic0,ic1)=arrayc(ie2,ic0,ic1)/areac
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 3d face-centered float data
c***********************************************************************
c
      subroutine cartwgtavgfaceflot3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      real
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2)
      double precision areaf,areac
      integer ie0,ic1,ic2,if0,if1,if2,ir1,ir2
c
c***********************************************************************
c
      areaf=dxf(1)*dxf(2)
      areac=dxc(1)*dxc(2)

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ie0=ifirstc0,ilastc0+1
               arrayc(ie0,ic1,ic2)=zero
            enddo
         enddo
      enddo

      do ir2=0,ratio(2)-1
         do ir1=0,ratio(1)-1
            do ic2=ifirstc2,ilastc2
               if2=ic2*ratio(2)+ir2
               do ic1=ifirstc1,ilastc1
                  if1=ic1*ratio(1)+ir1
                  do ie0=ifirstc0,ilastc0+1
                     if0=ie0*ratio(0)
                     arrayc(ie0,ic1,ic2)=arrayc(ie0,ic1,ic2)
     &                                +arrayf(if0,if1,if2)*areaf
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ie0=ifirstc0,ilastc0+1
               arrayc(ie0,ic1,ic2)=arrayc(ie0,ic1,ic2)/areac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgfaceflot3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      real
     &  arrayf(filo1:fihi1+1,
     &          filo2:fihi2,
     &          filo0:fihi0),
     &  arrayc(cilo1:cihi1+1,
     &          cilo2:cihi2,
     &          cilo0:cihi0)
      double precision areaf,areac
      integer ie1,ic2,ic0,if1,if2,if0,ir2,ir0
c
c***********************************************************************
c
      areaf=dxf(2)*dxf(0)
      areac=dxc(2)*dxc(0)

      do ic0=ifirstc0,ilastc0
         do ic2=ifirstc2,ilastc2
            do ie1=ifirstc1,ilastc1+1
               arrayc(ie1,ic2,ic0)=zero
            enddo
         enddo
      enddo

      do ir0=0,ratio(0)-1
         do ir2=0,ratio(2)-1
            do ic0=ifirstc0,ilastc0
               if0=ic0*ratio(0)+ir0
               do ic2=ifirstc2,ilastc2
                  if2=ic2*ratio(2)+ir2
                  do ie1=ifirstc1,ilastc1+1
                     if1=ie1*ratio(1)
                     arrayc(ie1,ic2,ic0)=arrayc(ie1,ic2,ic0)
     &                                +arrayf(if1,if2,if0)*areaf
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         do ic2=ifirstc2,ilastc2
            do ie1=ifirstc1,ilastc1+1
               arrayc(ie1,ic2,ic0)=arrayc(ie1,ic2,ic0)/areac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgfaceflot3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      real
     &  arrayf(filo2:fihi2+1,
     &          filo0:fihi0,
     &          filo1:fihi1),
     &  arrayc(cilo2:cihi2+1,
     &          cilo0:cihi0,
     &          cilo1:cihi1)
      double precision areaf,areac
      integer ie2,ic0,ic1,if2,if0,if1,ir0,ir1
c
c***********************************************************************
c
      areaf=dxf(0)*dxf(1)
      areac=dxc(0)*dxc(1)

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            do ie2=ifirstc2,ilastc2+1
               arrayc(ie2,ic0,ic1)=zero
            enddo
         enddo
      enddo

      do ir1=0,ratio(1)-1
         do ir0=0,ratio(0)-1
            do ic1=ifirstc1,ilastc1
               if1=ic1*ratio(1)+ir1
               do ic0=ifirstc0,ilastc0
                  if0=ic0*ratio(0)+ir0
                  do ie2=ifirstc2,ilastc2+1
                     if2=ie2*ratio(2)
                     arrayc(ie2,ic0,ic1)=arrayc(ie2,ic0,ic1)
     &                                +arrayf(if2,if0,if1)*areaf
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            do ie2=ifirstc2,ilastc2+1
               arrayc(ie2,ic0,ic1)=arrayc(ie2,ic0,ic1)/areac
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 3d face-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgfacecplx3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double complex
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2)
      double precision areaf,areac
      integer ie0,ic1,ic2,if0,if1,if2,ir1,ir2
c
c***********************************************************************
c
      areaf=dxf(1)*dxf(2)
      areac=dxc(1)*dxc(2)

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ie0=ifirstc0,ilastc0+1
               arrayc(ie0,ic1,ic2)=cmplx(zero,zero)
            enddo
         enddo
      enddo

      do ir2=0,ratio(2)-1
         do ir1=0,ratio(1)-1
            do ic2=ifirstc2,ilastc2
               if2=ic2*ratio(2)+ir2
               do ic1=ifirstc1,ilastc1
                  if1=ic1*ratio(1)+ir1
                  do ie0=ifirstc0,ilastc0+1
                     if0=ie0*ratio(0)
                     arrayc(ie0,ic1,ic2)=arrayc(ie0,ic1,ic2)
     &                                +arrayf(if0,if1,if2)*areaf
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ie0=ifirstc0,ilastc0+1
               arrayc(ie0,ic1,ic2)=arrayc(ie0,ic1,ic2)/areac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgfacecplx3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double complex
     &  arrayf(filo1:fihi1+1,
     &          filo2:fihi2,
     &          filo0:fihi0),
     &  arrayc(cilo1:cihi1+1,
     &          cilo2:cihi2,
     &          cilo0:cihi0)
      double precision areaf,areac
      integer ie1,ic2,ic0,if1,if2,if0,ir2,ir0
c
c***********************************************************************
c
      areaf=dxf(2)*dxf(0)
      areac=dxc(2)*dxc(0)

      do ic0=ifirstc0,ilastc0
         do ic2=ifirstc2,ilastc2
            do ie1=ifirstc1,ilastc1+1
               arrayc(ie1,ic2,ic0)=cmplx(zero,zero)
            enddo
         enddo
      enddo

      do ir0=0,ratio(0)-1
         do ir2=0,ratio(2)-1
            do ic0=ifirstc0,ilastc0
               if0=ic0*ratio(0)+ir0
               do ic2=ifirstc2,ilastc2
                  if2=ic2*ratio(2)+ir2
                  do ie1=ifirstc1,ilastc1+1
                     if1=ie1*ratio(1)
                     arrayc(ie1,ic2,ic0)=arrayc(ie1,ic2,ic0)
     &                                +arrayf(if1,if2,if0)*areaf
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         do ic2=ifirstc2,ilastc2
            do ie1=ifirstc1,ilastc1+1
               arrayc(ie1,ic2,ic0)=arrayc(ie1,ic2,ic0)/areac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgfacecplx3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double complex
     &  arrayf(filo2:fihi2+1,
     &          filo0:fihi0,
     &          filo1:fihi1),
     &  arrayc(cilo2:cihi2+1,
     &          cilo0:cihi0,
     &          cilo1:cihi1)
      double precision areaf,areac
      integer ie2,ic0,ic1,if2,if0,if1,ir0,ir1
c
c***********************************************************************
c
      areaf=dxf(0)*dxf(1)
      areac=dxc(0)*dxc(1)

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            do ie2=ifirstc2,ilastc2+1
               arrayc(ie2,ic0,ic1)=cmplx(zero,zero)
            enddo
         enddo
      enddo

      do ir1=0,ratio(1)-1
         do ir0=0,ratio(0)-1
            do ic1=ifirstc1,ilastc1
               if1=ic1*ratio(1)+ir1
               do ic0=ifirstc0,ilastc0
                  if0=ic0*ratio(0)+ir0
                  do ie2=ifirstc2,ilastc2+1
                     if2=ie2*ratio(2)
                     arrayc(ie2,ic0,ic1)=arrayc(ie2,ic0,ic1)
     &                                +arrayf(if2,if0,if1)*areaf
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            do ie2=ifirstc2,ilastc2+1
               arrayc(ie2,ic0,ic1)=arrayc(ie2,ic0,ic1)/areac
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 3d outerface double data
c***********************************************************************
c
      subroutine cartwgtavgoutfacedoub3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double precision
     &  arrayf(filo1:fihi1,
     &          filo2:fihi2),
     &  arrayc(cilo1:cihi1,
     &          cilo2:cihi2)
      double precision areaf,areac
      integer ic1,ic2,if1,if2,ir1,ir2
c
c***********************************************************************
c
      areaf=dxf(1)*dxf(2)
      areac=dxc(1)*dxc(2)

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            arrayc(ic1,ic2)=zero
         enddo
      enddo

      do ir2=0,ratio(2)-1
         do ir1=0,ratio(1)-1
            do ic2=ifirstc2,ilastc2
               if2=ic2*ratio(2)+ir2
               do ic1=ifirstc1,ilastc1
                  if1=ic1*ratio(1)+ir1
                  arrayc(ic1,ic2)=arrayc(ic1,ic2)
     &                              +arrayf(if1,if2)*areaf
               enddo
            enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            arrayc(ic1,ic2)=arrayc(ic1,ic2)/areac
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgoutfacedoub3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double precision
     &  arrayf(filo2:fihi2,
     &          filo0:fihi0),
     &  arrayc(cilo2:cihi2,
     &          cilo0:cihi0)
      double precision areaf,areac
      integer ic2,ic0,if2,if0,ir2,ir0
c
c***********************************************************************
c
      areaf=dxf(2)*dxf(0)
      areac=dxc(2)*dxc(0)

      do ic0=ifirstc0,ilastc0
         do ic2=ifirstc2,ilastc2
            arrayc(ic2,ic0)=zero
         enddo
      enddo

      do ir0=0,ratio(0)-1
         do ir2=0,ratio(2)-1
            do ic0=ifirstc0,ilastc0
               if0=ic0*ratio(0)+ir0
               do ic2=ifirstc2,ilastc2
                  if2=ic2*ratio(2)+ir2
                  arrayc(ic2,ic0)=arrayc(ic2,ic0)
     &                              +arrayf(if2,if0)*areaf
               enddo
            enddo
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         do ic2=ifirstc2,ilastc2
            arrayc(ic2,ic0)=arrayc(ic2,ic0)/areac
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgoutfacedoub3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double precision
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1)
      double precision areaf,areac
      integer ic0,ic1,if0,if1,ir0,ir1
c
c***********************************************************************
c
      areaf=dxf(0)*dxf(1)
      areac=dxc(0)*dxc(1)

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
     &                              +arrayf(if0,if1)*areaf
               enddo
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            arrayc(ic0,ic1)=arrayc(ic0,ic1)/areac
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 3d outerface float data
c***********************************************************************
c
      subroutine cartwgtavgoutfaceflot3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      real
     &  arrayf(filo1:fihi1,
     &          filo2:fihi2),
     &  arrayc(cilo1:cihi1,
     &          cilo2:cihi2)
      double precision areaf,areac
      integer ic1,ic2,if1,if2,ir1,ir2
c
c***********************************************************************
c
      areaf=dxf(1)*dxf(2)
      areac=dxc(1)*dxc(2)

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            arrayc(ic1,ic2)=zero
         enddo
      enddo

      do ir2=0,ratio(2)-1
         do ir1=0,ratio(1)-1
            do ic2=ifirstc2,ilastc2
               if2=ic2*ratio(2)+ir2
               do ic1=ifirstc1,ilastc1
                  if1=ic1*ratio(1)+ir1
                  arrayc(ic1,ic2)=arrayc(ic1,ic2)
     &                              +arrayf(if1,if2)*areaf
               enddo
            enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            arrayc(ic1,ic2)=arrayc(ic1,ic2)/areac
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgoutfaceflot3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      real
     &  arrayf(filo2:fihi2,
     &          filo0:fihi0),
     &  arrayc(cilo2:cihi2,
     &          cilo0:cihi0)
      double precision areaf,areac
      integer ic2,ic0,if2,if0,ir2,ir0
c
c***********************************************************************
c
      areaf=dxf(2)*dxf(0)
      areac=dxc(2)*dxc(0)

      do ic0=ifirstc0,ilastc0
         do ic2=ifirstc2,ilastc2
            arrayc(ic2,ic0)=zero
         enddo
      enddo

      do ir0=0,ratio(0)-1
         do ir2=0,ratio(2)-1
            do ic0=ifirstc0,ilastc0
               if0=ic0*ratio(0)+ir0
               do ic2=ifirstc2,ilastc2
                  if2=ic2*ratio(2)+ir2
                  arrayc(ic2,ic0)=arrayc(ic2,ic0)
     &                              +arrayf(if2,if0)*areaf
               enddo
            enddo
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         do ic2=ifirstc2,ilastc2
            arrayc(ic2,ic0)=arrayc(ic2,ic0)/areac
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgoutfaceflot3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      real
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1)
      double precision areaf,areac
      integer ic0,ic1,if0,if1,ir0,ir1
c
c***********************************************************************
c
      areaf=dxf(0)*dxf(1)
      areac=dxc(0)*dxc(1)

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
     &                              +arrayf(if0,if1)*areaf
               enddo
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            arrayc(ic0,ic1)=arrayc(ic0,ic1)/areac
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 3d outerface complex data
c***********************************************************************
c
      subroutine cartwgtavgoutfacecplx3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double complex
     &  arrayf(filo1:fihi1,
     &          filo2:fihi2),
     &  arrayc(cilo1:cihi1,
     &          cilo2:cihi2)
      double precision areaf,areac
      integer ic1,ic2,if1,if2,ir1,ir2
c
c***********************************************************************
c
      areaf=dxf(1)*dxf(2)
      areac=dxc(1)*dxc(2)

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            arrayc(ic1,ic2)=cmplx(zero,zero)
         enddo
      enddo

      do ir2=0,ratio(2)-1
         do ir1=0,ratio(1)-1
            do ic2=ifirstc2,ilastc2
               if2=ic2*ratio(2)+ir2
               do ic1=ifirstc1,ilastc1
                  if1=ic1*ratio(1)+ir1
                  arrayc(ic1,ic2)=arrayc(ic1,ic2)
     &                              +arrayf(if1,if2)*areaf
               enddo
            enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            arrayc(ic1,ic2)=arrayc(ic1,ic2)/areac
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgoutfacecplx3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double complex
     &  arrayf(filo2:fihi2,
     &          filo0:fihi0),
     &  arrayc(cilo2:cihi2,
     &          cilo0:cihi0)
      double precision areaf,areac
      integer ic2,ic0,if2,if0,ir2,ir0
c
c***********************************************************************
c
      areaf=dxf(2)*dxf(0)
      areac=dxc(2)*dxc(0)

      do ic0=ifirstc0,ilastc0
         do ic2=ifirstc2,ilastc2
            arrayc(ic2,ic0)=cmplx(zero,zero)
         enddo
      enddo

      do ir0=0,ratio(0)-1
         do ir2=0,ratio(2)-1
            do ic0=ifirstc0,ilastc0
               if0=ic0*ratio(0)+ir0
               do ic2=ifirstc2,ilastc2
                  if2=ic2*ratio(2)+ir2
                  arrayc(ic2,ic0)=arrayc(ic2,ic0)
     &                              +arrayf(if2,if0)*areaf
               enddo
            enddo
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         do ic2=ifirstc2,ilastc2
            arrayc(ic2,ic0)=arrayc(ic2,ic0)/areac
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgoutfacecplx3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double complex
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1)
      double precision areaf,areac
      integer ic0,ic1,if0,if1,ir0,ir1
c
c***********************************************************************
c
      areaf=dxf(0)*dxf(1)
      areac=dxc(0)*dxc(1)

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
     &                              +arrayf(if0,if1)*areaf
               enddo
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            arrayc(ic0,ic1)=arrayc(ic0,ic1)/areac
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 3d outerside double data
c***********************************************************************
c
      subroutine cartwgtavgoutsidedoub3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double precision
     &  arrayf(filo1:fihi1,
     &          filo2:fihi2),
     &  arrayc(cilo1:cihi1,
     &          cilo2:cihi2)
      double precision areaf,areac
      integer ic_outer,ic_inner,if_outer,if_inner,
     &  ir_outer,ir_inner
c
c***********************************************************************
c
      areaf=dxf(1)*dxf(2)
      areac=dxc(1)*dxc(2)


      do ic_outer=ifirstc2,ilastc2
         do ic_inner=ifirstc1,ilastc1

            arrayc(ic_inner,ic_outer)=zero
         enddo
      enddo


      do ic_outer=ifirstc2,ilastc2
         do ir_outer=0,ratio(2)-1
            if_outer=ic_outer*ratio(2)+ir_outer
            do ic_inner=ifirstc1,ilastc1
               do ir_inner=0,ratio(1)-1
               if_inner=ic_inner*ratio(1)+ir_inner
 
                  arrayc(ic_inner,ic_outer)=
     &               arrayc(ic_inner,ic_outer)
     &               +arrayf(if_inner,if_outer)*areaf
               enddo
            enddo
         enddo
      enddo


      do ic_outer=ifirstc2,ilastc2
         do ic_inner=ifirstc1,ilastc1
            arrayc(ic_inner,ic_outer)=
     &         arrayc(ic_inner,ic_outer)/areac
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgoutsidedoub3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double precision
     &  arrayf(filo0:fihi0,
     &          filo2:fihi2),
     &  arrayc(cilo0:cihi0,
     &          cilo2:cihi2)
      double precision areaf,areac
      integer ic_outer,ic_inner,if_outer,if_inner,
     &  ir_outer,ir_inner
c
c***********************************************************************
c
      areaf=dxf(0)*dxf(2)
      areac=dxc(0)*dxc(2)


      do ic_outer=ifirstc2,ilastc2
         do ic_inner=ifirstc0,ilastc0

            arrayc(ic_inner,ic_outer)=zero
         enddo
      enddo


      do ic_outer=ifirstc2,ilastc2
         do ir_outer=0,ratio(2)-1
            if_outer=ic_outer*ratio(2)+ir_outer
            do ic_inner=ifirstc0,ilastc0
               do ir_inner=0,ratio(0)-1
               if_inner=ic_inner*ratio(0)+ir_inner
 
                  arrayc(ic_inner,ic_outer)=
     &               arrayc(ic_inner,ic_outer)
     &               +arrayf(if_inner,if_outer)*areaf
               enddo
            enddo
         enddo
      enddo


      do ic_outer=ifirstc2,ilastc2
         do ic_inner=ifirstc0,ilastc0
            arrayc(ic_inner,ic_outer)=
     &         arrayc(ic_inner,ic_outer)/areac
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgoutsidedoub3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double precision
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1)
      double precision areaf,areac
      integer ic_outer,ic_inner,if_outer,if_inner,
     &  ir_outer,ir_inner
c
c***********************************************************************
c
      areaf=dxf(0)*dxf(1)
      areac=dxc(0)*dxc(1)


      do ic_outer=ifirstc1,ilastc1
         do ic_inner=ifirstc0,ilastc0

            arrayc(ic_inner,ic_outer)=zero
         enddo
      enddo


      do ic_outer=ifirstc1,ilastc1
         do ir_outer=0,ratio(1)-1
            if_outer=ic_outer*ratio(1)+ir_outer
            do ic_inner=ifirstc0,ilastc0
               do ir_inner=0,ratio(0)-1
               if_inner=ic_inner*ratio(0)+ir_inner
 
                  arrayc(ic_inner,ic_outer)=
     &               arrayc(ic_inner,ic_outer)
     &               +arrayf(if_inner,if_outer)*areaf
               enddo
            enddo
         enddo
      enddo


      do ic_outer=ifirstc1,ilastc1
         do ic_inner=ifirstc0,ilastc0
            arrayc(ic_inner,ic_outer)=
     &         arrayc(ic_inner,ic_outer)/areac
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 3d side-centered double data
c***********************************************************************
c
      subroutine cartwgtavgsidedoub3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double precision
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2)
      double precision areaf,areac
      integer ic0,ic1,ic2,if0,if1,if2,

     &  ir1,ir2
c
c***********************************************************************
c
      areaf=dxf(1)*dxf(2)
      areac=dxc(1)*dxc(2)


      do ic2=ifirstc2,ilastc2 

         do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0+1

               arrayc(ic0,ic1,ic2)=zero
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2 
         do ir2=0,ratio(2)-1
               if2=ic2*ratio(2)+ir2

         do ic1=ifirstc1,ilastc1
            do ir1=0,ratio(1)-1
                  if1=ic1*ratio(1)+ir1

            do ic0=ifirstc0,ilastc0+1
                  if0=ic0*ratio(0)
                     arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                                +arrayf(if0,if1,if2)*areaf
                  enddo
               enddo
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2 

         do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0+1
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/areac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgsidedoub3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double precision
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1,
     &          filo2:fihi2),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2)
      double precision areaf,areac
      integer ic0,ic1,ic2,if0,if1,if2,

     &  ir0,ir2
c
c***********************************************************************
c
      areaf=dxf(2)*dxf(0)
      areac=dxc(2)*dxc(0)


      do ic2=ifirstc2,ilastc2 

         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0

               arrayc(ic0,ic1,ic2)=zero
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2 
         do ir2=0,ratio(2)-1
               if2=ic2*ratio(2)+ir2

         do ic1=ifirstc1,ilastc1+1
                  if1=ic1*ratio(1)

            do ic0=ifirstc0,ilastc0
               do ir0=0,ratio(0)-1
                  if0=ic0*ratio(0)+ir0
                     arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                                +arrayf(if0,if1,if2)*areaf
                  enddo
               enddo
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2 

         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/areac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgsidedoub3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double precision
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2+1),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1)
      double precision areaf,areac
      integer ic0,ic1,ic2,if0,if1,if2,

     &  ir0,ir1
c
c***********************************************************************
c
      areaf=dxf(0)*dxf(1)
      areac=dxc(0)*dxc(1)


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0

               arrayc(ic0,ic1,ic2)=zero
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1
               if2=ic2*ratio(2)

         do ic1=ifirstc1,ilastc1
            do ir1=0,ratio(1)-1
                  if1=ic1*ratio(1)+ir1

            do ic0=ifirstc0,ilastc0
               do ir0=0,ratio(0)-1
                  if0=ic0*ratio(0)+ir0
                     arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                                +arrayf(if0,if1,if2)*areaf
                  enddo
               enddo
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/areac
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 3d side-centered float data
c***********************************************************************
c
      subroutine cartwgtavgsideflot3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      real
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2)
      double precision areaf,areac
      integer ic0,ic1,ic2,if0,if1,if2,

     &  ir1,ir2
c
c***********************************************************************
c
      areaf=dxf(1)*dxf(2)
      areac=dxc(1)*dxc(2)


      do ic2=ifirstc2,ilastc2 

         do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0+1

               arrayc(ic0,ic1,ic2)=zero
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2 
         do ir2=0,ratio(2)-1
               if2=ic2*ratio(2)+ir2

         do ic1=ifirstc1,ilastc1
            do ir1=0,ratio(1)-1
                  if1=ic1*ratio(1)+ir1

            do ic0=ifirstc0,ilastc0+1
                  if0=ic0*ratio(0)
                     arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                                +arrayf(if0,if1,if2)*areaf
                  enddo
               enddo
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2 

         do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0+1
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/areac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgsideflot3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      real
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1,
     &          filo2:fihi2),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2)
      double precision areaf,areac
      integer ic0,ic1,ic2,if0,if1,if2,

     &  ir0,ir2
c
c***********************************************************************
c
      areaf=dxf(2)*dxf(0)
      areac=dxc(2)*dxc(0)


      do ic2=ifirstc2,ilastc2 

         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0

               arrayc(ic0,ic1,ic2)=zero
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2 
         do ir2=0,ratio(2)-1
               if2=ic2*ratio(2)+ir2

         do ic1=ifirstc1,ilastc1+1
                  if1=ic1*ratio(1)

            do ic0=ifirstc0,ilastc0
               do ir0=0,ratio(0)-1
                  if0=ic0*ratio(0)+ir0
                     arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                                +arrayf(if0,if1,if2)*areaf
                  enddo
               enddo
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2 

         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/areac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgsideflot3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      real
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2+1),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1)
      double precision areaf,areac
      integer ic0,ic1,ic2,if0,if1,if2,

     &  ir0,ir1
c
c***********************************************************************
c
      areaf=dxf(0)*dxf(1)
      areac=dxc(0)*dxc(1)


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0

               arrayc(ic0,ic1,ic2)=zero
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1
               if2=ic2*ratio(2)

         do ic1=ifirstc1,ilastc1
            do ir1=0,ratio(1)-1
                  if1=ic1*ratio(1)+ir1

            do ic0=ifirstc0,ilastc0
               do ir0=0,ratio(0)-1
                  if0=ic0*ratio(0)+ir0
                     arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                                +arrayf(if0,if1,if2)*areaf
                  enddo
               enddo
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/areac
            enddo
         enddo
      enddo
c
      return
      end
c
c***********************************************************************
c Weighted averaging for 3d side-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgsidecplx3d0(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double complex
     &  arrayf(filo0:fihi0+1,
     &          filo1:fihi1,
     &          filo2:fihi2),
     &  arrayc(cilo0:cihi0+1,
     &          cilo1:cihi1,
     &          cilo2:cihi2)
      double precision areaf,areac
      integer ic0,ic1,ic2,if0,if1,if2,

     &  ir1,ir2
c
c***********************************************************************
c
      areaf=dxf(1)*dxf(2)
      areac=dxc(1)*dxc(2)


      do ic2=ifirstc2,ilastc2 

         do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0+1

               arrayc(ic0,ic1,ic2)=cmplx(zero,zero)
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2 
         do ir2=0,ratio(2)-1
               if2=ic2*ratio(2)+ir2

         do ic1=ifirstc1,ilastc1
            do ir1=0,ratio(1)-1
                  if1=ic1*ratio(1)+ir1

            do ic0=ifirstc0,ilastc0+1
                  if0=ic0*ratio(0)
                     arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                                +arrayf(if0,if1,if2)*areaf
                  enddo
               enddo
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2 

         do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0+1
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/areac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgsidecplx3d1(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double complex
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1+1,
     &          filo2:fihi2),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1+1,
     &          cilo2:cihi2)
      double precision areaf,areac
      integer ic0,ic1,ic2,if0,if1,if2,

     &  ir0,ir2
c
c***********************************************************************
c
      areaf=dxf(2)*dxf(0)
      areac=dxc(2)*dxc(0)


      do ic2=ifirstc2,ilastc2 

         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0

               arrayc(ic0,ic1,ic2)=cmplx(zero,zero)
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2 
         do ir2=0,ratio(2)-1
               if2=ic2*ratio(2)+ir2

         do ic1=ifirstc1,ilastc1+1
                  if1=ic1*ratio(1)

            do ic0=ifirstc0,ilastc0
               do ir0=0,ratio(0)-1
                  if0=ic0*ratio(0)+ir0
                     arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                                +arrayf(if0,if1,if2)*areaf
                  enddo
               enddo
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2 

         do ic1=ifirstc1,ilastc1+1

            do ic0=ifirstc0,ilastc0
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/areac
            enddo
         enddo
      enddo
c
      return
      end
c
      subroutine cartwgtavgsidecplx3d2(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
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
     &  dxf(0:3-1),
     &  dxc(0:3-1)
      double complex
     &  arrayf(filo0:fihi0,
     &          filo1:fihi1,
     &          filo2:fihi2+1),
     &  arrayc(cilo0:cihi0,
     &          cilo1:cihi1,
     &          cilo2:cihi2+1)
      double precision areaf,areac
      integer ic0,ic1,ic2,if0,if1,if2,

     &  ir0,ir1
c
c***********************************************************************
c
      areaf=dxf(0)*dxf(1)
      areac=dxc(0)*dxc(1)


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0

               arrayc(ic0,ic1,ic2)=cmplx(zero,zero)
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1
               if2=ic2*ratio(2)

         do ic1=ifirstc1,ilastc1
            do ir1=0,ratio(1)-1
                  if1=ic1*ratio(1)+ir1

            do ic0=ifirstc0,ilastc0
               do ir0=0,ratio(0)-1
                  if0=ic0*ratio(0)+ir0
                     arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                                +arrayf(if0,if1,if2)*areaf
                  enddo
               enddo
            enddo
         enddo
      enddo


      do ic2=ifirstc2,ilastc2+1

         do ic1=ifirstc1,ilastc1

            do ic0=ifirstc0,ilastc0
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/areac
            enddo
         enddo
      enddo
c
      return
      end
c
