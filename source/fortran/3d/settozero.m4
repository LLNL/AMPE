c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
c 
c***********************************************************************
c
       subroutine settozero(
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  arrayf)
c***********************************************************************
      implicit none
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      double precision
     &  arrayf(filo0:fihi0,
     &         filo1:fihi1,
     &         filo2:fihi2)
      integer if0,if1,if2
c
c***********************************************************************
c
c      print*,'settozero for index '
c      print*,ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
c
      do if2=ifirst2,ilast2
         do if1=ifirst1,ilast1
            do if0=ifirst0,ilast0
               arrayf(if0,if1,if2)=0.
            enddo
         enddo
      enddo
c
      return
      end
c
