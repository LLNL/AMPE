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
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  ilo0,ilo1,ihi0,ihi1,
     &  array)
c***********************************************************************
      implicit none
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  ilo0,ilo1,ihi0,ihi1
      double precision
     &  array(ilo0:ihi0,
     &        ilo1:ihi1)
      integer i0,i1
c
c***********************************************************************
c
      do i1=ifirst1,ilast1
         do i0=ifirst0,ilast0
            array(i0,i1)=0.
         enddo
      enddo
c
      return
      end
c
