c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
c Redistribution and use in source and binary forms, with or without 
c modification, are permitted provided that the following conditions are met:
c - Redistributions of source code must retain the above copyright notice,
c   this list of conditions and the disclaimer below.
c - Redistributions in binary form must reproduce the above copyright notice,
c   this list of conditions and the disclaimer (as noted below) in the
c   documentation and/or other materials provided with the distribution.
c - Neither the name of the LLNS/LLNL nor the names of its contributors may be
c   used to endorse or promote products derived from this software without
c   specific prior written permission.
c
c THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
c AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
c IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
c ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
c LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
c DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
c DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
c OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
c HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
c STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
c IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
c POSSIBILITY OF SUCH DAMAGE.
c 
c $Id:$

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
