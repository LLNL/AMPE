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
c        Weighted averaging for 3d cell-centered quaternion data
c        
      subroutine quatcoarsen(
     &   ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &   filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &   cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &   ratio,dxf,dxc,
     &   qfine,qcoarse,
     &   depth,
     &   symmetry_aware,
     &   iqrot_x, iqrot_y, iqrot_z,
     &   iqlo0, iqhi0, iqlo1, iqhi1, iqlo2, iqhi2 )
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c        
      integer
     &   ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &   filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &   cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &   depth,
     &   symmetry_aware,
     &   iqlo0, iqhi0, iqlo1, iqhi1, iqlo2, iqhi2
      integer ratio(0:3-1)
      double precision
     &   dxf(0:3-1),
     &   dxc(0:3-1)
      double precision
     &   qfine(filo0:fihi0,
     &   filo1:fihi1,
     &   filo2:fihi2,
     &   depth),
     &   qcoarse(cilo0:cihi0,
     &   cilo1:cihi1,
     &   cilo2:cihi2,
     &   depth)
      integer
     &   iqrot_x(iqlo0:iqhi0+1,iqlo1:iqhi1,iqlo2:iqhi2),
     &   iqrot_y(iqlo0:iqhi0,iqlo1:iqhi1+1,iqlo2:iqhi2),
     &   iqrot_z(iqlo0:iqhi0,iqlo1:iqhi1,iqlo2:iqhi2+1)
c
      double precision dVf,dVcinv
      integer ic0,ic1,ic2,if0,if1,if2,ir0,ir1,ir2
      double precision q0(depth)
      double precision q1(depth), q1_prime(depth)
      integer iq, m
c        
      iq = 0
c        
c***********************************************************************
c        
      dVf = dxf(0)*dxf(1)*dxf(2)
      dVcinv = 1.d0/(dxc(0)*dxc(1)*dxc(2))

      do ic2 = ifirstc2, ilastc2
         do ic1 = ifirstc1, ilastc1
            do ic0 = ifirstc0, ilastc0

               do m = 1, depth
                  qcoarse(ic0,ic1,ic2,m) = zero
                  if0 = ic0 * ratio(0)
                  if1 = ic1 * ratio(1)
                  if2 = ic2 * ratio(2)
                  q0(m) = qfine(if0,if1,if2,m)
               enddo

               do ir2 = 0, ratio(2) - 1
                  do ir1 = 0, ratio(1) - 1
                     do ir0 = 0, ratio(0) - 1

                        if0 = ic0 * ratio(0) + ir0
                        if1 = ic1 * ratio(1) + ir1
                        if2 = ic2 * ratio(2) + ir2

                        if ( symmetry_aware /= 0 )then
                           if ( ir0 == 0 .and. ir1 == 0 .and.
     &                        ir2 == 0 ) then
                              do m = 1, depth
                                 q1_prime(m) = qfine(if0,if1,if2,m)
                              enddo
                           else
                              do m = 1, depth
                                 q1(m) = qfine(if0,if1,if2,m)
                              enddo
                              call quatsymmrotate( q1, iq, q1_prime,
     &                           depth )
                           endif
                        else
                           do m = 1, depth
                              q1_prime(m) = qfine(if0,if1,if2,m)
                           enddo
                        endif

                        do m = 1, depth
                           qcoarse(ic0,ic1,ic2,m) =
     &                        qcoarse(ic0,ic1,ic2,m) +
     &                        q1_prime(m) * dVf
                        enddo

                     enddo
                  enddo
               enddo

               do m = 1, depth
                  q0(m) = qcoarse(ic0,ic1,ic2,m) * dVcinv
               enddo

               call quatnorm( q0, depth )

               do m = 1, depth
                  qcoarse(ic0,ic1,ic2,m) = q0(m)
               enddo

            enddo
         enddo
      enddo
c        
      return
      end
