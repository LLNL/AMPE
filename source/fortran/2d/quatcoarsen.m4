c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
c 
c        Weighted averaging for 2d cell-centered quaternion data
c        
      subroutine quatcoarsen(
     &   ifirstc0,ifirstc1,ilastc0,ilastc1,
     &   filo0,filo1,fihi0,fihi1,
     &   cilo0,cilo1,cihi0,cihi1,
     &   ratio,dxf,dxc,
     &   qfine,qcoarse,
     &   depth,
     &   symmetry_aware,
     &   iqrot_x, iqrot_y,
     &   iqlo0, iqhi0, iqlo1, iqhi1 )
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c        
      integer
     &   ifirstc0,ifirstc1,ilastc0,ilastc1,
     &   cilo0,cilo1,cihi0,cihi1,
     &   filo0,filo1,fihi0,fihi1,
     &   depth,
     &   symmetry_aware,
     &   iqlo0, iqhi0, iqlo1, iqhi1
      integer ratio(0:2-1)
      double precision
     &   dxf(0:2-1),
     &   dxc(0:2-1)
      double precision
     &   qfine(filo0:fihi0,
     &   filo1:fihi1,
     &   depth),
     &   qcoarse(cilo0:cihi0,
     &   cilo1:cihi1,
     &   depth)
      integer
     &   iqrot_x(iqlo0:iqhi0+1,iqlo1:iqhi1),
     &   iqrot_y(iqlo0:iqhi0,iqlo1+1:iqhi1)
c
      double precision dVf,dVcinv
      integer ic0,ic1,if0,if1,ir0,ir1
      double precision q0(depth)
      double precision q1(depth), q1_prime(depth)
      integer iq, m
c        
      iq = 0
c        
c***********************************************************************
c        
      dVf = dxf(0)*dxf(1)
      dVcinv = 1.d0/(dxc(0)*dxc(1))

      do ic1 = ifirstc1, ilastc1
         do ic0 = ifirstc0, ilastc0

            do m = 1, depth
               qcoarse(ic0,ic1,m) = zero
               if0 = ic0 * ratio(0)
               if1 = ic1 * ratio(1)
               q0(m) = qfine(if0,if1,m)
            enddo

            do ir1 = 0, ratio(1) - 1
               do ir0 = 0, ratio(0) - 1

                  if0 = ic0 * ratio(0) + ir0
                  if1 = ic1 * ratio(1) + ir1

                  if ( symmetry_aware /= 0 )then
                     if ( ir1 == 0 .and. ir0 == 0 ) then
                        do m = 1, depth
                           q1_prime(m) = qfine(if0,if1,m)
                        enddo
                     else
                        do m = 1, depth
                           q1(m) = qfine(if0,if1,m)
                        enddo
                        call quatsymmrotate( q1, iq, q1_prime, depth )
                     endif
                  else
                     do m = 1, depth
                        q1_prime(m) = qfine(if0,if1,m)
                     enddo
                  endif

                  do m = 1, depth
                     qcoarse(ic0,ic1,m) = qcoarse(ic0,ic1,m) +
     &                  q1_prime(m) * dVf
                  enddo

               enddo
            enddo

            do m = 1, depth
               q0(m) = qcoarse(ic0,ic1,m) * dVcinv
            enddo

            call quatnorm( q0, depth )

            do m = 1, depth
               qcoarse(ic0,ic1,m) = q0(m)
            enddo

         enddo
      enddo
c        
      return
      end
c        
