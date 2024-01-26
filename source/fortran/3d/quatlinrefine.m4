c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
c
c***********************************************************************
c Linear interpolation for 3d cell-centered double data
c***********************************************************************
c
      subroutine quatlinrefine3d(
     &   ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &   cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &   filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &   ratio,dxc,dxf,
     &   arrayc,arrayf,
     &   depth,
     &   iqrot_x, iqrot_y, iqrot_z,
     &   iqlo0, iqhi0, iqlo1, iqhi1, iqlo2, iqhi2 )
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.d0)
c        
      integer
     &   ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &   cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &   filo0,filo1,filo2,fihi0,fihi1,fihi2,depth,
     &   iqlo0, iqhi0, iqlo1, iqhi1, iqlo2, iqhi2
      integer ratio(0:3-1)
      double precision
     &   dxc(0:3-1),
     &   dxf(0:3-1)
      double precision
     &   arrayc(cilo0:cihi0,
     &   cilo1:cihi1,
     &   cilo2:cihi2,depth),
     &   arrayf(filo0:fihi0,
     &   filo1:fihi1,
     &   filo2:fihi2,depth)
      integer
     &   iqrot_x(iqlo0:iqhi0+1,iqlo1:iqhi1,iqlo2:iqhi2),
     &   iqrot_y(iqlo0:iqhi0,iqlo1:iqhi1+1,iqlo2:iqhi2),
     &   iqrot_z(iqlo0:iqhi0,iqlo1:iqhi1,iqlo2:iqhi2+1)
c        
      double precision deltax(0:15,0:3-1),x,y,z
      integer ic0,ic1,ic2,if0,if1,if2,ir0,ir1,ir2,iref0,iref1,iref2
      double precision q1(depth), q2(depth), q3(depth), q4(depth),
     &   q5(depth), q6(depth), q7(depth), q8(depth)
      double precision qref(depth)
      double precision q1_prime(depth), q2_prime(depth), 
     &   q3_prime(depth), q4_prime(depth), 
     &   q5_prime(depth), q6_prime(depth),
     &   q7_prime(depth), q8_prime(depth)
      integer iq,m
c        
      iq = 0
c        
c***********************************************************************
c        
      do ir0=0,ratio(0)-1
         deltax(ir0,0)=(dble(ir0)+half)*dxf(0)-dxc(0)*half
      enddo

      do ir1=0,ratio(1)-1
         deltax(ir1,1)=(dble(ir1)+half)*dxf(1)-dxc(1)*half
      enddo

      do ir2=0,ratio(2)-1
         deltax(ir2,2)=(dble(ir2)+half)*dxf(2)-dxc(2)*half
      enddo

      do if2=ifirstf2,ilastf2
         if (if2.lt.0) then
            ic2=(if2+1)/ratio(2)-1
         else
            ic2=if2/ratio(2)
         endif
         iref2=ic2
         ir2=if2-ic2*ratio(2)
         z=deltax(ir2,2)/dxc(2)
         if( z .lt. 0.d0 ) then
            ic2 = ic2-1
            z = z + one
         endif
         do if1=ifirstf1,ilastf1
            if (if1.lt.0) then
               ic1=(if1+1)/ratio(1)-1
            else
               ic1=if1/ratio(1)
            endif
            iref1=ic1
            ir1=if1-ic1*ratio(1)
            y=deltax(ir1,1)/dxc(1)
            if( y .lt. 0.d0 ) then
               ic1 = ic1-1
               y = y + one
            endif
            do if0=ifirstf0,ilastf0
               if (if0.lt.0) then
                  ic0=(if0+1)/ratio(0)-1
               else
                  ic0=if0/ratio(0)
               endif
               iref0=ic0
               ir0=if0-ic0*ratio(0)
               x=deltax(ir0,0)/dxc(0)
               if( x .lt. 0.d0 ) then
                  ic0 = ic0-1
                  x = x + one
               endif
               do m = 1, depth
                  qref(m) = arrayc(iref0,iref1,iref2,m)
                  q1(m) = arrayc(ic0,  ic1,  ic2,  m)
                  q2(m) = arrayc(ic0,  ic1+1,ic2,  m)
                  q3(m) = arrayc(ic0+1,ic1,  ic2,  m)
                  q4(m) = arrayc(ic0+1,ic1+1,ic2,  m)
                  q5(m) = arrayc(ic0,  ic1,  ic2+1,m)
                  q6(m) = arrayc(ic0,  ic1+1,ic2+1,m)
                  q7(m) = arrayc(ic0+1,ic1,  ic2+1,m)
                  q8(m) = arrayc(ic0+1,ic1+1,ic2+1,m)
               enddo
               call quatsymmrotate( q1, iq, q1_prime, depth )
               call quatsymmrotate( q2, iq, q2_prime, depth )
               call quatsymmrotate( q3, iq, q3_prime, depth )
               call quatsymmrotate( q4, iq, q4_prime, depth )
               call quatsymmrotate( q5, iq, q5_prime, depth )
               call quatsymmrotate( q6, iq, q6_prime, depth )
               call quatsymmrotate( q7, iq, q7_prime, depth )
               call quatsymmrotate( q8, iq, q8_prime, depth )
               do m = 1, depth
                  q1(m)=
     &               ( (q1_prime(m)
     &               +(q3_prime(m)-q1_prime(m))*x)*(one-y)
     &               +(q2_prime(m)
     &               +(q4_prime(m)-q2_prime(m))*x)*y )
     &               *(one-z)
     &               +( (q5_prime(m)
     &               +(q7_prime(m)-q5_prime(m))*x)
     &               *(one-y)
     &               +(q6_prime(m)
     &               +(q8_prime(m)-q6_prime(m))*x)*y)
     &               *z
               enddo
               call quatnorm( q1, depth )
               do m = 1, depth
                  arrayf(if0,if1,if2,m)=q1(m)
               enddo
            enddo
         enddo
      enddo
c        
      return
      end
c        
