c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
c 
      subroutine quatfindsymm( q1, q2, iq, q2_prime, qlen )

      implicit none

      integer iq, qlen
      double precision q1(qlen), q2(qlen)
      double precision q2_prime(qlen)

      if ( qlen == 4 ) then

         call quatfindsymm4( q1, q2, iq, q2_prime )

      else if ( qlen == 2 ) then

         call quatfindsymm2( q1, q2, iq, q2_prime )

      else if ( qlen == 1 ) then

         call quatfindsymm1( q1(1), q2(1), iq, q2_prime(1) )

      else
         
         print *, "Error in quatfindsymm"
         stop

      endif

      return
      end

c=======================================================================

      subroutine quatsymmrotate( q, iq, q_prime, qlen )

      implicit none

      integer iq, qlen
      double precision q(qlen)
      double precision q_prime(qlen)

      if ( qlen == 4 ) then

         call quatsymmrotate4( q, iq, q_prime )

      else if ( qlen == 2 ) then

         call quatsymmrotate2( q, iq, q_prime )

      else if ( qlen == 1 ) then

         call quatsymmrotate1( q(1), iq, q_prime(1) )

      else
         
         print *, "Error in quatsymmrotate, qlen = ", qlen
         stop

      endif

      return
      end

c=======================================================================

      subroutine quatfindsymm4( q1, q2, iq, q2_prime )

      implicit none

      integer NROTATIONS
      parameter( NROTATIONS=48 )
      integer QLEN
      parameter( QLEN=4 )

      integer iq
      double precision q1(QLEN), q2(QLEN)
      double precision q2_prime(QLEN)

      integer nn

      double precision qr(QLEN,NROTATIONS)
      save qr
      integer iq_qr_conj(NROTATIONS)
      save iq_qr_conj
      double precision dsq
      double precision min_dsq
      double precision min_q2_prime(QLEN)
      integer min_iq

      logical first_pass
      data first_pass / .true. /
      save first_pass

      double precision ZERO, ONE, PI
      double precision TWO_TIMES_SIN_OF_PI_OVER_16
      double precision SQ_TWO_TIMES_SIN_OF_PI_OVER_16
      parameter( ZERO=0.0d0 )
      parameter( ONE=1.0d0 )
      save PI
      save TWO_TIMES_SIN_OF_PI_OVER_16
      save SQ_TWO_TIMES_SIN_OF_PI_OVER_16

      if ( first_pass ) then

         PI = dacos( -1.0d0 )
         TWO_TIMES_SIN_OF_PI_OVER_16 = 2.0d0 * dsin( PI / 16.0d0 )
         SQ_TWO_TIMES_SIN_OF_PI_OVER_16 = 
     &      TWO_TIMES_SIN_OF_PI_OVER_16 * TWO_TIMES_SIN_OF_PI_OVER_16

c           Construct all qr
         call setqr( qr, iq_qr_conj )
         
         first_pass = .false.
      endif

      if ( iq == 0 .or. iq > NROTATIONS .or. iq < -NROTATIONS ) iq = 1
     
c        if iq is negative, use the conjugate of that rotation
      if ( iq < 0 ) iq = iq_qr_conj(-iq)

c      Test qr(:,iq) first
      if ( iq == 1 ) then
         call quatdiffsq4( q1, q2, dsq )
         call quatcopy4( q2, q2_prime )
      else
         call quatrotatediffsq4( q1, q2, qr(1,iq), q2_prime, dsq )
      endif
      if ( dsq .le. SQ_TWO_TIMES_SIN_OF_PI_OVER_16 ) return

      min_dsq = dsq
      min_iq = iq
      call quatcopy4( q2_prime, min_q2_prime )

      do nn = 1, NROTATIONS
         if ( nn == iq ) cycle
         if ( nn == 1 ) then
            call quatdiffsq4( q1, q2, dsq )
            call quatcopy4( q2, q2_prime )
         else
            call quatrotatediffsq4( q1, q2, qr(1,nn), q2_prime, dsq )
         endif
         if ( dsq < min_dsq ) then
            min_dsq = dsq
            min_iq = nn
            call quatcopy4( q2_prime, min_q2_prime )
         endif            
         if ( dsq <= SQ_TWO_TIMES_SIN_OF_PI_OVER_16 ) exit
      enddo

      call quatcopy4( min_q2_prime, q2_prime )
      iq = min_iq

      return
      end

c=======================================================================

      subroutine setqr( qr, iq_qr_conj )

      implicit none

      integer NROTATIONS
      parameter( NROTATIONS=48 )
      integer QLEN
      parameter( QLEN=4 )

      double precision ZERO, ONE
      parameter( ZERO=0.0d0 )
      parameter( ONE =1.0d0 )

      double precision qr(QLEN,NROTATIONS)
      integer iq_qr_conj(NROTATIONS)

c           \/\/ Reorder these in light of using fundamental_quat
c           \/\/ so that likely rotations are checked first

      call quatset( qr(1,1),   ONE, ZERO, ZERO, ZERO )
      call quatset( qr(1,2),  ZERO,  ONE, ZERO, ZERO )
      call quatset( qr(1,3),  ZERO, ZERO,  ONE, ZERO )
      call quatset( qr(1,4),  ZERO, ZERO, ZERO,  ONE )
      call quatset( qr(1,5),  -ONE, ZERO, ZERO, ZERO )
      call quatset( qr(1,6),  ZERO, -ONE, ZERO, ZERO )
      call quatset( qr(1,7),  ZERO, ZERO, -ONE, ZERO )
      call quatset( qr(1,8),  ZERO, ZERO, ZERO, -ONE )

      call quatset( qr(1,9),   ONE,  ONE, ZERO, ZERO )
      call quatset( qr(1,10),  ONE, ZERO,  ONE, ZERO )
      call quatset( qr(1,11),  ONE, ZERO, ZERO,  ONE )
      call quatset( qr(1,12), ZERO,  ONE,  ONE, ZERO )
      call quatset( qr(1,13), ZERO,  ONE, ZERO,  ONE )
      call quatset( qr(1,14), ZERO, ZERO,  ONE,  ONE )
      call quatset( qr(1,15), -ONE,  ONE, ZERO, ZERO )
      call quatset( qr(1,16), -ONE, ZERO,  ONE, ZERO )
      call quatset( qr(1,17), -ONE, ZERO, ZERO,  ONE )
      call quatset( qr(1,18), ZERO, -ONE,  ONE, ZERO )
      call quatset( qr(1,19), ZERO, -ONE, ZERO,  ONE )
      call quatset( qr(1,20), ZERO, ZERO, -ONE,  ONE )
      call quatset( qr(1,21),  ONE, -ONE, ZERO, ZERO )
      call quatset( qr(1,22),  ONE, ZERO, -ONE, ZERO )
      call quatset( qr(1,23),  ONE, ZERO, ZERO, -ONE )
      call quatset( qr(1,24), ZERO,  ONE, -ONE, ZERO )
      call quatset( qr(1,25), ZERO,  ONE, ZERO, -ONE )
      call quatset( qr(1,26), ZERO, ZERO,  ONE, -ONE )
      call quatset( qr(1,27), -ONE, -ONE, ZERO, ZERO )
      call quatset( qr(1,28), -ONE, ZERO, -ONE, ZERO )
      call quatset( qr(1,29), -ONE, ZERO, ZERO, -ONE )
      call quatset( qr(1,30), ZERO, -ONE, -ONE, ZERO )
      call quatset( qr(1,31), ZERO, -ONE, ZERO, -ONE )
      call quatset( qr(1,32), ZERO, ZERO, -ONE, -ONE )

      call quatset( qr(1,33),  ONE,  ONE,  ONE,  ONE )
      call quatset( qr(1,34), -ONE,  ONE,  ONE,  ONE )
      call quatset( qr(1,35),  ONE, -ONE,  ONE,  ONE )
      call quatset( qr(1,36),  ONE,  ONE, -ONE,  ONE )
      call quatset( qr(1,37),  ONE,  ONE,  ONE, -ONE )
      call quatset( qr(1,38), -ONE, -ONE,  ONE,  ONE )
      call quatset( qr(1,39), -ONE,  ONE, -ONE,  ONE )
      call quatset( qr(1,40), -ONE,  ONE,  ONE, -ONE )
      call quatset( qr(1,41),  ONE, -ONE, -ONE,  ONE )
      call quatset( qr(1,42),  ONE, -ONE,  ONE, -ONE )
      call quatset( qr(1,43),  ONE,  ONE, -ONE, -ONE )
      call quatset( qr(1,44),  ONE, -ONE, -ONE, -ONE )
      call quatset( qr(1,45), -ONE,  ONE, -ONE, -ONE )
      call quatset( qr(1,46), -ONE, -ONE,  ONE, -ONE )
      call quatset( qr(1,47), -ONE, -ONE, -ONE,  ONE )
      call quatset( qr(1,48), -ONE, -ONE, -ONE, -ONE )

c        identify index of conjugates
      iq_qr_conj(1) = 1
      iq_qr_conj(2) = 6
      iq_qr_conj(3) = 7
      iq_qr_conj(4) = 8
      iq_qr_conj(5) = 5
      iq_qr_conj(6) = 2
      iq_qr_conj(7) = 3
      iq_qr_conj(8) = 4
      iq_qr_conj(9) = 21
      iq_qr_conj(10) = 22
      iq_qr_conj(11) = 23
      iq_qr_conj(12) = 30
      iq_qr_conj(13) = 31
      iq_qr_conj(14) = 32
      iq_qr_conj(15) = 27
      iq_qr_conj(16) = 28
      iq_qr_conj(17) = 29
      iq_qr_conj(18) = 24
      iq_qr_conj(19) = 25
      iq_qr_conj(20) = 26
      iq_qr_conj(21) = 9
      iq_qr_conj(22) = 10
      iq_qr_conj(23) = 11
      iq_qr_conj(24) = 18
      iq_qr_conj(25) = 19
      iq_qr_conj(26) = 20
      iq_qr_conj(27) = 15
      iq_qr_conj(28) = 16
      iq_qr_conj(29) = 17
      iq_qr_conj(30) = 12
      iq_qr_conj(31) = 13
      iq_qr_conj(32) = 14
      iq_qr_conj(33) = 44
      iq_qr_conj(34) = 48
      iq_qr_conj(35) = 43
      iq_qr_conj(36) = 42
      iq_qr_conj(37) = 41
      iq_qr_conj(38) = 45
      iq_qr_conj(39) = 46
      iq_qr_conj(40) = 47
      iq_qr_conj(41) = 37
      iq_qr_conj(42) = 36
      iq_qr_conj(43) = 35
      iq_qr_conj(44) = 33
      iq_qr_conj(45) = 38
      iq_qr_conj(46) = 39
      iq_qr_conj(47) = 40
      iq_qr_conj(48) = 34
      
      return
      end

c=======================================================================

      subroutine quatsymmrotate4( q, iq, q_prime )

      implicit none

      integer NROTATIONS
      parameter( NROTATIONS=48 )
      integer QLEN
      parameter( QLEN=4 )

      integer iq
      double precision q(QLEN)
      double precision q_prime(QLEN)

      double precision qr(QLEN,NROTATIONS)
      save qr
      integer iq_qr_conj(NROTATIONS)
      save iq_qr_conj

      logical first_pass
      data first_pass / .true. /
      save first_pass

      if ( first_pass ) then

c           Construct all qr

c           \/\/ Reorder these in light of using fundamental_quat
c           \/\/ so that likely rotations are checked first

         call setqr( qr, iq_qr_conj )

         first_pass = .false.
      endif

      if ( iq == 0 .or. iq > NROTATIONS .or. iq < -NROTATIONS ) then
         print *, "Error in quatsymmrotate4, iq=", iq
         stop
      endif
     
c        if iq is negative, use the conjugate of that rotation
      if ( iq < 0 ) iq = iq_qr_conj(-iq)

      if ( iq == 1 ) then
         call quatcopy4( q, q_prime )
      else
         call quatmult4( q, qr(1,iq), q_prime )
      endif

      return
      end

c=======================================================================

      subroutine quatfindsymm2( q1, q2, iq, q2_prime )

      implicit none

      integer NROTATIONS
      parameter( NROTATIONS=4 )
      integer QLEN
      parameter( QLEN=2 )

      integer iq
      double precision q1(QLEN), q2(QLEN)
      double precision q2_prime(QLEN)

      integer nn

      double precision qr(QLEN,NROTATIONS)
      save qr
      double precision q0(QLEN)
      save q0
      integer iq_qr_conj(NROTATIONS)
      save iq_qr_conj
      double precision dsq
      double precision min_dsq
      double precision min_q2_prime(QLEN)
      integer min_iq

      logical first_pass
      data first_pass / .true. /
      save first_pass

      double precision ZERO, ONE, PI
      double precision TWO_TIMES_SIN_OF_PI_OVER_8
      double precision SQ_TWO_TIMES_SIN_OF_PI_OVER_8
      parameter( ZERO=0.0d0 )
      parameter( ONE=1.0d0 )
      save PI
      save TWO_TIMES_SIN_OF_PI_OVER_8
      save SQ_TWO_TIMES_SIN_OF_PI_OVER_8

      if ( first_pass ) then

         PI = dacos( -1.0d0 )
         TWO_TIMES_SIN_OF_PI_OVER_8 = 2.0d0 * dsin( PI / 8.0d0 )
         SQ_TWO_TIMES_SIN_OF_PI_OVER_8 = 
     &      TWO_TIMES_SIN_OF_PI_OVER_8 * TWO_TIMES_SIN_OF_PI_OVER_8

c           Construct all qr

         call quatset2( qr(1,1),   ONE, ZERO )
         call quatset2( qr(1,2),  ZERO,  ONE )
         call quatset2( qr(1,3),  -ONE, ZERO )
         call quatset2( qr(1,4),  ZERO, -ONE )

c           identify index of conjugates
         iq_qr_conj(1) = 1
         iq_qr_conj(2) = 4
         iq_qr_conj(3) = 3
         iq_qr_conj(4) = 2

         call quatset2( q0, ONE, ZERO )

         first_pass = .false.
      endif

      if ( iq == 0 .or. iq > NROTATIONS .or. iq < -NROTATIONS ) iq = 1
     
c        if iq is negative, use the conjugate of that rotation
      if ( iq < 0 ) iq = iq_qr_conj(-iq)

c      Test qr(:,iq) first
      if ( iq == 1 ) then
         call quatdiffsq2( q1, q2, dsq )
         call quatcopy2( q2, q2_prime )
      else
         call quatrotatediffsq2( q1, q2, qr(1,iq), q2_prime, dsq )
      endif
      if ( dsq .le. SQ_TWO_TIMES_SIN_OF_PI_OVER_8 ) return

      min_dsq = dsq
      min_iq = iq
      call quatcopy2( q2_prime, min_q2_prime )

      do nn = 1, NROTATIONS
         if ( nn == iq ) cycle
         if ( nn == 1 ) then
            call quatdiffsq2( q1, q2, dsq )
            call quatcopy2( q2, q2_prime )
         else
            call quatrotatediffsq2( q1, q2, qr(1,nn), q2_prime, dsq )
         endif
         if ( dsq < min_dsq ) then
            min_dsq = dsq
            min_iq = nn
            call quatcopy2( q2_prime, min_q2_prime )
         endif            
         if ( dsq <= SQ_TWO_TIMES_SIN_OF_PI_OVER_8 ) exit
      enddo

      call quatcopy2( min_q2_prime, q2_prime )
      iq = min_iq

      return
      end

c=======================================================================

      subroutine quatsymmrotate2( q, iq, q_prime )

      implicit none

      integer NROTATIONS
      parameter( NROTATIONS=4 )
      integer QLEN
      parameter( QLEN=2 )

      integer iq
      double precision q(QLEN)
      double precision q_prime(QLEN)

      double precision qr(QLEN,NROTATIONS)
      save qr
      integer iq_qr_conj(NROTATIONS)
      save iq_qr_conj

      logical first_pass
      data first_pass / .true. /
      save first_pass

      double precision ZERO, ONE, PI
      double precision TWO_TIMES_SIN_OF_PI_OVER_8
      double precision SQ_TWO_TIMES_SIN_OF_PI_OVER_8
      parameter( ZERO=0.0d0 )
      parameter( ONE=1.0d0 )
      save PI
      save TWO_TIMES_SIN_OF_PI_OVER_8
      save SQ_TWO_TIMES_SIN_OF_PI_OVER_8

      if ( first_pass ) then

         PI = dacos( -1.0d0 )
         TWO_TIMES_SIN_OF_PI_OVER_8 = 2.0d0 * dsin( PI / 8.0d0 )
         SQ_TWO_TIMES_SIN_OF_PI_OVER_8 = 
     &      TWO_TIMES_SIN_OF_PI_OVER_8 * TWO_TIMES_SIN_OF_PI_OVER_8

c           Construct all qr

         call quatset2( qr(1,1),   ONE, ZERO )
         call quatset2( qr(1,2),  ZERO,  ONE )
         call quatset2( qr(1,3),  -ONE, ZERO )
         call quatset2( qr(1,4),  ZERO, -ONE )

c           identify index of conjugates
         iq_qr_conj(1) = 1
         iq_qr_conj(2) = 4
         iq_qr_conj(3) = 3
         iq_qr_conj(4) = 2

         first_pass = .false.
      endif

      if ( iq == 0 .or. iq > NROTATIONS .or. iq < -NROTATIONS ) then
         print *, "Error in quatsymmrotate2"
         stop
      endif

c        if iq is negative, use the conjugate of that rotation
      if ( iq < 0 ) iq = iq_qr_conj(-iq)

      if ( iq == 1 ) then
         call quatcopy2( q, q_prime )
      else
         call quatmult2( q, qr(1,iq), q_prime )
      endif

      return
      end

c=======================================================================

      subroutine quatfindsymm1( q1, q2, iq, q2_prime )

      implicit none

      integer NROTATIONS
      parameter( NROTATIONS=9 )

      integer iq
      double precision q1, q2
      double precision q2_prime

      integer nn

      double precision qr(NROTATIONS)
      save qr
      integer iq_qr_conj(NROTATIONS)
      save iq_qr_conj
      double precision qdabs
      double precision min_qdabs
      double precision min_q2_prime
      integer min_iq

      logical first_pass
      data first_pass / .true. /
      save first_pass

      double precision PI, TWO_PI
      double precision PI_OVER_2, PI_OVER_4
      double precision THREE_PI_OVER_2
      save PI, TWO_PI, PI_OVER_2, PI_OVER_4, THREE_PI_OVER_2

      if ( first_pass ) then

         PI = dacos( -1.0d0 )
         PI_OVER_2 = 0.5d0 * PI
         PI_OVER_4 = 0.25d0 * PI
         THREE_PI_OVER_2 = 1.5d0 * PI
         TWO_PI = 2.0d0 * PI

c           Construct all qr
         qr(1) = 0.0
         qr(2) = PI_OVER_2
         qr(3) = -PI_OVER_2
         qr(4) = PI
         qr(5) = -PI
         qr(6) = THREE_PI_OVER_2
         qr(7) = -THREE_PI_OVER_2
         qr(8) = TWO_PI
         qr(9) = -TWO_PI

c           identify index of conjugates
         iq_qr_conj(1) = 1
         iq_qr_conj(2) = 3
         iq_qr_conj(3) = 2
         iq_qr_conj(4) = 5
         iq_qr_conj(5) = 4
         iq_qr_conj(6) = 7
         iq_qr_conj(7) = 6
         iq_qr_conj(8) = 9
         iq_qr_conj(9) = 8

         first_pass = .false.
      endif

      if ( iq == 0 .or. iq > NROTATIONS .or. iq < -NROTATIONS ) iq = 1
     
c        if iq is negative, use the conjugate of that rotation
      if ( iq < 0 ) iq = iq_qr_conj(-iq)

c        Test qr(:,iq) first
      if ( iq == 1 ) then
         call quatdiffabs1( q1, q2, qdabs )
         q2_prime = q2
      else
         call quatrotatediffabs1( q1, q2, qr(iq), q2_prime, qdabs )
      endif
      if ( qdabs <= PI_OVER_4 ) return

      min_qdabs = qdabs
      min_iq = iq
      min_q2_prime = q2_prime

      do nn = 1, NROTATIONS
         if ( nn == iq ) cycle
         if ( nn == 1 ) then
            call quatdiffabs1( q1, q2, qdabs )
            q2_prime = q2
         else
            call quatrotatediffabs1( q1, q2, qr(nn), q2_prime, qdabs )
         endif
         if ( qdabs < min_qdabs ) then
            min_qdabs = qdabs
            min_iq = nn
            min_q2_prime = q2_prime
         endif            
         if ( qdabs <= PI_OVER_4 ) exit
      enddo

      q2_prime = min_q2_prime
      iq = min_iq

      return
      end

c=======================================================================

      subroutine quatsymmrotate1( q, iq, q_prime )

      implicit none

      integer NROTATIONS
      parameter( NROTATIONS=9 )

      integer iq
      double precision q
      double precision q_prime

      double precision qr(NROTATIONS)
      save qr
      integer iq_qr_conj(NROTATIONS)
      save iq_qr_conj

      logical first_pass
      data first_pass / .true. /
      save first_pass

      double precision PI, TWO_PI
      double precision PI_OVER_2, PI_OVER_4
      double precision THREE_PI_OVER_2
      save PI, TWO_PI, PI_OVER_2, PI_OVER_4, THREE_PI_OVER_2

      if ( first_pass ) then

         PI = dacos( -1.0d0 )
         PI_OVER_2 = 0.5d0 * PI
         PI_OVER_4 = 0.25d0 * PI
         THREE_PI_OVER_2 = 1.5d0 * PI
         TWO_PI = 2.0d0 * PI

c           Construct all qr
         qr(1) = 0.0
         qr(2) = PI_OVER_2
         qr(3) = -PI_OVER_2
         qr(4) = PI
         qr(5) = -PI
         qr(6) = THREE_PI_OVER_2
         qr(7) = -THREE_PI_OVER_2
         qr(8) = TWO_PI
         qr(9) = -TWO_PI

c           identify index of conjugates
         iq_qr_conj(1) = 1
         iq_qr_conj(2) = 3
         iq_qr_conj(3) = 2
         iq_qr_conj(4) = 5
         iq_qr_conj(5) = 4
         iq_qr_conj(6) = 7
         iq_qr_conj(7) = 6
         iq_qr_conj(8) = 9
         iq_qr_conj(9) = 8

         first_pass = .false.
      endif

      if ( iq == 0 .or. iq > NROTATIONS .or. iq < -NROTATIONS ) then
         print *, "Error in quatsymmrotate1"
         stop
      endif

c        if iq is negative, use the conjugate of that rotation
      if ( iq < 0 ) iq = iq_qr_conj(-iq)

      q_prime = q + qr(iq)

      return
      end

c=======================================================================

      subroutine quatset( q, w, x, y, z )

      implicit none
      double precision q(4), w, x, y, z

      q(1) = w
      q(2) = x
      q(3) = y
      q(4) = z

      call quatnorm4( q )

      return
      end

c=======================================================================

      subroutine quatset2( q, x, y )

      implicit none
      double precision q(2), x, y

      q(1) = x
      q(2) = y

      call quatnorm2( q )

      return
      end

c=======================================================================

      subroutine quatcopy4( q1, q2 )

      implicit none
      double precision q1(4), q2(4)

      q2(1) = q1(1)
      q2(2) = q1(2)
      q2(3) = q1(3)
      q2(4) = q1(4)

      return
      end

c=======================================================================

      subroutine quatcopy2( q1, q2 )

      implicit none
      double precision q1(2), q2(2)

      q2(1) = q1(1)
      q2(2) = q1(2)

      return
      end

c=======================================================================

      subroutine quatswap( q1, q2 )

      implicit none
      double precision q1(4), q2(4)
      double precision qtmp

      qtmp = q1(1)
      q1(1) = q2(1)
      q2(1) = qtmp

      qtmp = q1(2)
      q1(2) = q2(2)
      q2(2) = qtmp

      qtmp = q1(3)
      q1(3) = q2(3)
      q2(3) = qtmp

      qtmp = q1(4)
      q1(4) = q2(4)
      q2(4) = qtmp

      return
      end
      
c=======================================================================

      subroutine quatswap2( q1, q2 )

      implicit none
      double precision q1(2), q2(2)
      double precision qtmp

      qtmp = q1(1)
      q1(1) = q2(1)
      q2(1) = qtmp

      qtmp = q1(2)
      q1(2) = q2(2)
      q2(2) = qtmp

      return
      end
      
c=======================================================================

      subroutine quatadd( q1, q2, q )

      implicit none
      double precision q1(4), q2(4)
      double precision q(4)

      q(1) = q1(1) + q2(1)
      q(2) = q1(2) + q2(2)
      q(3) = q1(3) + q2(3)
      q(4) = q1(4) + q2(4)

      return
      end

c=======================================================================

      subroutine quatadd2( q1, q2, q )

      implicit none
      double precision q1(2), q2(2)
      double precision q(2)

      q(1) = q1(1) + q2(1)
      q(2) = q1(2) + q2(2)

      return
      end

c=======================================================================

      subroutine quatsub4( q1, q2, q )

      implicit none
      double precision q1(4), q2(4)
      double precision q(4)

      q(1) = q1(1) - q2(1)
      q(2) = q1(2) - q2(2)
      q(3) = q1(3) - q2(3)
      q(4) = q1(4) - q2(4)

      return
      end

c=======================================================================

      subroutine quatsub2( q1, q2, q )

      implicit none
      double precision q1(2), q2(2)
      double precision q(2)

      q(1) = q1(1) - q2(1)
      q(2) = q1(2) - q2(2)

      return
      end

c=======================================================================

      subroutine quatmult4( q1, q2, q )

      implicit none
      double precision q1(4), q2(4)
      double precision q(4)

      q(1) = q1(1) * q2(1)
     &   - q1(2) * q2(2)
     &   - q1(3) * q2(3)
     &   - q1(4) * q2(4)

      q(2) = q1(1) * q2(2)
     &   + q1(2) * q2(1)
     &   + q1(3) * q2(4)
     &   - q1(4) * q2(3)

      q(3) = q1(1) * q2(3)
     &   + q1(3) * q2(1)
     &   + q1(4) * q2(2)
     &   - q1(2) * q2(4)

      q(4) = q1(1) * q2(4)
     &   + q1(4) * q2(1)
     &   + q1(2) * q2(3)
     &   - q1(3) * q2(2)

      return
      end

c=======================================================================

      subroutine quatmult2( q1, q2, q )

      implicit none
      double precision q1(2), q2(2)
      double precision q(2)

      q(1) = q1(1) * q2(1)
     &   - q1(2) * q2(2)

      q(2) = q1(1) * q2(2)
     &   + q1(2) * q2(1)

      return
      end

c=======================================================================

      subroutine quatconj( q1, q2 )

      implicit none
      double precision q1(4)
      double precision q2(4)

      q2(1) = q1(1)
      q2(2) = -q1(2)
      q2(3) = -q1(3)
      q2(4) = -q1(4)

      return
      end

c=======================================================================

      subroutine quatconj2( q1, q2 )

      implicit none
      double precision q1(2)
      double precision q2(2)

      q2(1) = q1(1)
      q2(2) = -q1(2)

      return
      end

c=======================================================================

      subroutine quatmagn4( q, m )

      implicit none
      double precision q(4)
      double precision m

      m = dsqrt( q(1)*q(1) + q(2)*q(2) + q(3)*q(3) + q(4)*q(4) )

      return
      end

c=======================================================================

      subroutine quatmagn2( q, m )

      implicit none
      double precision q(2)
      double precision m

      m = dsqrt( q(1)*q(1) + q(2)*q(2) )

      return
      end

c=======================================================================

      subroutine quatmagsq4( q, msq )

      implicit none
      double precision q(4)
      double precision msq

      msq = q(1)*q(1) + q(2)*q(2) + q(3)*q(3) + q(4)*q(4)

      return
      end

c=======================================================================

      subroutine quatmagsq2( q, msq )

      implicit none
      double precision q(2)
      double precision msq

      msq = q(1)*q(1) + q(2)*q(2)

      return
      end

c=======================================================================

      subroutine quatmaginv4( q, minv )

      implicit none
      double precision q(4)
      double precision minv
      double precision m

      minv = 0.d0

      call quatmagn4( q, m )
      if ( m .lt. 1.d-15 ) return
      minv = 1.d0 / m

c      call quatmagsq( q, msq )
c      if ( msq .lt. 1.d-30 ) return
c        minv = msq**(-0.5)
c      minv = dexp( -0.5d0 * dlog( msq ) )

      return
      end

c=======================================================================

      subroutine quatmaginv2( q, minv )

      implicit none
      double precision q(2)
      double precision minv
      double precision m

      minv = 0.d0

      call quatmagn2( q, m )
      if ( m .lt. 1.d-15 ) return
      minv = 1.d0 / m

c      call quatmagsq2( q, msq )
c      if ( msq .lt. 1.d-30 ) return
c        minv = msq**(-0.5)
c      minv = dexp( -0.5d0 * dlog( msq ) )

      return
      end

c=======================================================================

      subroutine quatnorm( q, qlen )

      implicit none
      integer qlen
      double precision q(qlen)

      if ( qlen == 4 ) then

         call quatnorm4( q )

      else if ( qlen == 2 ) then

         call quatnorm2( q )

      else if ( qlen /= 1 ) then

         print *, "Error in quatnorm"
         stop

      endif

      return
      end

c=======================================================================

      subroutine quatnorm4( q )

      implicit none
      double precision q(4)
      double precision minv

      call quatmaginv4( q, minv )

      q(1) = q(1) * minv
      q(2) = q(2) * minv
      q(3) = q(3) * minv
      q(4) = q(4) * minv

      return
      end

c=======================================================================

      subroutine quatnorm2( q )

      implicit none
      double precision q(2)
      double precision minv

      call quatmaginv2( q, minv )

      q(1) = q(1) * minv
      q(2) = q(2) * minv

      return
      end

c=======================================================================

c Make the value of q, which is periodic in TWO_PI, be between -PI and PI

      subroutine quatnorm1( q )

      implicit none
      double precision q

      logical first_pass
      data first_pass / .true. /
      double precision PI, TWO_PI
      save PI, TWO_PI, first_pass

      if ( first_pass ) then
         PI = dacos( -1.0d0 )
         TWO_PI = 2.0d0 * PI

         first_pass = .false.
      endif

      if ( q >= 0 ) then
         q = mod( q + PI, TWO_PI ) - PI
      else
         q = mod( q - PI, TWO_PI ) + PI
      endif

      end

c=======================================================================

      subroutine quatangl( d, theta )

      implicit none
      double precision d, theta

      if ( d > 2.0 ) then
         print *, 'd = ', d
      else
         theta = 4.0d0 * dasin( 0.5d0 * d )
      endif

      return
      end

c=======================================================================

      subroutine quatangl2( d, theta )

      implicit none
      double precision d, theta

      if ( d > 2.0 ) then
         print *, 'd = ', d
      else
         theta = 2.0d0 * dasin( 0.5d0 * d )
      endif

      return
      end

c=======================================================================

c        Note that [[d]] is the magnitude of the difference between a
c        _normalized_ q1 and a _normalized_ q2_prime.  This allows the
c        computation of an "angle" between the two.

c        This now works on the square of the magnitude of the difference
c        to avoid an extra sqrt call.

      subroutine quatdiffsq4( q1, q2, dsq )

      implicit none
      double precision q1(4), q2(4)
      double precision q1tmp(4)
      double precision q2tmp(4)
      double precision q_diff(4)
      double precision dsq
      integer n

      do n = 1, 4
         q2tmp(n) = q2(n)
         q1tmp(n) = q1(n)
      enddo
      call quatnorm4( q1tmp )
      call quatnorm4( q2tmp )
      call quatsub4( q2tmp, q1tmp, q_diff )
      call quatmagsq4( q_diff, dsq )
c      call quatangl( d, theta )

      return
      end

c=======================================================================

c        Note that [[d]] is the magnitude of the difference between a
c        _normalized_ q1 and a _normalized_ q2_prime.  This allows the
c        computation of an "angle" between the two.

c        This now works on the square of the magnitude of the difference
c        to avoid an extra sqrt call.

      subroutine quatdiffsq2( q1, q2, dsq )

      implicit none
      double precision q1(2), q2(2)
      double precision q1tmp(2)
      double precision q2tmp(2)
      double precision q_diff(2)
      double precision dsq
      integer n

      do n = 1, 2
         q2tmp(n) = q2(n)
         q1tmp(n) = q1(n)
      enddo
      call quatnorm2( q1tmp )
      call quatnorm2( q2tmp )
      call quatsub2( q2tmp, q1tmp, q_diff )
      call quatmagsq2( q_diff, dsq )
c      call quatangl2( d, theta )

      return
      end

c=======================================================================

      subroutine quatdiffabs1( q1, q2, q_diff )

      implicit none
      double precision q1, q2
      double precision q_diff

      q_diff = dabs( q2 - q1 )

      return
      end

c=======================================================================

c        Same as quatdiffsq above, but rotating q2 first

      subroutine quatrotatediffsq4( q1, q2, qr, q2_prime, dsq )

      implicit none
      double precision q1(4), q2(4), qr(4), q2_prime(4)
      double precision q1tmp(4)
      double precision q2tmp(4)
      double precision q_diff(4)
      double precision dsq
      integer n

      call quatmult4( q2, qr, q2_prime )
      do n = 1, 4
         q2tmp(n) = q2_prime(n)
         q1tmp(n) = q1(n)
      enddo
      call quatnorm4( q1tmp )
      call quatnorm4( q2tmp )
      call quatsub4( q2tmp, q1tmp, q_diff )
      call quatmagsq4( q_diff, dsq )
c      call quatangl( d, theta )

      return
      end

c=======================================================================

c        Same as quatdiffsq above, but rotating q2 first

      subroutine quatrotatediffsq2( q1, q2, qr, q2_prime, dsq )

      implicit none
      double precision q1(2), q2(2), qr(2), q2_prime(2)
      double precision q1tmp(2)
      double precision q2tmp(2)
      double precision q_diff(2)
      double precision dsq
      integer n

      call quatmult2( q2, qr, q2_prime )
      do n = 1, 2
         q2tmp(n) = q2_prime(n)
         q1tmp(n) = q1(n)
      enddo
      call quatnorm2( q1tmp )
      call quatnorm2( q2tmp )
      call quatsub2( q2tmp, q1tmp, q_diff )
      call quatmagsq2( q_diff, dsq )
c      call quatangl2( d, theta )

      return
      end

c=======================================================================

c        Same as quatdiffabs1 above, but rotating q2 first

      subroutine quatrotatediffabs1( q1, q2, qr, q2_prime, qdabs )

      implicit none
      double precision q1, q2, qr, q2_prime
      double precision qdabs

      q2_prime = q2 + qr
      qdabs = dabs( q2_prime - q1 )

      return
      end

c=======================================================================

c        compute rotation matrix corresponding to quaternion q

      subroutine quatrotationmatrix( q, matrixR )

      implicit none
      double precision q(4)
      double precision matrixR(3,3)

c first column
      matrixR(1,1)=q(1)*q(1)+q(2)*q(2)-q(3)*q(3)-q(4)*q(4)
      matrixR(2,1)=2.d0*(q(1)*q(4)+q(2)*q(3))
      matrixR(3,1)=2.d0*(q(2)*q(4)-q(1)*q(3))

c 2nd column
      matrixR(1,2)=2.d0*(q(2)*q(3)-q(1)*q(4))
      matrixR(2,2)=q(1)*q(1)-q(2)*q(2)+q(3)*q(3)-q(4)*q(4)
      matrixR(3,2)=2.d0*(q(1)*q(2)+q(3)*q(4))

c 3rd column
      matrixR(1,3)=2.d0*(q(1)*q(3)+q(2)*q(4))
      matrixR(2,3)=2.d0*(q(3)*q(4)-q(1)*q(2))
      matrixR(3,3)=q(1)*q(1)-q(2)*q(2)-q(3)*q(3)+q(4)*q(4)

      return
      end

c=======================================================================

c        compute derivative rotation matrix corresponding to quaternion q
c        with respect to q1

      subroutine derivquatrotationmatrixdq( q, matrixdRdq )

      implicit none
      double precision q(4)
      double precision matrixdRdq(3,3,4)

c d/dq1
      matrixdRdq(1,1,1)= 2.d0*q(1)
      matrixdRdq(2,1,1)= 2.d0*q(4)
      matrixdRdq(3,1,1)=-2.d0*q(3)

      matrixdRdq(1,2,1)=-2.d0*q(4)
      matrixdRdq(2,2,1)= 2.d0*q(1)
      matrixdRdq(3,2,1)= 2.d0*q(2)

      matrixdRdq(1,3,1)= 2.d0*q(3)
      matrixdRdq(2,3,1)=-2.d0*q(2)
      matrixdRdq(3,3,1)= 2.d0*q(1)

c d/dq2
      matrixdRdq(1,1,2)= 2.d0*q(2)
      matrixdRdq(2,1,2)= 2.d0*q(3)
      matrixdRdq(3,1,2)= 2.d0*q(4)

      matrixdRdq(1,2,2)= 2.d0*q(3)
      matrixdRdq(2,2,2)=-2.d0*q(2)
      matrixdRdq(3,2,2)= 2.d0*q(2)

      matrixdRdq(1,3,2)= 2.d0*q(4)
      matrixdRdq(2,3,2)=-2.d0*q(1)
      matrixdRdq(3,3,2)=-2.d0*q(2)

c d/dq3
      matrixdRdq(1,1,3)=-2.d0*q(3)
      matrixdRdq(2,1,3)= 2.d0*q(2)
      matrixdRdq(3,1,3)= 2.d0*q(1)

      matrixdRdq(1,2,3)= 2.d0*q(2)
      matrixdRdq(2,2,3)= 2.d0*q(3)
      matrixdRdq(3,2,3)= 2.d0*q(4)

      matrixdRdq(1,3,3)= 2.d0*q(1)
      matrixdRdq(2,3,3)= 2.d0*q(4)
      matrixdRdq(3,3,3)=-2.d0*q(3)

c d/dq4
      matrixdRdq(1,1,4)=-2.d0*q(4)
      matrixdRdq(2,1,4)= 2.d0*q(1)
      matrixdRdq(3,1,4)= 2.d0*q(2)

      matrixdRdq(1,2,4)=-2.d0*q(1)
      matrixdRdq(2,2,4)=-2.d0*q(4)
      matrixdRdq(3,2,4)= 2.d0*q(3)
      
      matrixdRdq(1,3,4)= 2.d0*q(2)
      matrixdRdq(2,3,4)= 2.d0*q(3)
      matrixdRdq(3,3,4)= 2.d0*q(4)

      return
      end


c=======================================================================

c        compute rotation matrix corresponding to complex q

      subroutine complexrotationmatrix( q, matrixR )

      implicit none
      double precision q(2)
      double precision matrixR(3,3)
      
c first column
      matrixR(1,1)= q(1)
      matrixR(2,1)= q(2)
      matrixR(3,1)= 0.d0

c 2nd column
      matrixR(1,2)=-q(2)
      matrixR(2,2)= q(1)
      matrixR(3,2)= 0.d0

c 3rd column
      matrixR(1,3)= 0.d0
      matrixR(2,3)= 0.d0
      matrixR(3,3)= 1.d0

      return
      end

c=======================================================================

c compute derivative rotation matrix corresponding to complex q

      subroutine derivcomplexrotationmatrixdq( q, matrixdRdq )

      implicit none
      double precision q(2)
      double precision matrixdRdq(3,3,2)

c d/dq1
      matrixdRdq(1,1,1)= 1.d0
      matrixdRdq(2,1,1)= 0.d0
      matrixdRdq(3,1,1)= 0.d0

      matrixdRdq(1,2,1)= 0.d0
      matrixdRdq(2,2,1)= 1.d0
      matrixdRdq(3,2,1)= 0.d0

      matrixdRdq(1,3,1)= 0.d0
      matrixdRdq(2,3,1)= 0.d0
      matrixdRdq(3,3,1)= 0.d0

c d/dq2
      matrixdRdq(1,1,2)= 0.d0
      matrixdRdq(2,1,2)= 1.d0
      matrixdRdq(3,1,2)= 0.d0

      matrixdRdq(1,2,2)=-1.d0
      matrixdRdq(2,2,2)= 0.d0
      matrixdRdq(3,2,2)= 0.d0

      matrixdRdq(1,3,2)= 0.d0
      matrixdRdq(2,3,2)= 0.d0
      matrixdRdq(3,3,2)= 0.d0

      return
      end

c=======================================================================

c        compute rotation tensor for cubic system
c        Aijkl=ai*aj*ak*al+...
c        where (a1,a2,a3)^T is a column of the rotation matrix

      function rotationtensorijkl( i, j, k, l, matrixR )

      implicit none
      double precision rotationtensorijkl
      integer          i,j,k,l
      double precision matrixR(3,3)

      rotationtensorijkl =
     &     matrixR(i,1)*matrixR(j,1)*matrixR(k,1)*matrixR(l,1)
     &    +matrixR(i,2)*matrixR(j,2)*matrixR(k,2)*matrixR(l,2)
     &    +matrixR(i,3)*matrixR(j,3)*matrixR(k,3)*matrixR(l,3)

      return
      end

c=======================================================================
c Valid values for floor_type are "max", "tanh" and "sqrt"

      function eval_grad_normi( grad_norm2, floor_type, 
     &                          floor_grad_norm2, 
     &                          max_grad_normi)

      implicit none
      double precision eval_grad_normi
      double precision grad_norm2, max_grad_normi
      character*(*) floor_type

      double precision gng2, tol_taylor2, floor_grad_norm2, grad_norm

      parameter(tol_taylor2=0.01d0)

      if ( floor_type(1:1) .eq. 'm' ) then
         if (grad_norm2 .gt. floor_grad_norm2) then
            eval_grad_normi = grad_norm2 ** (-0.5d0)
         else
            eval_grad_normi = max_grad_normi
         endif
      else if ( floor_type(1:1) .eq. 't' )then
            gng2 = grad_norm2 * max_grad_normi * max_grad_normi
            if (gng2 .gt. tol_taylor2) then
               grad_norm = dsqrt(grad_norm2)
               eval_grad_normi = tanh(max_grad_normi * grad_norm) 
     &                        / grad_norm
            else
               eval_grad_normi = max_grad_normi * 
     &              (1.d0 - gng2 * (5.d0 - 2.d0 * gng2) / 15.d0)
            endif
      else if ( floor_type(1:1) .eq. 's' )then
            eval_grad_normi = ( grad_norm2+floor_grad_norm2 )** (-0.5d0)
      else

         print *, "Error in eval_grad_normi: floor_type unknown"
         stop

      endif
      

      return
      end
