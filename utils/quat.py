# Copyright (c) 2018, Lawrence Livermore National Security, LLC and
# UT-Battelle, LLC.
# Produced at the Lawrence Livermore National Laboratory and
# the Oak Ridge National Laboratory
# LLNL-CODE-747500
# All rights reserved.
# This file is part of AMPE. 
# For details, see https://github.com/LLNL/AMPE
# Please also read AMPE/LICENSE.
# 
import numpy
import math
import random

def makeQuat( w=0, x=0, y=0, z=0 ) :
  q = numpy.array( [ w, x, y, z ], float )
  return q

def makeQuat2( x=0, y=0 ) :
  q = numpy.array( [ x, y ], float )
  return q

#=======================================================================

def makeNormalizedQuat( w=1, x=0, y=0, z=0 ) :
  q = numpy.array( [ w, x, y, z ], float )
  normalizeQuat( q )
  return q

def makeNormalizedQuat2( x=1, y=0 ) :
  q = numpy.array( [ x, y ], float )
  normalizeQuat( q )
  return q

#=======================================================================

def copyQuat( q1, q2 ) :
  assert( len( q1 ) == len( q2 ) )
  qlen = len( q1 )
  for j in range( qlen ) :
    q2[j] = q1[j]
  return q2

#=======================================================================

def duplicateQuat( q ) :
  qlen = len( q )
  assert( qlen == 4 or qlen == 2 )
  if ( qlen == 4 ) :
    q_copy = numpy.array( [ 1, 0, 0, 0 ], float )
  elif ( qlen == 2 ) :
    q_copy = numpy.array( [ 1, 0 ], float )

  for j in range( len( q ) ) :
    q_copy[j] = q[j]
  return q_copy

#=======================================================================

def addQuat( q1, q2, q ) :
  assert( len( q1 ) == len( q2 ) )
  assert( len( q ) == len( q1 ) )
  qlen = len( q )
  for j in range( qlen ) :
    q[j] = q1[j] + q2[j]

#=======================================================================

def subtractQuat( q1, q2, q ) :
  assert( len( q1 ) == len( q2 ) )
  assert( len( q ) == len( q1 ) )
  qlen = len( q )
  for j in range( qlen ) :
    q[j] = q1[j] - q2[j]

#=======================================================================

def multiplyQuat( q1, q2, q=0 ) :
  assert( len( q1 ) == len( q2 ) )
  qlen = len( q1 )
  assert( qlen == 2 or qlen == 4 )
  if ( len(q)>1 ) :
    assert( len( q ) == qlen )
    q_tmp = duplicateQuat( q )
  else :
    if ( qlen == 4 ) :
      q_tmp = makeQuat( 0, 0, 0, 0 )
    else :
      q_tmp = makeQuat2( 0, 0 )

  if ( qlen == 4 ) :
    q_tmp[0] = q1[0] * q2[0] \
               - q1[1] * q2[1] \
               - q1[2] * q2[2] \
               - q1[3] * q2[3]

    q_tmp[1] = q1[0] * q2[1] \
               + q1[1] * q2[0] \
               + q1[2] * q2[3] \
               - q1[3] * q2[2]

    q_tmp[2] = q1[0] * q2[2] \
               + q1[2] * q2[0] \
               + q1[3] * q2[1] \
               - q1[1] * q2[3]

    q_tmp[3] = q1[0] * q2[3] \
               + q1[3] * q2[0] \
               + q1[1] * q2[2] \
               - q1[2] * q2[1]
  
  else :
    q_tmp[0] = q1[0] * q2[0] \
               - q1[1] * q2[1]

    q_tmp[1] = q1[0] * q2[1] \
               + q1[1] * q2[0]

  if ( len(q)>1 ) :
    copyQuat( q_tmp, q )
  else :
    return q_tmp
  
#=======================================================================

def rotateQuat( q, q_rotation ) :
  q_prime = duplicateQuat( q )
  multiplyQuat( q, q_rotation, q_prime )
  return q_prime

#=======================================================================

def conjugateQuat( q1, q2=0 ) :
  qlen = len( q1 )

  ret_q = 0
  if ( q2 == 0 ) :
    ret_q = 1
    q2 = duplicateQuat( q )
  else :
    assert( len( q2 ) == qlen )
    
  q2[0] = q1[0]
  for j in range( 1, qlen ) :
    q2[j] = -q1[j]

  if ( ret_q ) :
    return q2

#=======================================================================

def quatMagnitude( q ) :
  qlen = len( q )

  sumsq = 0.0
  for j in range( qlen ) :
    sumsq += q[j] * q[j]

  return math.sqrt( sumsq )

#=======================================================================

def normalizeQuat( q ) :
  qlen = len( q )

  m = quatMagnitude( q )

  if ( m < 1.e-15 ) : return

  for j in range( qlen ) :
    q[j] = q[j] / m

#=======================================================================

def quatAngle( d ) :
  if ( d > 2.0 ) :
    print(  'd =', d )
    return -1

  else :
    return 4.0 * math.asin( 0.5 * d )

def quatAngle2( d ) :
  if ( d > 2.0 ) :
    print( 'd =', d )
    return -1

  else :
    return 2.0 * math.asin( 0.5 * d )

#=======================================================================

#        Note that [[d]] is the magnitude of the difference between a
#        _normalized_ q1 and a _normalized_ q2_prime.  This allows the
#        computation of an "angle" between the two.

def quatDiff( q1, q2, qr ) :
  assert( len( q1 ) == len( q2 ) )
  qlen = len( q1 )
  assert( qlen == 2 or qlen == 4 )

  q1tmp = duplicateQuat( q1 )
  normalizeQuat( q1tmp )

  if qlen == 2:
    q_diff = makeQuat2()
  else:
    q_diff = makeQuat()

  q2_prime = rotateQuat( q2, qr )
  q2tmp = duplicateQuat( q2_prime )
  normalizeQuat( q2tmp )

  subtractQuat( q2tmp, q1tmp, q_diff )
  d = quatMagnitude( q_diff )
  # theta = quatAngle( d )

  return d, q2_prime

def quatDiff1( q1, q2, qr ) :
  q1tmp = q1

  q2_prime = q2 + qr
  q_diff = math.fabs( q2_prime - q1 )

  return q_diff, q2_prime

#=======================================================================

qr4 = numpy.zeros( ( 48, 4 ), float )

qr4[ 0,:] = makeNormalizedQuat(  1,  0,  0,  0 )
qr4[ 1,:] = makeNormalizedQuat(  0,  1,  0,  0 )
qr4[ 2,:] = makeNormalizedQuat(  0,  0,  1,  0 )
qr4[ 3,:] = makeNormalizedQuat(  0,  0,  0,  1 )
qr4[ 4,:] = makeNormalizedQuat( -1,  0,  0,  0 )
qr4[ 5,:] = makeNormalizedQuat(  0, -1,  0,  0 )
qr4[ 6,:] = makeNormalizedQuat(  0,  0, -1,  0 )
qr4[ 7,:] = makeNormalizedQuat(  0,  0,  0, -1 )

qr4[ 8,:] = makeNormalizedQuat(  1,  1,  0,  0 )
qr4[ 9,:] = makeNormalizedQuat(  1,  0,  1,  0 )
qr4[10,:] = makeNormalizedQuat(  1,  0,  0,  1 )
qr4[11,:] = makeNormalizedQuat(  0,  1,  1,  0 )
qr4[12,:] = makeNormalizedQuat(  0,  1,  0,  1 )
qr4[13,:] = makeNormalizedQuat(  0,  0,  1,  1 )
qr4[14,:] = makeNormalizedQuat( -1,  1,  0,  0 )
qr4[15,:] = makeNormalizedQuat( -1,  0,  1,  0 )
qr4[16,:] = makeNormalizedQuat( -1,  0,  0,  1 )
qr4[17,:] = makeNormalizedQuat(  0, -1,  1,  0 )
qr4[18,:] = makeNormalizedQuat(  0, -1,  0,  1 )
qr4[19,:] = makeNormalizedQuat(  0,  0, -1,  1 )
qr4[20,:] = makeNormalizedQuat(  1, -1,  0,  0 )
qr4[21,:] = makeNormalizedQuat(  1,  0, -1,  0 )
qr4[22,:] = makeNormalizedQuat(  1,  0,  0, -1 )
qr4[23,:] = makeNormalizedQuat(  0,  1, -1,  0 )
qr4[24,:] = makeNormalizedQuat(  0,  1,  0, -1 )
qr4[25,:] = makeNormalizedQuat(  0,  0,  1, -1 )
qr4[26,:] = makeNormalizedQuat( -1, -1,  0,  0 )
qr4[27,:] = makeNormalizedQuat( -1,  0, -1,  0 )
qr4[28,:] = makeNormalizedQuat( -1,  0,  0, -1 )
qr4[29,:] = makeNormalizedQuat(  0, -1, -1,  0 )
qr4[30,:] = makeNormalizedQuat(  0, -1,  0, -1 )
qr4[31,:] = makeNormalizedQuat(  0,  0, -1, -1 )

qr4[32,:] = makeNormalizedQuat(  1,  1,  1,  1 )
qr4[33,:] = makeNormalizedQuat( -1,  1,  1,  1 )
qr4[34,:] = makeNormalizedQuat(  1, -1,  1,  1 )
qr4[35,:] = makeNormalizedQuat(  1,  1, -1,  1 )
qr4[36,:] = makeNormalizedQuat(  1,  1,  1, -1 )
qr4[37,:] = makeNormalizedQuat( -1, -1,  1,  1 )
qr4[38,:] = makeNormalizedQuat( -1,  1, -1,  1 )
qr4[39,:] = makeNormalizedQuat( -1,  1,  1, -1 )
qr4[40,:] = makeNormalizedQuat(  1, -1, -1,  1 )
qr4[41,:] = makeNormalizedQuat(  1, -1,  1, -1 )
qr4[42,:] = makeNormalizedQuat(  1,  1, -1, -1 )
qr4[43,:] = makeNormalizedQuat(  1, -1, -1, -1 )
qr4[44,:] = makeNormalizedQuat( -1,  1, -1, -1 )
qr4[45,:] = makeNormalizedQuat( -1, -1,  1, -1 )
qr4[46,:] = makeNormalizedQuat( -1, -1, -1,  1 )
qr4[47,:] = makeNormalizedQuat( -1, -1, -1, -1 )

iq_qr4_conj = numpy.zeros( 48, int )

iq_qr4_conj[0] = 0
iq_qr4_conj[1] = 5
iq_qr4_conj[2] = 6
iq_qr4_conj[3] = 7
iq_qr4_conj[4] = 4
iq_qr4_conj[5] = 1
iq_qr4_conj[6] = 2
iq_qr4_conj[7] = 3
iq_qr4_conj[8] = 20
iq_qr4_conj[9] = 21
iq_qr4_conj[10] = 22
iq_qr4_conj[11] = 29
iq_qr4_conj[12] = 30
iq_qr4_conj[13] = 31
iq_qr4_conj[14] = 26
iq_qr4_conj[15] = 27
iq_qr4_conj[16] = 28
iq_qr4_conj[17] = 23
iq_qr4_conj[18] = 24
iq_qr4_conj[19] = 25
iq_qr4_conj[20] = 8
iq_qr4_conj[21] = 9
iq_qr4_conj[22] = 10
iq_qr4_conj[23] = 17
iq_qr4_conj[24] = 18
iq_qr4_conj[25] = 19
iq_qr4_conj[26] = 14
iq_qr4_conj[27] = 15
iq_qr4_conj[28] = 16
iq_qr4_conj[29] = 11
iq_qr4_conj[30] = 12
iq_qr4_conj[31] = 13
iq_qr4_conj[32] = 43
iq_qr4_conj[33] = 47
iq_qr4_conj[34] = 42
iq_qr4_conj[35] = 41
iq_qr4_conj[36] = 40
iq_qr4_conj[37] = 44
iq_qr4_conj[38] = 45
iq_qr4_conj[39] = 46
iq_qr4_conj[40] = 38
iq_qr4_conj[41] = 35
iq_qr4_conj[42] = 34
iq_qr4_conj[43] = 32
iq_qr4_conj[44] = 37
iq_qr4_conj[45] = 38
iq_qr4_conj[46] = 39
iq_qr4_conj[47] = 33

#-----------------------------------------------------------------------

qr2 = numpy.zeros( ( 4, 2 ), float )

qr2[0,:] = makeNormalizedQuat2(  1,  0 )
qr2[1,:] = makeNormalizedQuat2(  0,  1 )
qr2[2,:] = makeNormalizedQuat2( -1,  0 )
qr2[3,:] = makeNormalizedQuat2(  0, -1 )

iq_qr2_conj = numpy.zeros( 4, int )

iq_qr2_conj[0] = 0
iq_qr2_conj[1] = 3
iq_qr2_conj[2] = 2
iq_qr2_conj[3] = 1

#-----------------------------------------------------------------------

qr1 = numpy.zeros( 9, float )

qr1[0] = 0.0
qr1[1] = 0.5 * math.pi
qr1[2] = -0.5 * math.pi
qr1[3] = math.pi
qr1[4] = -math.pi
qr1[5] = 1.5 * math.pi
qr1[6] = -1.5 * math.pi
qr1[7] = 2.0 * math.pi
qr1[8] = -2.0 * math.pi

iq_qr1_conj = numpy.zeros( 9, int )

iq_qr1_conj[0] = 0
iq_qr1_conj[1] = 2
iq_qr1_conj[2] = 1
iq_qr1_conj[3] = 4
iq_qr1_conj[4] = 3
iq_qr1_conj[5] = 6
iq_qr1_conj[6] = 5
iq_qr1_conj[7] = 8
iq_qr1_conj[8] = 7

#-----------------------------------------------------------------------

def getNumberSymmetryRotations( qlen ) :
  if ( qlen == 4 ) :
    return len( iq_qr4_conj )
  elif ( qlen == 2 ) :
    return len( iq_qr2_conj )
  else :
    return len( iq_qr1_conj )

def getQuatSymmetryRotation( n, qlen=4 ) :
  global qr1
  global qr2
  global qr4

  if ( qlen == 4 ) :
    return qr4[n,:]
  elif ( qlen == 2 ) :
    return qr2[n,:]
  else :
    return qr1[n]

def getQuatSymmetryRotationConjugateIdx( n, qlen=4 ) :
  global iq_qr1_conj
  global iq_qr2_conj
  global iq_qr4_conj

  if ( qlen == 4 ) :
    return iq_qr4_conj[n]
  elif ( qlen == 2 ) :
    return iq_qr2_conj[n]
  else :
    return iq_qr1_conj[n]

def getQuatSymmetryRotationConjugate( n, qlen=4 ) :
  global qr1
  global qr2
  global qr4
  global iq_qr1_conj
  global iq_qr2_conj
  global iq_qr4_conj

  if ( qlen == 4 ) :
    return qr4[iq_qr4_conj[n]]
  elif ( qlen == 2 ) :
    return qr2[iq_qr2_conj[n]]
  else :
    return qr1[iq_qr1_conj[n]]

#=======================================================================

def quatSymm( q1, q2, iq=0 ) :
  assert( len( q1 ) == len( q2 ) )
  qlen = len( q1 )
  assert( qlen == 2 or qlen == 4 )

  if ( qlen == 4 ) :
    global qr4
    global iq_qr4_conj
    qr = qr4
    iq_qr_conj = iq_qr4_conj
    CRIT_DIFF = 2.0 * math.sin( math.pi * 0.0625 )
    N_ROTATIONS = 48
  else :
    global qr2
    global iq_qr2_conj
    qr = qr2
    iq_qr_conj = iq_qr2_conj
    CRIT_DIFF = 2.0 * math.sin( math.pi / 8.0 )
    N_ROTATIONS = 4

  if ( iq > N_ROTATIONS or iq < -N_ROTATIONS ) : iq = 0

  if ( iq < 0 ) : iq = iq_qr_conj[-iq]

  # Test qr(:,iq) first
  diff, q2_prime = quatDiff( q1, q2, qr[iq,:] )

  min_diff = diff
  min_iq = iq
  min_q2_prime = duplicateQuat( q2_prime )
  if ( diff <= CRIT_DIFF ) : return min_q2_prime

  for nn in range( N_ROTATIONS ) :
    if ( nn == iq ) : continue

    diff, q2_prime = quatDiff( q1, q2, qr[nn,:] )

    if ( diff < min_diff ) :
      min_diff = diff
      min_iq = nn
      copyQuat( q2_prime, min_q2_prime )

    if ( diff <= CRIT_DIFF ) : break

  iq = min_iq
  return min_q2_prime

def quatSymm1( q1, q2, iq=0 ) :
  global qr1
  global iq_qr1_conj

  CRIT_DIFF = math.pi / 4.0
  N_ROTATIONS = 9

  if ( iq > N_ROTATIONS or iq < -N_ROTATIONS ) : iq = 0

  if ( iq < 0 ) : iq = iq_qr1_conj[-iq]

  # Test qr(:,iq) first
  diff, q2_prime = quatDiff1( q1, q2, qr1[iq] )
  if ( diff <= CRIT_DIFF ) : return

  min_diff = diff
  min_iq = iq
  min_q2_prime = q2_prime

  for nn in range( N_ROTATIONS ) :
    if ( nn == iq ) : continue

    diff, q2_prime = quatDiff1( q1, q2, qr1[nn] )

    if ( diff < min_diff ) :
      min_diff = diff
      min_iq = nn
      min_q2_prime = q2_prime

    if ( diff <= CRIT_DIFF ) : break

  iq = min_iq
  return min_q2_prime

#=======================================================================

def quatMakeFundamental( q ) :
  qref = makeNormalizedQuat() # (1,0,0,0)
  qsym = quatSymm( qref, q )
  return qsym

#=======================================================================

def quatToColor( q1,q2,q3,q4 ) :
  qq = makeQuat(q1,q2,q3,q4)
  qq = quatMakeFundamental( qq )
  return color4(128+127*qq[1],128+127*qq[2],128+127*qq[3],191+64*qq[0])

def getQuatRandom(t, qlen=4):
  if ( qlen == 4 ) :
    h1=0.1
    h2=0.5*h1
    h3=h2
    #K. Shoemake. "Uniform random rotations."
    #In D. Kirk, editor, Graphics Gems III, pages 124-132. Academic, New York, 1992.
    u1 = random.uniform(0, 1.)
    u2 = random.uniform(0, 1.)
    u3 = random.uniform(0, 1.)
    #print( u1,u2,u3

    #restricts possible orientations
    n  = math.floor( u1/h1 )
    u1 = n*h1
    n  = math.floor( u2/h2 )
    u2 = n*h2
    n  = math.floor( u3/h3 )
    u3 = n*h3

    w = math.sqrt(1.-u1)*math.sin(2.*math.pi*u2)
    x = math.sqrt(1.-u1)*math.cos(2.*math.pi*u2)
    y = math.sqrt(u1)*math.sin(2.*math.pi*u3)
    z = math.sqrt(u1)*math.cos(2.*math.pi*u3)

  else : #qlen=2
    w = math.cos( t )
    x = math.sin( t)
    y = 0.
    z = 0.

  return [w,x,y,z]
