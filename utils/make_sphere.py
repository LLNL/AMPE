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
# standard packages
import math
import random
import sys
import string
from optparse import OptionParser

# my required packages
import quat as Q

# other required packages
import numpy as N
import netCDF4 as nc4

#-----------------------------------------------------------------------
# command-line arguments and options

usage = "Usage: %prog [options] filename"

parser = OptionParser( usage = usage )

parser.add_option( "-1", "--one", action="store_true",
                   help="one sphere centered in the middle of the domain [default]" )
parser.add_option( "-2", "--two", action="store_true",
                   default=False,
                   help="two spheres centered in opposite quadrants/octants" )
parser.add_option( "--ratio", type="float", default=1.5,
                   help="ratio of sphere 1 to sphere 2 [default: %default]" )
parser.add_option( "-n", "--ncells", type="int",
                   help="number of cells in x,y,z (will all be equal)" )
parser.add_option( "-x", "--nx", type="int",
                   help="number of cells in x direction" )
parser.add_option( "-y", "--ny", type="int",
                   help="number of cells in y direction" )
parser.add_option( "-z", "--nz", type="int",
                   help="number of cells in z direction" )
parser.add_option( "--cx", type="int",
                   help="shift center by cx cells in x direction" )
parser.add_option( "--cy", type="int",
                   help="shift center by cy cells in y direction" )
parser.add_option( "--cz", type="int",
                   help="shift center by cz cells in z direction" )
parser.add_option( "-r", "--radius", type="int",
                   help="number of cells for sphere radius [default: nx/4]" )
parser.add_option( "-a", "--angle", type="float",
                   help="angular difference between interior and exterior" )
parser.add_option( "-i", "--phase-in", type="float", default=1.0,
                   help="phase in interior region [default: %default]" )
parser.add_option( "-o", "--phase-out", type="float", default=0.0,
                   help="phase in exterior region [default: %default]" )
parser.add_option( "--qlen", type="int", default=4,
                   help="Length of orientation vector" )
parser.add_option( "--quat-in", type="string", # something like "1,0,0,0",
                   help="quaternion in interior region" )
parser.add_option( "--quat-in-two", type="string", # something like "1,0,0,0",
                   help="quaternion in interior region of sphere 2" )
parser.add_option( "--quat-out", type="string",
                   help="quaternion in exterior region" )
parser.add_option( "--angle-in", type="float",
                   help="grain angle in interior region" )
parser.add_option( "--angle-in-two", type="float",
                   help="grain angle in interior region of sphere 2" )
parser.add_option( "--angle-out", type="float",
                   help="grain angle in exterior region" )
parser.add_option( "-c", "--concentration-in", type="float", default=0.1,
                   help="concentration in interior region [default: %default]" )
parser.add_option( "--concentration-out", type="float", default=0.06,
                   help="concentration in exterior region [default: %default]" )
parser.add_option( "-s", "--symmetry-test", action="store_true",
                   default=False,
                   help="rotate each octant by a random symmetry rotation" )
parser.add_option( "--random-symmetry", action="store_true",
                   default=False,
                   help="rotate each *CELL* by a random symmetry rotation" )
parser.add_option( "--random-spheres", type="int",
                   help="Use n random spheres" )

(options, args) = parser.parse_args()

if ( len( args ) < 1 ) :
  parser.error( "filename required argument missing" )
filename = args[0]

QLEN = 4
if ( not options.qlen is None ) : QLEN = options.qlen
if ( QLEN != 4 and QLEN != 2 and QLEN != 1 ) :
  print("Error: valid values of qlen are 1, 2, and 4")
  sys.exit(1)

n_spheres = "one"
if ( options.two ) :
  if ( options.one ) :
    print( "Error: cannot use options one and two together" )
    sys.exit(1)
  if ( not options.random_spheres is None ) :
    print( "Error: cannot use options two and random-spheres together" )
    sys.exit(1)
  n_spheres = "two"

if ( not options.random_spheres is None ) :
  if ( options.one ) :
    print( "Error: cannot use options one and random-spheres together" )
    sys.exit(1)
  n_spheres = "random"

nx = options.ncells
ny = options.ncells
nz = options.ncells

if ( not options.nx is None ) : nx = options.nx
if ( not options.ny is None ) : ny = options.ny
if ( not options.nz is None ) : nz = options.nz

if ( not ( nx and ny and nz ) ) :
  print( "Error: either -n or all of -nx -ny -nz are required" )
  sys.exit(1)

radii = []
cx = []
cy = []
cz = []
radius = options.radius
if radius is None :
  if ( n_spheres == "one" ) :
    radius = nx / 4
    print( "Radius not specified: using nx/4 = %d" % radius )
  else :
    radius = nx / 8
    print( "Radius not specified: using nx/8 = %d" % radius )

if ( n_spheres == "one" ) :
  print( "One sphere of radius %d" % radius)
  radii.append( radius )
  if( not options.cx is None ) :
    cx.append(options.cx)
  else :
    cx.append( nx / 2 )
  if( not options.cy is None ) :
    cy.append(options.cy)
  else :
    cy.append( ny / 2 )
  if( not options.cz is None ) :
    cz.append(options.cz)
  else:
    cz.append( nz / 2 )
  print( "Center: ", cx, cy, cz)
else :
  radii.append( radius )
  cx.append( nx / 4 )
  cx.append( 3 * nx / 4 )
  cy.append( ny / 4 )
  cy.append( 3 * ny / 4 )
  cz.append( nz / 2 )
  cz.append( nz / 2 )
  ratio = options.ratio
  radius_two = radius * ratio
  radii.append( radius_two )

phase_inside = options.phase_in
phase_outside = options.phase_out
conc_inside = options.concentration_in
conc_outside = options.concentration_out

use_simple_rotation = True

quat_inside = []
quat_outside = None

if ( not options.quat_in is None ) :
  if ( QLEN != 4 and QLEN != 2 ) :
    print( "Error: quat-in option is for QLEN=4 or 2" )
    sys.exit(1)

  use_simple_rotation = False

  q = map( float, options.quat_in.split( ',' ) )
  if ( QLEN == 4 ) :
    quat_inside.append( Q.makeNormalizedQuat( q[0], q[1], q[2], q[3] ) )
  else :
    quat_inside.append( Q.makeNormalizedQuat2( q[0], q[1] ) )
  
  if ( n_spheres == "one" and options.quat_out is None ) :
    print( "Error: must specify quat-in and quat-out together" )
    sys.exit(1)

  if ( n_spheres == "two" and
       ( options.quat_in_two is None or
         options.quat_out is None ) ) :
    print( "Error: must specify quat-in, quat-in-two, " \
          "and quat-out together" )
    sys.exit(1)

if ( not options.quat_in_two is None ) :
  if ( QLEN != 4 and QLEN != 2 ) :
    print( "Error: quat-in-two option is for QLEN=4 or 2" )
    sys.exit(1)

  use_simple_rotation = False

  q = map( float, options.quat_in_two.split(',') )
  if ( QLEN == 4 ) :
    quat_inside.append( Q.makeNormalizedQuat( q[0], q[1], q[2], q[3] ) )
  else :
    quat_inside.append( Q.makeNormalizedQuat2( q[0], q[1] ) )
  
  if ( n_spheres == "one" ) :
    print( "Warning: ignoring quat-in-two with only one sphere" )
    sys.exit(1)

  if ( n_spheres == "two" and
       ( options.quat_in is None or
         options.quat_out is None ) ) :
    print( "Error: must specify quat-in, quat-in-two, " \
          "and quat-out together" )
    sys.exit(1)
    
if ( not options.quat_out is None ) :
  if ( QLEN != 4 and QLEN != 2 ) :
    print( "Error: quat-out option is for QLEN=4 or 2" )
    sys.exit(1)

  use_simple_rotation = False

  q = map( float, string.split( options.quat_out, ',' ) )
  if ( QLEN == 4 ) :
    quat_outside = Q.makeNormalizedQuat( q[0], q[1], q[2], q[3] )
  else :
    quat_outside = Q.makeNormalizedQuat2( q[0], q[1] )
  
  if ( n_spheres == "one" and options.quat_in is None ) :
    print( "Error: must specify quat-in and quat-out together" )
    sys.exit(1)

  if ( n_spheres == "two" and
       ( options.quat_in_two is None or
         options.quat_in is None ) ) :
    print( "Error: must specify quat-in, quat-in-two, " \
          "and quat-out together" )
    sys.exit(1)

angle_inside = []
angle_outside = None

if ( not options.angle_in is None ) :
  if ( QLEN != 1 ) :
    print( "Error: angle-in option is for QLEN=1" )
    sys.exit(1)

  use_simple_rotation = False

  angle_inside.append( options.angle_in * math.pi / 180.0 )
  
  if ( n_spheres == "one" and options.angle_out is None ) :
    print( "Error: must specify angle-in and angle-out together" )
    sys.exit(1)

  if ( n_spheres == "two" and
       ( options.angle_in_two is None or
         options.angle_out is None ) ) :
    print( "Error: must specify angle-in, angle-in-two, " \
          "and angle-out together" )
    sys.exit(1)

if ( not options.angle_in_two is None ) :
  if ( QLEN != 1 ) :
    print( "Error: angle-in-two option is for QLEN=1" )
    sys.exit(1)

  use_simple_rotation = False

  angle_inside.append( options.angle_in_two * math.pi / 180.0 )
  
  if ( n_spheres == "one" ) :
    print( "Warning: ignoring angle-in-two with only one sphere" )
    sys.exit(1)

  if ( n_spheres == "two" and
       ( options.angle_in is None or
         options.angle_out is None ) ) :
    print( "Error: must specify angle-in, angle-in-two, " \
          "and angle-out together" )
    sys.exit(1)
    
if ( not options.angle_out is None ) :
  if ( QLEN != 1 ) :
    print( "Error: angle-out option is for QLEN=1" )
    sys.exit(1)

  use_simple_rotation = False

  angle_outside = options.angle_out * math.pi / 180.0
  
  if ( n_spheres == "one" and options.angle_in is None ) :
    print( "Error: must specify angle-in and angle-out together" )
    sys.exit(1)

  if ( n_spheres == "two" and
       ( options.angle_in_two is None or
         options.angle_in is None ) ) :
    print( "Error: must specify angle-in, angle-in-two, " \
          "and angle-out together" )
    sys.exit(1)

angle = options.angle
if ( use_simple_rotation ) :
  if ( QLEN == 1 ) :
    angle_outside = 0
  elif ( QLEN == 4 ) :
    quat_outside = Q.makeNormalizedQuat( 1, 0, 0, 0 )
  else :
    quat_outside = Q.makeNormalizedQuat2( 1, 0 )

  if ( angle is None ) :
    angle = 44.0
    print( "Angle not specified: using %f" % angle )

  if ( QLEN == 1 ) :
    angle_inside.append( angle * math.pi / 180.0 )
  else :
    q0 = math.cos( ( math.pi / 180.0 ) * ( angle / 2.0 ) )
    q1 = math.sqrt( 1 - q0 * q0 )
    if ( QLEN == 4 ) :
      quat_inside.append( Q.makeNormalizedQuat( q0, q1, 0, 0 ) )
    else :
      quat_inside.append( Q.makeNormalizedQuat2( q0, q1 ) )

  if ( n_spheres == "two" ) :
    if ( QLEN == 1 ) :
      angle_inside.append( -angle / 2.0 )
    else :
      q0 = math.cos( ( math.pi / 180.0 ) * ( -angle / 2.0 ) )
      q1 = math.sqrt( 1 - q0 * q0 )
      if ( QLEN == 4 ) :
        quat_inside.append( Q.makeNormalizedQuat( q0, q1, 0, 0 ) )
      else :
        quat_inside.append( Q.makeNormalizedQuat2( q0, q1 ) )

if ( n_spheres == "random" ) :
  radii = []
  angle_inside = []
  quat_inside = []
  cx = []
  cy = []
  cz = []
  
  min_radius = 2
  max_radius = min( nx/6, ny/6, nz/6 )
  use_z = True
  if ( nz < 12 ) :
    max_radius = min( nx/6, ny/6 )
    use_z = False

  print( options.random_spheres )
  for n in range( options.random_spheres ) :
    radius = min_radius + int( round( (max_radius - min_radius) * random.random() ) )
    radii.append( radius )

    cx.append( radius + int( round( (nx-2*radius) * random.random() ) ) )
    cy.append( radius + int( round( (ny-2*radius) * random.random() ) ) )
    if ( use_z ) :
      cz.append( radius + int( round( (nz-2*radius) * random.random() ) ) )
    else :
      cz.append( nz / 2 )
    if ( QLEN == 1 ) :
      angle_inside.append( math.pi - 2*math.pi*random.random() )
    elif ( QLEN == 2 ) :
      a = math.pi - 2*math.pi*random.random()
      q0 = math.sin( a )
      q1 = math.sqrt( 1 - q0 * q0 )
      quat_inside.append( Q.makeNormalizedQuat2( q0, q1 ) )
    else :
      raise Error

#-----------------------------------------------------------------------
# Open and define file

f = nc4.Dataset(filename, 'w', format='NETCDF4')

f.createDimension( 'x', nx )
f.createDimension( 'y', ny )
f.createDimension( 'z', nz )
f.createDimension( 'qlen', QLEN )

ncphase = f.createVariable( 'phase', 'd', ('z','y','x') )

ncquat = []
for qq in range( QLEN ) :
  q_comp = f.createVariable( 'quat%d' % (qq+1), 'd', ('z','y','x') )
  ncquat.append( q_comp )

ncconc = f.createVariable( 'concentration', 'd', ('z','y','x') )

phase = N.ones( (nz,ny,nx), N.float64 )
quat = N.zeros( (QLEN,nz,ny,nx), N.float64 )
conc = N.ones( (nz,ny,nx), N.float64 )

if ( n_spheres == "two" ) :
  r_sq_two = radius_two**2

  cx = nx / 2
  cy = ny / 2
  cz = nz / 2

  cx_one = 3 * cx / 4
  cy_one = 3 * cy / 4
  cz_one = cz

  cx_two = 5 * cx / 4
  cy_two = 5 * cy / 4
  cz_two = cz

#-----------------------------------------------------------------------
# Variables and functions for symmetry rotations

N_ROTATIONS = Q.getNumberSymmetryRotations( QLEN )
qrot = []
for n in range( 8 ) :
  nr = int( N_ROTATIONS * random.random() )
  qrot.append( Q.getQuatSymmetryRotation( nr, qlen=QLEN ) )

def rotateRandomly() :
  nr = int( N_ROTATIONS * random.random() )
  return Q.getQuatSymmetryRotation( nr, qlen=QLEN )

def rotateByOctant( x, y, z, cx, cy, cz ) :
  global qrot
  if ( x < cx ) :
    if ( y < cy ) :
      if ( z < cz ) :
        qr = qrot[0]
      else :
        qr = qrot[1]
    else :
      if ( z < cz ) :
        qr = qrot[2]
      else :
        qr = qrot[3]
  else :
    if ( y < cy ) :
      if ( z < cz ) :
        qr = qrot[4]
      else :
        qr = qrot[5]
    else :
      if ( z < cz ) :
        qr = qrot[6]
      else :
        qr = qrot[7]
  return qr

#-----------------------------------------------------------------------
# Fill data arrays

phase[:,:,:] = phase_outside
conc[:,:,:] = conc_outside
if ( QLEN == 1 ) :
  quat[:,:,:,:] = angle_outside
else :
  for qq in range( QLEN ) :
    quat[qq,:,:,:] = quat_outside[qq]

for n in range( len( radii ) ) :

  if ( QLEN == 1 ) :
    print( cx[n], cy[n], cz[n], radii[n], angle_inside[n] )
  else :
    print( cx[n], cy[n], cz[n], radii[n], quat_inside[n] )

  for k in range( nz ) :
    for j in range( ny ) :
      for i in range( nx ) :
        x = i + 0.5
        y = j + 0.5
        z = k + 0.5

        distance_sq = (x - cx[n])**2 + (y - cy[n])**2 + (z - cz[n])**2
        r_sq = radii[n] * radii[n]

        if ( distance_sq < r_sq ) :
          phase[k,j,i] = phase_inside

          if ( QLEN == 1 ) :
            quat[0,k,j,i] = angle_inside[n]
          else :
            for qq in range( QLEN ) :
              quat[qq,k,j,i] = quat_inside[n][qq]

          conc[k,j,i] = conc_inside

for k in range( nz ) :
  for j in range( ny ) :
    for i in range( nx ) :
      x = i + 0.5
      y = j + 0.5
      z = k + 0.5

      if ( QLEN == 1 ) :
        q = quat[0,k,j,i]
      else :
        q = quat[:,k,j,i] 

      qr = None
      if ( options.symmetry_test ) :
        qr = rotateByOctant( x, y, z, nx/2, ny/2, nz/2 )
      elif ( options.random_symmetry ) :
        qr = rotateRandomly()

      if ( not qr is None ) :
        if ( QLEN == 1 ) :
          q += qr
        else :
          q = Q.rotateQuat( q, qr )

      if ( QLEN == 1 ) :
        quat[0,k,j,i] = q
      else :
        for qq in range( QLEN ) :
          quat[qq,k,j,i] = q[qq]

#-----------------------------------------------------------------------
# Write data to file and close

ncphase[:,:,:]=phase

for n in range( QLEN ) :
  ncquat[n][:,:,:]=quat[n,:,:,:]

ncconc[:,:,:]= conc

f.close()
