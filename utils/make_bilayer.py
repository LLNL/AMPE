# Copyright (c) 2018, Lawrence Livermore National Security, LLC and
# UT-Battelle, LLC.
# Produced at the Lawrence Livermore National Laboratory and
# the Oak Ridge National Laboratory
# Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
# LLNL-CODE-747500
# All rights reserved.
# This file is part of AMPE. 
# For details, see https://github.com/LLNL/AMPE
# Please also read AMPE/LICENSE.
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# - Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the disclaimer below.
# - Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the disclaimer (as noted below) in the
#   documentation and/or other materials provided with the distribution.
# - Neither the name of the LLNS/LLNL nor the names of its contributors may be
#   used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
# LLC, UT BATTELLE, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
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
#from Scientific.IO import NetCDF
from scipy.io import netcdf as NetCDF

import netCDF4 as nc4

#-----------------------------------------------------------------------
# command-line arguments and options

usage = "Usage: %prog [options] filename"

parser = OptionParser( usage = usage )

parser.add_option( "-1", "--one", action="store_true",
                   help="one grain centered in the middle of the domain [default]" )
parser.add_option( "-2", "--two", action="store_true",
                   default=False,
                   help="two grains centered in opposite halves" )
parser.add_option( "--ratio", type="float", default=1.0,
                   help="ratio of grain 1 to grain 2 [default: %default]" )
parser.add_option( "-n", "--ncells", type="int",
                   help="number of cells in x,y,z (will all be equal)" )
parser.add_option( "-x", "--nx", type="int",
                   help="number of cells in x direction" )
parser.add_option( "-y", "--ny", type="int",
                   help="number of cells in y direction" )
parser.add_option( "-z", "--nz", type="int",
                   help="number of cells in z direction" )
parser.add_option( "-r", "--radius", type="int",
                   help="number of cells for grain half-width [default: nx/4]" )
parser.add_option( "-a", "--angle", type="float",
                   help="angular difference between interior and exterior" )
parser.add_option( "-d","--delta", type="int", default=0,
                   help="interface width (in number of cells) [default: 0.]" )
parser.add_option( "-i", "--phase-in", type="float", default=1.0,
                   help="phase in interior region [default: %default]" )
parser.add_option( "-o", "--phase-out", type="float", default=0.0,
                   help="phase in exterior region [default: %default]" )
parser.add_option( "--orient", type="string", # "x", "y", or "z"
                   help="direction of orientation of the bilayer" )
parser.add_option( "--qlen", type="int", default=4,
                   help="Length of orientation vector" )
parser.add_option( "--quat-in", type="string", # something like "1,0,0,0",
                   help="quaternion in interior region" )
parser.add_option( "--quat-in-two", type="string", # something like "1,0,0,0",
                   help="quaternion in interior region of grain 2" )
parser.add_option( "--quat-out", type="string",
                   help="quaternion in exterior region" )
parser.add_option( "--angle-in", type="float",
                   help="grain angle in interior region" )
parser.add_option( "--angle-in-two", type="float",
                   help="grain angle in interior region of grain 2" )
parser.add_option( "--angle-out", type="float",
                   help="grain angle in exterior region" )
parser.add_option( "-c", "--nomconc", type="string", # something like "0.1,0.2",
                   help="nominal concentration" )
parser.add_option( "--concentration-in", type="string",
                   help="concentration in interior region" )
parser.add_option( "--concentration-out", type="string", 
                   help="concentration in exterior region" )
parser.add_option( "--eta-in", type="float",
                   help="eta in interior region [default: %default]" )
parser.add_option( "--eta-out", type="float", default=0.,
                   help="eta in exterior region [default: %default]" )
parser.add_option( "-s", "--symmetry-test", action="store_true",
                   default=False,
                   help="rotate each octant by a random symmetry rotation" )
parser.add_option( "--random-symmetry", action="store_true",
                   default=False,
                   help="rotate each *CELL* by a random symmetry rotation" )
parser.add_option( "--temperature0", type="float",
                   help="temperature0" )
parser.add_option( "--gradxT", type="float",
                   help="gradxT" )
parser.add_option( "--gaussT", type="float",
                   help="gaussT" )
parser.add_option( "--double", action="store_true", dest="double_precision", default=False)
parser.add_option( "--centerx", type="float", 
                   help="x position of grain 1" )
parser.add_option( "--centery", type="float", 
                   help="y position of grain 1" )
parser.add_option( "--centerz", type="float", 
                   help="z position of grain 1" )

(options, args) = parser.parse_args()

if ( len( args ) < 1 ) :
  parser.error( "filename required argument missing" )
filename = args[0]

QLEN = options.qlen
if ( QLEN != 4 and QLEN != 2 and QLEN != 1 ) :
  print "Error: valid values of qlen are 1, 2, and 4"
  sys.exit(1)
else:
  print "QLEN=",QLEN
  
n_solid_layers = "one"
if ( options.two ) :
  if ( options.one ) :
    print "Error: cannot use options one and two together"
    sys.exit(1)
  n_solid_layers = "two"

double_precision = options.double_precision
if double_precision:
  print "use double precision..."

nx = options.ncells
ny = options.ncells
nz = options.ncells

if ( not options.nx is None ) : nx = options.nx
if ( not options.ny is None ) : ny = options.ny
if ( not options.nz is None ) : nz = options.nz

if ( not ( nx and ny and nz ) ) :
  print "Error: either -n or all of -x -y -z are required"
  sys.exit(1)

if ( options.two ) :
  if ( options.one ) :
    print "Error: cannot use options one and two together"
    sys.exit(1)

if ( not options.centerx is None ) : 
  centerx = options.centerx
else:
  centerx = nx / 2

if ( not options.centery is None ) :
  centery = options.centery
else:
  centery = ny / 2

if ( not options.centerz is None ) :
  centerz = options.centerz
else:
  centerz = nz / 2

nomconc = options.nomconc
temperature0 = options.temperature0
gradxT = options.gradxT
gaussT = options.gaussT

direction = 0
width = nx
if ( not options.orient is None ) :
  if ( options.orient == "x" ) :
    direction = 0
    width = nx
  elif ( options.orient == "y" ) :
    direction = 1
    width = ny
  elif ( options.orient == "z" ) :
    direction = 2
    width = nz
  else :
    print "Error: orient should be x, y, or z"
    sys.exit(1)

radius = options.radius
if radius is None :
  if ( n_solid_layers == "one" ) :
    radius = width / 4
    print "Radius not specified: using width/4 = %d" % radius
  else :
    radius = width / 8
    print "Radius not specified: using width/8 = %d" % radius

delta = options.delta
if delta is None :
  delta=0.
print 'interface width (in cells) =',delta

phase_inside = options.phase_in
phase_outside = options.phase_out

eta_inside = options.eta_in
eta_outside = options.eta_out

conc_inside   = options.concentration_in
conc_outside  = options.concentration_out

n_solid_layers = "one"
radius_two = 0.
if ( options.two ) :
  n_solid_layers = "two"
  ratio = options.ratio
  radius_two = radius * ratio

use_simple_rotation = True #if quaternions not specified

quat_inside     = None
quat_inside_two = None
quat_outside    = None

#set quat_inside
if ( not options.quat_in is None ) :
  #check for incompatible options first...
  if ( QLEN != 4 and QLEN != 2 ) :
    print "Error: quat-in option is for QLEN=4 or 2"
    sys.exit(1)
  if ( n_solid_layers == "one" and options.quat_out is None ) :
    print "Error: must specify quat-in and quat-out together"
    sys.exit(1)
  if ( n_solid_layers == "two" and
       ( options.quat_in_two is None or
         options.quat_out is None ) ) :
    print "Error: must specify quat-in, quat-in-two, " \
          "and quat-out together"
    sys.exit(1)

  use_simple_rotation = False

  q = map( float, string.split( options.quat_in, ',' ) )
  if ( QLEN == 4 ) :
    quat_inside = Q.makeNormalizedQuat( q[0], q[1], q[2], q[3] )
  else :
    quat_inside = Q.makeNormalizedQuat2( q[0], q[1] )


#set quat_inside_two
if ( not options.quat_in_two is None ) :
  #check for incompatible options first...
  if ( QLEN != 4 and QLEN != 2 ) :
    print "Error: quat-in-two option is for QLEN=4 or 2"
    sys.exit(1)
  if ( n_solid_layers == "two" and
       ( options.quat_in is None or
         options.quat_out is None ) ) :
    print "Error: must specify quat-in, quat-in-two, " \
          "and quat-out together"
    sys.exit(1)
  if ( n_solid_layers == "one" ) :
    print "Warning: ignoring quat-in-two with only one layer"
  else:
    q = map( float, string.split( options.quat_in_two, ',' ) )
    if ( QLEN == 4 ) :
      quat_inside_two = Q.makeNormalizedQuat( q[0], q[1], q[2], q[3] )
    else :
      quat_inside_two = Q.makeNormalizedQuat2( q[0], q[1] )


#set quat_outside
if ( not options.quat_out is None ) :
  #check for incompatible options first...
  if ( QLEN != 4 and QLEN != 2 ) :
    print "Error: quat-out option is for QLEN=4 or 2"
    sys.exit(1)
  if ( n_solid_layers == "one" and options.quat_in is None ) :
    print "Error: must specify quat-in and quat-out together"
    sys.exit(1)
  if ( n_solid_layers == "two" and
       ( options.quat_in_two is None or
         options.quat_in is None ) ) :
    print "Error: must specify quat-in, quat-in-two, " \
          "and quat-out together"
    sys.exit(1)

  use_simple_rotation = False

  q = map( float, string.split( options.quat_out, ',' ) )
  if ( QLEN == 4 ) :
    quat_outside = Q.makeNormalizedQuat( q[0], q[1], q[2], q[3] )
  else :
    quat_outside = Q.makeNormalizedQuat2( q[0], q[1] )
  


if ( not options.angle_in is None ) :
  if ( QLEN != 1 ) :
    print "Error: angle-in option is for QLEN=1"
    sys.exit(1)

  use_simple_rotation = False

  angle_inside = options.angle_in * math.pi / 180.0
  
  if ( n_solid_layers == "one" and options.angle_out is None ) :
    print "Error: must specify angle-in and angle-out together"
    sys.exit(1)

  if ( n_solid_layers == "two" and
       ( options.angle_in_two is None or
         options.angle_out is None ) ) :
    print "Error: must specify angle-in, angle-in-two, " \
          "and angle-out together"
    sys.exit(1)

if ( not options.angle_in_two is None ) :
  if ( QLEN != 1 ) :
    print "Error: angle-in-two option is for QLEN=1"
    sys.exit(1)

  use_simple_rotation = False

  angle_inside_two = options.angle_in_two * math.pi / 180.0
  
  if ( n_solid_layers == "one" ) :
    print "Warning: ignoring angle-in-two with only one layer"
    sys.exit(1)

  if ( n_solid_layers == "two" and
       ( options.angle_in is None or
         options.angle_out is None ) ) :
    print "Error: must specify angle-in, angle-in-two, " \
          "and angle-out together"
    sys.exit(1)
    
if ( not options.angle_out is None ) :
  if ( QLEN != 1 ) :
    print "Error: angle-out option is for QLEN=1"
    sys.exit(1)

  use_simple_rotation = False

  angle_outside = options.angle_out * math.pi / 180.0
  
  if ( n_solid_layers == "one" and options.angle_in is None ) :
    print "Error: must specify angle-in and angle-out together"
    sys.exit(1)

  if ( n_solid_layers == "two" and
       ( options.angle_in_two is None or
         options.angle_in is None ) ) :
    print "Error: must specify angle-in, angle-in-two, " \
          "and angle-out together"
    sys.exit(1)

#-----------------------------------------------------------------------
# specify quaternions using angles
angle = options.angle
if ( use_simple_rotation ) :
  if ( QLEN == 1 ) :
    angle_inside = 0
  elif ( QLEN == 4 ) :
    quat_inside = Q.makeNormalizedQuat( 1, 0, 0, 0 )
  else :
    quat_inside = Q.makeNormalizedQuat2( 1, 0 )

  if ( angle is None ) :
    angle = 44.0
    print "Angle not specified: using %f" % angle

  angle = angle * math.pi / 180.0

  if ( QLEN == 1 ) :
    angle_outside = angle
  else :
    if ( QLEN == 4 ) :
      q0 = math.cos( angle / 2.0 )
    else :
      q0 = math.cos( angle )
    q1 = math.sqrt( 1 - q0 * q0 )
    if ( QLEN == 4 ) :
      quat_outside = Q.makeNormalizedQuat( q0, q1, 0, 0 )
    else :
      quat_outside = Q.makeNormalizedQuat2( q0, q1 )
    print "Quat outside grain 1 =", quat_outside

  if ( n_solid_layers == "two" ) :
    if ( QLEN == 1 ) :
      angle_inside_two = angle_outside
      angle_outside = angle / 2.0
    else :
      quat_inside_two = quat_outside
      q0 = math.cos( angle / 4.0 )
      q1 = math.sqrt( 1 - q0 * q0 )
      if ( QLEN == 4 ) :
        quat_outside = Q.makeNormalizedQuat( q0, q1, 0, 0 )
      else :
        quat_outside = Q.makeNormalizedQuat2( q0, q1 )

#-----------------------------------------------------------------------
nspecies=0
if ( not ( nomconc is None ) ):
  c = map( float, string.split( options.nomconc, ',' ) )
  nspecies=len(c)
  print "Nominal composition=",c
if ( not ( conc_inside is None ) ):
  ci = map( float, string.split( options.concentration_in, ',' ) )
  if nspecies==0:
    nspecies=len(ci)
  print "Composition inside=",ci
else:
  ci = N.zeros( nspecies, N.float32 )
if ( not ( conc_outside is None ) ):
  co = map( float, string.split( options.concentration_out, ',' ) )
  print "Composition outside=",co
else:
  co = N.zeros( nspecies, N.float32 )

print "nspecies=",nspecies

#-----------------------------------------------------------------------
# Open and define file

#f = NetCDF.NetCDFFile( filename, 'w' )
f = nc4.Dataset(filename, 'w', format='NETCDF4')

f.createDimension( 'x', nx )
f.createDimension( 'y', ny )
f.createDimension( 'z', nz )
f.createDimension( 'qlen', QLEN )
f.createDimension( 'ns', nspecies )
ncquat = []
ncconc = []

if double_precision:
  ncphase       = f.createVariable( 'phase', 'd', ('z','y','x') )
  if( not eta_inside is None ):
    nceta         = f.createVariable( 'eta', 'd', ('z','y','x') )
  nctemperature = f.createVariable( 'temperature', 'd', ('z','y','x') )
  for n in range( QLEN ) :
    q_comp = f.createVariable( 'quat%d' % (n+1), 'd', ('z','y','x') )
    ncquat.append( q_comp )
  for s in range(nspecies):
    c_comp = f.createVariable( 'concentration%d' % s, 'd', ('z','y','x') )
    ncconc.append(c_comp)
else:
  ncphase       = f.createVariable( 'phase', 'f', ('z','y','x') )
  if( not eta_inside is None ):
    nceta         = f.createVariable( 'eta', 'f', ('z','y','x') )
  nctemperature = f.createVariable( 'temperature', 'f', ('z','y','x') )
  for n in range( QLEN ) :
    q_comp = f.createVariable( 'quat%d' % (n+1), 'f', ('z','y','x') )
    ncquat.append( q_comp )
  for s in range(nspecies):
    c_comp = f.createVariable( 'concentration%d' % s , 'f', ('z','y','x') )
    ncconc.append(c_comp)



if double_precision:
  phase = N.ones( (nz,ny,nx), N.float64 )
  if( not eta_inside is None ):
    eta   = N.ones( (nz,ny,nx), N.float64 )
  quat  = N.zeros( (QLEN,nz,ny,nx), N.float64 )
  conc  = N.ones( (nspecies,nz,ny,nx), N.float64 )
  temperature = N.zeros( (nz,ny,nx), N.float64 )
else:
  phase = N.ones( (nz,ny,nx), N.float32 )
  if( not eta_inside is None ):
    eta   = N.ones( (nz,ny,nx), N.float32 )
  quat  = N.zeros( (QLEN,nz,ny,nx), N.float32 )
  conc  = N.ones( (nspecies,nz,ny,nx), N.float32 )
  temperature = N.zeros( (nz,ny,nx), N.float32 )

r_sq = radius**2
r_sq_two = radius_two**2

cx = centerx
cy = centery
cz = centerz

cx_one = cx
cy_one = cy
cz_one = cz

if n_solid_layers == "two" :
  cx_one = cx / 2
  cy_one = cy / 2
  cx_two = 3 * cx / 2
  cy_two = 3 * cy / 2
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
vol=nx*ny*nz
if ( direction == 0 ) :
  vs=2*radius*ny*nz
elif ( direction == 1 ) :
  vs=2*radius*nx*nz
elif ( direction == 2 ) :
  vs=2*radius*ny*nx

if ( n_solid_layers == "two" ) :
  vs=2.*vs

vs=min(vs,vol)
vl=vol-vs
print "Volume           = ",vol
print "Volume of solid  = ",vs
print "Volume of liquid = ",vl

  
if ( not ( conc_inside is None ) ):
  if ( not ( nomconc is None ) and vl>0 ):
    for s in range(nspecies):
      co[s] = (c[s]*vol-ci[s]*vs)/vl
    print "Calculated composition outside=",co
if ( not ( options.concentration_out is None ) ):
  if ( not ( nomconc is None ) and vs>0 ):
    conc_inside = (c[0]*vol-conc_outside*vl)/vs
    print "Calculated composition inside=",conc_inside
if( ( conc_outside is None ) and ( conc_inside is None ) ):
  conc_inside = nomconc
  conc_outside = nomconc

if ( not ( conc_outside is None ) and not ( conc_inside is None ) ):
  for s in range(nspecies):
    print "Calculated nominal Composition=",(vl*co[s]+vs*ci[s])/vol

if( delta>0. ):
  invdelta = 1./delta

if ( QLEN > 1 ):
  print "Quat inside grain 1 =", quat_inside
  if ( n_solid_layers == "two" ) :
    print "Quat inside grain 2 =", quat_inside_two

for k in range( nz ) :
  for j in range( ny ) :
    for i in range( nx ) :
      x = i + 0.5
      y = j + 0.5
      z = k + 0.5

      if ( direction == 0 ) :
        distance_sq1 = (x - cx_one)**2
      elif ( direction == 1 ) :
        distance_sq1 = (y - cy_one)**2
      elif ( direction == 2 ) :
        distance_sq1 = (z - cz_one)**2

      d = distance_sq1 - r_sq
      if( delta>0. ):
        sd=math.sqrt(distance_sq1)-math.sqrt(r_sq)
        v = 0.5*(1.+math.tanh(-0.5*sd*invdelta))
      else :
        if( d>0. ):
          v=0.; #outside
        else:
          v=1.; #inside

      if ( not temperature0 is None ) :
        temperature[k,j,i] = temperature0
        if ( not gradxT is None ) :
          temperature[k,j,i] = temperature[k,j,i] + gradxT*(i-0.5*nx+0.5)
        if ( not gaussT is None ) :
          d2=(x - cx_one)**2
          temperature[k,j,i] = temperature[k,j,i] + gaussT*math.exp(-distance_sq1/(0.0625*nx*nx))
      
      #smooth interface
      phase[k,j,i] = v*phase_inside+(1.-v)*phase_outside
      
      if ( d<0. ) :  #inside grain 1
        if( not eta_inside is None ):
          eta[k,j,i]   = eta_inside
        if ( nspecies>0 ):
          c = ci
      else : #outside
        if( not eta_inside is None ):
          eta[k,j,i]   = eta_outside
        if ( nspecies>0 ):
          c = co
      if ( QLEN == 1 ) :
        q = angle_inside
      else :
        q = quat_inside

      if ( n_solid_layers == "two" ) :
        if ( direction == 0 ) :
          distance_sq2 = (x - cx_two)**2
        elif ( direction == 1 ) :
          distance_sq2 = (y - cy_two)**2
        elif ( direction == 2 ) :
          distance_sq2 = (z - cz_two)**2
        d = distance_sq2 - r_sq_two
        if( delta>0. ):
          v = 0.5*(1.+math.tanh(-0.5*d*invdelta))
        else :
          if( d>0. ):
            v=0.;
          else:
            v=1.;

        if ( d<delta ) : #inside grain 2
          phase[k,j,i] = v*phase_inside+(1.-v)*phase_outside
          if( not eta_inside is None ):
            eta[k,j,i]   = eta_inside
          c = ci
        
        if ( distance_sq2<distance_sq1 ) : 
          if ( QLEN == 1 ) :
            q = angle_inside_two
          else :
            q = quat_inside_two
        
      qr = None
      if ( options.symmetry_test ) :
        qr = rotateByOctant( x, y, z, cx, cy, cz )
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
        for n in range( QLEN ) :
          quat[n,k,j,i] = q[n]
      if ( nspecies>0 ):
        for s in range(nspecies):
          conc[s,k,j,i] = c[s]

#-----------------------------------------------------------------------
# Write data to file and close

ncphase[:,:,:]=phase
if( not eta_inside is None ):
  nceta[:,:,:]=eta
nctemperature[:,:,:]= temperature

for n in range( QLEN ) :
  ncquat[n][:,:,:]=quat[n,:,:,:]

if ( nspecies>0 ):
  for s in range(nspecies):
    ncconc[s][:,:,:]=conc[s,:,:,:]

f.close()
