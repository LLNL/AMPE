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
import Numeric as N
from Scientific.IO import NetCDF

#-----------------------------------------------------------------------
# command-line arguments and options

usage = "Usage: %prog [options] filename"

parser = OptionParser( usage = usage )

parser.add_option( "--ngrains", type="int",
                   help="number of grains to nucleate" )
parser.add_option( "-n", "--ncells", type="int",
                   help="number of cells in x,y,z (will all be equal)" )
parser.add_option( "-x", "--nx", type="int",
                   help="number of cells in x direction" )
parser.add_option( "-y", "--ny", type="int",
                   help="number of cells in y direction" )
parser.add_option( "-z", "--nz", type="int",
                   help="number of cells in z direction" )
parser.add_option( "-q", "--qlen", type="int",
                   help="number of component for q" )
parser.add_option( "-r", "--radius", type="int",
                   help="number of cells for sphere radius [default: nx/4]" )
parser.add_option( "-i", "--phase-in", type="float", default=1.0,
                   help="phase in interior region [default: %default]" )
parser.add_option( "-o", "--phase-out", type="float", default=0.0,
                   help="phase in exterior region [default: %default]" )
parser.add_option( "-c", "--concentration-in", type="float", default=0.1,
                   help="concentration in interior region [default: %default]" )
parser.add_option( "--concentration-out", type="float", default=0.06,
                   help="concentration in exterior region [default: %default]" )

(options, args) = parser.parse_args()

filename = args[0]

n_spheres = options.ngrains

nx = options.ncells
ny = options.ncells
nz = options.ncells

if ( options.nx ) : nx = options.nx
if ( options.ny ) : ny = options.ny
if ( options.nz ) : nz = options.nz

if ( not ( nx and ny and nz ) ) :
  print "Error: either -n or all of -nx -ny -nz are required"
  sys.exit(1)

qlen =4
if ( options.qlen ) : qlen = options.qlen

radius = options.radius
if radius is None :
  radius = nx / 10
  print "Radius not specified: using nx/n_spheres = %d" % radius

phase_inside  = options.phase_in
phase_outside = options.phase_out
conc_inside   = options.concentration_in
conc_outside  = options.concentration_out

quat_outside = Q.makeNormalizedQuat( 1, 0, 0, 0 )

# generate quaternions corresponding to 3 fcc orientations
quat_inside=[]
for g in range(n_spheres):

  q0=math.cos( math.pi/8. )
  q1=math.sin( math.pi/8. )
  
  q=[]
  q.append(q0)
  for i in range(3):
    if( (g+i)%3 == 0 ):
      q.append(q1)
    else:
      q.append(0.)

  print q
  quat_inside.append( Q.makeNormalizedQuat( q[0], q[1], q[2], q[3] ) )

#distance square function for periodic bc
def distance2(x1,y1,z1,x2,y2,z2):
  d2=(x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2
  for i in range(-1,2):
    dx2=(x1-x2+i*nx)**2
    for j in range(-1,2):
      dy2=(y1-y2+j*ny)**2
      for k in range(-1,2):
        dz2=(z1-z2+k*nz)**2
        l2=dx2+dy2+dz2
        if l2<d2:
          d2=l2
  return d2

#-----------------------------------------------------------------------
# simple algorithm to move nuclei away from each other using a dynamics
# based on forces between pairs.
# Forces depend on distances and go to zero at distance drmin
def imposeMinDistance2d(cx,cy,nx,ny,drmin):

  ncenters=len(cx)       

  fac=[0.,-1.,1.]

  step=1
  npairs=1000
  while npairs>0:
    print 'Step ',step
    npairs=0

    #reset forces
    forces=[]
    for index in range(ncenters):
      forces.append([0.,0.])
    
    for index0 in range(ncenters):
      x0=cx[index0]
      y0=cy[index0]
      for index1 in range(ncenters):
        if(index1!=index0):
          x1=cx[index1]
          y1=cy[index1]
      
          dr0=1000.
          dx0=nx
          dy0=ny
          #check minimage
          for i in range (3):
            dx=x1-x0+fac[i]*nx
            if dx<dr0:
              for j in range (3):
                dy=y1-y0+fac[j]*ny
                dr=math.sqrt(dx*dx+dy*dy)
                if dr<dr0:
                  dr0=dr
                  dx0=dx
                  dy0=dy
          #print dr0
          if dr0<drmin:
            factor=0.5
            if dr0>0.000001:
              af=factor*(drmin-dr0+0.0001)
              invdr0=1./dr0
              fx=af*dx0*invdr0
              fy=af*dy0*invdr0
            else:
              fx=factor*drmin*random.random()
              fy=factor*drmin*random.random()
            forces[index0][0]=forces[index0][0]-fx
            forces[index0][1]=forces[index0][1]-fy
            forces[index1][0]=forces[index1][0]+fx
            forces[index1][1]=forces[index1][1]+fy
              
            npairs=npairs+1
    print 'number of pairs=',npairs
    step=step+1
   
    # move cx and cy
    maxf=0.
    for index in range(ncenters):
      cx[index]=cx[index]+forces[index][0];
      cy[index]=cy[index]+forces[index][1];
      norm2f=forces[index][0]*forces[index][0]
      norm2f=norm2f+forces[index][1]*forces[index][1]
      
      if( norm2f>maxf ):
        maxf=norm2f
    print 'max force=',math.sqrt(maxf)

#-----------------------------------------------------------------------
# Open and define file

f = NetCDF.NetCDFFile( filename, 'w' )

f.createDimension( 'x', nx )
f.createDimension( 'y', ny )
f.createDimension( 'z', nz )
f.createDimension( 'qlen', qlen )

ncphase = f.createVariable( 'phase', N.Float, ('z','y','x') )
ncquat=[]
for m in range(qlen):
  name = "quat%d" % (m+1)
  print name
  ncquat.append( f.createVariable( name, N.Float, ('z','y','x') ) )
ncconc = f.createVariable( 'concentration', N.Float, ('z','y','x') )

phase = N.ones( (nz,ny,nx), N.Float )
quat = N.zeros( (qlen,nz,ny,nx), N.Float )
conc = N.ones( (nz,ny,nx), N.Float )

r_sq = radius**2

# generate random (x,y) positions
random.seed(21361)
cx=[]
cy=[]
cz=[]
for g in range(n_spheres):
  cx.append( random.randint(1,nx) )
for g in range(n_spheres):
  cy.append( random.randint(1,ny) )
for g in range(n_spheres):
  cz.append( nz/2 )

# move positions to avoid any overlaps between grains
imposeMinDistance2d(cx,cy,nx,ny,2.4*radius)

#-----------------------------------------------------------------------
# Fill data arrays
for k in range( nz ) :
  for j in range( ny ) :
    for i in range( nx ) :
      phase[k,j,i] = phase_outside
      conc[k,j,i] = conc_outside
      q = quat_outside
      for m in range(qlen):
        quat[m,k,j,i] = q[m]

for k in range( nz ) :
  for j in range( ny ) :
    for i in range( nx ) :
      x = i + 0.5
      y = j + 0.5
      z = k + 0.5
      for g in range(n_spheres):
      
        distance_sq = distance2(x,y,z,cx[g],cy[g],cz[g])
        if ( distance_sq < r_sq ) :
          phase[k,j,i] = phase_inside
          conc[k,j,i] = conc_inside
          q = quat_inside[g]
          for m in range(qlen):
            quat[m,k,j,i] = q[m]

#-----------------------------------------------------------------------
# Write data to file and close

ncphase.assignValue( phase )
for m in range(qlen):
  ncquat[m].assignValue( quat[m,...] )
ncconc.assignValue( conc )

f.close()


