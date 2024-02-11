# Copyright (c) 2018, Lawrence Livermore National Security, LLC and
# UT-Battelle, LLC.
# Produced at the Lawrence Livermore National Laboratory and
# Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
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
#import numpy as N
from Scientific.IO import NetCDF

#-----------------------------------------------------------------------
# command-line arguments and options

usage = "Usage: %prog [options] filename"

parser = OptionParser( usage = usage )

parser.add_option( "--ngrains", type="int", default=1,
                   help="number of grains to nucleate" )
parser.add_option( "-n", "--ncells", type="int",
                   help="number of cells in x,y,z (will all be equal)" )
parser.add_option( "-d", "--dimension", type="int",
                   help="dimension of subspace containing centers [default:3]" )
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
parser.add_option( "-c", "--concentration", type="float", default=0.06,
                   help="nominal concentration [default: %default]" )
parser.add_option( "--concentration-in",  type="float", default=0.10,
                   help="concentration in interior region [default: %default]" )

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

vol = nx * ny * nz

ndim=3
if ( options.dimension ) : ndim = options.dimension

qlen =4
if ( options.qlen ) : qlen = options.qlen

radius = options.radius
if radius is None :
  radius = nx / 10
  print "Radius not specified: using nx/n_spheres = %d" % radius

phase_inside  = options.phase_in
phase_outside = options.phase_out
conc_nominal  = options.concentration
conc_inside   = options.concentration_in

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
  for ir in range(-1,2):
    dx2=(x1-x2+ir*nx)**2
    if dx2<=d2:
      for jr in range(-1,2):
        dy2=(y1-y2+jr*ny)**2
        if dy2<=d2:
          for kr in range(-1,2):
            dz2=(z1-z2+kr*nz)**2
            l2=dx2+dy2+dz2
            if l2<d2:
              d2=l2
  return d2

#distance square function for periodic bc
def distance2_1d(x1,x2):
  d2=(x1 - x2)**2
  for ir in range(-1,2):
    l2=(x1-x2+ir*nx)**2
    if l2<d2:
      d2=l2
  return d2

#-----------------------------------------------------------------------
# simple algorithm to move nuclei away from each other using a dynamics
# based on forces between pairs.
# Forces depend on distances and go to zero at distance drmin
def imposeMinDistance(cx,cy,cz,nx,ny,nz,drmin):

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
      forces.append([0.,0.,0.])
    
    for index0 in range(ncenters):
      x0=cx[index0]
      y0=cy[index0]
      z0=cz[index0]
      for index1 in range(index0+1,ncenters):
        if(index1!=index0):
          x1=cx[index1]
          y1=cy[index1]
          z1=cz[index1]
      
          dr0=1000.
          dx0=nx
          dy0=ny
          dz0=nz
          #check minimage
          for i in range (3):
            dx=x1-x0+fac[i]*nx
            if dx<dr0:
              for j in range (3):
                dy=y1-y0+fac[j]*ny
                if dy<dr0:
                  for k in range (3):
                    dz=z1-z0+fac[k]*nz
                    dr=math.sqrt(dx*dx+dy*dy+dz*dz)
                    if dr<dr0:
                      dr0=dr
                      dx0=dx
                      dy0=dy
                      dz0=dz
          #print dr0
          if dr0<drmin:
            factor=0.5
            if dr0>0.000001:
              af=factor*(drmin-dr0+0.0001)
              invdr0=1./dr0
              fx=af*dx0*invdr0
              fy=af*dy0*invdr0
              fz=af*dz0*invdr0
            else:
              fx=factor*drmin*random.random()
              fy=factor*drmin*random.random()
              fz=factor*drmin*random.random()
            forces[index0][0]=forces[index0][0]-fx
            forces[index0][1]=forces[index0][1]-fy
            forces[index0][2]=forces[index0][2]-fz
            forces[index1][0]=forces[index1][0]+fx
            forces[index1][1]=forces[index1][1]+fy
            forces[index1][2]=forces[index1][2]+fz
              
            npairs=npairs+1
    print 'number of overlaping pairs=',npairs
    step=step+1
   
    # move cx and cy
    maxf=0.
    for index in range(ncenters):
      cx[index]=cx[index]+forces[index][0];
      cy[index]=cy[index]+forces[index][1];
      cz[index]=cz[index]+forces[index][2];
      norm2f=forces[index][0]*forces[index][0]
      norm2f=norm2f+forces[index][1]*forces[index][1]
      norm2f=norm2f+forces[index][2]*forces[index][2]
      
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

ncconc  = f.createVariable( 'concentration', N.Float, ('z','y','x') )
ncphase = f.createVariable( 'phase',         N.Float, ('z','y','x') )

ncquat=[]
for m in range(qlen):
  name = "quat%d" % (m+1)
  print name
  ncquat.append( f.createVariable( name, N.Float, ('z','y','x') ) )

conc  = N.ones( (nz,ny,nx), N.Float )
phase = N.ones( (nz,ny,nx), N.Float )
quat  = N.zeros( (qlen,nz,ny,nx), N.Float )

r_sq = radius**2

# generate random (x,y) positions
random.seed(21361)
cx=[]
cy=[]
cz=[]
if n_spheres>1:
  for g in range(n_spheres):
    cx.append( random.randint(1,nx) )
  for g in range(n_spheres):
    cy.append( random.randint(1,ny) )
  if ndim>2:
    for g in range(n_spheres):
      cz.append( random.randint(1,nz) )
  else:
    for g in range(n_spheres):
      cz.append( nz/2 )
else:
  cx.append( nx/2 )
  cy.append( ny/2 )
  cz.append( nz/2 )

# move positions to avoid any overlaps between grains
imposeMinDistance(cx,cy,cz,nx,ny,nz,3.*radius)

#-----------------------------------------------------------------------
# Fill data arrays
for k in range( nz ) :
  print 'layer ',k
  for j in range( ny ) :
    for i in range( nx ) :
      phase[k,j,i] = phase_outside
      q = quat_outside
      for m in range(qlen):
        quat[m,k,j,i] = q[m]

for g in range(n_spheres):
  print 'sphere ',g,', center: ',cx[g],cy[g],cz[g]
  q = quat_inside[g]
  for k in range( nz ) :
    z = k + 0.5
    dz2=distance2_1d(z,cz[g])
    if dz2<r_sq :
      for j in range( ny ) :
        y = j + 0.5
        dy2=distance2_1d(y,cy[g])
        if dy2<r_sq :
          for i in range( nx ) :
            x = i + 0.5
       
            distance_sq = distance2(x,y,z,cx[g],cy[g],cz[g])
            if ( distance_sq < r_sq ) :
              phase[k,j,i] = phase_inside
              for m in range(qlen):
                quat[m,k,j,i] = q[m]

#integral phase
int_phase=0.
for k in range( nz ) :
  for j in range( ny ) :
    for i in range( nx ) :
      int_phase=int_phase+phase[k,j,i]

#calculate concentration outside to satisfy nominal concentration
int_phase=int_phase/vol
conc_outside=(conc_nominal-conc_inside*int_phase)/(1.-int_phase)

print 'conc_outside=',conc_outside
print 'conc_inside =',conc_inside
for k in range( nz ) :
  print 'layer ',k
  for j in range( ny ) :
    for i in range( nx ) :
      conc[k,j,i]  = conc_outside

for g in range(n_spheres):
  print 'sphere ',g,', center: ',cx[g],cy[g],cz[g]
  q = quat_inside[g]
  for k in range( nz ) :
    z = k + 0.5
    dz2=distance2_1d(z,cz[g])
    if dz2<r_sq :
      for j in range( ny ) :
        y = j + 0.5
        dy2=distance2_1d(y,cy[g])
        if dy2<r_sq :
          for i in range( nx ) :
            x = i + 0.5
       
            distance_sq = distance2(x,y,z,cx[g],cy[g],cz[g])
            if ( distance_sq < r_sq ) :
              conc[k,j,i]  = conc_inside
              phase[k,j,i] = phase_inside
              for m in range(qlen):
                quat[m,k,j,i] = q[m]

#-----------------------------------------------------------------------
# Write data to file and close

print 'Write data to file'
ncconc.assignValue( conc )
ncphase.assignValue( phase )
for m in range(qlen):
  ncquat[m].assignValue( quat[m,...] )

f.close()


