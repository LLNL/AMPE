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
import sys
import string
from optparse import OptionParser

# other required packages
import numpy as N
import netCDF4 as nc4

from math import sqrt

print (sys.path)

#-----------------------------------------------------------------------
# command-line arguments and options

usage = "Usage: %prog [options] filename"

parser = OptionParser( usage = usage )

parser.add_option( "-x", "--nx", type="int",
                   help="number of cells in x direction" )
parser.add_option( "-y", "--ny", type="int",
                   help="number of cells in y direction" )
parser.add_option( "-z", "--nz", type="int", default=1,
                   help="number of cells in z direction" )
parser.add_option( "-r", "--radius", type="int",
                   help="number of cells for sphere radius [default: nx/4]" )

parser.add_option("--periodic", type="string", # something like "1,0,1" for periodic in x and z
                  default="1,1,1",
                  help="specify which dimensions are periodic") 

(options, args) = parser.parse_args()

filename = args[0]

nx = options.nx
ny = options.ny
nz = options.nz

if ( not ( nx and ny and nz ) ) :
  print ("Error: either -n or all of -nx -ny -nz are required")
  sys.exit(1)


radius = options.radius
if radius is None :
  radius = nx / 10
  print ("Radius not specified: using nx/n_spheres = {}".format(radius))

periodic = options.periodic.split( ',')

#-----------------------------------------------------------------------

#distance square function for periodic bc
if periodic[0]==1:
  rangex=range(-1,2)
else:
  rangex=range(1)

if periodic[1]==1:
  rangey=range(-1,2)
else:
  rangey=range(1)

if periodic[2]==1:
  rangez=range(-1,2)
else:
  rangez=range(1)

#distance square function for periodic bc
def distance2_1d_x(x1,x2):
  d2=(x1 - x2)**2
  for ir in rangex:
    l2=(x1-x2+ir*nx)**2
    if l2<d2:
      d2=l2
  return sqrt(d2)

def distance2_1d_y(y1,y2):
  d2=(y1 - y2)**2
  for ir in rangey:
    l2=(y1-y2+ir*ny)**2
    if l2<d2:
      d2=l2
  return sqrt(d2)

def distance2_1d_z(z1,z2):
  d2=(z1 - z2)**2
  for ir in rangez:
    l2=(z1-z2+ir*nz)**2
    if l2<d2:
      d2=l2
  return sqrt(d2)

#-----------------------------------------------------------------------
# Open and define file

f = nc4.Dataset(filename, 'w', format='NETCDF4') 

f.createDimension( 'x', nx )
f.createDimension( 'y', ny )
f.createDimension( 'z', nz )

ncconc = f.createVariable( 'concentration0',  'f', ('z','y','x') )
ncphase = f.createVariable( 'phase',  'f', ('z','y','x') )

conc  = N.ones( (nz,ny,nx), N.float32 )
phase = N.zeros( (nz,ny,nx), N.float32 )

#-----------------------------------------------------------------------
cx=nx/2
cy=ny/2
cz=nz/2

#-----------------------------------------------------------------------
# Fill data arrays
for k in range( nz ) :
  for j in range( ny ) :
    for i in range( nx ) :
      phase[k,j,i] = 0.

#fill phase value
for k in range( nz ) :
  z = k + 0.5
  dz2=distance2_1d_z(z,cz)
  if dz2<radius :
    for j in range( ny ) :
      y = j + 0.5
      dy2=distance2_1d_y(y,cy)
      if dy2<radius :
        for i in range( nx ) :
          x = i + 0.5

          #compute distance to center
          dx2=distance2_1d_x(x,cx)

          if( dx2<radius ):
            phase[k,j,i] = 1.

#fill conc values
for k in range( nz ) :
  for j in range( ny ) :
    for i in range( nx ) :
       conc[k,j,i] = 1.*phase[k,j,i]

#-----------------------------------------------------------------------
# Write data to file and close

print ("Write data to file")
ncphase[:,:,:]=phase
ncconc[:,:,:]=conc

f.close()
