# standard packages
import math
import random
import sys
import string
from optparse import OptionParser

# other required packages
import numpy as N
import netCDF4 as nc4

from math import pi

print( sys.path )

#-----------------------------------------------------------------------
# command-line arguments and options

usage = "Usage: %prog [options] filename"

parser = OptionParser( usage = usage )

parser.add_option( "-d", "--dimension", type="int", default=3,
                   help="dimension of subspace containing centers [default: %default]" )
parser.add_option( "-x", "--nx", type="int",
                   help="number of cells in x direction" )
parser.add_option( "-y", "--ny", type="int",
                   help="number of cells in y direction" )
parser.add_option( "-z", "--nz", type="int",
                   help="number of cells in z direction" )
parser.add_option( "--solid-fraction", type="float", default=1./60.,
                  help="solid-fraction at bottom")
parser.add_option( "--concL", type="float", help="conc en phase L")
parser.add_option( "--concA", type="float", help="conc en phase A")

(options, args) = parser.parse_args()

filename = args[0]

nx = options.nx
ny = options.ny
nz = options.nz

nn=[nx,ny,nz]

cl = options.concL
ca = options.concA

print(nn)

sf      = options.solid_fraction

if ( not ( nx and ny and nz ) ) :
  print( "Error: all of -nx -ny -nz are required")
  sys.exit(1)

ndim = options.dimension
if ndim < 3:
  nz=1

#-----------------------------------------------------------------------
# Open and define file

f = nc4.Dataset(filename, 'w', format='NETCDF4')

f.createDimension( 'x', nn[0] )
f.createDimension( 'y', nn[1] )
f.createDimension( 'z', nn[2] )

ncphase = f.createVariable( 'phase', 'f', ('z','y','x') )
phase = N.zeros( (nn[2],nn[1],nn[0]), N.float32 )

ncconc = f.createVariable( 'concentration', 'f', ('z','y','x') )
conc = N.zeros( (nn[2],nn[1],nn[0]), N.float32 )

#-----------------------------------------------------------------------

# Fill data arrays
index = [0,0,0]
invdelta=nn[1]

for j in range( nn[1] ) :
  #get a y in [0,1]
  xx = (j + 0.5)/(1.*nn[1])
  index[1] = j

  #d is negative for the lowest y
  #"sf" fraction of domain)
  d0 = (xx-sf)
  #print("d={}".format(d))

  for k in range( nn[2] ) :
    index[2] = k
    for i in range( nn[0] ) :
      index[0] = i
      phase[index[2],index[1],index[0]] = 0.5*(1.+math.tanh(-0.5*d0*invdelta))
      if( d0<0. ):
        conc[index[2],index[1],index[0]] = ca
      else:
        conc[index[2],index[1],index[0]] = cl

#-----------------------------------------------------------------------
# Write data to file and close

print("Write data to file")
ncphase[:,:,:]=phase[:,:,:]
ncconc[:,:,:]=conc

f.close()
