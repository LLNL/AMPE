# standard packages
import math
import sys
import string
from optparse import OptionParser

# other required packages
import numpy as N
import netCDF4 as nc4

from math import cos

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

(options, args) = parser.parse_args()

filename = args[0]

Lx = 100.
Ly = 100.
Lz = 100.
c0 = 0.5
epsilon = 0.01

nx = options.nx
ny = options.ny
nz = options.nz

nn=[nx,ny,nz]

print(nn)

if ( not ( nx and ny and nz ) ) :
  print( "Error: all of -nx -ny -nz are required")
  sys.exit(1)

ndim = options.dimension
if ndim < 3:
  nz=1

#-----------------------------------------------------------------------
# Open output file
f = nc4.Dataset(filename, 'w', format='NETCDF4')

f.createDimension( 'x', nn[0] )
f.createDimension( 'y', nn[1] )
f.createDimension( 'z', nn[2] )

s=0
ncconc = f.createVariable( 'concentration%d' % s , 'f', ('z','y','x') )

conc  = N.ones( (nn[2],nn[1],nn[0]), N.float32 )

#-----------------------------------------------------------------------

# Fill data arrays
hx = Lx/(1.*nx)
hy = Ly/(1.*ny)
hz = Lz/(1.*nz)

for i in range( nx ) :
  x = (i + 0.5)*hx

  for j in range( ny ) :
    #get a y in [0,1]
    y = (j + 0.5)*hy

    for k in range(nz):
      z = (k + 0.5)*hz

      t1 = cos(0.105*x)*cos(0.11*y)
      t2 = cos(0.13*x)*cos(0.087*y)
      t3 = cos(0.025*x-0.15*y)*cos(0.07*x-0.02*y)
      
      conc[k,j,i] = c0 + epsilon*(t1+t2*t2+t3)
      if k>0:
        conc[k,j,i] = conc[k,j,i] + epsilon*cos(0.11*z)*cos(0.105*x)

#-----------------------------------------------------------------------
# Write data to file and close

print("Write data to file")
ncconc[:,:,:]=conc

f.close()
