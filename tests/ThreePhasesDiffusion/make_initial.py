# standard packages
import math
import random
import sys
import string
from optparse import OptionParser

# other required packages
import numpy as N
import netCDF4 as nc4
from math import exp

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
parser.add_option( "--concB", type="float", help="conc en phase B")

(options, args) = parser.parse_args()

filename = args[0]

nx = options.nx
ny = options.ny
nz = options.nz

nn=[nx,ny,nz]

cb = options.concB

print(nn)

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

ncphase=[]
ncphase.append( f.createVariable( 'phase0', 'f', ('z','y','x') ))
ncphase.append( f.createVariable( 'phase1', 'f', ('z','y','x') ))
ncphase.append( f.createVariable( 'phase2', 'f', ('z','y','x') ))

phase = N.zeros( (3, nn[2],nn[1],nn[0]), N.float32 )

ncconc = f.createVariable( 'concentration', 'f', ('z','y','x') )
conc = N.zeros( (nn[2],nn[1],nn[0]), N.float32 )

#-----------------------------------------------------------------------
# Fill data arrays
for k in range( nn[2] ) :
  for j in range( nn[1] ) :
    for i in range( nn[0] ) :

      phase[2,k,j,i] = 1.
      conc[k,j,i] = cb*exp(-0.01*(i*i+j*j+k*k))

#-----------------------------------------------------------------------
# Write data to file and close
print("Write data to file")
for m in range(3):
  ncphase[m][:,:,:]=phase[m,:,:,:]
ncconc[:,:,:]=conc

f.close()
