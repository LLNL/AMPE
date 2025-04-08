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
parser.add_option( "--concB", type="float", help="conc en phase B")

parser.add_option( "--plane", type="int", default=1,
                  help="direction orthogonal to plane")

(options, args) = parser.parse_args()

filename = args[0]

plane = options.plane

nx = options.nx
ny = options.ny
nz = options.nz

nn=[nx,ny,nz]

cl = options.concL
ca = options.concA
cb = options.concB

print(nn)
dir0=0
dir1=1
dir2=2
if plane==0:
  dir1=0
  dir0=1

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

ncphase=[]
ncphase.append( f.createVariable( 'phase0', 'f', ('z','y','x') ))
ncphase.append( f.createVariable( 'phase1', 'f', ('z','y','x') ))
ncphase.append( f.createVariable( 'phase2', 'f', ('z','y','x') ))

phase = N.zeros( (3, nn[2],nn[1],nn[0]), N.float32 )

ncconc = f.createVariable( 'concentration', 'f', ('z','y','x') )
conc = N.zeros( (nn[2],nn[1],nn[0]), N.float32 )

#-----------------------------------------------------------------------

# Fill data arrays
index = [0,0,0]

for j in range( nn[dir1] ) :
  #get a y in [0,1]
  xx = (j + 0.5)/(1.*nn[dir1])
  index[dir1] = j

  #d is negative for the lowest y
  #"sf" fraction of domain)
  d0 = (xx-sf)
  #print("d={}".format(d))

  for k in range( nn[dir2] ) :
    index[dir2] = k
    for i in range( nn[dir0] ) :
      index[dir0] = i

      if( d0<0. ):
        if i<nn[dir0]/2:
          phase[1,index[2],index[1],index[0]] = 1.
          conc[index[2],index[1],index[0]] = ca
        else:
          phase[2,index[2],index[1],index[0]] = 1.
          conc[index[2],index[1],index[0]] = cb
      else:
        phase[0,index[2],index[1],index[0]] = 1.
        conc[index[2],index[1],index[0]] = cl

#-----------------------------------------------------------------------
# Write data to file and close

print("Write data to file")
for m in range(3):
  ncphase[m][:,:,:]=phase[m,:,:,:]
ncconc[:,:,:]=conc

f.close()
