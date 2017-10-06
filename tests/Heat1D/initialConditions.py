# standard packages
import math
import sys
import string
from optparse import OptionParser

# other required packages
import numpy as N
import netCDF4 as nc4

#-----------------------------------------------------------------------
# command-line arguments and options
usage = "Usage: %prog [options] filename"

parser = OptionParser( usage = usage )

parser.add_option( "-x", "--nx", type="int",
                   help="number of cells in x direction" )
parser.add_option( "-y", "--ny", type="int",
                   help="number of cells in y direction" )
parser.add_option( "-z", "--nz", type="int",
                   help="number of cells in z direction" )

(options, args) = parser.parse_args()

if ( len( args ) < 1 ) :
  parser.error( "filename required argument missing" )
filename = args[0]

nx = options.nx
ny = options.ny
nz = options.nz
if ( not ( nx and ny and nz ) ) :
  print "Error: either -n or all of -x -y -z are required"
  sys.exit(1)

#-----------------------------------------------------------------------
# Open and define file

f = nc4.Dataset(filename, 'w', format='NETCDF4')

f.createDimension( 'x', nx )
f.createDimension( 'y', ny )
f.createDimension( 'z', nz )

nctemperature = f.createVariable( 'temperature', 'f', ('z','y','x') )
temperature = N.zeros( (nz,ny,nx), N.float32 )

#-----------------------------------------------------------------------
# Fill data arrays
h=1./nx

for k in range( nz ) :
  for j in range( ny ) :
    for i in range( nx ) :
      x = (i + 0.5)*h
      y = (j + 0.5)*h
      z = (k + 0.5)*h

      temperature[k,j,i] = 2.*x+math.sin(2.*math.pi*x)+1.

#-----------------------------------------------------------------------
# Write data to file and close
nctemperature[:,:,:]= temperature

f.close()
