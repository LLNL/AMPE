# standard packages
import math
import random
import sys
import string
from optparse import OptionParser

# other required packages
import numpy as np
import netCDF4 as nc4

#-----------------------------------------------------------------------
# command-line arguments and options

usage = "Usage: python make_nuclei.py [options] filename"

parser = OptionParser( usage = usage )

parser.add_option( "-x", "--nx", type="int",
                   help="number of cells in x direction" )
parser.add_option( "-y", "--ny", type="int",
                   help="number of cells in y direction" )
parser.add_option( "-z", "--nz", type="int", default=1,
                   help="number of cells in z direction" )
parser.add_option( "-r", "--radius", type="int",
                   help="number of cells for grain radius" )
parser.add_option( "--concentration-in", type="string",
                   help="concentration in interior region" )
parser.add_option( "--concentration-out", type="string",
                   help="concentration in exterior region" )
parser.add_option( "--center", type="string", # something like "12.,14.,23"
                  help="position of grain")

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
  print ("Error: radius is required")
  sys.exit(1)

conc_inside   = options.concentration_in
conc_outside  = options.concentration_out

#-----------------------------------------------------------------------
nspecies=0

tmp = options.concentration_in.split( ',' )
ci = [float(i) for i in tmp]
if nspecies==0:
  nspecies=len(ci)
print ("Composition inside={}".format(ci))

tmp = options.concentration_out.split( ',' )
co = [float(i) for i in tmp]
print ("Composition outside={}".format(co))

print ("nspecies={}".format(nspecies))

#-----------------------------------------------------------------------
def distance2(x1,y1,z1,x2,y2,z2):
  d2=(x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2
  return d2

#-----------------------------------------------------------------------
# Open and define file

f = nc4.Dataset(filename, 'w', format='NETCDF4') 

f.createDimension( 'x', nx )
f.createDimension( 'y', ny )
f.createDimension( 'z', nz )
f.createDimension( 'ns', nspecies )

ncconc = []
for s in range(nspecies):
  c_comp = f.createVariable( 'concentration%d' % s , 'f', ('z','y','x') )
  ncconc.append(c_comp)
ncphase = f.createVariable( 'phase',         'f', ('z','y','x') )

conc  = np.ones( (nspecies,nz,ny,nx), float )
phase = np.zeros( (nz,ny,nx), float )

#-----------------------------------------------------------------------
if ( not (options.center is None) ):
  center = options.center.split( ',' )
else:
  print ("Error: center is required")
  sys.exit(1)

cx=eval(center[0])
cy=eval(center[1])
if len(center)>2:
  cz=eval(center[2])
else:
  cz=11

#-----------------------------------------------------------------------
# Fill data arrays
r_sq = radius**2
print ("grain, center: {},{},{}, radius: {}".format(cx,cy,cz,radius))
for k in range( nz ) :
  z = k + 0.5
  for j in range( ny ) :
    y = j + 0.5
    for i in range( nx ) :
      x = i + 0.5

      #compute distance to center
      distance_sq = distance2(x,y,z,cx,cy,cz)

      #compare with radius
      d = distance_sq - r_sq
      if( d<0. ):
        phase[k,j,i] = 1.

#fill conc values
if nspecies>0:
  print ("set composition for {} species".format(nspecies))
  for s in range(nspecies):
    for k in range( nz ) :
      print("k = {}".format(k))
      for j in range( ny ) :
        for i in range( nx ) :
          phi=phase[k,j,i]
          conc[s,k,j,i] = ci[s]*phi+co[s]*(1.-phi)

#-----------------------------------------------------------------------
# Write data to file and close

print ("Write data to file")
ncphase[:,:,:]=phase

if ( nspecies>0 ):
  for s in range(nspecies):
    ncconc[s][:,:,:]=conc[s,:,:,:]

f.close()
