# standard packages
import math
import sys
import string
from optparse import OptionParser

# my required packages
import quat as Q

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
parser.add_option( "-q", "--qlen", type="int", default=4,
                   help="number of component for q [default: %default]" )

(options, args) = parser.parse_args()

filename = args[0]

n_spheres = 4

nx = options.nx
ny = options.ny
nz = options.nz

QLEN = options.qlen
print ("qlen={}".format(QLEN))

radius = nx/2

#-----------------------------------------------------------------------
# generate 4 quaternions
quat_inside = []
q0 = 0.5
q1 = 0.5
q2 = 0.5
q3 = 0.5
quat_inside.append( Q.makeQuat( q0, q1, q2, q3 ) )

q0 = 0.5
q1 = 0.5
q2 = 0.5
q3 = 0.5
quat_inside.append( Q.makeQuat( q0, q1, q2, q3 ) )

q0 = -0.5
q1 = 0.5
q2 = 0.5
q3 = 0.5
quat_inside.append( Q.makeQuat( q0, q1, q2, q3 ) )

q0 = 0.5
q1 = -0.5
q2 = 0.5
q3 = 0.5
quat_inside.append( Q.makeQuat( q0, q1, q2, q3 ) )

#-----------------------------------------------------------------------
def distance2(x1,y1,z1,x2,y2,z2):
  d2=(x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2
  return d2

#-----------------------------------------------------------------------
def distancemax(x1,y1,z1,x2,y2,z2):
  d=(x1 - x2)
  if (y1 - y2)>d:
    d=(y1 - y2)
  if (z1 - z2)>d:
    d=(z1 - z2)
  return d

#-----------------------------------------------------------------------
# Open and define file
f = nc4.Dataset(filename, 'w', format='NETCDF4') 

f.createDimension( 'x', nx )
f.createDimension( 'y', ny )
f.createDimension( 'z', nz )
if QLEN>0:
  f.createDimension( 'qlen', QLEN )

ncphase = f.createVariable( 'phase',         'f', ('z','y','x') )

ncquat=[]
for m in range(QLEN):
  name = "quat%d" % (m+1)
  print (name)
  ncquat.append( f.createVariable( name, 'f', ('z','y','x') ) )

phase = N.zeros( (nz,ny,nx), N.float32 )
quat  = N.zeros( (QLEN,nz,ny,nx), N.float32 )

#-----------------------------------------------------------------------

# generate random (x,y) positions
cx=[]
cy=[]
cz=[]

# center lower left
cx.append(0.)
cy.append(0.)
cz.append(0.5)

# center upper left
cx.append(0.)
cy.append(ny)
cz.append(0.5)

# center lower right
cx.append(nx)
cy.append(0.)
cz.append(0.5)

# center upper right
cx.append(nx)
cy.append(ny)
cz.append(0.5)


#-----------------------------------------------------------------------
# Fill data arrays

#fill phase value
for g in range(n_spheres):
  print ("sphere {}, center: {},{},{}, radius: {}".format(g,cx[g],cy[g],cz[g],radius))
  for k in range( nz ) :
    z = k + 0.5
    for j in range( ny ) :
      y = j + 0.5
      for i in range( nx ) :
        x = i + 0.5

        #compute distance to center
        distance_sq = distance2(x,y,z,cx[g],cy[g],cz[g])

        #compare with radius
        d = N.sqrt(distance_sq) - radius
        if d<=0.:
          phase[k,j,i] = 1.

print ("Fill quaternion values...")
for i in range( nx ) :
  x = i + 0.5
  print ("Plane x={}".format(x))
  for j in range( ny ) :
    y = j + 0.5
    for k in range( nz ) :
      z = k + 0.5
      qi=[]
      for g in range(n_spheres):
        distance = distancemax(x,y,z,cx[g],cy[g],cz[g])
        if distance<=radius:
          q = quat_inside[g]
          for m in range(QLEN):
            quat[m,k,j,i] = q[m]
          break

#-----------------------------------------------------------------------
# Write data to file and close

print ("Write data to file")
ncphase[:,:,:]=phase
for m in range(QLEN):
  ncquat[m][:,:,:]=quat[m,:,:,:]

f.close()
