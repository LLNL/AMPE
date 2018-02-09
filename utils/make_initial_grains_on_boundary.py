# standard packages
import math
import random
import sys
import string
from optparse import OptionParser

# my required packages
import quat as Q

# other required packages
import numpy as N
#from Scientific.IO import NetCDF
#from scipy.io import netcdf as NetCDF
import netCDF4 as nc4

from math import pi

print sys.path

#-----------------------------------------------------------------------
# command-line arguments and options

usage = "Usage: %prog [options] filename"

parser = OptionParser( usage = usage )

parser.add_option( "--ngrains", type="int", default=1,
                   help="number of grains to nucleate [default: %default]" )
parser.add_option( "-d", "--dimension", type="int", default=3,
                   help="dimension of subspace containing centers [default: %default]" )
parser.add_option( "-x", "--nx", type="int",
                   help="number of cells in x direction" )
parser.add_option( "-y", "--ny", type="int",
                   help="number of cells in y direction" )
parser.add_option( "-z", "--nz", type="int",
                   help="number of cells in z direction" )
parser.add_option( "-q", "--qlen", type="int", default=0,
                   help="number of component for q [default: %default]" )
parser.add_option( "-c", "--nomconc", type="float", 
                   help="nominal concentration" )
parser.add_option( "--concentration-in", type="float", 
                   help="concentration in interior region" )
parser.add_option( "--concentration-out", type="float", 
                   help="concentration in outside region" )
parser.add_option( "--temperature", type="float",
                   help="temperature" )
parser.add_option( "-w", "--width", type="float",
                   help="interface width (in mesh points)", default=5.0 )
parser.add_option( "--solid-fraction", type="float", default=1./60.,
                  help="solid-fraction at bottom")
parser.add_option( "--noise", type="float", default=0.,
                  help="amplitude of noise in y interface")
parser.add_option( "--double", action="store_true", dest="double_precision", 
                  default=False)

(options, args) = parser.parse_args()

filename = args[0]

double_precision = options.double_precision
if double_precision:
  print "use double precision..."

nx = options.nx
ny = options.ny
nz = options.nz

ngrains = options.ngrains
sf      = options.solid_fraction
widthy  = options.width/ny
noise   = options.noise

if ( not ( nx and ny and nz ) ) :
  print "Error: all of -nx -ny -nz are required"
  sys.exit(1)

ndim = options.dimension
if ndim < 3:
  nz=1

QLEN = options.qlen
print "qlen=",QLEN

nomconc       = options.nomconc
conc_inside   = options.concentration_in
conc_outside  = options.concentration_out
if conc_inside is None :
  conc_inside = nomconc

temperature = options.temperature

# generate quaternions corresponding to random orientations
random.seed( 112345 )
quat_inside=[]
nangles = 20. #actually means fewer different angles if symmetry is on
h = 2.*math.pi/nangles

h1=0.1
h2=0.5*h1
h3=h2

#-----------------------------------------------------------------------
def setRandomQinGrains():
  print 'setRandomQinGrains...'
  if QLEN>0:
    
    for g in range(ngrains):

      #pick a discrete angle between 0 and 2*pi
      if ( ngrains == 1 ):
        n = 0
      else:
        n = random.randint(0,nangles-1)
      #print 'random number =',n
      t = n*h
      print 'angle=',t

      if ( QLEN == 1 ) :

        q0 = t
        q1 = 0.
        q2 = 0.
        q3 = 0.

      else :
        if ( QLEN == 4 ) :
          #K. Shoemake. "Uniform random rotations."
          #In D. Kirk, editor, Graphics Gems III, pages 124-132. Academic, New York, 1992. 
          u1 = random.uniform(0, 1.)
          u2 = random.uniform(0, 1.)
          u3 = random.uniform(0, 1.)
          #print u1,u2,u3

          #restricts possible orientations
          n  = math.floor( u1/h1 )
          u1 = n*h1
          n  = math.floor( u2/h2 )
          u2 = n*h2
          n  = math.floor( u3/h3 )
          u3 = n*h3

          w = math.sqrt(1.-u1)*math.sin(2.*math.pi*u2)
          x = math.sqrt(1.-u1)*math.cos(2.*math.pi*u2)
          y = math.sqrt(u1)*math.sin(2.*math.pi*u3)
          z = math.sqrt(u1)*math.cos(2.*math.pi*u3)
 
        else : #QLEN=2
          w = math.cos( t )
          x = math.sin( t)
          y = 0.
          z = 0.

        q0 = w
        q1 = x
        q2 = y
        q3 = z

      quat_inside.append( Q.makeQuat( q0, q1, q2, q3 ) )
      print '--- q=',quat_inside[g]

#-----------------------------------------------------------------------
# Open and define file

#f = NetCDF.NetCDFFile( filename, 'w' )
f = nc4.Dataset(filename, 'w', format='NETCDF4')

f.createDimension( 'x', nx )
f.createDimension( 'y', ny )
f.createDimension( 'z', nz )
if QLEN>0:
  f.createDimension( 'qlen', QLEN )

ncquat=[]

if double_precision:
  print 'Data in double precision...'
  if not(nomconc is None) or not(conc_inside is None):
    ncconc  = f.createVariable( 'concentration', 'd', ('z','y','x') )
  if not(temperature is None):
    nctemp  = f.createVariable( 'temperature', 'd', ('z','y','x') )
    temperature_field  = N.ones( (nz,ny,nx), N.float64 )
  ncphase = f.createVariable( 'phase',         'd', ('z','y','x') )

  for m in range(QLEN):
    name = "quat%d" % (m+1)
    print name
    ncquat.append( f.createVariable( name, 'd', ('z','y','x') ) )

  conc  = N.ones( (nz,ny,nx), N.float64 )
  phase = N.zeros( (nz,ny,nx), N.float64 )
  if QLEN>0:
    quat  = N.zeros( (QLEN,nz,ny,nx), N.float64 )
else:
  print 'Data in single precision...'
  if not(nomconc is None) or not(conc_inside is None):
    ncconc  = f.createVariable( 'concentration', 'f', ('z','y','x') )
  if not(temperature is None):
    nctemp  = f.createVariable( 'temperature', 'f', ('z','y','x') )
    temperature_field  = N.ones( (nz,ny,nx), N.float32 )
  ncphase = f.createVariable( 'phase',         'f', ('z','y','x') )

  for m in range(QLEN):
    name = "quat%d" % (m+1)
    print name
    ncquat.append( f.createVariable( name, 'f', ('z','y','x') ) )

  conc  = N.ones( (nz,ny,nx), N.float32 )
  phase = N.zeros( (nz,ny,nx), N.float32 )
  if QLEN>0:
    quat  = N.zeros( (QLEN,nz,ny,nx), N.float32 )

#-----------------------------------------------------------------------

# Fill data arrays

vol=nx*ny*nz
vs=0.
vl=0.;

shift=[]
for i in range( nx ) :
  shift.append(noise*random.uniform(-1., 1.))
  #shift.append(0.001*random.uniform(-1., 1.))
  #shift.append(0.0)
 
fraction=1./ngrains; 
for k in range( nz ) :
  z = k + 0.5
  for j in range( ny ) :
    #get a y in [0,1]
    y = (j + 0.5)/(1.*ny)
    
    #d is negative for the lowest y 
    #(1/60 fraction of domain)
    #d = (6.*y)-0.1
    d = (y-sf)
    #print 'd=',d
    for i in range( nx ) :
      x = i + 0.5
      
      d=d+shift[i]
      
      if( widthy>0. ):
        if( d<0.1 ):
          phase[k,j,i] = 0.5*(1.+N.tanh(-3.*d/widthy))
      else:
        if( d<0. ):
          phase[k,j,i] = 1.

      for g in range(1,ngrains):
        dx=abs(x-g*fraction*nx)
        if( dx<5 ):
          s=N.sin(0.5*pi*dx/5)
          phase[k,j,i]=phase[k,j,i]*s*s

      vs=vs+phase[k,j,i]

vl=vol-vs

#fill quat values
setRandomQinGrains()

offset=0.5*fraction*nx
if QLEN>0:
  gmin=0
  print 'Fill quaternion values...'
  for i in range( nx ) :
    x = i + 0.5

    #select quaternion based on x position
    gmin=0
    for g in range(ngrains):
      dx=abs(x-g*fraction*nx-offset)
      if dx<0.5*fraction*nx:
        gmin=g
    
    print 'Plane x=',x,', grain=',gmin,', q=',quat_inside[gmin]
    for j in range( ny ) :
      for k in range( nz ) :

        qi=quat_inside[gmin]

        for m in range(QLEN):
          quat[m,k,j,i] = qi[m]

#fill conc values
if nomconc is None:
  if not(conc_inside is None):
    nomconc=(conc_inside*vs+conc_outside*vl)/vol


if not(nomconc is None):
  conc_outside=(nomconc*vol-conc_inside*vs)/vl
  print 'Fill composition values'
  print 'conc_inside =',conc_inside
  print 'conc_outside=',conc_outside

  for x,y in N.nditer([conc,phase], op_flags=['readwrite']):
    x[...]= conc_inside*y+conc_outside*(1.-y)

#-----------------------------------------------------------------------
# Write data to file and close

print 'Write data to file'
if not(nomconc is None):
  ncconc[:,:,:]= conc
if not(temperature is None):
  nctemp[:,:,:]= temperature
ncphase[:,:,:]=phase
for m in range(QLEN):
  ncquat[m][:,:,:]=quat[m,:,:,:]

f.close()


