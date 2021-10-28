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
import netCDF4 as nc4

from math import pi

print( sys.path )

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
parser.add_option( "-c", "--nomconc", type="string",
                   help="nominal concentration" )
parser.add_option( "--concentration-in", type="string", 
                   help="concentration in interior region" )
parser.add_option( "--concentration-out", type="string", 
                   help="concentration in outside region" )
parser.add_option( "--temperature", type="float",
                   help="temperature" )
parser.add_option( "-w", "--width", type="float",
                   help="interface width (in mesh points)", default=5.0 )
parser.add_option( "--solid-fraction", type="float", default=1./60.,
                  help="solid-fraction at bottom")
parser.add_option( "--noise", type="float", default=0.,
                  help="amplitude of noise in y interface")
parser.add_option( "--plane", type="int", default=1,
                  help="direction orthogonal to plane")
parser.add_option( "--quatin", type="string", # something like "1,0,0,0",
                   help="quaternion in grain" )

(options, args) = parser.parse_args()

filename = args[0]


plane = options.plane

nx = options.nx
ny = options.ny
nz = options.nz

nn=[nx,ny,nz]

print(nn)
dir0=0
dir1=1
dir2=2
if plane==0:
  dir1=0
  dir0=1
if plane==2:
  dir1=2
  dir2=1


ngrains = options.ngrains
sf      = options.solid_fraction
widthy  = options.width/nn[1]
noise   = options.noise
quatin  = options.quatin

if ( not ( nx and ny and nz ) ) :
  print( "Error: all of -nx -ny -nz are required")
  sys.exit(1)

ndim = options.dimension
if ndim < 3:
  nz=1

QLEN = options.qlen
print( "qlen={}".format(QLEN))

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
  print("setRandomQinGrains...")
  if QLEN>0:
    
    for g in range(ngrains):

      #pick a discrete angle between 0 and 2*pi
      if ( ngrains == 1 ):
        n = 0
      else:
        n = random.randint(0,nangles-1)
      #print( 'random number =',n
      t = n*h
      print("angle={}".format(t))

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
          #print( u1,u2,u3

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
      print("--- q={}".format(quat_inside[g]))

#fill conc values
nspecies=0
if ( not ( nomconc is None ) ):
  c = map( float, string.split( options.nomconc, ',' ) )
  nspecies=len(c)
  print( "Nominal composition={}".format(c))
if not(conc_inside is None):
  ci = map( float, string.split( options.concentration_in, ',' ) )
  if nspecies==0:
    nspecies=len(ci)
  print( "Composition inside={}".format(ci))
else:
  ci = N.zeros( nspecies, N.float32 )

if ( not ( conc_outside is None ) ):
  co = map( float, string.split( options.concentration_out, ',' ) )
  print( "Composition outside={}".format(co))
else:
  co = N.zeros( nspecies, N.float32 )

print( "nspecies={}".format(nspecies))

if not(nomconc is None):
  for isp in range(nspecies):
    co[isp]=(c[isp]-ci[isp]*sf)/(1.-sf)
    print("Fill composition values")
    print("conc_inside ={}".format(ci[isp]))
    print("conc_outside={}".format(co[isp]))

#-----------------------------------------------------------------------
# Open and define file

#f = NetCDF.NetCDFFile( filename, 'w' )
f = nc4.Dataset(filename, 'w', format='NETCDF4')

f.createDimension( 'x', nn[0] )
f.createDimension( 'y', nn[1] )
f.createDimension( 'z', nn[2] )
if QLEN>0:
  f.createDimension( 'qlen', QLEN )
f.createDimension( 'ns', nspecies )

ncquat=[]
ncconc = []

print( 'Data in single precision...')
if not(temperature is None):
  nctemp  = f.createVariable( 'temperature', 'f', ('z','y','x') )
ncphase = f.createVariable( 'phase',         'f', ('z','y','x') )

for m in range(QLEN):
  name = "quat%d" % (m+1)
  print( name )
  ncquat.append( f.createVariable( name, 'f', ('z','y','x') ) )
for s in range(nspecies):
  c_comp = f.createVariable( 'concentration%d' % s , 'f', ('z','y','x') )
  ncconc.append(c_comp)


phase = N.zeros( (nn[2],nn[1],nn[0]), N.float32 )
if QLEN>0:
  quat  = N.zeros( (QLEN,nn[2],nn[1],nn[0]), N.float32 )
conc  = N.ones( (nspecies,nn[2],nn[1],nn[0]), N.float32 )

#-----------------------------------------------------------------------

# Fill data arrays
xx = [0., 0., 0.]
index = [0,0,0]
fraction=1./ngrains; 
for j in range( nn[dir1] ) :
  #get a y in [0,1]
  xx[dir1] = (j + 0.5)/(1.*nn[dir1])
  index[dir1] = j

  #d is negative for the lowest y
  #"sf" fraction of domain)
  d0 = (xx[dir1]-sf)
  #print("d={}".format(d))
  for k in range( nn[dir2] ) :
    xx[dir2] = k + 0.5
    index[dir2] = k
    for i in range( nn[dir0] ) :
      xx[dir0] = i + 0.5
      index[dir0] = i

      d=d0+noise*random.uniform(-1., 1.)
      
      if( widthy>0. ):
        if( d<0.1 ):
          phase[index[2],index[1],index[0]] = 0.5*(1.+N.tanh(-3.*d/widthy))
      else:
        if( d<0. ):
          phase[index[2],index[1],index[0]] = 1.

      for g in range(1,ngrains):
        dx=abs(xx[dir0]-g*fraction*nn[dir0])
        if( dx<5 ):
          s=N.sin(0.5*pi*dx/5)
          phase[index[2],index[1],index[0]]=phase[index[2],index[1],index[0]]*s*s



#fill quat values
if not quatin is None and ngrains==1:
  q = list(map( float, options.quatin.split(',') ))
  if ( QLEN == 4 ) :
    quat_inside.append( Q.makeNormalizedQuat( q[0], q[1], q[2], q[3] ) )
  else :
    quat_inside.append( Q.makeNormalizedQuat2( q[0], q[1] ) )
else:
  setRandomQinGrains()

offset=0.5*fraction*nn[0]
if QLEN>0:
  gmin=0
  print("Fill quaternion values...")
  for i in range( nn[0] ) :
    xx[0] = i + 0.5

    #select quaternion based on x position
    gmin=0
    for g in range(ngrains):
      dx=abs(xx[0]-g*fraction*nn[dir0]-offset)
      if dx<0.5*fraction*nn[0]:
        gmin=g
    
    print("Plane x={}, grain={}, q={}".format(xx,gmin,quat_inside[gmin]))
    for j in range( nn[1] ) :
      for k in range( nn[2] ) :

        qi=quat_inside[gmin]

        for m in range(QLEN):
          quat[m,k,j,i] = qi[m]

for isp in range(nspecies):
  for x,y in N.nditer([conc[isp,:,:,:],phase], op_flags=['readwrite']):
    x[...]= ci[isp]*y+co[isp]*(1.-y)

#-----------------------------------------------------------------------
# Write data to file and close

print("Write data to file")
if ( nspecies>0 ):
  for s in range(nspecies):
    ncconc[s][:,:,:]=conc[s,:,:,:]
if not(temperature is None):
  nctemp[:,:,:]= temperature
ncphase[:,:,:]=phase
for m in range(QLEN):
  ncquat[m][:,:,:]=quat[m,:,:,:]

f.close()


