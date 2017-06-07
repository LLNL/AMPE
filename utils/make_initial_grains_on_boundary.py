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
from Scientific.IO import NetCDF
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
parser.add_option( "-i", "--phase-in", type="float", default=1.0,
                   help="phase in interior region [default: %default]" )
parser.add_option( "-o", "--phase-out", type="float", default=0.0,
                   help="phase in exterior region [default: %default]" )
parser.add_option( "-c", "--nomconc", type="float", 
                   help="nominal concentration" )
parser.add_option( "--concentration-in", type="float", 
                   help="concentration in interior region" )
parser.add_option( "--concentration-out", type="float", 
                   help="concentration in outside region" )
parser.add_option( "--temperature", type="float",
                   help="temperature" )
parser.add_option( "-w", "--width", type="float",
                   help="interface width", default=0.0 )
parser.add_option( "--double", action="store_true", dest="double_precision", default=False)

(options, args) = parser.parse_args()

filename = args[0]

double_precision = options.double_precision
if double_precision:
  print "use double precision..."

nx = options.nx
ny = options.ny
nz = options.nz

ngrains = options.ngrains

width = 0.
if ( options.width ) : width = options.width

if ( not ( nx and ny and nz ) ) :
  print "Error: all of -nx -ny -nz are required"
  sys.exit(1)


ndim = options.dimension
if ndim < 3:
  nz=1

QLEN = options.qlen
print "qlen=",QLEN

phase_inside  = options.phase_in
phase_outside = options.phase_out

nomconc       = options.nomconc
conc_inside   = options.concentration_in
conc_outside  = options.concentration_out
if conc_inside is None :
  conc_inside = nomconc

temperature = options.temperature

# generate quaternions corresponding to random orientations
random.seed( 11234 )
quat_inside=[]
nangles = 20.
h = 2.*math.pi/nangles

h1=0.1
h2=0.5*h1
h3=h2

#-----------------------------------------------------------------------
def setRandomQinGrains():
  print 'setRandomQinGrains...'
  if QLEN>0:
    if ( QLEN == 4 ) :
      qref = Q.makeNormalizedQuat() # (1,0,0,0)
    else:
      qref = Q.makeNormalizedQuat2()
    
    for g in range(ngrains):

      #pich a discrete angle between -pi and pi
      t = math.pi * random.uniform(-1, 1.)
      n = math.floor( t/h )
      t = n*h

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

f = NetCDF.NetCDFFile( filename, 'w' )

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
print 'Fill Phase values'
if ( options.phase_out ) :
  for k in range( nz ) :
    for j in range( ny ) :
      for i in range( nx ) :
        phase[k,j,i] = phase_outside

vol=nx*ny*nz
vs=0.
vl=0.;

shift=[]
for i in range( nx ) :
  shift.append(0.001*random.uniform(-1., 1.))
  #shift.append(0.0)
  
for k in range( nz ) :
  z = k + 0.5
  for j in range( ny ) :
    y = j + 0.5
    d = (6.*y)/(1.*ny)-0.1
    #print 'd=',d
    for i in range( nx ) :
      x = i + 0.5
      
      d=d+shift[i]
      
      if( width>0. ):
        sq=abs(d)
        if( sq<0.15 ):
          phase[k,j,i] = 0.5*(1.+N.tanh(-0.5*d/width))
          vs=vs+phase[k,j,i]
      else:
        if( d<0.1 ):
          phase[k,j,i] = phase_inside
          vs=vs+phase[k,j,i]
      
      dx=abs(x-0.5*nx)
      if( dx<5 ):
        s=N.sin(0.5*pi*dx/5)
        phase[k,j,i]=phase[k,j,i]*s*s


vl=vol-vs

#fill quat values
setRandomQinGrains()

if QLEN>0:
  gmin=0
  print 'Fill quaternion values...'
  for i in range( nx ) :
    x = i + 0.5

    #select quaternion based on x position
    if( x<=0.5*nx or ngrains==1):
      gmin=0
    else:
      gmin=1
    
    print 'Plane x=',x
    for j in range( ny ) :
      y = j + 0.5
      for k in range( nz ) :
        z = k + 0.5

        qi=quat_inside[gmin]

        for m in range(QLEN):
          quat[m,k,j,i] = qi[m]

  for i in range( nx ) :
    for j in range( ny ) :
      for k in range( nz ) :
        for m in range(QLEN):
          quat[m,k,j,i] = quat[m,k,j,i]

#fill conc values
if nomconc is None:
  if not(conc_inside is None):
    nomconc=(conc_inside*vs+conc_outside*vl)/vol


if not(nomconc is None):
  conc_outside=(nomconc*vol-conc_inside*vs)/vl
  print 'Fill composition values'
  print 'conc_inside =',conc_inside
  print 'conc_outside=',conc_outside

  for g in range(ngrains):
    for k in range( nz ) :
      for j in range( ny ) :
        for i in range( nx ) :
          conc[k,j,i]  = conc_inside*phase[k,j,i]+conc_outside*(1.-phase[k,j,i])

#-----------------------------------------------------------------------
# Write data to file and close

print 'Write data to file'
if not(nomconc is None):
  ncconc.assignValue( conc )
print 'set T to ',temperature
if not(temperature is None):
  nctemp.assignValue( temperature )
ncphase.assignValue( phase )
for m in range(QLEN):
  ncquat[m].assignValue( quat[m,...] )

f.close()


