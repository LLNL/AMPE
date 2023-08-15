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
parser.add_option( "-w", "--width", type="float",
                   help="interface width (in mesh points)", default=5.0 )
parser.add_option( "--solid-fraction", type="float", default=1./60.,
                  help="solid-fraction at bottom")
parser.add_option( "--noise", type="float", default=0.,
                  help="amplitude of noise in y interface")
parser.add_option( "--fluctuation", type="float", default=0.,
                  help="amplitude of fluctuation in y interface")
parser.add_option( "--plane", type="int", default=1,
                  help="direction orthogonal to plane")
parser.add_option( "--quatin", type="string", # something like "1,0,0,0",
                   help="quaternion in grain" )
parser.add_option( "--periodicx", type="int", default=0,
                  help="periodic in direction x")
parser.add_option( "--periodicy", type="int", default=0,
                  help="periodic in direction y")
parser.add_option( "--periodicz", type="int", default=0,
                  help="periodic in direction z")
parser.add_option( "--smooth", type="int", default=1,
                  help="smooth quaternion field")

(options, args) = parser.parse_args()

filename = args[0]


plane = options.plane
smooth_quat = options.smooth

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
widthy  = options.width/nn[dir1]
noise   = options.noise
print( "noise={}".format(noise))
fluctuation = options.fluctuation
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

# generate quaternions corresponding to random orientations
random.seed( 112345 )
quat_inside=[]
nangles = 20. #actually means fewer different angles if symmetry is on
h = 2.*math.pi/nangles

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
        [w,x,y,z] = Q.getQuatRandom(t,QLEN)

        q0 = w
        q1 = x
        q2 = y
        q3 = z

      quat_inside.append( Q.makeQuat( q0, q1, q2, q3 ) )
      print("--- q={}".format(quat_inside[g]))

#-----------------------------------------------------------------------
def smoothQuat(quat,px,py,pz):
  print("smoothQuat...")
  if(px):
    ix=0
    nx=nn[0]
  else:
    ix=1
    nx=nn[0]-1
  if(py):
    iy=0
    ny=nn[1]
  else:
    iy=1
    ny=nn[1]-1
  if(pz or nn[2]==1):
    iz=0
    nz=nn[2]
  else:
    iz=1
    nz=nn[2]-1
  print("ix={}, nx={}".format(ix,nx))
  print("iy={}, ny={}".format(iy,ny))
  print("iz={}, nz={}".format(iz,nz))
  quat_new = N.zeros( (nn[2],nn[1],nn[0]), N.float32 )
  quat_new[:,:,:]=quat[:,:,:]
  for i in range(ix,nx):
    for j in range(iy,ny):
      for k in range(iz,nz):
        quat_new[k,j,i] = 0.25*quat[k,j,i] \
                        + 0.125*(quat[(k-1)%nz,j,i]) \
                        + 0.125*(quat[(k+1)%nz,j,i]) \
                        + 0.125*(quat[k,(j-1)%ny,i]) \
                        + 0.125*(quat[k,(j+1)%ny,i]) \
                        + 0.125*(quat[k,j,(i-1)%nx]) \
                        + 0.125*(quat[k,j,(i+1)%nx])
  quat[:,:,:]=quat_new[:,:,:]



#fill conc values
nspecies=0
if ( not ( nomconc is None ) ):
  c = list(map( float, options.nomconc.split(',' ) ))
  nspecies=len(c)
  print( "Nominal composition={}".format(c))
if not(conc_inside is None):
  ci = list(map( float, conc_inside.split(',' ) ))
  if nspecies==0:
    nspecies=len(ci)
  print( "Composition inside={}".format(ci))
else:
  ci = N.zeros( nspecies, N.float32 )

if ( not ( conc_outside is None ) ):
  co = list(map( float, options.concentration_out.split(',' ) ))
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
# Open output file
f = nc4.Dataset(filename, 'w', format='NETCDF4')

f.createDimension( 'x', nn[0] )
f.createDimension( 'y', nn[1] )
f.createDimension( 'z', nn[2] )
if QLEN>0:
  f.createDimension( 'qlen', QLEN )
f.createDimension( 'ns', nspecies )

ncquat=[]
ncconc = []

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

for k in range( nn[dir2] ) :
  xx[dir2] = (k + 0.5)/(1.*nn[dir2])
  index[dir2] = k

  for i in range( nn[dir0] ) :
    xx[dir0] = (i + 0.5)/(1.*nn[dir0])
    index[dir0] = i

    delta = 0.
    if noise>0.:
      delta = noise*random.uniform(-1., 1.)
    else:
      if fluctuation>0.:
        delta = fluctuation*math.sin(xx[dir0]*2.*math.pi)

    for j in range( nn[dir1] ) :
      #get a y in [0,1]
      xx[dir1] = (j + 0.5)/(1.*nn[dir1])
      index[dir1] = j

      #d is negative for the lowest y
      #"sf" fraction of domain)
      d = (xx[dir1]-sf)+delta
      #print("d={}".format(d))

      if( widthy>0. ):
        if( d<0.1 ):
          phase[index[2],index[1],index[0]] = 0.5*(1.+N.tanh(-3.*d/widthy))
      else:
        if( d<0. ):
          phase[index[2],index[1],index[0]] = 1.

      for g in range(1,ngrains):
        dx=abs(xx[dir0]-g*fraction)
        if( nn[dir0]*dx<5. ):
          s=N.sin(0.5*pi*dx*nn[dir0]/5)
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
  for i in range( nn[dir0] ) :
    xx[0] = i + 0.5
    index[dir0] = i

    #select quaternion based on x position
    gmin=0
    for g in range(ngrains):
      dx=abs(xx[0]-g*fraction*nn[dir0]-offset)
      if dx<0.5*fraction*nn[dir0]:
        gmin=g
    
    print("Plane x={}, grain={}, q={}".format(xx,gmin,quat_inside[gmin]))
    for j in range( nn[dir1] ) :
      index[dir1] = j
      for k in range( nn[dir2] ) :
        index[dir2] = k

        qi=quat_inside[gmin]

        for m in range(QLEN):
          quat[m,index[2],index[1],index[0]] = qi[m]

for isp in range(nspecies):
  for x,y in N.nditer([conc[isp,:,:,:],phase], op_flags=['readwrite']):
    x[...]= ci[isp]*y+co[isp]*(1.-y)

if smooth_quat:
  for it in range(3):
    for m in range(QLEN):
      smoothQuat(quat[m,:,:,:], options.periodicx, options.periodicy, options.periodicz)
else:
  #set q to random quaternion in liquid phase
  for i in range( nn[dir0] ) :
    index[dir0] = i
    for j in range( nn[dir1] ) :
      index[dir1] = j
      for k in range( nn[dir2] ) :
        index[dir2] = k
        if ( phase[index[2],index[1],index[0]]<0.1 ) :
          n = random.randint(0,nangles-1)
          #print( 'random number =',n
          t = n*h
          qq=Q.getQuatRandom(t,QLEN)
          for m in range(QLEN):
            quat[m,index[2],index[1],index[0]] = qq[m]

#-----------------------------------------------------------------------
# Write data to file and close

print("Write data to file")
if ( nspecies>0 ):
  for s in range(nspecies):
    ncconc[s][:,:,:]=conc[s,:,:,:]
ncphase[:,:,:]=phase
for m in range(QLEN):
  ncquat[m][:,:,:]=quat[m,:,:,:]

f.close()
