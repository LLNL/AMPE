# Copyright (c) 2018, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory
# Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
# LLNL-CODE-747500
# All rights reserved.
# This file is part of AMPE. 
# For details, see https://github.com/LLNL/AMPE
# Please also read AMPE/LICENSE.
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# - Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the disclaimer below.
# - Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the disclaimer (as noted below) in the
#   documentation and/or other materials provided with the distribution.
# - Neither the name of the LLNS/LLNL nor the names of its contributors may be
#   used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
# LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
# 
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
from scipy.io import netcdf as NetCDF

import netCDF4 as nc4

print sys.path

#-----------------------------------------------------------------------
# command-line arguments and options

usage = "Usage: %prog [options] filename"

parser = OptionParser( usage = usage )

parser.add_option( "--ngrains", type="int", default=1,
                   help="number of grains to nucleate [default: %default]" )
parser.add_option( "-n", "--ncells", type="int",
                   help="number of cells in x,y,z (will all be equal)" )
parser.add_option( "-x", "--nx", type="int",
                   help="number of cells in x direction" )
parser.add_option( "-y", "--ny", type="int",
                   help="number of cells in y direction" )
parser.add_option( "-z", "--nz", type="int", default=1,
                   help="number of cells in z direction" )
parser.add_option( "-q", "--qlen", type="int", default=0,
                   help="number of component for q [default: %default]" )
parser.add_option( "-r", "--radius", type="int",
                   help="number of cells for sphere radius [default: nx/4]" )
parser.add_option( "-i", "--phase-in", type="float", default=1.0,
                   help="phase in interior region [default: %default]" )
parser.add_option( "-o", "--phase-out", type="float", default=0.0,
                   help="phase in exterior region [default: %default]" )
parser.add_option( "-c", "--nomconc", type="string", # something like "0.1,0.2",
                   help="nominal concentration" )
parser.add_option( "--concentration-in", type="string",
                   help="concentration in interior region" )
parser.add_option( "--concentration-out", type="string",
                   help="concentration in exterior region" )
parser.add_option( "--temperature-in", type="float",
                   help="temperature in interior region" )
parser.add_option( "--temperature-out", type="float",
                   help="temperature in outside region" )
parser.add_option( "--quat-out", type="string", default=None,
                   help="quat in outside region" )
parser.add_option( "--jitter-factor", type="float", default=1.0, 
                   help="jitter factor" )
parser.add_option( "-s", "--crystal_sym", action="store_true", dest="crystal_sym", default=False)
parser.add_option( "-b", "--bain", action="store_true", dest="bain", default=False)
parser.add_option( "-3", "--three", action="store_true", dest="three", default=False)
parser.add_option( "-w", "--width", type="float",
                   help="interface width", default=0.0 )

parser.add_option( "--symmetry-test", action="store_true",
                   default=False,
                   help="rotate each octant by a random symmetry rotation" )
parser.add_option( "--double", action="store_true", dest="double_precision", default=False)
parser.add_option( "--center0", type="string", # something like "12.,14.,23"
                  help="position of center 0")
parser.add_option("--periodic", type="string", # something like "1,0,1" for periodic in x and z
                  default="1,1,1",
                  help="specify which dimensions are periodic") 
parser.add_option("--quat0", type="string",
                  help="specify value of quaternion in grain 0")

(options, args) = parser.parse_args()

filename = args[0]

n_spheres = options.ngrains
crystal_sym = options.crystal_sym
if crystal_sym:
  print "use crystal symmetry..."

double_precision = options.double_precision
if double_precision:
  print "use double precision..."

nx = options.ncells
ny = options.ncells
nz = options.ncells

if ( options.nx ) : nx = options.nx
if ( options.ny ) : ny = options.ny
if ( options.nz ) : nz = options.nz

width = 0.
if ( options.width ) : width = options.width

if ( not ( nx and ny and nz ) ) :
  print "Error: either -n or all of -nx -ny -nz are required"
  sys.exit(1)


QLEN = options.qlen
print "qlen=",QLEN

radius = options.radius
if radius is None :
  radius = nx / 10
  print "Radius not specified: using nx/n_spheres = %d" % radius

phase_inside  = options.phase_in
phase_outside = options.phase_out

nomconc       = options.nomconc
conc_inside   = options.concentration_in
conc_outside  = options.concentration_out
if conc_inside is None :
  conc_inside = nomconc

temperature_inside   = options.temperature_in
temperature_outside  = options.temperature_out

quat_outside  = options.quat_out
if( not(quat_outside is None) ):
  qout = map( float, string.split( options.quat_out, ',' ) )
  print 'qout=',qout

periodic = map( int, string.split( options.periodic, ',') )

#-----------------------------------------------------------------------
nspecies=0
if ( not ( nomconc is None ) ):
  c = map( float, string.split( options.nomconc, ',' ) )
  nspecies=len(c)
  print "Nominal composition=",c
if ( not ( conc_inside is None ) ):
  ci = map( float, string.split( options.concentration_in, ',' ) )
  if nspecies==0:
    nspecies=len(ci)
  print "Composition inside=",ci
else:
  ci = N.zeros( nspecies, N.float32 )
if ( not ( conc_outside is None ) ):
  co = map( float, string.split( options.concentration_out, ',' ) )
  print "Composition outside=",co
else:
  co = N.zeros( nspecies, N.float32 )

print "nspecies=",nspecies

# generate quaternions corresponding to random orientations
random.seed( 11234 )
quat_inside=[]
nangles = 20.
h = 2.*math.pi/nangles

h1=0.1
h2=0.5*h1
h3=h2

#-----------------------------------------------------------------------
def setRandomQinSpheres():
  print 'setRandomQinSpheres...'
  if QLEN>0:
    # crystal symmetry
    if ( QLEN == 4 ) :
      qref = Q.makeNormalizedQuat() # (1,0,0,0)
    else:
      qref = Q.makeNormalizedQuat2()
    for g in range(n_spheres):

      #pich a discrete angle between -pi and pi
      t = math.pi * random.uniform(-1, 1.)
      n = math.floor( t/h )
      t = n*h

      if ( QLEN == 1 ) :

        q0 = t
        q1 = 0.
        q2 = 0.
        q3 = 0.
        if crystal_sym :
          if t<0.:
            t = t+math.pi # between 0 and pi
          if t>0.5*math.pi:
            t = t-0.5*math.pi # between 0 and pi/2

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

        if crystal_sym :
          # crystal symmetry
          if ( QLEN == 4 ) :
            q  = Q.makeQuat( q0,q1,q2,q3 )
          else:
            q  = Q.makeQuat2( q0,q1 )
          
          print 'random q=',q
          print 'qref=',qref
          q = Q.quatSymm( qref, q )
          print 'rotated q=',q
          q = Q.quatSymm( qref, q )
          print 'q=',q
          q0 = q[0]
          q1 = q[1]
          if ( QLEN == 4 ) :
            q2 = q[2]
            q3 = q[3]
          if n_spheres==2:
            Q.copyQuat(q,qref)

      quat_inside.append( Q.makeQuat( q0, q1, q2, q3 ) )
      print '--- q=',quat_inside[g]

#distance square function for periodic bc
if periodic[0]==1:
  rangex=range(-1,2)
else:
  rangex=range(1)

if periodic[1]==1:
  rangey=range(-1,2)
else:
  rangey=range(1)

if periodic[2]==1:
  rangez=range(-1,2)
else:
  rangez=range(1)

def distance2(x1,y1,z1,x2,y2,z2):
  d2=(x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2
  for ir in rangex:
    dx2=(x1-x2+ir*nx)**2
    if dx2<=d2:
      for jr in rangey:
        dy2=(y1-y2+jr*ny)**2
        if dy2<=d2:
          for kr in rangez:
            dz2=(z1-z2+kr*nz)**2
            l2=dx2+dy2+dz2
            if l2<d2:
              d2=l2
  return d2

#distance square function for periodic bc
def distance2_1d_x(x1,x2):
  d2=(x1 - x2)**2
  for ir in rangex:
    l2=(x1-x2+ir*nx)**2
    if l2<d2:
      d2=l2
  return d2

def distance2_1d_y(y1,y2):
  d2=(y1 - y2)**2
  for ir in rangey:
    l2=(y1-y2+ir*ny)**2
    if l2<d2:
      d2=l2
  return d2

def distance2_1d_z(z1,z2):
  d2=(z1 - z2)**2
  for ir in rangez:
    l2=(z1-z2+ir*nz)**2
    if l2<d2:
      d2=l2
  return d2

#-----------------------------------------------------------------------
def setBainQinSpheres():
  print 'setBainQinSpheres...'
  if QLEN>0:
    q0=math.cos( math.pi/8. )
    q1=math.sin( math.pi/8. )
    for g in range(n_spheres):
      q=[]
      q.append(q0)
      for i in range(3): 
        if( (g+i)%3 == 0 ):
          q.append(q1)
        else:
          q.append(0.)

      quat_inside.append( Q.makeNormalizedQuat( q[0], q[1], q[2], q[3] ) )
      print '--- q=',quat_inside[g]

#-----------------------------------------------------------------------
# Variables and functions for symmetry rotations

N_ROTATIONS = Q.getNumberSymmetryRotations( QLEN )
qrot = []
for n in range( 8 ) :
  nr = int( N_ROTATIONS * random.random() )
  qrot.append( Q.getQuatSymmetryRotation( nr, qlen=QLEN ) )

def rotateRandomly() :
  nr = int( N_ROTATIONS * random.random() )
  return Q.getQuatSymmetryRotation( nr, qlen=QLEN )

def rotateByOctant( x, y, z, cx, cy, cz ) :
  global qrot
  if ( x < cx ) :
    if ( y < cy ) :
      if ( z < cz ) :
        qr = qrot[0]
      else :
        qr = qrot[1]
    else :
      if ( z < cz ) :
        qr = qrot[2]
      else :
        qr = qrot[3]
  else :
    if ( y < cy ) :
      if ( z < cz ) :
        qr = qrot[4]
      else :
        qr = qrot[5]
    else :
      if ( z < cz ) :
        qr = qrot[6]
      else :
        qr = qrot[7]
  return qr

#-----------------------------------------------------------------------
# simple algorithm to move nuclei away from each other using a dynamics
# based on forces between pairs.
# Forces depend on distances and go to zero at distance drmin
def imposeMinDistance(cx,cy,cz,nx,ny,nz,drmin):

  ncenters=len(cx)       

  fac=[0.,-1.,1.]

  step=1
  npairs=1000
  while npairs>0:
    print 'Step ',step,
    npairs=0

    #reset forces
    forces=[]
    for index in range(ncenters):
      forces.append([0.,0.,0.])
    
    for index0 in range(ncenters):
      x0=cx[index0]
      y0=cy[index0]
      z0=cz[index0]
      for index1 in range(index0+1,ncenters):
        if(index1!=index0):
          x1=cx[index1]
          y1=cy[index1]
          z1=cz[index1]
      
          dr0=1000.
          dx0=nx
          dy0=ny
          dz0=nz
          #check minimage
          for i in range (3):
            dx=x1-x0+fac[i]*nx
            if dx<dr0:
              for j in range (3):
                dy=y1-y0+fac[j]*ny
                if dy<dr0:
                  for k in range (3):
                    dz=z1-z0+fac[k]*nz
                    dr=math.sqrt(dx*dx+dy*dy+dz*dz)
                    if dr<dr0:
                      dr0=dr
                      dx0=dx
                      dy0=dy
                      dz0=dz
          #print dr0
          if dr0<drmin:
            factor=0.5
            if dr0>0.000001:
              af=factor*(drmin-dr0+0.0001)
              invdr0=1./dr0
              fx=af*dx0*invdr0
              fy=af*dy0*invdr0
              fz=af*dz0*invdr0
            else:
              fx=factor*drmin*random.random()
              fy=factor*drmin*random.random()
              fz=factor*drmin*random.random()
            forces[index0][0]=forces[index0][0]-fx
            forces[index0][1]=forces[index0][1]-fy
            forces[index0][2]=forces[index0][2]-fz
            forces[index1][0]=forces[index1][0]+fx
            forces[index1][1]=forces[index1][1]+fy
            forces[index1][2]=forces[index1][2]+fz
              
            npairs=npairs+1
    print '--- number of overlaping pairs=',npairs
    step=step+1
   
    # move cx and cy
    maxf=0.
    for index in range(ncenters):
      cx[index]=cx[index]+forces[index][0];
      cy[index]=cy[index]+forces[index][1];
      cz[index]=cz[index]+forces[index][2];
      norm2f=forces[index][0]*forces[index][0]
      norm2f=norm2f+forces[index][1]*forces[index][1]
      norm2f=norm2f+forces[index][2]*forces[index][2]
      
      if( norm2f>maxf ):
        maxf=norm2f
    print 'max force=',math.sqrt(maxf)

#-----------------------------------------------------------------------
# Open and define file

#f = NetCDF.NetCDFFile( filename, 'w' )
f = nc4.Dataset(filename, 'w', format='NETCDF4') 

f.createDimension( 'x', nx )
f.createDimension( 'y', ny )
f.createDimension( 'z', nz )
f.createDimension( 'ns', nspecies )
if QLEN>0:
  f.createDimension( 'qlen', QLEN )

ncquat=[]
ncconc = []

if double_precision:
  print 'Data in double precision...'
  for s in range(nspecies):
    c_comp = f.createVariable( 'concentration%d' % s, 'd', ('z','y','x') )
    ncconc.append(c_comp)
  if not(temperature_outside is None) or not(temperature_inside is None):
    nctemp  = f.createVariable( 'temperature', 'd', ('z','y','x') )
    temperature  = N.ones( (nz,ny,nx), N.float64 )
  ncphase = f.createVariable( 'phase',         'd', ('z','y','x') )
  if( options.three ):
    nceta = f.createVariable( 'eta',         'd', ('z','y','x') )

  for m in range(QLEN):
    name = "quat%d" % (m+1)
    print name
    ncquat.append( f.createVariable( name, 'd', ('z','y','x') ) )

  conc  = N.ones( (nspecies,nz,ny,nx), N.float64 )
  phase = N.zeros( (nz,ny,nx), N.float64 )
  if QLEN>0:
    quat  = N.zeros( (QLEN,nz,ny,nx), N.float64 )
  if( options.three ):
    eta = N.zeros( (nz,ny,nx), N.float64 )
else:
  print 'Data in single precision...'
  for s in range(nspecies):
    c_comp = f.createVariable( 'concentration%d' % s , 'f', ('z','y','x') )
    ncconc.append(c_comp)
  if not(temperature_outside is None) or not(temperature_inside is None):
    nctemp  = f.createVariable( 'temperature', 'f', ('z','y','x') )
    temperature  = N.ones( (nz,ny,nx), N.float32 )
  ncphase = f.createVariable( 'phase',         'f', ('z','y','x') )
  if( options.three ):
    nceta = f.createVariable( 'eta',         'f', ('z','y','x') )

  for m in range(QLEN):
    name = "quat%d" % (m+1)
    print name
    ncquat.append( f.createVariable( name, 'f', ('z','y','x') ) )

  conc  = N.ones( (nspecies,nz,ny,nx), N.float32 )
  phase = N.zeros( (nz,ny,nx), N.float32 )
  if QLEN>0:
    quat  = N.zeros( (QLEN,nz,ny,nx), N.float32 )
  if( options.three ):
    eta = N.zeros( (nz,ny,nx), N.float32 )

#-----------------------------------------------------------------------
def generateRandomCenters(cx,cy,cz,n):
  for g in range(n):
    cx.append( random.randint(1,nx) )
  for g in range(n):
    cy.append( random.randint(1,ny) )
  for g in range(n):
    cz.append( random.randint(1,nz) )
  else:
    for g in range(n):
      cz.append( 1 )


#-----------------------------------------------------------------------
def generateJitterCenters2D(cx,cy,cz,ncx,ncy,jitter):
  print 'Generate jitter centers in 2D...'
  hx=float(nx)/float(ncx)
  hy=float(ny)/float(ncy)
  for gx in range(ncx):
    for gy in range(ncy):
      jx=(2.*random.random()-1.0)*0.5*jitter
      jy=(2.*random.random()-1.0)*0.5*jitter
      cx.append( math.floor((gx+0.5)*hx+jx) )
      cy.append( math.floor((gy+0.5)*hy+jy) )
      cz.append( 1 )
        
#-----------------------------------------------------------------------
def generateJitterCenters3D(cx,cy,cz,ncx,ncy,ncz,jitter):
  print 'Generate jitter centers in 3D...'
  hx=float(nx)/float(ncx)
  hy=float(ny)/float(ncy)
  hz=float(nz)/float(ncz)
  for gx in range(ncx):
    for gy in range(ncy):
      for gz in range(ncz):
        jx=(2.*random.random()-1.0)*0.5*jitter
        jy=(2.*random.random()-1.0)*0.5*jitter
        jz=(2.*random.random()-1.0)*0.5*jitter
        cx.append( math.floor((gx+0.5)*hx+jx) )
        cy.append( math.floor((gy+0.5)*hy+jy) )
        cz.append( math.floor((gz+0.5)*hz+jz) )
        
#-----------------------------------------------------------------------

mind=4.*radius
halfmind2=0.25*mind*mind

# generate random (x,y) positions
random.seed(21361)
cx=[]
cy=[]
cz=[]
r=[]

if n_spheres==1:
  r.append(radius)
  if ( not (options.center0 is None) ):
    center = map( float, string.split( options.center0, ',' ) )
    cx.append( center[0] )
    cy.append( center[1] )
    if len(center)>2:
      cz.append( center[2] )
    else:
      cz.append(1)
  else:
    cx.append( nx/2 )
    cy.append( ny/2 )
    if nz>1:
      cz.append( nz/2 )
    else:
      cz.append(1)
else:
  if n_spheres==2:
    print 'nspheres=',n_spheres
    cx.append( nx/3 )
    cx.append( 2*nx/3 )
    cy.append( ny/3 )
    cy.append( 2*ny/3 )
    cz.append( nz/2 )
    cz.append( nz/2 )
    r.append(radius)
    r.append(radius* 1.5)
    
  else:
    if nz==1:
      rtn=int(math.floor(math.sqrt(n_spheres)))
      nnn=rtn*rtn
    if nz>1:
      rtn=int(round(math.pow(n_spheres,1./3.)))
      nnn=rtn*rtn*rtn
    if( nnn==n_spheres ):
      hx=float(nx)/float(rtn)
      print 'hx=',hx
      #jitter= hx-mind
      jitter= options.jitter_factor*hx
      print 'jitter=',jitter
      if nz==1:
        generateJitterCenters2D(cx,cy,cz,rtn,rtn,jitter)
      if nz>1:
        generateJitterCenters3D(cx,cy,cz,rtn,rtn,rtn,jitter)
    else:
      generateRandomCenters(cx,cy,cz,n_spheres)
    for k in range(n_spheres):
      r.append(radius)

    # move positions to avoid any overlaps between grains
    #imposeMinDistance(cx,cy,cz,nx,ny,nz,2.4*radius)
    imposeMinDistance(cx,cy,cz,nx,ny,nz,mind)

#-----------------------------------------------------------------------
# Fill data arrays
print 'Fill Phase values'
if ( options.phase_out ) :
  for k in range( nz ) :
    for j in range( ny ) :
      for i in range( nx ) :
        phase[k,j,i] = phase_outside

if( options.three ):
  for k in range( nz ) :
    for j in range( ny ) :
      for i in range( nx ) :
        phase[k,j,i] = 1.
        eta[k,j,i] = 0.



vol=nx*ny*nz
vs=0.
vl=0.;

#fill phase value and compute volume solid
for g in range(n_spheres):
  r_sq = r[g]**2
  threshold = (r[g]+width)**2
  print 'sphere ',g,', center: ',cx[g],cy[g],cz[g],', radius: ',r[g]
  for k in range( nz ) :
    z = k + 0.5
    dz2=distance2_1d_z(z,cz[g])
    if dz2<threshold :
      for j in range( ny ) :
        y = j + 0.5
        dy2=distance2_1d_y(y,cy[g])
        if dy2<threshold :
          for i in range( nx ) :
            x = i + 0.5
       
            distance_sq = distance2(x,y,z,cx[g],cy[g],cz[g])
            d = distance_sq - r_sq
            if( width>0. ):
              sq=N.sqrt(abs(d))
              if( sq<2.*width or d<0. ):
                phase[k,j,i] = phase_inside*0.5*(1.+N.tanh(0.5*sq/width))
                vs=vs+phase[k,j,i]
                if( options.three ):
                  eta[k,j,i] = g%2
            else:
              if( d<0. ):
                phase[k,j,i] = phase_inside
                vs=vs+phase[k,j,i]
                if( options.three ):
                  eta[k,j,i] = g%2

vl=vol-vs

#fill quat values
if ( not (options.quat0 is None) ):
  q = map( float, string.split( options.quat0, ',' ) )
  quat_inside.append( q )
else:
  if options.bain :
    setBainQinSpheres()
  else:
    setRandomQinSpheres()

if QLEN>0:
  gmin=0
  print 'Fill quaternion values...'
  for i in range( nx ) :
    x = i + 0.5
    print 'Plane x=',x
    for j in range( ny ) :
      y = j + 0.5
      for k in range( nz ) :
        z = k + 0.5
        d2min=distance2(x,y,z,cx[gmin],cy[gmin],cz[gmin])
        if( not(quat_outside is None)):
          if( d2min<radius*radius ):
            qi=quat_inside[gmin]
          else:
            qi=[]
            for m in range(QLEN):
              qi.append(qout[m])
        else:
          if d2min>halfmind2:
            for g in range(n_spheres):
              distance_sq = distance2(x,y,z,cx[g],cy[g],cz[g])
              if distance_sq<d2min:
                d2min=distance_sq
                gmin=g
                if d2min<halfmind2:
                  break
          qi=quat_inside[gmin]

        qr = None
        if ( options.symmetry_test ) :
          qr = rotateByOctant( x, y, z, nx/2, ny/2, nz/2 )
        if ( not qr is None ) :
          if ( QLEN == 1 ) :
            qi += qr
          else :
            qi = Q.rotateQuat( qi, qr )

        for m in range(QLEN):
          quat[m,k,j,i] = qi[m]

#fill conc values
if nspecies>0:
  if ( not ( conc_inside is None ) ):
    if ( not ( nomconc is None ) and vl>0 ):
      for s in range(nspecies):
        co[s] = (c[s]*vol-ci[s]*vs)/vl
      print "Calculated composition outside=",co
  if ( not ( options.concentration_out is None ) ):
    if ( not ( nomconc is None ) and vs>0 ):
      conc_inside = (c[0]*vol-conc_outside*vl)/vs
      print "Calculated composition inside=",conc_inside
  if( ( conc_outside is None ) and ( conc_inside is None ) ):
    conc_inside = nomconc
    conc_outside = nomconc

  if ( not ( conc_outside is None ) and not ( conc_inside is None ) ):
    for s in range(nspecies):
      print "Calculated nominal Composition=",(vl*co[s]+vs*ci[s])/vol

  print 'set composition...'
  for k in range( nz ) :
    for j in range( ny ) :
      for i in range( nx ) :
        for s in range(nspecies):
          conc[s,k,j,i]  = ci[s]*phase[k,j,i]+co[s]*(1.-phase[k,j,i])

if not(temperature_inside is None):
  print 'Fill temperature values'
  print 'temperature_inside =',temperature_inside
  print 'temperature_outside=',temperature_outside
  for k in range( nz ) :
    for j in range( ny ) :
      for i in range( nx ) :
        temperature[k,j,i]  = temperature_outside

  for g in range(n_spheres):
    r_sq = r[g]**2
    print 'temperature, grain ',g
    for k in range( nz ) :
      z = k + 0.5
      dz2=distance2_1d_z(z,cz[g])
      if dz2<r_sq :
        for j in range( ny ) :
          y = j + 0.5
          dy2=distance2_1d_y(y,cy[g])
          if dy2<r_sq :
            for i in range( nx ) :
              x = i + 0.5
         
              distance_sq = distance2(x,y,z,cx[g],cy[g],cz[g])
              d = distance_sq - r_sq
              if( d<0. ):
                temperature[k,j,i]  = temperature_inside

#-----------------------------------------------------------------------
# Write data to file and close

print 'Write data to file'
if not(temperature_inside is None):
  nctemp[:,:,:]= temperature
ncphase[:,:,:]=phase
if( options.three ):
  nceta[:,:,:]=eta
for m in range(QLEN):
  ncquat[m][:,:,:]=quat[m,:,:,:]

if ( nspecies>0 ):
  for s in range(nspecies):
    ncconc[s][:,:,:]=conc[s,:,:,:]

f.close()


