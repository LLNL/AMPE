# Copyright (c) 2018, Lawrence Livermore National Security, LLC and
# UT-Battelle, LLC.
# Produced at the Lawrence Livermore National Laboratory and
# the Oak Ridge National Laboratory
# LLNL-CODE-747500
# All rights reserved.
# This file is part of AMPE. 
# For details, see https://github.com/LLNL/AMPE
# Please also read AMPE/LICENSE.
# 
# standard packages
import math
import random
import sys
import string
import csv
from optparse import OptionParser

# other required packages
import numpy as N
import netCDF4 as nc4

print (sys.path)

#-----------------------------------------------------------------------
# command-line arguments and options

usage = "Usage: %prog [options] filename"

parser = OptionParser( usage = usage )

parser.add_option( "--spheres", type="string", help="csv file with list of spheres" )
parser.add_option( "-n", "--ncells", type="int",
                   help="number of cells in x,y,z (will all be equal)" )
parser.add_option( "-x", "--nx", type="int",
                   help="number of cells in x direction" )
parser.add_option( "-y", "--ny", type="int",
                   help="number of cells in y direction" )
parser.add_option( "-z", "--nz", type="int", default=1,
                   help="number of cells in z direction" )
parser.add_option( "-c", "--nomconc", type="string", # something like "0.1,0.2",
                   help="nominal concentration" )
parser.add_option( "--concentration-in", type="string",
                   help="concentration in interior region" )
parser.add_option( "--concentration-out", type="string",
                   help="concentration in exterior region" )
parser.add_option( "-w", "--width", type="float",
                   help="interface width", default=0.0 )
parser.add_option("--periodic", type="string", # something like "1,0,1" for periodic in x and z
                  default="1,1,1",
                  help="specify which dimensions are periodic") 

(options, args) = parser.parse_args()

filename = args[0]

#read spheres coordinates and radius
spheres_filename =  options.spheres
nspheres = 0

nx = options.ncells
ny = options.ncells
nz = options.ncells

radius = []
centers = []
with open(spheres_filename, mode ='r') as file:    
  csvFile = csv.reader(file)
  for line in csvFile:
    print(line)
    nspheres = nspheres+1
    radius.append(eval(line[3]))
    cx = eval(line[0])
    cy = eval(line[1])
    cz = eval(line[2])
    if nz==1:
      cz = 1    
    centers.append([cx,cy,cz])

print("Centers:")
print(centers)
print("Radius:")
print(radius)

if ( options.nx ) : nx = options.nx
if ( options.ny ) : ny = options.ny
if ( options.nz ) : nz = options.nz

width = 0.
if ( options.width ) : width = options.width

if ( not ( nx and ny and nz ) ) :
  print ("Error: either -n or all of -nx -ny -nz are required")
  sys.exit(1)

nomconc       = options.nomconc
conc_inside   = options.concentration_in
conc_outside  = options.concentration_out
if conc_inside is None :
  conc_inside = nomconc

periodic = options.periodic.split( ',')

#-----------------------------------------------------------------------
nspecies=0
if ( not ( nomconc is None ) ):
  tmp = options.nomconc.split( ',' )
  c = [float(i) for i in tmp]
  nspecies=len(set(c))
  print ("Nominal composition={}".format(c))
if ( not ( conc_inside is None ) ):
  tmp = options.concentration_in.split( ',' )
  ci = [float(i) for i in tmp]
  if nspecies==0:
    nspecies=len(ci)
  print ("Composition inside={}".format(ci))
else:
  ci = N.zeros( nspecies, N.float32 )
if ( not ( conc_outside is None ) ):
  tmp = options.concentration_out.split( ',' )
  co = [float(i) for i in tmp]
  print ("Composition outside={}".format(co))
else:
  co = N.zeros( nspecies, N.float32 )

print ("nspecies={}".format(nspecies))


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
# Open and define file

f = nc4.Dataset(filename, 'w', format='NETCDF4') 

f.createDimension( 'x', nx )
f.createDimension( 'y', ny )
f.createDimension( 'z', nz )
f.createDimension( 'ns', nspecies )

ncconc = []
ncphase = []
for s in range(nspecies):
  c_comp = f.createVariable( 'concentration%d' % s , 'f', ('z','y','x') )
  ncconc.append(c_comp)
for p in range(nspheres):
  phase = f.createVariable( 'phase%d' % p,         'f', ('z','y','x') )
  ncphase.append(phase)

conc  = N.zeros( (nspecies,nz,ny,nx), N.float32 )
phase = N.zeros( (nspheres,nz,ny,nx), N.float32 )

#-----------------------------------------------------------------------

vol=nx*ny*nz
vs=0.
vl=0.;

#fill phase value and compute volume solid
for g in range(nspheres):
  center = centers[g]
  r_sq = radius[g]**2
  threshold = (radius[g]+5.*width)**2
  print ("sphere {}, center: {},{},{}, radius: {}".format(g,center[0],center[1],center[2],radius[g]))
  for k in range( nz ) :
    z = k + 0.5
    dz2=distance2_1d_z(z,center[2])
    print("z = {}".format(z))
    if dz2<threshold or nz==1:
      for j in range( ny ) :
        y = j + 0.5
        dy2=distance2_1d_y(y,center[1])
        if dy2<threshold :
          for i in range( nx ) :
            x = i + 0.5

            #compute distance to center
            distance_sq = distance2(x,y,z,center[0],center[1],center[2])
            #compare with rads[i]ius
            d = N.sqrt(distance_sq) - N.sqrt(r_sq)
            if( d<0. ):
              phase[g,k,j,i] = 1.
            if( width>0. ):
              if( abs(d)<8.*width ):
                phase[g,k,j,i] = 0.5*(1.+N.tanh(-1.*d/(2.*width)))

#fill conc values
if nspecies>0:
  if ( not ( conc_inside is None ) ):
    if ( not ( nomconc is None ) and vl>0 ):
      for s in range(nspecies):
        co[s] = (c[s]*vol-ci[s]*vs)/vl
      print ("Calculated composition outside={}".format(co))
  if ( not ( options.concentration_out is None ) ):
    if ( not ( nomconc is None ) and vs>0 ):
      conc_inside = (c[0]*vol-conc_outside*vl)/vs
      print ("Calculated composition inside={}".format(conc_inside))
  if( ( conc_outside is None ) and ( conc_inside is None ) ):
    conc_inside = nomconc
    conc_outside = nomconc

  if ( not ( conc_outside is None ) and not ( conc_inside is None ) ):
    for s in range(nspecies):
      print ("Calculated nominal Composition={}".format((vl*co[s]+vs*ci[s])/vol))

  print ("set composition for {} species".format(nspecies))
  for k in range( nz ) :
    for j in range( ny ) :
      for i in range( nx ) :
        for s in range(nspecies):
          phi=0.
          for p in range(nspheres):
            phi=phi+phase[p,k,j,i]
          conc[s,k,j,i] = ci[s]*phi+co[s]*(1.-phi)

#-----------------------------------------------------------------------
# Write data to file and close

print ("Write data to file")
if ( nspecies>0 ):
  for s in range(nspecies):
    ncconc[s][:,:,:]=conc[s,:,:,:]
for p in range(nspheres):
  ncphase[p][:,:,:]=phase[p,:,:,:]

f.close()

