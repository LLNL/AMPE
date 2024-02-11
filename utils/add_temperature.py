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
import sys
import string
from optparse import OptionParser

# other required packages
import Numeric as N
from Scientific.IO import NetCDF

#-----------------------------------------------------------------------
# command-line arguments and options

usage = "Usage: %prog [options] filename"

parser = OptionParser( usage = usage )

(options, args) = parser.parse_args()

if ( len( args ) < 1 ) :
  parser.error( "filename required argument missing" )
filename = args[0]

#-----------------------------------------------------------------------
# Open and define file

f = NetCDF.NetCDFFile( filename, 'a' )

print f.dimensions.keys()

nx = f.dimensions['x']
ny = f.dimensions['y']
nz = f.dimensions['z']

print nx, ny, nz

nctemp = f.createVariable( 'temperature', N.Float, ('z','y','x') )

temp = N.ones( (nz,ny,nx), N.Float )

#-----------------------------------------------------------------------
# Fill data arrays

tavg = 837.15
delt = 40.0
dt = delt / nx
t0 = tavg - 0.5 * delt

for i in range( nx ) :
  x = i + 0.5
  temp[:,:,i] = t0 + x * dt

#-----------------------------------------------------------------------
# Write data to file and close

nctemp.assignValue( temp )

f.close()
