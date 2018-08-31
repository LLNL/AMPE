# Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
# Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
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
# LLC, UT BATTELLE, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
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
