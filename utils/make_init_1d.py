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
import sys
from optparse import OptionParser

# my required packages
import quat as Q

# other required packages
import Numeric as N
from Scientific.IO import NetCDF

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

(options, args) = parser.parse_args()

if ( len( args ) < 1 ) :
  parser.error( "filename required argument missing" )
filename = args[0]

QLEN = 1

nx = None
ny = 4
nz = 1

if ( not options.nx is None ) : nx = options.nx
if ( not options.ny is None ) : ny = options.ny
if ( not options.nz is None ) : nz = options.nz

if ( not nx ) :
  print "Error: --nx=[NCELLS in x direction] is required"
  sys.exit(1)

#-----------------------------------------------------------------------
# Open and define file

f = NetCDF.NetCDFFile( filename, 'w' )

f.createDimension( 'x', nx )
f.createDimension( 'y', ny )
f.createDimension( 'z', nz )
f.createDimension( 'qlen', QLEN )

ncphase = f.createVariable( 'phase', N.Float, ('z','y','x') )

ncquat = []
for n in range( QLEN ) :
  q_comp = f.createVariable( 'quat%d' % (n+1), N.Float, ('z','y','x') )
  ncquat.append( q_comp )

ncconc = f.createVariable( 'concentration', N.Float, ('z','y','x') )

nctemp = f.createVariable( 'temperature', N.Float, ('z','y','x') )

phase = N.ones( (nz,ny,nx), N.Float )
quat = N.zeros( (QLEN,nz,ny,nx), N.Float )
conc = N.ones( (nz,ny,nx), N.Float )
temp = N.ones( (nz,ny,nx), N.Float )

cx = nx / 2

#-----------------------------------------------------------------------
# Fill data arrays

phase[:,:,:] = 1.0
quat[:,:,:,:] = 0.0

conc_center = 0.02
del_conc = 0.04

temp_center = 860.0
del_temp = -60.0

for i in range( nx ) :
  x = i + 0.5
  distance = abs(x - cx)

  conc[:,:,i] = conc_center + del_conc * distance / cx
  temp[:,:,i] = temp_center + del_temp * distance / cx

#-----------------------------------------------------------------------
# Write data to file and close

ncphase.assignValue( phase )

for n in range( QLEN ) :
  ncquat[n].assignValue( quat[n,...] )

ncconc.assignValue( conc )

nctemp.assignValue( temp )

f.close()
