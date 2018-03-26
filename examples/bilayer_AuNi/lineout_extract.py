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
#usage:
#  visit -nowin -cli -verbose -s lineout_extract.py > lineout_extract.dat
import sys, string, os

maindir =os.getcwd()

# Adapt these parameters
visitdir="/v.bilayer_AuNi_restart/"
filename=maindir+"/bilayer_AuNi_restart.log"

#end points coordinates in [um]
p0 = ( 0.0,0.)
p1 = ( 3.2,0.)

myfile=open(filename,'r')
L=myfile.readlines()

#use dt as set for visit output frequency
dt=4.
t=dt
frames=["00000"]
search_string = 'cycle'
for line in L:
  num_matches = string.count(line, search_string)
  if num_matches:
    w=string.split(line)
    tc=eval(w[6])
    if tc>t:
      t=t+dt
      cycle=w[2]
      while len(cycle)<5:
        cycle="0"+cycle
      frames.append(cycle)

#print frames

#load data
datadir=maindir+visitdir
DB=[]
for frame in frames:
  DB.append(datadir+"visit_dump."+frame+"/summary.samrai")

Var = "concentration0"

#plot data
for db in DB:
  OpenDatabase(db)

  AddPlot("Pseudocolor", Var)
  DrawPlots()

  # Do a lineout on variable to produce curve.
  Lineout(p0, p1, ("default"))

# extract data into (x,y) format
SetActiveWindow(2)
nplots=len(DB)
for j in range(nplots):
  SetActivePlots(j)
  vals = GetPlotInformation()["Curve"]

  Query("Time")
  t0 = GetQueryOutputValue()

  print "#time: %g" % t0

  # Write data as "x  1.-y"
  for i in range(len(vals) / 2):
    print "%g  %g" % (p0[0]+vals[2*i], 1.-vals[2*i+1])
  print "\n"

sys.exit()
