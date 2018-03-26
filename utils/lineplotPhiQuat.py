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
#  visit -cli -verbose -nowin -s lineplotPhiQuat.py
#in main visit data directory
import sys, string, os

visitdir =os.getcwd()
filename=visitdir+"/dumps.visit"

myfile=open(filename,'r')
L=myfile.readlines()

frames=[]
for line in L:
  words=string.split(line)
  frames.append(words[0])

print "#",frames

#load data
DB=[]
for frame in frames:
  DB.append(visitdir+"/"+frame)

Vars=[]
Vars.append("phase")
Vars.append("concentration0")
Vars.append("q0")
Vars.append("q1")
Vars.append("q2")
Vars.append("q3")
Vars.append("temperature")

#plot data
for db in DB:
  OpenDatabase(db)

  vars2remove=[]
  for var in Vars:
    ierr = AddPlot("Pseudocolor", var)
    if ierr==1:
      DrawPlots()

      Query("SpatialExtents")
      pxy = GetQueryOutputValue()
      p0 = ( pxy[0], 0. )
      p1 = ( pxy[1], 0. )

      # Do a lineout on variable to produce curve.
      Lineout(p0, p1, ("default"))
    else:
      print '#remove var ',var,' (not found)...'
      vars2remove.append(var)
  for var in vars2remove:
    Vars.remove(var)

# extract data into (x,y) format
SetActiveWindow(2)
n=len(Vars)
nplots=n*len(DB)
for k in range(len(DB)):
  for j in range(n):
    SetActivePlots(k*n+j)
    vals = GetPlotInformation()["Curve"]

    Query("Time")
    t0 = GetQueryOutputValue()

    filename = "%s_%s.dat" % (Vars[j],t0)
    output = open(filename,'w')
    output.write("# %s --- time: %g" % (Vars[j],t0))
    print "# %s --- time: %g" % (Vars[j],t0)

    # Write data as "x  1.-y"
    for i in range(len(vals) / 2):
      #print "%g  %g" % (p0[0]+vals[2*i], vals[2*i+1])
      output.write("%g  %g\n" % (p0[0]+vals[2*i], vals[2*i+1]))
    #print "\n"
    output.write("\n")

sys.exit()
