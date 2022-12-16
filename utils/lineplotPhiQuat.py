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
#usage:
#  visit -cli -verbose -nowin -s lineplotPhiQuat.py
#in main visit data directory
import sys, string, os

visitdir =os.getcwd()
filename=visitdir+"/dumps.visit"

myfile=open(filename,'r')
lines=myfile.readlines()

frames=[]
for line in lines:
  words=line.split()
  frames.append(words[0])

print("# {}".format(frames))

#load data
DB=[]
for frame in frames:
  DB.append(visitdir+"/"+frame)

Vars=[]
Vars.append("phase0")
Vars.append("concentration0")
Vars.append("concentration1")
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
      print("remove var {} (not found)...".format(var))
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
    print("# {} --- time: {}".format(Vars[j],t0))

    # Write data as "x  y"
    for i in range(len(vals) / 2):
      output.write("%g  %g\n" % (p0[0]+vals[2*i], vals[2*i+1]))
    output.write("\n")

sys.exit()
