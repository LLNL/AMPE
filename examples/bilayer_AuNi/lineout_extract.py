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
