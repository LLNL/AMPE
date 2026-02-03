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
#  visit -cli -verbose -nowin -o dumps.visit -s plotPhase.py
import os
import sys

DeleteAllPlots()

if( len(sys.argv)>1 ):
  db = sys.argv[1]
else:
  db = "dumps.visit"
OpenDatabase(db)

HideActivePlots()
flag = 1
i = 0
while(flag):
  fieldname = "phase"+str(i)
  print("Try to add field {}".format(fieldname))
  ierr = AddPlot( "Pseudocolor", fieldname )
  if ierr==0:
    flag = 0
  else:
    i=i+1

annot_atts = AnnotationAttributes()
#default value=1
annot_atts.SetLegendInfoFlag( 0 )
annot_atts.SetDatabaseInfoFlag( 0 )
annot_atts.SetUserInfoFlag( 0 )
axes = annot_atts.GetAxes2D()
xa = axes.GetXAxis()
ya = axes.GetYAxis()
xt = AxisTitles()
yt = AxisTitles()
xt.SetVisible( 0 )
yt.SetVisible( 0 )
xa.SetTitle( xt )
ya.SetTitle( yt )
SetAnnotationAttributes( annot_atts )

DrawPlots()

Query("SpatialExtents")
pxy = GetQueryOutputValue()
ll = ( pxy[0], pxy[2] )
ur = ( pxy[1], pxy[3])

Query("MinMax")
mmval=GetQueryOutputValue()
print ("MinMax={}".format(mmval))

p = PseudocolorAttributes()
# Set the min/max values
p.min, p.minFlag = 0., 1
p.max, p.maxFlag = 1., 1
SetPlotOptions(p)

v0 = View2DAttributes()
v0.SetFullFrameActivationMode( 0 )
v0.SetWindowCoords( (ll[0],ur[0],ll[1],ur[1]) )
SetView2D( v0 )

DrawPlots()

swa = SaveWindowAttributes()
swa.family = swa.PNG
#swa.family = swa.JPEG
swa.family = 0
swa.width = 1280
swa.height = 1280
SetSaveWindowAttributes( swa )

#plot all frame
N = GetDatabaseNStates()
for i in range(N):
  SetTimeSliderState( i )
  swa.fileName = "phase_%04d" % i
  SetSaveWindowAttributes( swa )
  SaveWindow()
  DrawPlots()

sys.exit()
