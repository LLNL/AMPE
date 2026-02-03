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
#  visit -cli -verbose -nowin -o dumps.visit -s plotComposition0.py
#note on Windows:
#  lauch Xming
import os
import sys

DeleteAllPlots()

if( len(sys.argv)>1 ):
  db = sys.argv[1]
else:
  db = "dumps.visit"
OpenDatabase(db)

HideActivePlots()

#DefineScalarExpression( "concentration","1.-concentration0" )
AddPlot( "Pseudocolor", "concentration0" )

annot_atts = AnnotationAttributes()
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
#remove labels
xl = AxisLabels()
xl.SetVisible( 0 )
xa.SetLabel( xl )
ya.SetLabel( xl )

SetAnnotationAttributes( annot_atts )

DrawPlots()

Query("SpatialExtents")
pxy = GetQueryOutputValue()
ll = ( pxy[0], pxy[2] )
ur = ( pxy[1], pxy[3])

p = PseudocolorAttributes()
N = GetDatabaseNStates()
maxval=0.
minval=1.
for i in range(0,N):
  SetTimeSliderState( i )
  Query("MinMax")
  mmval=GetQueryOutputValue()
  minval=min(minval,mmval[0])
  maxval=max(maxval,mmval[1])

# Set the min/max values
print("MinMax={},{}".format(minval,maxval))
p.min, p.minFlag = minval, 1
p.max, p.maxFlag = maxval, 1

SetPlotOptions(p)

v0 = View2DAttributes()
v0.SetFullFrameActivationMode( 0 )
v0.SetWindowCoords( (ll[0],ur[0],ll[1],ur[1]) )
SetView2D( v0 )

DrawPlots()

swa = SaveWindowAttributes()
swa.family = swa.PNG
swa.family = 0
swa.width = 1280
swa.height = 1280
SetSaveWindowAttributes( swa )

#plot frames
for i in range(0,N):
  SetTimeSliderState( i )
  swa.fileName = "composition_%04d" % i
  SetSaveWindowAttributes( swa )
  SaveWindow()
  DrawPlots()

sys.exit()
