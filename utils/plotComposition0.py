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
#  module load visit/2.13.3
#  export PYTHONPATH=$VISIT_DIR/2.13.3/linux-x86_64/lib/site-packages:$PYTHONPATH
#  python plotComposition0.py dumps.visit 
#note on Windows:
#  lauch Xming
import os
import sys

_visit_path = "/software/dev_tools/swtree/cs400_centos7.2_pe2016-08/visit/2.13.3/centos7.5_binary/visit2_13_3.linux-x86_64/2.13.3/linux-x86_64/lib"
sys.path.insert(0, _visit_path)
import visit

visit.AddArgument("-nowin")
visit.Launch()
visit.DeleteAllPlots()

if( len(sys.argv)>1 ):
  db = sys.argv[1]
else:
  db = "dumps.visit"
visit.OpenDatabase(db)

visit.HideActivePlots()

visit.DefineScalarExpression( "concentration",
"1.-concentration0" )
visit.AddPlot( "Pseudocolor", "concentration" )

annot_atts = visit.AnnotationAttributes()
#default value=1
annot_atts.SetLegendInfoFlag( 0 )
annot_atts.SetDatabaseInfoFlag( 0 )
annot_atts.SetUserInfoFlag( 0 )
axes = annot_atts.GetAxes2D()
xa = axes.GetXAxis()
ya = axes.GetYAxis()
xt = visit.AxisTitles()
yt = visit.AxisTitles()
xt.SetVisible( 0 )
yt.SetVisible( 0 )
xa.SetTitle( xt )
ya.SetTitle( yt )
#remove labels
xl = visit.AxisLabels()
xl.SetVisible( 0 )
xa.SetLabel( xl )
ya.SetLabel( xl )

visit.SetAnnotationAttributes( annot_atts )

visit.DrawPlots()

visit.Query("SpatialExtents")
pxy = visit.GetQueryOutputValue()
ll = ( pxy[0], pxy[2] )
ur = ( pxy[1], pxy[3])

p = visit.PseudocolorAttributes()
N = visit.GetDatabaseNStates()
maxval=0.
minval=1.
for i in range(0,N):
  visit.SetTimeSliderState( i )
  visit.Query("MinMax")
  mmval=visit.GetQueryOutputValue()
  minval=min(minval,mmval[0])
  maxval=max(maxval,mmval[1])

minval=0.5
maxval=0.9

# Set the min/max values
print("MinMax={},{}".format(minval,maxval))
p.min, p.minFlag = minval, 1
p.max, p.maxFlag = maxval, 1

visit.SetPlotOptions(p)

v0 = visit.View2DAttributes()
v0.SetFullFrameActivationMode( 0 )
v0.SetWindowCoords( (ll[0],ur[0],ll[1],ur[1]) )
visit.SetView2D( v0 )

visit.DrawPlots()

swa = visit.SaveWindowAttributes()
swa.family = swa.PNG
swa.family = 0
swa.width = 1280
swa.height = 1280
visit.SetSaveWindowAttributes( swa )

#plot frames
for i in range(0,N):
  visit.SetTimeSliderState( i )
  swa.fileName = "composition_%04d" % i
  visit.SetSaveWindowAttributes( swa )
  visit.SaveWindow()
  visit.DrawPlots()

sys.exit()
