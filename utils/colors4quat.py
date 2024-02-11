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
#  visit -cli -verbose -nowin -s colors4quat.py 
#in main visit data directory
import os
import sys

visitdir =os.getcwd()
filename=visitdir+"/dumps.visit"

OpenDatabase( filename )

DeleteAllPlots()

DefineVectorExpression( "quat",
"color4(128+127*q1,128+127*q2,128+127*q3,191+64*q0)" )

AddPlot( "Truecolor", "quat" )

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
SetAnnotationAttributes( annot_atts )

DrawPlots()

Query("SpatialExtents")
pxy = GetQueryOutputValue()
ll = ( pxy[0], pxy[2] )
ur = ( pxy[1], pxy[3])

v0 = View2DAttributes()
v0.SetFullFrameActivationMode( 0 )
v0.SetWindowCoords( (ll[0],ur[0],ll[1],ur[1]) )
SetView2D( v0 )

swa = SaveWindowAttributes()
swa.family = swa.PNG
swa.family = 0
swa.width = 1280
swa.height = 1280
SetSaveWindowAttributes( swa )

N = GetDatabaseNStates()

#for i in range( 1 ) :
for i in range( N ) :
   SetTimeSliderState( i )
   swa.fileName = "q_%04d" % i
   SetSaveWindowAttributes( swa )
   SaveWindow()
   DrawPlots()
