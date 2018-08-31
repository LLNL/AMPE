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
#usage:
#  python plotPhase.py filename
import os
import sys
_visit_path = "/usr/gapps/visit/current/linux-x86_64/lib"
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
visit.AddPlot( "Pseudocolor", "phase" )

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
visit.SetAnnotationAttributes( annot_atts )

visit.DrawPlots()

visit.Query("SpatialExtents")
pxy = visit.GetQueryOutputValue()
ll = ( pxy[0], pxy[2] )
ur = ( pxy[1], pxy[3])

visit.Query("MinMax")
mmval=visit.GetQueryOutputValue()
print 'MinMax=',mmval

p = visit.PseudocolorAttributes()
# Set the min/max values
p.min, p.minFlag = 0., 1
p.max, p.maxFlag = 1., 1
visit.SetPlotOptions(p)

v0 = visit.View2DAttributes()
v0.SetFullFrameActivationMode( 0 )
v0.SetWindowCoords( (ll[0],ur[0],ll[1],ur[1]) )
visit.SetView2D( v0 )

visit.DrawPlots()

swa = visit.SaveWindowAttributes()
swa.family = swa.PNG
#swa.family = swa.JPEG
swa.family = 0
swa.width = 1280
swa.height = 1280
visit.SetSaveWindowAttributes( swa )

#plot all frame
N = visit.GetDatabaseNStates()
for i in range(N):
  visit.SetTimeSliderState( i )
  swa.fileName = "phase_%04d" % i
  visit.SetSaveWindowAttributes( swa )
  visit.SaveWindow()
  visit.DrawPlots()

sys.exit()

