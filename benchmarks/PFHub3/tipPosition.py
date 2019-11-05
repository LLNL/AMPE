#usage:
#  visit -cli -verbose -nowin -s tipPosition.py
#in main visit data directory
import sys, string, os

visitdir =os.getcwd()
filename=visitdir+"/dumps.visit"

myfile=open(filename,'r')
lines=myfile.readlines()

frames=[]
counter=0
for line in lines:
  words=line.split()
  if counter%5 == 0:
    frames.append(words[0])
  counter=counter+1

print( "#{}".format(frames))

#load data
DB=[]
for frame in frames:
  DB.append(visitdir+"/"+frame)

Var="phase"

#extract coordinates of end points in x direction
OpenDatabase(DB[0])
ierr = AddPlot("Pseudocolor", Var)
if ierr==1:
  DrawPlots()
  Query("SpatialExtents")
  pxy = GetQueryOutputValue()
  p0 = ( pxy[0], 0. ) 
  p1 = ( pxy[1], 0. ) 

for db in DB:
  OpenDatabase(db)

  ierr = AddPlot("Pseudocolor", Var)
  if ierr==1:
    DrawPlots()

    # Do a lineout on variable to produce curve.
    Lineout(p0, p1, ("default"))

# extract data into CSV format
SetActiveWindow(2)
output = open("tip_position.csv",'w')

for k in range(len(DB)):
  SetActivePlots(k)
  vals = GetPlotInformation()["Curve"]

  Query("Time")
  t0 = GetQueryOutputValue()
  print "# time: %g" % (t0)

  xm=0.
  xp=0.
  ym=0.
  yp=0.
  for i in range(len(vals) / 2):
    if vals[2*i+1]>0.5:
      xm=p0[0]+vals[2*i]
      ym=vals[2*i+1]
    if vals[2*i+1]<0.5:
      xp=p0[0]+vals[2*i]
      yp=vals[2*i+1]
      slope=(yp-ym)/(xp-xm)
      x0=(0.5-ym)/slope+xm
      output.write("%g, %g\n" % (t0,x0))
      break

sys.exit()
