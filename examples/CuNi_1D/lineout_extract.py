#usage:
#  visit -nowin -cli -verbose -s lineout_extract.py > lineout_extract.dat
import sys, string, os

maindir =os.getcwd()

# Adapt these parameters
visitdir="/v.bilayer/"
filename=maindir+"/bilayer.log"

#end points coordinates in [um]
p0 = ( 0.0,0.)
p1 = ( 2.0,0.)

t=3.e-5

frames=[]
myfile=open(filename,'r')
lines=myfile.readlines()

#figure out dimension of computational domain
found=0
for line in lines:
  if string.count(line, 'x_lo'):
    w=string.split(line)[2]
    x_lo=eval(w.replace(',',''))
    found=found+1
  if string.count(line, 'x_up'):
    w=string.split(line)[2]
    x_up=eval(w.replace(',',''))
    found=found+1
  if found>1:
    break;

for i,line1 in enumerate(lines):
  num_matches = string.count(line1, 'cycle')
  if num_matches:
    w=string.split(line1)
    tc1=eval(w[6])
    if tc1>t:
      cycle=w[2]
      while len(cycle)<5:
        cycle="0"+cycle
      frames.append(cycle)

      #now loop over following lines
      count=0
      f1=0.
      for line2 in lines[i:]:
        num_matches = string.count(line2, 'cycle')
        if num_matches:
          w=string.split(line2)
          tc2=eval(w[6])
        num_matches = string.count(line2, 'fraction')
        if num_matches:
          w=string.split(line2)
          f2=eval(w[6])
          count=count+1
          if count>1:
            print '#velocity=',(x_up-x_lo)*(f2-f1)/(tc2-tc1),'um/s'
            print '#interface at x=',x_lo+(x_up-x_lo)*f1
            break
          else:
            f1=f2
      break

#print frames

#load data
datadir=maindir+visitdir
DB=[]
for frame in frames:
  DB.append(datadir+"visit_dump."+frame+"/summary.samrai")

Vars=[]
Vars.append("phase")
Vars.append("concentration0")
Vars.append("conc_l0")
Vars.append("conc_a0")
Vars.append("driving_force")

#plot data
OpenDatabase(DB[0])
for Var in Vars:

  AddPlot("Pseudocolor", Var)

  DrawPlots()

  # Do a lineout on variable to produce curve.
  Lineout(p0, p1, ("default"))

# extract data into (x,y) format
SetActiveWindow(2)
nplots=len(Vars)
set_of_vals=[]
for j in range(nplots):
  SetActivePlots(j)
  set_of_vals.append(GetPlotInformation()["Curve"])

  Query("Time")
  t0 = GetQueryOutputValue()

print "#time: %g" % t0

# Write data as "x  y1 y2 y3 ..."
tol=1.e-6
for i in range(len(set_of_vals[0])/2):
  vals = set_of_vals[0]
  if vals[2*i+1]>tol and vals[2*i+1]<1-tol:
    print "%g" % (p0[0]+vals[2*i]),
    for vals in set_of_vals:
      print "  %g" %  (vals[2*i+1]),
    print '\n',

sys.exit()
