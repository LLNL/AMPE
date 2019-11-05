#usage:
#  visit -cli -verbose -nowin -s phase_field_1500.py
#in main visit data directory
import sys, string, os

visitdir =os.getcwd()
filename=visitdir+"/dumps.visit"

myfile=open(filename,'r')
lines=myfile.readlines()

for line in lines:
  words=line.split()
  frame=words[0]

#load data
DB=visitdir+"/"+frame
print( "#{}".format(DB))


Var="phase"

OpenDatabase(DB)

DeleteAllPlots()
AddPlot("Contour", "phase", 1, 1)

# change Contour plot attributes
c=ContourAttributes()
c.contourNLevels=1
#c.contourMethod=Value
c.contourValue=(0.5)
SetPlotOptions(c)

DrawPlots()

e=ExportDBAttributes()
e.db_type = "XYZ"
e.variables = ("x", "y")
e.filename = "phase_contour"
ExportDatabase(e)

sys.exit()

