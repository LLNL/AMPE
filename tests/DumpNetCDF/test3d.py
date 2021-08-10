#!/usr/bin/env python
import sys
import subprocess
import os

print("Test DumpNetCDF...")

#prepare initial conditions file
subprocess.call(["python3", "../../utils/make_bilayer.py", "-2",
  "-x", "32", "-y", "16", "-z", "16", "-r", "5",
  "-c", "0.05", "--concentration-in", "0.04",
  "--quat-in", "1,0,0,0", "--quat-in-two", "0.95,0.31225,0.,0.",
  "--quat-out", "0,0,1,0", "--qlen", "4",
  "initial.nc"])

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]

#run AMPE
print("Run AMPE...")
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

#analyse AMPE standard output
lines=output.split(b'\n')

end_reached = False
count=0
for line in lines:
  if line.count(b'cycle'):
    print(line)
    count=count+1
  if line.count(b'Wrote'):
    end_reached = True

if count!=10:
  print("1st run: count!=10")
  sys.exit(1)

if not end_reached:
  print("NetCDF dump not reached")
  sys.exit(1)

print("Rename NetCDF file...")
os.rename('final.nc','initial.nc')

print("Run AMPE again...")
output = subprocess.check_output(command,shell=True)

#analyse AMPE standard output
lines=output.split(b'\n')

end_reached = False
count=0
for line in lines:
  if line.count(b'cycle'):
    print (line)
    count=count+1
  if line.count(b'Wrote'):
    end_reached = True

if count!=10:
  print("2nd run: count!=10")
  sys.exit(1)

if end_reached:
  sys.exit(0)
else:
  print("NetCDF dump not reached")
  sys.exit(1)
