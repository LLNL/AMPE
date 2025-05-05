#!/usr/bin/env python
import sys
import subprocess
import os

print("Test porosity measure in sintering...")

mpicmd  = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe     = sys.argv[4]
inp     = sys.argv[5]
datadir = sys.argv[6]

#make symbolic link to calphad data
data = "3spheres.csv"
src = datadir+'/'+data
try:
  os.symlink(src, data)
except FileExistsError:
  os.remove(data)
  os.symlink(src, data)

#prepare initial conditions file
initfilename="3spheres.nc"
subprocess.call(["python3", "../../utils/make_multi_spheres.py",
  "--nx", "96", "--ny", "96", "--nz", "1",
  "--concentration-A", "1.,0.", "--concentration-B", "0.,1.",
  "--concentration-out", "0.,0.",
  "--spheres", data,
  initfilename])

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

os.remove(initfilename)
os.unlink(data)

#analyse AMPE standard output
lines=output.split(b'\n')
densities=[]

end_reached = False
end_time = 0.006
for line in lines:

  if line.count(b'cycle'):
    print(line)
    words=line.split()
    time=eval(words[6])
    if time>end_time:
      end_reached = True

  if line.count(b'Density'):
    print(line)
    words=line.split()
    density=eval(words[2])
    densities.append(density)

mind=1.
maxd=0.
for d in densities:
  if d<mind:
    mind = d
  if d>maxd:
    maxd = d

expected_value=0.936
if abs(maxd-expected_value)>0.001:
  print("Expected max density = {}, found {}".format(expected_value,maxd))
  sys.exit(1)

expected_value=0.895
if abs(mind-expected_value)>0.001:
  print("Expected min density = {}, found {}".format(expected_value,mind))
  sys.exit(1)

if end_reached:
  sys.exit(0)
else:
  print("End time not reached!")
  sys.exit(1)
