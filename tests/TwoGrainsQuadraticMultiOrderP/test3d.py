#!/usr/bin/env python
import sys
import subprocess
import os

print("Test Two grains with Quadratic energies...")

mpicmd  = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe     = sys.argv[4]
inp     = sys.argv[5]
datadir = sys.argv[6]

#make symbolic link to calphad data
data = "2spheres.csv"
src = datadir+'/'+data
try:
  os.symlink(src, data)
except FileExistsError:
  os.remove(data)
  os.symlink(src, data)

#prepare initial conditions file
initfilename="2spheres.nc"
subprocess.call(["python3", "../../utils/make_multi_spheres.py",
  "--nx", "64", "--ny", "64", "--nz", "64",
  "--concentration-A", "0.1", "--concentration-out", "0.06",
  "--spheres", data,
  initfilename])

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

os.remove(initfilename)
os.unlink(data)

#analyse AMPE standard output
lines=output.split(b'\n')
volfractions=[]

end_reached = False
end_time = 0.08
for line in lines:
  if (line.count(b'phase 0') or line.count(b'phase 1') ) and line.count(b'Volume'):
    print(line)
    words=line.split()
    volume=eval(words[6])
    volfractions.append(volume)

  if line.count(b'cycle'):
    print(line)
    words=line.split()
    time=eval(words[6])
    if time > end_time:
      end_reached = True

minv=1.
maxv=0.
for v in volfractions:
  if v<minv:
    minv = v
  if v>maxv:
    maxv = v

expected_value=0.099
if abs(maxv-expected_value)>0.003:
  print("Expected maxv = {}, found {}".format(expected_value,maxv))
  sys.exit(1)

expected_value=0.009
if abs(minv-expected_value)>0.001:
  print("Expected minv = {}, found {}".format(expected_value,minv))
  sys.exit(1)

if end_reached:
  sys.exit(0)
else:
  print("End time not reached!")
  sys.exit(1)
