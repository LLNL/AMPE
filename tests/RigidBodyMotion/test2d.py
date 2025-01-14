#!/usr/bin/env python
import sys
import subprocess
import os

print("Test rigid body motion...")

mpicmd  = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe     = sys.argv[4]
inp     = sys.argv[5]
datadir = sys.argv[6]

#make symbolic link to calphad data
data = "1sphere.csv"
src = datadir+'/'+data
try:
  os.symlink(src, data)
except FileExistsError:
  os.remove(data)
  os.symlink(src, data)

#prepare initial conditions file
initfilename="1sphere.nc"
subprocess.call(["python3", "../../utils/make_multi_spheres.py",
  "--nx", "40", "--ny", "40", "--nz", "1",
  "--concentration-A", "1.", "--concentration-out", "0.",
  "--spheres", data,
  initfilename])

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

os.unlink(data)
os.remove(initfilename)

#analyse AMPE standard output
lines=output.split(b'\n')

volume=0.
end_reached = False
end_time = 0.02
for line in lines:

  if line.count(b'cycle'):
    print(line)
    words=line.split()
    time=eval(words[6])
    if time>end_time:
      end_reached = True

  if end_reached:
    if line.count(b'phase 0') and line.count(b'Volume'):
      print(line)
      words=line.split()
      volume=eval(words[6])

if end_reached:
  expected_value=0.147
  if abs(volume-expected_value)>0.001:
    print("Expected volume = {}, found {}".format(expected_value,volume))
    sys.exit(1)

if end_reached:
  sys.exit(0)
else:
  print("End time not reached!")
  sys.exit(1)
