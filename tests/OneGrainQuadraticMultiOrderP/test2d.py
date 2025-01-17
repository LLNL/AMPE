#!/usr/bin/env python
import sys
import subprocess
import os

print("Test One grains with multiple order parameters...")

mpicmd  = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe     = sys.argv[4]
inp     = sys.argv[5]
datadir = sys.argv[6]

#make symbolic link to input data
data = "spheres.csv"
src = datadir+'/'+data
print("Create symlink {}".format(src))
try:
  os.symlink(src, data)
except FileExistsError:
  os.remove(data)
  os.symlink(src, data)

#prepare initial conditions file
initfilename="sphere.nc"
subprocess.call(["python3", "../../utils/make_multi_spheres.py",
  "--nx", "64", "--ny", "64", "--nz", "1",
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
volumes=[]

end_reached = False
target_sf = 0.205
for line in lines:

  if line.count(b'cycle'):
    print(line)
    words=line.split()
    time=eval(words[6])
    if time>0.25:
      end_reached = True

  if line.count(b'fraction')  and line.count(b'phase 0'):
    print(line)
    if end_reached:
      words=line.split()
      sfraction=eval(words[6])
      if abs(sfraction-target_sf)>1.e-2:
        print("Wrong solid fraction:")
        print("found {}, expected {}".format(sfraction, target_sf))
        sys.exit(1)

if end_reached:
  sys.exit(0)
else:
  print("End time not reached")
  sys.exit(1)
