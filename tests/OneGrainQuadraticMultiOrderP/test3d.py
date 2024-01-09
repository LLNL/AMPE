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
if not os.path.exists(data):
  src = datadir+'/'+data
  print("Create symlink {}".format(src))
  os.symlink(src, data)

#prepare initial conditions file
initfilename="sphere.nc"
subprocess.call(["python3", "../../utils/make_multi_spheres.py",
  "--nx", "48", "--ny", "48", "--nz", "48",
  "--concentration-in", "0.1", "--concentration-out", "0.06",
  "--spheres", data,
  initfilename])

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

#analyse AMPE standard output
lines=output.split(b'\n')
volumes=[]

end_reached = False
for line in lines:

  if line.count(b'cycle'):
    print(line)
    words=line.split()
    time=eval(words[6])
    if time>0.15:
      end_reached = True

  if line.count(b'fraction'):
    print(line)
    if end_reached:
      words=line.split()
      sfraction=eval(words[6])
      if abs(sfraction-0.04)>1.e-2:
        print("Wrong solid fraction")
        sys.exit(1)

os.remove(initfilename)
os.unlink(data)

if end_reached:
  sys.exit(0)
else:
  sys.exit(1)
