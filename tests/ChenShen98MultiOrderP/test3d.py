#!/usr/bin/env python
import sys
import subprocess
import os

print("Test ChenShen98 multi-order parameters...")

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]
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
initfilename = "32x32x32.nc"
subprocess.call(["python3", "../../utils/make_multi_spheres.py",
  "-x", "32", "-y", "32", "-z", "32",
  "--spheres", data,
  initfilename])

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

os.remove(initfilename)
os.unlink(data)

#analyse AMPE standard output
lines=output.split(b'\n')

end_reached = False
for line in lines:
  if line.count(b'cycle'):
    print(line)
    words=line.split()
    if eval(words[6])>120.:
      end_reached = True
  if line.count(b'fraction') and line.count(b'phase 0'):
    print(line)
    if end_reached:
      words=line.split()
      if abs(eval(words[6])-0.040)>1.e-2:
        print("Wrong solid fraction")
        sys.exit(1)

if end_reached:
  sys.exit(0)
else:
  print("End time not reached!")
  sys.exit(1)

