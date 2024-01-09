#!/usr/bin/env python
import sys
import subprocess
import os

print("Test multiple order parameters in 1D...")

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]
datadir = sys.argv[6]

#make symbolic link to input data
data = "sphere.csv"
if not os.path.exists(data):
  src = datadir+'/'+data
  print("Create symlink {}".format(src))
  os.symlink(src, data)

#prepare initial conditions file
initfilename="512x4.nc"
subprocess.call(["python3", "../../utils/make_multi_spheres.py",
  "--nx", "512", "--ny", "4", "--nz", "1",
  "--concentration-in", "0.136", "--concentration-out", "0.147",
  "--spheres", data,
  initfilename])

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

input_file=open(inp,'r')
lines=input_file.readlines()

#analyse AMPE standard output
lines=output.split(b'\n')

fs = 0.
time = 0.
for line in lines: ## loop over lines of file
  if line.count(b'cycle #'):
    words = line.split()
    time = eval(words[6])

  if line.count(b'Volume fraction'):
    words = line.split()
    fs = eval(words[6])

target_vol_fraction = 0.513
print("Solid fraction : {}".format(fs))
if abs(fs-target_vol_fraction)>0.005:
  print("Wrong solid fraction {}".format(fs))
  sys.exit(1)

target_time = 3.e-2
print("Time reached : {}".format(time))
if time<target_time:
  print("Target time {} not reached!".format(target_time))
  sys.exit(1)

os.remove(initfilename)

sys.exit(0)

