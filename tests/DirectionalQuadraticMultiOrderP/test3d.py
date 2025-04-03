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
src = datadir+'/'+data
print("Create symlink {}".format(src))
try:
  os.symlink(src, data)
except FileExistsError:
  os.remove(data)
  os.symlink(src, data)

#prepare initial conditions file
initfilename="512x4x4.nc"
subprocess.call(["python3", "../../utils/make_multi_spheres.py",
  "--nx", "512", "--ny", "4", "--nz", "4",
  "--concentration-A", "0.136", "--concentration-out", "0.147",
  "--spheres", data,
  initfilename])

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
print(command)
output = subprocess.check_output(command,shell=True)

os.remove(initfilename)
os.unlink(data)

#analyse AMPE standard output
lines=output.split(b'\n')

fs = 0.
time = 0.
for line in lines: ## loop over lines of file
  if line.count(b'cycle #'):
    words = line.split()
    time = eval(words[6])

  if line.count(b'Volume fraction') and line.count(b'phase 0'):
    words = line.split()
    fs = eval(words[6])

  if line.count(b'Max') and line.count(b'concentration'):
    words = line.split()
    cmax = eval(words[3])

target_cmax = 0.163
if abs(cmax-target_cmax)>0.002:
  print("Wrong cmax {}, expected {}".format(cmax,target_cmax))
  sys.exit(1)

target_vol_fraction = 0.53
print("Solid fraction : {}".format(fs))
if abs(fs-target_vol_fraction)>0.01:
  print("Wrong solid fraction {}, expected {}".format(fs,target_vol_fraction))
  sys.exit(1)

target_time = 2.e-2
print("Time reached : {}".format(time))
if time<target_time:
  print("Target time {} not reached!".format(target_time))
  sys.exit(1)

sys.exit(0)
