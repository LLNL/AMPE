#!/usr/bin/env python
import sys
import subprocess
import os

print("Test Gaussian Temperature profile...")

#prepare initial conditions file
initfilename="64x64.nc"
subprocess.call(["python3", "../../utils/make_nuclei.py",
  "--nx", "64", "--ny", "64", "--nz", "1", "-r", "20",
  "--concentration-in", "0.11", "--concentration-out", "0.28",
  initfilename])

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]
thermdatadir = sys.argv[6]

#make symbolic link to calphad database
calphad_data = "calphadAuNi.dat"
src = thermdatadir+'/'+calphad_data
try:
  os.symlink(src, calphad_data)
except FileExistsError:
  os.remove(calphad_data)
  os.symlink(src, calphad_data)

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

os.remove(initfilename)
os.unlink(calphad_data)

#analyse AMPE standard output
lines=output.split(b'\n')

for line in lines:
  if line.count(b'cycle'):
    print(line)

#solid fraction
solid_fraction = 0.
for line in lines:
  if line.count(b'fraction'):
    print(line)
    words=line.split()
    solid_fraction = eval(words[6])

tol = 1.e-3
expected_solid_fraction = 0.31
if abs(solid_fraction-expected_solid_fraction)>tol:
  print("solid fraction = {}".format(solid_fraction))
  print("Final solid fraction differs from expected value")
  sys.exit(1)

#final time
final_time = 0.
for line in lines:
  if line.count(b'cycle'):
    words=line.split()
    final_time = eval(words[6])
print("final time = {}".format(final_time))

tol = 1.e-3
expected_final_time = 0.04
if abs(final_time-expected_final_time)>tol:
  print(final_time)
  print("Final time differs from expected value")
  sys.exit(1)

#max. temperature
maxT = 0.
for line in lines:
  if line.count(b'Temperature') and line.count(b'Max.'):
    print(line)
    words=line.split()
    maxT = eval(words[3])

tol = 1.e-1
expected_max_T = 1423. + final_time * 500.
if abs(maxT-expected_max_T)>tol:
  print("max. T = {}".format(maxT))
  print("expected value = {}".format(expected_max_T))
  print("Final Max. temperature differs from expected value")
  sys.exit(1)

#min. temperature
minT = 0.
for line in lines:
  if line.count(b'Temperature') and line.count(b'Min.'):
    print(line)
    words=line.split()
    minT = eval(words[3])

expected_min_T = 1423.
if abs(minT-expected_min_T)>tol:
  print("Final Min. temperature differs from expected value")
  sys.exit(1)

#number of cycles
for line in lines:
  if line.count(b'cycle'):
    words=line.split()
    cycle = eval(words[2])

expected_ncycles = 259
if cycle>expected_ncycles+10:
  print("Number of steps larger than expected")
  sys.exit(1)

sys.exit(0)
