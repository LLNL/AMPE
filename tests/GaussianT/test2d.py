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

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

#analyse AMPE standard output
lines=output.split(b'\n')

#solid fraction
solid_fraction = 0.
for line in lines:
  if line.count(b'fraction'):
    print(line)
    words=line.split()
    solid_fraction = eval(words[6])

tol = 1.e-3
expected_solid_fraction = 0.309
if abs(solid_fraction-expected_solid_fraction)>tol:
  printf("Final solid fraction differs from expected value")
  sys.exit(1)

#max. temperature
maxT = 0.
for line in lines:
  if line.count(b'Temperature') and line.count(b'Max.'):
    print(line)
    words=line.split()
    maxT = eval(words[3])

tol = 2.e-2
expected_max_T = 1443.01
if abs(maxT-expected_max_T)>tol:
  printf("Final Max. temperature differs from expected value")
  sys.exit(1)

#min. temperature
minT = 0.
for line in lines:
  if line.count(b'Temperature') and line.count(b'Min.'):
    print(line)
    words=line.split()
    maxT = eval(words[3])

expected_max_T = 1423.
if abs(maxT-expected_max_T)>tol:
  printf("Final Min. temperature differs from expected value")
  sys.exit(1)

#number of cycles
for line in lines:
  if line.count(b'cycle'):
    words=line.split()
    cycle = eval(words[2])

expected_ncycles = 259
if cycle>expected_ncycles+10:
  printf("Number of steps larger than expected")
  sys.exit(1)

os.remove(initfilename)

sys.exit(0)
