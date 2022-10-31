#!/usr/bin/env python
import sys
import subprocess
import os

print("Test Binary2Ph1Sl...")

#prepare initial conditions file
subprocess.call(["python3", "../../tests/2Ph1Sl/make_initial.py", "-d", "2",
  "-x", "16", "-y", "64", "-z", "4", "--solid-fraction", "0.2",
  "--concL", "0.825", "--concA", "0.68",
  "test.nc"])

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]
thermdatadir = sys.argv[6]

#make symbolic link to calphad database
calphad_data = "calphadAlCuLTheta.json"
if not os.path.exists(calphad_data):
  src = thermdatadir+'/'+calphad_data
  os.symlink(src, calphad_data)

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

#analyse AMPE standard output
lines=output.split(b'\n')

os.remove("test.nc")

time = 0.
f0=0.
for line in lines:
  if line.count(b'cycle'):
    print(line)
    words=line.split()
    time = eval(words[6])
  if line.count(b'Volume fraction'):
    print(line)
    words=line.split()
    f0= eval(words[6])

#check phase fraction
tol=0.01
if abs(f0-0.25)>tol:
  print("Final phase fraction = {} is incorrect".format(f0))
  sys.exit(1)

#check target time is reached
if time<5.e-5:
  print("Final time not reached")
  sys.exit(1)

sys.exit(0)
