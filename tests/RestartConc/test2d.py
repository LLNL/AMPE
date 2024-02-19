#!/usr/bin/env python
import sys
import subprocess
import os

print("Test ThreePhases...")

#prepare initial conditions file
subprocess.call(["python3", "../../tests/ThreePhasesCALPHAD/make_initial.py", "-d", "2",
  "-x", "32", "-y", "32", "-z", "1", "--solid-fraction", "0.5",
  "--concL", "0.5309", "--concA", "0.7686", "--concB", "0.2314",
  "test.nc"])

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]
thermdatadir = sys.argv[6]

#make symbolic link to calphad database
calphad_data = "calphad3phases.json"
if not os.path.exists(calphad_data):
  src = thermdatadir+'/'+calphad_data
  os.symlink(src, calphad_data)

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

#second run (restart)
command = "{} {} {} r.test 375".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

#analyse AMPE standard output
lines=output.split(b'\n')

os.remove("test.nc")

time = 0.
f0=0.
f1=0.
f2=0.
for line in lines:
  if line.count(b'cycle'):
    print(line)
    words=line.split()
    time = eval(words[6])
  if line.count(b'Volume fraction'):
    if line.count(b'phase 0'):
      print(line)
      words=line.split()
      f0= eval(words[6])
    if line.count(b'phase 1'):
      print(line)
      words=line.split()
      f1= eval(words[6])
    if line.count(b'phase 2'):
      print(line)
      words=line.split()
      f2= eval(words[6])

#check phase fractions
tol=0.01
if abs(f0-0.51)>tol:
  print("Final f0 = {} is incorrect".format(f0))
  sys.exit(1)
if abs(f1-0.24)>tol:
  print("Final f1 = {} is incorrect".format(f1))
  sys.exit(1)
if abs(f2-0.24)>tol:
  print("Final f2 = {} is incorrect".format(f2))
  sys.exit(1)

#check target time is reached
if time<130.:
  print("Final time not reached")
  sys.exit(1)

sys.exit(0)
