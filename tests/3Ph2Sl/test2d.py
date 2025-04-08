#!/usr/bin/env python
import sys
import subprocess
import os

print("Test Binary3Ph2Sl...")

#prepare initial conditions file
subprocess.call(["python3", "../../tests/FolchPlappCALPHAD/make_initial.py", "-d", "2",
  "-x", "32", "-y", "32", "-z", "1", "--solid-fraction", "0.5",
  "--concL", "0.82", "--concA", "0.97", "--concB", "0.68",
  "test.nc"])

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]
thermdatadir = sys.argv[6]

#make symbolic link to calphad database
calphad_data = "calphadAlCuLFccTheta.json"
src = thermdatadir+'/'+calphad_data
try:
  os.symlink(src, calphad_data)
except FileExistsError:
  os.remove(calphad_data)
  os.symlink(src, calphad_data)

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

os.remove("test.nc")
os.unlink(calphad_data)

#analyse AMPE standard output
lines=output.split(b'\n')

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
if abs(f0-0.53)>tol:
  print("Final f0 = {} is incorrect".format(f0))
  sys.exit(1)
if abs(f1-0.22)>tol:
  print("Final f1 = {} is incorrect".format(f1))
  sys.exit(1)
if abs(f2-0.25)>tol:
  print("Final f2 = {} is incorrect".format(f2))
  sys.exit(1)

#check target time is reached
if time<2.2:
  print("Final time not reached")
  sys.exit(1)

sys.exit(0)
