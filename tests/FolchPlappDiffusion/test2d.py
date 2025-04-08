#!/usr/bin/env python
import sys
import subprocess
import os

print("Test FolchPlappDiffusion...")

#prepare initial conditions file
subprocess.call(["python3", "../../tests/FolchPlappDiffusion/make_initial.py", "-d", "2",
  "-x", "32", "-y", "32", "-z", "1",
  "--concB", "0.8",
  "test.nc"])

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]
thermdatadir = sys.argv[6]

#make symbolic link to calphad database
calphad_data = "calphad3phases.json"
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
cmax=20.
for line in lines:
  if line.count(b'cycle'):
    print(line)
    words=line.split()
    time = eval(words[6])
  if line.count(b'Max.') and  line.count(b'concentration'):
    print(line)
    words=line.split()
    cmax = eval(words[3])

#check target time is reached
if time<1000.:
  print("Final time not reached")
  sys.exit(1)
if cmax>0.1:
  print("cmax did not decrease as expected")
  sys.exit(1)

sys.exit(0)
