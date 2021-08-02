#!/usr/bin/env python
import sys
import subprocess
import os

print("Test ThreePhases...")

#prepare initial conditions file
subprocess.call(["python3", "../../tests/ThreePhases/make_initial.py", "-d", "3",
  "-x", "64", "-y", "64", "-z", "32", "--solid-fraction", "0.5",
  "test.nc"])

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

#analyse AMPE standard output
lines=output.split(b'\n')

os.remove("test.nc")
mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]

time = 0.
for line in lines:
  if line.count(b'cycle'):
    print(line)
    words=line.split()
    time = eval(words[6])

if time<500.:
  print("Final time not reached")
  sys.exit(1)

sys.exit(0)
