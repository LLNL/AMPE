#!/usr/bin/env python
import sys
import subprocess
import os

print("Test AlCuMovingFrame...")

#prepare initial conditions file
subprocess.call(["python3", "../../utils/make_bilayer.py",
  "-r", "512", "-d", "1",
  "-x", "1024", "-y", "1", "-z", "1", "--centerx", "0.",
  "1024.nc"])

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

#analyse AMPE standard output
lines=output.split(b'\n')

os.remove("1024.nc")

previous_time = -1.

tol = 1.e-3

sfraction = 0.5
c0=0.
time=0.

for line in lines:
  if line.count(b'cycle'):
    print(line)
    words=line.split()
    time = eval(words[6])

  if line.count(b'concentration') and line.count(b'Max.'):
    words=line.split()
    c0 = eval(words[3]);

  if line.count(b'fraction'):
    print(line)
    words=line.split()
    previous_sfraction = sfraction
    sfraction=eval(words[6])

    delta_sfraction = (sfraction-previous_sfraction)
    print("delta_sfraction={}".format(delta_sfraction))

    if time>0.003:
      if abs(delta_sfraction)>tol:
        print("Solid fraction is still changing!")
        sys.exit(1)

maxc = 0.0211
print("Max. concentration: {}".format(c0))
if abs(c0-maxc)>1.e-4:
  print("Max. concentration incorrect, expected {}!!!".format(maxc))
  sys.exit(1)

sys.exit(0)

