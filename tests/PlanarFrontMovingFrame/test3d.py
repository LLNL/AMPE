#!/usr/bin/env python
import sys
import subprocess
import os

print("Test PlanarFrontMovingFrame...")

#prepare initial conditions file
subprocess.call(["python3", "../../utils/make_initial_grains_on_boundary.py",
  "--ngrains", "1", "-d", "3",
  "-x", "256", "-y", "32", "-z", "32", "--solid-fraction", "0.5",
  "--width", "1", "--plane", "0",
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

previous_time = -1.

tol = 1.e-3

ly = 400. #length of domain
sfraction = 0.5

for line in lines:
  num_matches = line.count(b'cycle')
  if num_matches:
    print(line)
    words=line.split()
    time = eval(words[6])

  num_matches = line.count(b'fraction')
  if num_matches:
    print(line)
    words=line.split()
    previous_sfraction = sfraction
    sfraction=eval(words[6])

    delta_sfraction = (sfraction-previous_sfraction)
    print("delta_sfraction={}".format(delta_sfraction))

    if time>50:
      if abs(delta_sfraction)>tol:
        print("Solid fraction is still changing!")
        sys.exit(1)

sys.exit(0)

