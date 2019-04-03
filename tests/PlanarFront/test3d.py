#!/usr/bin/env python
import sys
import subprocess
import os

print("Test PlanarFront...")

#prepare initial conditions file
subprocess.call(["python3", "../../utils/make_initial_grains_on_boundary.py",
  "--ngrains", "1", "-d", "3",
  "-x", "16", "-y", "256", "-z", "16", "--solid-fraction", "0.05",
  "--width", "1",
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

previous_sfraction = -1.
previous_time = -1.

tol_percent = 0.025 #2.5%
tol_velocity = 6.5*tol_percent

ly = 400. #length of domain
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
    sfraction=eval(words[6])

    if previous_sfraction>0.:
      velocity = ly*(sfraction-previous_sfraction)/(time-previous_time)
      print("velocity={}".format(velocity))

      if time>20:
        if abs(velocity-6.5)>tol_velocity:
          print("Wrong velocity!")
          sys.exit(1)

    previous_sfraction = sfraction
    previous_time      = time


sys.exit(0)

