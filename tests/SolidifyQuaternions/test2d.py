#!/usr/bin/env python
import sys
import subprocess
import os

print("Test solidify quaternions...")

#prepare initial conditions file
initfilename="test.nc"
subprocess.call(["python3", "../../utils/make_initial_grains_on_boundary.py",
  "-x", "64", "-y", "32", "-z", "1", "--solid-fraction", "0.25",
  "--smooth", "0", "--qlen", "4", "--ngrains", "2",
  initfilename])

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

#analyse AMPE standard output
lines=output.split(b'\n')

end_reached = False
volumes=[]

target_fs = 0.42
for line in lines:
  if line.count(b'cycle'):
    print(line)
    words=line.split()
    time=eval(words[6])
    if time>1.:
      end_reached = True
      print("End time reached")

  if end_reached:
    if line.count(b'grain') and line.count(b'Volume'):
      print(line)
      words=line.split()
      volume=eval(words[5])
      volumes.append(volume)

  if end_reached:
    if line.count(b'fraction'):
      print(line)
      words=line.split()
      sfraction=eval(words[6])
      if abs(sfraction-target_fs)>1.e-2:
        print("Wrong solid fraction, expected {}!".format(target_fs))
        sys.exit(1)

os.remove(initfilename)

if end_reached:
  if len(volumes) != 2:
    print("Expected two grains, found {}".format(len(volumes)))
    sys.exit(1)

if end_reached:
  sys.exit(0)
else:
  print("End time not reached.")
  sys.exit(1)

