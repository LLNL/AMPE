#!/usr/bin/env python
import sys
import subprocess
import os

print("Test Dendrite...")

#prepare initial conditions file
initfilename="3d.nc"
subprocess.call(["python3", "../../utils/make_nuclei.py",
  "--nx", "60", "--ny", "60", "--nz", "60", "-r", "10",
  "--center0", "0, 0, 0",
  "-w", "1.4",
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
for line in lines:
  num_matches = line.count(b'cycle')
  if num_matches:
    print(line)
    words=line.split()
    time=eval(words[6])
    if time>40.:
      end_reached = True

  num_matches = line.count(b'fraction')
  if num_matches:
    print(line)
    if end_reached:
      words=line.split()
      sfraction=eval(words[6])
      if abs(sfraction-0.15)>1.e-2:
        print("Wrong solid fraction")
        sys.exit(1)

os.remove(initfilename)

if end_reached:
  sys.exit(0)
else:
  sys.exit(1)

