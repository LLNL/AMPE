#!/usr/bin/env python
import sys
import subprocess
import os

print("Test One grain with Quadratic energies...")

#prepare initial conditions file
initfilename="nuclei.nc"
subprocess.call(["python3", "../../utils/make_nuclei.py",
  "--nx", "48", "--ny", "48", "--nz", "48", "-r", "8",
  "--concentration-in", "0.1", "--concentration-out", "0.06",
  "--ngrains", "1",
  initfilename])

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
print(command)
output = subprocess.check_output(command,shell=True)

#analyse AMPE standard output
lines=output.split(b'\n')
volumes=[]

end_reached = False
for line in lines:

  if line.count(b'cycle'):
    print(line)
    words=line.split()
    time=eval(words[6])
    if time>0.15:
      end_reached = True

  if line.count(b'fraction'):
    print(line)
    if end_reached:
      words=line.split()
      sfraction=eval(words[6])
      if abs(sfraction-0.14)>1.e-2:
        print("Wrong solid fraction")
        sys.exit(1)

os.remove(initfilename)

if end_reached:
  sys.exit(0)
else:
  sys.exit(1)

