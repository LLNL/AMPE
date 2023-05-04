#!/usr/bin/env python
import sys
import subprocess
import os

print("Test 1 grain Quadratic energies...")

#prepare initial conditions file
initfilename="3d.nc"
subprocess.call(["python3", "../../tests/ConservedVolume/make_square.py",
  "--nx", "32", "--ny", "32", "--nz", "32", "-r", "8",
  initfilename])

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

#analyse AMPE standard output
lines=output.split(b'\n')
sfractions=[]

end_reached = False
for line in lines:
  if line.count(b'cycle'):
    print(line)
    words=line.split()
    time=eval(words[6])
    if time>0.008:
      end_reached = True

  if line.count(b'fraction'):
    print(line)
    words=line.split()
    sfractions.append(eval(words[6]))

for i in range(4, len(sfractions)):
  print("fraction at step {} = {}".format(i,sfractions[i]))
  diff=sfractions[i-1]-sfractions[i]
  if abs(diff)>1.e-4:
    print("change in solid fraction exceeds tol")
    sys.exit(1)


os.remove(initfilename)

if end_reached:
  sys.exit(0)
else:
  print("end time not reached")
  sys.exit(1)

