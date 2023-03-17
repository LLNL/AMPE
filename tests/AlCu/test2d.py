#!/usr/bin/env python
import sys
import subprocess
import os

print("Test AlCu...")

#prepare initial conditions file
initfilename="160x160.nc"
subprocess.call(["python3", "../../utils/make_nuclei.py",
  "--nx", "160", "--ny", "160", "--nz", "1", "-r", "15",
  "--center0", "0, 0, 0",
  "--concentration-in", "0.003", "--concentration-out", "0.02",
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
first_concentration=-1.
for line in lines:
  if line.count(b'Integral'):
    words=line.split()
    concentration=eval(words[3])
    if first_concentration<0.:
      first_concentration=concentration
    if abs(concentration-first_concentration)>1.e-4:
      sys.exit(1)

  if line.count(b'cycle'):
    print(line)
    words=line.split()
    time=eval(words[6])
    if time>0.00695:
      end_reached = True
      dt=eval(words[10])
      if (dt-1.e-5)<0.:
        print("Wrong dt: too small")
        sys.exit(1)

  if line.count(b'fraction'):
    print(line)
    if end_reached:
      words=line.split()
      sfraction=eval(words[6])
      if abs(sfraction-0.21)>2.e-3:
        print("Wrong solid fraction")
        sys.exit(1)

os.remove(initfilename)

if end_reached:
  sys.exit(0)
else:
  sys.exit(1)

