#!/usr/bin/env python
import sys
import subprocess
import os

print("Test Cahn-Hilliard solver...")

#prepare initial conditions file
initfilename="64x64.nc"
subprocess.call(["python3", "../../tests/CahnHilliard/make_initial.py",
  "--nx", "64", "--ny", "64", "--nz", "1",
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
oldconcentration=-1.
end_time = 300.
for line in lines:
  num_matches = line.count(b'Integral')
  if num_matches:
    words=line.split()
    concentration=eval(words[3])
    print("Integral concentration = {}".format(concentration))
    if oldconcentration<0.:
      oldconcentration=concentration
    if abs(concentration-oldconcentration) > 1.e-5*concentration:
      print("Two consecutive integral of composition differ by more than relative tolerance")
      sys.exit(1)

  num_matches = line.count(b'cycle')
  if num_matches:
    print(line)
    words=line.split()
    time=eval(words[6])
    if time>end_time:
      end_reached = True

  num_matches = line.count(b'fraction')
  if num_matches:
    print(line)
    if end_reached:
      words=line.split()
      sfraction=eval(words[6])
      if abs(sfraction-0.32)>1.e-2:
        print("Wrong solid fraction")
        sys.exit(1)

os.remove(initfilename)

if end_reached:
  sys.exit(0)
else:
  print("End time {} not reached".format(end_time))
  sys.exit(1)
