#!/usr/bin/env python
import sys
import subprocess
import os

print("Test parabolic multi-order three phases...")

#prepare initial conditions file
initfilename="160x120.nc"
subprocess.call(["python3", "../../tests/ParabolicMultiOrderThreePhases/make_initial.py",
  "--nx", "160", "--ny", "120", "--nz", "1",
  "--concL", "0.5", "--concB", "0.75", "--concA", "0.25",
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
end_time = 3.e-3
final_lfraction = 0.647
lfraction_checked = False

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
    if time>end_time:
      end_reached = True
      dt=eval(words[10])
      if (dt-1.e-6)<0.:
        print("Wrong dt: too small")
        sys.exit(1)

  if line.count(b'fraction of phase 2'):
    print(line)
    if end_reached:
      lfraction_checked = True
      words=line.split()
      lfraction=eval(words[6])
      print("Final liquid fraction: {}".format(lfraction))
      if abs(lfraction-final_lfraction)>1.e-3:
        print("Wrong liquid fraction")
        sys.exit(1)

  if line.count(b'fraction of phase 0'):
    print(line)
    fraction0=eval(words[6])
  if line.count(b'fraction of phase 1'):
    print(line)
    fraction1=eval(words[6])

if abs(fraction0-fraction1)>1.e-3:
  print("fraction0 = {}".format(fraction0))
  print("fraction1 = {}".format(fraction1))
  print("Phase 0 and 1 should have the same fraction of the domain")
  sys.exit(1)



os.remove(initfilename)

if end_reached and lfraction_checked:
  sys.exit(0)
else:
  print("End time not reached...")
  sys.exit(1)
