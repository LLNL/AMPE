#!/usr/bin/env python
import sys
import subprocess
import os

print("Test AlCu parabolic...")

#prepare initial conditions file
initfilename="64x64x64.nc"
subprocess.call(["python3", "../../utils/make_nuclei.py",
  "--nx", "64", "--ny", "64", "--nz", "64", "-r", "10",
  "--center0", "0, 0, 0",
  "--concentration-in", "0.03", "--concentration-out", "0.08",
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
end_time = 0.0004
first_concentration=-1.
target_solid_fraction = 0.192494
for line in lines:
  if line.count(b'Integral'):
    words=line.split()
    concentration=eval(words[3])
    if first_concentration<0.:
      first_concentration=concentration
    #check conservation of mass
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

  if line.count(b'fraction'):
    print(line)
    if end_reached:
      words=line.split()
      sfraction=eval(words[6])
      if abs(sfraction-target_solid_fraction)>2.e-4:
        print("Wrong solid fraction")
        print("{} instead of {}".format(sfraction,target_solid_fraction))
        sys.exit(1)

os.remove(initfilename)

if end_reached:
  sys.exit(0)
else:
  print("Target end time {} not reached!".format(end_time))
  sys.exit(1)
