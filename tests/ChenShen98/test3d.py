#!/usr/bin/env python
import sys
import subprocess

print("Test ChenShen98...")

#prepare initial conditions file
subprocess.call(["python3", "../../utils/make_sphere.py",
  "-x", "32", "-y", "32", "-z", "32", "-r", "25",
  "--cx", "0", "--cy", "0", "--cz", "0", "--qlen", "1",
  "32x32x32.nc"])

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
    if eval(words[6])>120.:
      end_reached = True
      dt=words[10]
      if (eval(dt)-0.65)>1.e-2:
        print("Wrong dt")
        sys.exit(1)
  num_matches = line.count(b'fraction')
  if num_matches:
    print(line)
    if end_reached:
      words=line.split()
      if abs(eval(words[6])-0.035)>1.e-3:
        print("Wrong solid fraction")
        sys.exit(1)

if end_reached:
  sys.exit(0)
else:
  sys.exit(1)

