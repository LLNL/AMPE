#!/usr/bin/env python
import sys
import subprocess
import os

print("Test Four grain orientations...")

#prepare initial conditions file
initfilename="test.nc"
subprocess.call(["python3", "../../utils/make4corners.py",
  "-x", "64", "-y", "64", "-z", "1",
  initfilename])

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
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
    if time>0.3:
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
      if abs(sfraction-0.93)>1.e-2:
        print("Wrong solid fraction")
        sys.exit(1)

tol=0.0001
for v in volumes:
  if abs(v-0.0028)>tol and abs(v-0.0087)>tol:
    print("Expected volume of 0.0028 or 0.0087, found {}".format(v))
    sys.exit(1)

os.remove(initfilename)

if end_reached:
  sys.exit(0)
else:
  sys.exit(1)

