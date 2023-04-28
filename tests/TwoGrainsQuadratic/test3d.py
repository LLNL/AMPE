#!/usr/bin/env python
import sys
import subprocess
import os

print("Test Two grains with Quadratic energies...")

#prepare initial conditions file
initfilename="2spheres.nc"
subprocess.call(["python3", "../../utils/make_nuclei.py",
  "--nx", "64", "--ny", "64", "--nz", "48", "-r", "8",
  "--concentration-in", "0.1", "--concentration-out", "0.06",
  "--ngrains", "2", "-q", "4", 
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
  if line.count(b'grain') and line.count(b'Volume'):
    print(line)
    words=line.split()
    volume=eval(words[5])
    volumes.append(volume)

  if line.count(b'cycle'):
    print(line)
    words=line.split()
    time=eval(words[6])
    if time>0.08:
      end_reached = True

  if line.count(b'fraction'):
    print(line)
    if end_reached:
      words=line.split()
      sfraction=eval(words[6])
      if abs(sfraction-0.13)>1.e-2:
        print("Wrong solid fraction")
        sys.exit(1)

minv=1000.
maxv=0.
for v in volumes:
  if v<minv:
    minv = v
  if v>maxv:
    maxv = v

expected_value=2.0
if abs(maxv-expected_value)>0.01:
  print("Expected maxv = {}, found {}".format(expected_value,maxv))
  sys.exit(1)

expected_value=0.179
if abs(minv-expected_value)>0.001:
  print("Expected minv = {}, found {}".format(expected_value,minv))
  sys.exit(1)


os.remove(initfilename)

if end_reached:
  sys.exit(0)
else:
  sys.exit(1)

