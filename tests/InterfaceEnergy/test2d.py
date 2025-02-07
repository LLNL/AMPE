#!/usr/bin/env python
import sys
import subprocess
from math import sqrt
from math import pi

print("Test interface energy...")

subprocess.call(["python3", "../../utils/make_sphere.py",
  "-x", "64", "-y", "64", "-z", "1", "-r", "32",
  "--cx", "0", "--cy", "0", "--qlen", "1",
  "test.nc"])

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

#analyse AMPE standard output
lines=output.split(b'\n')

domainl=3.2
tol_interface_energy=2.e-2
tol_bulk_energy=2.e-2
end_time=4.
end_sf=0.63
end_reached = False
time=0.

for line in lines:
  if line.count(b'cycle'):
    print(line)
    words=line.split()
    time=eval(words[6])
    print("time: ",time)
    if eval(words[6])>end_time:
      end_reached = True
  if line.count(b'bulk'):
    words=line.split()
    ebulk=eval(words[4])
    print("ebulk: ",ebulk)
  if line.count(b'interface'):
    words=line.split()
    eint=eval(words[4])
    print("eint: ",eint)
  if line.count(b'fraction'):
    words=line.split()
    sf=eval(words[6])
    print("solid fraction: ",sf)

    area=domainl*domainl*sf
    print("area: ",area)

    liquid_area=domainl*domainl*(1.-sf)
    print("liquid_area: ", liquid_area)

    circle_radius = sqrt(4.*area/pi)
    print("circle_radius: ", circle_radius)

    arc_length = 0.5*pi*circle_radius
    print("arc_length: ", arc_length)

    #start testing when interface profile reach "steady" shape
    if time>2:
      if(abs(eint-arc_length)>tol_interface_energy):
        print("Interface energy not matching arc length!")
        sys.exit(1)
      if(abs(ebulk-liquid_area)>tol_bulk_energy):
        print("Bulk energy not matching liquid area!")
        sys.exit(1)

    if end_reached:
      if abs(sf-end_sf)>1.e-2:
        print("Wrong solid fraction")
        sys.exit(1)

if end_reached:
  sys.exit(0)
else:
  print("End time not reached!")
  sys.exit(1)
