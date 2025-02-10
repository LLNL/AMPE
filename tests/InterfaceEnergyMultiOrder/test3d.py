#!/usr/bin/env python
import sys
import subprocess
from math import pi

print("Test interface energy...")

subprocess.call(["python3", "../../utils/make_sphere.py",
  "-x", "32", "-y", "32", "-z", "32", "-r", "28",
  "--cx", "0", "--cy", "0", "--cz", "0", "--qlen", "4",
  "test.nc"])

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

#analyse AMPE standard output
lines=output.split(b'\n')

domainl=1.6
tol_interface_energy=5.e-2
tol_bulk_energy=5.e-2
end_time=0.6
end_sf=0.28
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
  if line.count(b'fraction of phase 0'):
    words=line.split()
    sf=eval(words[6])
    print("solid fraction: ",sf)

    volume=domainl*domainl*domainl*sf
    print("volume: ",volume)

    liquid_volume=domainl*domainl*domainl*(1.-sf)
    print("liquid_volume: ", liquid_volume)

    circle_radius = (3.*8.*volume/(4.*pi))**(1./3.)
    print("circle_radius: ", circle_radius)

    eight_sphere_vol=0.125*4.*pi*circle_radius*circle_radius*circle_radius/3.
    print("eight_sphere_vol: ", eight_sphere_vol)

    # 1/8 sphere surface
    surface = 0.5*pi*circle_radius*circle_radius
    print("surface: ", surface)

    #start testing when interface profile reach "steady" shape
    if time>0.4:
      if(abs(eint-surface)>tol_interface_energy):
        print("Interface energy not matching surface area!")
        sys.exit(1)
      if(abs(10.*ebulk-liquid_volume)>tol_bulk_energy):
        print("Bulk energy not matching liquid volume!")
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
