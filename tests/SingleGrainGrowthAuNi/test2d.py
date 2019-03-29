#!/usr/bin/env python
import sys
import subprocess
import os

print("Test SingleGrainGrowthAuNi...")

#prepare initial conditions file
initfilename="64x64.nc"
subprocess.call(["python3", "../../utils/make_nuclei.py",
  "--nx", "64", "--ny", "64", "--nz", "1", "-r", "16",
  "--center0", "0, 0, 0",
  "-c", "0.25", "--concentration-in", "0.096",
  initfilename])

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]
thermdatadir = sys.argv[6]

#make symbolic link to calphad database
calphad_data = "calphadAuNi.dat"
if not os.path.exists(calphad_data):
  src = thermdatadir+'/'+calphad_data
  os.symlink(src, calphad_data)

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

#analyse AMPE standard output
lines=output.split(b'\n')

end_reached = False
oldconcentration=-1.
for line in lines:
  num_matches = line.count(b'Integral')
  if num_matches:
    words=line.split()
    concentration=eval(words[3])
    if oldconcentration<0.:
      oldconcentration=concentration
    if abs(concentration-oldconcentration)>1.e-4:
      sys.exit(1)

  num_matches = line.count(b'cycle')
  if num_matches:
    print(line)
    words=line.split()
    time=eval(words[6])
    if time>0.3:
      end_reached = True
      dt=eval(words[10])
      if (dt-0.0012)<0.:
        print("Wrong dt")
        sys.exit(1)
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
  sys.exit(1)

