#!/usr/bin/env python
import sys
import subprocess
import os

print("Test FeNbNiWithNoise...")

#prepare initial conditions file
initfilename="test.nc"
subprocess.call(["python3", "../../utils/make_initial_grains_on_boundary.py",
  "--ngrains", "1", "-d", "2", "--fluctuation", "0.01",
  "-x", "32", "-y", "64", "-z", "1", "--solid-fraction", "0.2",
  "--width", "3",
  initfilename])

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]
thermdatadir = sys.argv[6]

#make symbolic link to calphad database
calphad_data = "calphadFeNbNi_Mathon_et_al.dat"
src = thermdatadir+'/'+calphad_data
try:
  os.symlink(src, calphad_data)
except FileExistsError:
  os.remove(calphad_data)
  os.symlink(src, calphad_data)

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

os.remove(initfilename)
os.unlink(calphad_data)

#analyse AMPE standard output
lines=output.split(b'\n')

end_reached = False
oldconcentration0=-1.
oldconcentration1=-1.

end_time=2.e-6
end_dt=8.e-9

nghosts=-1
for line in lines:
  if line.count(b'ghost') and line.count(b'cells'):
    print(line)
    words=line.split()
    nghosts=eval(words[2])
    print("Number of ghost cells: {}".format(nghosts))
    if nghosts!=2:
      sys.exit(1)
    else:
      break

if nghosts<0:
  sys.exit(1)

for line in lines:
  if line.count(b'Integral concentration 0'):
    print(line)
    words=line.split()
    concentration=eval(words[3])
    if oldconcentration0<0.:
      oldconcentration0=concentration
    if abs(concentration-oldconcentration0)>1.e-4:
      print("Concentration changing...")
      sys.exit(1)

  if line.count(b'Integral concentration 1'):
    print(line)
    words=line.split()
    concentration=eval(words[3])
    if oldconcentration1<0.:
      oldconcentration1=concentration
    if abs(concentration-oldconcentration1)>1.e-4:
      print("Concentration changing...")
      sys.exit(1)

  if line.count(b'cycle'):
    print(line)
    words=line.split()
    time=eval(words[6])
    if time>end_time:
      print(time)
      end_reached = True
      dt = eval(words[10])
      if dt<end_dt:
        print("dt = {} not reached!".format(end_dt))
        sys.exit(1)

  if line.count(b'fraction'):
    print(line)
    if end_reached:
      words=line.split()
      sfraction=eval(words[6])
      if abs(sfraction-0.54)>1.e-2:
        print("Wrong solid fraction")
        sys.exit(1)

if end_reached:
  sys.exit(0)
else:
  print("End time not reached")
  sys.exit(1)
