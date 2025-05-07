#!/usr/bin/env python
import sys
import subprocess
import os

print("Test AlCu calphad multi-order...")

#prepare initial conditions file
initfilename="40x40x40.nc"
subprocess.call(["python3", "../../utils/make_nuclei.py",
  "--nx", "40", "--ny", "40", "--nz", "40", "-r", "15",
  "--center0", "0, 0, 0",
  "--concentration-in", "0.997", "--concentration-out", "0.983",
  initfilename])

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]
thermdatadir = sys.argv[6]

#make symbolic link to calphad database
calphad_data = "calphadAlCuLFccTheta.json"
src = thermdatadir+'/'+calphad_data
try:
  os.symlink(src, calphad_data)
except FileExistsError:
  os.remove(calphad_data)
  os.symlink(src, calphad_data)

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

#analyse AMPE standard output
lines=output.split(b'\n')

end_reached = False
end_time = 4.e-5
final_sfraction = 0.114
sfraction_checked = False

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
      if (dt-2.e-7)<0.:
        print("Wrong dt: too small")
        sys.exit(1)

  if line.count(b'fraction of phase 0'):
    print(line)
    if end_reached:
      sfraction_checked = True
      words=line.split()
      sfraction=eval(words[6])
      print("Final solid fraction: {}".format(sfraction))
      if abs(sfraction-final_sfraction)>1.e-3:
        print("Wrong solid fraction")
        sys.exit(1)

os.remove(initfilename)

if end_reached and sfraction_checked:
  sys.exit(0)
else:
  print("End time not reached...")
  sys.exit(1)
