#!/usr/bin/env python
import sys
import subprocess
import os

print("Test Two grains sintering...")

mpicmd  = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe     = sys.argv[4]
inp     = sys.argv[5]
datadir = sys.argv[6]

#make symbolic link to calphad data
data = "2spheres.csv"
src = datadir+'/'+data
try:
  os.symlink(src, data)
except FileExistsError:
  os.remove(data)
  os.symlink(src, data)

#prepare initial conditions file
initfilename="2spheres.nc"
subprocess.call(["python3", "../../utils/make_multi_spheres.py",
  "--nx", "120", "--ny", "200", "--nz", "120",
  "--concentration-A", "1.,0.", "--concentration-B", "0.,1.", "--concentration-out", "0.,0.",
  "--spheres", data,
  initfilename])

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

os.remove(initfilename)
os.unlink(data)

#analyse AMPE standard output
lines=output.split(b'\n')
volfractions=[]

end_reached = False
end_time = 0.01
integral = -1.
integral0 = -1.
tol = 1.e-6
for line in lines:

  if line.count(b'cycle'):
    print(line)
    words=line.split()
    time=eval(words[6])
    if time>end_time:
      end_reached = True

  if line.count(b'Integral'):
    print(line)
    words=line.split()
    integral = eval(words[3])
    if integral0<0.:
      integral0 = integral
    else:
      if( abs(integral-integral0)>tol*integral0):
        print("Integral of composition is {}, expected {}".format(integral,integral0))
        sys.exit(1)

if end_reached:
  sys.exit(0)
else:
  print("End time not reached!")
  sys.exit(1)
