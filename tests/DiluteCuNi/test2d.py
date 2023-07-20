#!/usr/bin/env python
import sys
import subprocess
import os

print("Test Dilute_CuNi...")

#prepare initial conditions file
initfilename="1000.nc"
subprocess.call(["python3", "../../utils/make_bilayer.py",
  "--nx", "1000", "--ny", "1", "--nz", "1", "-r", "20",
  "-d", "1", "--centerx", "0.",
  initfilename])

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

input_file=open(inp,'r')
lines=input_file.readlines()

#determine domain size
lx=0.
for line in lines:
  if line.count('x_lo'):
    words=line.split()
    low=eval(words[2].split(',')[0])
  if line.count('x_up'):
    words=line.split()
    up=eval(words[2].split(',')[0])
    lx=(up-low)*1.e-6
    print("Lx = {}".format(lx))
    break

#analyse AMPE standard output
lines=output.split(b'\n')

vol=-1
time=-1.
target_velocity=0.06
tol=0.01
count=0
l=len(lines)  ## no of lines in file
print("Time         Velocity")
for line in range(l): ## loop over lines of file 
  if lines[line].count(b'Volume fraction'):
    time_old=time
    #search for time in previous lines
    for line2 in range(line-1,line-1000,-1):
      if lines[line2].count(b'cycle'):
        w=lines[line2].split()
        time=eval(w[6])
        #print time
        break
    words=lines[line].split()
    vol_old=vol
    vol=eval(words[6])
    if vol_old>0.:
      velocity=(vol-vol_old)*lx/(time-time_old)
      print("{}, {}".format(time,velocity))
      if(count>2):
        print("velocity = {}".format(velocity))
        print("target velocity = {}".format(target_velocity))
        if(abs(velocity-target_velocity)>tol):
          sys.exit(1)
      count=count+1

target_time = 7.5e-6
if time<target_time:
  print("Target time {} not reached!".format(target_time))
  sys.exit(1)

os.remove(initfilename)

sys.exit(0)
