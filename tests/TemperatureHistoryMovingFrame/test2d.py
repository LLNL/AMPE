#!/usr/bin/env python
import sys
import subprocess
import os

print("Test moving frame with gradient given by T history...")

#prepare initial conditions file
initfilename="500.nc"
subprocess.call(["python3", "../../utils/make_bilayer.py",
  "--nx", "500", "--ny", "1", "--nz", "1", "-r", "250",
  "-d", "1", "--centerx", "0.",
  initfilename])

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]
T_data_dir = sys.argv[6]

#make symbolic link to T data
T_data = "temperature.csv"
src = T_data_dir+'/'+T_data
try:
  os.symlink(src, T_data)
except FileExistsError:
  os.remove(T_data)
  os.symlink(src, T_data)

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

os.remove(initfilename)
os.unlink(T_data)

input_file=open(inp,'r')
lines=input_file.readlines()

#analyse AMPE standard output
lines=output.split(b'\n')

time=-1.
fs=-1.
for line in lines: ## loop over lines of file 
  if line.count(b'cycle'):
    w=line.split()
    time=eval(w[6])
    print("Time: {}".format(time))
  if line.count(b'Volume fraction'):
    w=line.split()
    fs=eval(w[6])
    print("Volume fraction: {}".format(fs))

target_fs = 0.50
if abs(target_fs-fs)>0.005:
  print("Target fs not reached")
  sys.exit(1)

target_time = 0.022
if time<target_time:
  print("Target time not reached")
  sys.exit(1)

sys.exit(0)
