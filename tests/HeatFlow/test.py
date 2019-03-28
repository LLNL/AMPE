#!/usr/bin/env python
import sys
import subprocess

print("Test HeatFlux...")

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
exe = sys.argv[4]
inp = sys.argv[5]

#run AMPE
command = "{} {} {}".format(mpicmd,exe,inp)
output = subprocess.check_output(command,shell=True)

#analyse AMPE standard output
lines=output.split(b'\n')

end_reached = False
for line in lines:
  num_matches = line.count(b'cycle')
  if num_matches:
    print(line)
    words=line.split()
    time=eval(words[6])
  num_matches = line.count(b'Average')
  if num_matches:
    print(line)
    words=line.split()
    temperature=eval(words[3])
    if abs(10.*time-temperature)>1.e-4*time:
      print("Wrong temperature")
      print("time={}".format(time))
      print("temperature={}".format(temperature))
      sys.exit(1)

