#!/usr/bin/env python
import sys
import subprocess

print("Test TwoBilayers...")

#prepare initial conditions file
subprocess.call(["python3", "../../utils/make_bilayer.py", "-2",
  "-x", "256", "-y", "8", "-z", "8", "-r", "50", "-d", "1",
  "-c", "0.25", "--concentration-in", "0.2",
  "--quat-in", "1,0,0,0", "--quat-in-two", "0.95,0.31225,0.,0.",
  "--quat-out", "0,0,1,0", "--qlen", "4",
  "initial.nc"])

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
    #check for end time reached
    if eval(words[6])>0.05:
      end_reached = True

  num_matches = line.count(b'fraction')
  if num_matches:
    print(line)
    if end_reached:
      words=line.split()
      if abs(eval(words[6])-0.91)>1.e-2:
        print("Wrong solid fraction")
        sys.exit(1)

  num_matches = line.count(b'qint')
  if num_matches:
    if end_reached:
      words=line.split()
      qe = eval(words[4])
      print("Q energy = {}".format(qe))
      if abs(qe-0.000764)>1.e-5:
        print("Wrong qint energy")
        sys.exit(1)

  num_matches = line.count(b'orient')
  if num_matches:
    if end_reached:
      words=line.split()
      oe = eval(words[4])
      print("Orient enegy = {}".format(oe))
      if abs(oe-0.00025)>1.e-5:
        print("Wrong orient energy")
        sys.exit(1)

if end_reached:
  sys.exit(0)
else:
  sys.exit(1)

