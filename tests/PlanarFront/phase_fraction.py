#usage:
#  python phase_fraction.py ampe.log > volfraction.dat
#  
import sys, string
from math import pi
myfile=open(sys.argv[1],'r')
lines=myfile.readlines()

print '#phase fraction vs. time [s]'
for line in lines: ## loop over lines of file 
  num_matches = string.count(line, 'cycle')

  if num_matches :
    w=string.split(line)
    time=eval(w[6])
    #print time

  num_matches = string.count(line, 'fraction')

  if num_matches :
    w=string.split(line)
    fraction=eval(w[6])
    print time, fraction

