
import sys, string
input_   =open(sys.argv[1],'r')

lines = input_.readlines()

shift=-960.*960.*0.25

time=-1.
for line in lines:
  if line.count('cycle'):
    words=line.split()
    time = words[6]
  if line.count('Total') and line.count('energy'):
    words=line.split()
    energy=eval(words[4])+shift
    print( "{}, {}".format(time, energy))

