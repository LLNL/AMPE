
import sys, string
input_   =open(sys.argv[1],'r')

lines = input_.readlines()

count=0
for line in lines:
  count=count+1
  if count<3:
    continue
  words=line.split()
  print( "{}, {}".format(words[1], words[2]))

