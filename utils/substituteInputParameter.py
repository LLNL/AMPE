#substitute parameter in input file
#usage:
#   python substituteInputParameter.py input_file \
#      parameter_name parameter_value > new_input_file
import sys, string

lines = open(sys.argv[1],'r')

parameter_name = sys.argv[2]
parameter_val  = sys.argv[3]
newline = " "

# loop over lines of file
for line in lines:
  words=line.split()
  if len(words)>0:
    if words[0]==parameter_name:
      words[2]=parameter_val
      newline = newline.join(words)
      print("  {}".format(newline))
    else:
      print(line, end = '')
  else:
    print(line, end = '')

