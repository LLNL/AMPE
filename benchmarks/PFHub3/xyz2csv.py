import sys, string

input_   =open(sys.argv[1],'r')

#function to order data points
def comparedata(x,y):
  if y>x:
    return -y
  else:
    return x

lines = input_.readlines()

#build list with data points
data=[]
for line in lines[2:]:
  words=line.split()
  data.append((eval(words[1]),eval(words[2])))

#sort data
data.sort(key = lambda x: comparedata(x[0],x[1]))

#print data in csv format
for d in data:
  print("{}, {}".format(d[0], d[1]))
