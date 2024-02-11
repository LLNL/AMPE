# Copyright (c) 2018, Lawrence Livermore National Security, LLC and
# UT-Battelle, LLC.
# Produced at the Lawrence Livermore National Laboratory and
# the Oak Ridge National Laboratory
# LLNL-CODE-747500
# All rights reserved.
# This file is part of AMPE. 
# For details, see https://github.com/LLNL/AMPE
# Please also read AMPE/LICENSE.
# 
import sys, string
from math import pi
import matplotlib.pyplot as plt

myfile=open(sys.argv[1],'r')
lines=myfile.readlines()
nlines=len(lines)  ## no of lines in file

time=0.
cell_length=12.8
cell_volume=cell_length*cell_length*cell_length

diameters=[]
times=[]

for line in range(nlines): ## loop over lines of file 
  num_matches1 = string.count(lines[line], 'Volume')

  if num_matches1:
    for line2 in range(line,line-25,-1):
      num_matches2 = string.count(lines[line2], 'cycle')
      if num_matches2:
        w=string.split(lines[line2])
        time=eval(w[6])
        break
    w=string.split(lines[line])
    vol=eval(w[6])*cell_volume
    d=2.*(3.*vol/(4.*pi))**(1./3.) # diameter for 3d sphere
    print time,d
    times.append(time)
    diameters.append(d)

# plot results
plt.plot(times,diameters,'r.')
plt.ylabel('Diameter (um)')
plt.xlabel('time (s)')
plt.axis([0.,1.,0.,8.])

plt.show()
#plt.savefig('growth.png', dpi=100)
