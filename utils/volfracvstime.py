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
#usage:
#  python volfracvstime.py ampe.log > solid_fraction.csv
#  
import sys, string
from math import pi
myfile=open(sys.argv[1],'r')
L=myfile.readlines()
l=len(L)  ## no of lines in file

search_string1 = 'Volume fraction'
search_string2 = 'cycle'
time=0.
print('#Volume fraction vs. time [s]')
for line in range(l): ## loop over lines of file 
  num_matches1 = L[line].count(search_string1)

  if num_matches1 :
    #search for time in previous lines
    for line2 in range(line-1,line-1000,-1):
      num_matches2 = L[line2].count(search_string2)
      if num_matches2:
        w=L[line2].split()
        time=eval(w[6])
        #print time
        break
    words=L[line].split()
    vol=eval(words[6])
    print("{}, {}".format(time,vol))
