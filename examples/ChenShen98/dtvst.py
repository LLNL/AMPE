# Copyright (c) 2018, Lawrence Livermore National Security, LLC and
# UT-Battelle, LLC.
# Produced at the Lawrence Livermore National Laboratory and
# the Oak Ridge National Laboratory
# Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
# LLNL-CODE-747500
# All rights reserved.
# This file is part of AMPE. 
# For details, see https://github.com/LLNL/AMPE
# Please also read AMPE/LICENSE.
# 
#usage:
#  python dtvst.py ampe.log > volfraction.dat
#  
import sys, string
from math import pi
myfile=open(sys.argv[1],'r')
lines=myfile.readlines()

search_string = 'cycle'
print '#dt fraction vs. time [s]'
for line in lines: ## loop over lines of file 
  num_matches = string.count(line, search_string)

  if num_matches :
    w=string.split(line)
    time=eval(w[6])
    #print time
    dt=eval(w[10])
    print time,dt
