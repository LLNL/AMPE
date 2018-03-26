# Copyright (c) 2018, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory
# Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
# LLNL-CODE-747500
# All rights reserved.
# This file is part of AMPE. 
# For details, see https://github.com/LLNL/AMPE
# Please also read AMPE/LICENSE.
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# - Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the disclaimer below.
# - Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the disclaimer (as noted below) in the
#   documentation and/or other materials provided with the distribution.
# - Neither the name of the LLNS/LLNL nor the names of its contributors may be
#   used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
# LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
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
