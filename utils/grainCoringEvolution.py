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
#plot coring vs. time for multiple grains
#usage: python coringvstime.py ampe.log
#
import sys, string
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy

from math import pi, sqrt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc

rc('text', usetex=True)
#rc('font', family='serif')
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
rc('font', **font)

myfile=open(sys.argv[1],'r')
L=myfile.readlines()

search_string1 = 'Concentration deviation'
search_string2 = 'Simulation'
search_string3 = 'time'
search_string4 = 'cycle'
search_string5 = 'Volume of grain'

# loop over lines of file to collect data
times=[]
vols={}
cexs=[]
time='0'
time_found=1
maxt=0.
for line in L:
  num_matches1 = string.count(line, search_string1)
  num_matches2 = string.count(line, search_string2)
  num_matches3 = string.count(line, search_string3)
  num_matches4 = string.count(line, search_string4)
  num_matches5 = string.count(line, search_string5)
  
  if num_matches2 & num_matches3:
    w=string.split(line)
    time_found=1
  else:
    if num_matches4:
      w=string.split(line)
      time=w[6]
      time_found=1
  
  if num_matches5:
    w=string.split(line)
    
    volt=eval(w[5])
    gid=eval(w[3])
    vols[gid]=volt

  if num_matches1:
    if time_found:
      times.append(eval(time))
      maxt=eval(time)
      time_found=0
      
    w=string.split(line)
    
    gid   =eval(w[4])
    #coring=sqrt(eval(w[6])*vols[gid])
    coring=sqrt(eval(w[6]))
    #print gid
    if gid>=len(cexs):
      cexs.append([])
    cexs[gid].append(coring)

fig = plt.figure(1, figsize=(7.5,5.5))
axScatter = plt.subplot(111)
tlim=1.05*maxt
axScatter.set_xlim([0,tlim])

print 'Read',len(times),'time steps...'
print times

#number of time steps to skip between prints
inc=len(times)/50
print 'Skip',inc,'steps for prints...'

alltimes=[]
allcexs=[]
lastcex=[]
for cex in cexs:
  #print len(cex)
  for i in range(0,len(cex),inc):
    alltimes.append(times[i])
    allcexs.append(cex[i])
    if i+1>len(times)-inc:
      lastcex.append(cex[i])
    
axScatter.scatter(alltimes, allcexs, c=allcexs)
axScatter.set_xlabel('time (s)')
axScatter.set_ylabel('coring (mole fraction of Ni)')
#axScatter.set_xticks(np.arange(0.,5.,1.))
axScatter.set_xticks(np.arange(0.,32.,5.))

# create new axes on the right and on the top of the current axes
# The first argument of the new_vertical(new_horizontal) method is
# the height (width) of the axes to be created in inches.
divider = make_axes_locatable(axScatter)
axHisty = divider.append_axes("right", 1.2, pad=0.2, sharey=axScatter)

plt.setp(axHisty.get_yticklabels(), visible=False)

# setup bins:
xymax = np.max(np.fabs(lastcex))
binwidth = xymax/25
lim = ( int(xymax/binwidth) + 1) * binwidth

bins = np.arange(0, lim + binwidth, binwidth)
#print bins
N, bins, patches = axHisty.hist(lastcex, bins=bins, orientation='horizontal', color='r')

maxN=max(N)
maxN=maxN+1
if maxN%2>0:
  maxN=maxN+1

maxbin=bins[-1]
for bin, thispatch in zip(bins,patches):
  color = cm.jet(bin/maxbin)
  thispatch.set_facecolor(color)

# the yaxis of axHisty is shared with axScatter,
# thus there is no need to manually adjust the xlim and ylim of these
# axis.
for tl in axHisty.get_yticklabels():
    tl.set_visible(False)
axHisty.set_xticks([0, maxN/2, maxN])

plt.draw()
#plt.show()
#plt.ylim([0,1.05*xymax])
plt.ylim([0,0.085])

plt.savefig('coring.png')
