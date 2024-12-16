import sys, string
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from math import pi, ceil
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc

rc('text', usetex=True)
#rc('font', family='serif')
#rc('xtick', labelsize=16) 
#rc('ytick', labelsize=16) 
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 20}
rc('font', **font)

myfile=open(sys.argv[1],'r')
lines=myfile.readlines()

# loop over lines of file 
times=[]
vols=[]
time='0'
time_found=1
maxt=0.
maxv=0.
for line in lines:
  num_matches2 = line.count('Simulation')
  num_matches3 = line.count('time')

  if num_matches2 & num_matches3:
    w=line.split()
    time_found=1
  else:
    if line.count('cycle'):
      w=line.split()
      if w[0]=='cycle':
        time=w[6]
        print(time)
        time_found=1

  if line.count('Volume fraction of phase'):
    if time_found:
      times.append(eval(time))
      maxt=eval(time)
      time_found=0
      
    w=line.split()
    
    volt=eval(w[6])
    maxv=max(maxv,volt)
    gid=eval(w[4])
    #print gid
    if gid>=len(vols):
      vols.append([])
    vols[gid].append(volt)
    #print time+' '+volt


fig = plt.figure(1, figsize=(8.,6.))
axScatter = plt.subplot(111)
tlim=1.05*maxt
axScatter.set_xlim([0,tlim])

print("#Read {} time steps...".format(len(times)))

#number of time steps to skip between prints
inc=int(len(times)/50)
inc=max(1,inc)

#extract data for plotting by sub-sampling and thersholding
alltimes=[]
allvols=[]
threshold=1.e-6
interval=100.
for vol in vols:
  oldtime=-10000.
  for i in range(0,len(vol)):
    if vol[i]>threshold:
      if times[i]>oldtime+interval:
        oldtime=times[i]
        alltimes.append(times[i])
        allvols.append(vol[i])

#build array with last measured volumes for bar chart
lastvol=[]
for vol in vols:
  v=vol[-1]
  if v>threshold:
    lastvol.append(vol[-1])

for vol in vols:
  print(' ')
  for i in range(0,len(vol),inc):
    print(times[i], vol[i])

#define colors
allcolors=[]
for vol in allvols:
  allcolors.append(vol/maxv)

axScatter.scatter(alltimes, allvols, c=allcolors)
axScatter.set_xlabel('time (s)')
axScatter.set_ylabel(r"grain size ($\mu m^3$)")

# create new axes on the right and on the top of the current axes
# The first argument of the new_vertical(new_horizontal) method is
# the height (width) of the axes to be created in inches.
divider = make_axes_locatable(axScatter)
axHisty = divider.append_axes("right", 1.2, pad=0.2, sharey=axScatter)

plt.setp(axHisty.get_yticklabels(), visible=False)

# setup bins:
xymax = np.max(np.fabs(lastvol))
binwidth = xymax/25
lim = ( int(xymax/binwidth) + 1) * binwidth

bins = np.arange(0, lim + binwidth, binwidth)
#print bins
N, bins, patches = axHisty.hist(lastvol, bins=bins, orientation='horizontal', color='r')

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

print(lastvol)

plt.draw()
#plt.show()
plt.ylim([threshold,1.05*xymax])

plt.savefig('volumes.png')
