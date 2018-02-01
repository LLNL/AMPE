import numpy as np
import matplotlib.pyplot as plt
import glob

filenames=sorted(glob.glob('concentration0_*.dat'))
colors=['b','g','r','c','m','y','k']
color=0
for fname in filenames:
    data=np.loadtxt(fname)
    X=data[:,0]
    Y=data[:,1]
    plt.plot(X,Y,colors[color])
    color=color+1
    if color>6:
      color=0

#save plot into png file    
plt.savefig('figure.png')

