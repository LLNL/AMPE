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
from math import floor
myfile=open(sys.argv[1],'r')
nreduce=eval(sys.argv[2])
L=myfile.readlines()

search_string1 = 'Volume of grain'
search_string2 = 'Simulation'
search_string3 = 'time'
search_string4 = 'cycle'

vols={}
vol=[]
time=0.
# loop over lines of file 
for line in L:
  num_matches1 = string.count(line, search_string1)
  num_matches2 = string.count(line, search_string2)
  num_matches3 = string.count(line, search_string3)
  num_matches4 = string.count(line, search_string4)
  
  if num_matches2 & num_matches3:
    w=string.split(line)
    time_found=1
  else:
    if num_matches4:
      w=string.split(line)
      time_found=1
      time=w[6]
      #print time

  if num_matches1:
    w=string.split(line)
    
    volume=w[5]
    gid=eval(w[3])
    #print gid
    if gid>=len(vol):
      vol.append([])
    vol[gid]=volume
  else:
    if len(vol)>0:
      vols[time]=vol
      vol=[] 

maxV=0.
for key in vols.keys():
  vol=vols[key]
  for vv in vol:
    if eval(vv)>maxV:
      maxV=eval(vv)
print '#max. V=',maxV

nintervals=20
h=1.01*maxV/nintervals

count=0
for key in sorted(vols.keys(), key=float):
  if( count % nreduce )==0:
    vol=vols[key]
    print '#time=',key
    print '#Number of grains=',len(vol)

    dist=[]
    for i in range(nintervals):
      dist.append(0.)

    for vv in vol:
      v=int(floor(eval(vv)/h))
      dist[v]=dist[v]+1

    for i in range(len(dist)):
      print (i+0.5)*h, dist[i]/len(vol)
    print '\n'
  count=count+1
