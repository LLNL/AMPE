import sys

beginval = eval(sys.argv[1])
endval   = eval(sys.argv[2])

end_time = 5000.
nslots   = 10
dt = end_time/nslots

print ('type = "slope"')

slope = (endval-beginval)/end_time

for i in range(nslots+1):
  #print ('time_'+str(i),'=', i*dt, ',', beginval+i*dt*slope)
  itime = 'time_'+str(i)
  time = round(i*dt, 2)
  val = beginval+i*dt*slope
  print( itime+" = {}, {:.5g}".format(time, val) )
