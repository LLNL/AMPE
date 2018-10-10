import sys

beginval = eval(sys.argv[1])
endval   = eval(sys.argv[2])

end_time = 5000.
nslots   = 100
dt = end_time/nslots

print 'type = "value"'

slope = (endval-beginval)/end_time

for i in range(nslots+1):
  print 'time_'+str(i),'=', i*dt, ',', beginval+i*dt*slope

