import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#
def poly2(x, A, B, C):
  y = A*x*x+B*x+C
  return y

ra = np.genfromtxt('CvsTliquid.csv', delimiter=',',dtype=None, names=True)
rb = np.genfromtxt('CvsTsolid.csv', delimiter=',',dtype=None, names=True)

aa = ra.view(np.float64).reshape(len(ra), -1)
ab = rb.view(np.float64).reshape(len(rb), -1)

xref=900.
x=aa[:,0]-xref
print(x)
y=aa[:,1]

parameters, covariance = curve_fit(poly2, x,y)
fit_A=parameters[0]
fit_B=parameters[1]
fit_C=parameters[2]
print("polynomial {}*x^2 + {}*x + {}".format(fit_A,fit_B,fit_C))

fit_y = poly2(x, fit_A, fit_B, fit_C)
plt.plot(x+xref,y,'o')
plt.plot(x+xref,fit_y,'-')

x=ab[:,0]-xref
print(x)
y=ab[:,1]

parameters, covariance = curve_fit(poly2, x,y)
fit_A=parameters[0]
fit_B=parameters[1]
fit_C=parameters[2]
print("polynomial {}*x^2 + {}*x + {}".format(fit_A,fit_B,fit_C))

fit_y = poly2(x, fit_A, fit_B, fit_C)
plt.plot(x+xref,y,'o')
plt.plot(x+xref,fit_y,'-')

plt.savefig('fit_ceq.png', dpi=100)

