#Fit a quadratic polynomial through the data points in file 'zeta.csv'
#Printout polynomial coefficients, as well as plot in 'fit.png'
#Reference xref may need to be adjusted
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#
def poly2(x, A, B, C):
  y = A*x*x+B*x+C
  return y

ra = np.genfromtxt('zeta.csv', delimiter=',',dtype=None, names=True)

aa = ra.view(np.float64).reshape(len(ra), -1)

xref=800.
x=aa[:,0]-xref
print(x)
y=aa[:,1]

parameters, covariance = curve_fit(poly2, x,y)
fit_A=parameters[0]
fit_B=parameters[1]
fit_C=parameters[2]
print("polynomial {}*x^2 + {}*x + {}".format(fit_A,fit_B,fit_C))

fit_y = poly2(x, fit_A, fit_B, fit_C)
plt.plot(x,y,'o')
plt.plot(x,fit_y,'-')

plt.savefig('fit.png', dpi=100)

