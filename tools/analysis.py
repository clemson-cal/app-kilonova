import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

a = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5])
b = np.array([6.74e9, 13.40e9, 20.07e9, 26.66e9, 33.53e9, 39.92e9, 46.52e9])

x = np.array([b[1], b[2], b[3], b[4], b[5], b[6]])
y = np.array([(b[1]-b[0])/0.5, (b[2]-b[1])/0.5, (b[3]-b[2])/0.5, (b[4]-b[3])/0.5, (b[5]-b[4])/0.5, (b[6]-b[5])/0.5])

popt, pcov = curve_fit(lambda fx,a,b: a*fx**b, x, y)
power_y = popt[0]*x**popt[1]
plt.plot(x, power_y, '--', label = "Fit equation: $y = %.6f \cdot x^{%.6f}$" %(popt[0], popt[1]))

plt.plot(x, y, 'o', label ='')
plt.xlabel(r'Distance (cm)')
plt.ylabel(r'Speed (cm/s)')
plt.title(r'Speed vs Distance')
plt.legend()
plt.show()
