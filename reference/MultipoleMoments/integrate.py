from numpy import genfromtxt
from scipy.integrate import simps

import numpy as np

data = genfromtxt('RhoTildaL', delimiter=None)

r = data[:990, 2]
a = data[:990, 5]
b = data[:990, 6]

c = a + b * 1.0j


l = 4
c = 4 * np.pi * r ** l * r ** 2 * c
c[-1] = 0.0
res = simps(c, r)
print(res)

r = data[:990, 2]
a = data[:990, 3]
b = data[:990, 4]

c = a + b * 1.0j


l = 0
c = 4 * np.pi * r ** l * r ** 2 * c
c[-1] = 0.0
res = simps(c, r)
print(res)