import numpy as np
from sympy.physics.wigner import gaunt

# On the computation of the integrated products of three
# spherical harmonic

#
# Rules
#
# 1. Triangle relalation
# 2. Even rule
# 3. Sum of m must be 0
#

l1 = 0
l2 = 0
l3 = 0

m1 = 0
m2 = 0
m3 = 0

res = gaunt(l1, l2, l3, m1, m2, m3)

print(res)

print(res.n(64))
