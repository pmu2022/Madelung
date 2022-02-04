import numpy as np
from scipy.integrate import quad, tplquad
from scipy.special import sph_harm
from sympy.physics.wigner import gaunt

#
# G(r,E) = \sum_{L_p,L_pp} Y_{L_p} G_{L_p,L_pp}(r,E) Y_{L_pp}
#
# Easiest would be to have channels
#




# theta: azimuthal angle
# phi: polar angle


l = 2
m = 0

l_p = 2
m_p = 2

l_pp = 2
m_pp = -2

#
# Numerical solution
#

prefactor = np.sqrt(4 * np.pi) / (2 * l + 1)


def multi(phi, theta, r, l, m):
    ylm_p = sph_harm(m_p, l_p, theta, phi)
    ylm_pp = np.conj(sph_harm(m_pp, l_pp, theta, phi))
    ylm = sph_harm(m, l, theta, phi)
    return np.real(r**l * np.exp(-r) * ylm * r**2 * np.sin(theta) * ylm_p * ylm_pp)


real_qlm, _ = tplquad(
    multi,
    0,
    1,
    lambda theta: 0,
    lambda theta: np.pi,
    lambda theta, phi: 0,
    lambda theta, phi: 2 * np.pi,
    args=(l, m),
    epsabs=1e-8,
)


def multi(phi, theta, r, l, m):
    ylm_p = sph_harm(m_p, l_p, theta, phi)
    ylm_pp = np.conj(sph_harm(m_pp, l_pp, theta, phi))
    ylm = sph_harm(m, l, theta, phi)
    return np.imag(r**l * np.exp(-r) * ylm * r**2 * np.sin(theta) * ylm_p * ylm_pp)


imag_qlm, _ = tplquad(
    multi,
    0,
    1,
    lambda theta: 0,
    lambda theta: np.pi,
    lambda theta, phi: 0,
    lambda theta, phi: 2 * np.pi,
    args=(l, m),
    epsabs=1e-8,
)

res = (real_qlm + 1j * imag_qlm) * prefactor

print(res)
