import numpy as np
from scipy.integrate import quad, tplquad
from scipy.special import sph_harm
from sympy.physics.wigner import gaunt

# theta: azimuthal angle
# phi: polar angle

# Higher order moments

l_p = 0
m_p = 0

l = 2
m = 0

# Prefactors are the same for both cases

prefactor = np.sqrt(4 * np.pi) / (2 * l + 1)

#
# Direct solution with Gaunt
#


def density(r):
    return r**l * np.exp(-r) * r**2

real_direct_qlm, _ = quad(density, 0, 1)

gaunt_number = np.float64(gaunt(l, l_p, l_p, m, m_p, m_p))

real_direct_qlm = real_direct_qlm * prefactor * gaunt_number

print(real_direct_qlm)

#
# Numerical solution
#


def multi(phi, theta, r, l, m):
    ylm_p = sph_harm(m_p, l_p, theta, phi)
    ylm = sph_harm(m, l, theta, phi)
    return np.real(
        r**l * np.exp(-r) * ylm * r**2 * np.sin(theta) * ylm_p * np.conj(ylm_p)
    )


real_qlm, _ = tplquad(
    multi,
    0,
    1,
    lambda theta: 0,
    lambda theta: np.pi,
    lambda theta, phi: 0,
    lambda theta, phi: 2 * np.pi,
    args=(l, m),
    epsabs=1e-12,
)


def multi(phi, theta, r, l, m):
    ylm_p = sph_harm(m_p, l_p, theta, phi)
    ylm = sph_harm(m, l, theta, phi)
    return np.imag(
        r**l * np.exp(-r) * ylm * r**2 * np.sin(theta) * ylm_p * np.conj(ylm_p)
    )


imag_qlm, _ = tplquad(
    multi,
    0,
    1,
    lambda theta: 0,
    lambda theta: np.pi,
    lambda theta, phi: 0,
    lambda theta, phi: 2 * np.pi,
    args=(l, m),
    epsabs=1e-12,
)

res = (real_qlm + 1j * imag_qlm) * prefactor

print(res)
