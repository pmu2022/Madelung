import numpy as np
from scipy.integrate import quad, tplquad
from scipy.special import sph_harm

# regular integration


def density(r):
    return r**2 * np.exp(-r) * 4 * np.pi


q_00, _ = quad(density, 0, 1, limit=150, epsabs=1e-14)


print("q_00")
print(q_00)

# theta: azimuthal angle
# phi: polar angle

for l in range(0,3):

    m = 0

    prefactor = np.sqrt(4 * np.pi) / (2 * l + 1)


    def multi(phi, theta, r, l, m):
        ylm = np.real(sph_harm(m, l, theta, phi))
        return r**l * np.exp(-r) * ylm * r**2 * np.sin(theta)


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


    res = real_ * prefactor

    print(res)


    def multi(phi, theta, r, l, m):
        ylm = np.imag(sph_harm(m, l, theta, phi))
        return r**l * np.exp(-r) * ylm * r**2 * np.sin(theta)


    res, _ = tplquad(
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

    res = res * prefactor

    print(res)
