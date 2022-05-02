import numpy as np
import numpy
from scipy.special import perm
import sympy

def factorial_fraction(n, i):
    # calculate n! / (n-i)!
    res = 1
    for i in xrange(i):
        res *= n - i
    return res

def _M(i,j,k,l,n10,n01,n11,n00, n, f0):
    res = np.power(1 - 1./n/f0, n10 - i) * np.power(1 - 1./n/f0, n01 - j) * np.power(1 - 2./n/f0, n11 - k)
    res *= perm(n10, i) * perm(n01, j) * perm(n11, k) * perm(n00, l)
    # res *= factorial_fraction(n10, i) * factorial_fraction(n01, j) * factorial_fraction(n11, k) * factorial_fraction(n00, l)
    res /= np.power(n, i + j + k + l)
    # res /= perm(n, i+j+k+l)  # used in Ben's rare LD scripts
    return res


def M(i,j,k,l,n_obs,n,f0):
    # n_obs.shape = (# observations, 4)
    # n_obs = [n10, n01, n11, n00]  Make sure the order is correct!
    return _M(i,j,k,l,n_obs[:,0],n_obs[:,1],n_obs[:,2],n_obs[:,3],n,f0)

def D(n00, n10, n01, n11, n, f0):
    n_obs = np.hstack([n10,n01,n11,n00])
    return M(0,0,1,1,n_obs,n,f0) - M(0,0,1,1,n_obs,n,f0)

def _D2(n00, n10, n01, n11, n, f0):
    n_obs = np.hstack([n10,n01,n11,n00])
    return D2(n_obs, n, f0)

def D2(n_obs, n, f0):
    return M(0,0,2,2,n_obs,n,f0) - 2*M(1,1,1,1,n_obs,n,f0) + M(2,2,0,0,n_obs,n,f0)

# def LD_denom_manual(n00, n10, n01, n11, n, f0):
def LD_denom_manual(n_obs, n, f0):
    # manually typing all the terms, following appendix
    # n_obs = np.hstack([n00,n10,n01,n11])
    res = M(2,2,0,0,n_obs,n,f0)
    res += M(1,2,0,1,n_obs,n,f0)
    res += M(2,1,1,0,n_obs,n,f0)
    res += M(1,1,1,1,n_obs,n,f0)
    res += M(2,1,0,1,n_obs,n,f0)
    res += M(1,1,0,2,n_obs,n,f0)
    res += M(2,0,1,1,n_obs,n,f0)
    res += M(1,0,1,2,n_obs,n,f0)
    res += M(1,2,1,0,n_obs,n,f0)
    res += M(0,2,1,1,n_obs,n,f0)
    res += M(1,1,2,0,n_obs,n,f0)
    res += M(0,1,2,1,n_obs,n,f0)
    res += M(1,1,1,1,n_obs,n,f0)
    res += M(0,1,1,2,n_obs,n,f0)
    res += M(1,0,2,1,n_obs,n,f0)
    res += M(0,0,2,2,n_obs,n,f0)
    return res

def any_poly(mono_coeffs, n_obs, n, f0):
    # compute estimator for any polynomial of fAb, faB, fAB, fab
    # let each monomial be C * fAb**i * faB**j * fAB**k * fab**l
    # mono_coeffs will be a list of ((i,j,k,l), C)
    res = 0
    for powers, coeff in mono_coeffs:
        res += coeff * M(powers[0], powers[1], powers[2], powers[3], n_obs, n, f0)
    return res

def poly_to_mono_coeffs(p):
    return [(m, p.coeff_monomial(m)) for m in p.monoms()]

def generate_LE_poly():
    x, y, z, w = sympy.symbols('f_Ab,f_aB,f_AB,f_ab')
    exp = ((x+z)*(y+w)*(y+z)*(x+w))**2
    denom = sympy.Poly(exp, x, y, z, w)
    denom_monos = poly_to_mono_coeffs(denom)

    exp = x*y*z*w
    numer = sympy.Poly(exp, x, y, z, w)
    numer_monos = poly_to_mono_coeffs(numer)
    return numer_monos, denom_monos

def generate_LD_poly():
    x, y, z, w = sympy.symbols('f_Ab,f_aB,f_AB,f_ab')
    exp = (x+z)*(y+w)*(y+z)*(x+w)
    denom = sympy.Poly(exp, x, y, z, w)
    denom_monos = poly_to_mono_coeffs(denom)

    exp = (z*w - x*y)**2
    numer = sympy.Poly(exp, x, y, z, w)
    numer_monos = poly_to_mono_coeffs(numer)
    return numer_monos, denom_monos

def calculate_LE(n_obs, f0):
    numer_monos, denom_monos = generate_LE_poly()
    n_tots = np.sum(n_obs, axis=1)
    numer = any_poly(numer_monos, n_obs, n_tots, f0).mean()
    denom = any_poly(denom_monos, n_obs, n_tots, f0).mean()
    return numer / denom

def calculate_LD(n_obs, f0):
    numer_monos, denom_monos = generate_LD_poly()
    n_tots = np.sum(n_obs, axis=1)
    numer = any_poly(numer_monos, n_obs, n_tots, f0).mean()
    denom = any_poly(denom_monos, n_obs, n_tots, f0).mean()
    return numer / denom

def calculate_LD_Good2022(n_obs, n, f0):
    # ns = joint_passed_sites.sum(axis=2)

    # n11s = ((genotypes_1[None, :, :]) * (genotypes_2[:, None, :]) * joint_passed_sites).sum(axis=2)
    # n10s = (genotypes_1[None, :, :] * (1 - genotypes_2[:, None, :]) * joint_passed_sites).sum(axis=2)
    # n01s = ((1 - genotypes_1[None, :, :]) * (genotypes_2[:, None, :]) * joint_passed_sites).sum(axis=2)
    # n00s = ((1 - genotypes_1[None, :, :]) * (1 - genotypes_2[:, None, :]) * joint_passed_sites).sum(axis=2)
    n10s = n_obs[:, 0]
    n01s = n_obs[:, 1]
    n11s = n_obs[:, 2]
    n00s = n_obs[:, 3]
    ntots = n_obs.sum(axis=1)
    fstar = f0

    # Now do version based on counts
    # First calculate numerator
    rsquared_numerators = n11s * (n11s - 1) * n00s * (n00s - 1) * numpy.power(1 - 2.0 / ntots / fstar,
                                                                              n11s - 2) * numpy.power(
        1 - 1.0 / ntots / fstar, n10s - 0 + n01s - 0)
    rsquared_numerators += -2 * n10s * n01s * n11s * n00s * numpy.power(1 - 2.0 / ntots / fstar,
                                                                        n11s - 1) * numpy.power(1 - 1.0 / ntots / fstar,
                                                                                                n10s - 1 + n01s - 1)
    rsquared_numerators += n10s * (n10s - 1) * n01s * (n01s - 1) * numpy.power(1 - 2.0 / ntots / fstar,
                                                                               n11s - 0) * numpy.power(
        1 - 1.0 / ntots / fstar, n10s - 2 + n01s - 2)

    # Divide by sample size
    rsquared_numerators = rsquared_numerators * (ntots > 3.5) * 1.0 / (
                ntots * (ntots - 1) * (ntots - 2) * (ntots - 3) + 10 * (ntots < 3.5))

    # 1
    rsquared_denominators = n10s * (n10s - 1) * n01s * (n01s - 1) * numpy.power(1 - 2.0 / ntots / fstar,
                                                                                n11s - 0) * numpy.power(
        1 - 1.0 / ntots / fstar, n10s - 2 + n01s - 2)
    # 2
    rsquared_denominators += n10s * n01s * (n01s - 1) * n00s * numpy.power(1 - 2.0 / ntots / fstar,
                                                                           n11s - 0) * numpy.power(
        1 - 1.0 / ntots / fstar, n10s - 1 + n01s - 2)
    # 3
    rsquared_denominators += n10s * (n10s - 1) * n01s * n11s * numpy.power(1 - 2.0 / ntots / fstar,
                                                                           n11s - 1) * numpy.power(
        1 - 1.0 / ntots / fstar, n10s - 2 + n01s - 1)
    # 4
    rsquared_denominators += n10s * n01s * n11s * n00s * numpy.power(1 - 2.0 / ntots / fstar, n11s - 1) * numpy.power(
        1 - 1.0 / ntots / fstar, n10s - 1 + n01s - 1)
    # 5
    rsquared_denominators += n10s * (n10s - 1) * n01s * n00s * numpy.power(1 - 2.0 / ntots / fstar,
                                                                           n11s - 0) * numpy.power(
        1 - 1.0 / ntots / fstar, n10s - 2 + n01s - 1)
    # 6
    rsquared_denominators += n10s * n01s * n00s * (n00s - 1) * numpy.power(1 - 2.0 / ntots / fstar,
                                                                           n11s - 0) * numpy.power(
        1 - 1.0 / ntots / fstar, n10s - 1 + n01s - 1)
    # 7
    rsquared_denominators += n10s * (n10s - 1) * n11s * n00s * numpy.power(1 - 2.0 / ntots / fstar,
                                                                           n11s - 1) * numpy.power(
        1 - 1.0 / ntots / fstar, n10s - 2 + n01s - 0)
    # 8
    rsquared_denominators += n10s * n11s * n00s * (n00s - 1) * numpy.power(1 - 2.0 / ntots / fstar,
                                                                           n11s - 1) * numpy.power(
        1 - 1.0 / ntots / fstar, n10s - 1 + n01s - 0)
    # 9
    rsquared_denominators += n10s * n01s * (n01s - 1) * n11s * numpy.power(1 - 2.0 / ntots / fstar,
                                                                           n11s - 1) * numpy.power(
        1 - 1.0 / ntots / fstar, n10s - 1 + n01s - 2)
    # 10
    rsquared_denominators += n01s * (n01s - 1) * n11s * n00s * numpy.power(1 - 2.0 / ntots / fstar,
                                                                           n11s - 1) * numpy.power(
        1 - 1.0 / ntots / fstar, n10s - 0 + n01s - 2)
    # 11
    rsquared_denominators += n10s * n01s * n11s * (n11s - 1) * numpy.power(1 - 2.0 / ntots / fstar,
                                                                           n11s - 2) * numpy.power(
        1 - 1.0 / ntots / fstar, n10s - 1 + n01s - 1)
    # 12
    rsquared_denominators += n01s * n11s * (n11s - 1) * n00s * numpy.power(1 - 2.0 / ntots / fstar,
                                                                           n11s - 2) * numpy.power(
        1 - 1.0 / ntots / fstar, n10s - 0 + n01s - 1)
    # 13
    rsquared_denominators += n10s * n01s * n11s * n00s * numpy.power(1 - 2.0 / ntots / fstar, n11s - 1) * numpy.power(
        1 - 1.0 / ntots / fstar, n10s - 1 + n01s - 1)
    # 14
    rsquared_denominators += n01s * n11s * n00s * (n00s - 1) * numpy.power(1 - 2.0 / ntots / fstar,
                                                                           n11s - 1) * numpy.power(
        1 - 1.0 / ntots / fstar, n10s - 0 + n01s - 1)
    # 15
    rsquared_denominators += n10s * n11s * (n11s - 1) * n00s * numpy.power(1 - 2.0 / ntots / fstar,
                                                                           n11s - 2) * numpy.power(
        1 - 1.0 / ntots / fstar, n10s - 1 + n01s - 0)
    # 16
    rsquared_denominators += n11s * (n11s - 1) * n00s * (n00s - 1) * numpy.power(1 - 2.0 / ntots / fstar,
                                                                                 n11s - 2) * numpy.power(
        1 - 1.0 / ntots / fstar, n10s - 0 + n01s - 0)

    # divide by sample size
    rsquared_denominators = rsquared_denominators * (ntots > 3.5) * 1.0 / (
                ntots * (ntots - 1) * (ntots - 2) * (ntots - 3) + 10 * (ntots < 3.5))
    return rsquared_numerators, rsquared_denominators

