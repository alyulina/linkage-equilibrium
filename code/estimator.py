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


def sterling_factorial(n):
    # return log(n!) using sterling's approximation
    # if n=0, return 0 (because 0! = 1)
    n = np.array(n)
    res = np.zeros(n.shape)
    res[n<=0] = 0
    n_mask = n[n>0]
    res[n>0] = n_mask * np.log(n_mask) - n_mask + np.log(2*np.pi*n_mask) / 2. + np.log(1 + 1. / 12 / n_mask)
    return res


def _M(i,j,k,l,n10,n01,n11,n00, n, f0):
    res = np.power(1 - 1./n/f0, n10 - i) * np.power(1 - 1./n/f0, n01 - j) * np.power(1 - 2./n/f0, n11 - k)
    res *= perm(n10, i) / np.power(n, i)
    res *= perm(n01, j) / np.power(n, j)
    res *= perm(n11, k) / np.power(n, k)
    res *= perm(n00, l) / np.power(n, l)
    # res *= factorial_fraction(n10, i) * factorial_fraction(n01, j) * factorial_fraction(n11, k) * factorial_fraction(n00, l)
    # res /= np.power(n, i + j + k + l)
    # res /= perm(n, i+j+k+l)  # used in Ben's rare LD scripts
    return res


def M(i,j,k,l,n_obs,n,f0):
    # n_obs.shape = (# observations, 4)
    # n_obs = [n10, n01, n11, n00]  Make sure the order is correct!
    return _M(i,j,k,l,n_obs[:,0],n_obs[:,1],n_obs[:,2],n_obs[:,3],n,f0)


def _M_half(i,j,k,l,n10,n01,n11,n00, n, f0):
    # The half rare version of M
    # assuming that n11+n01 is fixed to be nB, and all counts add to n
    # working with logs to avoid overflow

    nB = n01 + n11
    res = np.log(1 - 1./n/f0) * (n11 - k) + np.log(1 - 1./n/f0) * (n10 - i)
    res += np.log(perm(n10, i))  # this will give -np.inf if n10 < i; but that's okay
    res += np.log(perm(n01, j))
    res += np.log(perm(n11, k))
    res += np.log(perm(n00, l))
    res += sterling_factorial(nB - k - j) + sterling_factorial(n - nB - i - l) - sterling_factorial(n)
    # get rid of the small factor due to fB**nB
    fB_star = nB / n
    res -= np.log(fB_star) * nB + np.log(1-fB_star)*(n-nB)
    res = np.exp(res)

    # remove the entries where counts are smaller than i,j,k,l
    mask = (n10 < i) | (n01 < j) | (n11 < k) | (n00 < l)
    res[mask] = 0
    return res

def M_half(i,j,k,l,n_obs,n,f0):
    # n_obs.shape = (# observations, 4)
    # n_obs = [n10, n01, n11, n00]  Make sure the order is correct!
    return _M_half(i,j,k,l,n_obs[:,0],n_obs[:,1],n_obs[:,2],n_obs[:,3],n,f0)

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
    return [(m, int(p.coeff_monomial(m))) for m in p.monoms()]

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
    numer = any_poly(numer_monos, n_obs, n_tots, f0)
    denom = any_poly(denom_monos, n_obs, n_tots, f0)
    return numer, denom

def calculate_LD(n_obs, f0):
    numer_monos, denom_monos = generate_LD_poly()
    n_tots = np.sum(n_obs, axis=1)
    numer = any_poly(numer_monos, n_obs, n_tots, f0)
    denom = any_poly(denom_monos, n_obs, n_tots, f0)  # not yet averaged over pairs
    return numer, denom


def calculate_LE_half(n_obs, f0):
    n_tots = np.sum(n_obs, axis=1)
    nB = n_obs[:, 1] + n_obs[:, 2]
    nb = n_tots - nB
    fB = np.mean(nB / n_tots.astype(float))
    numer = M_half(1, 0, 1, 0, n_obs, n_tots, f0)
    denom = M_half(2, 0, 0, 0, n_obs, n_tots, f0) * nb / nB + 2 * M_half(1, 0, 1, 0, n_obs, n_tots, f0) + \
            M_half(0, 0, 2, 0, n_obs, n_tots, f0) * nB / nb
    denom = denom * fB * (1-fB)
    return numer, denom

def calculate_LE_numer_alternative(n_obs, f0, num_reps=1):
    """
    Alternative to the moment method, this function computes the numerator of LE using an estimator based on
    indicators of configurations of n_obs
    :param n_obs:
    :param f0:
    :param num_reps: how many subsampling reps to perform; if zero, use exact hypergeometric probability weights
    (but Sterling's approximation for factorial)
    :return: A single number for the numerator at this frequency scale
    """
    subsample_size = int(1 / f0)
    if subsample_size <= 1:
        return 0
    n_obs = n_obs.astype(int)
    ntots = n_obs.sum(axis=1)
    numers = []
    if num_reps == 0:
        n00s = n_obs[:, 3]
        factorials = sterling_factorial(n00s) - sterling_factorial(n00s + 3 - subsample_size) \
                    - sterling_factorial(ntots) + sterling_factorial(ntots - subsample_size)
        numers = n_obs[:, 0] * n_obs[:, 1] * n_obs[:, 2] * perm(subsample_size, subsample_size - 3) / (subsample_size**3)
        numers *= np.exp(factorials)
    for i in range(num_reps):
        sub_n10s = np.random.hypergeometric(n_obs[:, 0], ntots - n_obs[:, 0], subsample_size)
        sub_n01s = np.random.hypergeometric(n_obs[:, 1], ntots - n_obs[:, 1], subsample_size)
        sub_n11s = np.random.hypergeometric(n_obs[:, 2], ntots - n_obs[:, 2], subsample_size)
        passed = (sub_n11s==1) & (sub_n10s==1) & (sub_n01s==1)  # the indicator results
        numers.append(np.mean(passed) / (subsample_size**3))  # already averaging over site pairs
    return np.mean(numers)

def calculate_LE_denom_single_site_moment(nAs, ntots, f0):
    """
    Alternative way of calculating the denominator, using only single-site SFS
    Here, estimator is based on the moment method, and res should give the same mean as
    <fA^2 (1-fA)^2>

    TODO: might need to change nAs from the true SFS, instead of the one from site pairs
    :param nAs: single mutant counts at sites
    :param ntots: total number of samples at each site
    :param f0: frequency scale
    :return: a single number for the denominator of LE
    """
    res = np.power(1 - 1. / ntots / f0, nAs - 2) * perm(nAs, 2) / np.power(ntots, 2)
    res -= 2 * np.power(1 - 1. / ntots / f0, nAs - 3) * perm(nAs, 3) / np.power(ntots, 3)
    res += np.power(1 - 1. / ntots / f0, nAs - 4) * perm(nAs, 4) / np.power(ntots, 4)
    return np.mean(res) ** 2  # both A and B sites


def calculate_LE_denom_single_site_indicator(nAs, ntots, f0, num_reps=1):
    """
    Similar to calculate_LE_denom_single_site_moment, but using the indicator method
    """
    subsample_size = int(1 / f0)
    if subsample_size <= 1:
        return 0
    Mtmp = []
    nAs = nAs.astype(int)
    ntots = ntots.astype(int)
    for k in range(num_reps):
        subsamples = np.random.hypergeometric(nAs, ntots - nAs, subsample_size)
        M = np.mean(subsamples==2) * 2 / (subsample_size**2)
        M -= 2 * np.mean(subsamples==3) * 6 / (subsample_size**3)
        M += np.mean(subsamples==4) * 24 / (subsample_size**4)
        Mtmp.append(M)
    return np.mean(Mtmp)**2  # squaring for both sites in the denominator

# n10s = n_obs[:, 0]
# n01s = n_obs[:, 1]
# n11s = n_obs[:, 2]
# n00s = n_obs[:, 3]
def calculate_LD_Good2022(n_obs, f0):
    return _calculate_LD_Good2022(n_obs[:, 2], n_obs[:, 0], n_obs[:, 1], n_obs[:, 3], n_obs.sum(axis=1), f0)


def _calculate_LD_Good2022(n11s, n10s, n01s, n00s, ntots, fstar):
    """
    The manual version used in Good2022 paper
    :return:
    """
    # Now do version based on counts
    # First calculate numerator
    rsquared_numerators = n11s*(n11s-1)*n00s*(n00s-1)*numpy.power(1-2.0/ntots/fstar,n11s-2)*numpy.power(1-1.0/ntots/fstar,n10s-0+n01s-0)
    rsquared_numerators += -2*n10s*n01s*n11s*n00s*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-1)
    rsquared_numerators += n10s*(n10s-1)*n01s*(n01s-1)*numpy.power(1-2.0/ntots/fstar,n11s-0)*numpy.power(1-1.0/ntots/fstar,n10s-2+n01s-2)

    # Divide by sample size
    rsquared_numerators = rsquared_numerators*(ntots>3.5)*1.0/(ntots*(ntots-1)*(ntots-2)*(ntots-3)+10*(ntots<3.5))

    #1
    rsquared_denominators = n10s*(n10s-1)*n01s*(n01s-1)*numpy.power(1-2.0/ntots/fstar,n11s-0)*numpy.power(1-1.0/ntots/fstar,n10s-2+n01s-2)
    #2
    rsquared_denominators += n10s*n01s*(n01s-1)*n00s*numpy.power(1-2.0/ntots/fstar,n11s-0)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-2)
    #3
    rsquared_denominators += n10s*(n10s-1)*n01s*n11s*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-2+n01s-1)
    #4
    rsquared_denominators += n10s*n01s*n11s*n00s*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-1)
    #5
    rsquared_denominators += n10s*(n10s-1)*n01s*n00s*numpy.power(1-2.0/ntots/fstar,n11s-0)*numpy.power(1-1.0/ntots/fstar,n10s-2+n01s-1)
    #6
    rsquared_denominators += n10s*n01s*n00s*(n00s-1)*numpy.power(1-2.0/ntots/fstar,n11s-0)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-1)
    #7
    rsquared_denominators += n10s*(n10s-1)*n11s*n00s*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-2+n01s-0)
    #8
    rsquared_denominators += n10s*n11s*n00s*(n00s-1)*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-0)
    #9
    rsquared_denominators += n10s*n01s*(n01s-1)*n11s*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-2)
    #10
    rsquared_denominators += n01s*(n01s-1)*n11s*n00s*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-0+n01s-2)
    #11
    rsquared_denominators += n10s*n01s*n11s*(n11s-1)*numpy.power(1-2.0/ntots/fstar,n11s-2)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-1)
    #12
    rsquared_denominators += n01s*n11s*(n11s-1)*n00s*numpy.power(1-2.0/ntots/fstar,n11s-2)*numpy.power(1-1.0/ntots/fstar,n10s-0+n01s-1)
    #13
    rsquared_denominators += n10s*n01s*n11s*n00s*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-1)
    #14
    rsquared_denominators += n01s*n11s*n00s*(n00s-1)*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-0+n01s-1)
    #15
    rsquared_denominators += n10s*n11s*(n11s-1)*n00s*numpy.power(1-2.0/ntots/fstar,n11s-2)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-0)
    #16
    rsquared_denominators += n11s*(n11s-1)*n00s*(n00s-1)*numpy.power(1-2.0/ntots/fstar,n11s-2)*numpy.power(1-1.0/ntots/fstar,n10s-0+n01s-0)

    # divide by sample size
    rsquared_denominators = rsquared_denominators*(ntots>3.5)*1.0/(ntots*(ntots-1)*(ntots-2)*(ntots-3)+10*(ntots<3.5))
    return rsquared_numerators, rsquared_denominators

