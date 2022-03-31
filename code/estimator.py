import numpy as np

def factorial_fraction(n, i):
    # calculate n! / (n-i)!
    res = 1
    for i in xrange(i):
        res *= n - i
    return res

def M(i,j,k,l,n00,n10,n01,n11, n, f0):
    res = np.powers(1 - 1./n/f0, n10 - i) * np.powers(1 - 1./n/f0, n01 - j) * np.powers(1 - 2./n/f0, n11 - k)
    res *= factorial_fraction(n10, i) * factorial_fraction(n01, j) * factorial_fraction(n11, k) * factorial_fraction(n00, l)
    res /= np.powers(n, i + j + k + l)
    return res
