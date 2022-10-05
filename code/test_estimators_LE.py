import numpy as np
import estimator

n = 4000
# f0 = 3e2
f0 = 0.1
# fs = [0.05, 0.05, 0.02, 1-0.02-0.1]
f = 0.2
fAb,faB,fAB = f, f, f
fab = 1 - fAb-faB-fAB
fs = [fAb, faB, fAB, fab]
fA = fAb+fAB
fB = faB+fAB

test_data = np.random.multinomial(n, fs, size=1000000)
print(test_data.shape)


def calculate_sympy(dat, n, f0):
    numer_monos, denom_monos = estimator.generate_LE_poly()
    denom = estimator.any_poly(denom_monos, dat, n, f0).mean()
    numer = estimator.any_poly(numer_monos, dat, n, f0).mean()
    numer_monos, denom_monos = estimator.generate_LD_poly()
    LD_numer = estimator.any_poly(numer_monos, dat, n, f0).mean()
    LD_denom = estimator.any_poly(denom_monos, dat, n, f0).mean()
    return numer, denom, LD_numer, LD_denom


auto_LE_numer, auto_LE_denom, auto_LD_numer, auto_LD_denom = calculate_sympy(test_data, n, f0)

true_numer = (fAB*fab*fAb*faB) * np.exp(-(fA+fB)/f0)
true_denom = (fA * fB * (1-fA) * (1-fB))**2 * np.exp(-(fA+fB)/f0)

true_LD_denom = (fA * fB * (1-fA) * (1-fB)) * np.exp(-(fA+fB)/f0)
true_LD_numer = (fAB*fab - fAb*faB)**2 * np.exp(-(fA+fB)/f0)

print(auto_LE_numer, true_numer)
print(auto_LE_denom, true_denom)
# print(auto_LD_numer, true_LD_numer)
# print(auto_LD_denom, true_LD_denom)
print(auto_LE_numer/auto_LE_denom, true_numer/true_denom)
