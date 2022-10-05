import numpy as np
import estimator

n = 4000
f0 = 3e2
# fs = [0.05, 0.05, 0.02, 1-0.02-0.1]
f = 0.2
fAb,faB,fAB = f, f, f
fab = 1 - fAb-faB-fAB
fs = [fAb, faB, fAB, fab]
fA = fAb+fAB
fB = faB+fAB

test_data = np.random.multinomial(n, fs, size=1000000)
print(test_data.shape)

def calculate_manual(dat, n, f0):
    numer = estimator.D2(dat, n, f0).mean()
    denom = estimator.LD_denom_manual(dat, n, f0).mean()
    return numer, denom

def calculate_sympy(dat, n, f0):
    numer_monos, denom_monos = estimator.generate_LD_poly()
    numer = estimator.any_poly(numer_monos, dat, n, f0).mean()
    denom = estimator.any_poly(denom_monos, dat, n, f0).mean()
    return numer, denom


manual_LD_numer, manual_LD_denom = calculate_manual(test_data, n, f0)

auto_LD_numer, auto_LD_denom = calculate_sympy(test_data, n, f0)

true_numer = (fAB*fab - fAb*faB)**2 * np.exp(-(fA+fB)/f0)
true_denom = fA * fB * (1-fA) * (1-fB) * np.exp(-(fA+fB)/f0)

numer, denom = estimator.calculate_LD_Good2022_(test_data, f0)
numer = numer.mean()
denom = denom.mean()

print(manual_LD_numer, auto_LD_numer, numer, true_numer)
print(manual_LD_denom, auto_LD_denom, denom, true_denom)
print(manual_LD_numer/manual_LD_denom, auto_LD_numer/auto_LD_denom, numer/denom, true_numer/true_denom)
