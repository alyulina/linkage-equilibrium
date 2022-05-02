import numpy as np
import estimator

n = 1000
f0 = 3e2
# fs = [0.05, 0.05, 0.02, 1-0.02-0.1]
f = 0.2
fAb,faB,fAB = f, f, f
fab = 1 - fAb-faB-fAB
fs = [fAb, faB, fAB, fab]
fA = fAb+fAB
fB = faB+fAB

test_data = np.random.multinomial(n, fs, size=5000)
print(test_data.shape)

manual_LD_numer = estimator.D2(test_data, n, f0).mean()
manual_LD_denom = estimator.LD_denom_manual(test_data, n, f0).mean()

numer_monos, denom_monos = estimator.generate_LD_poly()
auto_LD_numer = estimator.any_poly(numer_monos, test_data, n, f0).mean()
auto_LD_denom = estimator.any_poly(denom_monos, test_data, n, f0).mean()

true_numer = (fAB*fab - fAb*faB)**2 * np.exp(-(fA+fB)/f0)
true_denom = fA * fB * (1-fA) * (1-fB) * np.exp(-(fA+fB)/f0)

numer, denom = estimator.calculate_LD_Good2022(test_data, n, f0)
numer = numer.mean()
denom = denom.mean()

print(manual_LD_numer, auto_LD_numer, numer, true_numer)
print(manual_LD_denom, auto_LD_denom, denom, true_denom)
print(manual_LD_numer/manual_LD_denom, auto_LD_numer/auto_LD_denom, numer/denom, true_numer/true_denom)
