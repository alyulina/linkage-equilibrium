import parameters
import config
import ld_theory

import gzip
import numpy as np
from math import log10
import os

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-p", "--parameters", type="string",
                  help="name pf the parameters regime from parameters.py",
                  dest="regime")
parser.add_option("--path", type="string",
                  help="root path to the data files",
                  dest="dat_path")

(options, args) = parser.parse_args()
regime = options.regime
dat_path = options.dat_path

params = parameters.params[regime]  # list of tuples of parameters
f0s = np.hstack([np.logspace(-4,-0.5,20), np.array([0.001, 0.003, 0.01, 0.03, 0.1])])
rs = np.array([params[idx][6] for idx in range(len(params))])
fB_stars = [0.05, 0.1, 0.2, 0.3, 0.5]

numers = np.empty((len(rs), len(fB_stars), len(f0s)))
denoms = np.empty((len(rs), len(fB_stars), len(f0s)))

for i in range(len(rs)):
    # LEs = {fstar: [] for fstar in fstars}
    # denominatorsquareds = {fstar: [] for fstar in fstars}
    # bare_numeratorsquareds = {fstar: [] for fstar in fstars}

    # gammas = []

    # N = params[type][param_idx][2]
    # r = params[type][param_idx][6]
    # gamma = 2 * N * r
    r = rs[i]

    #filename = os.path.join(dat_path, 'output_%s_%d.txt.gz' % (regime, i))
    filename = os.path.join(dat_path, 'fs_%s_%d.npy'%(regime,i)) # faster with npy
    if filename.endswith('npy'):
        all_fs = np.load(filename)
        f11s = all_fs[:, 0]
        f10s = all_fs[:, 1]
        f01s = all_fs[:, 2]
        f00s = 1-f11s-f10s-f01s
        fAs = f11s+f10s
        fBs = f11s+f01s
    else:
        print("Processing %s" % filename)
        file = gzip.GzipFile(filename, "r")
        f11s = []
        f10s = []
        f01s = []

        for line in file:
            if line.startswith('//'):
                continue
            items = line.split()

            f11 = float(items[0])
            f10 = float(items[1])
            f01 = float(items[2])

            f11s.append(f11)
            f10s.append(f10)
            f01s.append(f01)
        file.close()

        f11s = np.array(f11s)
        f10s = np.array(f10s)
        f01s = np.array(f01s)
        f00s = 1 - f11s - f10s - f01s
        fAs = f11s + f10s
        fBs = f11s + f01s

    print("Computing LE for %s" % filename)
    numerator = f11s * f00s * f10s * f01s
    sampling_variances = (fAs * (1 - fAs) * fBs * (1 - fBs)) ** 2

    for j, fB_star in enumerate(fB_stars):
        # pick a narror window for fB conditioning
        B_mask = (fBs < fB_star * 1.001) & (fBs > fB_star * 0.999)
        numer_masked = numerator[B_mask]
        var_masked = sampling_variances[B_mask]
        for k, f0 in enumerate(f0s):
            W_numer = np.exp(-f11s[B_mask] / f0 / fBs[B_mask] - f10s[B_mask] / f0 / (1-fBs[B_mask]))
            Ws = np.exp(-fAs[B_mask] / f0)
            LE_numerator = (numer_masked * Ws).mean()
            LE_denominator = (var_masked * Ws).mean()

            numers[i, j, k] = LE_numerator
            denoms[i, j, k] = LE_denominator

# saving all the intermediate files
savepath = os.path.join(config.CACHED_DIR, 'half_LE')
np.save(os.path.join(savepath, '%s_rs' % regime), rs)
np.save(os.path.join(savepath, '%s_fstars' % regime), f0s)
np.save(os.path.join(savepath, '%s_fB_stars' % regime), fB_stars)
np.save(os.path.join(savepath, '%s_numers' % regime), numers)
np.save(os.path.join(savepath, '%s_denoms' % regime), denoms)
