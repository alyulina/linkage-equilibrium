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
fstars = np.hstack([np.logspace(-4,-0.5,20), np.array([0.001, 0.003, 0.01, 0.03, 0.1])])
rs = np.array([params[idx][6] for idx in range(len(params))])

numers = np.empty((len(rs), len(fstars)))
denoms = np.empty((len(rs), len(fstars)))
normalization = np.empty((len(rs), len(fstars)))
LDs = np.empty((len(rs), len(fstars)))

for i in range(len(rs)):
    # LEs = {fstar: [] for fstar in fstars}
    # denominatorsquareds = {fstar: [] for fstar in fstars}
    # bare_numeratorsquareds = {fstar: [] for fstar in fstars}

    # gammas = []

    # N = params[type][param_idx][2]
    # r = params[type][param_idx][6]
    # gamma = 2 * N * r
    r = rs[i]

    filename = os.path.join(dat_path, 'output_%s_%d.txt.gz' % (regime, i))
    #filename = os.path.join(dat_path, 'output_%s_%d_large.txt.gz' % (regime, i))
    #filename = os.path.join(dat_path, 'output_%s_%d_mega.txt.gz' % (regime, i))
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

    # ZHIRU: can consider adding multinomial sampling for a figure about finite sample

    numerator = f11s * f00s * f10s * f01s
    sampling_variances = (fAs * (1 - fAs) * fBs * (1 - fBs)) ** 2
    norm = (f10s * f01s) ** 2  # our approximation for the denom

    LD_numerator_uw = (f11s * f00s - f01s * f10s) ** 2
    LD_sampling_variances = fAs * (1 - fAs) * fBs * (1 - fBs)

    for j in range(len(fstars)):
        fstar = fstars[j]
        Hs = np.exp(-fAs / fstar - fBs / fstar)

        LE_numerator = (numerator * Hs).mean()
        LE_denominator = (sampling_variances * Hs).mean()
        LE_norm = (norm * Hs).mean()
        numers[i, j] = LE_numerator
        denoms[i, j] = LE_denominator
        normalization[i, j] = LE_norm

        # compute LD for comparison too
        LD_numerator = (LD_numerator_uw * Hs).mean()
        LD_denominator = (LD_sampling_variances * Hs).mean()
        LDs[i, j] = LD_numerator / LD_denominator

# saving all the intermediate files
savepath = os.path.join(config.CACHED_DIR, 'simulated')
np.save(os.path.join(savepath, '%s_rs' % regime), rs)
np.save(os.path.join(savepath, '%s_fstars' % regime), fstars)
np.save(os.path.join(savepath, '%s_numers' % regime), numers)
np.save(os.path.join(savepath, '%s_denoms' % regime), denoms)
np.save(os.path.join(savepath, '%s_norms' % regime), normalization)
np.save(os.path.join(savepath, '%s_LDs' % regime), LDs)
