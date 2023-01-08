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
rs = np.array([params[idx][6] for idx in range(len(params))])
ns = np.array([30, 100, 300, 1000, 3000, 10000])
FGs = []

for i in range(len(rs)):
    r = rs[i]

    filename = os.path.join(dat_path, 'output_%s_%d.txt.gz' % (regime, i))
    print("Processing %s" % filename)
    file = gzip.GzipFile(filename, "r")
    n11s = [[] for n in ns]
    n10s = [[] for n in ns]
    n01s = [[] for n in ns]
    n00s = []

    for line in file:
        if line.startswith('//'):
            continue
        items = line.split()

        f11 = float(items[0])
        f10 = float(items[1])
        f01 = float(items[2])

        for idx, n in enumerate(ns):
            nsample = np.random.binomial(n, [f11, f10, f01, 0])
            n11s[idx].append(nsample[0])
            n10s[idx].append(nsample[1])
            n01s[idx].append(nsample[2])

    file.close()
    for idx, n in enumerate(ns):
        # multinomial sampling because four gamete is discrete
        n11s[idx] = np.array(n11s[idx])
        n10s[idx] = np.array(n10s[idx])
        n01s[idx] = np.array(n01s[idx])
        ntots = np.ones_like(n11s[idx])*n
        n00s.append(ntots-n10s[idx]-n11s[idx]-n01s[idx])

    four_gametes = [np.sum((n11s[i]>0)*(n10s[i]>0)*(n01s[i]>0)*(n00s[i]>0)) / float(len(n11s[i])) for i in range(len(ns))]
    FGs.append(four_gametes)

FGs = np.array(FGs)
savepath = os.path.join(config.CACHED_DIR, 'four_gamete')
np.save(os.path.join(savepath, '%s_FGs' % regime), FGs)
np.save(os.path.join(savepath, '%s_rs' % regime), rs)
np.save(os.path.join(savepath, '%s_ns' % regime), ns)
