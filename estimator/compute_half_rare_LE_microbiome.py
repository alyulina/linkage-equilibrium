"""
Similar to compute_rare_LE_microbiome.py, but for Lambda_1 (the version that
conditions the frequency of the B allele).
"""

import sys
import numpy
import os
import logging
import argparse
from datetime import datetime

from estimator import calculate_LE_half

parser = argparse.ArgumentParser()
parser.add_argument('--path', type=str, required=True, 
                    help='Path to the npy file with all pairs of sites')
parser.add_argument('--savepath', type=str,
                    help='Path to where outputs should be saved')
parser.add_argument('--debug', action='store_true', 
                    help='will terminate after parsing the filename')
parser.add_argument('--pairtype', type=int, choices=range(3), help='the type of pair to process')
args = parser.parse_args()

path = args.path
accession = path.split('/')[-1].split('.')[0]

now = datetime.now()
dt_string = now.strftime("%H:%M:%S")
print("Started processing {} at {}".format(accession, dt_string))

if args.pairtype is None:
    typename = 'all'
else:
    typename = ['syn', 'syn_non_syn', 'non_syn'][args.pairtype]

stats = 'half_LE'
if args.savepath:
    savepath = os.path.join(args.savepath, stats, typename, accession)
else:
    savepath = os.path.join('./cached', stats, typename, accession)
print("saving to {}".format(savepath))
os.makedirs(savepath, exist_ok=True)

stat_func = calculate_LE_half

if args.debug:
    print("Debug mode: exiting")
    quit()

"""
Loading the ns for a species
"""
if path.endswith('npy'):
    ns = numpy.load(path).astype(float)
    print("Done loading")
    n11s = ns[:, 0]
    n10s = ns[:, 1]
    n01s = ns[:, 2]
    n00s = ns[:, 3]
    ells = ns[:, 4]
    types = ns[:, 5]
    ntots = n11s+n10s+n01s+n00s
elif path.endswith('.gz'):
    import gzip
    file = gzip.GzipFile(path, "r")
    file.readline()  # header
    n11s = []
    n10s = []
    n01s = []
    n00s = []
    ells = []
    print("Loading data...")
    for line in file:
        items = line.split()
        n11 = int(items[0])
        n10 = int(items[1])
        n01 = int(items[2])
        n00 = int(items[3])
        ell = int(items[4])

        n11s.append(n11)
        n10s.append(n10)
        n01s.append(n01)
        n00s.append(n00)
        ells.append(ell)
    file.close()
    n11s = numpy.array(n11s) * 1.0
    n10s = numpy.array(n10s) * 1.0
    n01s = numpy.array(n01s) * 1.0
    n00s = numpy.array(n00s) * 1.0
    ntots = n11s + n10s + n01s + n00s
    ells = numpy.array(ells)
    types = numpy.ones(ells.shape)
else:
    raise RuntimeError("site count file format is not recognized")

max_ntot = ntots.max()

# keeping only well covered sites
good_idxs = (ntots>0.95*max_ntot)
print("max_ntot", max_ntot)
print("Kept", good_idxs.sum()*1.0/len(good_idxs))


# filter pairs according to pair type
if args.pairtype is None:
    mask = numpy.ones(ntots.shape).astype(bool)
else:
    mask = (types==args.pairtype)
good_idxs = good_idxs & mask

n11s = n11s[good_idxs]
n10s = n10s[good_idxs]
n01s = n01s[good_idxs]
n00s = n00s[good_idxs]
n_obs = numpy.vstack([n10s, n01s, n11s, n00s]).T
print("Total data shape: {}, {} bytes".format(n_obs.shape, n_obs.nbytes))
ntots = ntots[good_idxs]
ells = ells[good_idxs]
nB = n11s + n01s

# set genome wide pairs to 1e7 so that later code can work
mask = (ells > 1e4) & (ells < 1e8)
ells[mask] = 1e7

# for conditioning the frequency of B allele
nB_stars = [1000, 1500, 2000]

"""
Fix a few ells, scan a range of f0 and nB_star
"""
f0s = numpy.hstack([numpy.logspace(-3,-0.5,20),[1e02]])
capped_f0s = numpy.clip(f0s,0,1.5)

ellranges = [(1e04,1e08),(1e03,2e03),(1e02,3e02)]
colors = ['#1f77b4','#2ca02c','#ff7f0e']
labels =['Genome avg','1000<$\ell$<2000','100<$\ell$<300']

for B_idx in range(len(nB_stars)):
    nB_star = nB_stars[B_idx]
    nB_mask = (nB > nB_star * 0.9) & (nB < nB_star * 1.1)
    for idx in range(0,len(ellranges)):
        ellmin,ellmax = ellranges[idx]
        color = colors[idx]
        label=labels[idx]

        bins = numpy.arange(0,ntots.max())-0.5
        xs = numpy.arange(0,ntots.max()-1)

        good_idxs = (ells>ellmin)&(ells<ellmax)&(nB_mask)

        LEs = []
        all_numers = []
        all_denoms = []
        for f0 in f0s:
            print(f0)
            numers, denoms = stat_func(n_obs[good_idxs, :],f0)
            if len(numers) == 0:
                numer = 0
                denom = 0
            else:
                numer = numers.mean()
                denom = denoms.mean()
            all_numers.append(numer)
            all_denoms.append(denom)

        all_numers = numpy.array(all_numers)
        all_denoms = numpy.array(all_denoms)
        LEs = all_numers / all_denoms
        numpy.save(os.path.join(savepath, "f_scan_stats_{}_{}".format(idx, B_idx)), LEs)
        numpy.save(os.path.join(savepath, "f_scan_stats_numer_{}_{}".format(idx, B_idx)), all_numers)
        numpy.save(os.path.join(savepath, "f_scan_stats_denom_{}_{}".format(idx, B_idx)), all_denoms)

numpy.save(os.path.join(savepath, "f_scan_f0s".format(idx)), capped_f0s)
now = datetime.now()
dt_string = now.strftime("%H:%M:%S")
print("Finished computing f0 scan at {}".format(dt_string))

"""
Fix a few f0, scan a range of distances
"""
raw_theory_ells = numpy.logspace(0,3,200)*3  # Specify the range of ells
theory_ells = numpy.hstack([raw_theory_ells,[6e03]])  # for plotting on x axis
big_theory_ells = numpy.hstack([raw_theory_ells,[1e07]])  # for filtering data; 1e7 is chosen to represent genome wide pairs

fstars = [1e02,1e-01,3e-02, 1e-2, 3e-3, 1e-3]  # Specify the f0s to calculate


for idx in range(0,len(fstars)):
    fstar = fstars[idx]
    capped_fstar = min([fstar,1.0])
    print("Processing fstar =", fstar, capped_fstar)

    # numer, denom = calculate_sigmasquared(n11s,n10s,n01s,n00s,ntots,fstar)
    # numer, denom = calculate_LD(n_obs,fstar)
    # numers, denoms = stat_func(n_obs, fstar)

    for B_idx in range(len(nB_stars)):
        nB_star = nB_stars[B_idx]
        nB_mask = (nB > nB_star * 0.9) & (nB < nB_star * 1.1)

        avg_ells2 = []
        all_numers = []
        all_denoms = []
        for l in big_theory_ells:
            lmin = l*numpy.power(10,-0.1)
            lmax = l*numpy.power(10,0.1)
            good_idxs = (ells>=lmin)&(ells<=lmax)&(nB_mask)
            if good_idxs.sum() < 10:
                # fewer than 10 points!
                continue
            
            numers, denoms = stat_func(n_obs[good_idxs, :],fstar)
            curr_numer = numers.mean()
            curr_denom = denoms.mean()

            local_ells = ells[good_idxs]
            # Harmonic mean
            avg_ell = 1.0/(1.0/local_ells).mean()  # TODO: zhiru askes why this mean
            # regular mean
            # avg_ell = local_ells.mean()

            avg_ells2.append(avg_ell)
            all_numers.append(curr_numer)
            all_denoms.append(curr_denom)

        avg_ells = numpy.array(avg_ells2)
        avg_ells[-1]=6e03
        all_numers = numpy.array(all_numers)
        all_denoms = numpy.array(all_denoms)
        LEs = all_numers / all_denoms

        numpy.save(os.path.join(savepath, "ell_scan_ells_{}_{}".format(idx, B_idx)), avg_ells)  # the genome wide distance is set to be 6e03 here
        numpy.save(os.path.join(savepath, "ell_scan_stats_{}_{}".format(idx, B_idx)), LEs)
        numpy.save(os.path.join(savepath, "ell_scan_stats_numer_{}_{}".format(idx, B_idx)), all_numers)
        numpy.save(os.path.join(savepath, "ell_scan_stats_denom_{}_{}".format(idx, B_idx)), all_denoms)

now = datetime.now()
dt_string = now.strftime("%H:%M:%S")
print("Finished ell scan at {}".format(dt_string))
