import sys
import numpy
import os
import logging
import argparse
from datetime import datetime
from numpy.random import randint, shuffle, poisson, binomial, choice, hypergeometric

from estimator import calculate_LD_Good2022, calculate_LD, calculate_LE

parser = argparse.ArgumentParser()
parser.add_argument('--path', type=str, required=True, 
                    help='Path to the npy file with all pairs of sites')
parser.add_argument('stats', choices=['LD', 'LE', 'LD_manual'])
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

if args.savepath:
    savepath = os.path.join(args.savepath, args.stats, typename, accession)
else:
    savepath = os.path.join('./cached', args.stats, typename, accession)
os.makedirs(savepath, exist_ok=True)

if args.stats=='LD':
    stat_func = calculate_LD
elif args.stats=='LE':
    stat_func = calculate_LE
else:
    stat_func = calculate_LD_Good2022

if args.debug:
    print("Debug mode: exiting")
    quit()

"""
Loading the ns for a species
"""
ns = numpy.load(path).astype(float)
print("Done loading")
n11s = ns[:, 0]
n10s = ns[:, 1]
n01s = ns[:, 2]
n00s = ns[:, 3]
ells = ns[:, 4]
types = ns[:, 5]
ntots = n11s+n10s+n01s+n00s

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

# set genome wide pairs to 1e7 so that later code can work
mask = (ells > 1e4) & (ells < 1e8)
ells[mask] = 1e7

"""
Fix a few ells, scan a range of f0
"""
# f0s = numpy.hstack([numpy.logspace(-1.5,-0.5,20),[1e02]])
f0s = numpy.hstack([numpy.logspace(-3,-0.5,20),[1e02]])
# f0s = numpy.hstack([numpy.logspace(-3,-0.5,1),[1e02]])
capped_f0s = numpy.clip(f0s,0,1.5)

ellranges = [(1e04,1e08),(1e03,2e03),(1e02,3e02)]
colors = ['#1f77b4','#2ca02c','#ff7f0e']
labels =['Genome avg','1000<$\ell$<2000','100<$\ell$<300']

for idx in range(0,len(ellranges)):
    ellmin,ellmax = ellranges[idx]
    color = colors[idx]
    label=labels[idx]

    bins = numpy.arange(0,ntots.max())-0.5
    xs = numpy.arange(0,ntots.max()-1)

    good_idxs = (ells>ellmin)*(ells<ellmax)

    sigmasquareds = []
    for f0 in f0s:
        print(f0)
        # numer, denom = calculate_LD(n_obs[good_idxs, :],f0)
        # numer, denom = calculate_LD_Good2022(n11s[good_idxs],n10s[good_idxs],n01s[good_idxs],n00s[good_idxs],ntots[good_idxs],f0)
        numer, denom = stat_func(n_obs[good_idxs, :],f0)
        sigmasquared = numer.mean() / denom.mean()
        sigmasquareds.append(sigmasquared)

    sigmasquareds = numpy.array(sigmasquareds)
    numpy.save(os.path.join(savepath, "f_scan_f0s_{}".format(idx)), capped_f0s)
    numpy.save(os.path.join(savepath, "f_scan_stats_{}".format(idx)), sigmasquareds)

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
    numer, denom = stat_func(n_obs, fstar)

    sigmasquareds2 = []
    avg_ells2 = []
    for l in big_theory_ells:
        lmin = l*numpy.power(10,-0.1)
        lmax = l*numpy.power(10,0.1)
        good_idxs = (ells>=lmin)*(ells<=lmax)
        if good_idxs.sum() < 10:
            # fewer than 10 points!
            continue

        sigmasquared = (numer*(ells>=lmin)*(ells<=lmax)).mean() / (denom*(ells>=lmin)*(ells<=lmax)).mean()

        local_ells = ells[(ells>=lmin)*(ells<=lmax)]
        # Harmonic mean
        avg_ell = 1.0/(1.0/local_ells).mean()  # TODO: zhiru askes why this mean
        # regular mean
        # avg_ell = local_ells.mean()

        avg_ells2.append(avg_ell)
        sigmasquareds2.append(sigmasquared)

    avg_ells = numpy.array(avg_ells2)
    avg_ells[-1]=6e03
    sigmasquareds = numpy.array(sigmasquareds2)

    numpy.save(os.path.join(savepath, "ell_scan_ells_{}".format(idx)), avg_ells)  # the genome wide distance is set to be 6e03 here
    numpy.save(os.path.join(savepath, "ell_scan_stats_{}".format(idx)), sigmasquareds)

now = datetime.now()
dt_string = now.strftime("%H:%M:%S")
print("Finished ell scan at {}".format(dt_string))
