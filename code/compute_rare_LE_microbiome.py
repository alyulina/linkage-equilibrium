import sys
import numpy
import os
import gzip
from scipy.special import gammaln

import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
import matplotlib.pyplot as plt
from datetime import datetime

import matplotlib.gridspec as gridspec
from numpy.random import randint, shuffle, poisson, binomial, choice, hypergeometric
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe

from estimator import calculate_LD_Good2022, calculate_LD, calculate_LE

def calculate_sigmasquared(n11s,n10s,n01s,n00s,ntots,fstar):

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

now = datetime.now()
dt_string = now.strftime("%H:%M:%S")
print("Started at {}".format(dt_string))

debug = True

mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['legend.fontsize'] = 'small'

theory_xs = numpy.logspace(-3,3.5,100)

# Set up figure
f, axes = plt.subplots(2,2,squeeze=False)
f.set_size_inches(6,4.5)

ldl_axis = axes[0][0]
ldl_axis.set_ylabel("LD, $\\sigma_d^2(f_0)$")
ldl_axis.set_xlabel("Distance between SNVs, $\ell$")
# ldl_axis.set_ylim([6e-04,6e-01])
ldl_axis.set_ylim([1e-01,1e1])
ldl_axis.set_xlim([2,8e03])

# ldl_axis.fill_between([4e03,1e04],[1,1],[1e-04,1e-04],color='0.8')

ldl_axis.loglog([1e-06],[1e-06],'k.')

ldf_axis = axes[0][1]
#ldf_axis.set_ylabel("$\\sigma_d^2")
ldf_axis.set_xlabel("Frequency scale, $f_0$")
# ldf_axis.set_ylim([1e-04,6e-01])
ldf_axis.set_ylim([1e-01,1e1])
ldf_axis.set_xlim([0.8e-03,2])
ldf_axis.loglog([1e-06],[1e-06],'k.')
# ldf_axis.fill_between([1.0,2.0],[1,1],[1e-04,1e-04],color='0.8')

sfs_axis = axes[1][0]  # TODO: remove this axis
sfs_axis.set_ylabel("Fraction SNV pairs")
sfs_axis.set_xlabel("Single mutant count, $n_A$")

n11_axis = axes[1][1]  # TODO: zhiru: remove this
n11_axis.set_xlabel("Double mutant count, $n_{AB}$")

plt.subplots_adjust(hspace=0.35,wspace=0.25)


"""
Loading the ns for a species
"""
species_name = "E_Rectale"
ns = numpy.load('/Volumes/Botein/uhgg/E_rectale_site_pairs_083022.npy').astype(float)
print("Done loading")
n11s = ns[:, 0]
n10s = ns[:, 1]
n01s = ns[:, 2]
n00s = ns[:, 3]
ells = ns[:, 4]
ntots = n11s+n10s+n01s+n00s

max_ntot = ntots.max()

# keeping only well covered sites
good_idxs = (ntots>0.95*max_ntot)
print("max_ntot", max_ntot)
print("Kept", good_idxs.sum()*1.0/len(good_idxs))

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


def plot_SFS(ax, n11s, n10s, n01s, n00s, ntots):
    min_ntot = ntots.min()

    f11s = n11s*1.0/ntots
    f10s = n10s*1.0/ntots
    f01s = n01s*1.0/ntots
    f00s = n00s*1.0/ntots

    fAs = f11s+f10s
    fBs = f11s+f01s

    nAs = n11s+n10s
    nBs = n11s+n10s


    macs = numpy.fmin(nAs,ntots-nAs)
    bins = numpy.arange(0,max_ntot)+0.5
    h,bin_edges = numpy.histogram(macs,bins=bins)
    ks = numpy.arange(1,max_ntot)

    minor_h = (h+h[::-1])*1.0
    minor_h /= minor_h.sum()
    minor_h *= 2

    ax.loglog(ks,minor_h,label='Observed')
    ax.set_xlim([1,min_ntot/2.0])

    nstar = long(min_ntot/2.0)
    ylims = ax.get_ylim()
    bottom = ylims[0]
    top = ylims[1]

    neutral_h = 1.0/(ks*(max_ntot-ks))

    neutral_h = neutral_h*minor_h[nstar-1]/neutral_h[nstar-1]
    ax.loglog(ks,neutral_h,'k:',label='$1/n_A(n-n_A)$')

    ax.set_ylim([bottom,top])

    ax.legend(loc='lower left',frameon=False,numpoints=1)


def plot_double_mutant_freq(ax, n11s, n10s, n01s, n00s, ntots):
    min_ntot = ntots.min()

    f11s = n11s*1.0/ntots
    f10s = n10s*1.0/ntots
    f01s = n01s*1.0/ntots
    f00s = n00s*1.0/ntots

    fAs = f11s+f10s
    fBs = f11s+f01s

    nAs = n11s+n10s
    nBs = n11s+n10s

    fmin = 0.13
    fmax = 0.17

    ellranges = [(1e04,1e08),(1e03,2e03),(1e02,3e02)]
    colors = ['#1f77b4','#2ca02c','#ff7f0e']
    labels =['Genome avg','1000<$\ell$<2000','100<$\ell$<300']

    for idx in xrange(0,len(ellranges)):
        ellmin,ellmax = ellranges[idx]
        color = colors[idx]
        step = 10
        bins = numpy.arange(0,ntots.max(), 20)-0.5
        xs = (bins[1:]+bins[:-1])/2
        # xs = numpy.arange(0,ntots.max(), 10)

        good_idxs = (ells>ellmin)*(ells<ellmax)*(fAs<fmax)*(fBs<fmax)*(fAs>fmin)*(fBs>fmin)

        h,bin_edges = numpy.histogram(n11s[good_idxs],bins=bins)
        h = h*1.0/h.sum()
        ax.semilogy(xs,h,'.-',markersize=3,color=color)

    ax.set_xlim([-0.5,900])

    ylims = ax.get_ylim()
    bottom = ylims[0]
    top = ylims[1]
    ax.plot([min_ntot*fmin*fmin,max_ntot*fmin*fmin],[bottom,top],'k:')
    ax.plot([min_ntot*fmin,min_ntot*fmin],[bottom,top],'k:')

    ax.set_ylim([bottom,top])
    print("Done plotting nAB")

# plt.savefig('%s_example.pdf' % species_name,bbox_inches='tight')


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

for idx in xrange(0,len(ellranges)):
    ellmin,ellmax = ellranges[idx]
    color = colors[idx]
    label=labels[idx]

    bins = numpy.arange(0,ntots.max())-0.5
    xs = numpy.arange(0,ntots.max()-1)

    good_idxs = (ells>ellmin)*(ells<ellmax)

    sigmasquareds = []
    for f0 in f0s:
        print f0
        # numer, denom = calculate_sigmasquared(n11s[good_idxs],n10s[good_idxs],n01s[good_idxs],n00s[good_idxs],ntots[good_idxs],f0)
        # numer, denom = calculate_LD(n_obs[good_idxs, :],f0)
        numer, denom = calculate_LE(n_obs[good_idxs, :],f0)
        # numer, denom = calculate_LD_Good2022(n11s[good_idxs],n10s[good_idxs],n01s[good_idxs],n00s[good_idxs],ntots[good_idxs],f0)
        sigmasquared = numer.mean() / denom.mean()
        sigmasquareds.append(sigmasquared)

    sigmasquareds = numpy.array(sigmasquareds)
    line, = ldf_axis.loglog(capped_f0s[:-1],sigmasquareds[:-1],'-',markersize=3,color=color,label=label)
    ldf_axis.loglog(capped_f0s[-2:],sigmasquareds[-2:],'.:',color=color,markersize=3)

    numpy.save("./cached/f_scan_f0s_{}".format(idx), capped_f0s)
    numpy.save("./cached/f_scan_LEs_{}".format(idx), sigmasquareds)

ldf_axis.legend(loc='upper left',frameon=False,numpoints=1)
print "Done plotting f0 scan!"
plt.savefig('%s_LE_example.pdf' % species_name,bbox_inches='tight')

"""
Fix a few f0, scan a range of distances
"""
raw_theory_ells = numpy.logspace(0,3,200)*3  # Specify the range of ells
# raw_theory_ells = numpy.logspace(0,3,20)*3  # Specify the range of ells
theory_ells = numpy.hstack([raw_theory_ells,[6e03]])  # for plotting on x axis
big_theory_ells = numpy.hstack([raw_theory_ells,[1e07]])  # for filtering data; 1e7 is chosen to represent genome wide pairs

fstars = [1e02,1e-01,3e-02, 1e-2, 3e-3, 1e-3]  # Specify the f0s to calculate
fstar_colors = ['#045a8d','#2b8cbe','#74a9cf', '#9bbbd2', '#c1cbd2', '#d2d2d2']
fstar_labels=['$f_0=\\infty$','$f_0=0.1$', '$f_0=0.03$', '$f_0=0.01$', '$f_0=0.003$', '$f_0=0.001$']

for idx in xrange(0,len(fstars)):
    fstar = fstars[idx]
    color=fstar_colors[idx]
    label=fstar_labels[idx]

    capped_fstar = min([fstar,1.0])
    print "Processing fstar =", fstar, capped_fstar
    # numer, denom = calculate_sigmasquared(n11s,n10s,n01s,n00s,ntots,fstar)
    # numer, denom = calculate_LD(n_obs,fstar)
    numer, denom = calculate_LE(n_obs,fstar)
    # numer, denom = calculate_LD_Good2022(n11s,n10s,n01s,n00s,ntots,fstar)

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

    line, = ldl_axis.loglog(avg_ells[:-1],sigmasquareds[:-1],'-',label=label,color=color)
    ldl_axis.loglog(avg_ells[-2:],sigmasquareds[-2:],'.:',color=color,markersize=3)

    numpy.save("./cached/ell_scan_ells_{}".format(idx), avg_ells)
    numpy.save("./cached/ell_scan_LEs_{}".format(idx), sigmasquareds)

ldl_axis.text(6e03,1.7e-04,'Genome\navg', horizontalalignment='center',fontsize='7')

ldl_axis.legend(loc='lower left',frameon=False,numpoints=1)

plt.savefig('%s_LE_example_2.pdf' % species_name,bbox_inches='tight')

now = datetime.now()
dt_string = now.strftime("%H:%M")
print("Finished at {}".format(dt_string))
