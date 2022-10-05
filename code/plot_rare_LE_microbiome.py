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

mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['legend.fontsize'] = 'small'

theory_xs = numpy.logspace(-3,3.5,100)

# Set up figure
f, axes = plt.subplots(1,2,squeeze=False)
f.set_size_inches(6,2.)

ldl_axis = axes[0][0]
ldl_axis.set_ylabel("LE, $\\Lambda(f_0)$")
ldl_axis.set_xlabel("Distance between SNVs, $\ell$")
# ldl_axis.set_ylim([6e-04,6e-01])
ldl_axis.set_ylim([0.3e-1,3e1])
ldl_axis.set_xlim([2,8e03])

ldl_axis.fill_between([4e03,1e04],[100,100],[1e-04,1e-04],color='0.8')

ldl_axis.loglog([1e-06],[1e-06],'k.')

ldf_axis = axes[0][1]
#ldf_axis.set_ylabel("$\\sigma_d^2")
ldf_axis.set_xlabel("Frequency scale, $f_0$")
# ldf_axis.set_ylim([1e-04,6e-01])
ldf_axis.set_ylim([0.3e-1,3e1])
ldf_axis.set_xlim([0.8e-03,2])
ldf_axis.loglog([1e-06],[1e-06],'k.')
ldf_axis.fill_between([1.0,2.0],[100,100],[1e-04,1e-04],color='0.8')

plt.subplots_adjust(hspace=0.35,wspace=0.25)

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
    capped_f0s = numpy.load("./cached/f_scan_f0s_{}.npy".format(idx))
    sigmasquareds = numpy.load("./cached/f_scan_LEs_{}.npy".format(idx))
    line, = ldf_axis.loglog(capped_f0s[:-1],sigmasquareds[:-1],'-',markersize=3,color=color,label=label)
    ldf_axis.loglog(capped_f0s[-2:],sigmasquareds[-2:],'.:',color=color,markersize=3)

ldf_axis.legend(loc='lower center', bbox_to_anchor=(0.5, 1), ncol=2, frameon=False,numpoints=1)
print "Done plotting f0 scan!"

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

    avg_ells = numpy.load("./cached/ell_scan_ells_{}.npy".format(idx))
    sigmasquareds = numpy.load("./cached/ell_scan_LEs_{}.npy".format(idx))
    line, = ldl_axis.loglog(avg_ells[:-1],sigmasquareds[:-1],'-',label=label,color=color)
    ldl_axis.loglog(avg_ells[-2:],sigmasquareds[-2:],'.:',color=color,markersize=3)

# ldl_axis.text(6e03,1.7e-04,'Genome\navg', horizontalalignment='center',fontsize='7')

ldl_axis.legend(loc='lower center', bbox_to_anchor=(0.5, 1), ncol=3, frameon=False,numpoints=1)

plt.savefig('LE_example.pdf', bbox_inches='tight')
