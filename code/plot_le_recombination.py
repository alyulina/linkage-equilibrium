# modifying plot_recombination_figure.py

import parameters
import ld_theory

import gzip
import numpy as np
from math import log10

import matplotlib
matplotlib.use('Agg') # uncomment if running on the cluster
import matplotlib.pyplot as plt

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-p", "--parameters", type="string",
                  help="name pf the parameters regime from parameters.py",
                  dest="regime")

(options, args) = parser.parse_args()
regime = options.regime

plt.figure(1,figsize=(6,4.5))
f = plt.gcf()

eta_axis = plt.gca()
eta_axis.set_ylabel("$\\Lambda(f_0)$", fontsize='small')
eta_axis.set_xlabel("$2 N R f_0$", fontsize='small')

theory_xs = np.logspace(-6,5,25)
eta_axis.loglog(theory_xs, theory_xs, color="#d1d1d1", label='$\\Lambda(f_0)=2 N R f_0$')
eta_axis.loglog(theory_xs, np.ones_like(theory_xs), color="#d1d1d1")

vmin=-3
vmax=-1
cmap='viridis'

jet = cm = plt.get_cmap(cmap) 
cNorm  = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=jet)

n=1e05 # actual population size
fstars = np.array([0.001,0.003,0.01,0.03,0.1])
params = parameters.params
for type,symbol,counts_symbol in zip([regime],['o'],['s']):
    LEs = {fstar:[] for fstar in fstars}
    denominatorsquareds = {fstar:[] for fstar in fstars}
    bare_numeratorsquareds = {fstar:[] for fstar in fstars}
    
    gammas = []
    
    for param_idx in xrange(0,len(params[type])):
    
        N = params[type][param_idx][2]
        r = params[type][param_idx][6]
        gamma = 2*N*r
        
        gammas.append(gamma)
    
        # filename = 'output/output_%s_%d.txt.gz' % (type,param_idx)
        filename = '/scratch/users/alyulina/recombination/output/output_%s_%d.txt.gz' % (type,param_idx) # uncomment if running on the cluster
        file = gzip.GzipFile(filename,"r")
        f11s = []
        f10s = []
        f01s = []
        
        # for counts (regenerated on the fly)
        n11s = []
        n10s = []
        n01s = []
        print "Loading %s" % filename
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
            
        
            n11 = long(f11*n)
            n10 = long(f10*n)
            n01 = long(f01*n)
            
            n11s.append(n11)
            n10s.append(n10)
            n01s.append(n01)
            
        file.close()
   
        f11s = np.array(f11s)
        f10s = np.array(f10s)
        f01s = np.array(f01s)
        f00s = 1-f11s-f10s-f01s
        fAs = f11s+f10s
        fBs = f11s+f01s

        n11s = np.array(n11s)
        n10s = np.array(n10s)
        n01s = np.array(n01s)
        ns = np.ones_like(n11s)*n
        n00s = ns-n10s-n11s-n01s
        
        for fstar in fstars:
        
            numerator = f11s*f00s*f10s*f01s
            sampling_variances = (fAs*(1-fAs)*fBs*(1-fBs))**2

            Hs = np.exp(-fAs/fstar-fBs/fstar)
    
            LE_numerator = (numerator*Hs).mean()
            LE_denominator = (sampling_variances*Hs).mean()
            LE = LE_numerator/LE_denominator
            LEs[fstar].append(LE)

    gammas = np.array(gammas)
    for fstar in fstars:
        LEs[fstar] = np.array(LEs[fstar])
        denominatorsquareds[fstar] = np.array(denominatorsquareds[fstar])
        bare_numeratorsquareds[fstar] = np.array(bare_numeratorsquareds[fstar])
        
    for fstar in fstars:
        plt.figure(1)
        collapse_xs = gammas *fstar
        collapse_ys = LEs[fstar]

        # saving data! f* \t x \t y \n
        with open('LE_%s_points.txt'%regime, 'a+') as o:
        	o.writelines(['{:.1g}'.format(fstar) + '\t' + str(x[0]) + '\t' + str(x[1]) + '\n' for x in zip(collapse_xs, collapse_ys)])
        
        colorVal = scalarMap.to_rgba(log10(fstar))

        l = '$f_0=$' + str(fstar)
        line, = eta_axis.loglog(collapse_xs,collapse_ys,symbol,markersize=5,color=colorVal,alpha=0.7,markeredgewidth=0,label=l)

eta_axis.set_xticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4])
eta_axis.set_xticklabels(['$10^{-5}$', '$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$10^{0}$', '$10^{1}$', '$10^{2}$', '$10^{3}$', '$10^{4}$'])
eta_axis.set_yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
eta_axis.set_yticklabels(['$10^{-5}$', '$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$10^{0}$'])
eta_axis.set_ylim([5e-6,2])
eta_axis.set_xlim([5e-6,2e4])
eta_axis.minorticks_off()
eta_axis.tick_params(axis='both', labelsize=8)

eta_axis.legend(frameon=False,fontsize=8,loc='lower right')

plt.savefig('LE_%s.png'%regime,dpi=600,bbox_inches='tight') 
