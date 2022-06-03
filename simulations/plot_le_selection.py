# modifying plot_selection_figure.py

import sys
import numpy
import pylab
import gzip
import parameters
from math import log10,log
import ld_theory
from numpy.random import multinomial

import matplotlib
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from numpy.random import randint, shuffle, poisson, binomial, choice, hypergeometric
import matplotlib.patheffects as pe

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-p", "--parameters", type="string",
                  help="name pf the parameters regime from parameters.py",
                  dest="regime")

(options, args) = parser.parse_args()
regime = options.regime

mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'


pylab.figure(1,figsize=(6,4.5))
f = pylab.gcf()

eta_axis = pylab.gca()
eta_axis.set_ylabel("$\\Lambda(f_0)$")
eta_axis.set_xlabel("$2 N s f_0$")

theory_xs = numpy.logspace(-6,5,25)

vmin=-3
vmax=-1
cmap='viridis'

jet = cm = pylab.get_cmap(cmap) 
cNorm  = colors.Normalize(vmin=vmin, vmax=vmax)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

n=1e05 # actual population size
#fstars = numpy.logspace(-3,-1,20)
fstars = numpy.array([0.001,0.003,0.01,0.03,0.1])
params = parameters.params
for type,symbol,counts_symbol in zip([regime],['o'],['s']):
    LEs = {fstar:[] for fstar in fstars}
    denominatorsquareds = {fstar:[] for fstar in fstars}
    bare_numeratorsquareds = {fstar:[] for fstar in fstars}
    
    gammas = []
    
    for param_idx in xrange(0,len(params[type])):

        N = params[type][param_idx][2]
        sA = params[type][param_idx][3]
        sB = params[type][param_idx][4]
        eps = params[type][param_idx][5]
        gamma = 2*(N*eps+N*sA+N*sB)
        
        gammas.append(gamma)
    
        filename = 'output/output_%s_%d.txt.gz' % (type,param_idx)
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
   
        f11s = numpy.array(f11s)
        f10s = numpy.array(f10s)
        f01s = numpy.array(f01s)
        f00s = 1-f11s-f10s-f01s
        fAs = f11s+f10s
        fBs = f11s+f01s

        n11s = numpy.array(n11s)
        n10s = numpy.array(n10s)
        n01s = numpy.array(n01s)
        ns = numpy.ones_like(n11s)*n
        n00s = ns-n10s-n11s-n01s
        
        for fstar in fstars:
        
            numerator = f11s*f00s*f10s*f01s
            sampling_variances = (fAs*(1-fAs)*fBs*(1-fBs))**2

            Hs = numpy.exp(-fAs/fstar-fBs/fstar)
    
            LE_numerator = (numerator*Hs).mean()
            LE_denominator = (sampling_variances*Hs).mean()
            LE = LE_numerator/LE_denominator
            LEs[fstar].append(LE)

    gammas = numpy.array(gammas)
    for fstar in fstars:
        LEs[fstar] = numpy.array(LEs[fstar])
        denominatorsquareds[fstar] = numpy.array(denominatorsquareds[fstar])
        bare_numeratorsquareds[fstar] = numpy.array(bare_numeratorsquareds[fstar])
        
    for fstar in fstars:
        pylab.figure(1)
        collapse_xs = gammas*fstar
        collapse_ys = LEs[fstar]
        
        colorVal = scalarMap.to_rgba(log10(fstar))

        l = '$f_0=$' + str(fstar)
        line, = eta_axis.semilogx(collapse_xs,collapse_ys,symbol,markersize=5,color=colorVal,alpha=0.7,markeredgewidth=0,label=l)


eta_axis.set_xticks([1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3])
eta_axis.set_xticklabels(['$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$10^{0}$', '$10^{1}$', '$10^{2}$', '$10^{3}$'])
#eta_axis.set_yticks([0, 0.5e-6, 1e-6, 1.5e-6, 2e-6])
#eta_axis.set_yticklabels(['$0.0$', '$0.5x10^{-6}$', '$1.0x10^{-6}$', '$1.5x10^{-6}$', '$2.0x10^{-6}$'])
#eta_axis.set_ylim([5e-6,2])
eta_axis.set_xlim([1e-3,1e3])
eta_axis.minorticks_off()


eta_axis.legend(frameon=False,loc='upper right')

pylab.savefig('LE_%s.png'%regime,dpi=600,bbox_inches='tight')
