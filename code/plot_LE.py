import sys
import numpy
import pylab
import gzip
import parameters
from math import log10,log
import ld_theory
from numpy.random import multinomial

import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from numpy.random import randint, shuffle, poisson, binomial, choice, hypergeometric
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe

mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

theory_xs = numpy.logspace(-3,3.5,50)

# Set up figure
pylab.figure(1,figsize=(3.42,2.5))
f = pylab.gcf()

eta_axis = pylab.gca()
eta_axis.set_ylabel("$\\Lambda(f_0)$")
eta_axis.set_xlabel("$2 N R f_0$")
eta_axis.loglog(theory_xs,numpy.ones_like(theory_xs),'k:')
# eta_axis.set_xlim([theory_xs[0],theory_xs[-1]])

theory_xs = numpy.logspace(-3,0,25)
eta_axis.loglog(theory_xs, theory_xs, 'k--', label='y=x')

vmin=-3
vmax=-1
cmap='jet_r'

jet = cm = pylab.get_cmap(cmap) 
cNorm  = colors.Normalize(vmin=vmin, vmax=vmax)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

# Plot DATA

n=1e05 # actual population size
#fstars = numpy.logspace(-3,-1,20)
fstars = numpy.array([0.001,0.01,0.03,0.1])
params = parameters.params
for type,symbol,counts_symbol in zip(['r'],['o'],['s']):
    LEs = {fstar:[] for fstar in fstars}
    denominatorsquareds = {fstar:[] for fstar in fstars}
    bare_numeratorsquareds = {fstar:[] for fstar in fstars}
    
    sigmasquareds_counts = {fstar:[] for fstar in fstars}
    
    gammas = []
    
    for param_idx in xrange(0,len(params[type])):
    
        N = params[type][param_idx][2]
        r = params[type][param_idx][6]
        gamma = 2*N*r
        theta = 1.0/log(N)
    
        # if gamma<5:
        #     print gamma
        #     pass
        #     continue
        
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

            # denominatorsquareds[fstar].append(sigmasquared_denominator/theta**2)
            #
            # bare_numeratorsquareds[fstar].append(( numpy.square(bare_Ds)*Hs).mean()/theta**2)

    gammas = numpy.array(gammas)
    for fstar in fstars:
        LEs[fstar] = numpy.array(LEs[fstar])
        denominatorsquareds[fstar] = numpy.array(denominatorsquareds[fstar])
        bare_numeratorsquareds[fstar] = numpy.array(bare_numeratorsquareds[fstar])
        
        sigmasquareds_counts[fstar] = numpy.array(sigmasquareds_counts[fstar])
        
    for fstar in fstars:
        pylab.figure(1)
        collapse_xs = gammas *fstar
        collapse_ys = LEs[fstar]
        
        colorVal = scalarMap.to_rgba(log10(fstar))
        
        line, = eta_axis.loglog(collapse_xs,collapse_ys,symbol,markersize=3,color=colorVal)

    # theory_ys = (1 - 1 / theory_xs) / (1 + 8 / theory_xs)
    # eta_axis.semilogx(theory_xs,theory_ys,'-',color='k',zorder=0)


eta_axis.set_ylim([1e-4,2])
eta_axis.set_xlim([1e-4,3e03])

eta_axis.legend(frameon=False,loc='lower left',numpoints=1,scatterpoints=1)

pylab.savefig('LE.pdf',bbox_inches='tight')
