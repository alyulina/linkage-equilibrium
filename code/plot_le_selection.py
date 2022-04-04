# Zhiru probably modified plot_recombination_figure.py

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

mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

theory_xs = numpy.logspace(-3,2.5,50)

# Set up figure
pylab.figure(1,figsize=(6,4.5))
f = pylab.gcf()

eta_axis = pylab.gca()
eta_axis.set_ylabel("$\\Lambda(f_0)$")
eta_axis.set_xlabel("$2 N s f_0$")
# eta_axis.loglog(theory_xs,numpy.ones_like(theory_xs),'k:')
# eta_axis.set_xlim([theory_xs[0],theory_xs[-1]])

theory_xs = numpy.logspace(-6,5,25)
# eta_axis.loglog(theory_xs, theory_xs, color="#d1d1d1", label='$\\Lambda(f_0)=2 N R f_0$')

vmin=-3
vmax=-1
cmap='viridis'

jet = cm = pylab.get_cmap(cmap) 
cNorm  = colors.Normalize(vmin=vmin, vmax=vmax)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

# Plot DATA

n=1e05 # actual population size
#fstars = numpy.logspace(-3,-1,20)
fstars = numpy.array([0.001,0.003,0.01,0.03,0.1])
params = parameters.params
for type,symbol,counts_symbol in zip(['selAB'],['o'],['s']):
    LEs = {fstar:[] for fstar in fstars}
    denominatorsquareds = {fstar:[] for fstar in fstars}
    bare_numeratorsquareds = {fstar:[] for fstar in fstars}
    
    # sigmasquareds_counts = {fstar:[] for fstar in fstars}
    
    gammas = []
    
    for param_idx in xrange(0,len(params[type])):

        N = params[type][param_idx][2]
        sA = params[type][param_idx][3]
        sB = params[type][param_idx][4]
        eps = params[type][param_idx][5]
        gamma = 2*(N*eps+N*sA+N*sB)
        theta = 1.0/log(N)
    
        thetaA = ld_theory.calculate_theta(N,sA) # what is this??
        thetaB = ld_theory.calculate_theta(N,sB) # what is this??
        # print N,sA,sB,eps,gamma
    
        # if gamma<20:
            # print gamma
            # pass
            # continue
        
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

            #denominatorsquareds[fstar].append(sigmasquared_denominator/thetaA/thetaB)
            
            #bare_numeratorsquareds[fstar].append(( numpy.square(bare_Ds)*Hs).mean()/thetaA/thetaB)

    gammas = numpy.array(gammas)
    for fstar in fstars:
        LEs[fstar] = numpy.array(LEs[fstar])
        denominatorsquareds[fstar] = numpy.array(denominatorsquareds[fstar])
        bare_numeratorsquareds[fstar] = numpy.array(bare_numeratorsquareds[fstar])
        
        # sigmasquareds_counts[fstar] = numpy.array(sigmasquareds_counts[fstar])
        
    for fstar in fstars:
        pylab.figure(1)
        collapse_xs = gammas*fstar
        collapse_ys = LEs[fstar]
        
        colorVal = scalarMap.to_rgba(log10(fstar))

        l = '$f_0=$' + str(fstar)
        #line, = eta_axis.loglog(collapse_xs,collapse_ys,symbol,markersize=5,color=colorVal,alpha=0.7,markeredgewidth=0,label=l)
        line, = eta_axis.plot(collapse_xs,collapse_ys,symbol,markersize=5,color=colorVal,alpha=0.7,markeredgewidth=0,label=l)

    # theory_ys = (1 - 1 / theory_xs) / (1 + 8 / theory_xs)
    # eta_axis.semilogx(theory_xs,theory_ys,'-',color='k',zorder=0)

#eta_axis.set_xticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4])
#eta_axis.set_xticklabels(['$10^{-5}$', '$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$10^{0}$', '$10^{1}$', '$10^{2}$', '$10^{3}$', '$10^{4}$'])
#eta_axis.set_yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
#eta_axis.set_yticklabels(['$10^{-5}$', '$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$10^{0}$'])
#eta_axis.set_ylim([5e-6,2])
#eta_axis.set_xlim([5e-6,2e4])
#eta_axis.minorticks_off()

# this is not working
#cax = f.add_axes([0.95, 0.95, 0.62, 0.02])
#cbar = f.colorbar(line,cax=cax,orientation='vertical',ticks=[-3,-2,-1])
#cbar.ax.tick_params(labelsize=8) 
#f.text(0.94,0.94,'$\log_{10} f_0$')

#eta_axis.legend(frameon=False,loc='upper right',numpoints=1,scatterpoints=1)
eta_axis.legend(frameon=False,loc='lower right')

pylab.savefig('LE_sAB.png',dpi=600,bbox_inches='tight')
