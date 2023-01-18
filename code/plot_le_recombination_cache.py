import parameters
import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('Agg')
import config

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-p", "--parameters", type="string",
                  help="name pf the parameters regime from parameters.py",
                  dest="regime")
(options, args) = parser.parse_args()
regime = options.regime

# loading intermediate files
savepath = os.path.join(config.CACHED_DIR, 'sherlock', 'simulated')
rs = np.load(os.path.join(savepath, '%s_rs.npy' % regime))
fstars = np.load(os.path.join(savepath, '%s_fstars.npy' % regime))
numers = np.load(os.path.join(savepath, '%s_numers.npy' % regime))
denoms = np.load(os.path.join(savepath, '%s_denoms.npy' % regime))
# norms = np.load(os.path.join(savepath, '%s_large_norms' % regime), denoms)

# set up figure
mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['legend.fontsize'] = 'small'

fig, axes = plt.subplots(1, 3, figsize=(6.5, 1.5))
plt.subplots_adjust(wspace=0.28)
collapsed_ax, r_ax, f_ax = axes
vmin=-4
vmax=-0.5
cmap='viridis'
jet = cm = plt.get_cmap(cmap)
cNorm  = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
scalarMap = mpl.cm.ScalarMappable(norm=cNorm, cmap=jet)

# Plotting collapsed scan
freq_pairs = list(enumerate(fstars))
for j, fstar in freq_pairs[:-5]:
    # the last five are 1e-3, 3e-3, 1e-2, 3e-2, 1e-1
    xs = []
    ys = []
    for i, r in enumerate(rs):
        N = parameters.params[regime][i][2]
        rho = 2 * N * r
        LE = numers[i, j] / denoms[i, j]
        xs.append(rho)
        ys.append(LE)
    colorVal = scalarMap.to_rgba(np.log10(fstar))
    xs = np.array(xs)
    ys = np.array(ys)
    line, = collapsed_ax.loglog(xs*fstar, ys, 'o', markersize=2, color=colorVal, alpha=0.7,
                            markeredgewidth=0)

# Plotting R scan
for j, fstar in freq_pairs[-5:][::-1]:
    # the last five are 1e-3, 3e-3, 1e-2, 3e-2, 1e-1
    xs = []
    ys = []
    for i, r in enumerate(rs):
        N = parameters.params[regime][i][2]
        rho = 2 * N * r
        LE = numers[i, j] / denoms[i, j]
        xs.append(rho)
        ys.append(LE)
    colorVal = scalarMap.to_rgba(np.log10(fstar))
    xs = np.array(xs)
    ys = np.array(ys)
    l = '$f_0=$' + str(fstar)
    line, = r_ax.loglog(xs, ys, 'o', markersize=2, color=colorVal, alpha=0.7,
                                markeredgewidth=0, label=l)

# Plotting f0 scan
r_pairs = list(enumerate(rs))
good_r_pairs = [r_pairs[ind] for ind in [8, 14, 20]]
colors = ['#1f77b4','#2ca02c','#ff7f0e']
for idx, pair in enumerate(good_r_pairs):
    i, r = pair
    xs = []
    ys = []
    N = parameters.params[regime][i][2]
    rho = 2 * N * r
    for j, fstar in freq_pairs[:-5]:
        LE = numers[i, j] / denoms[i, j]
        xs.append(fstar)
        ys.append(LE)
    xs = np.array(xs)
    ys = np.array(ys)
    color = colors[idx]
    l = '$\\rho={:.1e}$'.format(rho)
    line, = f_ax.loglog(xs, ys, 'o', markersize=2, color=color, alpha=0.7,
                        markeredgewidth=0, label=l)

theory_xs = np.logspace(-6,5,25)
collapsed_ax.loglog(theory_xs, theory_xs, color="#d1d1d1", label='$\\Lambda(f_0)=2 N R f_0$', zorder=1)
collapsed_ax.loglog(theory_xs, np.ones_like(theory_xs), color="#d1d1d1", zorder=1, linestyle=':')
r_ax.loglog(theory_xs, np.ones_like(theory_xs), color="#d1d1d1", zorder=1, linestyle=':')
f_ax.loglog(theory_xs, np.ones_like(theory_xs), color="#d1d1d1", zorder=1, linestyle=':')

# set up figure axis
# collapsed_ax.set_xticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4])
# collapsed_ax.set_xticklabels(
#     ['$10^{-5}$', '$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$10^{0}$', '$10^{1}$', '$10^{2}$', '$10^{3}$',
#      '$10^{4}$'])
# collapsed_ax.set_yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
# collapsed_ax.set_yticklabels(['$10^{-5}$', '$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$10^{0}$'])
collapsed_ax.set_ylim([5e-6, 3])
collapsed_ax.set_xlim([5e-6, 2e4])
# collapsed_ax.minorticks_off()
collapsed_ax.tick_params(axis='both')
collapsed_ax.legend(frameon=False, loc='lower right')
collapsed_ax.set_ylabel(r"$\Lambda$")
collapsed_ax.set_xlabel(r"$2NRf_0$")

r_ax.set_ylim([5e-6, 3])
r_ax.set_xlim([1e-2, 1e5])
r_ax.legend(frameon=False, loc='lower right')
r_ax.set_xlabel("$2NR$")

f_ax.set_ylim([5e-6, 3])
f_ax.set_xlim([0.8e-4, 0.5])
f_ax.legend(frameon=False, loc='lower right')
f_ax.set_xlabel("$f_0$")

plt.savefig('LE_%s.pdf' % regime, dpi=600, bbox_inches='tight')