# ======================================================================
# LevelSetWear, figure 13 of the article (narrow band level set update)
#
# Written by Iñigo Llavori
#
# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
# ======================================================================
# Narrow band level set update against the height function on a
# self-consistent Winkler wear model (rigid cylinder R = 6 mm, gross slip,
# Archard), with the four operations of Section 3.4, kernel of radius
# le = 10 um and grid spacing le/16 = 0.625 um.
#   (a) worn profiles at several numbers of cycles
#   (b) level set - height function at the last one, in % of the maximum depth
# Reads the CSV files of this folder.
#   python plot.py      -> levelsetUpdate.pdf / levelsetUpdate.png
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
META = lambda e: ({'Author': 'Iñigo Llavori', 'Creator': 'Iñigo Llavori'}
                  if e == 'pdf' else {'Author': 'Iñigo Llavori'})


def load(name):
    """Columns of a CSV file of this folder, {header: array}."""
    with open(os.path.join(HERE, name), encoding='utf-8') as f:
        lines = [l for l in f if not l.startswith('#')]
    d = np.loadtxt(lines[1:], delimiter=',', ndmin=2)
    return {c: d[:, i] for i, c in enumerate(lines[0].strip().split(','))}


H, LS = load('height_function.csv'), load('level_set.csv')
prof = {}                                # cycles: x, y (um) of both methods
for n in np.unique(H['cycles']):
    mh, ml = H['cycles'] == n, LS['cycles'] == n
    prof[int(n)] = (H['x_mm'][mh], H['y_um'][mh], LS['x_mm'][ml], LS['y_um'][ml])
xh, yh, xl, yl = prof[max(prof)]
d = yl - np.interp(xl, xh, yh)

plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Times New Roman'],
                     'mathtext.fontset': 'stix', 'font.size': 10, 'axes.labelsize': 10,
                     'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 8})
W = 252/72.27
fig, ax = plt.subplots(2, 1, figsize=(W, 1.05*W), sharex=True,
                       gridspec_kw=dict(height_ratios=[2.2, 1]))
fig.subplots_adjust(left=0.18, right=0.97, bottom=0.12, top=0.98, hspace=0.08)
shades = ['#7fd0d2', '#00a0a5', '#006d70', '#003d3d']   # light to dark with the cycles
for (n, (xh_, yh_, xl_, yl_)), col in zip(sorted(prof.items()), shades):
    # level set as lines, the height function (exact reference) as open triangles
    ax[0].plot(xl_, yl_, color=col, lw=1.0, label='%d cycles' % n)
    ax[0].plot(xh_[::30], yh_[::30], '^', color=col, mfc='none', ms=3.2, mew=0.6)
ax[0].plot([], [], 'k-', lw=1.0, label='level set')
ax[0].plot([], [], 'k^', mfc='none', ms=3.2, mew=0.6, label='height function')
ax[0].legend(loc='lower center', frameon=False, ncol=2, handlelength=1.4, columnspacing=0.8)
ax[0].set_ylabel(r'$y$ [$\mu$m]')
ax[0].set_ylim(-5.5, 0.4)
hmax = -np.nanmin(yh)                    # maximum depth of the scar (um)
ax[1].plot(xl, 100*d/hmax, color='k', lw=0.7)
ax[1].axhline(0, color='0.25', lw=0.5, zorder=0)
ax[1].set_ylabel('difference [%]')
ax[1].set_xlabel(r'$x$ [mm]')
ax[1].set_xlim(-0.3, 0.3)
for a_, lab, yt in ((ax[0], '(a)', 0.86), (ax[1], '(b)', 0.94)):
    a_.text(0.98, yt, lab, transform=a_.transAxes, ha='right', va='top', fontsize=7.5,
            bbox=dict(facecolor='white', edgecolor='none', pad=0.6))
for e in ('pdf', 'png'):
    fig.savefig(os.path.join(HERE, 'levelsetUpdate.' + e), dpi=300, metadata=META(e))
