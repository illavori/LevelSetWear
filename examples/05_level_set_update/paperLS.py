# ======================================================================
# LevelSetWear, verification of the level set update (example 05)
#
# Written by Iñigo Llavori
#
# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
# ======================================================================
# Figure of the paper (Section 4.3): narrow band level set update against
# the height function on the self-consistent Winkler wear model of
# protoNBLS.py (rigid cylinder R = 6 mm, gross slip, Archard), with the four
# operations of Section 3.4 and dx = le/16 = 0.625 um.
#   (a) worn profiles at several numbers of cycles
#   (b) level set - height function at the last one, in % of the maximum depth
#   python paperLS.py      -> levelsetUpdate.pdf / .png
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# PDF metadata: the author of the paper
META = lambda e: ({'Author': 'Iñigo Llavori', 'Creator': 'Iñigo Llavori'}
                  if e == 'pdf' else {'Author': 'Iñigo Llavori'})

import protoNBLS as pn

HERE = os.path.dirname(os.path.abspath(__file__))
pn.NSTEP = 180                           # 18 000 cycles, as the test of McColl et al.
SAVE = (45, 90, 135, 180)                # steps of 100 cycles: quarters of the test

h, ls = pn.Height(), pn.NBLS()
prof = {}
for k in range(pn.NSTEP):
    h.step()
    ls.step(k)
    if k + 1 in SAVE:
        prof[(k + 1)*pn.DN] = (h.x.copy(), h.y.copy(), ls.x.copy(), ls.profile())
xh, yh, xl, yl = prof[max(prof)]
d = yl - np.interp(xl, xh, yh)
hmax = -np.nanmin(yh)                    # maximum depth of the scar (mm)
scar = np.interp(xl, xh, yh) < -0.01*hmax   # worn part of the flat
e = 100*d/hmax
print('max depth %.4f um, scar width %.3f mm' % (1e3*hmax, np.ptp(xl[scar])))
print('over the scar: difference mean %+.3f %%, std %.3f %%, max |.| %.3f %% of the depth'
      % (np.nanmean(e[scar]), np.nanstd(e[scar]), np.nanmax(np.abs(e[scar]))))
c = np.abs(xl) < 0.5*np.ptp(xl[scar])/2
print('central half of the scar: max |.| %.3f %%' % np.nanmax(np.abs(e[c])))
print('worn area difference %.4f %%'
      % (100*(np.trapz(yl, xl) - np.trapz(yh, xh))/np.trapz(yh, xh)))

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
    ax[0].plot(xl_, 1e3*yl_, color=col, lw=1.0, label='%d cycles' % n)
    ax[0].plot(xh_[::30], 1e3*yh_[::30], '^', color=col, mfc='none', ms=3.2, mew=0.6)
ax[0].plot([], [], 'k-', lw=1.0, label='level set')
ax[0].plot([], [], 'k^', mfc='none', ms=3.2, mew=0.6, label='height function')
ax[0].legend(loc='lower center', frameon=False, ncol=2, handlelength=1.4, columnspacing=0.8)
ax[0].set_ylabel(r'$y$ [$\mu$m]')
ax[0].set_ylim(-5.5, 0.4)
ax[1].plot(xl, e, color='k', lw=0.7)
ax[1].axhline(0, color='0.25', lw=0.5, zorder=0)
ax[1].set_ylabel('difference [%]')
ax[1].set_xlabel(r'$x$ [mm]')
ax[1].set_xlim(-0.3, 0.3)
print('whole profile: max |.| %.3f %% of the depth' % np.nanmax(np.abs(e)))
for a_, lab, yt in ((ax[0], '(a)', 0.86), (ax[1], '(b)', 0.94)):
    a_.text(0.98, yt, lab, transform=a_.transAxes, ha='right', va='top', fontsize=7.5,
            bbox=dict(facecolor='white', edgecolor='none', pad=0.6))
for e_ in ('pdf', 'png'):
    fig.savefig(os.path.join(HERE, 'levelsetUpdate.' + e_), dpi=300, metadata=META(e_))
