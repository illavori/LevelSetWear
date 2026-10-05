# ======================================================================
# LevelSetWear, figure 15 of the article (block-on-ring, Cruzado 2010)
#
# Written by Iñigo Llavori
#
# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
# ======================================================================
# Worn surface of the block of the block-on-ring test of Cruzado et al.
# (2010), 65 N, after 800, 1600, 2400 and 3600 m of sliding. Level set
# model (lines) and mesh motion model on the same mesh (crosses), with the
# three profiles measured after 3600 m (grey), digitised from Cruzado,
# Zabala, Urchegui and Gomez, Revista de Metalurgia 46 (2010) 106-114,
# Fig. 12.
# Reads the CSV files of this folder.
#   python plot.py      -> bor.pdf / bor.png
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


plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Times New Roman'],
                     'mathtext.fontset': 'stix', 'font.size': 10,
                     'axes.labelsize': 10, 'axes.titlesize': 10,
                     'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 8})
DS = (800, 1600, 2400, 3600)
SHADES = ['#7fd0d2', '#00a0a5', '#006d70', '#003d3d']      # light to dark with sliding

me = load('measured_cruzado2010.csv')
ls, um = load('level_set.csv'), load('mesh_motion.csv')
W1 = 252/72.27
fig, ax = plt.subplots(figsize=(W1, 0.72*W1))
fig.subplots_adjust(left=0.17, right=0.97, bottom=0.17, top=0.97)
for k in (1, 2, 3):
    ax.plot(me['x_mm'], me['test%d_um' % k], color='0.6', lw=0.8)
hs = []
for d, col in zip(DS, SHADES):
    ml, mu = ls['sliding_m'] == d, um['sliding_m'] == d
    hs += ax.plot(ls['x_mm'][ml], ls['y_um'][ml], color=col, lw=1.3)
    ax.plot(um['x_mm'][mu], um['y_um'][mu], 'x', color=col, ms=3, mew=0.7)
st1 = [ax.plot([], [], color='0.6', lw=0.8)[0], ax.plot([], [], color='k', lw=1.3)[0],
       ax.plot([], [], 'kx', ms=3, mew=0.7)[0]]
ax.add_artist(ax.legend(hs, ['%d' % d for d in DS], title='sliding [m]', frameon=False,
                        loc='lower right', handlelength=1.5, borderaxespad=0.2))
ax.legend(st1, ['measured', 'level set', 'mesh motion'], frameon=False,
          loc='lower left', handlelength=1.8, borderaxespad=0.2)
ax.set_xlim(-2.0, 2.0)
ax.set_ylim(-75, 8)
ax.set_xlabel('$x$ [mm]')
ax.set_ylabel(r'worn surface of the block [$\mu$m]')
for e in ('pdf', 'png'):
    fig.savefig(os.path.join(HERE, 'bor.' + e), dpi=300, metadata=META(e))
