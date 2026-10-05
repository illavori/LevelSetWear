# ======================================================================
# LevelSetWear, figure 17 of the article (Newton iterations, fretting)
#
# Written by Iñigo Llavori
#
# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
# ======================================================================
# Newton iterations per increment along the fretting test of McColl et
# al. (18 000 cycles, one point per wear block of 100 cycles), level set
# model and mesh motion model.
# Reads the CSV files of this folder.
#   python plot.py      -> iterations.pdf / iterations.png
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
                     'axes.labelsize': 10, 'xtick.labelsize': 9,
                     'ytick.labelsize': 9, 'legend.fontsize': 8})
W1 = 252/72.27
fig, ax = plt.subplots(figsize=(W1, 0.72*W1))
fig.subplots_adjust(left=0.15, right=0.97, bottom=0.17, top=0.97)
for name, col, lab in (('level_set.csv', '#003d3d', 'level set'),
                       ('mesh_motion.csv', '#00a0a5', 'mesh motion')):
    d = load(name)
    ax.plot(d['cycles'], d['iterations_per_increment'], color=col, lw=1.1, label=lab)
ax.set_xlim(0, 18000)
ax.set_ylim(0, 5)
ax.set_yticks([0, 1, 2, 3, 4, 5])
ax.set_xlabel('cycles')
ax.set_ylabel('iterations per increment')
ax.legend(frameon=False, loc='lower right', handlelength=1.8, borderaxespad=0.3)
for e in ('pdf', 'png'):
    fig.savefig(os.path.join(HERE, 'iterations.' + e), dpi=300, metadata=META(e))
