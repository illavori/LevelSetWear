# ======================================================================
# LevelSetWear, figure 12 of the article (von Mises stress, Mindlin)
#
# Written by Iñigo Llavori
#
# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
# ======================================================================
# Von Mises stress below a cylinder on a flat in loading (Q = Q*,
# Mindlin-Cattaneo), (a) standard Abaqus contact, (b) level set model on
# the same mesh and (c) their difference, stacked in one column.
# Reads the CSV file of this folder.
#   python plot.py      -> stressMC.pdf / stressMC.png
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LinearSegmentedColormap

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
teal_div = LinearSegmentedColormap.from_list(
    'MUdiv', [(0.0, '#4d4d4d'), (0.25, '#bdbdbd'), (0.5, '#ffffff'),
              (0.67, '#7fd0d2'), (0.83, '#00a0a5'), (1.0, '#005a5d')])
lev = np.linspace(0, 0.7, 22)
dm = 2.5
d_ = load('von_mises_loading.csv')
x, y = d_['x_over_a'], d_['y_over_a']
va, vu = d_['vm_standard_over_p0'], d_['vm_level_set_over_p0']
d = 100*(vu - va)                         # difference in % of p0
tri = mtri.Triangulation(x, y)

W1 = 252/72.27
pw1 = 2.50
YD = 1.5                                  # depth shown (in a)
ph1 = pw1*YD/5
GAP = 0.08
HF1 = 0.40 + 3*ph1 + 2*GAP + 0.08
IN = lambda l, b, w, h: [l/W1, b/HF1, w/W1, h/HF1]
fig = plt.figure(figsize=(W1, HF1))
bot = lambda k: 0.40 + (2 - k)*(ph1 + GAP)
for ir, (v, lab) in enumerate(((va, '(a) standard'), (vu, '(b) level set'), (d, '(c) difference'))):
    ax = fig.add_axes(IN(0.42, bot(ir), pw1, ph1))
    if ir < 2:
        cf = ax.tricontourf(tri, v, levels=lev, cmap='jet')
    else:
        cd = ax.tricontourf(tri, np.clip(v, -dm, dm),
                            levels=np.r_[np.linspace(-dm, -0.25, 10), np.linspace(0.25, dm, 10)],
                            cmap=teal_div)
    ax.text(-2.45, -YD + 0.05, lab, ha='left', va='bottom', fontsize=7.5,
            bbox=dict(facecolor='white', edgecolor='none', pad=0.8))
    ax.set_aspect('equal')
    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-YD, 0)
    ax.set_yticks([-1, 0])
    ax.set_ylabel('$y/a$')
    if ir == 2:
        ax.set_xlabel('$x/a$')
    else:
        ax.set_xticklabels([])
cax = fig.add_axes(IN(0.42 + pw1 + 0.07, bot(1), 0.08, 2*ph1 + GAP))
cb = fig.colorbar(cf, cax=cax, ticks=np.arange(0, 0.71, 0.1))
cb.set_label(r'$\sigma_{\mathrm{vM}}/p_0$')
cax2 = fig.add_axes(IN(0.42 + pw1 + 0.07, bot(2), 0.08, ph1))
cb2 = fig.colorbar(cd, cax=cax2, ticks=[-2, 0, 2])
cb2.set_label(r'$\Delta\sigma_{\mathrm{vM}}/p_0$ [%]')
for e in ('pdf', 'png'):
    fig.savefig(os.path.join(HERE, 'stressMC.' + e), dpi=300, metadata=META(e))
