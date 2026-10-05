# ======================================================================
# LevelSetWear, figure 11 of the article (surface stresses, Mindlin)
#
# Written by Iñigo Llavori
#
# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
# ======================================================================
# Surface stresses over a fretting half cycle of a cylinder on a flat,
# (a) loading to Q* (Mindlin-Cattaneo) and (b) unloading to Q = 0
# (Mindlin-Deresiewicz), analytical solution against standard Abaqus
# contact and the level set model on the same mesh.
# Reads the CSV files of this folder.
#   python plot.py      -> surfaceMD.pdf / surfaceMD.png
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
an = load('analytical.csv')
zb = load('zone_boundaries.csv')
c, cp = zb['c_over_a'][0], zb['cprime_over_a'][0]
cs, cu = load('contact_standard.csv'), load('contact_level_set.csv')
ss, su = load('sxx_standard.csv'), load('sxx_level_set.csv')
xs = an['x_over_a']

# three shades of teal, dark to light (grey is kept for measured data)
cols = {'sxx': '#004b4e', 'p': '#1a8a8c', 'q': '#6cc3c4'}
mk = dict(ms=3, mew=0.6)

# one column, loading above unloading, shared x axis
WC = 252/72.27
fig, axs = plt.subplots(2, 1, figsize=(WC, 1.36*WC), sharex=True)
fig.subplots_adjust(left=0.15, right=0.97, bottom=0.09, top=0.955, hspace=0.16)
cases = ((1, 'loading', '(a) loading, $Q = Q^*$'),
         (2, 'unloading', '(b) unloading, $Q = 0$'))
# stick zone boundary c (loading) and counter-slip boundary c' (unloading)
marks = ([(-c, '$-c$'), (c, '$c$')],
         [(-cp, "$-c'$"), (-c, '$-c$'), (c, '$c$'), (cp, "$c'$")])
for ax, (case, st, title), mks in zip(axs, cases, marks):
    # Cartesian axes
    ax.axhline(0, color='0.25', lw=0.5, zorder=0)
    ax.axvline(0, color='0.25', lw=0.5, zorder=0)
    for xv, lab in mks:
        ax.axvline(xv, color='0.6', lw=0.5, ls=':', zorder=0)
        ax.text(xv, 1.01, lab, ha='center', va='bottom', fontsize=8, color='0.3',
                transform=ax.get_xaxis_transform())
    ax.plot(xs, an['sxx_over_p0_' + st], color=cols['sxx'], lw=1.0, label=r'$\sigma_{xx}$')
    ax.plot(xs, an['p_over_p0'], color=cols['p'], lw=1.0, label=r'$|\sigma_{yy}|=p$')
    ax.plot(xs, an['q_over_p0_' + st], color=cols['q'], lw=1.0, label=r'$\sigma_{xy}=q$')
    ma, ms_ = cs['case'] == case, ss['case'] == case
    for x, v, col in ((cs['x_over_a'][ma], cs['p_over_p0'][ma], 'p'),
                      (cs['x_over_a'][ma], cs['q_over_p0'][ma], 'q'),
                      (ss['x_over_a'][ms_], ss['sxx_over_p0'][ms_], 'sxx')):
        ax.plot(x, v, 'o', mfc='none', mec=cols[col], **mk)
    mu_, ms_ = cu['case'] == case, su['case'] == case
    for x, v, col in ((cu['x_over_a'][mu_], cu['p_over_p0'][mu_], 'p'),
                      (cu['x_over_a'][mu_], cu['q_over_p0'][mu_], 'q'),
                      (su['x_over_a'][ms_], su['sxx_over_p0'][ms_], 'sxx')):
        ax.plot(x, v, 'x', color=cols[col], **mk)
    # panel label inside the top right corner, load case in the caption
    ax.text(0.98, 0.95, title[:3], transform=ax.transAxes, ha='right', va='top',
            fontsize=7.5, bbox=dict(facecolor='white', edgecolor='none', pad=0.6))
    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.55, 1.08)   # free band at the bottom for the legend
    ax.set_yticks([-1, -0.5, 0, 0.5, 1])
    ax.set_ylabel(r'$\sigma_{ij}/p_0$')
axs[1].set_xlabel('$x/a$')
# legends inside the panels, in the free band below the curves:
# the stress components in (a), the models in (b)
LEG = dict(loc='lower center', frameon=True, facecolor='white', edgecolor='none',
           framealpha=1, fancybox=False, ncol=3, handlelength=1.4,
           handletextpad=0.4, columnspacing=1.2, fontsize=7.5, borderaxespad=0.3)
h, l = axs[0].get_legend_handles_labels()
axs[0].legend(h, l, **LEG)
axs[1].legend([plt.Line2D([], [], color='k', lw=1.0),
               plt.Line2D([], [], ls='none', marker='o', mfc='none', mec='k', **mk),
               plt.Line2D([], [], ls='none', marker='x', color='k', **mk)],
              ['analytical', 'standard', 'level set'], **LEG)
for e in ('pdf', 'png'):
    fig.savefig(os.path.join(HERE, 'surfaceMD.' + e), dpi=300, metadata=META(e))
