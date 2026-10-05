# ======================================================================
# LevelSetWear, figure 9 of the article (plate with a hole, Kirsch)
#
# Written by Iñigo Llavori
#
# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
# ======================================================================
# Plate with a hole on a fixed grid cut by the level set phi = r0 - r.
#   (a) sigma_xx/S of the level set model, 256 x 256 grid, quarter plate
#   (b) error against the exact solution, in % of S, scale clipped at +-0.05 %
#   (c) energy norm error against the element size
#   (d) error of the cut elements against their material fraction
# Reads the CSV files of this folder.
#   python plot.py      -> kirsch.pdf / kirsch.png
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.path import Path

HERE = os.path.dirname(os.path.abspath(__file__))
META = lambda e: ({'Author': 'Iñigo Llavori', 'Creator': 'Iñigo Llavori'}
                  if e == 'pdf' else {'Author': 'Iñigo Llavori'})
A, L = 1.0, 4.0                          # hole radius, size of the quarter plate (mm)


def load(name):
    """Columns of a CSV file of this folder, {header: array}."""
    with open(os.path.join(HERE, name), encoding='utf-8') as f:
        lines = [l for l in f if not l.startswith('#')]
    d = np.loadtxt(lines[1:], delimiter=',', ndmin=2)
    return {c: d[:, i] for i, c in enumerate(lines[0].strip().split(','))}


el = load('elements_N256.csv')
x, y, sxx, err = el['x_mm'], el['y_mm'], el['sxx_over_S'], el['error_pct_of_S']
zc = load('zero_contour_N256.csv')
xc, yc = zc['x_mm'], zc['y_mm']


def extend(xx, yy, v, xs_, ys_):
    """Adds the zero contour and the symmetry lines with the value of the
    nearest element, so that the contour plot reaches the boundaries."""
    bx = np.r_[xs_, np.zeros(40), np.linspace(A, L, 40), np.linspace(0, L, 60), np.full(60, L)]
    by = np.r_[ys_, np.linspace(A, L, 40), np.zeros(40), np.full(60, L), np.linspace(0, L, 60)]
    j = np.argmin((xx[None, :] - bx[:, None])**2 + (yy[None, :] - by[:, None])**2, axis=1)
    return np.r_[xx, bx], np.r_[yy, by], np.r_[v, v[j]]


def masked(xx, yy, poly):
    t = mtri.Triangulation(xx, yy)
    xm = xx[t.triangles].mean(axis=1)
    ym = yy[t.triangles].mean(axis=1)
    t.set_mask(poly.contains_points(np.c_[xm, ym]))
    return t


plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Times New Roman'],
                     'mathtext.fontset': 'stix', 'font.size': 10, 'axes.labelsize': 10,
                     'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 8})
W = 252/72.27
HF = 5.15
IN = lambda l, b, w, hh: [l/W, b/HF, w/W, hh/HF]
pw = 1.38
fig = plt.figure(figsize=(W, HF))
B0 = 3.12
axA = fig.add_axes(IN(0.24, B0, pw, pw))
axB = fig.add_axes(IN(0.24 + pw + 0.28, B0, pw, pw))
axC = fig.add_axes(IN(0.50, 1.90, W - 0.62, 1.05))
axD = fig.add_axes(IN(0.50, 0.40, W - 0.62, 1.05))
void = Path(np.c_[np.r_[0, xc, 0], np.r_[0, yc, 0]])
lev = np.linspace(0, 3.2, 17)
xe, ye, ve = extend(x, y, sxx, xc, yc)
c = axA.tricontourf(masked(xe, ye, void), ve, levels=lev, cmap='jet', extend='both')
axA.fill(np.r_[0, xc, 0], np.r_[0, yc, 0], color='0.82', lw=0, zorder=2)
div = LinearSegmentedColormap.from_list('MUdiv', [(0.0, '#4d4d4d'), (0.25, '#bdbdbd'),
                                                  (0.5, '#ffffff'), (0.67, '#7fd0d2'),
                                                  (0.83, '#00a0a5'), (1.0, '#005a5d')])
# error scale clipped at +-dm % (larger errors saturate), white band +-dm/10 %
dm = 0.05
xe, ye, ee = extend(x, y, err, xc, yc)
cd = axB.tricontourf(masked(xe, ye, void), ee,
                     levels=np.r_[np.linspace(-dm, -dm/10, 10), np.linspace(dm/10, dm, 10)],
                     cmap=div, extend='both')
for ax in (axA, axB):
    ax.plot(xc, yc, color=(126 / 255, 7 / 255, 42 / 255), lw=1.2, zorder=4)
    ax.set_aspect('equal')
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    # no axes: the frame is kept as the outline of the quarter plate (0 to 4 r0)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.6)
axA.annotate(r'$\phi = 0$', xy=(np.cos(0.9), np.sin(0.9)), xytext=(0.08, 0.3), fontsize=7, va='center',
             color=(126 / 255, 7 / 255, 42 / 255),
             arrowprops=dict(arrowstyle='-', color=(126 / 255, 7 / 255, 42 / 255), lw=0.5),
             zorder=5)
for ax, lab in ((axA, '(a)'), (axB, '(b)')):
    ax.text(3.9, 3.9, lab, ha='right', va='top', fontsize=7.5, zorder=6,
            bbox=dict(facecolor='white', edgecolor='none', pad=0.6))
cb = fig.colorbar(c, cax=fig.add_axes(IN(0.24, B0 + pw + 0.10, pw, 0.07)),
                  orientation='horizontal', ticks=[0, 1, 2, 3])
cb.set_label(r'$\sigma_{xx}/S$', labelpad=5)
cb.ax.xaxis.set_label_position('top')
cb.ax.xaxis.set_ticks_position('top')
cb2 = fig.colorbar(cd, cax=fig.add_axes(IN(0.24 + pw + 0.28, B0 + pw + 0.10, pw, 0.07)),
                   orientation='horizontal', ticks=[-dm, 0, dm])
cb2.set_label(r'error in $\sigma_{xx}/S$ [%]', labelpad=5)
cb2.ax.xaxis.set_label_position('top')
cb2.ax.xaxis.set_ticks_position('top')
# rules of the cut elements: 0 cell rule (present), 1 material polygon, 2 points
STYLES = ((1, 'exact (material polygon)',
           dict(color='0.55', ls='--', marker='s', mfc='none', ms=3.5, lw=1.0, zorder=2)),
          (0, 'cell rule, $n = 4$ (present)',
           dict(color='#005a5d', marker='o', mfc='none', ms=4, lw=1.0, zorder=3)),
          (2, 'point rule',
           dict(color='#00a0a5', marker='x', ms=4, lw=0.9, zorder=2)))
# (c) convergence of the three integration rules of the cut elements
en = load('energy_error.csv')
for r, lab, st in STYLES:
    m = en['rule'] == r
    axC.loglog(en['le_over_r0'][m], en['energy_error_pct'][m], label=lab, **st)
# slope triangle below the curves, legs labelled as in the usual convergence plots
x0, x1, y0 = 0.05, 0.2, 0.16
axC.plot([x0, x1, x1, x0], [y0, y0, y0*x1/x0, y0], color='0.3', lw=0.6)
axC.text(np.sqrt(x0*x1), y0/1.12, '1', fontsize=8, ha='center', va='top')
axC.text(x1*1.08, y0*np.sqrt(x1/x0), '1', fontsize=8, va='center')
axC.set_xlabel(r'element size $\ell_e/r_0$', labelpad=1)
axC.set_xlim(0.003, 0.33)
axC.set_ylim(0.07, 9)
axC.set_xticks([0.005, 0.01, 0.02, 0.05, 0.1, 0.2])
axC.set_xticklabels(['0.005', '0.01', '0.02', '0.05', '0.1', '0.2'])
axC.set_yticks([0.1, 0.3, 1, 3])
axC.set_yticklabels(['0.1', '0.3', '1', '3'])
axC.minorticks_off()
axC.set_ylabel('energy error [%]')
axC.legend(frameon=False, loc='upper left', bbox_to_anchor=(0.07, 1.0), fontsize=7.5, handlelength=2.2)
axC.text(0.03, 0.95, '(c)', transform=axC.transAxes, ha='left', va='top', fontsize=7.5)
# (d) error of the cut elements against their material fraction (binned rms)
ce = load('cut_element_error.csv')
for r, lab, st in STYLES:
    m = ce['rule'] == r
    axD.loglog(ce['material_fraction'][m], ce['rms_error_pct_of_S'][m], label=lab, **st)
# stresses of cut elements with less than 1 % of material are not used
axD.axvspan(0.002, 0.01, color='0.9', lw=0, zorder=0)
axD.text(0.0045, 1.35, 'not used', ha='center', va='bottom', fontsize=7.5, color='0.35')
axD.set_xlabel('material fraction of the cut element', labelpad=1)
axD.set_ylabel(r'error in $\sigma_{xx}$ [% of $S$]')
axD.set_xlim(0.002, 1)
axD.set_xticks([0.003, 0.01, 0.03, 0.1, 0.3, 1])
axD.set_xticklabels(['0.003', '0.01', '0.03', '0.1', '0.3', '1'])
axD.set_yticks([1, 3, 10, 30, 100])
axD.set_yticklabels(['1', '3', '10', '30', '100'])
axD.minorticks_off()
axD.text(0.97, 0.95, '(d)', transform=axD.transAxes, ha='right', va='top', fontsize=7.5)
for e_ in ('pdf', 'png'):
    fig.savefig(os.path.join(HERE, 'kirsch.' + e_), dpi=300, metadata=META(e_))
