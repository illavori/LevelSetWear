# ======================================================================
# LevelSetWear, data of Fig. 9 from the plate with a hole runs (example 01)
#
# Written by Inigo Llavori
#
# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
# ======================================================================
# Turns the text files of extractKirsch.py into the CSV files read by
# figures/fig09_kirsch/plot.py.
#   energy_error.csv       energy norm error against the element size,
#                          sqrt(|U - U_h|/U), with U_h = F.u/2 (nodal forces
#                          of makeKirsch.py and nodal displacements) and U the
#                          exact strain energy (fine midpoint rule, Richardson)
#   cut_element_error.csv  error in sigma_xx of the cut elements against their
#                          material fraction, N = 64, 128, 256 pooled, rms in
#                          bins. The element stress (average over its material
#                          part) is compared with the exact stress averaged
#                          over the same part (24 x 24 sub-sampling)
#   elements_N256.csv      sigma_xx/S and its error at the elements of the
#                          256 x 256 grid with material
#   zero_contour_N256.csv  zero contour of phi seen by the cut elements
# It also writes kirsch.txt, kirsch_p0.txt and kirsch_pt.txt (strain energy
# per grid, as reference/kirsch.txt) in the data folder.
#
#   python postKirsch.py [--data DIR] [--out DIR]
#     --data  folder of the <job>_U.txt and <job>_S.txt files (default: here)
#     --out   folder of the CSV files (default: ../../figures/fig09_kirsch)
# Grids or rules without results are skipped. Panels (c) and (d) need the
# three rules, (a) and (b) the 256 x 256 grid with the cell rule.
import os
import sys

import numpy as np

import makeKirsch as K

HERE = os.path.dirname(os.path.abspath(__file__))
RULES = ((4, 'cell rule n = 4 (present)'), (0, 'exact (material polygon)'), (-4, 'point rule'))
CUTGRIDS = (64, 128, 256)                # grids pooled in cut_element_error.csv
NFIG = 256                               # grid of panels (a) and (b)
KS = 24                                  # sub-sampling of the cut elements
FBINS = np.array([0.001, 0.01, 0.03, 0.1, 0.3, 1.0])
CITE = ['# LevelSetWear, data of Fig. 9 of',
        '# I. Llavori, Wear simulation on a fixed mesh with level sets',
        '# and immersed frictional contact,',
        '# submitted for publication. Mondragon Unibertsitatea.',
        '# https://github.com/illavori/LevelSetWear       BSD 3-Clause']
RULELINE = ('# rule of the cut elements: 0 cell rule n = 4 (present), '
            '1 exact (material polygon), 2 point rule')


def energyDensity(sxx, syy, sxy):
    """Plane strain strain energy density from the in-plane stresses."""
    return ((1 - K.NU**2)*(sxx**2 + syy**2) - 2*K.NU*(1 + K.NU)*sxx*syy
            + 2*(1 + K.NU)*sxy**2)/(2*K.E)


def exactEnergy(n):
    """Exact strain energy of the quarter plate outside the hole (per unit
    thickness), midpoint rule on an n x n grid, the cells crossed by the
    hole weighted by their material fraction (8 x 8 sub-sampling)."""
    h = K.L/n
    xc = (np.arange(n) + 0.5)*h
    X, Y = np.meshgrid(xc, xc)
    R = np.hypot(X, Y)
    w = np.where(R > K.A, 1.0, 0.0)
    band = np.abs(R - K.A) < h
    k = 8
    g = ((np.arange(k) + 0.5)/k - 0.5)*h
    GX, GY = np.meshgrid(g, g)
    xb, yb = X[band], Y[band]
    fr = (np.hypot(xb[:, None] + GX.ravel(), yb[:, None] + GY.ravel()) > K.A).mean(axis=1)
    w[band] = fr
    Xs, Ys = np.where(R > 0.5*K.A, X, 1.0), np.where(R > 0.5*K.A, Y, 1.0)
    u = energyDensity(*K.kirsch(Xs, Ys))
    return (u*w).sum()*h*h


def strainEnergy(fname, n):
    """U_h = F.u/2 with the nodal displacements of the odb."""
    d = np.loadtxt(fname)
    U = {int(r[0]): r[1:] for r in d}
    nodes = K.mesh(n)[0]
    return 0.5*sum(f @ U[k] for e in ('X', 'Y') for k, f in K.nodalLoads(nodes, e).items())


def energyErrors(data):
    """Rows (rule index, N, le/r0, error in %) and the files kirsch<tag>.txt."""
    U = (9*exactEnergy(3000) - 4*exactEnergy(2000))/5    # Richardson, O(h^2)
    rows = []
    for ir, (rule, lab) in enumerate(RULES):
        done = []
        for n in K.GRIDS:
            f = os.path.join(data, K.jobName(n, rule) + '_U.txt')
            if os.path.exists(f):
                Uh = strainEnergy(f, n)
                done.append('%d %.6e %.10e %.6e' % (n, K.L/n, Uh, np.sqrt(abs(U - Uh)/U)))
        if not done:
            print('no results of the %s' % lab)
            continue
        with open(os.path.join(data, 'kirsch%s.txt' % K.TAG[rule]), 'w') as f:
            f.write('# exact strain energy %.10e N mm/mm\n'
                    '# grid  l_e (mm)  U_h  energy norm error (rel.)\n' % U)
            f.write('\n'.join(done) + '\n')
        # the CSV is written from the rounded values of the text file
        kd = np.loadtxt(done, ndmin=2)
        for r in kd:
            rows.append((ir, r[0], r[1]/K.A, 100*r[3]))
            print('%-27s N %4d  error %.3f %%' % (lab, r[0], 100*r[3]))
    return rows


def exactAverage(x, y, h, k=KS):
    """Exact sxx averaged over the material part (r > r0) of the element of
    centre (x, y) and size h, k x k sub-sampling. None if no sample is in
    the material."""
    g = ((np.arange(k) + 0.5)/k - 0.5)*h
    GX, GY = np.meshgrid(g, g)
    px, py = x + GX.ravel(), y + GY.ravel()
    q = np.hypot(px, py) > K.A
    if not q.any():
        return None, 0.0
    return K.kirsch(px[q], py[q])[0].mean(), q.mean()


def cutErrors(fname, n):
    """(material fraction, error in % of S) of the elements cut by the hole.
    The fraction is the geometric one, the same for the three rules."""
    d = np.loadtxt(fname)
    h = K.L/n
    x, y, sxx = d[:, 1], d[:, 2], d[:, 3]
    r = [np.hypot(x + sx*h/2, y + sy*h/2) for sx in (-1, 1) for sy in (-1, 1)]
    cut = (np.min(r, axis=0) < K.A) & (np.max(r, axis=0) > K.A)
    out = []
    for i in np.nonzero(cut)[0]:
        ex, f = exactAverage(x[i], y[i], h)
        if ex is None:
            continue
        # rounded to 5 decimals as in the paper (bin edges)
        out.append((float('%.5f' % f), float('%.5f' % (100*(sxx[i] - ex)/K.S))))
    return np.array(out)


def cutRows(data):
    rows = []
    for ir, (rule, lab) in enumerate(RULES):
        e = [cutErrors(os.path.join(data, K.jobName(n, rule) + '_S.txt'), n) for n in CUTGRIDS
             if os.path.exists(os.path.join(data, K.jobName(n, rule) + '_S.txt'))]
        if not e:
            print('no ghost results of the %s' % lab)
            continue
        e = np.vstack(e)
        for lo, hi in zip(FBINS[:-1], FBINS[1:]):
            m = (e[:, 0] >= lo) & (e[:, 0] < hi)
            if m.any():
                rows.append((ir, np.sqrt(lo*hi), np.sqrt((e[m, 1]**2).mean())))
    return rows


def elementRows(fname, n):
    """Centroid, sxx/S and error in % of S of the elements with material.
    Exact stress at the centroid of uncut elements, averaged over the
    material part of cut ones."""
    d = np.loadtxt(fname)
    m = d[:, 7] > 1e-6
    x, y, sxx, fr = d[m, 1], d[m, 2], d[m, 3], d[m, 7]
    h = K.L/n
    ex = K.kirsch(x, y)[0]
    cut = fr < 0.999
    for i in np.nonzero(cut)[0]:
        a, _ = exactAverage(x[i], y[i], h)
        if a is None:
            # sliver thinner than the sampling grid: the point of the edge of
            # the hole closest to the element centre, just outside it
            t = 1.0001*K.A/np.hypot(x[i], y[i])
            a = K.kirsch(np.array([x[i]*t]), np.array([y[i]*t]))[0].mean()
        ex[i] = a
    return x, y, sxx/K.S, 100*(sxx - ex)/K.S


def chords(n):
    """Zero contour of phi = r0 - r as seen by the cut elements: crossings
    with the grid lines, ordered by angle."""
    h = K.L/n
    xn = np.linspace(0, K.L, n + 1)
    pts = []
    for c in xn:
        for vert in (True, False):
            f = K.A - np.hypot(c, xn) if vert else K.A - np.hypot(xn, c)
            s = np.nonzero((f[:-1] < 0) != (f[1:] < 0))[0]
            for j in s:
                t = f[j]/(f[j] - f[j + 1])
                v = xn[j] + t*h
                pts.append((c, v) if vert else (v, c))
            for j in np.nonzero(f == 0)[0]:
                pts.append((c, xn[j]) if vert else (xn[j], c))
    pts = np.array(sorted(set(pts), key=lambda p: np.arctan2(p[1], p[0])))
    return pts[:, 0], pts[:, 1]


def writeCsv(out, name, head, cols, rows, fmt=repr):
    with open(os.path.join(out, name), 'w', newline='\n') as f:
        f.write('\n'.join(CITE + head) + '\n' + ','.join(cols) + '\n')
        for r in rows:
            f.write(','.join(fmt(v) for v in r) + '\n')
    print('written %s' % os.path.join(out, name))


def main(argv):
    data = HERE
    out = os.path.join(HERE, '..', '..', 'figures', 'fig09_kirsch')
    if '--data' in argv:
        data = argv[argv.index('--data') + 1]
    if '--out' in argv:
        out = argv[argv.index('--out') + 1]
    os.makedirs(out, exist_ok=True)
    num = lambda v: repr(float(v))
    rows = energyErrors(data)
    if rows:
        writeCsv(out, 'energy_error.csv',
                 ['# energy norm error against the element size le, hole radius r0 = 1 mm', RULELINE],
                 ['rule', 'N', 'le_over_r0', 'energy_error_pct'],
                 [(r[0], num(r[1]), num(r[2]), num(r[3])) for r in rows], fmt=str)
    rows = cutRows(data)
    if rows:
        writeCsv(out, 'cut_element_error.csv',
                 ['# error in sigma_xx of the cut elements against their material fraction, '
                  'N = %s pooled' % ', '.join(str(n) for n in CUTGRIDS),
                  '# rms error in bins of material fraction with edges %s, geometric bin centre'
                  % ' '.join('%g' % b for b in FBINS), RULELINE],
                 ['rule', 'material_fraction', 'rms_error_pct_of_S'],
                 [(r[0], num(r[1]), num(r[2])) for r in rows], fmt=str)
    f = os.path.join(data, K.jobName(NFIG, 4) + '_S.txt')
    if os.path.exists(f):
        writeCsv(out, 'elements_N%d.csv' % NFIG,
                 ['# Plate with a hole (Kirsch), hole radius r0 = 1 mm, quarter plate 4 mm x 4 mm, '
                  'remote stress S = 100 MPa',
                  '# level set model on a fixed %d x %d grid, elements with material '
                  '(fraction > 1e-6), centroid coordinates' % (NFIG, NFIG),
                  '# sxx_over_S = sigma_xx/S of the level set model',
                  '# error_pct_of_S = 100 (sigma_xx - exact)/S, exact averaged over the '
                  'material part of the cut elements'],
                 ['x_mm', 'y_mm', 'sxx_over_S', 'error_pct_of_S'],
                 zip(*elementRows(f, NFIG)), fmt=lambda v: '%.8g' % v)
    else:
        print('no ghost results of the %d x %d grid with the cell rule' % (NFIG, NFIG))
    writeCsv(out, 'zero_contour_N%d.csv' % NFIG,
             ['# zero contour of phi = r0 - r as seen by the cut elements of the '
              '%d x %d grid' % (NFIG, NFIG),
              '# (crossings of the zero level set with the grid lines, ordered by angle)'],
             ['x_mm', 'y_mm'], zip(*chords(NFIG)), fmt=lambda v: '%.12g' % v)


if __name__ == '__main__':
    main(sys.argv[1:])
