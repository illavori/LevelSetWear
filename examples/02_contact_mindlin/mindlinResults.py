# ======================================================================
# LevelSetWear, results of the Mindlin contact example (example 02)
#
# Written by Inigo Llavori
#
# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
# ======================================================================
# Cylinder on a flat over a fretting half cycle, loading to Q* (Mindlin-
# Cattaneo) and unloading to Q = 0 (Mindlin-Deresiewicz). From the results
# of the jobs of this example it
#   - writes the data of Fig. 11 (surface stresses) and Fig. 12 (von Mises
#     stress) of the article, the CSV files read by
#     figures/fig11_surface_mindlin/plot.py and figures/fig12_stress_mindlin/plot.py
#   - prints the numbers of Table 1 (error of the surface stresses of the
#     aligned level set model and of the standard contact against the
#     analytical solution) and of Table 2 (immersed surfaces)
#
# Steps
#   1. run the jobs (run.bat in aligned/, standard_abaqus/, immersed/ and,
#      for the other immersed cases of Table 2, immersed/makeImmersed.py)
#   2. extract the odb results (Abaqus Python)
#        cd aligned          & abaqus python ..\extractOdb.py uel_flat_MD.odb
#        cd standard_abaqus  & abaqus python ..\extractOdb.py std_flat_MD.odb
#   3. python mindlinResults.py              -> CSV files in figures/fig11_*, fig12_*
#      python mindlinResults.py --out DIR    -> CSV files in DIR instead
#      then python plot.py in each figure folder
# The contact output of the level set jobs is contactOut.txt in the folder
# of each job. When a job has not been run, its reference result
# (reference/ folders) is used for the tables, and the script says so.
# Figures 11 and 12 need the odb results of step 2.
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, '..', '..', 'figures')

# ---------------------------------------------------------------- parameters
# plane strain cylinder (R1) on a flat, both steel, normal load P per unit
# length, tangential load Q* = QRATIO mu P
E, NU, MU = 200000.0, 0.3, 0.6        # MPa, -, -
P, QRATIO, R1 = 120.0, 0.5, 6.0       # N/mm, -, mm
Es = E/(2*(1 - NU**2))                # contact modulus of two equal solids
a = np.sqrt(4*P*R1/(np.pi*Es))        # half contact width (Hertz)
p0 = 2*P/(np.pi*a)                    # peak pressure
Qs = QRATIO*MU*P
c = a*np.sqrt(1 - QRATIO)             # stick zone in loading (Mindlin-Cattaneo)
cp = a*np.sqrt(1 - Qs/(2*MU*P))       # counter-slip boundary at Q = 0 (Mindlin-Deresiewicz)

# immersed cases of Table 2: name, depth of the surface in the flat and in
# the cylinder (fractions of the surface element size le)
IMMERSED = (('F01', 0.01, 0.0), ('F25', 0.25, 0.0), ('F50', 0.50, 0.0), ('F75', 0.75, 0.0),
            ('F99', 0.99, 0.0), ('F150', 1.50, 0.0), ('F50C50', 0.50, 0.50))

# ---------------------------------------------------------------- analytical
# a traction A sqrt(1 - x^2/r^2) on the half-plane and its surface sigma_xx
xs = np.linspace(-1.6, 1.6, 3201)*a
root = lambda x, r: np.sqrt(np.clip(1 - (x/r)**2, 0, None))
out = lambda x, r: np.where(np.abs(x) > r, np.sign(x)*np.sqrt(np.clip((x/r)**2 - 1, 0, None)), 0.0)
T = lambda A, r: A*root(xs, r)
S = lambda A, r: -2*A*(xs/r - out(xs, r))
p_an = p0*root(xs, a)
qL = T(MU*p0, a) - T(MU*p0*c/a, c)
sL = -p_an + S(MU*p0, a) - S(MU*p0*c/a, c)
qU = qL - 2*(T(MU*p0, a) - T(MU*p0*cp/a, cp))
sU = sL - 2*(S(MU*p0, a) - S(MU*p0*cp/a, cp))


# ---------------------------------------------------------------- readers
def find(*names):
    """First existing file of the list, None if there is none."""
    for n in names:
        if os.path.isfile(os.path.join(HERE, n)):
            return os.path.join(HERE, n)
    return None


def contactOut(path, stepno):
    """x, p, q of the surface points of the flat in contactOut.txt (UEL)."""
    rows, st = [], None
    for line in open(path):
        t = line.split()
        if t[0] == 'STEP':
            st = int(t[1])
        elif st == stepno:
            rows.append([float(v) for v in t[1:4]])
    u = np.array(rows)
    return u[:, 0], u[:, 1], u[:, 2]


def stdContact(path):
    """x, CPRESS, CSHEAR1 of the surface nodes (standard model)."""
    d = np.loadtxt(path)
    return d[:, 0], d[:, 1], d[:, 2]


def surfaceSxx(path):
    """sigma_xx at the surface nodes from the nodal displacements of the
    first row of elements (bilinear quadrilaterals, plane strain), averaged
    over the elements that share the node."""
    lam, G = E*NU/((1 + NU)*(1 - 2*NU)), E/(2*(1 + NU))
    XN, EN = np.array([-1, 1, 1, -1.0]), np.array([-1, -1, 1, 1.0])
    acc = {}
    lines = open(path).read().split('\n')
    for i, l in enumerate(lines):
        if not l.startswith('E'):
            continue
        d = np.array([[float(v) for v in lines[i + k].split()] for k in range(1, 5)])
        lab, xy, uv = d[:, 0].astype(int), d[:, 1:3], d[:, 3:5]
        for k in range(4):
            if abs(xy[k, 1]) > 1e-9:
                continue
            dN = 0.25*np.array([XN*(1 + EN*EN[k]), EN*(1 + XN*XN[k])])
            g = np.linalg.solve(dN @ xy, dN) @ uv
            acc.setdefault(lab[k], []).append((xy[k, 0], (lam + 2*G)*g[0, 0] + lam*g[1, 1]))
    r = np.array([[v[0][0], np.mean([t[1] for t in v])] for v in acc.values()])
    r = r[np.argsort(r[:, 0])]
    return r[:, 0], r[:, 1]


def vonMises(path):
    """{element label: (xc, yc, von Mises)} of the element-average stresses."""
    d = np.loadtxt(path)
    s11, s22, s33, s12 = d[:, 3], d[:, 4], d[:, 5], d[:, 6]
    v = np.sqrt(0.5*((s11 - s22)**2 + (s22 - s33)**2 + (s33 - s11)**2) + 3*s12**2)
    return {int(l): (x, y, w) for l, x, y, w in zip(d[:, 0], d[:, 1], d[:, 2], v)}


def below(A, U):
    """Elements of both models below the contact, |x| < 2.5a, -2.5a < y < 0."""
    return [l for l in U if l in A and U[l][1] < 0 and abs(U[l][0]) < 2.5*a and U[l][1] > -2.5*a]


# ---------------------------------------------------------------- CSV output
CITE = ('# LevelSetWear, data of Fig. %d of\n'
        '# I. Llavori, Wear simulation on a fixed mesh with level sets\n'
        '# and immersed frictional contact,\n'
        '# submitted for publication. Mondragon Unibertsitatea.\n'
        '# https://github.com/illavori/LevelSetWear       BSD 3-Clause\n'
        '# cylinder on flat, half-contact width a = %.6f mm, peak pressure p0 = %.2f MPa, '
        'friction coefficient %g, Q*/(mu P) = %g\n')
CASES = '# case: 1 loading to Q* (Mindlin-Cattaneo), 2 unloading to Q = 0 (Mindlin-Deresiewicz)\n'


def writeCsv(folder, name, fig, comments, header, rows):
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(folder, name), 'w', encoding='utf-8', newline='\n') as f:
        f.write(CITE % (fig, a, p0, MU, QRATIO) + comments + header + '\n')
        for r in rows:
            f.write(','.join(str(int(v)) if isinstance(v, int) else repr(float(v)) for v in r) + '\n')
    print('  written', os.path.normpath(os.path.join(folder, name)))


def figureData(outDir):
    # odb results of extractOdb.py
    need = ['aligned/uel_flat_MD_%s_%s.txt' % (k, s) for k in 'US' for s in ('TANGENTIAL', 'UNLOAD')]
    need += ['standard_abaqus/std_flat_MD_%s_%s.txt' % (k, s)
             for k in ('U', 'S', 'contact') for s in ('TANGENTIAL', 'UNLOAD')]
    miss = [n for n in need if find(n) is None]
    co = find('aligned/contactOut.txt')
    if miss or co is None:
        print('Figures 11 and 12 not written, missing (run the jobs and extractOdb.py):')
        for n in miss + ([] if co else ['aligned/contactOut.txt']):
            print('   ', n)
        return None
    f11 = outDir or os.path.join(FIGS, 'fig11_surface_mindlin')
    f12 = outDir or os.path.join(FIGS, 'fig12_stress_mindlin')
    data = {}
    rows = {'cu': [], 'cs': [], 'su': [], 'ss': []}
    for case, (st, sn) in enumerate((('TANGENTIAL', 3), ('UNLOAD', 4)), start=1):
        xa, pa, qa = stdContact(find('standard_abaqus/std_flat_MD_contact_%s.txt' % st))
        xu, pu, qu = contactOut(co, sn)
        xsa, sa = surfaceSxx(find('standard_abaqus/std_flat_MD_U_%s.txt' % st))
        xsu, su = surfaceSxx(find('aligned/uel_flat_MD_U_%s.txt' % st))
        data[st] = dict(std=(xa, pa, qa, xsa, sa), uel=(xu, pu, qu, xsu, su))
        rows['cs'] += [(case, x/a, p/p0, q/p0) for x, p, q in zip(xa, pa, qa)]
        rows['cu'] += [(case, x/a, p/p0, q/p0) for x, p, q in zip(xu, pu, qu)]
        rows['ss'] += [(case, x/a, s/p0) for x, s in zip(xsa, sa)]
        rows['su'] += [(case, x/a, s/p0) for x, s in zip(xsu, su)]

    print('Data of Fig. 11')
    writeCsv(f11, 'analytical.csv', 11,
             '# analytical surface stresses, loading to Q* (Mindlin-Cattaneo) and unloading '
             'to Q = 0 (Mindlin-Deresiewicz)\n'
             '# p = |sigma_yy| contact pressure, q = sigma_xy shear traction, sxx = sigma_xx '
             'at the surface\n',
             'x_over_a,p_over_p0,q_over_p0_loading,sxx_over_p0_loading,q_over_p0_unloading,'
             'sxx_over_p0_unloading',
             zip(xs/a, p_an/p0, qL/p0, sL/p0, qU/p0, sU/p0))
    writeCsv(f11, 'contact_level_set.csv', 11,
             '# level set model (UEL), contact pressure p and shear traction q on the '
             'immersed surface\n' + CASES, 'case,x_over_a,p_over_p0,q_over_p0', rows['cu'])
    writeCsv(f11, 'contact_standard.csv', 11,
             '# standard Abaqus contact (penalty), contact pressure p and shear traction q '
             'at the surface nodes\n' + CASES, 'case,x_over_a,p_over_p0,q_over_p0', rows['cs'])
    writeCsv(f11, 'sxx_level_set.csv', 11,
             '# level set model (UEL), sigma_xx at the surface nodes (from the nodal '
             'displacements)\n' + CASES, 'case,x_over_a,sxx_over_p0', rows['su'])
    writeCsv(f11, 'sxx_standard.csv', 11,
             '# standard Abaqus contact, sigma_xx at the surface nodes (from the nodal '
             'displacements)\n' + CASES, 'case,x_over_a,sxx_over_p0', rows['ss'])
    writeCsv(f11, 'zone_boundaries.csv', 11,
             "# stick zone boundary c (loading) and counter-slip boundary c' (unloading)\n",
             'c_over_a,cprime_over_a', [(c/a, cp/a)])

    print('Data of Fig. 12')
    A = vonMises(find('standard_abaqus/std_flat_MD_S_TANGENTIAL.txt'))
    U = vonMises(find('aligned/uel_flat_MD_S_TANGENTIAL.txt'))
    com = below(A, U)
    writeCsv(f12, 'von_mises_loading.csv', 12,
             '# von Mises stress at the element centroids in loading (Q = Q*, '
             'Mindlin-Cattaneo),\n'
             '# standard Abaqus contact and level set model (UEL) on the same mesh, '
             '|x| < 2.5 a, -2.5 a < y < 0\n',
             'x_over_a,y_over_a,vm_standard_over_p0,vm_level_set_over_p0',
             [(U[l][0]/a, U[l][1]/a, A[l][2]/p0, U[l][2]/p0) for l in com])
    return data


# ---------------------------------------------------------------- tables
pm = lambda e: '%+.1f +- %.1f' % (e.mean(), e.std())


def table1(data):
    """Error against the analytical solution at the same x, in % of p0,
    over |x| <= a for p and q and |x| <= 1.5a for sigma_xx."""
    print('\nTable 1. Error of the surface stresses against the analytical solution,')
    print('mean +- standard deviation in % of p0')
    print('%-20s %-10s %-14s %-14s %-14s' % ('', 'Model', 'p', 'q', 'sxx'))
    for st, qan, san, lab in (('TANGENTIAL', qL, sL, 'Loading, Q = Q*'),
                              ('UNLOAD', qU, sU, 'Unloading, Q = 0')):
        for kind, name in (('std', 'Standard'), ('uel', 'Level set')):
            xp, pp, qq, xss, ss = data[st][kind]
            m = np.abs(xp) <= a
            ep = 100*(pp[m] - np.interp(xp[m], xs, p_an))/p0
            eq = 100*(qq[m] - np.interp(xp[m], xs, qan))/p0
            ms = np.abs(xss) <= 1.5*a
            es = 100*(ss[ms] - np.interp(xss[ms], xs, san))/p0
            print('%-20s %-10s %-14s %-14s %-14s' % (lab, name, pm(ep), pm(eq), pm(es)))
            lab = ''
    # von Mises stress of the two models, element by element (Fig. 12)
    print('\nVon Mises stress, level set - standard at the same element (% of p0)')
    for st in ('TANGENTIAL', 'UNLOAD'):
        A = vonMises(find('standard_abaqus/std_flat_MD_S_%s.txt' % st))
        U = vonMises(find('aligned/uel_flat_MD_S_%s.txt' % st))
        d = np.array([100*(U[l][2] - A[l][2])/p0 for l in below(A, U)])
        print('%-10s %+.2f +- %.2f  (%d elements)' % (st, d.mean(), d.std(), len(d)))


def analyticalAt(x, stepno):
    p = p0*root(x, a)
    q = MU*p0*root(x, a) - MU*p0*c/a*root(x, c)
    if stepno == 4:
        q = q - 2*(MU*p0*root(x, a) - MU*p0*cp/a*root(x, cp))
    return p, q


def table2():
    """Error of p and q against the analytical solution with the surfaces
    immersed, over |x| < a, and the largest difference node by node from
    the aligned configuration, in % of p0."""
    ref = find('aligned/contactOut.txt', 'aligned/reference/contactOut.txt')
    cases = [('aligned', 0.0, 0.0, ref)]
    for name, dF, dC in IMMERSED:
        f = find('immersed/%s/contactOut.txt' % name,
                 *(['immersed/contactOut.txt'] if name == 'F50C50' else []),
                 'immersed/reference/contactOut_imm_%s.txt' % name)
        cases.append((name, dF, dC, f))
    print('\nTable 2. Error of the surface stresses against the analytical solution with')
    print('the surfaces immersed, mean +- standard deviation in % of p0 over |x| < a,')
    print('and max |immersed - aligned| node by node (% of p0)')
    print('%-6s %-6s | %-14s %-14s | %-14s %-14s | %-11s %-11s | %s'
          % ('dyF/le', 'dyC/le', 'p loading', 'q loading', 'p unloading', 'q unloading',
             'max dp', 'max dq', 'source'))
    for name, dF, dC, f in cases:
        if f is None:
            print('%-6.2f %-6.2f   (no result, run the job)' % (dF, dC))
            continue
        cols, dmax = [], [0.0, 0.0]
        for sn in (3, 4):
            x0, p0_, q0_ = contactOut(ref, sn)
            x, p, q = contactOut(f, sn)
            m = np.abs(x0) < a
            pa, qa = analyticalAt(x0[m], sn)
            cols += [pm(100*(p[m] - pa)/p0), pm(100*(q[m] - qa)/p0)]
            dmax[0] = max(dmax[0], 100*np.abs(p[m] - p0_[m]).max()/p0)
            dmax[1] = max(dmax[1], 100*np.abs(q[m] - q0_[m]).max()/p0)
        print('%-6.2f %-6.2f | %-14s %-14s | %-14s %-14s | %-11.2f %-11.2f | %s'
              % ((dF, dC) + tuple(cols) + tuple(dmax) + (os.path.relpath(f, HERE),)))


if __name__ == '__main__':
    outDir = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else None
    print('a = %.6f mm, p0 = %.2f MPa, c/a = %.4f, c\'/a = %.4f' % (a, p0, c/a, cp/a))
    data = figureData(outDir)
    if data is not None:
        table1(data)
    table2()
