# ======================================================================
# LevelSetWear, input files of the plate with a hole (example 01)
#
# Written by Inigo Llavori
#
# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
# ======================================================================
# Infinite plate with a circular hole under uniaxial tension (Kirsch), the
# verification of the cut elements (Section 4.1, Fig. 9 of the article).
# Quarter plate [0, L] x [0, L], hole of radius r0 = A at the origin, remote
# stress S in x. Symmetry on x = 0 (ux = 0) and y = 0 (uy = 0), exact
# Kirsch tractions on x = L and y = L (consistent nodal forces). The hole is
# the level set phi = r0 - r (circle option of uelWearBulk.for, NSB < 0 in
# the header of surf1.txt) on a regular grid of N x N elements, plane strain.
# Grids up to N = 256 carry a ghost mesh of CPE4 elements that shows the
# stresses of the cut elements in Abaqus/Viewer (and in the odb).
#
# Integration rule of the cut elements (PROPS(6) = NFIX of the UEL)
#   --rule 4    cell rule, 4 x 4 points (present)   kirsch_<N>.inp
#   --rule 0    exact, material polygon             kirsch_p0_<N>.inp
#   --rule -4   point rule (verification only)      kirsch_pt_<N>.inp
#   --rule all  the three of them
#
#   python makeKirsch.py                     -> kirsch_64.inp, surf1.txt, surf2.txt
#   python makeKirsch.py 16 32 64 --rule all -> the three rules on three grids
#   python makeKirsch.py --paper             -> the 21 inputs of Fig. 9
#                                               (N = 16 to 1024, three rules)
# Then run each job with the UEL, for instance
#   abq2020 job=kirsch_64 user=uelWearBulk.for interactive
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
E, NU, EPSV = 200000.0, 0.3, 1e-6        # Young modulus (MPa), Poisson, void stiffness factor
A, L, S = 1.0, 4.0, 100.0                # hole radius, plate size (mm), remote stress (MPa)
GHOST = 100000                           # label offset of the ghost elements
NGHOST = 256                             # ghost mesh up to this grid
GRIDS = (16, 32, 64, 128, 256, 512, 1024)
TAG = {4: '', 0: '_p0', -4: '_pt'}
RULE = {4: '', 0: ' (exact rule, material polygon)', -4: ' (point rule)'}
CITE = ('LevelSetWear, written by Inigo Llavori. Lines after the data are not read.\n'
        'If you use it, please cite I. Llavori, Wear simulation on a fixed mesh with level sets\n'
        'and immersed frictional contact, submitted for publication. Mondragon Unibertsitatea.\n'
        'https://github.com/illavori/LevelSetWear       BSD 3-Clause\n')
HEAD = ['** ' + '=' * 70,
        '** LevelSetWear, input file of an example',
        '**',
        '** Written by Inigo Llavori',
        '**',
        '** I. Llavori, Wear simulation on a fixed mesh with level sets',
        '** and immersed frictional contact,',
        '** submitted for publication. Mondragon Unibertsitatea.',
        '** https://github.com/illavori/LevelSetWear       BSD 3-Clause',
        '** ' + '=' * 70]


def kirsch(x, y):
    """Exact stresses (sxx, syy, sxy) of the infinite plate with a hole."""
    r2 = x**2 + y**2
    t = np.arctan2(y, x)
    q2, q4 = A**2/r2, A**4/r2**2
    srr = S/2*(1 - q2) + S/2*(1 - 4*q2 + 3*q4)*np.cos(2*t)
    stt = S/2*(1 + q2) - S/2*(1 + 3*q4)*np.cos(2*t)
    srt = -S/2*(1 + 2*q2 - 3*q4)*np.sin(2*t)
    c, s = np.cos(t), np.sin(t)
    sxx = srr*c*c + stt*s*s - 2*srt*s*c
    syy = srr*s*s + stt*c*c + 2*srt*s*c
    sxy = (srr - stt)*s*c + srt*(c*c - s*s)
    return sxx, syy, sxy


def mesh(n):
    """Nodes {label: (x, y)}, elements {label: 4 nodes} and symmetry sets."""
    nid = lambda i, j: j*(n + 1) + i + 1
    h = L/n
    nodes = {nid(i, j): (i*h, j*h) for j in range(n + 1) for i in range(n + 1)}
    elems = {}
    e = 1
    for j in range(n):
        for i in range(n):
            elems[e] = (nid(i, j), nid(i+1, j), nid(i+1, j+1), nid(i, j+1))
            e += 1
    sets = {'XSYM': [nid(0, j) for j in range(n + 1)],
            'YSYM': [nid(i, 0) for i in range(n + 1)]}
    return nodes, elems, sets


def nodalLoads(nodes, edge):
    """Consistent nodal forces of the exact traction on the edge x = L
    (edge 'X') or y = L (edge 'Y'), 3-point Gauss per segment."""
    k = 0 if edge == 'X' else 1
    on = sorted((n for n, c in nodes.items() if abs(c[k] - L) < 1e-12),
                key=lambda n: nodes[n][1 - k])
    F = {n: np.zeros(2) for n in on}
    gp, gw = np.polynomial.legendre.leggauss(3)
    for n1, n2 in zip(on[:-1], on[1:]):
        s1, s2 = nodes[n1][1 - k], nodes[n2][1 - k]
        ln = s2 - s1
        for p, w in zip(gp, gw):
            sp = 0.5*(s1 + s2) + 0.5*ln*p
            x, y = (L, sp) if edge == 'X' else (sp, L)
            sxx, syy, sxy = kirsch(x, y)
            t = np.array([sxx, sxy]) if edge == 'X' else np.array([sxy, syy])
            N1, N2 = 0.5*(1 - p), 0.5*(1 + p)
            F[n1] += N1*t*0.5*ln*w
            F[n2] += N2*t*0.5*ln*w
    return F


def jobName(n, rule):
    return 'kirsch%s_%d' % (TAG[rule], n)


def writeInp(n, rule, folder=HERE):
    name = jobName(n, rule)
    nodes, elems, sets = mesh(n)
    L_ = HEAD + ['*HEADING',
                 'Plate with a hole (Kirsch), fixed %d x %d grid, cut elements%s' % (n, n, RULE[rule]),
                 '*NODE, NSET=NALL']
    L_ += ['%d, %.15g, %.15g' % (k, x, y) for k, (x, y) in sorted(nodes.items())]
    L_ += ['*USER ELEMENT, NODES=4, TYPE=U1, PROPERTIES=6, COORDINATES=2,'
           ' VARIABLES=%d' % (5 + abs(rule)**2), '1, 2', '*ELEMENT, TYPE=U1, ELSET=EALL']
    L_ += ['%d, %d, %d, %d, %d' % ((e,) + c) for e, c in sorted(elems.items())]
    for s in ('XSYM', 'YSYM'):
        ns = sets[s]
        L_ += ['*NSET, NSET=%s' % s]
        L_ += [', '.join(str(k) for k in ns[i:i+12]) for i in range(0, len(ns), 12)]
    L_ += ['*UEL PROPERTY, ELSET=EALL', '%.12g, %.12g, 1., 1., %.12g, %d' % (E, NU, EPSV, rule)]
    ghost = n <= NGHOST
    if ghost:
        L_ += ['*ELEMENT, TYPE=CPE4, ELSET=GHOST']
        L_ += ['%d, %d, %d, %d, %d' % ((e + GHOST,) + c) for e, c in sorted(elems.items())]
        L_ += ['*SOLID SECTION, ELSET=GHOST, MATERIAL=GHOST', '1.',
               '*MATERIAL, NAME=GHOST', '*USER MATERIAL, CONSTANTS=1', '%d' % GHOST,
               '*DEPVAR', '5', '1, S11', '2, S22', '3, S33', '4, S12', '5, VFRAC']
    L_ += ['*STEP', '*STATIC', '*BOUNDARY', 'XSYM, 1, 1', 'YSYM, 2, 2', '*CLOAD']
    F = {}
    for edge in ('X', 'Y'):
        for k, f in nodalLoads(nodes, edge).items():
            F[k] = F.get(k, np.zeros(2)) + f
    for k in sorted(F):
        for d in (0, 1):
            if abs(F[k][d]) > 0:
                L_ += ['%d, %d, %.15g' % (k, d + 1, F[k][d])]
    L_ += ['*ENERGY PRINT', '*OUTPUT, FIELD', '*NODE OUTPUT', 'U']
    if ghost:
        L_ += ['*ELEMENT OUTPUT, ELSET=GHOST', 'SDV']
    L_ += ['*END STEP']
    with open(os.path.join(folder, name + '.inp'), 'w', newline='\r\n') as f:
        f.write('\n'.join(L_) + '\n')
    return name


def writeSurfaces(folder=HERE):
    """surf1.txt: the hole (2 points, NSB = -1 selects the circle with centre
    (0, 0) and radius r0). surf2.txt: a second body far away (not used)."""
    with open(os.path.join(folder, 'surf1.txt'), 'w', newline='\r\n') as f:
        f.write('2 -1\n0 0\n%.15g 0\n' % A + CITE)
    with open(os.path.join(folder, 'surf2.txt'), 'w', newline='\r\n') as f:
        f.write('2 1\n-10 100\n10 100\n' + CITE)


def main(argv):
    rules = [4]
    argv = list(argv)
    if '--rule' in argv:
        i = argv.index('--rule')
        r = argv[i + 1]
        del argv[i:i + 2]
        rules = [4, 0, -4] if r == 'all' else [int(r)]
        if any(q not in TAG for q in rules):
            sys.exit('rule must be 4, 0, -4 or all')
    grids = [int(a) for a in argv if a.isdigit()] or [64]
    if '--paper' in argv:
        rules, grids = [4, 0, -4], list(GRIDS)
    writeSurfaces()
    for rule in rules:
        for n in grids:
            print('written %s.inp' % writeInp(n, rule))


if __name__ == '__main__':
    main(sys.argv[1:])
