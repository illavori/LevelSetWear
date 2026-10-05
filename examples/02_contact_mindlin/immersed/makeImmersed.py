# ======================================================================
# LevelSetWear, input files of the immersed Mindlin contact (example 02)
#
# Written by Inigo Llavori
#
# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
# ======================================================================
# The immersed cases of Table 2 of the article. The surface of the flat
# lies at a depth dyF inside the flat and the surface of the cylinder at a
# depth dyC inside the cylinder (fractions of the surface element size
# le = 0.01 mm), so that the surfaces cut the elements instead of lying on
# the nodes. The mesh does not change: the nodes of the cylinder are
# lowered by (dyF + dyC) le, so that both true surfaces still touch at
# x = 0, and the surface tables surf1.txt (flat) and surf2.txt (cylinder)
# are lowered by dyF le.
#   case    dyF/le  dyC/le
#   F01     0.01    0       nearly on the nodes, 99 % of material in the cut row
#   F25     0.25    0
#   F50     0.50    0
#   F75     0.75    0
#   F99     0.99    0       1 % of material in the cut row
#   F150    1.50    0       second row, the first row of the flat is void
#   F50C50  0.50    0.50    both surfaces at mid-height (the inputs of this folder)
# The inputs are built from those of this example, aligned/uel_flat_MD.inp
# (nodes on the surfaces) and uel_imm_F50C50.inp (immersed contact options),
# and are identical to the inputs of the article.
#   python makeImmersed.py F25        -> F25/uel_imm_F25.inp, surf1.txt,
#                                        surf2.txt, v2opts.txt, run.bat
#   python makeImmersed.py all        -> the seven cases
# Then run.bat in the folder of each case, and python ..\mindlinResults.py
# for the numbers of Table 2.
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CASES = {'F01': (0.01, 0.0), 'F25': (0.25, 0.0), 'F50': (0.50, 0.0), 'F75': (0.75, 0.0),
         'F99': (0.99, 0.0), 'F150': (1.50, 0.0), 'F50C50': (0.50, 0.50)}
HSURF = 0.01            # surface element size le (mm)
NSUB = 8                # table points per surface element
R1 = 6.0                # cylinder radius (mm)
NOFF_FLAT = 20000       # nodes of the flat are numbered from here
CITE = ('LevelSetWear, written by Inigo Llavori. Lines after the data are not read.\n'
        'If you use it, please cite I. Llavori, Wear simulation on a fixed mesh with level sets\n'
        'and immersed frictional contact, submitted for publication. Mondragon Unibertsitatea.\n'
        'https://github.com/illavori/LevelSetWear       BSD 3-Clause\n')
RUN = ('@echo off\n'
       'rem LevelSetWear, written by Inigo Llavori. If you use it, please cite\n'
       'rem   I. Llavori, Wear simulation on a fixed mesh with level sets\n'
       'rem   and immersed frictional contact,\n'
       'rem   submitted for publication. Mondragon Unibertsitatea.\n'
       'rem   https://github.com/illavori/LevelSetWear       BSD 3-Clause\n'
       'rem Mindlin contact, {title}\n'
       'rem Replace abq2020 by the command of your Abaqus installation.\n'
       'call abq2020 job={job} user=..\\..\\..\\..\\src\\uelWearLS.for cpus=1 interactive '
       'ask_delete=OFF\n')


def title(dF, dC):
    if dF == dC == 0.5:
        return 'surfaces immersed at half an element'
    if dC > 0:
        return 'both surfaces immersed in their first row'
    if dF > 1:
        return 'surface of the flat immersed at %.2f le, in its second row' % dF
    return 'surface of the flat immersed at %.2f le in its first row' % dF


def nodeBlock(lines):
    """Index range of the node lines after *NODE, NSET=NALL."""
    i0 = lines.index('*NODE, NSET=NALL') + 1
    i1 = i0
    while not lines[i1].startswith('*'):
        i1 += 1
    return i0, i1


def readTable(name):
    """Vertices of a surface table (every NSUB-th point)."""
    lines = open(os.path.join(HERE, '..', 'aligned', name)).read().splitlines()
    n = int(lines[0].split()[0])
    pts = [tuple(float(v) for v in l.split()) for l in lines[1:n + 1]]
    return pts[::NSUB]


def surfTable(pts, yfun):
    """Surface table through the surface nodes (exactly) with NSUB-1 extra
    points per interval on the true surface yfun."""
    t = []
    for (x0, y0), (x1, _) in zip(pts[:-1], pts[1:]):
        t.append((x0, y0))
        t += [(x0 + (x1 - x0)*k/NSUB, yfun(x0 + (x1 - x0)*k/NSUB)) for k in range(1, NSUB)]
    t.append(pts[-1])
    return t


def write(path, lines):
    with open(path, 'w', newline='\r\n') as f:
        f.write('\n'.join(lines) + '\n')


def make(case):
    dF, dC = CASES[case][0]*HSURF, CASES[case][1]*HSURF
    out = os.path.join(HERE, case)
    os.makedirs(out, exist_ok=True)
    job = 'uel_imm_' + case
    # nodes of the cylinder on the surfaces (aligned input), lowered by dF + dC
    al = open(os.path.join(HERE, '..', 'aligned', 'uel_flat_MD.inp')).read().splitlines()
    i0, i1 = nodeBlock(al)
    cyl = {}
    for l in al[i0:i1]:
        n, x, y = l.split(',')
        if int(n) < NOFF_FLAT:
            cyl[int(n)] = (float(x), float(y))
    # the immersed input of this folder, with the new nodes and title
    L = open(os.path.join(HERE, 'uel_imm_F50C50.inp')).read().splitlines()
    i0, i1 = nodeBlock(L)
    for i in range(i0, i1):
        n = int(L[i].split(',')[0])
        if n < NOFF_FLAT:
            x, y = cyl[n]
            L[i] = '%d, %.15g, %.15g' % (n, x, y - 0.0 - dF - dC)
    k = L.index('*HEADING') + 1
    L[k] = 'Mindlin contact, cylinder on flat, %s, immersed contact element' % title(*CASES[case])
    write(os.path.join(out, job + '.inp'), L)
    # surface tables: flat at y = -dF, cylinder through its lowered nodes
    # plus dC, on the circle lowered by dF
    flat = [(x, y - dF) for x, y in readTable('surf1.txt')]
    cylt = [(x, y - 0.0 - dF - dC + dC) for x, y in readTable('surf2.txt')]
    tables = {1: surfTable(flat, lambda x: -dF),
              2: surfTable(cylt, lambda x: R1 - 0.0 - dF - math.sqrt(R1**2 - x**2))}
    for ib, pts in tables.items():
        write(os.path.join(out, 'surf%d.txt' % ib),
              ['%d %d' % (len(pts), NSUB)] + ['%.15g %.15g' % p for p in pts]
              + CITE.splitlines())
    with open(os.path.join(HERE, 'v2opts.txt')) as f, \
            open(os.path.join(out, 'v2opts.txt'), 'w', newline='\r\n') as g:
        g.write(f.read())
    with open(os.path.join(out, 'run.bat'), 'w', newline='\r\n') as f:
        f.write(RUN.format(title=title(*CASES[case]), job=job))
    print('%-7s dyF = %.2f le, dyC = %.2f le  -> %s' % (case, CASES[case][0], CASES[case][1],
                                                         os.path.relpath(out, HERE)))


if __name__ == '__main__':
    args = sys.argv[1:]
    if not args or any(a not in CASES and a != 'all' for a in args):
        print('python makeImmersed.py CASE [CASE ...] | all,  CASE in', ', '.join(CASES))
        sys.exit(1)
    for case in (CASES if 'all' in args else args):
        make(case)
