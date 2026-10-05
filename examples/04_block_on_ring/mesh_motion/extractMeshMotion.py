# ======================================================================
# LevelSetWear, results of the block-on-ring test, mesh-motion model
# (example 04)
#
# Written by Inigo Llavori
#
# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
# ======================================================================
# Reads the output of the job BOR20S_3600 (run.bat) and writes the
# mesh-motion data of Fig. 15 of the article, the worn surface of the
# block after 800, 1600, 2400 and 3600 m of sliding (every second surface
# node, as plotted). It also prints the mesh-motion column of Table 3 for
# this test.
#   Results/wearFlatSurfaceCycle_<n>.txt   surface of the block at the end
#       of wear block n, 80 m each (written from the odb by
#       abq2020 python odbSurfaces.py)
#   BOR20S_3600.sta    Newton iterations of each increment
#   BOR20S_3600.dat    CPU time (or BOR20S_3600_time.txt, its last lines)
# Step 1 is CONTACT, step 2 FRETTING_CYCLE (first block, no wear), and
# step s >= 3 is the wear step CYCLE_(s - 2).
#   python extractMeshMotion.py             job run in this folder
#   python extractMeshMotion.py FOLDER      job folder elsewhere
#   python extractMeshMotion.py --figures
#       overwrites ../../../figures/fig15_block_on_ring/mesh_motion.csv,
#       then python plot.py in that folder redraws the figure
# Without --figures the CSV file is written as fig15_mesh_motion.csv to the
# job folder, or to the folder given with --out DIR. Pure Python 3 with numpy.
import argparse
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
JOB = 'BOR20S_3600'
FIRST = 3                       # first wear step (CYCLE_1)
DS = 80                         # m of sliding per wear block
SLIDING = (800, 1600, 2400, 3600)
HEADER = """# LevelSetWear, data of Fig. 15 of
# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
# block-on-ring test of Cruzado et al. (2010), 65 N, worn surface of the block (wear depth, negative down)
# mesh motion model (UMESHMOTION) on the same mesh, one analysis, every second surface node (as plotted)
sliding_m,x_mm,y_um
"""


def surface(res, n):
    """x, y (mm) of the surface of the block after wear block n, sorted by x"""
    a = np.loadtxt(os.path.join(res, 'wearFlatSurfaceCycle_%d.txt' % n), delimiter=',')
    return a[np.argsort(a[:, 1])][:, 1:3]


def sta(path):
    """{step: (Newton iterations, converged increments)} of an Abaqus .sta"""
    it, n = {}, {}
    for l in open(path, errors='replace'):
        w = l.split()
        if len(w) > 8 and w[0].isdigit() and w[1].isdigit():
            s = int(w[0])
            it[s] = it.get(s, 0) + int(w[5])
            if 'U' not in w[2]:                 # 'nU' marks a cutback
                n[s] = n.get(s, 0) + 1
    return {s: (it[s], n.get(s, 0)) for s in it}


def cpuTime(folder, job):
    """last TOTAL CPU TIME (s) of the .dat file, or of <job>_time.txt"""
    for name in (job + '.dat', job + '_time.txt'):
        path = os.path.join(folder, name)
        if os.path.isfile(path):
            t = re.findall(r'TOTAL CPU TIME \(SEC\)\s*=\s*([0-9.Ee+-]+)', open(path).read())
            if t:
                return float(t[-1])
    return None


def scar(x, w, frac=0.01):
    """depth, and width of the zone where the wear w exceeds frac times the
    depth, interpolated between nodes and between the outermost nodes"""
    d = w.max()
    t = frac*d
    on = np.nonzero(w > t)[0]
    i, j = on[0], on[-1]
    xl = x[i] if i == 0 else x[i - 1] + (t - w[i - 1])*(x[i] - x[i - 1])/(w[i] - w[i - 1])
    xr = x[j] if j == len(x) - 1 else x[j] + (t - w[j])*(x[j + 1] - x[j])/(w[j + 1] - w[j])
    return d, xr - xl, x[j] - x[i]


def main():
    ap = argparse.ArgumentParser(
        description='Mesh-motion results of Fig. 15 and Table 3 (block-on-ring)')
    ap.add_argument('folder', nargs='?', default=HERE,
                    help='job folder (default this folder)')
    ap.add_argument('--out', help='folder for the CSV file (default the job folder)')
    ap.add_argument('--figures', action='store_true',
                    help='overwrite mesh_motion.csv in ../../../figures')
    opt = ap.parse_args()
    folder = opt.folder
    if not os.path.isdir(folder):
        folder = os.path.join(HERE, folder)
    folder = os.path.abspath(folder)
    out = os.path.abspath(opt.out) if opt.out else folder
    res = os.path.join(folder, 'Results')
    jump = int(open(os.path.join(folder, 'cycleJump')).read().split()[0])
    lines = []
    for s in SLIDING:
        for x, y in surface(res, s//DS*jump)[::2]:
            lines.append('%d,%r,%r\n' % (s, float(x), float(1000*y)))
    out = os.path.join(out, 'fig15_mesh_motion.csv')
    if opt.figures:
        out = os.path.join(HERE, '..', '..', '..', 'figures', 'fig15_block_on_ring',
                           'mesh_motion.csv')
    with open(out, 'w', newline='\n') as f:
        f.write(HEADER + ''.join(lines))
    print('written', os.path.normpath(out))

    last = max(int(re.search(r'_(\d+)\.txt$', f).group(1)) for f in os.listdir(res)
               if re.match(r'wearFlatSurfaceCycle_\d+\.txt$', f))
    a = surface(res, last)
    d, w, wn = scar(a[:, 0], np.maximum(-1000*a[:, 1], 0.0))
    print('mesh motion, block-on-ring, %d m (wear block %d)' % (last//jump*DS, last//jump))
    print('  depth  %.1f um' % d)
    w5 = scar(a[:, 0], np.maximum(-1000*a[:, 1], 0.0), 0.05)[1]
    print('  width  %.2f mm (5 %% of the depth, Table 3), %.2f mm (1 %%)' % (w5, w))
    o = np.argsort(a[:, 0])
    print('  worn area  %.3f mm2' % (1e-3*np.trapz(np.maximum(-1000*a[o, 1], 0.0), a[o, 0])))
    path = os.path.join(folder, JOB + '.sta')
    if os.path.isfile(path):
        r = sta(path)
        it = sum(v[0] for s, v in r.items() if s >= FIRST)
        n = sum(v[1] for s, v in r.items() if s >= FIRST)
        print('  iterations per increment  %.2f (%d iterations, %d increments of the wear steps)'
              % (it/n, it, n))
    t = cpuTime(folder, JOB)
    if t is not None:
        print('  CPU time  %.1f h (%.0f s)' % (t/3600, t))


if __name__ == '__main__':
    main()
