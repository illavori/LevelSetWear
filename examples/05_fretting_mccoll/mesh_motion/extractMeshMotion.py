# ======================================================================
# LevelSetWear, results of the fretting test of McColl et al.,
# mesh-motion model (example 05)
#
# Written by Inigo Llavori
#
# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
# ======================================================================
# Reads the output of the job MCKS_18000 (run.bat) and writes the
# mesh-motion data of Fig. 16 of the article, the worn surface of the
# flat after 1000, 5000, 10 000 and 18 000 cycles (every third surface
# node, as plotted), and of Fig. 17, the Newton iterations per increment
# of each wear block. It also prints the mesh-motion column of Table 3
# for this test.
#   ContactOut_slave.txt  one line per call of the subroutine for a node
#   ContactOut_mast.txt   of the flat (slave) or of the cylinder (master),
#       node, CPRESS, CSHEAR, MASTER, CSLIP, DSLIP, WDIST, x, y
#   initialFlatSurface.txt, initialCylinderSurface.txt  unworn surfaces
#   cycleJump, cylinderNodes.txt   cycles per wear block, cylinder nodes
#   MCKS_18000.sta    Newton iterations of each increment
#   MCKS_18000.dat    CPU time (or MCKS_18000_time.txt, its last lines)
# Step 1 is CONTACT, step 2 FRETTING_CYCLE (first cycle, no wear), and
# step s >= 3 is the wear step CYCLE_(s - 2), that ends cycle 100 (s - 2).
#
# The surfaces in the odb are loaded and deformed, and in ALE the
# displacement includes the mesh motion, so the wear is summed from the
# depths WDIST that the subroutine applies. The ContactOut files are
# split into wear blocks of 82 call blocks (41 increments, 2 calls per
# increment). Abaqus calls the subroutine twice per increment with the
# same values and the wear is applied once, so a line equal to the
# previous line of the same node is counted once. The wear of each node is
# placed at its mean x over the calls, and the wear depth is measured
# against the unworn surface interpolated there.
#   python extractMeshMotion.py             job run in this folder
#   python extractMeshMotion.py FOLDER      job folder elsewhere
#   python extractMeshMotion.py --figures
#       overwrites mesh_motion.csv in ../../../figures/fig16_fretting_mccoll
#       and ../../../figures/fig17_iterations, then python plot.py there
#       redraws them
# Without --figures the CSV files are written as fig16_mesh_motion.csv and
# fig17_mesh_motion.csv to the job folder, or to the folder given with
# --out DIR. Pure Python 3 with numpy.
import argparse
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
JOB = 'MCKS_18000'
FIRST = 3                       # first wear step (CYCLE_1)
NBLOCK = 82                     # call blocks per wear block (41 increments x 2 calls)
CYCLES = (1000, 5000, 10000, 18000)
CITE = """# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
"""
HEADER16 = ('# LevelSetWear, data of Fig. 16 of\n' + CITE +
            '# fretting test of McColl et al. (2004), 185 N, worn surface of the flat '
            '(wear depth, negative down)\n'
            '# mesh motion model (UMESHMOTION) on the same mesh, single analysis, '
            'every third surface node (as plotted)\n'
            'cycles,x_mm,y_um\n')
HEADER17 = ('# LevelSetWear, data of Fig. 17 of\n' + CITE +
            '# Newton iterations per converged increment along the fretting test of '
            'McColl et al. (2004), 185 N,\n'
            '# one point per wear block of 100 cycles (from the Abaqus .sta files)\n'
            '# mesh motion model (UMESHMOTION), single analysis\n'
            'cycles,iterations_per_increment\n')


def first(folder, name):
    """first value of a small input file"""
    return int(open(os.path.join(folder, name)).read().split()[0])


def blocks(folder, body):
    """[(cycle, {node: wear (mm)}, {node: [sum of x, calls]})], one entry
    per wear block, from the ContactOut file of the body"""
    dn = first(folder, 'cycleJump')
    off = first(folder, 'cylinderNodes.txt') if body == 'flat' else 0
    src = os.path.join(folder, 'ContactOut_%s.txt' % ('slave' if body == 'flat' else 'mast'))
    out, lab0, block = [], None, -1
    w, pos, last = {}, {}, {}
    for l in open(src):
        v = l.split()
        if len(v) < 9:
            continue
        lab = int(v[0])
        if lab0 is None:
            lab0 = lab
        if lab == lab0:                         # a new call block
            block += 1
            if block % NBLOCK == 0:             # a new wear block
                w, pos, last = {}, {}, {}
                out.append(((block//NBLOCK + 1)*dn, w, pos))
        n = lab - off
        if last.get(n) == l:                    # second call, same values
            continue
        last[n] = l
        w[n] = w.get(n, 0.0) + float(v[6])
        q = pos.setdefault(n, [0.0, 0])
        q[0] += float(v[7])
        q[1] += 1
    return out


def profile(folder, data, cycle, body):
    """x (mm), worn surface y (mm) of the body after `cycle` cycles"""
    ini = np.loadtxt(os.path.join(folder, 'initial%sSurface.txt'
                                  % ('Flat' if body == 'flat' else 'Cylinder')), delimiter=',')
    ini = ini[np.argsort(ini[:, 1])]
    tot, pos = {}, {}
    for j, w, p in data:
        if j <= cycle:
            for n, d in w.items():
                tot[n] = tot.get(n, 0.0) + d
            for n, q in p.items():
                r = pos.setdefault(n, [0.0, 0])
                r[0] += q[0]
                r[1] += q[1]
    if not pos:                                 # no wear block yet
        return ini[:, 1:3]
    sgn = -1.0 if body == 'flat' else 1.0
    x = np.array(sorted(q[0]/q[1] for q in pos.values()))
    n = sorted(pos, key=lambda k: pos[k][0]/pos[k][1])
    y0 = np.interp(x, ini[:, 1], ini[:, 2])
    return np.column_stack([x, y0 + sgn*np.array([tot.get(k, 0.0) for k in n])])


def wear(folder, data, cycle, body):
    """x (mm) and wear depth (um) of the body after `cycle` cycles"""
    u0 = profile(folder, data, 0, body)
    u = profile(folder, data, cycle, body)
    return np.column_stack([u[:, 0], 1000*np.abs(u[:, 1] - np.interp(u[:, 0], u0[:, 0], u0[:, 1]))])


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


def write(path, header, lines):
    with open(path, 'w', newline='\n') as f:
        f.write(header + ''.join(lines))
    print('written', os.path.normpath(path))


def main():
    ap = argparse.ArgumentParser(
        description='Mesh-motion results of Figs. 16 and 17 and Table 3 (fretting)')
    ap.add_argument('folder', nargs='?', default=HERE,
                    help='job folder (default this folder)')
    ap.add_argument('--out', help='folder for the CSV files (default the job folder)')
    ap.add_argument('--figures', action='store_true',
                    help='overwrite mesh_motion.csv in ../../../figures')
    opt = ap.parse_args()
    folder = opt.folder
    if not os.path.isdir(folder):
        folder = os.path.join(HERE, folder)
    folder = os.path.abspath(folder)
    out = os.path.abspath(opt.out) if opt.out else folder
    figs = opt.figures
    FIG = os.path.join(HERE, '..', '..', '..', 'figures')

    # Fig. 16, worn surface of the flat
    flat = blocks(folder, 'flat')
    lines = []
    for n in CYCLES:
        a = wear(folder, flat, n, 'flat')[::3]
        for x, y in zip(a[:, 0], -a[:, 1]):
            lines.append('%d,%r,%r\n' % (n, float(x), float(y)))
    write(os.path.join(FIG, 'fig16_fretting_mccoll', 'mesh_motion.csv') if figs
          else os.path.join(out, 'fig16_mesh_motion.csv'), HEADER16, lines)

    # Fig. 17, iterations per converged increment of each wear block
    path = os.path.join(folder, JOB + '.sta')
    r = sta(path) if os.path.isfile(path) else None
    dn = first(folder, 'cycleJump')
    if r:
        lines = ['%d,%r\n' % ((s - FIRST + 1)*dn, r[s][0]/max(r[s][1], 1))
                 for s in sorted(r) if s >= FIRST]
        write(os.path.join(FIG, 'fig17_iterations', 'mesh_motion.csv') if figs
              else os.path.join(out, 'fig17_mesh_motion.csv'), HEADER17, lines)
    else:
        print('no %s.sta in %s, Fig. 17 not written' % (JOB, folder))

    # Table 3
    last = flat[-1][0]
    print('mesh motion, fretting, %d cycles (%d wear blocks)' % (last, len(flat)))
    a = wear(folder, flat, last, 'flat')
    d, w, wn = scar(a[:, 0], a[:, 1])
    print('  depth  %.2f um' % d)
    print('  width  %.2f mm (1 %% of the depth, interpolated), %.2f mm between nodes'
          % (w, wn))
    o = np.argsort(a[:, 0])
    print('  worn area  %.2f x 1e-3 mm2' % np.trapz(a[o, 1], a[o, 0]))
    a = wear(folder, blocks(folder, 'cylinder'), last, 'cylinder')
    print('  cylinder depth  %.2f um' % a[:, 1].max())
    if r:
        it = sum(v[0] for s, v in r.items() if s >= FIRST)
        n = sum(v[1] for s, v in r.items() if s >= FIRST)
        print('  iterations per increment  %.2f (%d iterations, %d increments of the wear steps)'
              % (it/n, it, n))
    t = cpuTime(folder, JOB)
    if t is not None:
        print('  CPU time  %.1f h (%.0f s)' % (t/3600, t))


if __name__ == '__main__':
    main()
