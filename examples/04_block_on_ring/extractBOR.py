# ======================================================================
# LevelSetWear, results of the block-on-ring test (example 04)
#
# Written by Inigo Llavori
#
# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
# ======================================================================
# Reads the output of the job BOR20LS_3600 (run.bat) and writes the level
# set data of Fig. 15 of the article, the worn surface of the block after
# 800, 1600, 2400 and 3600 m of sliding. It also prints the level set
# column of Table 3 for this test.
#   wearProfile.txt     worn surface of each body at the end of each step
#   BOR20LS_3600.sta    Newton iterations of each increment
#   BOR20LS_3600.dat    CPU time (or BOR20LS_3600_time.txt, its last lines)
# Step 1 is APPROACH, step 2 NORMAL, and step s >= 3 ends the wear block
# s - 2, of 80 m each (45 blocks, 3600 m).
#   python extractBOR.py               job run in this folder
#   python extractBOR.py reference     output shipped with the repository
#   python extractBOR.py reference --figures
#       overwrites ../../figures/fig15_block_on_ring/level_set.csv, then
#       python plot.py in that folder redraws the figure
# Without --figures the CSV file is written as fig15_level_set.csv to the
# job folder, or to the folder given with --out DIR. Pure Python 3 with numpy.
import argparse
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
JOB = 'BOR20LS_3600'
FIRST = 3                       # first wear step (CYCLE_1)
DS = 80                         # m of sliding per wear block
SLIDING = (800, 1600, 2400, 3600)
HEADER = """# LevelSetWear, data of Fig. 15 of
# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
# block-on-ring test of Cruzado et al. (2010), 65 N, worn surface of the block (wear depth, negative down)
# level set model (UEL v5), element size l_e = 20 um
sliding_m,x_mm,y_um
"""


def profiles(path):
    """{step: {body: array of x, y (mm), sorted by x}} of wearProfile.txt"""
    out, st, b = {}, None, None
    for l in open(path):
        w = l.split()
        if not w:
            continue
        if w[0] == 'STEP':
            st = int(w[1])
            out[st] = {}
        elif w[0] == 'BODY':
            b = int(w[1])
            out[st][b] = []
        else:
            out[st][b].append([float(w[0]), float(w[1])])
    for d in out.values():
        for k in d:
            a = np.array(d[k])
            d[k] = a[np.argsort(a[:, 0], kind='stable')]
    return out


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



def measuredBOR():
    """Measured column of Table 3, block-on-ring: depth, width (5 % of the
    depth, between the digitised points) and worn area of the three tests,
    the area measured from the base line of each profile outside the scar
    (mean of |x| > 1.6 mm on each side)."""
    p = os.path.join(HERE, '..', '..', 'figures', 'fig15_block_on_ring',
                     'measured_cruzado2010.csv')
    lines = [l for l in open(p, encoding='utf-8') if not l.startswith('#')]
    d = np.loadtxt(lines[1:], delimiter=',', ndmin=2)
    x = d[:, 0]
    print('measured, block-on-ring (Cruzado et al.)')
    for k in (1, 2, 3):
        w = -d[:, k]
        on = np.nonzero(w > 0.05*w.max())[0]
        base = np.interp(x, [x.min(), x.max()], [w[x < -1.6].mean(), w[x > 1.6].mean()])
        print('  test %d  depth %.1f um, width %.2f mm, worn area %.3f mm2'
              % (k, w.max(), x[on[-1]] - x[on[0]], 1e-3*np.trapz(w - base, x)))

def main():
    ap = argparse.ArgumentParser(
        description='Level set results of Fig. 15 and Table 3 (block-on-ring)')
    ap.add_argument('folder', nargs='?', default=HERE,
                    help='job folder (default this folder, or reference)')
    ap.add_argument('--out', help='folder for the CSV files (default the job folder)')
    ap.add_argument('--figures', action='store_true',
                    help='overwrite level_set.csv in ../../figures')
    opt = ap.parse_args()
    folder = opt.folder
    if not os.path.isdir(folder):
        folder = os.path.join(HERE, folder)
    folder = os.path.abspath(folder)
    out = os.path.abspath(opt.out) if opt.out else folder
    P = profiles(os.path.join(folder, 'wearProfile.txt'))
    lines = []
    for s in SLIDING:
        for x, y in P[s//DS + FIRST - 1][1]:          # body 1 is the block
            lines.append('%d,%r,%r\n' % (s, x, 1000*y))
    out = os.path.join(out, 'fig15_level_set.csv')
    if opt.figures:
        out = os.path.join(HERE, '..', '..', 'figures', 'fig15_block_on_ring', 'level_set.csv')
    with open(out, 'w', newline='\n') as f:
        f.write(HEADER + ''.join(lines))
    print('written', os.path.normpath(out))

    last = max(P)
    a = P[last][1]
    d, w, wn = scar(a[:, 0], np.maximum(-1000*a[:, 1], 0.0))
    print('level set, block-on-ring, %d m (step %d)' % ((last - FIRST + 1)*DS, last))
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
    measuredBOR()


if __name__ == '__main__':
    main()
