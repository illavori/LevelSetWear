# ======================================================================
# LevelSetWear, results of the fretting test of McColl et al. (example 04)
#
# Written by Inigo Llavori
#
# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
# ======================================================================
# Reads the output of the job MCKLS5_18000 (run.bat) and writes the level
# set data of Fig. 16 of the article, the worn surface of the flat after
# 1000, 5000, 10 000 and 18 000 cycles, and of Fig. 17, the Newton
# iterations per increment of each wear block. It also prints the level
# set column of Table 3 for this test.
#   wearProfile.txt     worn surface of each body at the end of each step
#   MCKLS5_18000.sta    Newton iterations of each increment
#   MCKLS5_18000.dat    CPU time (or MCKLS5_18000_time.txt, its last lines)
# Step 1 is APPROACH, step 2 NORMAL, step 3 CYCLE_0 (no wear), and step
# s >= 4 ends the wear block s - 3, that is cycle 100 (s - 3). The wear is
# measured against the surface of step 1. Body 1 is the flat and body 2
# the cylinder.
#   python extractMcColl.py               job run in this folder
#   python extractMcColl.py reference     output shipped with the repository
#   python extractMcColl.py reference --figures
#       overwrites level_set.csv in ../../figures/fig16_fretting_mccoll and
#       ../../figures/fig17_iterations, then python plot.py there redraws them
# Without --figures the CSV files are written as fig16_level_set.csv and
# fig17_level_set.csv to the job folder, or to the folder given with
# --out DIR. Pure Python 3 with numpy.
import argparse
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
JOB = 'MCKLS5_18000'
FIRST = 4                       # first wear step (CYCLE_1)
DN = 100                        # cycles per wear block
CYCLES = (1000, 5000, 10000, 18000)
CITE = """# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
"""
HEADER16 = ('# LevelSetWear, data of Fig. 16 of\n' + CITE +
            '# fretting test of McColl et al. (2004), 185 N, worn surface of the flat '
            '(wear depth, negative down)\n'
            '# level set model (UEL v5, dual mortar contact)\n'
            'cycles,x_mm,y_um\n')
HEADER17 = ('# LevelSetWear, data of Fig. 17 of\n' + CITE +
            '# Newton iterations per converged increment along the fretting test of '
            'McColl et al. (2004), 185 N,\n'
            '# one point per wear block of 100 cycles (from the Abaqus .sta files)\n'
            '# level set model (UEL v5, dual mortar contact)\n'
            'cycles,iterations_per_increment\n')


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


def wear(xy, ref):
    """wear depth (um) against the unworn surface ref, on the points of xy"""
    return 1000*np.abs(np.interp(xy[:, 0], ref[:, 0], ref[:, 1]) - xy[:, 1])


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



def measuredMcColl():
    """Measured worn area of the flat (Table 3), from the base line of the
    profile outside the scar (mean of |x| > 0.3 mm on each side). The
    measured depth and width of Table 3 are those of McColl et al."""
    p = os.path.join(HERE, '..', '..', 'figures', 'fig16_fretting_mccoll',
                     'measured_mccoll2004.csv')
    lines = [l for l in open(p, encoding='utf-8') if not l.startswith('#')]
    d = np.loadtxt(lines[1:], delimiter=',', ndmin=2)
    x, w = d[:, 0], -d[:, 1]
    base = np.interp(x, [x.min(), x.max()], [w[x < -0.3].mean(), w[x > 0.3].mean()])
    print('measured, fretting (McColl et al.)  worn area %.2f x 1e-3 mm2'
          % np.trapz(w - base, x))

def main():
    ap = argparse.ArgumentParser(
        description='Level set results of Figs. 16 and 17 and Table 3 (fretting)')
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
    figs = opt.figures
    FIG = os.path.join(HERE, '..', '..', 'figures')

    # Fig. 16, worn surface of the flat
    P = profiles(os.path.join(folder, 'wearProfile.txt'))
    lines = []
    for n in CYCLES:
        a = P[FIRST - 1 + n//DN][1]
        for x, y in zip(a[:, 0], 0.0 - wear(a, P[1][1])):
            lines.append('%d,%r,%r\n' % (n, x, y))
    write(os.path.join(FIG, 'fig16_fretting_mccoll', 'level_set.csv') if figs
          else os.path.join(out, 'fig16_level_set.csv'), HEADER16, lines)

    # Fig. 17, iterations per converged increment of each wear block
    path = os.path.join(folder, JOB + '.sta')
    r = sta(path) if os.path.isfile(path) else None
    if r:
        lines = ['%d,%r\n' % ((s - FIRST + 1)*DN, r[s][0]/max(r[s][1], 1))
                 for s in sorted(r) if s >= FIRST]
        write(os.path.join(FIG, 'fig17_iterations', 'level_set.csv') if figs
              else os.path.join(out, 'fig17_level_set.csv'), HEADER17, lines)
    else:
        print('no %s.sta in %s, Fig. 17 not written' % (JOB, folder))

    # Table 3
    last = max(P)
    print('level set, fretting, %d cycles (step %d)' % ((last - FIRST + 1)*DN, last))
    a = P[last][1]
    d, w, wn = scar(a[:, 0], wear(a, P[1][1]))
    print('  depth  %.2f um' % d)
    print('  width  %.2f mm (1 %% of the depth, interpolated), %.2f mm between nodes'
          % (w, wn))
    ww = wear(a, P[1][1])
    o = np.argsort(a[:, 0])
    print('  worn area  %.2f x 1e-3 mm2' % np.trapz(ww[o], a[o, 0]))
    a = P[last][2]
    print('  cylinder depth  %.2f um' % wear(a, P[1][2]).max())
    if r:
        it = sum(v[0] for s, v in r.items() if s >= FIRST)
        n = sum(v[1] for s, v in r.items() if s >= FIRST)
        print('  iterations per increment  %.2f (%d iterations, %d increments of the wear steps)'
              % (it/n, it, n))
    t = cpuTime(folder, JOB)
    if t is not None:
        print('  CPU time  %.1f h (%.0f s)' % (t/3600, t))
    measuredMcColl()


if __name__ == '__main__':
    main()
