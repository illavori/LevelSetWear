# ======================================================================
# LevelSetWear, results of the plate with a hole (example 01)
#
# Written by Inigo Llavori
#
# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
# ======================================================================
# Abaqus Python script (Python 2.7). Reads the odb of a plate with a hole
# job (makeKirsch.py) and writes, next to it, two text files for
# postKirsch.py
#   <job>_U.txt  nodal displacements of the last frame (strain energy)
#                label, U1, U2
#   <job>_S.txt  results of the cut elements, read from the ghost mesh
#                (grids up to N = 256 only). One line per element
#                label, x and y of the centroid, material-averaged S11 S22
#                S33 S12 (SDV1-4), material fraction (SDV5)
#   abq2020 python extractKirsch.py kirsch_64.odb [kirsch_p0_64.odb ...]
import os
import sys

from odbAccess import openOdb

GHOST = 100000                           # label offset of the ghost elements
SDV = ('SDV_S11', 'SDV_S22', 'SDV_S33', 'SDV_S12', 'SDV_VFRAC')


def extract(odbName):
    job = os.path.splitext(odbName)[0]
    odb = openOdb(odbName, readOnly=True)
    fr = odb.steps.values()[-1].frames[-1]
    f = open(job + '_U.txt', 'w')
    for v in fr.fieldOutputs['U'].values:
        f.write('%d %.9e %.9e\n' % (v.nodeLabel, v.data[0], v.data[1]))
    f.close()
    print('%s_U.txt' % job)
    if SDV[0] not in fr.fieldOutputs.keys():
        odb.close()
        return
    inst = odb.rootAssembly.instances['PART-1-1']
    X = dict((n.label, n.coordinates) for n in inst.nodes)
    C = {}
    for e in inst.elements:
        xs = [X[n] for n in e.connectivity]
        C[e.label] = (sum(p[0] for p in xs)/len(xs), sum(p[1] for p in xs)/len(xs))
    # average of the integration points of each ghost element (all equal)
    comp = {}
    for k in SDV:
        for v in fr.fieldOutputs[k].values:
            a = comp.setdefault(v.elementLabel, {}).setdefault(k, [0.0, 0])
            a[0] += v.data
            a[1] += 1
    f = open(job + '_S.txt', 'w')
    for lab, d in comp.items():
        m = [d[k][0]/d[k][1] for k in SDV]
        xc, yc = C[lab]
        f.write('%d %.8e %.8e %.8e %.8e %.8e %.8e %.6f\n'
                % ((lab - GHOST, xc, yc) + tuple(m)))
    f.close()
    print('%s_S.txt' % job)
    odb.close()


for name in sys.argv[1:]:
    if name.endswith('.odb'):
        extract(name)
