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
# Abaqus Python script. Reads the odb of the level set model (aligned/)
# or of the standard Abaqus model (standard_abaqus/) and writes, next to
# the odb, the results that mindlinResults.py turns into the data of
# Figs. 11 and 12 and Table 1 of the article. Last frame of the steps
# TANGENTIAL (loading to Q*) and UNLOAD (unloading to Q = 0).
#   <job>_U_<STEP>.txt        elements of the first row of the flat below
#                             the contact (|x| < 0.2 mm, -0.01 < y < 0),
#                             "E label" and 4 lines "node x y ux uy"
#   <job>_S_<STEP>.txt        element-average stresses of the flat,
#                             "label xc yc s11 s22 s33 s12 frac". In the
#                             level set model these are the material
#                             averages of the cut elements (SDV1-5 of the
#                             ghost elements), with the UEL label.
#   <job>_contact_<STEP>.txt  standard model only, "x CPRESS CSHEAR1" of
#                             the surface nodes of the flat (y = 0)
# The kind of model is found from the odb (SDV_S11 present: level set).
# Run it after the job, from the folder of the job, for instance
#   cd aligned          & abaqus python ..\extractOdb.py uel_flat_MD.odb
#   cd standard_abaqus  & abaqus python ..\extractOdb.py std_flat_MD.odb
# (abaqus is the command of your Abaqus installation, e.g. abq2020).
import os
import sys

from odbAccess import openOdb

STEPS = ('TANGENTIAL', 'UNLOAD')
GHOST = 100000          # ghost element label = UEL element label + GHOST
NOFF_FLAT = 20000       # nodes of the flat are numbered from here
YMIN = -0.01            # first row of elements below the surface (mm)
XMAX = 0.2              # half width of the strip of the sigma_xx output (mm)

odbName = sys.argv[-1]
base = os.path.splitext(odbName)[0]
odb = openOdb(odbName, readOnly=True)
inst = odb.rootAssembly.instances['PART-1-1']
X = dict((n.label, n.coordinates) for n in inst.nodes)
C = {}
for e in inst.elements:
    xs = [X[n] for n in e.connectivity]
    C[e.label] = (sum(p[0] for p in xs)/len(xs), sum(p[1] for p in xs)/len(xs))
uel = 'SDV_S11' in odb.steps[STEPS[0]].frames[-1].fieldOutputs.keys()
print('%s: %s model' % (odbName, 'level set' if uel else 'standard'))

for stp in STEPS:
    fr = odb.steps[stp].frames[-1]

    # nodal displacements of the first row of elements of the flat
    U = dict((v.nodeLabel, v.data) for v in fr.fieldOutputs['U'].values)
    out = open('%s_U_%s.txt' % (base, stp), 'w')
    for e in inst.elements:
        # standard model: CPE4. Level set model: its ghost CPE4 elements
        # (label + GHOST) share nodes and connectivity with the UEL elements
        c = e.connectivity
        if len(c) != 4:
            continue
        yc = sum(X[n][1] for n in c)/4.0
        xc = sum(X[n][0] for n in c)/4.0
        if not (YMIN < yc < 0.0) or abs(xc) > XMAX:
            continue
        out.write('E %d\n' % e.label)
        for n in c:
            out.write('%d %.10e %.10e %.10e %.10e\n' % (n, X[n][0], X[n][1], U[n][0], U[n][1]))
    out.close()

    # element-average stresses
    if not uel:
        acc = {}
        for v in fr.fieldOutputs['S'].values:
            a = acc.setdefault(v.elementLabel, [0.0]*5)
            d = v.data
            a[0] += d[0]; a[1] += d[1]; a[2] += d[2]; a[3] += d[3]; a[4] += 1
        rows = [(l, a[0]/a[4], a[1]/a[4], a[2]/a[4], a[3]/a[4], 1.0) for l, a in acc.items()]
    else:
        keys = ('SDV_S11', 'SDV_S22', 'SDV_S33', 'SDV_S12', 'SDV_VFRAC')
        comp = {}
        for k in keys:
            for v in fr.fieldOutputs[k].values:
                a = comp.setdefault(v.elementLabel, {}).setdefault(k, [0.0, 0])
                a[0] += v.data; a[1] += 1
        rows = []
        for l, d in comp.items():
            m = [d[k][0]/d[k][1] for k in keys]
            rows.append((l - GHOST, m[0], m[1], m[2], m[3], m[4]))
    out = open('%s_S_%s.txt' % (base, stp), 'w')
    for r in rows:
        lab = r[0]
        xc, yc = C[lab + GHOST if uel else lab]
        out.write('%d %.8e %.8e %.8e %.8e %.8e %.8e %.6f\n' % ((lab, xc, yc) + tuple(r[1:])))
    out.close()

    # contact pressure and shear of the surface nodes of the flat (standard)
    if not uel:
        kp = [k for k in fr.fieldOutputs.keys() if k.startswith('CPRESS')][0]
        ks = [k for k in fr.fieldOutputs.keys() if k.startswith('CSHEAR1')][0]
        p = dict((v.nodeLabel, v.data) for v in fr.fieldOutputs[kp].values)
        q = dict((v.nodeLabel, v.data) for v in fr.fieldOutputs[ks].values)
        out = open('%s_contact_%s.txt' % (base, stp), 'w')
        for n in sorted(p, key=lambda n: X[n][0]):
            if n < NOFF_FLAT or abs(X[n][1]) > 1e-9:
                continue
            out.write('%.10e %.10e %.10e\n' % (X[n][0], p[n], q.get(n, 0.0)))
        out.close()
    print('  %s written' % stp)
odb.close()
