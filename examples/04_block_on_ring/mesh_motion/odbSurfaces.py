# ======================================================================
# LevelSetWear, worn surfaces of the mesh-motion model (example 04)
#
# Written by Inigo Llavori
#
# I. Llavori, Wear simulation on a fixed mesh with level sets
# and immersed frictional contact,
# submitted for publication. Mondragon Unibertsitatea.
# https://github.com/illavori/LevelSetWear       BSD 3-Clause
# ======================================================================
# Abaqus Python script. Reads the odb of the job BOR20S_3600 (run.bat)
# and writes, for every wear step CYCLE_k, the coordinates of the surface
# nodes of the block (node set ADAPTN) and of the ring (ADAPTMN) at the
# last frame of the step, in the folder Results of the job folder
#   Results/wearFlatSurfaceCycle_<n>.txt       block, "label, x, y" (mm)
#   Results/wearCylinderSurfaceCycle_<n>.txt   ring
# with n = k times the value of cycleJump (wear block k). The surfaces
# are those of the loaded model (COORD), as in Fig. 15 of the article.
# Then python extractMeshMotion.py writes the data of the figure.
#   abq2020 python odbSurfaces.py               job run in this folder
#   abq2020 python odbSurfaces.py FOLDER        job folder elsewhere
# Abaqus fails with paths longer than 255 characters. Run it from a
# short folder if needed.
import os
import sys

from odbAccess import openOdb

JOB = 'BOR20S_3600'
args = [a for a in sys.argv[1:] if not a.lower().endswith('.py')]
folder = os.path.abspath(args[0] if args else os.path.dirname(os.path.abspath(__file__)))
jump = int(open(os.path.join(folder, 'cycleJump')).read().split()[0])
out = os.path.join(folder, 'Results')
if not os.path.isdir(out):
    os.makedirs(out)
odb = openOdb(os.path.join(folder, JOB + '.odb'), readOnly=True)
flat = odb.rootAssembly.nodeSets['ADAPTN']
cyl = odb.rootAssembly.nodeSets['ADAPTMN']
n = 0
for name in odb.steps.keys():
    if not name.upper().startswith('CYCLE_'):
        continue
    k = int(name.split('_')[1])
    coord = odb.steps[name].frames[-1].fieldOutputs['COORD']
    for region, tag in ((flat, 'Flat'), (cyl, 'Cylinder')):
        f = open(os.path.join(out, 'wear%sSurfaceCycle_%d.txt' % (tag, k*jump)), 'w')
        for v in coord.getSubset(region=region).values:
            f.write('%d,   %.10g,   %.10g\n' % (v.nodeLabel, v.data[0], v.data[1]))
        f.close()
    n += 1
odb.close()
print('%d wear steps written to %s' % (n, out))
