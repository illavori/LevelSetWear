LevelSetWear, example 01, plate with a hole (Kirsch)

Written by Inigo Llavori. If you use it, please cite
  I. Llavori, Wear simulation on a fixed mesh with level sets
  and immersed frictional contact,
  submitted for publication. Mondragon Unibertsitatea.
  https://github.com/illavori/LevelSetWear       BSD 3-Clause

Quarter plate with a circular hole under uniaxial tension on a fixed grid
cut by the level set phi = r0 - r (Section 4.1 and Fig. 9 of the article).
run.bat runs the 64 x 64 grid with the cell rule. reference/kirsch.txt
holds the strain energy of every grid with the cell rule.

Fig. 9 from scratch (replace abq2020 by the command of your Abaqus)

1. Input files of the 21 runs (grids N = 16 to 1024, three integration
   rules of the cut elements, cell rule, material polygon and point rule)
     python makeKirsch.py --paper
   Single grids and rules are also possible, for instance
     python makeKirsch.py 128 --rule 0      -> kirsch_p0_128.inp
   Grids up to N = 256 carry the ghost mesh with the stresses of the cut
   elements.

2. Each job with the UEL, one at a time
     abq2020 job=kirsch_p0_128 user=uelWearBulk.for interactive

3. The results of each odb (Abaqus Python)
     abq2020 python extractKirsch.py kirsch_16.odb kirsch_32.odb ...
   writes <job>_U.txt (displacements) and, with the ghost mesh,
   <job>_S.txt (stresses and material fraction of the cut elements).

4. The CSV files of figures/fig09_kirsch
     python postKirsch.py
   It also writes kirsch.txt, kirsch_p0.txt and kirsch_pt.txt (strain
   energy per grid, compare kirsch.txt with reference/kirsch.txt).

5. The figure
     cd ../../figures/fig09_kirsch
     python plot.py      -> kirsch.pdf / kirsch.png

The inputs of step 1 are the ones of the article. On the runs of the
article, steps 3 to 5 give the CSV files and the figure of the repository
bit for bit.
