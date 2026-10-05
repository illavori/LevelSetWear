# LevelSetWear

Wear simulation on a fixed finite element mesh with level sets, implemented
as user elements (`UEL`) for Abaqus/Standard.

Conventional finite element wear models move the nodes of the mesh to follow
the worn surface, so every variable with history has to be transferred to the
new positions after each wear update. In LevelSetWear the mesh never moves.
The worn surfaces of the two bodies in contact are level sets, the elements
cut by them are integrated with fixed points whose weights follow the material
left in their cells, and the contact between the worn surfaces is enforced by
an immersed contact element with a dual mortar law. The integration points stay
with the material while the surfaces wear.

## Contents

| Folder | Content |
|---|---|
| `src/uelWearLS.for` | The user subroutines: cut bulk element U1, contact element U2, ghost `UMAT` for visualisation and `UEXTERNALDB` for the wear update |
| `src/umeshArchard2D.for` | The reference mesh-motion model (`UMESHMOTION`) used for the comparison in the article |
| `docs/LevelSetWear_Manual.pdf` | User manual, to use the code. Model, element properties, state variables, input files, usage and examples |
| `examples/` | The input files of every analysis of the article, with the scripts that extract their results and the results of the article in `reference/` |
| `figures/` | The data and the plotting script of every figure of the article |

Every result of the article can be reproduced from this repository. Each
example runs its analyses, its scripts extract the results into the CSV
files of `figures/`, and `python plot.py` in each figure folder draws the
figure. On the runs of the article the scripts give the CSV files of
`figures/` bit for bit.

## Examples

| Example | What it shows | Results of the article |
|---|---|---|
| `01_plate_with_hole` | Cut elements with the fixed cell rule against the solution of Kirsch. `makeKirsch.py` writes the 21 inputs of the article (7 grids, 3 integration rules), `extractKirsch.py` (Abaqus Python) and `postKirsch.py` extract the results, see its `README.txt` | Section 4.1, Fig. 9 |
| `02_contact_mindlin` | Immersed contact against the Mindlin–Cattaneo and Mindlin–Deresiewicz solutions, surfaces aligned with the mesh, immersed in it (`immersed/makeImmersed.py` writes the seven cases) and the standard Abaqus contact. `extractOdb.py` (Abaqus Python) and `mindlinResults.py` extract the results | Section 4.2, Figs. 11 and 12, Tables 1 and 2 |
| `03_block_on_ring` | Sliding wear, block-on-ring test of Cruzado et al. Level set model (`extractBOR.py`) and mesh-motion model in `mesh_motion/` (`odbSurfaces.py`, Abaqus Python, and `extractMeshMotion.py`) | Section 5.1, Fig. 15, Table 3 |
| `04_fretting_mccoll` | Fretting wear of both bodies, test of McColl et al. Level set model (`extractMcColl.py`) and mesh-motion model in `mesh_motion/` (`extractMeshMotion.py`) | Section 5.2, Figs. 16 and 17, Table 3 |
| `05_level_set_update` | Narrow band level set update against the exact height-function solution (Python, no Abaqus). `python paperLS.py` draws the figure of the article | Section 4.3, Fig. 13 |

The level set scripts of examples 03 and 04 also read the shipped results,
`python extractBOR.py reference`, so that the figures can be redrawn
without running Abaqus. Each script prints the numbers of the tables of the
article and has its usage in its header.

## Figures

Each folder in `figures/` reproduces one figure of the article from the data
in its CSV files. Run `python plot.py` in the folder (NumPy and Matplotlib).

| Folder | Figure of the article |
|---|---|
| `fig09_kirsch` | Plate with a hole on the fixed grid, stress, errors and convergence |
| `fig11_surface_mindlin` | Surface stresses in loading and unloading, analytical, standard contact and level set |
| `fig12_stress_mindlin` | Von Mises stress below the contact, standard and level set models |
| `fig13_level_set_update` | Level set update against the height function |
| `fig15_block_on_ring` | Block-on-ring test of Cruzado et al., worn profiles against the measurements |
| `fig16_fretting_mccoll` | Fretting test of McColl et al., worn profiles against the measurement |
| `fig17_iterations` | Newton iterations per increment in the fretting test |

## Quick start

Requirements: Abaqus/Standard with a Fortran compiler configured for user
subroutines (tested with Abaqus 2020 and Intel Fortran 2021) on Windows.
Example 05 only needs Python 3 with NumPy, SciPy and Matplotlib.

```bat
cd examples\02_contact_mindlin\aligned
run.bat
```

Each `run.bat` calls `abq2020`; replace it with the command of your Abaqus
installation. The subroutines share data between elements through `COMMON`
blocks, so every job must run with `cpus=1`. Only the equation solver can run
in parallel (`standard_parallel = SOLVER` in `abaqus_v6.env`), which gives
identical results but no speed-up for these 2D models. Three-dimensional models
would need the code to be reorganised so that the elements run truly in
parallel. Compare `contactOut.txt` or `wearProfile.txt` with the files
in `reference/`. The block-on-ring (about 18.5 h of CPU time for the level
set model and 25 h for the mesh-motion model) and fretting analyses (5.7 h
and 4.7 h) are long. Abaqus fails with paths longer than 255 characters, so
copy an example to a short folder if needed.

## Citation

If you use this code, please cite

> I. Llavori, Wear simulation on a fixed mesh with level sets and immersed frictional contact,
> submitted for publication.

## License

BSD 3-Clause, see [LICENSE](LICENSE).

## Contact

Iñigo Llavori, Department of Mechanical Engineering, Mondragon Unibertsitatea
(illavori@mondragon.edu).
