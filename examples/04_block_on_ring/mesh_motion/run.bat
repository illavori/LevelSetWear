@echo off
rem LevelSetWear, written by Inigo Llavori. If you use it, please cite
rem   I. Llavori, Wear simulation on a fixed mesh with level sets
rem   and immersed frictional contact,
rem   submitted for publication. Mondragon Unibertsitatea.
rem   https://github.com/illavori/LevelSetWear       BSD 3-Clause
rem Block-on-ring, 3600 m in 45 wear blocks, mesh-motion model (UMESHMOTION)
rem Replace abq2020 by the command of your Abaqus installation.
call abq2020 job=BOR20S_3600 user=..\..\..\src\umeshArchard2D.for cpus=1 interactive ask_delete=OFF
rem Results, Fig. 15 and Table 3 (see extractMeshMotion.py)
rem   abq2020 python odbSurfaces.py
rem   python extractMeshMotion.py
