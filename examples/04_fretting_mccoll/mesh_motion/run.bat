@echo off
rem LevelSetWear, written by Inigo Llavori. If you use it, please cite
rem   I. Llavori, Wear simulation on a fixed mesh with level sets
rem   and immersed frictional contact,
rem   submitted for publication. Mondragon Unibertsitatea.
rem   https://github.com/illavori/LevelSetWear       BSD 3-Clause
rem Fretting of McColl et al., 18000 cycles in 180 wear blocks of 100 cycles, mesh-motion model (UMESHMOTION)
rem Replace abq2020 by the command of your Abaqus installation.
call abq2020 job=MCKS_18000 user=..\..\..\src\umeshArchard2D.for cpus=1 interactive ask_delete=OFF
rem Results, Figs. 16 and 17 and Table 3 (see extractMeshMotion.py)
rem   python extractMeshMotion.py
