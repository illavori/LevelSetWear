@echo off
rem LevelSetWear, written by Inigo Llavori. If you use it, please cite
rem   I. Llavori, Wear simulation on a fixed mesh with level sets
rem   and immersed frictional contact,
rem   submitted for publication. Mondragon Unibertsitatea.
rem   https://github.com/illavori/LevelSetWear       BSD 3-Clause
rem Plate with a hole, 64 x 64 grid, cut elements with the cell rule
rem Replace abq2020 by the command of your Abaqus installation.
call abq2020 job=kirsch_64 user=uelWearBulk.for cpus=1 interactive ask_delete=OFF
