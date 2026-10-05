@echo off
rem LevelSetWear, written by Inigo Llavori. If you use it, please cite
rem   I. Llavori, Wear simulation on a fixed mesh with level sets
rem   and immersed frictional contact,
rem   submitted for publication. Mondragon Unibertsitatea.
rem   https://github.com/illavori/LevelSetWear       BSD 3-Clause
rem Mindlin-Cattaneo and Mindlin-Deresiewicz, surfaces on the nodes
rem Replace abq2020 by the command of your Abaqus installation.
call abq2020 job=uel_flat_MD user=..\..\..\src\uelWearLS.for cpus=1 interactive ask_delete=OFF
