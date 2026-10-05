@echo off
rem LevelSetWear, written by Inigo Llavori. If you use it, please cite
rem   I. Llavori, Wear simulation on a fixed mesh with level sets
rem   and immersed frictional contact,
rem   submitted for publication. Mondragon Unibertsitatea.
rem   https://github.com/illavori/LevelSetWear       BSD 3-Clause
rem Fretting of McColl et al., 18000 cycles, level set model
rem Replace abq2020 by the command of your Abaqus installation.
call abq2020 job=MCKLS5_18000 user=..\..\src\uelWearLS.for cpus=1 interactive ask_delete=OFF
