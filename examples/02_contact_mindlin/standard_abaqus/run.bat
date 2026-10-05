@echo off
rem LevelSetWear, written by Inigo Llavori. If you use it, please cite
rem   I. Llavori, Wear simulation on a fixed mesh with level sets
rem   and immersed frictional contact,
rem   submitted for publication. Mondragon Unibertsitatea.
rem   https://github.com/illavori/LevelSetWear       BSD 3-Clause
rem Reference model with the standard contact of Abaqus
call abq2020 job=std_flat_MD cpus=1 interactive ask_delete=OFF
