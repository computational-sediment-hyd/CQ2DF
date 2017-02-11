@echo off

gfortran -std=f2003 -shared ^
Cross-sectionQuasi2DFlow.f90 ^
-static -o Cross-sectionQuasi2DFlow.dll 

del *.mod
