# Copyright (C) 2017 computational-sediment-hyd https://github.com/computational-sediment-hyd
# Released under the MIT license
# http://opensource.org/licenses/mit-license.php

## @package CQ2DF.py
#  Documentation for this module.
# coding: utf-8
# using python 3.5

import numpy as np
from ctypes import *

## Documentation for a function.
#
#  read section data : tab separated
def inputTXT(inputfile):
    with open(inputfile, "r") as f:
        try:
            next(f) # skip-header
            y, z=[], []
            while True:
                data = next(f).split()  # read section header
                y.append(float(data[0]))
                z.append(float(data[1]))
        except StopIteration:
            pass   
    return y, z

## Documentation for a function.
#
# split section
def search(y, z, yval):
    imin=0
    imax=len(y)
    imid = int(np.round(0.5*(imin+imax)))
    
    while (imax-imin) > 1 :        
        if y[imid] < yval :
            imin = imid
        else:
            imax = imid
        
        imid = int(np.round(0.5*(imin+imax)))            
    zval = z[imin]+(yval-y[imin])*(z[imax]-z[imin])/(y[imax]-y[imin])
    return zval

## Documentation for a function.
#
# split section
def split(y,z, dy):
    nn = np.floor((y[-1] - y[0])/dy)
    ysplit = np.arange(int(nn))
    ysplit = y[0] + ysplit * dy + dy / 2
    zsplit = np.array([search(y,z,v) for v in ysplit])
    
    return ysplit, zsplit


## Documentation for a function.
#
# set cal Area
def setCalArea(y, z, dy, dz, waterLevel):

    ysplit, zsplit = split(y, z, dy)
    
    j=0
    while waterLevel-zsplit[j] < dz :
        j += 1
    jstart = j

    j=len(ysplit)-1
    while waterLevel-zsplit[j] < dz :
        j -= 1
    jend = j

    ny = jend - jstart + 1
    nz = int(np.ceil((waterLevel - min(zsplit[jstart:jend])) / dz))

    yOrg = ysplit[jstart] - dy / 2
    zOrg = waterLevel - nz * dz

    yOrgC = yOrg + dy / 2
    zOrgC = zOrg + dz / 2
    
#set Volume of Fluid
    VoF = np.array([[ 0 for j in range(ny)] for k in range(nz)])

    j, k = 0, 0
    for jj in range(jstart, jend+1):
        kk = 0
        while zsplit[jj] > (zOrgC + dz*kk) :
            kk += 1            
        VoF[kk:nz, j] = 1
        j +=1

# set Zb
    zb = np.array([0.0 for j in range(ny)])
    zb[0:ny] = zsplit[jstart:jend+1]

    return yOrg, zOrg, ny, nz, VoF, zb


## Documentation for a function.
#
# create cell condition
#  - if fluid CellCnd = 0001 elseif out of Area CellCnd = 0000
#  - if left next to cell   == bank CellCnd=1000
#  - if bottom next to cell == bank CellCnd=0100
#  - if right next to cell  == bank CellCnd=0010
# example 
# Cell Codition Fluid & left and bottom next to cell == bank CellCnd = 1101 
def setCellCnd(VoF,ny,nz):
    cellCnd = np.array([[ 0 for j in range(ny)] for k in range(nz)])

    for j in range(ny):
        for k in range(nz):
            if VoF[k,j] == 1:
                cellCnd[k,j] = 1
            #left-side
                if j == 0 :
                    cellCnd[k,j] += 1000
                else :
                    if VoF[k,j-1] == 0 : cellCnd[k,j] += 1000
            #right-side
                if j == ny-1 :
                    cellCnd[k,j] += 10
                else :
                    if VoF[k,j+1] == 0 : cellCnd[k,j] += 10
            #bottom
                if k == 0 :
                    cellCnd[k,j] += 100
                else :
                    if VoF[k-1,j] == 0 : cellCnd[k,j] += 100            

            else :
                cellCnd[k,j] == 0
    return cellCnd

## Documentation for a function.
#
#  Load DLL
def LoadDLL():
    fort = np.ctypeslib.load_library("Cross-sectionQuasi2DFlow.dll",".")

    fort.initilize_.argtypes = [
                            POINTER(c_double)
                            ,POINTER(c_double)
                            ,POINTER(c_double)
                            ,POINTER(c_int)
                            ,POINTER(c_int)
                            ,POINTER(c_double)
                            ,POINTER(c_double)
                            ]

    fort.calsectionprofile_.argtypes = [np.ctypeslib.ndpointer(dtype=np.int)
                                    , np.ctypeslib.ndpointer(dtype=np.float64)]


    return fort


# Main

# Set Calculation Condition with python
y, z = inputTXT('section.txt')

f = open('condition.txt', "r")
dy = float(f.readline().split()[1])
dz = float(f.readline().split()[1])
waterLevel = float(f.readline().split()[1])
ib = float(f.readline().split()[1])
ks = float(f.readline().split()[1])
kappa = float(f.readline().split()[1])
f.close

# Set Calculation Area
yOrg, zOrg, ny, nz, VoF, zb = setCalArea(y, z, dy, dz, waterLevel)

# Set Cell Condition
cellCnd = setCellCnd(VoF,ny,nz)

print('Origin(y,z)=',yOrg,',',zOrg)
print('Number of mesh(y,z)=',ny,',',nz)

## Calculation with Fortran
ib_, ks_, kappa_ = c_double(ib),c_double(ks),c_double(kappa)
ny_, nz_, dy_, dz_ = c_int(ny),c_int(nz),c_double(dy),c_double(dz)

fort = LoadDLL()

# Constructor
fort.initilize_(
            byref(ib_),byref(ks_),byref(kappa_)
            ,byref(ny_),byref(nz_),byref(dy_),byref(dz_)
            )
  
# Calculation
un = np.array([[ 0.0 for j in range(ny)] for k in range(nz)], dtype=np.float64)
cellCnd_ = cellCnd.astype(np.int)
fort.calsectionprofile_(cellCnd_, un)

# Destructor
fort.finalize_()

print('Calculation is finished')

# make gragh
print('making gragh figure....')

#import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
#import seaborn as sns

#sns.set(context='notebook', font_scale=1.5, style="whitegrid", palette="deep")

fig, ax = plt.figure(figsize=(8, 3),dpi=300), plt.subplot()
ax.plot(y,z,'k-')
ax.plot([y[0], y[-1]],[waterLevel, waterLevel],'b-')

yp, zp = np.arange(ny), np.arange(nz)
yp, zp = yOrg+yp*dy+0.5*dy, zOrg+zp*dz+0.5*dz

X, Y = np.meshgrid(yp, zp)

# masked-9999 
unmask = np.ma.masked_values(un, -9999.)

#norm = mpl.colors.Normalize(vmin=0.0, vmax=3.0)
fig = ax.contourf(X, Y, unmask, 15, cmap = cm.jet)
plt.colorbar(fig)

#ax.set_xlim(20,120)
#ax.set_ylim(4., 8)

plt.tight_layout()
plt.savefig('out.svg', bbox_inches='tight')

print('End')

