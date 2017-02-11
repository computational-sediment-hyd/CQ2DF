# CQ2DF - Cross section Quasi 2D Flow Model

====

CQ2DF is numerical simulation model that calculates simply flow velocity distribution in cross section of river.

## Description

 - CQ2DF is numerical simulation model that calculates the flow velocity distribution by Reynolds equation in x direction ignoring the differential term of x by assuming a uniform channel

 - Parameters are only river bed slope and equivalent roughness.

 - A Cartesian grid is used. 

## Demo

![example](/out.png "example")

## Requirement

- Windows 10
- Windows 7
- Python 3.5 or more
- gfortran 6.0 (Fortran 2003)


## Usage

- making DLL
```cmd
gfortran -std=f2003 -shared Cross-sectionQuasi2DFlow.f90 -static -o Cross-sectionQuasi2DFlow.dll
```
or
```cmd
build.bat
```

- execute
```cmd
python CQ2DF.py
```

## Licence

[MIT](/LICENCE)

## Author

[computational-sediment-hyd](https://github.com/computational-sediment-hyd)


