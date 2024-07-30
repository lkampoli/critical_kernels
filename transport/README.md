# CO2-6T_HARM

This code computes the transport coefficients in a CO2/O2/CO/O/C mixture
with a six-temperature (T, Tv1, Tv2, Tv3, TvO2, TvCO) model using an 
harmonic oscillator approximation.

The output is given by:
  - shear viscosity coefficient [1]
  - bulk viscosity coefficient  [1]
  - heat conductivity coefficient [1]
  - vibrational heat conductivity coefficients (CO2, symmetric, bending, asymmetric mode, O2 and CO) ([5] otional, otherwise only [1])
  - mass diffusion coefficients [25]
  - thermal diffusion coefficients [5]

## How to compile and run?

Go to ```src``` folder and execute:
```
	./clean.sh
	./compile.sh
```

On successful compilation the executable file ```co2_harm_6T.x```
is created. To run the computation launch:
```
	./co2_harm_6T.x
```

Note that in order to successfully compile the source is required:
 - a fortran compilar (gfortran, intel) and correnspondingly,
 - openblas or mkl library

Note also that typically gfortran is slower than ifort.

The simplest option is to use gfortran and openblas.
More info about how to install openblas can be found at:
```
	https://github.com/xianyi/OpenBLAS
```

On linux machines (ubuntu), this is achieved by:
```
	sudp apt-get install gfortran libopenblas-dev
```

## Post-processing

The output file (database) is defined in the source file ```co2-o2_transp.f90```
at the line 77:
```
open (2,file='shear_viscosity_lite.dat',status='unknown',action='write',position='append')
```

Note the following observations:
 - to keep the size of the database file manageable it is recommended to
   create one file for each transport properties (this may not apply for
   large storage and RAM machines) with meaningful names.
 - Any output is "appended" which means that consecutive runs will write
   on the same file. So, be sure to delete any previous files before any
   run or disable the positioning flag in the aforementioned line.
 - Compile as explained above after any modification to sources files,
   otherwise no change will be applied and the executables files will be
   the same. 

An example file has been generated 
```
	shear_viscosity_lite.dat
```
and renamed *.csv after substituting spaces with commas. This is not
strictly necessary but could be convenient if working with pandas.
To reduce the size of the file (~114 MB) on the repository, it has 
been compressed. Thus, it should be uncompressed by doing:
```
	tar -xvfz shear_viscosity_lite.dat.tar.gz
```

## TODO

Several modification to the ```co2-o2_transp.f90``` file should be done.

 1. reduce the precision of the molar factions, x when printing (2 or 3) [X]
 2. modify the loops of the temperatures to improve the sampling space   [X]
 3. possibly introduce the pressure dependency among input features      [ ]
