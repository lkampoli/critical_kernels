# change here as well ...
module purge 
module load PrgEnv-cray/8.4.0
module load openblas/0.3.24

rm *.o *.mod gmon.out *.x

FC=ifx
FC=gfortran
FC=ftn

$FC -c CONSTANT.f90 \
       Specific_heat.f90 \
       Omega_integrals.f90 \
       Bracket_integrals.f90 \
       invers.f90 \
       co2-o2_transp.f90

# using intel MKL -qmkl
#$FC *.o -pg -o co2_harm_6T.x

# using OpenBlas -lopenblas
$FC *.o -pg -o co2_harm_6T.x

# Cray
#$FC *.o -o co2_harm_6T.x

time ./co2_harm_6T.x
