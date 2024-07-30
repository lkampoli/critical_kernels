I have identified some critical kernels in both kinetic and transport modules.

They compile and run with intel, gnu and cray compilers, but ...

There are some issues in the kinetic part, specifically in brent module, which
produce an error (non initialized field) in some cases. As workaround, it was
commented in kvt_fho.f90 and this makes cray to complain about it. To be fixed,
but not relevant for the sake of optimization.

Morever, to compile with ifort the Makefile is preferable, for the kinetic
module, otherwise the compile.sh is ok.
