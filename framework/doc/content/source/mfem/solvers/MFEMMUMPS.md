# MFEMMUMPS

!if! function=hasCapability('mfem')

## Overview

Defines and builds an `mfem::MUMPSSolver` to use as a direct solver to solve the MFEM equation system.

By default, a complex-valued system is factorized as its equivalent real-valued system, which has
twice as many unknowns. Setting [!param](/Solvers/MFEMMUMPS/complex_arithmetic) to `true` instead
factorizes the complex-valued system directly using `mfem::ComplexMUMPSSolver`. This requires a
complex weak form assembled with `LEGACY` assembly.

!syntax parameters /Solvers/MFEMMUMPS

!syntax inputs /Solvers/MFEMMUMPS

!syntax children /Solvers/MFEMMUMPS

!if-end!

!else
!include mfem/mfem_warning.md
