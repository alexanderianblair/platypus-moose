# MFEMOrthoSolver

!if! function=hasCapability('mfem')

## Overview

Defines and builds an `mfem::OrthoSolver`, which applies the solver named in
[!param](/Solvers/MFEMOrthoSolver/solver) on the subspace orthogonal to the constant vector. The
mean value is subtracted from both the right hand side passed to the wrapped solver and from the
solution it returns.

This makes singular systems whose nullspace is spanned by the constant vector solvable, as arise
for example from problems with pure Neumann boundary conditions. The zero-mean solution is selected
from the family of solutions differing by a constant.

Note that `mfem::OrthoSolver` overwrites the `iterative_mode` of the solver it wraps, so
[!param](/Solvers/MFEMOrthoSolver/use_initial_guess) must be set on this object rather than on the
wrapped solver. Note also that the orthogonalisation is performed on host, so this solver is not
supported on non-CPU devices.

## Example Input File Syntax

!listing test/tests/mfem/solvers/ortho_neumann_diffusion.i block=Solvers

!syntax parameters /Solvers/MFEMOrthoSolver

!syntax inputs /Solvers/MFEMOrthoSolver

!syntax children /Solvers/MFEMOrthoSolver

!if-end!

!else
!include mfem/mfem_warning.md
