# MFEMScalarVariable

!if! function=hasCapability('mfem')

## Overview

`MFEMScalarVariable` defines a global scalar unknown of an MFEM problem: a variable holding a
single degree of freedom for the whole problem, rather than a field discretized on a finite
element space. It is the MFEM analogue of a MOOSE scalar variable.

MFEM provides no finite element space with a single global degree of freedom, so unlike
[MFEMVariable.md] this variable is not backed by an `mfem::ParGridFunction`, and it takes no
`fespace` parameter. Its degree of freedom is owned by rank 0 and occupies a single-entry block
of the block system appended after the blocks of the field variables, while its value is stored
redundantly on every rank so that it may be read anywhere.

Scalar variables are coupled to field variables by [MFEMIntegralConstraint.md], which uses one to
hold the multiplier of a weakly enforced integral constraint. Because no MFEM kernel currently
contributes to a scalar variable's equation, a scalar variable is only useful in a problem that
also has an integral constraint naming it. Its solved value can be reported with
[MFEMScalarVariableValue.md].

## Example Input File Syntax

!listing test/tests/mfem/constraints/closed_coil_net_current.i block=Variables Constraints

!syntax parameters /Variables/MFEMScalarVariable

!syntax inputs /Variables/MFEMScalarVariable

!syntax children /Variables/MFEMScalarVariable

!if-end!

!else
!include mfem/mfem_warning.md
