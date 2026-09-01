# MFEMVectorFromScalarsAux

!if! function=hasCapability('mfem')

## Overview

AuxKernel that assembles a vector MFEM auxiliary variable from one scalar
coefficient per component,

!equation
\vec u = \sum_i u_i \hat e_i

where each $u_i$ is a scalar coefficient, which may be a scalar variable, a
function, a material property or a literal value. This is the inverse of
[MFEMInnerProductAux.md], which extracts a single component of a vector.

It is needed wherever a vector field has to be rebuilt from components that were
produced separately. Transfers to and from MFEM applications only support scalar
variables, so a vector field arriving from another application arrives one
component at a time and must be reassembled before it can be used as a vector,
for instance as the displacement consumed by [MFEMMesh.md].

The number of coefficients supplied must match the vector dimension of the
target auxiliary variable.

## Example Input File Syntax

!listing modules/solid_mechanics/test/tests/mfem_lorentz_force/lorentz_force_coil_displaced.i block=/AuxKernels/assemble_displacement

!syntax parameters /AuxKernels/MFEMVectorFromScalarsAux

!syntax inputs /AuxKernels/MFEMVectorFromScalarsAux

!syntax children /AuxKernels/MFEMVectorFromScalarsAux

!if-end!

!else
!include mfem/mfem_warning.md
