# MFEMMesh

!if! function=hasCapability('mfem')

## Overview

`MFEMMesh` is responsible for building an `mfem::ParMesh` object from the provided mesh input file
for use in an `MFEMProblem`. Exodus files are supported, along with other mesh formats listed
 [here](https://mfem.org/mesh-formats/).

As MOOSE checks for the existence of a `libMesh` MOOSE mesh at various points during setup,
`MFEMMesh` currently builds a dummy MOOSE mesh of a single point alongside the MFEM mesh. This dummy
mesh should not be used in an `MFEMProblem`; all MFEM objects should access the `mfem::ParMesh` via
the `getMFEMParMesh()` accessor as needed.

## Mesh Displacement

Setting `displacement` to the name of a vector MFEM variable makes the mesh move by that field
after every solve, so that subsequent solves are assembled on the deformed geometry.

By default the displacement is treated as an *increment*: it is added to the current node
positions each time the mesh is displaced, so a mesh displaced on every step of a transient
accumulates successive displacements. Set `displacement_is_total = true` when the variable holds
the *total* displacement measured from the undeformed mesh instead. The undeformed node positions
are then stored the first time the mesh is displaced and the mesh is repositioned relative to
them, which makes repeated displacement idempotent. That is what a coupled solve driven by a
total displacement field needs, since such a mesh is displaced once per solve and once per fixed
point iteration within a solve.

## Example Input File Syntax

!listing test/tests/mfem/kernels/diffusion.i block=Problem Mesh

And displacing the mesh by a total displacement received from a coupled application:

!listing modules/solid_mechanics/test/tests/mfem_lorentz_force/lorentz_force_coil_displaced.i block=Mesh

!syntax parameters /Mesh/MFEMMesh

!syntax inputs /Mesh/MFEMMesh

!syntax children /Mesh/MFEMMesh

!if-end!

!else
!include mfem/mfem_warning.md
