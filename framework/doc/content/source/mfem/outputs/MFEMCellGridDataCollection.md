# MFEMCellGridDataCollection

!if! function=hasCapability('mfem')

## Overview

`MFEMDataCollection` controlling output of data as VTK cell grids (`.dg` files) for visualisation
in ParaView, which opens them with the "Composite CellGrid Reader". A cell grid stores each field
as the coefficients of a finite element basis on every cell, rather than as values sampled at
points. This means that:

- discontinuous (L2) variables are represented without duplicating points;
- curved (second-order) meshes and variables of order up to two are represented exactly;
- lowest-order Nédélec (ND) and Raviart-Thomas (RT) variables on tetrahedra, hexahedra and wedges
  are stored in the native VTK `HCURL` and `HDIV` function spaces, preserving their tangential or
  normal continuity, when [!param](/Outputs/MFEMCellGridDataCollection/native_vector_spaces) is
  true. Other ND and RT variables are interpolated to discontinuous nodal vector fields.

The VTK cell-grid bases have a polynomial order of at most two (one for edges and pyramids).
Meshes or variables of higher order are written on uniformly subdivided elements, which
approximates them with an error that decreases with
[!param](/Outputs/MFEMCellGridDataCollection/refinements).

The file opened in ParaView is `<file_base>/Run0/Run0_<cycle>.dg`. It is a small JSON index listing
one data file per process, `<file_base>/Run0/Cycle<cycle>/proc<rank>.dg`, which holds the data in
the chosen [!param](/Outputs/MFEMCellGridDataCollection/encoding). The compact binary `MESSAGEPACK`
encoding is recommended; `JSON` is human-readable. Each output step is a separate file, since
ParaView's cell-grid reader does not yet support time series; the simulation time is stored in
each file as metadata.

Cell grids require VTK 9.3 or later (ParaView 5.12 or later); VTK 9.7 or later is recommended.
Writing them does not require VTK.

## Example Input File Syntax

!listing test/tests/mfem/kernels/diffusion.i block=Outputs/CellGridDataCollection

!syntax parameters /Outputs/MFEMCellGridDataCollection

!syntax inputs /Outputs/MFEMCellGridDataCollection

!syntax children /Outputs/MFEMCellGridDataCollection

!if-end!

!else
!include mfem/mfem_warning.md
