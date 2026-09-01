# MFEMCohomologyCutIC

!if! function=hasCapability('mfem')

## Overview

An `MFEMCohomologyCutIC` sets a first order $H(\mathrm{curl})$ conforming (Nedelec) variable to a
cohomology basis cochain computed by [Gmsh](https://gmsh.info)'s `Cohomology` command, scaled by a
user-specified amplitude.

On a domain that is not simply connected, such as the air region surrounding a toroidal conductor,
a curl-free field can have a non-zero circulation about a loop linking the hole. Such a field cannot
be written as the gradient of a single-valued scalar potential, and is instead represented by a
generator of the first cohomology space $H^1$ of the domain. The resulting field is curl-free
throughout the domain, and its circulation about any loop linking the hole is the amplitude given,
so it can be used directly as the source field of a global current constraint in an $H$-$\phi$
formulation, or of a global loop voltage constraint in an $A$-$V$ formulation.

This provides an alternative to building a cut surface in the geometry and using an
[MFEMCutTransitionSubMesh.md] to define a transition region, as described in
[ClosedCoilMagnetostatic.md]. The cut is instead computed by Gmsh, using the coreduction algorithm
of [M. Pellikka, S. Suuriniemi, L. Kettunen and C. Geuzaine, SIAM Journal on Scientific Computing,
35, no. 5 (2013): B1195-B1214.](https://doi.org/10.1137/130906556), and no cut surface need be
present in the geometry.

## Producing the cochain in Gmsh

Gmsh stores the computed basis cochains as physical groups of line elements in the mesh file. To
request the generators of $H^1$ of the region with physical tag 1, add to the `.geo` file

```
Cohomology(1) {{1}, {}};
```

or, using the Gmsh Python API,

```python
gmsh.model.mesh.addHomologyRequest("Cohomology", domainTags=[1], subdomainTags=[], dims=[1])
gmsh.model.mesh.computeHomology()
```

Gmsh names each resulting group after the space and the domain it was computed on, so the first
basis cochain of $H^1$ of the domain of physical group 1 is named `H^1{1}1`. That name is what
should be given to the `cut_name` parameter.

A group of line elements lies two dimensions below a three dimensional mesh, so it takes no part in
the mesh topology and has no place in the element or boundary arrays. MFEM's Gmsh reader keeps such
groups in `mfem::Mesh::subdim_entities`, and [MFEMMesh.md] takes them from the serial mesh as it is
read, before it is partitioned. The cochain therefore always comes from the file the mesh itself was
built from.

## Limitations

- The mesh format must be one whose reader keeps groups lying below the mesh boundary, which at
  present means a Gmsh mesh. Both the 2.2 and 4.1 formats work, in ASCII or binary.
- The variable must be defined on a `FIRST` order `ND` finite element space, since the cochain
  supplies a single value per mesh edge.
- The mesh must not be refined after being read, as refinement replaces the edges the cochain is
  defined on. Refine the geometry in Gmsh and recompute the cohomology there instead of setting
  `serial_refine`, `uniform_refine` or `parallel_refine` on the `[Mesh]` block.
- The field is only curl-free on the domain the cochain was computed for, so kernels using it as a
  source should be restricted to that domain.
- The variable must be defined on the parent mesh rather than on a submesh. The sign of an edge's
  degree of freedom follows the order of that edge's vertices in the local numbering of the mesh;
  an `mfem::ParMesh` takes that numbering from the serial mesh, so ranks sharing an edge agree on
  it, but an `mfem::ParSubMesh` builds its own and they need not. To use the cochain on a submesh,
  set it on a variable of the parent mesh and move it across with an [MFEMSubMeshTransfer.md], as
  in the magnetostatic example below.

## Example Input File Syntax

!listing test/tests/mfem/cohomology/cohomology_cut_source.i block=ICs

The example below solves for the magnetic field about a toroidal conductor carrying a prescribed
current, with the global current constraint carried by the imported cochain rather than by a cut
surface built in the geometry:

!listing test/tests/mfem/cohomology/cohomology_cut_magnetostatic.i

!syntax parameters /ICs/MFEMCohomologyCutIC

!syntax inputs /ICs/MFEMCohomologyCutIC

!syntax children /ICs/MFEMCohomologyCutIC

!if-end!

!else
!include mfem/mfem_warning.md
