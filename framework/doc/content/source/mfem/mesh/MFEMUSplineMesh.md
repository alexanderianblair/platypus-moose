# MFEMUSplineMesh

!if! function=hasCapability('mfem')

## Overview

`MFEMUSplineMesh` reads a U-spline in Bezier extraction form from a BEXT 2.0 file, as written by
Coreform Cubit with `export uspline ... bext`, for isogeometric analysis (IGA) with
[MFEMUSplineFESpace.md].

A BEXT file holds the spline control points and their weights, and for each Bezier element the
spline functions supported on it together with their extraction operator, which expresses each of
those functions in terms of the Bernstein polynomials of the element. The rational spline
functions are

!equation
R_A = \frac{w_A N_A}{\sum_B w_B N_B}, \qquad N_A |_e = \sum_j C^e_{Aj} B_j,

where $w_A$ are the spline weights, $C^e$ is the extraction operator of element $e$, and $B_j$ are
its tensor product Bernstein polynomials.

The MFEM mesh holds the Bezier elements as a NURBS mesh with one single-element patch per element.
Each patch is the rational Bezier representation of the spline geometry on its element, so the
geometry, including curved boundaries, is exact. When the mesh is read, adjacent elements are
checked to agree on the spline combinations forming their shared Bezier control points.

The following restrictions apply:

- Elements must be quadrilaterals or hexahedra of a single polynomial degree, the same in every
  parametric direction. Two-dimensional meshes must lie in the $z=0$ plane.
- Refinement and reordering of the mesh are not available, as they would invalidate the
  extraction operator.
- BEXT files hold no sidesets or subdomains. All elements have attribute 1, and all boundary
  elements have boundary attribute 1.
- Recovery is not supported.

## Example Input File Syntax

!listing test/tests/mfem/uspline/uspline_diffusion.i block=Problem Mesh

!syntax parameters /Mesh/MFEMUSplineMesh

!syntax inputs /Mesh/MFEMUSplineMesh

!syntax children /Mesh/MFEMUSplineMesh

!if-end!

!else
!include mfem/mfem_warning.md
