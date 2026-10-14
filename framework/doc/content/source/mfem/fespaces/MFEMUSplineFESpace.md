# MFEMUSplineFESpace

!if! function=hasCapability('mfem')

## Overview

Builds a finite element space of the rational U-spline functions of an [MFEMUSplineMesh.md], for
isogeometric analysis (IGA). The degrees of freedom of the space are the spline functions
themselves, so the space has as many degrees of freedom per component as the U-spline has control
points, and its basis functions have the smoothness of the U-spline across element boundaries. The
order of the space is the degree of the spline. Vector-valued spaces, for example for the
displacement in linear elasticity, are built by setting `vdim`.

On each Bezier element, the space evaluates the spline functions supported on the element through
their extraction operator. Since the same functions represent the geometry, the space contains all
linear functions, and the solution of a patch test is exact.

Spline functions do not interpolate the solution at nodes, so [MFEMScalarDirichletBC.md] and
[MFEMVectorDirichletBC.md] set the coefficients of the spline functions with nonzero traces on
the boundary from local $L^2$ projections of the boundary values onto those traces, averaged over
the boundary elements each function's trace is supported on. Boundary values that lie in the trace
space of the U-spline are reproduced exactly.

Every Bezier element must support exactly $(p+1)^d$ spline functions, where $p$ is the degree and
$d$ the dimension, as MFEM projects coefficients, such as initial conditions, onto NURBS spaces
through a square local system on each element. The space cannot be built on a submesh.

## Example Input File Syntax

!listing test/tests/mfem/uspline/uspline_diffusion.i block=Mesh FESpaces Variables

!syntax parameters /FESpaces/MFEMUSplineFESpace

!syntax inputs /FESpaces/MFEMUSplineFESpace

!syntax children /FESpaces/MFEMUSplineFESpace

!if-end!

!else
!include mfem/mfem_warning.md
