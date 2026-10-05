# MFEMMixedNormalCrossGradKernel

!if! function=hasCapability('mfem')

## Overview

Adds the domain integrator for integrating the mixed bilinear form

!equation
(k \hat n \times \vec \nabla_\Gamma u, \vec v)_\Gamma \,\,\, \forall \vec v \in V

on a two-dimensional surface $\Gamma$ embedded in three dimensions, where $u \in H^1(\Gamma)$,
$\vec v \in H(\mathrm{curl})$ or $H(\mathrm{div})$, $\hat n$ is the unit normal to each element of
$\Gamma$, and $k$ is a scalar coefficient. The normal is oriented by the vertex ordering of the
element, which for a surface extracted from a sideset points out of the elements the sideset was
defined on.

The surface $\Gamma$ is the domain of the variable defined on it, which is typically defined on
an [MFEMBoundarySubMesh.md], while $\vec v$ may be defined on the parent mesh. In that case, the
form is integrated over the elements of the submesh.

This term arises in the T-A formulation of thin superconducting tapes, where the sheet current
density in a tape of thickness $d$ is $d \vec\nabla_\Gamma T \times \hat n$ for a scalar current
vector potential $T$ defined on the tape surface. Setting `transpose = true` gives the term
$(k \vec u, \hat n \times \vec\nabla_\Gamma v)_\Gamma$, which equals
$-(k \hat n \cdot \vec\nabla \times \vec u, v)_\Gamma$ for constant $k$ when $v$ vanishes on the
boundary of $\Gamma$.

## Example Input File Syntax

!listing test/tests/mfem/submeshes/ta_tape.i block=Kernels

!syntax parameters /Kernels/MFEMMixedNormalCrossGradKernel

!syntax inputs /Kernels/MFEMMixedNormalCrossGradKernel

!syntax children /Kernels/MFEMMixedNormalCrossGradKernel

!if-end!

!else
!include mfem/mfem_warning.md
