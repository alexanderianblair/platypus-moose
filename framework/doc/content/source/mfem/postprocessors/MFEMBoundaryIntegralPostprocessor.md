# MFEMBoundaryIntegralPostprocessor

!if! function=hasCapability('mfem')

## Overview

Postprocessor for calculating the integral of a scalar coefficient over a set of
boundaries,

!equation
(k, 1)_{\partial\Omega}

where $k$ is a scalar coefficient and $\partial\Omega$ the user-specified boundaries,
which may lie in the interior of the mesh. The coefficient may depend on variables, as
through an [MFEMParsedFunction.md]. On an interior boundary an $H^1$ variable takes its
trace on the boundary, whichever element it is evaluated from. The integral is computed
with a quadrature rule of order `quadrature_order` on each boundary element.

## Example Input File Syntax

!listing mfem/interfacekernels/jump_mass.i block=Postprocessors

!syntax parameters /Postprocessors/MFEMBoundaryIntegralPostprocessor

!syntax inputs /Postprocessors/MFEMBoundaryIntegralPostprocessor

!syntax children /Postprocessors/MFEMBoundaryIntegralPostprocessor

!if-end!

!else
!include mfem/mfem_warning.md
