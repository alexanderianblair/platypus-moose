# MFEMVectorFEInnerProductIntegralPostprocessor

!if! function=hasCapability('mfem')

## Overview

Postprocessor for calculating the integral of the dot product between two vector variables over a
user-specified subdomain of the mesh:

!equation
(k \vec u, \vec v)_\Omega

where $\vec u$ is the primal variable and $\vec v$ the dual variable, $k$ is an optional scalar
coefficient, and $\Omega$ is the user-specified mesh subdomain. The primal variable may be in
$H(\mathrm{curl})$ or $H(\mathrm{div})$, or in a vector space built from scalar elements, such as
$[L^2]^d$ or $[H^1]^d$; the dual variable may be any vector variable.

## Example Input File Syntax

!listing mfem/submeshes/av_magnetostatic.i block=Postprocessors/CoilPower

!syntax parameters /Postprocessors/MFEMVectorFEInnerProductIntegralPostprocessor

!syntax inputs /Postprocessors/MFEMVectorFEInnerProductIntegralPostprocessor

!syntax children /Postprocessors/MFEMVectorFEInnerProductIntegralPostprocessor

!if-end!

!else
!include mfem/mfem_warning.md
