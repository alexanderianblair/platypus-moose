# MFEMNLDiffusionKernel

!if! function=hasCapability('mfem')

## Overview

Adds the domain integrator for integrating the non-linear form

!equation
(k(u)\vec\nabla u, \vec\nabla v)_\Omega \,\,\, \forall v \in V

where $u, v \in H^1$ and $k(u)$ is a scalar non-linear diffusivity coefficient that may depend on
the trial variable $u$.

The above terms arises from the weak form of the non-linear operator

!equation
- \vec\nabla \cdot \left( k(u) \vec\nabla u \right)

## Coupling to other variables

The diffusivity may also depend on other variables solved for alongside $u$. Naming those
variables in `coupled_variables`, along with one derivative coefficient each in
`dk_dcoupled_coefficients`, adds the off-diagonal Jacobian contribution

!equation
\left(\frac{\partial k}{\partial c}\phi_c \vec\nabla u, \vec\nabla v\right)_\Omega

for each coupled variable $c$, where $\phi_c$ are the basis functions of the space $c$ is
discretised on. Only variables that are solved for need to be named: a dependence of $k$ on any
other variable is already carried by the coefficient itself when the residual is evaluated.

## Example Input File Syntax

!listing mfem/nonlinear/nldiffusion_common.i block=/Kernels

!syntax parameters /Kernels/MFEMNLDiffusionKernel

!syntax inputs /Kernels/MFEMNLDiffusionKernel

!syntax children /Kernels/MFEMNLDiffusionKernel

!if-end!

!else
!include mfem/mfem_warning.md
