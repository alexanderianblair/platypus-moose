# MFEMNLCurlCurlKernel

!if! function=hasCapability('mfem')

## Overview

Adds the domain integrator for integrating the non-linear form

!equation
(k(|\vec\nabla \times \vec u|)\vec\nabla \times \vec u, \vec\nabla \times \vec v)_\Omega \,\,\,
\forall \vec v \in V

where $\vec u, \vec v \in H(\mathrm{curl})$ and $k$ is a scalar coefficient that may depend on the
magnitude of the curl of the trial variable $\vec u$.

This term arises from the weak form of the curl curl operator

!equation
\vec\nabla \times \left(k(|\vec\nabla \times \vec u|) \vec\nabla \times \vec u\right)

The coefficient $k$ may also depend on scalar variables $c_i$, such as temperature, listed in
[!param](/Kernels/MFEMNLCurlCurlKernel/coupled_variables). If a $c_i$ is solved for in the same
equation system as $\vec u$, the Jacobian includes the off-diagonal contribution

!equation
\left(\frac{\partial k}{\partial c_i} \delta c_i \vec\nabla \times \vec u,
\vec\nabla \times \vec v\right)_\Omega

formed from [!param](/Kernels/MFEMNLCurlCurlKernel/dk_dcoupled_coefficients).

## Example Input File Syntax

!listing mfem/submeshes/nl_hphi_magnetodynamic.i block=/Kernels

The magnetic field and temperature of a conductor with temperature-dependent resistivity can be
solved for in one equation system:

!listing mfem/nonlinear/nl_electrothermal.i block=/Kernels

!syntax parameters /Kernels/MFEMNLCurlCurlKernel

!syntax inputs /Kernels/MFEMNLCurlCurlKernel

!syntax children /Kernels/MFEMNLCurlCurlKernel

!if-end!

!else
!include mfem/mfem_warning.md
