# MFEMNLJouleHeatingKernel

!if! function=hasCapability('mfem')

## Overview

Adds the domain integrator for the nonlinear Joule heating source

!equation
(k |\vec\nabla \times \vec u|^2, v)_\Omega \,\,\, \forall v \in V

where $v$ is a scalar test function, $\vec u \in H(\mathrm{curl})$ is the coupled variable set by
[!param](/Kernels/MFEMNLJouleHeatingKernel/curl_variable) whose curl is the current density, and
$k$ is the resistivity. In the H formulation, $\vec u$ is the magnetic field and the source is the
heating rate $\vec E \cdot \vec J = k |\vec J|^2$ of a conductor obeying
$\vec E = k \vec J$.

The resistivity may depend on $|\vec\nabla \times \vec u|$ and on scalar variables $c_i$ listed in
[!param](/Kernels/MFEMNLJouleHeatingKernel/coupled_variables), which may include the kernel
variable itself. The source is evaluated in the residual at each nonlinear iteration. If
$\vec u$ or any $c_i$ is solved for in the same equation system as the kernel variable, the
Jacobian includes the derivatives

!equation
\frac{\partial (k |\vec J|^2)}{\partial \vec J} = \left(2 k + |\vec J|
\frac{\partial k}{\partial |\vec J|}\right) \vec J, \quad
\frac{\partial (k |\vec J|^2)}{\partial c_i} = \frac{\partial k}{\partial c_i} |\vec J|^2

formed from [!param](/Kernels/MFEMNLJouleHeatingKernel/curlu_dk_dcurlu_coefficient) and
[!param](/Kernels/MFEMNLJouleHeatingKernel/dk_dcoupled_coefficients).

## Example Input File Syntax

!listing mfem/nonlinear/nl_electrothermal.i block=/Kernels

!syntax parameters /Kernels/MFEMNLJouleHeatingKernel

!syntax inputs /Kernels/MFEMNLJouleHeatingKernel

!syntax children /Kernels/MFEMNLJouleHeatingKernel

!if-end!

!else
!include mfem/mfem_warning.md
