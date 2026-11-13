# MFEMJumpMassInterfaceKernel

!if! function=hasCapability('mfem')

## Overview

Adds the bilinear form

!equation
(k [u], [v])_\Gamma \,\,\, \forall v \in V

over an interior boundary $\Gamma$, where $[u] = u - u_n$ is the jump from the primary
variable $u$ to the neighbor variable $u_n$ across $\Gamma$, $[v]$ the jump in the test
function, and $k$ a scalar coefficient. Added to a diffusion problem, it imposes a flux
$k [u]$ from the primary side of $\Gamma$ to the neighbor side, as for an electrical or
thermal contact conductance $k$.

An $H^1$ variable is continuous across every face of the mesh, so a field that jumps
across $\Gamma$ is represented by two variables, one for each side. Each is defined over
the whole mesh, and

- governed by kernels restricted to the blocks on its own side, and
- fixed elsewhere, other than on its interfaces with its own side, by an
  [MFEMScalarEssentialConstraint.md] with `constrain_subdomain_boundary = false` and
  `constrain_exterior_boundary = true`.

On $\Gamma$ the two variables have degrees of freedom of their own, which this kernel
couples. Where the blocks meet outside $\Gamma$ nothing couples them, so no flux crosses.

The boundary given must contain each face of $\Gamma$ once. The kernel is not yet
supported in time-dependent or complex problems.

## Example Input File Syntax

!listing mfem/interfacekernels/jump_mass.i block=InterfaceKernels Constraints

!syntax parameters /InterfaceKernels/MFEMJumpMassInterfaceKernel

!syntax inputs /InterfaceKernels/MFEMJumpMassInterfaceKernel

!syntax children /InterfaceKernels/MFEMJumpMassInterfaceKernel

!if-end!

!else
!include mfem/mfem_warning.md
