# MFEMNetCurrentIntegralConstraint

!if! function=hasCapability('mfem')

## Summary

Weakly constrains the net current passing through an interior cut surface of a topologically
closed conductor, solving for the loop voltage that drives it.

## Overview

This is the current-driven counterpart of the loop voltage constraint described in
[ClosedCoilMagnetostatic.md]. There, the loop voltage $\mathcal V$ is prescribed strongly and the
resulting current falls out of the solution; here the net current $I$ is prescribed and the loop
voltage is solved for, as the multiplier of an [MFEMIntegralConstraint.md].

The conductor $\Omega_c$ is made simply connected by a cut surface $\Gamma_c$, and a one element
wide transition subdomain $\Omega_t$ on one side of it is built with an
[MFEMCutTransitionSubMesh.md]. The externally applied electric field is represented in $\Omega_t$
by

\begin{equation}
\vec E_{ext} = -\lambda \vec \nabla w
\end{equation}

where $w$ is the unit cut function given by `cut_function`, taking the value one on $\Gamma_c$ and
zero on the remainder of $\partial \Omega_t$, and $\lambda$ is the loop voltage held by
`scalar_variable`. With $\phi$ the induced scalar potential given by `variable` and $\sigma$ the
conductivity given by `coefficient`, the total current density in $\Omega_t$ is
$\vec J = -\sigma \vec \nabla (\phi + \lambda w)$.

Taking $K$ to be the stiffness matrix $(\sigma \vec \nabla \cdot, \vec \nabla \cdot)_{\Omega_t}$
assembled over `block`, the constraint contributes the coupling vector $c = K w$ and the diagonal
$d = w^T K w$. Its row is therefore

\begin{equation}
c^T \phi + d \lambda
= (\sigma \vec \nabla (\phi + \lambda w), \vec \nabla w)_{\Omega_t}
= -(\vec J, \vec \nabla w)_{\Omega_t}
= - \int_{\Gamma_c} \vec J \cdot \hat n \, dS
= -I
\end{equation}

where the last two steps integrate by parts over $\Omega_t$ and use $\vec \nabla \cdot \vec J = 0$
together with $\vec J \cdot \hat n = 0$ on the insulated conductor surface, leaving only the cut
surface, on which $w = 1$, contributing to the boundary term.

Restricting $K$ to the transition subdomain with `block` is essential: it is what makes
$\vec \nabla w$ extended by zero a field with non-zero circulation around the conductor, rather
than the gradient of a single-valued potential, which would have no circulation at all.

The prescribed `current` is counted positive in the direction of the normal $\hat n$ of the cut
surface pointing out of the transition subdomain. Which of the two possible directions this is
depends on the side of the cut the transition subdomain was built on.

!alert note
The coupling vector is simultaneously the source term
$\lambda (\sigma \vec \nabla w, \vec \nabla \phi')_{\Omega_t}$ that the loop voltage contributes to
the equation for $\phi$. A problem using this constraint must therefore *not* also supply that
source as an `MFEMMixedGradGradKernel`, as the voltage-driven input in
[ClosedCoilMagnetostatic.md] does.

## Example Input File Syntax

!listing test/tests/mfem/constraints/closed_coil_net_current.i block=SubMeshes Variables AuxVariables ICs Kernels Constraints

!syntax parameters /Constraints/MFEMNetCurrentIntegralConstraint

!syntax inputs /Constraints/MFEMNetCurrentIntegralConstraint

!syntax children /Constraints/MFEMNetCurrentIntegralConstraint

!if-end!

!else
!include mfem/mfem_warning.md
