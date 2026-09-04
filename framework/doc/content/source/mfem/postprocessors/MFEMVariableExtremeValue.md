# MFEMVariableExtremeValue

!if! function=hasCapability('mfem')

## Overview

Postprocessor for finding the maximum or minimum value taken by a scalar MFEM
variable, optionally restricted to a subset of the subdomains of the mesh.

!equation
\max_{i} u_i \qquad \mathrm{or} \qquad \min_{i} u_i

where $u_i$ are the degrees of freedom of $u \in H^1 \lor L^2$ that belong to
the selected subdomains. For the nodal bases used by $H^1$ and $L^2$ spaces
these are the nodal values of the variable, so the result is the extreme nodal
value rather than the extremum of the interpolant.

Which extremum is returned is set by [!param](/Postprocessors/MFEMVariableExtremeValue/value_type).

Restricting to subdomains makes this usable for the flux extrema that define a
plasma equilibrium: the flux on the magnetic axis is the maximum over the
vessel, and for a limited plasma the flux on the plasma boundary is the maximum
over the limiter. Both are used in [this way](syntax/MFEM/FreeBoundaryEquilibrium.md).

## Example Input File Syntax

!listing mfem/gradshafranov/freeboundary.i block=Postprocessors

!syntax parameters /Postprocessors/MFEMVariableExtremeValue

!syntax inputs /Postprocessors/MFEMVariableExtremeValue

!syntax children /Postprocessors/MFEMVariableExtremeValue

!if-end!

!else
!include mfem/mfem_warning.md
