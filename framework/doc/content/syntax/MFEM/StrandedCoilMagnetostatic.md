# Magnetostatic Problem on a Homogenised Stranded Conductor

## Summary

Solves for the magnetic field around a topologically closed toroidal coil carrying a specified
 total current, where the coil is treated as a homogenised stranded conductor rather than a solid
 one. The direction the current follows around the coil is obtained from the coil geometry alone,
 with no vector function for the current density given in the input.

## Description

[ClosedCoilMagnetostatic.md] solves for the field around a *solid* toroidal conductor driven by a
 global loop voltage $\mathcal V$. There the current density follows Ohm's law, $\vec J = \sigma
 \vec E$, and for a torus of major radius $R$ this gives

\begin{equation}
    \vec J = \frac{\sigma \mathcal V}{2 \pi r} \hat \phi,
\end{equation}

so that $|\vec J|$ varies as $1/r$ over the cross-section, being largest on the inboard side.

A coil wound from many fine strands behaves differently. The strands are individually insulated, so
 no current redistribution across the cross-section is possible, and eddy currents within each
 strand are negligible. Homogenising over the turns, the current density is instead taken to have
 constant magnitude across the cross-section,

\begin{equation}
    \vec J = \frac{I}{S} \hat t,
\end{equation}

where $I$ is the total current through the coil, $S$ its cross-sectional area, and $\hat t$ the unit
 vector tangent to the turns. For a coil of uniform cross-section this remains a valid source
 current; in the toroidal case $\hat t = \hat \phi$ and

\begin{equation}
    \vec \nabla \cdot \left( f(r, z) \hat \phi \right)
    = \frac{1}{r} \partial_\phi f
    = 0,
\end{equation}

so $\vec \nabla \cdot \vec J = 0$ is satisfied as required for a magnetostatic source.

## Obtaining the Current Direction from the Geometry

Rather than prescribe $\hat t$ as a function of position, which would require knowing the path of
 the coil analytically, it is recovered from the coil geometry. Solving the electrokinetic problem
 of [ClosedCoilMagnetostatic.md] on the closed conductor produces an electric field $\vec E$ that is
 everywhere tangential to the coil and satisfies $\vec E \cdot \hat n = 0$ on $\partial \Omega_c$,
 so

\begin{equation}
    \hat t = \frac{\vec E}{|\vec E|}.
\end{equation}

Since only the direction of $\vec E$ is used, the conductivity and loop voltage assigned in that
 problem are immaterial; the magnitude of the loop voltage divides out, and reversing its sign
 reverses $\vec E$ and $S$ below together, leaving $\vec J$ unchanged. The sense in which the
 current circulates is instead fixed by the orientation chosen for the cut surface $\Gamma_c$, which
 is a property of the coil geometry.

The cross-sectional area follows from the same field, as the flux of the unit tangent through the
 cut surface,

\begin{equation}
    S = \int_{\Gamma_c} \hat t \cdot d \vec S,
\end{equation}

which is evaluated with an [MFEMVectorBoundaryFluxIntegralPostprocessor.md] scaled by the
 coefficient $1/|\vec E|$. Measuring $S$ this way returns the area of the discretised cut surface
 rather than that of the underlying geometry, which is what makes the total current through the
 discretised coil come out equal to $I$.

This is shown explicitly in the example file:

!listing test/tests/mfem/submeshes/stranded_coil_source.i

## Assembling the Current Density

The tangential field and the measured area are passed to the magnetostatic problem, where the
 source current density is assembled as the scaled field

\begin{equation}
    \vec J = \frac{I}{S |\vec E|} \vec E,
\end{equation}

using an `MFEMMixedVectorMassKernel` with $\vec E$ as its trial variable and the scalar coefficient
 $I / (S |\vec E|)$, built with an [MFEMParsedFunction.md]. The area enters that expression as a
 postprocessor value; postprocessor values are available as spatially uniform scalar coefficients,
 including values transferred in from a subapp.

The magnetostatic problem itself is then as in [ClosedCoilMagnetostatic.md],

\begin{equation}
(\mu^{-1} \vec \nabla \times \vec A, \vec \nabla \times \vec A')_\Omega
- (\vec H \times \vec n, \vec A')_{\partial \Omega}
+ (\vec J, \vec A')_{\Omega_c} = 0 \,\,\, \forall \vec A' \in V_A
\end{equation}

with the same caveat that the problem is singular in non-conductive regions, so that either
`singular = true` should be passed to the `HypreAMS` preconditioner or a small mass term added for
stability.

The total current is recovered by integrating $\vec J$ over a measurement plane cutting the coil
 elsewhere than $\Gamma_c$, which checks that the homogenised current density is normalised
 correctly. The sign of that integral depends on the orientation of the measurement plane relative
 to that of the cut surface.

## Example File

The full homogenised stranded coil magnetostatic example detailed above can be found below:

!listing test/tests/mfem/submeshes/stranded_coil_magnetostatic.i
