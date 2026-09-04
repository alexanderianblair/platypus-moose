# Grad-Shafranov Problem

## Summary

Solves the Grad-Shafranov equation for the poloidal flux of an axisymmetric
plasma equilibrium on the $(R, z)$ half-plane. The particular case solved is a
Solov'ev equilibrium, for which the source terms are constant and a closed-form
solution is available, so the discrete solution can be verified directly.

## Description

The magnetic field of an axisymmetric plasma may be written in terms of the
poloidal flux function $\psi(R, z)$ and the poloidal current function
$F(\psi) = R B_\phi$ as

\begin{equation}
\vec B = \frac{1}{R} \vec \nabla \psi \times \hat e_\phi + \frac{F(\psi)}{R} \hat e_\phi .
\end{equation}

Requiring force balance $\vec J \times \vec B = \vec \nabla p$ gives the
Grad-Shafranov equation,

\begin{equation}
\Delta^* \psi = -\mu_0 R^2 p'(\psi) - F(\psi) F'(\psi),
\qquad
\Delta^* \psi \equiv R \frac{\partial}{\partial R}
  \left( \frac{1}{R} \frac{\partial \psi}{\partial R} \right)
  + \frac{\partial^2 \psi}{\partial z^2} ,
\end{equation}

where $p$ is the plasma pressure and a prime denotes differentiation with
respect to $\psi$.

The elliptic operator can be written in divergence form as
$\Delta^* \psi = R \vec \nabla \cdot \left( \frac{1}{R} \vec \nabla \psi \right)$,
so dividing through by $R$ and testing against $v$ gives the weak form solved
here,

\begin{equation}
\left( \frac{1}{R} \vec \nabla \psi, \vec \nabla v \right)_\Omega
  = \left( \mu_0 R p'(\psi) + \frac{F(\psi) F'(\psi)}{R}, v \right)_\Omega
  \,\,\, \forall v \in V ,
\end{equation}

with $\psi, v \in H^1(\Omega)$ and $\Omega$ a rectangle in the $(R, z)$ plane.
Note that this is a diffusion problem with the axisymmetric metric weight
$1/R$, so it is assembled with [MFEMDiffusionKernel.md] and
[MFEMDomainLFKernel.md] without any Grad-Shafranov specific object.

### Solov'ev equilibrium

Following Solov'ev, the free profile functions are chosen so that
$\mu_0 p'$ and $F F'$ are constants, which makes the equation linear. This
example takes $\mu_0 p' = 8$ and $F F' = 2$, so that
$\Delta^* \psi = -8 R^2 - 2$ and the right hand side above is $8 R + 2/R$.

Since $\Delta^* R^4 = 8 R^2$ and $\Delta^* z^2 = 2$, while $R^2$ and $1$ both
lie in the kernel of $\Delta^*$, that choice admits the closed-form solution

\begin{equation}
\psi(R, z) = -R^4 + 2.18 R^2 - 0.8281 - z^2 .
\end{equation}

The two free constants have been chosen to place the $\psi = 0$ separatrix at
$R = 0.7$ and $R = 1.3$ on the midplane. The resulting plasma has its magnetic
axis at $(R, z) = (1.044, 0)$ with $\psi_\mathrm{axis} = 0.36$, a half-height of
$0.6$, and hence an elongation of $\kappa = 2$. The pressure
$p = 8 \psi / \mu_0$ is positive inside the separatrix and peaked on axis.

The computational domain is the box $R \in [0.6, 1.4]$, $z \in [-0.7, 0.7]$,
which strictly contains that plasma, and the analytic $\psi$ is imposed as a
Dirichlet condition on all four sides. This is a fixed-boundary calculation: the
plasma shape is set by the imposed flux rather than by external coils.

### Poloidal field

The poloidal field follows from the solution as
$\vec B_\mathrm{pol} = \frac{1}{R} \vec \nabla \psi \times \hat e_\phi$, so its
magnitude is $\left|\vec \nabla \psi\right| / R$. Because
[MFEMVariable.md] automatically declares a `<name>_grad_mag` coefficient for
scalar $H^1$ variables, this is assembled with an [MFEMParsedFunction.md] and
projected onto a discontinuous auxiliary variable with
[MFEMScalarProjectionAux.md].

### Note on the radial coordinate

This example writes the $1/R$ weight out explicitly as a `ParsedFunction` of
`x` rather than using [MFEMCoordinateTransformations.md]. The coefficients
declared by that object are built on `mfem::CylindricalRadialCoefficient`, which
returns $\sqrt{x^2 + y^2}$, and so gives the cylindrical radius only on a three
dimensional mesh whose symmetry axis is $z$. On a two dimensional $(R, z)$
half-plane, where MOOSE's `x` is the major radius and `y` is the height, it
would instead return $\sqrt{R^2 + z^2}$.

Here the plasma boundary is imposed through the Dirichlet data. For the
counterpart in which it is determined by the solution, see the
[free-boundary example](syntax/MFEM/FreeBoundaryEquilibrium.md).

## Verification

Solving with first order $H^1$ elements on uniformly refined meshes gives the
following $L^2$ errors against the analytic solution, reported by
[MFEML2Error.md]:

| `uniform_refine` | Elements | $\left\Vert \psi - \psi_h \right\Vert_{L^2}$ | Ratio |
| :--- | :--- | :--- | :--- |
| 0 | 28 | 4.4711e-2 | |
| 1 | 112 | 1.1262e-2 | 3.97 |
| 2 | 448 | 2.8208e-3 | 3.99 |
| 3 | 1792 | 7.0552e-4 | 4.00 |

confirming the expected second order convergence. Because the analytic solution
is a quartic polynomial it lies exactly in a fourth order $H^1$ space; solving
with `fec_order = FOURTH` reduces the $L^2$ error to $1.9 \times 10^{-15}$,
which verifies the source term and the $1/R$ weight themselves rather than only
their convergence rate.

## Example File

!listing test/tests/mfem/gradshafranov/gradshafranov.i
