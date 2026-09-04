# Free-Boundary Grad-Shafranov Problem

## Summary

Solves the Grad-Shafranov equation for an axisymmetric plasma equilibrium in
which the plasma boundary is not prescribed. The only boundary condition is
$\psi = 0$ on the far edge of a computational box that also contains the
poloidal field coils, and the plasma shape follows from the balance between the
plasma current and the coil currents. This is the free-boundary counterpart of
the fixed-boundary [Solov'ev equilibrium](syntax/MFEM/GradShafranov.md).

## Description

The weak form is the one derived for the [fixed-boundary case](syntax/MFEM/GradShafranov.md),
written in terms of the toroidal current density $J_\phi$ rather than the
pressure and poloidal current profiles separately:

\begin{equation}
\left( \frac{1}{R} \vec \nabla \psi, \vec \nabla v \right)_\Omega
  = \left( \mu_0 J_\phi, v \right)_\Omega \,\,\, \forall v \in V ,
\qquad \psi = 0 \,\, \mathrm{on} \,\, \partial\Omega ,
\end{equation}

with $\psi, v \in H^1(\Omega)$. What makes the problem free-boundary is the
current density. Inside the plasma it is taken to be the Luxon-Brown profile

\begin{equation}
J_\phi(R, \psi) = \lambda
  \left[ \beta_0 \frac{R}{R_0} + (1 - \beta_0) \frac{R_0}{R} \right]
  \left( 1 - \psi_N^2 \right)^2 ,
\qquad
\psi_N = \frac{\psi_\mathrm{axis} - \psi}{\psi_\mathrm{axis} - \psi_\mathrm{boundary}} ,
\end{equation}

and zero outside it, where the bracket interpolates between the $R$ and $1/R$
weightings that the $\mu_0 R^2 p'(\psi)$ and $F F'(\psi)$ terms of the
Grad-Shafranov source carry. The normalised flux $\psi_N$ runs from $0$ on the
magnetic axis to $1$ on the plasma boundary, so clamping $(1 - \psi_N^2)$ at
zero is what confines the current to the closed flux surfaces.

Both normalising fluxes are measured from the solution itself with
[MFEMVariableExtremeValue.md]:

- $\psi_\mathrm{axis}$, the flux on the magnetic axis, is the largest flux
  anywhere inside the vessel;
- $\psi_\mathrm{boundary}$, the flux on the plasma boundary, is the largest flux
  on the limiter. It is the last flux surface that closes without striking the
  wall, so the plasma is bounded by wherever it makes contact.

The equation is therefore nonlinear: where current flows at all depends on the
solution. Because MFEM-MOOSE rebuilds the linear forms on every implicit solve,
adding an [MFEMTimeDerivativeMassKernel.md] and stepping in pseudo-time
evaluates the source at the previous iterate, which is a Picard iteration with
the step size $\Delta t$ setting how strongly each update is damped. The
pseudo-time derivative vanishes at convergence, so the fixed point reached is a
solution of the Grad-Shafranov system itself.

### Geometry

The computational box is $R \in [0.1, 2.5]$, $z \in [-1.5, 1.5]$ in metres,
meshed with named subdomains for

- `plasma_region`, $R \in [0.6, 1.6]$, $z \in [-0.7, 0.7]$, where the plasma
  current may flow;
- `limiter`, the surrounding frame out to $R \in [0.5, 1.7]$, $z \in [-0.8, 0.8]$,
  which bounds the plasma;
- `vertical_field_coils`, a pair at $R \in [1.9, 2.1]$, $|z| \in [0.4, 0.6]$,
  each carrying $-0.48$ MA. Their current runs anti-parallel to the plasma
  current, so the vertical field they produce at the plasma, crossed into the
  plasma current, supplies the inward force that balances the hoop force. A
  tokamak plasma has no equilibrium without them: with these coils switched off
  the Picard iteration collapses to zero current;
- `shaping_coils`, a pair at $R \in [0.9, 1.3]$, $|z| \in [1.0, 1.2]$, each
  carrying $+0.25$ MA parallel to the plasma current, which pulls the flux
  surfaces vertically and elongates the plasma.

The coil currents are prescribed. Codes that solve the inverse problem instead
choose them to produce a requested plasma shape.

### Result

Starting from a seed flux that only has to give the first iteration a well
defined axis and boundary, the iteration converges to

| Quantity | Value |
| :--- | :--- |
| Magnetic axis flux $\psi_\mathrm{axis}$ | 0.38353 Wb |
| Plasma boundary flux $\psi_\mathrm{boundary}$ | 0.18566 Wb |
| Total plasma current $I_p$ | 0.9916 MA |
| Midplane plasma extent | $R \in [0.87, 1.60]$ m |
| Magnetic axis position | $R = 1.27$ m |
| Elongation $\kappa$ | 1.37 |

The plasma is limited on the outboard side, touching the limiter at $R = 1.6$ m,
and its magnetic axis sits outboard of the geometric centre of the plasma at
$R = 1.24$ m. That outward displacement is the Shafranov shift, and like the
plasma boundary it is an output of the solve rather than an input.

## Verification

The pseudo-time derivative damps each Picard update but is absent from the
converged equations, so the equilibrium reached must not depend on the step
size. Running with $\Delta t = 1$ and $\Delta t = 5$ gives

| $\Delta t$ | $\psi_\mathrm{axis}$ | $\psi_\mathrm{boundary}$ |
| :--- | :--- | :--- |
| 1 | 0.38353321245824 | 0.18565655682035 |
| 5 | 0.38353321246051 | 0.18565655681891 |

agreeing to eleven significant figures, and both test entries are run against
gold files that record this. The equilibrium is likewise unchanged to the
digits printed when the calculation is partitioned over 1, 2 or 4 MPI ranks.

## Limitations

Two simplifications are worth naming.

- The plasma boundary is taken to be the limiter contact. A diverted plasma is
  instead bounded by the flux at an X-point, which would need a saddle point
  search over $\psi$ rather than a maximum over the limiter subdomain, so the
  shaping coil current here is kept below the value at which an X-point enters
  the vessel.
- The current scale $\lambda$ is prescribed and the resulting $I_p$ reported.
  Codes normally do the reverse, rescaling $\lambda$ each iteration so that
  $\int J_\phi \, \mathrm{d}R \, \mathrm{d}z$ matches a requested plasma current.

## Example File

!listing test/tests/mfem/gradshafranov/freeboundary.i
