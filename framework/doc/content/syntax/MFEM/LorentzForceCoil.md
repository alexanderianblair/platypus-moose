# Lorentz Force Driven Coil Deformation

## Summary

Solves a transient, two-application problem in which MFEM computes the Lorentz
force density acting on a current carrying coil bar and a libMesh based
[solid mechanics](modules/solid_mechanics/index.md optional=True) parent
application uses that force at every timestep to deform the coil. Both the
electromagnetic and the mechanical response are verified against analytic
solutions.

## Description

The coil is represented by a straight conductor of length $L$, height $h$ and
thickness $w$ occupying $0 \le x \le L$, $0 \le y \le h$, $0 \le z \le w$. It is
energised by a terminal voltage that ramps linearly in time and it sits in a
uniform applied magnetic flux density $\vec B = B_0 \hat e_y$.

### Electromagnetic sub-application

The MFEM sub-application solves the steady current conservation problem for the
electric potential $V$ in the conductor,

\begin{equation}
\begin{split}
-\vec\nabla \cdot \left(\sigma \vec\nabla V\right) = 0 \,\,\,&\mathrm{on}\,\, \Omega \\
V = V_0 t \,\,\, &\mathrm{on}\,\, \Gamma_1 \\
V = 0 \,\,\, &\mathrm{on}\,\, \Gamma_2
\end{split}
\end{equation}

where $\Gamma_1$ and $\Gamma_2$ are the two end faces of the bar and $\sigma$ is
the electrical conductivity. The Dirichlet value on $\Gamma_1$ is a function of
time, so the problem is re-solved at every timestep of the parent application.
The potential is discretised with $H^1$ conforming nodal elements, and its
gradient is projected onto $H(\mathrm{curl})$ conforming Nédélec elements to
obtain the current density

\begin{equation}
\vec J = -\sigma \vec\nabla V
\end{equation}

The Lorentz force density is then formed as the cross product

\begin{equation}
\vec f = \vec J \times \vec B
\end{equation}

and projected onto a vector $L^2$ auxvariable using
[MFEMCrossProductAux.md]. Because transfers of MFEM vector variables to libMesh
are not supported, the three Cartesian components of $\vec f$ are extracted onto
scalar $L^2$ auxvariables with [MFEMInnerProductAux.md] before being handed to
the parent application by a
[MultiAppMFEMTolibMeshShapeEvaluationTransfer.md].

The end faces are the only boundaries carrying current, so the exact potential is
$V = V_0 t \left(1 - x/L\right)$ and both the current density and the Lorentz
force density are uniform over the conductor and ramp linearly in time,

\begin{equation}
\vec J = \frac{\sigma V_0 t}{L} \hat e_x, \qquad
\vec f = \frac{\sigma V_0 B_0 t}{L} \hat e_z
\end{equation}

### Mechanical parent application

The parent application is a transient libMesh solid mechanics problem for the
displacement $\vec u$ of the same coil bar, solved quasi-statically at each
timestep,

\begin{equation}
-\partial_j \sigma_{ij} = f_i
\end{equation}

with the isotropic linear elastic stress/strain relation of
[Linear Elasticity](syntax/MFEM/LinearElasticity.md). The transferred force
components enter the residual through [CoupledForce.md] kernels, one per
displacement component. Second order Lagrange displacements are used so that the
exact solution of the verification case below lies in the finite element space.

## Verification

### Uniaxial strain

In `uniaxial_strain.i` rollers on the four faces normal to $x$ and $y$ suppress
all lateral strain, the $z = 0$ face is held and the $z = w$ face is traction
free. The problem then reduces to the one dimensional equation

\begin{equation}
\left(\lambda + 2\mu\right)\frac{d^2 u_z}{dz^2} + f_z(t) = 0
\end{equation}

whose solution

\begin{equation}
u_z(z, t) = \frac{f_z(t)}{\lambda + 2\mu}\left(w z - \frac{z^2}{2}\right),
\qquad u_x = u_y = 0
\end{equation}

is quadratic in $z$ and therefore lies in the second order Lagrange space. The
computed displacement reproduces it to solver tolerance at every timestep, so
the $L^2$ error postprocessor and the relative error in the displacement of the
traction free surface are both zero to within the solver tolerance.

### Cantilever

In `cantilever.i` the coil bar is clamped at $x = 0$ and free elsewhere, so it
deflects under the transverse Lorentz load. The bar is slender, $L/w = 20$, and
the tip deflection is compared against the Euler-Bernoulli result for a
cantilever under a uniformly distributed load,

\begin{equation}
\delta(t) = \frac{q(t) L^4}{8 E I}, \qquad q(t) = f_z(t)\, h\, w, \qquad
I = \frac{h w^3}{12}
\end{equation}

This is a slender beam approximation rather than an exact solution of the three
dimensional problem, so unlike the uniaxial strain case the agreement is close
rather than exact. On the deliberately coarse mesh used here the computed tip
deflection is 0.98% below the Euler-Bernoulli value; refining the mesh brings the
three dimensional result to roughly 0.5% below it, the remainder being the
finite thickness correction to slender beam theory. Note that the anticlastic
restraint of the fully clamped end stiffens the bar, which outweighs the
softening due to shear deformation at this aspect ratio.

## Example Files

The MFEM sub-application solving for the current density and Lorentz force:

!listing modules/solid_mechanics/test/tests/mfem_lorentz_force/lorentz_force_coil.i

The mechanics parent application used for the uniaxial strain verification:

!listing modules/solid_mechanics/test/tests/mfem_lorentz_force/uniaxial_strain.i

The mechanics parent application used for the cantilever verification:

!listing modules/solid_mechanics/test/tests/mfem_lorentz_force/cantilever.i
