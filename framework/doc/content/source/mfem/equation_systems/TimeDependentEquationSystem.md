# TimeDependentEquationSystem

!if! function=hasCapability('mfem')

For transient problems, time derivatives $\dot{u}$ of the trial variables $u$ will also be
present in the weak form.

!equation
\left(\mathcal{T}(\dot{u}), v\right)_{\Omega} + {\left(\mathcal{L}(u), v\right)_{\Omega}
=\left(f,v\right)_{\Omega}\,\,\,\forall v \in V}

Contributions to $\left(\mathcal{T}(\phi_j), \varphi_i\right)_{\Omega}$ are given by
time derivative kernels such as
 [MFEMTimeDerivativeMassKernel.md].

Writing $\mathcal{M}$ for the mass operator contributed by time derivative kernels and
$f(u, t) = f - \mathcal{L}(u)$ for the remaining spatial terms, the semi-discrete system solved is

!equation
\mathcal{M}\dot{u} = f_E(u, t) + f_I(u, t)

subject to essential constraints $u = g(t)$ on constrained degrees of freedom. The spatial terms
are split into those $f_E$ from kernels and integrated boundary conditions with `implicit = false`,
and those $f_I$ from the rest. The split is used only by implicit-explicit (IMEX) schemes; all other
schemes treat $f = f_E + f_I$ as a whole.

The time integration schemes selected in [MFEMTransient.md] advance $u$ through two kinds of stage,
each evaluated at a stage base state $u$ assembled from the solution at the start of the step and
earlier stages.

An implicit stage with stage coefficient $\gamma$ solves for the stage state $w$ satisfying

!equation
\left([\mathcal{M}+\gamma\mathcal{L}](w), v\right)_{\Omega}
=\left([\gamma f + \mathcal{M}(u)],v\right)_{\Omega}\,\,\,\forall v \in V

with $w = g(t)$ on constrained degrees of freedom, so that essential data are imposed exactly on
every stage state. For backwards Euler, $\gamma = \delta t$ and $u = u(t)$.

An explicit stage solves the mass system for the stage slope $k$ satisfying

!equation
\left(\mathcal{M}(k), v\right)_{\Omega}
=\left([f - \mathcal{L}(u)],v\right)_{\Omega}\,\,\,\forall v \in V

with $k$ on constrained degrees of freedom set to the rate of change of the essential data,
approximated by fourth order central differences in time. In the explicit part of an IMEX scheme, $k$ is instead
zero on constrained degrees of freedom, since the essential data are imposed by its implicit stages.

At the end of each step, the essential data at the new time are imposed on the constrained degrees
of freedom, which schemes that are not stiffly accurate would otherwise only approximate.

!if-end!

!else
!include mfem/mfem_warning.md
