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

Transient problems in MOOSE's wrapping of MFEM are advanced with an explicit or diagonally
implicit Runge-Kutta scheme, selected in the [TimeIntegrators](syntax/TimeIntegrators/index.md)
block and defaulting to backwards Euler.

An $s$ stage scheme with Butcher tableau $(a, b, c)$ advances the state from $u_n = u(t)$ to
$u_{n+1} = u_n + \delta t \sum_i b_i \dot{u}_i$ through stage states $u^{(i)}$ evaluated at the
stage times $t + c_i \delta t$. Denoting the base state of stage $i$ by
$y^{(i)} = u_n + \delta t \sum_{j<i} a_{ij} \dot{u}_j$, a stage with a nonzero diagonal
coefficient $\gamma = a_{ii} \delta t$ is solved for the stage state,

!equation
\left([\mathcal{T}+\gamma\mathcal{L}](u^{(i)}), v\right)_{\Omega}
=\left([\gamma f + \mathcal{T}(y^{(i)})],v\right)_{\Omega}\,\,\,\forall v \in V

from which the stage slope follows as
$\dot{u}_i = (u^{(i)} - y^{(i)})/\gamma$. Backwards Euler is the single stage case
$a_{11} = b_1 = c_1 = 1$, for which $y^{(1)} = u_n$ and $u_{n+1} = u^{(1)}$.

Solving for the stage state rather than the stage slope lets essential boundary conditions be
imposed directly, as values of $u^{(i)}$ at the stage time.

A stage with a vanishing diagonal coefficient is explicit, and is evaluated by solving the mass
system

!equation
\left(\mathcal{T}(\dot{u}_i), v\right)_{\Omega}
=\left([f - \mathcal{L}(y^{(i)})],v\right)_{\Omega}\,\,\,\forall v \in V

for the stage slope directly. The unknown of such a stage is a time derivative, so the values
imposed on essentially constrained degrees of freedom are the time derivative of the essential
data, obtained by central differencing it about the stage time.

!if-end!

!else
!include mfem/mfem_warning.md
