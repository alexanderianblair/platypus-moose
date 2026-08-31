# EquationSystem

!if! function=hasCapability('mfem')

The EquationSystem is responsible for defining and assembling the weak form of the PDE into an
[`mfem::Operator`](https://docs.mfem.org/html/classmfem_1_1Operator.html) used to solve an iteration
of the FE problem. This operator is passed to a linear or nonlinear solver as appropriate in a
[ProblemOperator.md], which handles the update of the
state of all variables (including any required nonlinear iterations).

!equation
{\left(\mathcal{L}(u), v\right)_{\Omega}=\left(f,v\right)_{\Omega}\,\,\,\forall v \in V}

Discretizing the trial variable $u$ by $u\approx u_h = \sum_i u_i \phi_i(\vec x)$,
and approximating the test space $V$ by a finite dimensional subspace spanned by
the basis $\{\varphi_i\}$, the weak form becomes

!equation
{\sum_{j}u_j\left(\mathcal{L}(\phi_j),\varphi_i\right)_{\Omega}=\left(f,\varphi_i\right)_{\Omega}}

after variation over the test variable $v$. This can be expressed in matrix form by
defining $A_{ij} = \left(\mathcal{L}(\phi_j), \varphi_i\right)_{\Omega}$ and
${b_i=\left(f,\varphi_i\right)_{\Omega}}$.

!equation
{\sum_j A_{ij} u_j = b_i}

[MFEMKernels](source/mfem/kernels/MFEMKernel.md) contribute domain integrators to the weak form, and
[MFEMIntegratedBCs](source/mfem/bcs/MFEMIntegratedBC.md) contribute boundary integrators to the weak
form. [`mfem::BilinearFormIntegrators`](https://mfem.org/bilininteg/) add contributions to
$A_{ij}(\varphi_i, \phi_j)$ and [`mfem::LinearFormIntegrators`](https://mfem.org/lininteg/) add
contributions to $b_i(\varphi_i)$ when assembled.

[MFEMKernels](source/mfem/kernels/MFEMKernel.md) can also contribute to domain integrators for
non-linear actions. This allows to form the residual $\mathcal{L}(u)$ for non-linear Newton's
method as shown below

!equation
{\mathbf{J}\left(\vec{u}_n\right) \delta \vec{u}_{n+1}=-\vec{R}\left(\vec{u}_n\right)}

!equation
{\vec{u}_{n+1}=\vec{u}_n+\delta \vec{u}_{n+1}}

where $\mathbf{J}$ is the Jacobian, and $\delta \vec{u}$ is the incremental solution.

## Systems of several variables

When several variables are solved for together, the residual and the Jacobian are block
structured, with one row per test variable and one column per trial variable,

!equation
{\sum_j \mathbf{J}_{ij}\left(\vec{u}_n\right)\delta \vec{u}_{n+1,j} = -\vec{R}_i\left(\vec{u}_n\right)}

The blocks $\mathbf{J}_{ij}$ receive contributions from two places. Kernels and integrated boundary
conditions whose weak form does not depend on the solution are assembled once per solve, into a
bilinear form for $i = j$ and a mixed bilinear form otherwise. Everything else is assembled into a
single [`mfem::ParBlockNonlinearForm`](https://docs.mfem.org/html/classmfem_1_1ParBlockNonlinearForm.html)
spanning all trial variables, which supplies both the nonlinear part of $\vec{R}$ and every
$\partial\vec{R}_i/\partial\vec{u}_j$, at the current iterate. Off-diagonal blocks may therefore
carry nonlinear terms.

A kernel contributes to the nonlinear form by supplying one of

- a `mfem::NonlinearFormIntegrator` from `createNLIntegrator()`, whose action on the DoFs of the
  kernel's own test variable gives the residual, and whose gradient gives the $(i,i)$ block; or
- a `mfem::BilinearFormIntegrator` from `createNLMixedIntegrator()`, whose action on the DoFs of
  the kernel's trial variable gives a residual that is linear in those DoFs but whose coefficients
  depend on the solution, and which gives the corresponding $(i,j)$ block.

The second of these is what routes a kernel with solution-dependent coefficients away from the
bilinear forms. This matters for correctness rather than for the convergence rate: a bilinear form
is assembled once per solve, so the coefficients of a term left there are held at the values they
took when the system was formed, and Newton converges to the solution of that frozen problem
instead of the intended one.

Kernels supplying either of the above may additionally name, through `getCoupledVariableNames()`,
the variables their coefficients depend on, and return a `mfem::BilinearFormIntegrator` from
`createOffDiagJacobianIntegrator()` giving the derivative of their residual with respect to each.
Those integrators populate the remaining off-diagonal blocks. Supplying them is optional: without
them the Jacobian is approximate and Newton converges more slowly, but to the same solution.

All variables solved for by a nonlinear equation system must be defined on the same mesh.

!if-end!

!else
!include mfem/mfem_warning.md
