# Solves du/dt - k d^2u/dx^2 + c u = s on (0, 1), with a time-dependent Dirichlet condition at
# x = 0 and a convective condition -k du/dx = h (u - T_inf) at x = 1 whose far-field data depend
# on time. The exact solution is quadratic in x, so is represented exactly by second order
# elements and the error measured is that of the time integration alone. The diffusivity is small
# enough for explicit schemes to be stable at the timesteps used.
k = 0.05
c = 1
h = 1

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/ref-segment.mesh
  uniform_refine = 1
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = SECOND
  []
[]

[Variables]
  [u]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[Functions]
  [exact]
    type = ParsedFunction
    expression = 'cos(t) * (1 + x^2) + x * sin(t)'
  []
  [source]
    type = ParsedFunction
    expression = '-sin(t) * (1 + x^2) + x * cos(t) - 2 * k * cos(t) + c * (cos(t) * (1 + x^2) + x * sin(t))'
    symbol_names = 'k c'
    symbol_values = '${k} ${c}'
  []
  [T_inf]
    # u(1, t) + k du/dx(1, t) / h
    type = ParsedFunction
    expression = '(2 * cos(t) + sin(t)) * (1 + k / h)'
    symbol_names = 'k h'
    symbol_values = '${k} ${h}'
  []
[]

[ICs]
  [u_ic]
    type = MFEMScalarIC
    variable = u
    coefficient = exact
  []
[]

[Kernels]
  [dudt]
    type = MFEMTimeDerivativeMassKernel
    variable = u
  []
  [diffusion]
    type = MFEMDiffusionKernel
    variable = u
    coefficient = ${k}
  []
  [reaction]
    type = MFEMMassKernel
    variable = u
    coefficient = ${c}
  []
  [source]
    type = MFEMDomainLFKernel
    variable = u
    coefficient = source
  []
[]

[BCs]
  [left]
    type = MFEMScalarDirichletBC
    variable = u
    boundary = 1
    coefficient = exact
  []
  [right]
    type = MFEMConvectiveHeatFluxBC
    variable = u
    boundary = 2
    T_infinity = T_inf
    heat_transfer_coefficient = ${h}
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
  []
  [main]
    type = MFEMHyprePCG
    preconditioner = boomeramg
    l_tol = 1e-16
  []
[]

[Executioner]
  type = MFEMTransient
  device = cpu
  dt = 0.1
  end_time = 1
[]

[Postprocessors]
  [l2_error]
    type = MFEML2Error
    variable = u
    function = exact
  []
[]

[Outputs]
  csv = true
[]
