# Manufactured solution u = x*cos(t) for the nonlinear diffusion problem
# du/dt - d/dx(k(u) du/dx) = f with k(u) = 1 + u^2.
#
# As in mms.i the solution is linear in space, so the semi-discrete solution reproduces it and
# the L2 error reported at the end of the run is purely the time integration error. The essential
# data are time dependent, and every implicit stage is a nonlinear solve in its own right at its
# own stage time; an explicit stage instead evaluates the nonlinear residual at the stage base
# state and moves it to the right hand side of a linear mass solve.

[Mesh]
  [gen]
    type = MFEMGeneratedMeshGenerator
    dim = 1
    nx = 8
  []
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [h1]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
  []
[]

[Variables]
  [u]
    type = MFEMVariable
    fespace = h1
  []
[]

[Functions]
  [exact]
    type = MFEMParsedFunction
    expression = 'x*cos(t)'
  []
  # f = du/dt - d/dx(k(u) du/dx) for u = x*cos(t) and k(u) = 1 + u^2
  [source]
    type = MFEMParsedFunction
    expression = '-x*sin(t) - 2*x*cos(t)^3'
  []
  [k]
    type = MFEMParsedFunction
    expression = '1 + u*u'
    symbol_names = 'u'
    symbol_values = 'u'
  []
  [dk_du]
    type = MFEMParsedFunction
    expression = '2*u'
    symbol_names = 'u'
    symbol_values = 'u'
  []
[]

[ICs]
  [u0]
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
  [source]
    type = MFEMDomainLFKernel
    variable = u
    coefficient = source
  []
  [nldiff]
    type = MFEMNLDiffusionKernel
    variable = u
    k_coefficient = k
    dk_du_coefficient = dk_du
  []
[]

[BCs]
  [ends]
    type = MFEMScalarDirichletBC
    variable = u
    boundary = 'left right'
    coefficient = exact
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
  []
  [lin]
    type = MFEMHyprePCG
    preconditioner = boomeramg
    l_tol = 1e-12
    l_max_its = 1000
  []
  [newton]
    type = MFEMNewtonNonlinearSolver
    max_its = 100
    abs_tol = 1.0e-12
    rel_tol = 1.0e-10
  []
[]

[TimeIntegrators]
  [ti]
    type = MFEMRungeKuttaTimeIntegrator
    scheme = sdirk22
  []
[]

[Executioner]
  type = MFEMTransient
  device = cpu
  dt = 0.25
  start_time = 0.0
  end_time = 1.0
[]

[Postprocessors]
  [error]
    type = MFEML2Error
    variable = u
    function = exact
    execute_on = 'FINAL'
  []
[]

[Outputs]
  [out]
    type = CSV
    file_base = nonlinear
    execute_on = 'FINAL'
  []
[]
