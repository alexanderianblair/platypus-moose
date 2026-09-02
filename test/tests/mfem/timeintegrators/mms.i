# Manufactured solution u = x*cos(t) for du/dt - laplacian(u) = f.
#
# The solution is linear in space, so it lies in the first order H1 space exactly and the
# semi-discrete Galerkin solution reproduces it. The L2 error reported at the end of the run is
# therefore purely the time integration error, and refining dt exposes the order of the scheme.
# The essential data are time dependent, so the boundary treatment of each stage is exercised too.

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
  # f = du/dt - laplacian(u)
  [source]
    type = MFEMParsedFunction
    expression = '-x*sin(t)'
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
  [diff]
    type = MFEMDiffusionKernel
    variable = u
  []
  [source]
    type = MFEMDomainLFKernel
    variable = u
    coefficient = source
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
  [main]
    type = MFEMHypreGMRES
    preconditioner = boomeramg
    l_tol = 1e-14
    l_max_its = 1000
  []
[]

[TimeIntegrators]
  [ti]
    type = MFEMRungeKuttaTimeIntegrator
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
    file_base = mms
    execute_on = 'FINAL'
  []
[]
