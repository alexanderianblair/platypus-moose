# Mixed nonlinear diffusion problem, coupling two variables through their diffusivities:
#
#   -div((1 + u v) grad(u)) = 0
#   -div((2 + u^2) grad(v)) = 0
#
# Each equation's diffusivity depends on the other variable, so both off-diagonal Jacobian blocks
# carry nonlinear contributions. The v equation's diffusivity does not depend on v, so its diagonal
# block has no contribution beyond the diffusion term itself.

[Mesh]
  type = MFEMMesh
  file = ../mesh/square.e
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
  []
[]

[Variables]
  [u]
    type = MFEMVariable
    fespace = H1FESpace
  []
  [v]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[Functions]
  [u_initial]
    type = ParsedFunction
    expression = '2 * y + 1'
  []
  [v_initial]
    type = ParsedFunction
    expression = '2 - y'
  []
  [k_u]
    type = MFEMParsedFunction
    expression = '1 + u * v'
    symbol_names = 'u v'
    symbol_values = 'u v'
  []
  [dk_u_du]
    type = MFEMParsedFunction
    expression = 'v'
    symbol_names = 'v'
    symbol_values = 'v'
  []
  [dk_u_dv]
    type = MFEMParsedFunction
    expression = 'u'
    symbol_names = 'u'
    symbol_values = 'u'
  []
  [k_v]
    type = MFEMParsedFunction
    expression = '2 + u * u'
    symbol_names = 'u'
    symbol_values = 'u'
  []
  [dk_v_du]
    type = MFEMParsedFunction
    expression = '2 * u'
    symbol_names = 'u'
    symbol_values = 'u'
  []
[]

[ICs]
  [u_ic]
    type = MFEMScalarIC
    coefficient = u_initial
    variable = u
  []
  [v_ic]
    type = MFEMScalarIC
    coefficient = v_initial
    variable = v
  []
[]

[BCs]
  [u_top]
    type = MFEMScalarDirichletBC
    variable = u
    boundary = 'top'
    coefficient = 3.0
  []
  [u_bottom]
    type = MFEMScalarDirichletBC
    variable = u
    boundary = 'bottom'
    coefficient = 1.0
  []
  [v_top]
    type = MFEMScalarDirichletBC
    variable = v
    boundary = 'top'
    coefficient = 1.0
  []
  [v_bottom]
    type = MFEMScalarDirichletBC
    variable = v
    boundary = 'bottom'
    coefficient = 2.0
  []
[]

[Kernels]
  active = 'nl_u nl_v'
  [nl_u]
    type = MFEMNLDiffusionKernel
    variable = u
    k_coefficient = k_u
    dk_du_coefficient = dk_u_du
    coupled_variables = 'v'
    dk_dcoupled_coefficients = 'dk_u_dv'
  []
  [nl_v]
    type = MFEMNLDiffusionKernel
    variable = v
    k_coefficient = k_v
    dk_du_coefficient = 0.0
    coupled_variables = 'u'
    dk_dcoupled_coefficients = 'dk_v_du'
  []
  # The same two weak forms written as bilinear forms with solution-dependent coefficients.
  # Declaring the variables those coefficients depend on assembles them into the nonlinear form,
  # so that they are re-evaluated at each iterate. No derivative of the coefficients is supplied,
  # so Newton is left with a Picard-like Jacobian and converges linearly, but to the same solution.
  [picard_u]
    type = MFEMMixedGradGradKernel
    variable = u
    coefficient = k_u
    coupled_variables = 'u v'
  []
  [picard_v]
    type = MFEMMixedGradGradKernel
    variable = v
    coefficient = k_v
    coupled_variables = 'u'
  []
[]

[Solvers]
  [lin]
    type = MFEMMUMPS
    print_level = 0
  []
  [newton]
    type = MFEMNewtonNonlinearSolver
    # Newton converges quadratically in three iterations with the full Jacobian. The budget is
    # kept tight so that dropping the off-diagonal blocks, which leaves convergence linear, is
    # detected rather than merely slower.
    max_its = 6
    abs_tol = 1.0e-10
    rel_tol = 1.0e-9
    print_level = 1
  []
[]

[VectorPostprocessors]
  [u_sample]
    type = MFEMVariablePointValueSampler
    variable = u
    points = '0.5 0.25 0
              0.5 0.50 0
              0.5 0.75 0'
    execute_on = 'timestep_end'
  []
  [v_sample]
    type = MFEMVariablePointValueSampler
    variable = v
    points = '0.5 0.25 0
              0.5 0.50 0
              0.5 0.75 0'
    execute_on = 'timestep_end'
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
  assembly_level = legacy
[]

[Outputs]
  [csv]
    type = CSV
    file_base = NLMixedDiffusion
    execute_on = 'timestep_end'
  []
[]
