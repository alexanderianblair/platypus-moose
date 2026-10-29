# Monolithic electro-thermal problem in the H formulation: eddy currents decay through a conductor
# whose resistivity depends on both the current density and the temperature, heating it through
# Joule losses. The magnetic field and temperature are solved for in one equation system, so the
# Jacobian includes the off-diagonal couplings between them.

rho0 = 1.0
alpha = 0.5
e_c = 1.0
j_c = 1.0
heat_capacity = 1.0
conductivity = 1.0

[Mesh]
  [gen]
    type = MFEMGeneratedMeshGenerator
    dim = 3
    nx = 2
    ny = 2
    nz = 2
  []
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
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
  []
[]

[Variables]
  [h_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [temperature]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[Functions]
  [initial_h_field]
    type = ParsedVectorFunction
    expression_x = 'sin(pi * y)'
    expression_y = 'sin(pi * z)'
    expression_z = 'sin(pi * x)'
  []
  [initial_temperature]
    type = ParsedFunction
    expression = '1 + x'
  []
  # Ohmic resistivity increasing linearly with temperature, in parallel with a power law of
  # index 3 in the current density
  [resistivity]
    type = MFEMParsedFunction
    expression = 'rho0 * (1 + alpha * T) + (e_c / j_c) * (j / j_c)^2'
    symbol_names = 'j T rho0 alpha e_c j_c'
    symbol_values = 'h_field_curl_mag temperature ${rho0} ${alpha} ${e_c} ${j_c}'
  []
  [j_dresistivity_dj]
    type = MFEMParsedFunction
    expression = '2 * (e_c / j_c) * (j / j_c)^2'
    symbol_names = 'j e_c j_c'
    symbol_values = 'h_field_curl_mag ${e_c} ${j_c}'
  []
[]

[ICs]
  [h_field_ic]
    type = MFEMVectorIC
    variable = h_field
    vector_coefficient = initial_h_field
  []
  [temperature_ic]
    type = MFEMScalarIC
    variable = temperature
    coefficient = initial_temperature
  []
[]

[Kernels]
  [dBdt]
    type = MFEMTimeDerivativeVectorFEMassKernel
    variable = h_field
  []
  [curlE]
    type = MFEMNLCurlCurlKernel
    variable = h_field
    k_coefficient = resistivity
    curlu_dk_dcurlu_coefficient = j_dresistivity_dj
    coupled_variables = temperature
    dk_dcoupled_coefficients = ${fparse rho0 * alpha}
  []
  [dTdt]
    type = MFEMTimeDerivativeMassKernel
    variable = temperature
    coefficient = ${heat_capacity}
  []
  [conduction]
    type = MFEMDiffusionKernel
    variable = temperature
    coefficient = ${conductivity}
  []
  [joule_heating]
    type = MFEMNLJouleHeatingKernel
    variable = temperature
    curl_variable = h_field
    k_coefficient = resistivity
    curlu_dk_dcurlu_coefficient = j_dresistivity_dj
    coupled_variables = temperature
    dk_dcoupled_coefficients = ${fparse rho0 * alpha}
  []
[]

[Solvers]
  active = 'newton mumps'
  [newton]
    type = MFEMNewtonNonlinearSolver
    max_its = 10
    abs_tol = 1e-12
    rel_tol = 1e-10
    print_level = 1
  []
  [mumps]
    type = MFEMMUMPS
    print_level = 0
  []
  [petsc_nl]
    type = MFEMPetscNonlinearSolver
    max_its = 10
    abs_tol = 1e-12
    rel_tol = 1e-10
    print_level = 1
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
    petsc_options_value = 'lu mumps'
  []
[]

[Executioner]
  type = MFEMTransient
  dt = 0.25
  num_steps = 2
[]

[Postprocessors]
  [average_temperature]
    type = MFEMElementAverageValue
    variable = temperature
  []
  [magnetic_energy]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    coefficient = 0.5
    primal_variable = h_field
    dual_variable = h_field
  []
[]

[Outputs]
  csv = true
[]
