# Checks that an auxkernel reading an MFEM variable through a coefficient parameter
# executes after the auxkernel that updates that variable on the same execute flag.
#
# e_field = (t, 0, 0) is projected at TIMESTEP_END, and lorentz_force = e_field x b
# with b = (0, 1, 0) reads it through the coefficient parameter first_source_vec.
# The cross product auxkernel is listed first, so only the dependency on e_field
# orders it correctly; read in the wrong order it would lag one step behind.

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/ref-cube.mesh
[]

[Problem]
  type = MFEMProblem
  solve = false
[]

[FESpaces]
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
  []
  [L2VectorFESpace]
    type = MFEMVectorFESpace
    fec_type = L2
    fec_order = CONSTANT
  []
[]

[AuxVariables]
  [e_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [lorentz_force]
    type = MFEMVariable
    fespace = L2VectorFESpace
  []
[]

[Functions]
  [e_function]
    type = ParsedVectorFunction
    expression_x = 't'
    expression_y = '0'
    expression_z = '0'
  []
  [expected_force]
    type = ParsedVectorFunction
    expression_x = '0'
    expression_y = '0'
    expression_z = 't'
  []
[]

[AuxKernels]
  [cross]
    type = MFEMCrossProductAux
    variable = lorentz_force
    first_source_vec = e_field
    second_source_vec = '0 1 0'
    execute_on = TIMESTEP_END
  []
  [e_field]
    type = MFEMVectorProjectionAux
    variable = e_field
    vector_coefficient = e_function
    execute_on = TIMESTEP_END
  []
[]

[Postprocessors]
  [force_error]
    type = MFEMVectorL2Error
    variable = lorentz_force
    function = expected_force
    execute_on = TIMESTEP_END
  []
[]

[Executioner]
  type = MFEMTransient
  device = cpu
  dt = 1
  num_steps = 2
[]

[Outputs]
  [csv]
    type = CSV
    file_base = OutputData/CoefficientDependency
  []
[]
