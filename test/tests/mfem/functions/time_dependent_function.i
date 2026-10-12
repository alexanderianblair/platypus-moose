# Projects the function f(t) = t onto an L2 variable without solving any equation
# system, so that the projection only sees the current time if MFEMProblem passes it
# to the coefficients. The start time is nonzero so that INITIAL is checked too.

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/ref-cube.mesh
[]

[Problem]
  type = MFEMProblem
  solve = false
[]

[FESpaces]
  [L2FESpace]
    type = MFEMScalarFESpace
    fec_type = L2
    fec_order = CONSTANT
  []
[]

[AuxVariables]
  [time_field]
    type = MFEMVariable
    fespace = L2FESpace
  []
[]

[Functions]
  [time_function]
    type = ParsedFunction
    expression = 't'
  []
[]

[AuxKernels]
  [project]
    type = MFEMScalarProjectionAux
    variable = time_field
    coefficient = time_function
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Postprocessors]
  [projected_time]
    type = MFEMElementAverageValue
    variable = time_field
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = MFEMTransient
  device = cpu
  start_time = 1
  dt = 0.5
  num_steps = 2
[]

[Outputs]
  [csv]
    type = CSV
    file_base = OutputData/TimeDependentFunction
  []
[]
