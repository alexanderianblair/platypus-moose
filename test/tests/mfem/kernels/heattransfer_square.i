[Mesh]
  type = MFEMMesh
  file = square.e
  dim = 2
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
  [temperature]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[AuxVariables]
  inactive = average_temperature
  [average_temperature]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[AuxKernels]
  inactive = average_field
  [average_field]
    type = MFEMScalarTimeAverageAux
    variable = average_temperature
    source = temperature
  []
[]

[Kernels]
  [diff]
    type = MFEMDiffusionKernel
    variable = temperature
  []
  [dT_dt]
    type = MFEMTimeDerivativeMassKernel
    variable = temperature
  []
[]

[BCs]
  active = 'bottom top_dirichlet'
  [bottom]
    type = MFEMScalarDirichletBC
    variable = temperature
    boundary = 'bottom'
    coefficient = 1.0
  []
  [top_convective]
    type = MFEMConvectiveHeatFluxBC
    variable = temperature
    boundary = '2'
    T_infinity = .5
    heat_transfer_coefficient = 5
  []
  [top_dirichlet]
    type = MFEMScalarDirichletBC
    variable = temperature
    boundary = 'top'
  []
[]

[Preconditioner]
  [boomeramg]
    type = MFEMHypreBoomerAMG
    print_level = 0
  []
  [jacobi]
    type = MFEMOperatorJacobiSmoother
  []
[]

[Solver]
  type = MFEMHypreGMRES
  preconditioner = boomeramg
  print_level = 1
  l_tol = 1e-16
  l_max_its = 1000
[]

[Executioner]
  type = MFEMTransient
  device = cpu
  assembly_level = legacy
  dt = 0.25
  start_time = 0.0
  end_time = 0.5
  print_level = 1
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/HeatTransferSquare
    vtk_format = ASCII
  []
[]
