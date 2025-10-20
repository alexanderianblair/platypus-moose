# Magnetostatic problem solved on a closed conductor subject to
# global loop voltage constraint.

[Mesh]
  type = MFEMMesh
  file = ../mesh/stacked_embedded_concentric_torus.e
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
  []
  [HDivFESpace]
    type = MFEMVectorFESpace
    fec_type = RT
    fec_order = CONSTANT
  []
[]

[Variables]
  [a_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[AuxVariables]
  [b_field]
    type = MFEMVariable
    fespace = HDivFESpace
  []
  [e_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [top_coil_e_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [bottom_coil_e_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
[]

[AuxKernels]
  [update_total_coil_e_field]
    type = MFEMSumAux
    variable = e_field
    source_variables = 'top_coil_e_field bottom_coil_e_field'
    scale_factors = '1.0 1.0'
    execute_on = INITIAL
    execution_order_group = 2
  []
  [curl]
    type = MFEMCurlAux
    variable = b_field
    source = a_field
    scale_factor = 1.0
    execute_on = TIMESTEP_END
    execution_order_group = 3
  []
[]

[Functions]
  [exact_a_field]
    type = ParsedVectorFunction
    expression_x = '0'
    expression_y = '0'
    expression_z = '0'
  []
[]

[BCs]
  [tangential_a_bdr]
    type = MFEMVectorTangentialDirichletBC
    variable = a_field
    vector_coefficient = exact_a_field
    boundary = 'Exterior'
  []
[]

[FunctorMaterials]
  [Vacuum]
    type = MFEMGenericFunctorMaterial
    prop_names = reluctivity
    prop_values = '${fparse (1.0e7)/(4*pi)}'
  []
  [Conductor]
    type = MFEMGenericFunctorMaterial
    prop_names = conductivity
    prop_values = 5.96e7
    block = 'TorusCore1 TorusSheath1 TorusCore2 TorusSheath2'
  []
[]

[Kernels]
  [mass]
    type = MFEMVectorFEMassKernel
    variable = a_field
    coefficient = 1.0
  []
  [curlcurl]
    type = MFEMCurlCurlKernel
    variable = a_field
    coefficient = reluctivity
  []
  [source]
    type = MFEMMixedVectorMassKernel
    variable = a_field
    trial_variable = e_field
    coefficient = conductivity
    block = 'TorusCore1 TorusSheath1 TorusCore2 TorusSheath2'
  []
[]

[Preconditioner]
  [ams]
    type = MFEMHypreAMS
    fespace = HCurlFESpace
  []
[]

[Solver]
  type = MFEMHyprePCG
  preconditioner = ams
  l_tol = 1e-14
  l_max_its = 1000
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[MultiApps]
  [top_coil]
    type = FullSolveMultiApp
    input_files = stacked_cut_closed_coil_1.i
    execute_on = INITIAL
  []
  [bottom_coil]
    type = FullSolveMultiApp
    input_files = stacked_cut_closed_coil_2.i
    execute_on = INITIAL
  []
[]

[Transfers]
  [from_top_coil]
    type = MultiAppMFEMCopyTransfer
    source_variable = e_field
    variable = top_coil_e_field
    from_multi_app = top_coil
  []
  [from_bottom_coil]
    type = MultiAppMFEMCopyTransfer
    source_variable = e_field
    variable = bottom_coil_e_field
    from_multi_app = bottom_coil
  []
[]

# [Postprocessors]
#   [TotalCurrent]
#     type = MFEMBoundaryNetFluxPostprocessor
#     variable = e_field
#     boundary = 'MeasurementPlane'
#   []
# []

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/MagnetostaticStackedClosedCoilsClosedCoil
    vtk_format = ASCII
  []
  [ReportedCurrent]
    type = CSV
    file_base = OutputData/MagnetostaticStackedClosedCoilsCSV
  []
[]
