# Definite Maxwell problem solved with Nedelec elements of the first kind
# based on MFEM Example 3.

[Mesh]
  type = MFEMMesh
  file = ../mesh/coil.gen
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [HCurlFESpace]
    type = MFEMGenericFESpace
    fec_name = ND@Gi_3D_P1
  []
  [HDivFESpace]
    type = MFEMGenericFESpace
    fec_name = RT@Gi_3D_P0
  []
[]

[Variables]
  [a_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [current_density]
    type = MFEMVariable
    fespace = HCurlFESpace
  []  
[]

[AuxVariables]
  [b_field]
    type = MFEMVariable
    fespace = HDivFESpace
  []
[]

[AuxKernels]
  [curl]
    type = MFEMCurlAux
    variable = b_field
    source = a_field
    scale_factor = 1.0
    execute_on = TIMESTEP_END
  []
[]

[Functions]
  [exact_a_field]
    type = ParsedVectorFunction
    expression_x = '0'
    expression_y = '0'
    expression_z = '0'

    symbol_names = kappa
    symbol_values = 3.1415926535
  []
[]

[BCs]
  [tangential_a_bdr]
    type = MFEMVectorFunctorTangentialDirichletBC
    variable = a_field
    vector_coefficient = exact_a_field
    boundary = '1 2 4'
  []
[]

[FunctorMaterials]
  [Vacuum]
    type = MFEMGenericConstantFunctorMaterial
    prop_names = reluctivity
    prop_values = 1.0
  []
[]

[Kernels]
  [curlcurl]
    type = MFEMCurlCurlKernel
    variable = a_field
    coefficient = reluctivity
  []
  [source]
    type = MFEMVectorFEDomainLFKernel
    variable = a_field
    vector_coefficient = current_density
    block = '1'
  []
[]

[Preconditioner]
  [ams]
    type = MFEMHypreAMS
    fespace = HCurlFESpace
    # singular = true
  []
[]

[Solver]
  type = MFEMHypreGMRES
  preconditioner = ams
  l_tol = 1e-9
  l_max_its = 1000
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[MultiApps]
  [subapp]
    type = FullSolveMultiApp
    input_files = coil_open_coil_source.i
    execute_on = INITIAL
  []
[]

[Transfers]
  [from_sub]
    type = MultiAppMFEMCopyTransfer
    source_variable = current_density
    variable = current_density
    from_multi_app = subapp
  []
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/CurlCurl
    vtk_format = ASCII
  []
[]
