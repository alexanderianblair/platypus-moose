# Read the H^1 basis cochain Gmsh computed on the air region surrounding a toroidal
# conductor, and check that the resulting Nedelec field is curl-free there.

[Mesh]
  type = MFEMMesh
  file = ../mesh/torus_in_box.msh
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
  [HDivFESpace]
    type = MFEMVectorFESpace
    fec_type = RT
    fec_order = CONSTANT
  []
[]

[AuxVariables]
  [cut_field]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [cut_curl]
    type = MFEMVariable
    fespace = HDivFESpace
  []
[]

[ICs]
  [cut]
    type = MFEMCohomologyCutIC
    variable = cut_field
    cut_name = 'H^1{1}1'
    amplitude = 1.0
  []
[]

[AuxKernels]
  [curl]
    type = MFEMCurlAux
    variable = cut_curl
    source = cut_field
    execute_on = 'INITIAL'
  []
[]

[Postprocessors]
  [FieldSquaredInAir]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = cut_field
    dual_variable = cut_field
    block = Air
    execute_on = 'INITIAL'
  []
  [CurlSquaredInAir]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = cut_curl
    dual_variable = cut_curl
    block = Air
    execute_on = 'INITIAL'
  []
  [CurlSquaredInCoil]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = cut_curl
    dual_variable = cut_curl
    block = Coil
    execute_on = 'INITIAL'
  []
[]

[Executioner]
  type = MFEMSteady
[]

[Outputs]
  [ReportedPostprocessors]
    type = CSV
    file_base = OutputData/CohomologyCutSource
  []
[]
