# Integrates inner products of vector variables over the unit cube, with primal
# variables in both an ND space and an [L2]^3 space. The fields u = (1, 2, 3) and
# v = (x, 2y, 3z) are represented exactly, so the integrals are exact:
#
#   (u, v) = 7,   (2 v, u) = 14,   (v, v) = 14/3.

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
    fec_order = FIRST
  []
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = FIRST
  []
[]

[AuxVariables]
  [u]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [v]
    type = MFEMVariable
    fespace = L2VectorFESpace
  []
  [scalar]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[Functions]
  [u_function]
    type = ParsedVectorFunction
    expression_x = '1'
    expression_y = '2'
    expression_z = '3'
  []
  [v_function]
    type = ParsedVectorFunction
    expression_x = 'x'
    expression_y = '2 * y'
    expression_z = '3 * z'
  []
[]

[ICs]
  [u]
    type = MFEMVectorIC
    variable = u
    vector_coefficient = u_function
  []
  [v]
    type = MFEMVectorIC
    variable = v
    vector_coefficient = v_function
  []
[]

[Postprocessors]
  [nd_primal]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = u
    dual_variable = v
  []
  [l2_primal]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = v
    dual_variable = u
    coefficient = 2
  []
  [l2_primal_l2_dual]
    type = MFEMVectorFEInnerProductIntegralPostprocessor
    primal_variable = v
    dual_variable = v
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[Outputs]
  [csv]
    type = CSV
    file_base = OutputData/VectorInnerProduct
  []
[]
