# Isogeometric analysis (IGA) of linear elasticity on a three-dimensional U-spline mesh of a block
# with a cylindrical hole, with a vector-valued U-spline space for the displacement. The
# displacement is prescribed on the whole boundary to be the linear field u, whose constant stress
# satisfies the equilibrium equations without body forces. The U-spline space contains every
# linear field, so the error is that of the Gauss quadrature of the rational basis functions alone.

[Mesh]
  type = MFEMUSplineMesh
  file = ../mesh/BlockWithHole_Patch9.bxt
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [USplineFESpace]
    type = MFEMUSplineFESpace
    vdim = 3
  []
[]

[Variables]
  [displacement]
    type = MFEMVariable
    fespace = USplineFESpace
  []
[]

[Functions]
  [exact_displacement]
    type = ParsedVectorFunction
    expression_x = '0.01 + 0.02*x - 0.01*y + 0.03*z'
    expression_y = '-0.02 + 0.01*x + 0.04*y - 0.02*z'
    expression_z = '0.03 - 0.03*x + 0.02*y + 0.01*z'
  []
[]

[BCs]
  [boundary]
    type = MFEMVectorDirichletBC
    variable = displacement
    vector_coefficient = exact_displacement
  []
[]

[Kernels]
  [elasticity]
    type = MFEMLinearElasticityKernel
    variable = displacement
    lambda = 2.0
    mu = 1.0
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
  []
  [main]
    type = MFEMHyprePCG
    preconditioner = boomeramg
    l_tol = 1e-16
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[Postprocessors]
  [l2_error]
    type = MFEMVectorL2Error
    variable = displacement
    function = exact_displacement
  []
[]

[Outputs]
  csv = true
  file_base = OutputData/USplineElasticity
[]
