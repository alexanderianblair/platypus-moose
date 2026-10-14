# Isogeometric analysis (IGA) of the Laplace equation on a U-spline mesh read from a BEXT file,
#
#   -div(grad(u)) = 0    in the domain,
#               u = g    on the boundary,
#
# where g is the exact solution. The default exact solution is linear, which makes this a patch
# test: the U-spline space contains the geometry map, and hence every linear function. On the
# default mesh, whose weights are all one, the basis functions are polynomials, the default
# quadrature is exact, and so is the discrete solution. On meshes with rational basis functions,
# such as those with exactly represented circular holes, Gauss quadrature of the rational
# integrands is inexact, so the error of the patch test is the quadrature error, which vanishes
# as the quadrature order is raised.

[Mesh]
  type = MFEMUSplineMesh
  file = ../mesh/25_quad.bxt
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [USplineFESpace]
    type = MFEMUSplineFESpace
  []
[]

[Variables]
  [u]
    type = MFEMVariable
    fespace = USplineFESpace
  []
[]

[Functions]
  [exact_solution]
    type = ParsedFunction
    expression = '1 + 2*x + 3*y + 4*z'
  []
[]

[BCs]
  [boundary]
    type = MFEMScalarDirichletBC
    variable = u
    coefficient = exact_solution
  []
[]

[Kernels]
  [diffusion]
    type = MFEMDiffusionKernel
    variable = u
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
    type = MFEML2Error
    variable = u
    function = exact_solution
  []
[]

[Outputs]
  csv = true
  file_base = OutputData/USplinePatchTest
[]
