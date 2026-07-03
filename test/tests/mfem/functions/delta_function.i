[Mesh]
  type = MFEMMesh
  file = ../mesh/square_tri6.e
[]

[Problem]
  type = MFEMProblem
[]

[FESpaces]
  [H1FESpace]
    type = MFEMScalarFESpace
    fec_type = H1
    fec_order = SECOND
  []
[]

[Variables]
  [variable]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[BCs]
  [bc]
    type = MFEMScalarDirichletBC
    variable = variable
  []
[]

[Functions]
  [delta_func]
    type = MFEMDeltaFunction
    position_x = 0.5
    position_y = 0.5
  []
[]

[Kernels]
  [diff]
    type = MFEMDiffusionKernel
    variable = variable
  []
  [point_source]
    type = MFEMDomainLFKernel
    variable = variable
    coefficient = delta_func
  []
[]

[Preconditioner]
  [boomeramg]
    type = MFEMHypreBoomerAMG
  []
[]

[Solver]
  type = MFEMHyprePCG
  preconditioner = boomeramg
  l_tol = 1e-16
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

# [Postprocessors]
#   [error]
#     type = MFEML2Error
#     variable = variable
#     function = solution
#   []
# []

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/DeltaFunction
    vtk_format = ASCII
  []
[]
