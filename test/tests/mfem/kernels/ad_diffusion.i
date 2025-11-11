[Mesh]
  type = MFEMMesh
  file = ../mesh/square.e
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

[Functions]
  [ic_function]
    type = ParsedFunction
    expression = '0.01*sin(pi*x)'
  []
[]

[ICs]
  [ic]
    type = MFEMScalarIC
    variable = concentration
    coefficient = ic_function
  []
  [ic_bnd]
    type = MFEMScalarBoundaryIC
    variable = concentration
    coefficient = ic_function
  []  
[]

[Variables]
  [concentration]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[Solver]
  type = MFEMCGSolver
  l_tol = 1e-4
  l_abs_tol = 0.0
  l_max_its = 500
  print_level = 2
[]

[Executioner]
  type = MFEMSteady
  device = cpu
  use_ad = true
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/ADDiffusion
    vtk_format = ASCII
  []
[]
