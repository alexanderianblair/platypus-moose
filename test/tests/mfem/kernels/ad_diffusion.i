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
    expression = 'log(cos(a * pi * (x-x0)) / cos(a * pi * (y-y0))) / a'
    symbol_names = 'a x0 y0'
    symbol_values = '1e-2 0. 0.'
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
