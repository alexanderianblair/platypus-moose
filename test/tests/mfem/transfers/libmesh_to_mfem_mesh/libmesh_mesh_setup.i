[Mesh]
  [./gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 10
    ny = 10
    nz = 10
  []
  [shift_zero_block]
    type = RenameBlockGenerator
    old_block = '0'
    new_block = '1'
    input = gmg
  []  
[]

[Problem]
  type = FEProblem
  solve = false
[]

[Executioner]
  type = Steady
[]


[MultiApps]
  [sub]
    type = FullSolveMultiApp
    input_files = mfem_diffusion.i
    execute_on = 'INITIAL'
  []
[]