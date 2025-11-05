[Mesh]
  [./gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 10
    ny = 10
  []
  [shift_zero_block]
    type = RenameBlockGenerator
    old_block = '0'
    new_block = '1'
    input = gmg
  []
  # [shift_zero_boundary]
  #   type = RenameBoundaryGenerator
  #   old_boundary = '0'
  #   new_boundary = '1'
  #   input = shift_zero_block
  # []
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