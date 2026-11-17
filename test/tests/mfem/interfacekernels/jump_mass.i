# Diffusion across an interface of finite conductance. The unit square blocks 'one',
# 0 < x < 1, and 'two', 1 < x < 2, meet on the interface x = 1, across which the flux
# is k [u], with [u] = u_one - u_two the jump in u, and k = 1 the conductance. With
# u = 1 at x = 0 and u = 0 at x = 2 and unit diffusivity, the flux is q = 1 / (2 + 1/k)
# = 1/3 throughout, and
#
#   u_one = 1 - x / 3,   u_two = (2 - x) / 3,   [u] = 1/3,
#
# which first order elements represent exactly. Each side's variable is defined over
# the whole mesh, governed by a kernel on its own block and fixed at zero on the other
# but for the interface.

[Mesh]
  type = MFEMFileMesh
  file = two_block_interface.e
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

[Variables]
  [u_one]
    type = MFEMVariable
    fespace = H1FESpace
  []
  [u_two]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[Functions]
  [jump]
    type = MFEMParsedFunction
    expression = 'u - v'
    symbol_names = 'u v'
    symbol_values = 'u_one u_two'
  []
[]

[Kernels]
  [diffusion_one]
    type = MFEMDiffusionKernel
    variable = u_one
    block = one
  []
  [diffusion_two]
    type = MFEMDiffusionKernel
    variable = u_two
    block = two
  []
[]

[InterfaceKernels]
  [conductance]
    type = MFEMJumpMassInterfaceKernel
    variable = u_one
    neighbor_var = u_two
    boundary = interface
    coefficient = 1
  []
[]

[Constraints]
  [u_one_in_two]
    type = MFEMScalarEssentialConstraint
    variable = u_one
    block = two
    constrain_subdomain_boundary = false
    constrain_exterior_boundary = true
  []
  [u_two_in_one]
    type = MFEMScalarEssentialConstraint
    variable = u_two
    block = one
    constrain_subdomain_boundary = false
    constrain_exterior_boundary = true
  []
[]

[BCs]
  [left]
    type = MFEMScalarDirichletBC
    variable = u_one
    boundary = left
    coefficient = 1
  []
  [right]
    type = MFEMScalarDirichletBC
    variable = u_two
    boundary = right
    coefficient = 0
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
    print_level = 0
  []
  [pcg]
    type = MFEMHyprePCG
    preconditioner = boomeramg
    l_tol = 1e-14
    print_level = 0
  []
[]

[Executioner]
  type = MFEMSteady
[]

[Postprocessors]
  # The integral of [u] over the interface, of unit length: 1/3.
  [interface_jump]
    type = MFEMBoundaryIntegralPostprocessor
    coefficient = jump
    boundary = interface
  []
[]

[VectorPostprocessors]
  # u_one = 5/6 and u_two = 1/6 at the centres of their blocks, and on the interface
  # 2/3 and 1/3.
  [point_sample]
    type = MFEMVariablePointValueSampler
    variables = 'u_one u_two'
    points = '0.5 0.5 0  1 0.5 0  1.5 0.5 0'
  []
[]

[Outputs]
  csv = true
  file_base = OutputData/JumpMass
[]
