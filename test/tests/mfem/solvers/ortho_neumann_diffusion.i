# Diffusion solve on a singular, pure Neumann system using mfem::OrthoSolver.
#
# PDE:  -Laplacian(u) = f  on [0,1]^2
# Exact solution: u = cos(pi x)cos(pi y)
# Forcing:        f = 2 pi^2 cos(pi x)cos(pi y)
# BCs:            homogeneous Neumann on all four sides, satisfied naturally since
#                 grad(u).n = 0 there, so no [BCs] block is needed.
#
# With no essential boundary conditions the operator is singular, its nullspace spanned by
# the constant vector. MFEMOrthoSolver subtracts the mean from both the right hand side
# handed to the wrapped CG solve and the solution it returns, picking out the zero-mean
# solution. The exact solution integrates to zero over the domain, and the node
# distribution is symmetric about x = 0.5 and y = 0.5, so cos(pi x)cos(pi y) sums to zero
# over the nodes as well; the two normalisations therefore agree and the L2 error below
# measures discretisation error alone.

[Mesh]
  type = MFEMFileMesh
  file = ../mesh/square.e
  uniform_refine = 2
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
  [concentration]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[Functions]
  [u_exact]
    type = ParsedFunction
    expression = cos(pi*x)*cos(pi*y)
  []
  [forcing]
    type = ParsedFunction
    expression = 2*pi^2*cos(pi*x)*cos(pi*y)
  []
[]

[Kernels]
  [diff]
    type = MFEMDiffusionKernel
    variable = concentration
  []
  [rhs]
    type = MFEMDomainLFKernel
    variable = concentration
    coefficient = forcing
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
    print_level = 0
  []
  [cg]
    type = MFEMCGSolver
    preconditioner = boomeramg
    l_tol = 1e-12
    l_max_its = 300
    print_level = -1
  []
  [ortho]
    type = MFEMOrthoSolver
    solver = cg
  []
[]

[Executioner]
  type = MFEMSteady
  device = cpu
[]

[Postprocessors]
  [l2_error]
    type = MFEML2Error
    variable = concentration
    function = u_exact
  []
[]

[Outputs]
  [csv]
    type = CSV
    file_base = OutputData/OrthoNeumannDiffusion
  []
[]
