!include nldiffusion_common.i

# Nonlinear diffusion solved without ever assembling the Jacobian as a matrix: the Newton update
# is driven by GMRES acting on the matrix-free Jacobian action, preconditioned by its diagonal.
# Intended to be run with a non-legacy Executioner/assembly_level.

[Solvers]
  [jacobi]
    type = MFEMOperatorJacobiSmoother
  []
  [lin]
    type = MFEMGMRESSolver
    preconditioner = jacobi
    print_level = 1
    l_tol = 1e-12
    l_max_its = 1000
  []
  [native_mfem_nl]
    type = MFEMNewtonNonlinearSolver
    max_its = 100
    abs_tol = 1.0e-10
    rel_tol = 1.0e-9
    print_level = 1
  []
[]
