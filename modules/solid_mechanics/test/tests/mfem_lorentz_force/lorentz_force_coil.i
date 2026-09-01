# MFEM sub-application computing the Lorentz force density acting on a
# current-carrying coil bar held in a uniform applied magnetic field.
#
# The conductor occupies 0 <= x <= 1, 0 <= y <= 0.1, 0 <= z <= 0.05. A potential
# difference is applied between the two end faces, so the steady current problem
#
#   div(sigma grad(V)) = 0,   V = V0 * t on x = 0,   V = 0 on x = L
#
# has the exact solution V = V0 * t * (1 - x / L), giving the uniform, linearly
# ramped current density
#
#   J = -sigma grad(V) = (sigma * V0 * t / L) e_x .
#
# Crossing that with the uniform applied field B = B0 e_y gives a Lorentz force
# density that is uniform in space and ramps linearly in time,
#
#   f = J x B = (sigma * V0 * B0 * t / L) e_z .
#
# The three components of f are projected onto scalar L2 auxvariables so that
# they can be picked up by the libMesh mechanics parent application, which does
# not support transfers of MFEM vector variables.

conductivity = 1.0 # Electrical conductivity of the conductor, S/m
coil_voltage = 1.0 # Terminal voltage ramp rate applied across the coil bar, V/s
applied_field = 1.0 # Uniform applied magnetic flux density along y, T

[Mesh]
  type = MFEMMesh
  file = coil_bar.mesh
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
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
  []
  [L2VectorFESpace]
    type = MFEMVectorFESpace
    fec_type = L2
    fec_order = CONSTANT
  []
  [L2FESpace]
    type = MFEMScalarFESpace
    fec_type = L2
    fec_order = CONSTANT
  []
[]

[Variables]
  [electric_potential]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[AuxVariables]
  [current_density]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [lorentz_force]
    type = MFEMVariable
    fespace = L2VectorFESpace
  []
  [lorentz_force_x]
    type = MFEMVariable
    fespace = L2FESpace
  []
  [lorentz_force_y]
    type = MFEMVariable
    fespace = L2FESpace
  []
  [lorentz_force_z]
    type = MFEMVariable
    fespace = L2FESpace
  []
[]

[Functions]
  [terminal_voltage]
    type = ParsedFunction
    expression = '${coil_voltage} * t'
  []
  [applied_b_field]
    type = ParsedVectorFunction
    expression_x = '0'
    expression_y = '${applied_field}'
    expression_z = '0'
  []
[]

[FunctorMaterials]
  [conductor]
    type = MFEMGenericFunctorMaterial
    prop_names = conductivity
    prop_values = ${conductivity}
  []
[]

[Kernels]
  [conduction]
    type = MFEMDiffusionKernel
    variable = electric_potential
    coefficient = conductivity
  []
[]

[BCs]
  [high_terminal]
    type = MFEMScalarDirichletBC
    variable = electric_potential
    boundary = '1'
    coefficient = terminal_voltage
  []
  [low_terminal]
    type = MFEMScalarDirichletBC
    variable = electric_potential
    boundary = '2'
    coefficient = 0.0
  []
[]

# These auxkernels form a chain, each one consuming the variable updated by the
# previous one. The cross and inner products couple to their sources through
# coefficients rather than through declared variable dependencies, so they are
# listed here in the order they must be executed in.
[AuxKernels]
  [current]
    type = MFEMGradAux
    variable = current_density
    source = electric_potential
    scale_factor = ${fparse -1.0 * conductivity}
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [lorentz]
    type = MFEMCrossProductAux
    variable = lorentz_force
    first_source_vec = current_density
    second_source_vec = applied_b_field
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [lorentz_x]
    type = MFEMInnerProductAux
    variable = lorentz_force_x
    first_source_vec = lorentz_force
    second_source_vec = '1 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [lorentz_y]
    type = MFEMInnerProductAux
    variable = lorentz_force_y
    first_source_vec = lorentz_force
    second_source_vec = '0 1 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [lorentz_z]
    type = MFEMInnerProductAux
    variable = lorentz_force_z
    first_source_vec = lorentz_force
    second_source_vec = '0 0 1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Solvers]
  [boomeramg]
    type = MFEMHypreBoomerAMG
  []
  [main]
    type = MFEMHyprePCG
    preconditioner = boomeramg
    l_tol = 1e-14
  []
[]

# The timestep is set by the parent application through the TransientMultiApp.
[Executioner]
  type = MFEMTransient
  dt = 0.1
  num_steps = 10
  device = cpu
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/LorentzForceCoil
    vtk_format = ASCII
  []
[]
