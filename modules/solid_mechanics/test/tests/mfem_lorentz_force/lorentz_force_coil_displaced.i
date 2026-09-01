# MFEM sub-application for the two-way coupled form of the Lorentz force example.
#
# This differs from lorentz_force_coil.i in two ways. The mesh is displaced by the
# deformation computed in the parent application, so the conduction problem is
# re-solved on the deformed conductor at every fixed point iteration; and the
# applied magnetic flux density varies through the thickness of the coil bar,
#
#   B(z) = applied_field * (1 + z / field_length) e_y
#
# so a conductor that moves in z samples a different field and feels a different
# Lorentz force. That is the feedback the one-way examples cannot represent.
#
# The parent sends the three displacement components as separate scalars, because
# transfers of MFEM vector variables are not supported. They are reassembled into
# the vector field that MFEMMesh consumes with an MFEMVectorFromScalarsAux. The
# displacement is a total displacement measured from the undeformed mesh, so
# Mesh/displacement_is_total is set; this keeps repeated displacement idempotent
# across fixed point iterations rather than accumulating.

conductivity = 1.0 # Electrical conductivity of the conductor, S/m
coil_voltage = 1.0 # Terminal voltage ramp rate applied across the coil bar, V/s
applied_field = 1.0 # Applied magnetic flux density at z = 0, T
field_length = 0.05 # Length scale over which the applied field varies in z, m

[Mesh]
  type = MFEMMesh
  file = coil_bar.mesh
  displacement = mesh_displacement
  displacement_is_total = true
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
  [H1VectorFESpace]
    type = MFEMVectorFESpace
    fec_type = H1
    fec_order = FIRST
    range_dim = 3
    ordering = vdim
  []
  [HCurlFESpace]
    type = MFEMVectorFESpace
    fec_type = ND
    fec_order = FIRST
  []
  [L2VectorFESpace]
    type = MFEMVectorFESpace
    fec_type = L2
    fec_order = FIRST
  []
  [L2FESpace]
    type = MFEMScalarFESpace
    fec_type = L2
    fec_order = FIRST
  []
[]

[Variables]
  [electric_potential]
    type = MFEMVariable
    fespace = H1FESpace
  []
[]

[AuxVariables]
  # Displacement components received from the mechanics parent application
  [mesh_disp_x]
    type = MFEMVariable
    fespace = H1FESpace
  []
  [mesh_disp_y]
    type = MFEMVariable
    fespace = H1FESpace
  []
  [mesh_disp_z]
    type = MFEMVariable
    fespace = H1FESpace
  []
  [mesh_displacement]
    type = MFEMVariable
    fespace = H1VectorFESpace
  []
  [current_density]
    type = MFEMVariable
    fespace = HCurlFESpace
  []
  [current_density_x]
    type = MFEMVariable
    fespace = L2FESpace
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
  # Evaluated on the displaced mesh, so a conductor that moves in z sees a
  # different field than it did in its undeformed position.
  [applied_b_field]
    type = ParsedVectorFunction
    expression_x = '0'
    expression_y = '${applied_field} * (1 + z / ${field_length})'
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

[AuxKernels]
  # Reassemble the vector the mesh is displaced by from the transferred components
  # Runs on TIMESTEP_BEGIN, before the solve that displaces the mesh, so the mesh
  # is moved by the displacement just received rather than by the previous one.
  # The parent samples this application on its own displaced mesh, so the two
  # must agree: were the mesh to lag, points on the moving surface would fall
  # outside it and the transfer would return its out-of-mesh value.
  [assemble_displacement]
    type = MFEMVectorFromScalarsAux
    variable = mesh_displacement
    component_coefficients = 'mesh_disp_x mesh_disp_y mesh_disp_z'
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [current]
    type = MFEMGradAux
    variable = current_density
    source = electric_potential
    scale_factor = ${fparse -1.0 * conductivity}
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [current_x]
    type = MFEMInnerProductAux
    variable = current_density_x
    first_source_vec = current_density
    second_source_vec = '1 0 0'
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

[Executioner]
  type = MFEMTransient
  dt = 0.1
  num_steps = 10
  device = cpu
[]

[Outputs]
  [ParaViewDataCollection]
    type = MFEMParaViewDataCollection
    file_base = OutputData/LorentzForceCoilDisplaced
    vtk_format = ASCII
  []
[]
